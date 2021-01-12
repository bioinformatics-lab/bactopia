#!/usr/bin/env bash
# update-docker
#
# Automate the building of Bactopia related Docker containers
REPOSITORY="quay.io"
VERSION=1.5.5
CONTAINER_VERSION="${VERSION%.*}.x"

function docker_build {
    recipe=$1
    image=$2
    latest=${3:-0}

    echo "Working on ${recipe}"
    docker build --rm -t ${image} -f ${recipe} .
    docker push ${image}

    if [[ "${latest}" != "0" ]]; then
        docker tag ${image} ${latest}
        docker push ${latest}
    fi
}

if [[ $# == 0 ]]; then
    echo ""
    echo "build-containers.sh BACTOPIA_DIR OUTPUT_DIR"
    echo ""
    echo "Example Command"
    echo "build-containers.sh /home/bactopia/bactopia"
    echo ""
    exit
fi

BACTOPIA_DIR=$1
if [ -z  ${BACTOPIA_DIR} ]; then
    echo "Got ${#} arguement"
    echo "Must give the path to Bactopia repository"
    exit 1
fi

MAJOR_VERSION=${3:-"0"}

# Build Bactopia containers
docker_build Dockerfile bactopia/bactopia:${VERSION} ${REPOSITORY}/bactopia/bactopia:latest

# Build Docker
docker_build Dockerfile bactopia/bactopia:${CONTAINER_VERSION} bactopia/bactopia:latest
for recipe in $(ls "${BACTOPIA_DIR}/containers/docker" | grep ".Dockerfile"); do
    recipe_path="${BACTOPIA_DIR}/containers/docker/${recipe}"
    recipe_name=$(echo ${recipe} | sed 's/.Dockerfile//')
    recipe_image="${REPOSITORY}/bactopia/${recipe_name}:${CONTAINER_VERSION}"
    docker_build ${recipe_path} ${recipe_image}
done

# Build Bactopia Tools containers
for tool in $(ls "${BACTOPIA_DIR}/tools"); do
    recipe_path="${BACTOPIA_DIR}/tools/${tool}"
    if [ -f "${BACTOPIA_DIR}/tools/${tool}/environment-linux.yml" ]; then
        docker_file="${recipe_path}/Dockerfile"
        docker_image="${REPOSITORY}/bactopia/tools-${tool}:${CONTAINER_VERSION}"
        docker_build ${docker_file} ${docker_image}
    fi
done
