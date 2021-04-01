# estimate_genome_size process testing:

This process estimate a genome size using Mash.

## About testing this process:

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`

## How to test it:

$ nextflow run estimate_genome_size.nf -entry test -params-file test_params.yaml -profile test
