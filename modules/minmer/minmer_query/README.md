# minmer_query process testing:

This process querys minmer sketches against pre-computed RefSeq (Mash, k=21) and GenBank (Sourmash, k=21,31,51)

## About testing this process:

Using DSL2 each module can be tested separately, using a test workflow inside the process.nf file, testing requires 3 itens:  
- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`

## How to test it:

$ nextflow run minmer_query.nf -entry test -params-file test_params.yaml -profile test
