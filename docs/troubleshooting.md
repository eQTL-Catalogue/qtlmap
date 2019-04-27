# kerimoff/qtlmap: Troubleshooting

## Input files not found
Each mandatory argument (input files) accepts only one file. 
If any of the mandatory files will be missing the pipeline will not work and show the corresponding message
```
ERROR ~ Cannot find any genotype vcf file: testdat/genotype_samples.vcf.gzm
```
Note that the pathes can be surronded with quotes but it is not necessary, because pipeline does not support multiple files for any of the input files.


## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/kerimoff/qtlmap) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
