FROM bioconductor/bioconductor_docker:RELEASE_3_13
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for susieR fine mapping"

RUN R -e "BiocManager::install(c('remotes', 'dplyr', 'optparse', 'readr', 'purrr', 'tidyr','assertthat', 'matrixstats', 'expm', 'data.table','stringr','SummarizedExperiment','Rsamtools','susieR', 'zoo'))"
RUN R -e "remotes::install_github('kauralasoo/eQTLUtils@v20.04.1')"