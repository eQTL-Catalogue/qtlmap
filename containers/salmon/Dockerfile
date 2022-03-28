FROM nfcore/base:1.8
LABEL authors="Nurlan Kerimov, Kaur Alasoo" \
      description="Docker image containing all requirements for the eQTL-Catalogue/rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/eqtl-catalogue-22.03.01/bin:$PATH

