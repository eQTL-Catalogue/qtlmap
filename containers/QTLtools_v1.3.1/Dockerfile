FROM nfcore/base
LABEL authors="Nurlan Kerimov" \
      description="Docker image containing all requirements for QTLTools v1.3.1"

RUN apt-get update && apt-get install -y \
   build-essential \
   libboost-dev \
   libgsl0-dev \
   gfortran \
   zlib1g-dev \
   libcurl4-openssl-dev \
   libblas-dev \
   libboost-all-dev \
   libbz2-dev 

COPY bin/QTLtools /usr/bin/
