FROM rocker/rstudio:4

RUN apt-get update && \
    apt-get install -y libz-dev

RUN install2.r \
    BiocManager

RUN R -q -e 'BiocManager::install("biomaRt")'
