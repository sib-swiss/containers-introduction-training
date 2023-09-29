FROM rocker/rstudio:4

RUN apt-get update && \
    apt-get install -y libz-dev

# install R packages on CRAN with install2.r
RUN install2.r \
    --error \
    optparse \
    BiocManager

# install bioconductor packages
RUN R -q -e 'BiocManager::install("DESeq2") '
