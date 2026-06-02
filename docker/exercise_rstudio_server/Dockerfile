# Use Rocker Rstudio v4 as base image
FROM rocker/rstudio:4

# Install zlib development library (required for biomaRt)
RUN apt-get update && \
    apt-get install -y libz-dev

# Install BiocManager
RUN install2.r \
    BiocManager

# Install biomaRt using BiocManager
RUN R -q -e 'BiocManager::install("biomaRt")'
