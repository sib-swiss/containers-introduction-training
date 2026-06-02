# Use Rocker Rstudio v4 as base image
FROM rocker/rstudio:4

# Add metadata to image
LABEL org.opencontainers.image.created="2024-10-12" \
    org.opencontainers.image.authors="Geert van Geest" \
    org.opencontainers.image.description="Container with DESeq2 and friends"

# Install system dependencies (required for DEseq2 and IHW)
RUN apt update && apt install -y \
    fonts-noto \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    libxml2-dev \
    tzdata \
    && apt clean && apt autoremove && apt purge \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN install2.r \
    BiocManager

# Install DESeq2 and friends using BiocManager
RUN R -q -e "BiocManager::install(c( \
  'DESeq2', \
  'optparse', \
  'apeglm', \
  'IHW', \
  'limma', \
  'data.table', \
  'ggrepel', \
  'pheatmap', \
  'stringr' \
  ))"

# Set working directory
WORKDIR /opt

# Copy R script inside container
COPY test_deseq2.R .

# Add /opt to PATH so we can execute script directly
ENV PATH=/opt:$PATH
