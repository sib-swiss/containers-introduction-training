# # Use Rocker R with Ubuntu 22 as base image
FROM rocker/r2u:jammy

# Add metadata to image
LABEL org.opencontainers.image.created="2024-10-12" \
    org.opencontainers.image.authors="Geert van Geest" \
    org.opencontainers.image.description="Container with DESeq2 and friends"

# Install DESeq2 and friends
RUN install2.r \
    DESeq2 \
    optparse \
    apeglm \
    IHW \
    limma \
    data.table \
    ggrepel \
    pheatmap \
    stringr 

# Set working directory
WORKDIR /opt

# Copy R script inside container
COPY test_deseq2.R .

# Add /opt to PATH so we can execute script directly
ENV PATH=/opt:$PATH

# Note that if you want to be able to combine the two, both ENTRYPOINT and CMD
# need to written in the exec form
ENTRYPOINT ["test_deseq2.R"]

# Default option (if positional arguments are not specified)
CMD ["--rows", "100"]
