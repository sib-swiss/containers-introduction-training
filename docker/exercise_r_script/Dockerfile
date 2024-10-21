FROM rocker/r2u:jammy

LABEL org.opencontainers.image.created="2024-10-12" \
    org.opencontainers.image.authors="Geert van Geest" \
    org.opencontainers.image.description="Container with DESeq2 and friends"

# install DESeq2 and friends
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

WORKDIR /opt

COPY test_deseq2.R .

ENV PATH=/opt:$PATH

# note that if you want to be able to combine the two
# both ENTRYPOINT and CMD need to written in the exec form
ENTRYPOINT ["test_deseq2.R"]

# default option (if positional arguments are not specified)
CMD ["--rows", "100"]
