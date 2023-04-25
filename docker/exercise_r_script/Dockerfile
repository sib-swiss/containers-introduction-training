FROM r-base:4.2.3

# LABEL org.opencontainers.image.created="2023-04-12" \
#     org.opencontainers.image.authors="Geert van Geest" \
#     org.opencontainers.image.description="Container to search ensembl datasets with biomart"

RUN apt-get update \
    && apt-get install -y \
    r-cran-optparse \
    r-bioc-biomart

WORKDIR /opt

COPY search_biomart_datasets.R .

ENV PATH=/opt:$PATH

# note that if you want to be able to combine the two
# both ENTRYPOINT and CMD need to written in the exec form
ENTRYPOINT ["search_biomart_datasets.R"]

# default option (if positional arguments are not specified)
CMD ["--pattern", "mouse"]
