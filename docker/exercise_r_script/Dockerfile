FROM r-base:4.2.3

RUN apt-get update
RUN apt-get install -y \
    r-cran-optparse \
    r-bioc-biomart

COPY search_biomart_datasets.R /opt 

ENV PATH=/opt:$PATH

# note that if you want to be able to combine the two
# both ENTRYPOINT and CMD need to written in the exec form
ENTRYPOINT ["search_biomart_datasets.R"]

# default option (if positional arguments are not specified)
CMD ["--pattern", "mouse"]
