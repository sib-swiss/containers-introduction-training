FROM eddelbuettel/r2u:jammy

# RUN install2.r biomaRt optparse
RUN R -q -e 'install.packages("biomaRt", "optparse")'

WORKDIR /opt

COPY search_biomart_datasets.R .

ENV PATH=/opt:$PATH

# note that if you want to be able to combine the two
# both ENTRYPOINT and CMD need to written in the exec form
ENTRYPOINT ["search_biomart_datasets.R"]

# default option (if positional arguments are not specified)
CMD ["--pattern", "mouse"]
