FROM rocker/r-ubuntu

RUN apt-get update \
&& apt-get install -y --no-install-recommends \
r-cran-rmarkdown \
r-cran-dt \
r-bioc-biostrings \
&& rm -rf /var/lib/apt/lists/*