FROM rocker/r-ver:4
# RUN R -q -e 'install.packages("curl")'
# FROM rocker/rstudio:4.2

# RUN R -q -e 'install.packages(c("rmarkdown", "DT"))'

RUN install2.r \
    rmarkdown \
    DT \
    Biostrings \
    && rm -rf /tmp/downloaded_packages \
    && strip /usr/local/lib/R/site-library/*/libs/*.so
# RUN install2.r \
# rmarkdown \
# DT \
# shiny \
# knitr \
# htmltools \
# readr \
# readxl \
# ggplot2 \
# Biostrings