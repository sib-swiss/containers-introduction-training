# FROM rocker/r-bspm:focal
FROM rocker/r-bspm:testing

RUN install2.r --error --skipinstalled \
    rmarkdown \
    DT \
    shiny \
    knitr \
    htmltools \
    readr \
    readxl \
    ggplot2 \
    Biostrings \
    && rm -rf /var/lib/apt/lists/*
# # 9.8s
# RUN R -q -e 'install.packages(c("rmarkdown", "DT", "Biostrings"))' \
#     && rm -rf /var/lib/apt/lists/*

# FROM rocker/r-bspm:testing
# RUN R -q -e 'install.packages("curl")' \
#     && rm -rf /var/lib/apt/lists/*