FROM rocker/r-ver:4

RUN echo 'options(BioC_mirror = "https://packagemanager.posit.co/bioconductor")' >>"${R_HOME}/etc/Rprofile.site"
RUN install2.r --error \
    rmarkdown \
    DT \
    BiocManager

RUN $R_HOME/site-library/littler/examples/installBioc.r --error Biostrings 

# RUN R -q -e 'install.packages(c("rmarkdown", "Biostrings"))'
