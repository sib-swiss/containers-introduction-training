# Use Jupyter with Python 3.13 as base image
FROM jupyter/base-notebook:x86_64-ubuntu-22.04

# Add metadata to image
LABEL org.opencontainers.image.created="2026-06-02" \
    org.opencontainers.image.authors="Antonin Thiébaut" \
    org.opencontainers.image.description="Container with PyDESeq2 and JupyterLab"

# Install pydeseq2 and friends
RUN pip install pydeseq2

# Set working directory
WORKDIR /opt

# Copy Python script inside container
COPY test_pydeseq2.py .

# Add /opt to PATH so we can execute script directly
ENV PATH=/opt:$PATH
