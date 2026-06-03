# Use Python 3.13 as base image
FROM python:3.13

# Add metadata to image
LABEL org.opencontainers.image.created="2026-06-02" \
    org.opencontainers.image.authors="Antonin Thiébaut" \
    org.opencontainers.image.description="Container with PyDESeq2 and friends"

# Install pydeseq2 and friends
RUN pip install pydeseq2

# Set working directory
WORKDIR /opt

# Copy Python script inside container
COPY test_pydeseq2.py .

# Add /opt to PATH so we can execute script directly
ENV PATH=/opt:$PATH

# Specify entrypoint
# Note that if you want to be able to combine the two, both ENTRYPOINT and CMD
# need to written in the exec form
ENTRYPOINT ["test_pydeseq2.py"]

# Default option (if positional arguments are not specified)
CMD ["--rows", "100"]
