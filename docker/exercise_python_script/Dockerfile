FROM python:3.9.16

LABEL org.opencontainers.image.created="2022-04-12" \
    org.opencontainers.image.authors="Geert van Geest" \
    org.opencontainers.image.description="Great container for getting all dates in a week! \
    You will never use a calender again"

RUN pip install pandas 

WORKDIR /opt

COPY daterange.py .

ENV PATH=/opt:$PATH

# note that if you want to be able to combine the two
# both ENTRYPOINT and CMD need to written in the exec form
ENTRYPOINT ["daterange.py"]

# default option (if positional arguments are not specified)
CMD ["--date", "20220226"]