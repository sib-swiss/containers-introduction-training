FROM python:3.9.16

RUN pip install pandas 

COPY daterange.py /opt 

ENV PATH=/opt:$PATH

# note that if you want to be able to combine the two
# both ENTRYPOINT and CMD need to written in the exec form
ENTRYPOINT ["daterange.py"]

# default option (if positional arguments are not specified)
CMD ["--date", "20220226"]
