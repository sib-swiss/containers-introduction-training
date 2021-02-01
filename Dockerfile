FROM ubuntu
RUN apt-get update
RUN apt-get install figlet
CMD figlet 'My image works!'
# CMD ["/bin/bash", "-c", "figlet My image works!"]
