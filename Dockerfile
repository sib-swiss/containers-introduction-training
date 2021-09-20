FROM ubuntu
RUN apt-get update
RUN apt-get install figlet
CMD ["/bin/sh", "-c", "figlet 'My image works!'"]
