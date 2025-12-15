FROM ubuntu:20.04

RUN apt-get update -y && \
    apt-get install -y gcc make libopenblas-dev libopenblas0 && \
    rm -rf /var/lib/apt/lists/*


WORKDIR /app
COPY . .

ENV HOSTNAME=ambre
RUN mkdir bin && \
  make all && \
  make run
