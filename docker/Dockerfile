FROM ubuntu:latest
RUN apt-get update && apt-get install -y \
    vim \
    build-essential \
    cmake \
    libopenmpi-dev \
    openmpi-bin

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
