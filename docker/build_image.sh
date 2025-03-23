#!/bin/sh

if [ "$1" = "cuda" ]; then
    docker build -f Dockerfile.cuda -t mole_sim_cuda .
else
    docker build -f Dockerfile.cuda -t mole_sim_cuda .
    #docker build -t mole_sim .
fi
