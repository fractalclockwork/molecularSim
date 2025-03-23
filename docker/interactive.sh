#!/bin/sh

if [ "$1" = "cuda" ]; then
    docker run --rm -it --gpus all -v $(pwd)/..:/repo mole_sim_cuda
else
    docker run --rm -it -v $(pwd)/..:/repo mole_sim_cuda
    #docker run --rm -it -v $(pwd)/..:/repo mole_sim
fi
