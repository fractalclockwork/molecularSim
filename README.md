# Molecular Simulation of Fluids 

Based on text and source from:
Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation, and Parallel Computing, 2nd Edition by Richard J. Sadus

https://www.elsevier.com/books-and-journals/book-companion/9780323853989

This project containerizes the build and (MPI and CUDA) runtime environment.

To use:
```
$ cd docker
docker$ ./build_image.sh && ./interactive.sh # No GPU 
docker$ ./build_image.sh && ./interactive.sh cuda # Requires nvidia-container-toolkit 
...
root@6e11566453a4:/repo# mkdir -p build && cd build && cmake .. && make
...
root@6e11566453a4:/repo/build# ls bin/
mcCU_simulation  mcMPI_simulation  mc_simulation  mdCU_simulation  mdMPI_simulation  md_simulation
```

| Program | Description                         | Notes                               |
|---------|-------------------------------------|-------------------------------------|
| mc      | Serial Monte Carlo program         | None                                | 
| mcMPI   | MPI Monte Carlo program            | None                                |
| mcCU    | CUDA Monte Carlo program           | None                                |
| mcmd    | Combined serial MC and MD program  | None                                |
| md      | Serial molecular dynamics program  | None                                |
| mdMPI   | MPI molecular dynamics program     | None                                |
| mdCU    | CUDA molecular dynamics program    | None                                |
