# Molecular Simulation of Fluids 

Based on text and source from:
Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation, and Parallel Computing, 2nd Edition by Richard J. Sadus

https://www.elsevier.com/books-and-journals/book-companion/9780323853989

This project containerizes the build and (MPI and CUDA) runtime environment.

To use:
```
$ cd docker
docker$ ./build_image.sh && ./interactive.sh
...
root@6e11566453a4:/repo# mkdir -p build && cd build && cmake .. && make
...
root@6e11566453a4:/repo/build# ls bin/
mcCU_simulation  mcMPI_simulation  mc_simulation  mdCU_simulation  mdMPI_simulation  md_simulation
```

\[
\begin{array}{|l|l|l|}
\hline
\textbf{Program} & \textbf{Description} & \textbf{Notes} \\
\hline
\text{mc} & \text{Serial Monte Carlo program} & \text{Updated to use cmake and common data.} \\
\hline
\text{mcMPI} & \text{MPI Monte Carlo program} & \text{None} \\
\hline
\text{mcCU} & \text{CUDA Monte Carlo program} & \text{None} \\
\hline
\text{mcmd} & \text{Combined serial MC and MD program} & \text{None} \\
\hline
\text{md} & \text{Serial molecular dynamics program} & \text{None} \\
\hline
\text{mdMPI} & \text{MPI molecular dynamics program} & \text{None} \\
\hline
\text{mdCU} & \text{CUDA molecular dynamics program} & \text{None} \\
\hline
\end{array}
\]



