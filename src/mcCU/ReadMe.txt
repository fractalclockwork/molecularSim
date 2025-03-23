This directory contains the source, header and data files required for
the MC CUDA program. A simple compilation is:

nvcc --fmad=false -o mcCU.exe *.cpp *.cu

The use of --fmad=false is important for the correct execution of the CUDA code. 
Further details are found is the software user's manual, which is part of the book.
Check the book's website for any updates.
