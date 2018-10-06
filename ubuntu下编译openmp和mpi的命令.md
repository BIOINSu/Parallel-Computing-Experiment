## OpenMP

g++ -c test.cpp -o test.o -fopenmp

g++ test.o -o test -fopenmp -lpthread



## MPI

mpicxx -o test test.cpp 

mpirun -np 8 ./test 