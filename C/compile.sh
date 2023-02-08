export OMP_NUM_THREADS=10
gcc main.c rede.c calc.c mtwister.c CM.c SBM.c LCM.c -o main -lm -O3 -fopenmp
./main 1000 2