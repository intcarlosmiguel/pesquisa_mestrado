export OMP_NUM_THREADS=10
gcc main.c rede.c calc.c mtwister.c infect.c CM.c SBM.c LCM.c -o main -lm -O3 -fopenmp
for counter in {0..200}; do  time ./main 42 100 $counter; done
#time ./main 42 100 0