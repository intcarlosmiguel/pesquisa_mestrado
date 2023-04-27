export OMP_NUM_THREADS=10
gcc main.c rede.c calc.c mtwister.c infect.c estagio.c CM.c SBM.c LCM.c -o main -lm -O3 -fopenmp
time for counter in {1..80}; do ./main 42 500 $counter; done
#time ./main 42 1
#for i in {6..100}; do
#  for j in {0..100}; do
#    ./main 42 100 $j $i
#  done
#don