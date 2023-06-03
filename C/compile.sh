export OMP_NUM_THREADS=10
gcc main.c rede.c calc.c mtwister.c infect.c CM.c SBM.c LCM.c -o main -lm -O3 -fopenmp
./main 42
#for i in {6..100}; do
#  for j in {0..100}; do
#    ./main 42 100 $j $i
#  done
#don