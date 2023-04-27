#export OMP_NUM_THREADS=10
gcc -fpic -shared calc.c mtwister.c infect.c -o main.so -lm -O3 #-fopenmp
julia main.jl
#time for counter in {0..100}; do ./main 42 500 $counter; done
#time ./main 42 1 100
#for i in {6..100}; do
#  for j in {0..100}; do
#    ./main 42 100 $j $i
#  done
#done