export OMP_NUM_THREADS=10
gcc -I/home/miguel/Downloads/igraph-0.10.4/build/include -I/home/miguel/Downloads/igraph-0.10.4/include main.c bib/rede.c bib/calc.c bib/mtwister.c bib/infect.c bib/CM.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -O3 -fopenmp
./main 42
#for i in {6..100}; do
#  for j in {0..100}; do
#    ./main 42 100 $j $i
#  done
#don