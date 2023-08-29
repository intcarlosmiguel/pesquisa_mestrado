export OMP_NUM_THREADS=10
gcc -I/home/miguel/Downloads/igraph-0.10.4/build/include -I/home/miguel/Downloads/igraph-0.10.4/include main.c bib/rede.c bib/calc.c bib/mtwister.c bib/infect.c bib/CM.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

for i in {1..100}; do
    ./main $i 2029 700 2 0 $i &
    if (($i % 10 == 0)); then
        wait
    fi
done
for i in {1..100}; do
    ./main $i 2029 700 2 50 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
for i in {1..100}; do
    ./main $i 2029 700 2 100 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done