gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

seed=203
for i in {91..100}; do
    ./main $seed 10000 200 2 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 3 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
