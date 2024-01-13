gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3


for i in {1..100}; do
    ./main $i 2028 800 11 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 12 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 13 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 14 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 