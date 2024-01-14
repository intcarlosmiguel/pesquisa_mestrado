gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

for i in {41..100}; do
    ./main $i 10000 200 4 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 10000 200 3 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 10000 200 2 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 10000 200 1 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 10000 200 0 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 