gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

seed=203
for i in {21..100}; do
    ./main $seed 10000 300 7 100 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 0 100 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 1 100 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 

for i in {1..100}; do
    ./main $seed 10000 300 14 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 13 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 12 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 11 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 9 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 8 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 7 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 0 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 300 1 100 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 