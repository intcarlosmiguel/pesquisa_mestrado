gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

seed=33
for i in {1..100}; do
    ./main $seed 10000 300 14 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 13 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 12 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 11 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 9 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 8 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
for i in {1..100}; do
    ./main $seed 10000 200 7 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
clear 
./main 4 10000 500 6 0 0 1
./main 43 10000 500 6 25 0 1
./main 433 10000 500 6 50 0 1
./main 4333 10000 500 6 75 0 1
./main 43333 10000 500 6 100 0 1