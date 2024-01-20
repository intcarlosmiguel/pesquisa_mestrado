gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

for i in {1..100}; do
    ./main $i 7189 200 14 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 13 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 12 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 11 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 10 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 0 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 1 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 9 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 8 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 7 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 200 6 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 