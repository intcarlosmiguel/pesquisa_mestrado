gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

for i in {1..100}; do
    ./main $i 7189 300 5 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 7189 300 5 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
