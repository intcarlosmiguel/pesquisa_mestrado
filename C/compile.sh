gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

for i in {61..100}; do
    ./main $i 10000 300 4 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
