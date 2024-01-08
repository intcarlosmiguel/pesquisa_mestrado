gcc -I/home/miguel/Downloads/igraph-0.10.4/build/include -I/home/miguel/Downloads/igraph-0.10.4/include main.c bib/rede.c bib/calc.c bib/mtwister.c bib/infect.c bib/CM.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3
for i in {1..100}; do
    ./main $i 7189 500 0 0 $i 1 &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 7189 500 0 0 $i 0 &
    if (($i % 10 == 0)); then
        wait
    fi
done