gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++

seed=54
for i in {2..20}; do
    ./main $seed 10000 300 $i 100 1
    ((seed += 100))
done