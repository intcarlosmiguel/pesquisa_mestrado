gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++

seed=90
for i in {0..100}; do
    ./main $seed $i &
    if (($i % 10 == 0)); then
        wait
    fi
    ((seed += 1))
done
