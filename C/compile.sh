gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=6
for ((i = 0; i <= 100; i += 25)); do
    ./main $seed 10000 400 21 $i 1
    ((seed += 400))
    ./main 88 10000 400 23 $i 1
    ((seed += 400))
    ./main 888 10000 400 25 $i 1
    ((seed += 400))
    ./main 8888 10000 400 27 $i 1
    ((seed += 400))
    ./main 88888 10000 400 29 $i 1
    ((seed += 400))
done