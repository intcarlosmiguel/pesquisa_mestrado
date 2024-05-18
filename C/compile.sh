gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
seed=62
for ((i = 50; i <= 100; i += 25)); do
    ./main $seed 10000 400 22 $i 1
    ((seed += 400))
    ./main $seed 10000 400 24 $i 1
    ((seed += 400))
    ./main $seed 10000 400 26 $i 1
    ((seed += 400))
    ./main $seed 10000 400 28 $i 1
    ((seed += 400))
    ./main $seed 10000 400 30 $i 1
    ((seed += 400))
done