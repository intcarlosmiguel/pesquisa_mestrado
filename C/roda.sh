gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
echo "Compilou!"
seed=2
for ((i = 0; i <= 100; i += 1)); do

    ./main $seed 10000 100 $i
    ((seed += 400))
done