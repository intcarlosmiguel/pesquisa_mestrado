

gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
echo "Compilou!"
seed=2
for ((i = 25; i <= 100; i += 25)); do

    ./main $seed 10000 400 0 $i 0
    ((seed += 400))
done