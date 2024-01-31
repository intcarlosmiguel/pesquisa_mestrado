gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3

./main 5 10000 1000 2 25 0 0
./main 5 10000 1000 2 50 0 0
./main 5 10000 1000 2 75 0 0
./main 5 10000 1000 2 100 0 0