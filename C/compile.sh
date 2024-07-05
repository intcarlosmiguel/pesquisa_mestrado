gcc -I/home/carlos/igraph/build/include -I/home/carlos/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++
./main 2 10000 400 16 0 1
./main 4 10000 400 17 0 1
./main 8 10000 400 18 0 1
./main 10 10000 400 19 0 1
./main 12 10000 400 20 0 1