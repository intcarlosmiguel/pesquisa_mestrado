gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++


./main 1 10000 400 25 0 0
./main 11 10000 400 25 25 0
./main 111 10000 400 25 50 0
./main 1111 10000 400 25 75 0
./main 11111 10000 400 25 100 0