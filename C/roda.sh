gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++ -Wall -Werror

./main 42323 10000 200 0 100 0
./main 400 10000 200 1 100 0
./main 410 10000 200 5 100 0
./main 420 10000 200 7 100 0
./main 430 10000 200 8 100 0
./main 440 10000 200 9 100 0
./main 450 10000 200 10 100 0
./main 460 10000 200 11 100 0
./main 470 10000 200 12 100 0
./main 480 10000 200 13 100 0
./main 490 10000 200 14 100 0
./main 500 10000 200 15 100 0
./main 540 10000 200 2 100 0
./main 640 10000 200 3 100 0
./main 4640 10000 200 4 100 0
./main 4740 10000 200 6 100 0