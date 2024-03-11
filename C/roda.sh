gcc -I/home/miguel/igraph/build/include -I/home/miguel/igraph/include main.c bib/calc.c bib/mtwister.c bib/infect.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3 -lstdc++ -Wall -Werror


time ./main 47 10000 200 12 0 1
time ./main 48 10000 200 13 0 1
time ./main 49 10000 200 14 0 1
time ./main 50 10000 200 15 0 1
time ./main 54 10000 200 2 0 1
time ./main 64 10000 200 3 0 1
time ./main 464 10000 200 4 0 1
time ./main 474 10000 200 6 0 1
time ./main 4982 10000 200 0 0 1
time ./main 42323 10000 200 0 100 1
time ./main 400 10000 200 1 100 1
time ./main 410 10000 200 5 100 1
time ./main 420 10000 200 7 100 1
time ./main 430 10000 200 8 100 1
time ./main 440 10000 200 9 100 1
time ./main 450 10000 200 10 100 1
time ./main 460 10000 200 11 100 1
time ./main 470 10000 200 12 100 1
time ./main 480 10000 200 13 100 1
time ./main 490 10000 200 14 100 1
time ./main 500 10000 200 15 100 1
time ./main 540 10000 200 2 100 1
time ./main 640 10000 200 3 100 1
time ./main 4640 10000 200 4 100 1
time ./main 4740 10000 200 6 100 1