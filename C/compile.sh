gcc -I/home/miguel/Downloads/igraph-0.10.4/build/include -I/home/miguel/Downloads/igraph-0.10.4/include main.c bib/rede.c bib/calc.c bib/mtwister.c bib/infect.c bib/CM.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3


./main 0 7189 1 0 100 0
