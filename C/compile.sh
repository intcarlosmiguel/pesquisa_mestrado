gcc -I/home/miguel/Downloads/igraph-0.10.4/build/include -I/home/miguel/Downloads/igraph-0.10.4/include main.c bib/rede.c bib/calc.c bib/mtwister.c bib/infect.c bib/CM.c bib/SBM.c bib/LCM.c -o main -lm -Ibib -ligraph -fopenmp -O3


for i in {1..100}; do
    ./main $i 2028 800 5 0 $i &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 5 50 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 5 100 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 6 0 $i &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 6 50 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 6 100 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 7 0 $i &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 7 50 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 7 100 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 8 0 $i &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 8 50 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 8 100 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 9 0 $i &
    if (($i % 10 == 0)); then
        wait
    fi
done
clear 
for i in {1..100}; do
    ./main $i 2028 800 9 50 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear
for i in {1..100}; do
    ./main $i 2028 800 9 100 $i & 
    if (($i % 10 == 0)); then
        wait
    fi
done
clear