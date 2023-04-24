gcc main.c rede.c calc.c mtwister.c infect.c CM.c SBM.c LCM.c -o main -lm -O3 -fopenmp
for i in {0..90..10}; do
  gnome-terminal --title="Meu Programa" -- bash -c "time ./main 42 500 $i; exec bash"
done