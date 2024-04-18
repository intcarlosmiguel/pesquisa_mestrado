#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <igraph.h>

#include <omp.h>

#include "bib/mtwister.h"
#include "bib/calc.h"
//#include "bib/CM.h"
#include "bib/SBM.h"
#include "bib/LCM.h"
//#include "bib/define.h"
#include "bib/infect.h"



int main(int argc,char *argv[ ]){
    clock_t inicio = clock();
    //int seed = atoi(argv[1]);
    //double prob = atof(argv[2])/100;
    //generate_local_configuration_model(prob ,100,seed);

    int seed = atoi(argv[1]);
    double N = atof(argv[2]);
    uint16_t redes = atoi(argv[3]);
    uint8_t vacina = atoi(argv[4]);
    double prob = atof(argv[5])/100;
    bool weight = atoi(argv[6]);
    generate_infect(N,prob,seed, redes,vacina,weight);
    clock_t fim = clock();

    // Calcula a diferença e converte para segundos
    double tempo_gasto = (double)(fim - inicio) / CLOCKS_PER_SEC;

    // Exibe o tempo de execução
    printf(" Tempo gasto: %f minutos\n", tempo_gasto/60/10);

    
}
