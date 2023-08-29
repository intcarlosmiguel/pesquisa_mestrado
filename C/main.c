#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <igraph.h>

#include <omp.h>

#include "bib/mtwister.h"
#include "bib/calc.h"
#include "bib/rede.h"
#include "bib/CM.h"
#include "bib/SBM.h"
#include "bib/LCM.h"
#include "bib/infect.h"

int main(int argc,char *argv[ ]){
    //igraph_t Grafo = local_configuration_model(2029, 1.0,738278);
    //generate_configuration_model((double) 1.0,500);
    //generate_SBM_p_model(1,1,0);
    //generate_local_configuration_model(0,1000);
    int seed = atoi(argv[1]);
    int N = atoi(argv[2]);
    int redes = atoi(argv[3]);
    int vacina = atoi(argv[4]);
    int prob = atoi(argv[5]);
    int freq = atof(argv[6]);
    generate_infect(N,(double)prob/100,seed, redes,(double) freq/100,vacina);
    /* int T = atoi(argv[1]);
    ///int seed = atoi(argv[1]);
    //int redes = atoi(argv[2]);
    //int vacina = atoi(argv[3]);
    //double prop = atoi(argv[4]);
    //generate_infect(2029,seed, redes,prop/100,vacina);
    generate_infect(2029,seed, redes,prop/100,vacina);
    //generate_configuration_model((double) 0.0,1);
    //local_configuration_model(7286,1.0,232,2);
    //teste_resto_ligacoes();
    /*int seed = atoi(argv[1]);
    int redes = atoi(argv[2]);
    int vacina = atoi(argv[3]);*/
    
}
