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
    //igraph_t Grafo = local_configuration_model(7286, 0.0,0);
    //generate_configuration_model((double) 1.0,500);
    //generate_SBM_p_model(1,1,0);
    //generate_local_configuration_model(0,1000);
    /* int T = atoi(argv[1]);
    int model = atoi(argv[2]);
    switch (model){
    case 0: // Modelo de Configuração
        
        //generate_configuration_model((double) 0,T);
        generate_configuration_model((double) 1,T);
        //generate_configuration_model((double) 1.0,T);
        break;
    case 1: // Modelo SBM
        generate_SBM_p_model(T,model,0);
        break;
    case 2: // Modelo SBM probabilidade

        //for(int i = 0;i<= 100;i++) generate_SBM_p_model(T,model,(double) i/100);
        generate_SBM_p_model(T,model,0);
        break;

    case 3:

        generate_local_configuration_model(1,T,0);
        //generate_local_configuration_model(0.5,T,1);
        //generate_local_configuration_model(1,T,1);
        break;
        
    default:
        break;
    } */
    int seed = atoi(argv[1]);
    generate_infect(2029,seed, 1,0.0,-1);

    //generate_configuration_model((double) 0.0,1);
    //local_configuration_model(7286,1.0,232,2);
    //teste_resto_ligacoes();
    /*int seed = atoi(argv[1]);
    int redes = atoi(argv[2]);
    int vacina = atoi(argv[3]);
    
    for(int i = 50;i < 51; i++)generate_infect(seed+i, redes,(double)i/100,vacina);*/
}
