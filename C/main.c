#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"
#include "CM.h"
#include "SBM.h"
#include "LCM.h"
#include "infect.h"

void teste_resto_ligacoes(){
    struct Graph G;
    for(int i = 0;i < 1; i++) G = local_configuration_model(2025,0,8,0);
}

int main(int argc,char *argv[ ]){
    //generate_configuration_model((double) 0.5,1000);
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
    local_configuration_model(160029,0.0,seed,2);
    //teste_resto_ligacoes();
    /*int seed = atoi(argv[1]);
    int redes = atoi(argv[2]);
    int vacina = atoi(argv[3]);
    for(int i = 50;i < 51; i++)generate_infect(seed+i, redes,(double)i/100,vacina);*/
}
