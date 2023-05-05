#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"
//#include "CM.h"
//#include "SBM.h"
#include "LCM.h"
#include "infect.h"
//#include "estagio.h"

int main(/*void (*ptr_funcao)(int** viz,int N)*/){
    //generate_configuration_model((double) 0,1);
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
    //struct Graph G = local_configuration_model(2029,0.0,1);
    //printf("%d\n",G.viz[0][0]);
    //ptr_funcao(G.viz,2029);
    //if(numero == 90) for(int i =numero;i<=numero+10;i++) generate_infect(seed, redes,(double)i/100);
    //else for(int i =numero;i<numero+10;i++) generate_infect(seed, redes,(double)i/100);
    //for(int i =0;i<=100;i++) generate_infect(seed, redes,(double)i/100);
    //for(int i =0;i<=100;i++) generate_infect(seed, redes,(double)i/100,tau);
    generate_infect(42, 1,0.0);
}
