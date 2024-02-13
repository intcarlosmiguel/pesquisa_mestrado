#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <igraph.h>

#include <omp.h>

#include "bib/mtwister.h"
#include "bib/calc.h"
//#include "bib/rede.h"
//#include "bib/CM.h"
#include "bib/SBM.h"
#include "bib/LCM.h"
#include "bib/infect.h"

int main(int argc,char *argv[ ]){
    //double avg;
    //for(int i =0; i< 100; i++)SBM_p(7189,i,0.0);
    //double prob = atof(argv[1]);
    //generate_local_configuration_model(0.0 ,100,43322);
    //generate_local_configuration_model(0.25,100,4242);
    //generate_local_configuration_model(0.5,100,4242);
    //generate_local_configuration_model(0.75,100,4242);
    //generate_local_configuration_model(1.,1,424242);

    int seed = atoi(argv[1]);
    double N = atof(argv[2]);
    int redes = atoi(argv[3]);
    int vacina = atoi(argv[4]);
    double prob = atof(argv[5])/100;
    double fracao_vacinados = atof(argv[6])/100;
    bool weight = atoi(argv[7]);
    generate_infect(N,prob,seed, redes,fracao_vacinados,vacina,weight);
    
}
