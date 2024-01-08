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
    //SBM_p(7189,true,0);
    //double prob = atof(argv[1]);
    //generate_local_configuration_model(0 ,1,422);
    //generate_local_configuration_model(0.5,1000,4242);
    //generate_local_configuration_model(1.,1,424242);

    int seed = atoi(argv[1]);
    int N = atoi(argv[2]);
    int redes = atoi(argv[3]);
    int vacina = atoi(argv[4]);
    double prob = atof(argv[5]);
    double freq = atof(argv[6]);
    bool w = (bool)atoi(argv[7]);
    generate_infect(N,prob/100,seed, redes,freq/100,vacina,w);
}
