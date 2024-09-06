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
    /*int seed = atoi(argv[1]);
    int N = atoi(argv[2]);
    int redes = atoi(argv[3]);
    double prob = atof(argv[4])/100;
    generate_local_configuration_model(N,prob ,redes,seed);*/


    int seed = atoi(argv[1]);
    double N = atof(argv[2]);
    uint16_t redes = atoi(argv[3]);
    uint8_t vacina = atoi(argv[4]);
    double prob = atof(argv[5])/100;
    bool weight = atoi(argv[6]);
    generate_infect(N,prob,seed, redes,vacina,weight);
    
}
