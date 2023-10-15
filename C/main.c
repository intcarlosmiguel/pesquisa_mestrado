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
    //generate_local_configuration_model(0.,100);
    //generate_local_configuration_model(0.5,100);
    //generate_local_configuration_model(1.,100);

    int seed = atoi(argv[1]);
    int N = atoi(argv[2]);
    int redes = atoi(argv[3]);
    int vacina = atoi(argv[4]);
    int prob = atoi(argv[5]);
    int freq = atof(argv[6]);
    generate_infect(N,(double)prob/100,seed, redes,(double) freq/100,vacina);
    
}
