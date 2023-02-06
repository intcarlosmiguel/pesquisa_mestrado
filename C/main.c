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

int main(int argc,char *argv[ ]){
    int T = atoi(argv[1]);
    int model = atoi(argv[2]);
    switch (model){
    case 0:
        generate_configuration_model((double) 0,T);
        //generate_configuration_model((double) 0.5,T);
        //generate_configuration_model((double) 1.0,T);
        break;
    case 1:
        generate_SBM_p_model(T,model);
        break;
    case 2:
        generate_SBM_p_model(T,model);
        break;
    default:
        break;
    }
    
}
