#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mtwister.h"
#include "calc.h"
#include <math.h>

struct Graph{
    int **mat;
    int **viz;
    int **n_existir;
    int Nodes;
    int edges;
    int existir;
};