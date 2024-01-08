#ifndef SBM_H
#define SBM_H

struct Graph;
void generate_SBM_p_model(int T,int model,double p);
igraph_t SBM_p(int N,int seed,double p,bool weight,double* avg);

#endif