#ifndef LCM_H
#define LCM_H
struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
void generate_local_configuration_model(double p, int T,int seed);
igraph_t local_configuration_model(int N, double p,int seed,bool weight,double* avg_degree);


#endif