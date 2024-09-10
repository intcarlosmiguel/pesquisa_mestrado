#ifndef LCM_H
#define LCM_H
struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
void generate_local_configuration_model(int N,double p, int redes,int seed);
igraph_t local_configuration_model(int N, double p,int seed,const bool weight,double *avg,igraph_vector_int_t* centralidade,uint8_t estrategy,bool calcula,igraph_matrix_t *W,double* perca);


#endif