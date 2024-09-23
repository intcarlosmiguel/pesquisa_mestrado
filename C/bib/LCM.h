#ifndef LCM_H
#define LCM_H
struct Graph{
    int **viz;
    double **W;
    int Nodes;
    double edges;
    double avg;
    int *faixas;
};
void generate_local_configuration_model(int N,double p, int redes,int seed);
void local_configuration_model(struct Graph* Grafo,int N, double p,int seed,const bool weight,igraph_vector_int_t* centralidade,uint8_t estrategy,bool calcula,double* perca);


#endif