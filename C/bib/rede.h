#ifndef REDE_H
#define REDE_H

struct Graph;
void result(struct Graph G,double* resultados);
void create_network(struct Graph G,double p);
double shortest_length(int** viz,int N,int site,int* diametro);
double* list_clustering(int** viz,int N);

#endif