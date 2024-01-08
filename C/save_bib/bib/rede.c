#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mtwister.h"
#include "calc.h"
#include <math.h>

struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
struct Graph G;


double* degree_list(int** viz,int N){
    double* degree = (double*) malloc(sizeof(double*)*N);
    for (int i = 0; i < N; i++) degree[i] = (double) viz[i][0];
    return degree;
}

double global_clustering(int** viz,int N){
    int c = 0;
    int n = 0;
    for (int i = 0; i < N; i++){
        c += combination(viz[i][0],2);
        if(combination(viz[i][0],2)!=0) for (int first = 1; first < viz[i][0]-1; first++) for (int second = first+1; second < viz[i][0]; second++) for(int vizinhos = 1; vizinhos<viz[second][0]; second++) if(viz[i][first] == viz[second][vizinhos]) {n++;break;}
    }
    return(double) n/c;
}

double avarage_clustering(int** viz,int N){
    int c = 0;
    double n = 0;
    double local = 0;
    for (int i = 0; i < N; i++){
        c = combination(viz[i][0],2);
        if(c!=0)
            for (int first = 1; first <=viz[i][0]; first++)
                for (int second = 1; second <=viz[i][0]; second++)
                    if(second!= first)
                        for(int vizinhos = 1; vizinhos<=viz[viz[i][second]][0]; vizinhos++) 
                            if(viz[i][first] == viz[viz[i][second]][vizinhos]) 
                                {n++;break;}
        if(c!=0)local += (double) n/(2*c);
        n = 0;
    }
    return(double) local/N;
}

double* list_clustering(int** viz,int N){
    double c = 0;
    double n = 0;
    double* clustering = (double*) malloc(N*sizeof(double));
    for (int i = 0; i < N; i++){
        c = combination(viz[i][0],2);
        if(c!=0) 
            for (int first = 1; first <=viz[i][0]; first++) 
                for (int second = 1; second <=viz[i][0]; second++) 
                    if(second!= first)
                        for(int vizinhos = 1; vizinhos<=viz[viz[i][second]][0]; vizinhos++) 
                            if(viz[i][first] == viz[viz[i][second]][vizinhos]) 
                                {n++;break;}
        clustering[i] = (c!=0)? n/(2*c): 0;
        n = 0;
    }
    return clustering;
}

int find_cluster(int** viz,int N,int site){
    int tam = 1;
    int infinity = N+2;
    int* lista = (int*) malloc(1*sizeof(int));
    int* distance = (int*) malloc(N*sizeof(int));
    int current = 0;

    for (int i = 0; i < N; i++) distance[i] = infinity;
    lista[0] = site;
    distance[site] = 0;
    while ((current<tam) && (tam != N)){

        int sitio = lista[current];
        int vizinhos = viz[sitio][0];
        for(int j = 1; j<=vizinhos;j++){
            int vizinho = viz[sitio][j];
            if((distance[vizinho]>distance[sitio]+1)){
                distance[vizinho] = distance[sitio]+1;
                tam++;
                lista = (int*) realloc(lista,tam*sizeof(int));
                lista[tam-1] = vizinho;
            }
        }
        current++;
    }
    free(distance);
    free(lista);
    return tam;
}

int* find_largest_cluster(struct Graph G,int* big){
    int maior = 0;
    int* resultado = (int*) malloc(sizeof(int)*G.Nodes);
    for(int i = 0; i < G.Nodes; i++){
        resultado[i] = find_cluster(G.viz,G.Nodes,i);
        if(resultado[i] > maior) maior = resultado[i];
    }
    int* resultado2 = malloc(sizeof(int)*maior);
    int s = 0;
    for(int i = 0; i < G.Nodes; i++){
        if(resultado[i] == maior) {
            resultado2[s] = i;
            s++;
        }
    }
    free(resultado);
    *big = maior;
    return resultado2;
}

double shortest_length(struct Graph G,int site,int* diametro){

    double l = 0;
    int tam = 1;
    int infinity = G.Nodes+2;
    int* lista = (int*) malloc(1*sizeof(int));
    int* distance = (int*) malloc(G.Nodes*sizeof(int));
    int current = 0;

    for (int i = 0; i < G.Nodes; i++) distance[i] = infinity;
    lista[0] = site;
    distance[site] = 0;
    *diametro = 0;
    if(G.viz[site][0] == 0) return 0;
    while ((current<tam) && (tam != G.Nodes)){

        int sitio = lista[current];
        int vizinhos = G.viz[sitio][0];
        for(int j = 1; j<=vizinhos;j++){
            int vizinho = G.viz[sitio][j];
            if((distance[vizinho]>distance[sitio]+1)){
                distance[vizinho] = distance[sitio]+1;
                tam++;
                lista = (int*) realloc(lista,tam*sizeof(int));
                lista[tam-1] = vizinho;
                if(distance[sitio]+1 > *diametro) *diametro = distance[sitio]+1;
            }
        }
        current++;
    }
    //printf("=======================\n");
    
    for (int i = 0; i < G.Nodes; i++)if(distance[i] != infinity) l += distance[i];
    free(distance);
    free(lista);
    //print_vetor(distance,tam);
    return (double) l/(tam-1);
}

double av_path_length(struct Graph G,int* diametro){
    double l = 0;
    //int big = 0;
    //int* sites = find_largest_cluster(G,&big);
    for(int i = 0; i < G.Nodes; i++) l += shortest_length(G,i,diametro);
    //free(sites);
    return l/G.Nodes;
}

void degree_distribution(struct Graph G,double* media1,double* median1,double* std1){

    double media = 0;
    double media2 = 0;
    int* degree_ = (int*) malloc(sizeof(int*)*G.Nodes);
    int median = (G.Nodes+1)/2 - 1;
    media = (double)2*(G.edges)/G.Nodes;

    for (int i = 0; i < G.Nodes; i++){
        degree_[i] = G.viz[i][0];
        media2 += pow(degree_[i]-media,2);
    }
    
    degree_ = bubble_sort(degree_,G.Nodes);

    *media1 = media;
    *std1 = pow(media2/(G.Nodes-1),0.5);
    *median1 = degree_[median];

    free(degree_);
}

void result(struct Graph G,double* resultados){
    int l;
    resultados[5] = av_path_length(G,&l);
    resultados[6] = (double) l;
    double *deg = degree_list(G.viz,G.Nodes);
    double *clustering = list_clustering(G.viz,G.Nodes);

    resultados[3] = 0;
    for (int i = 0; i < G.Nodes; i++) resultados[3] += clustering[i];
    resultados[3] /= G.Nodes;

    resultados[4] = correlation(deg,clustering,G.Nodes);
    degree_distribution(G,&resultados[0],&resultados[1],&resultados[2]);

    free(deg);
    free(clustering);
}
 void create_network(struct Graph G,double p){
    char arquivo[100];
    sprintf(arquivo, "./rede.txt");
    FILE* file;
    file = fopen(arquivo,"w");
    for (int i = 0; i < G.Nodes; i++) for (int j = 0; j < G.viz[i][0]; j++) fprintf(file,"%d %d\n",i,G.viz[i][j+1]);
    
    fclose(file);
} 