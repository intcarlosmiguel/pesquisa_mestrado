#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "mtwister.h"
#include "calc.h"
#include "rede.h"
#include <math.h>
#include <igraph.h>

struct Graph{
    int **viz;
    int Nodes;
    int edges;
};

int* get_degree(int N){
    FILE* file;
    file = fopen("./dados/graus.txt","r");
    int* degree = (int*) malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&degree[i]));
    fclose(file);
    return degree;
}

struct Graph conf_model_p(struct Graph G,int* degree,double p,int* shuff,int n,int site,igraph_vector_int_t* vector){    

    for (int j = 0; j < n; j++){

        if(degree[site] == 0) break;

        int vizinho = shuff[j];

        if(genrand64_real2()<=p){
            int rep = 0;
            for (int k = 0; k < G.viz[site][0]; k++) if((vizinho == G.viz[site][k+1]) || (vizinho == -G.viz[site][k+1])) {rep = 1; break;}
            if((rep == 1) || (degree[vizinho] == 0)) continue;

            G.edges += 1;

            degree[site]--;
            degree[vizinho]--;

            G.viz[site][0]++;
            G.viz[vizinho][0]++;
            G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
            G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
            G.viz[site][G.viz[site][0]] = vizinho;
            G.viz[vizinho][G.viz[vizinho][0]] = site;
            igraph_vector_int_push_back(vector, site);
            igraph_vector_int_push_back(vector, vizinho);
        }
        else{
            int rep = 0;
            for (int k = 0; k < G.viz[site][0]; k++) if(vizinho == -G.viz[site][k+1]) {rep = 1; break;}
            if(rep == 1) continue;
            G.viz[site][0]++;
            G.viz[vizinho][0]++;
            G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
            G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
            G.viz[site][G.viz[site][0]] = -vizinho;
            G.viz[vizinho][G.viz[vizinho][0]] = -site;
        }
    }
    
    return G;
}

struct Graph add_edge(struct Graph G,int* degree,int* shuff, int n,int site,double p,igraph_vector_int_t* vector){
    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(degree[site] == 0) break;
        int rep = 0;
        //if(site == 2) printf("Vizinho: %d\n",vizinho);
        for (int j = 0; j < G.viz[site][0]; j++) if(vizinho == G.viz[site][j+1]) {rep = 1; break;}

        if((rep == 1) || (degree[vizinho] == 0)) continue;
        
        G.edges += 1;
        degree[site]--;
        degree[vizinho]--;
        
        G.viz[site][0]++;
        G.viz[vizinho][0]++;
        G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
        G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
        G.viz[site][G.viz[site][0]] = vizinho;
        G.viz[vizinho][G.viz[vizinho][0]] = site;
        igraph_vector_int_push_back(vector, site);
        igraph_vector_int_push_back(vector, vizinho);
    }
    if(p >0){
        //printf("Entrou===========================================\n");
        int *vizinhos = (int*) malloc(0*sizeof(int));
        int n = 0;
        for (int i = 0; i < G.viz[site][0]; i++){
            int vizinho = G.viz[site][i+1];
            if(vizinho < 0) continue;
            if(degree[vizinho] !=0){
                vizinhos = (int*) realloc(vizinhos,(n+1)*sizeof(int));
                vizinhos[n] = vizinho;
                n++;
            }
        }
        if(n > 1){
            for (int i = 0; i < n-1; i++){
                int vizinho1 = vizinhos[i];
                int *shuff = (int*) malloc((n-(i+1))*sizeof(int));
                for (int j = i+1; j < n; j++) shuff[j - (i+1)] = vizinhos[j];
                shuff = randomize(shuff,n-(i+1),i);
                G = conf_model_p(G,degree,p,shuff,n-(i+1),vizinho1,vector);
                
                free(shuff);
            }
        }
        free(vizinhos);
    }
    
    return G;
}

igraph_t configuration_model(int N, double p,int seed){

    struct Graph G;
    char filename[200];
    
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){ 
        G.viz[j] =(int*) malloc(1*sizeof(int));
        G.viz[j][0] = 0;
    }
    G.edges = 0;
    
    int *shuff = (int*) malloc(N*sizeof(int));

    char file[200] = "./dados/degree.txt";
    sprintf(file,"./dados/degree.txt");
    int* degree = (int*) calloc(N,sizeof(int));
    load_file(file,degree,sizeof(degree[0]));
    degree = bubble_sort(degree,N);
    
    igraph_t Grafo;
    igraph_vector_int_t vector;
    igraph_vector_int_init(&vector, 0);
    igraph_empty(&Grafo, N, IGRAPH_UNDIRECTED);

    for (int i = 0; i < N; i++){
        init_genrand64(seed);

        for (int j = 0; j < N; j++) shuff[j] = j;
        swap(&shuff[i],&shuff[N-1]);
        shuff = randomize(shuff,N-1,seed+i);
        
        G = add_edge(G,degree,shuff,N-1,i,p,&vector);
        
    }
    
    igraph_add_edges(&Grafo, &vector, NULL);
    //create_network(G,0.0);
    
    free(shuff);
    igraph_vector_int_destroy(&vector);
    for(int j = 0; j < N; j++)free(G.viz[j]);
        
    free(G.viz);
    return Grafo;
}



void get_params(double* resultados,igraph_t* Grafo ){
    igraph_vector_int_t v;
    igraph_vector_t local_clustering;
    double clustering;
    double diametro;
    double path;
    int N = igraph_vcount(Grafo);
    
    igraph_vector_int_init(&v, N);
    igraph_vector_init(&local_clustering, 0);
    
    
    igraph_degree(Grafo, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);

    igraph_transitivity_local_undirected(Grafo,&local_clustering,igraph_vss_all(),IGRAPH_TRANSITIVITY_ZERO);
    double media = (double)igraph_vector_int_sum(&v)/N;
    int mediana = (N%2 == 0) ? VECTOR(v)[(int)N/2] : VECTOR(v)[(int)(N-1)/2];

    igraph_vector_int_add_constant(&v,-media);
    igraph_vector_int_mul(&v,&v);

    double std = (double)igraph_vector_int_sum(&v)/N;
    std = sqrt(std);

    igraph_transitivity_avglocal_undirected(Grafo,&clustering,IGRAPH_TRANSITIVITY_ZERO);

    igraph_diameter(Grafo, &diametro, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);
    igraph_average_path_length(Grafo,&path,NULL,IGRAPH_UNDIRECTED, 1);

    resultados[0] += media;
    resultados[1] += mediana;
    resultados[2] += std;
    resultados[3] += clustering;
    resultados[4] += path;
    resultados[5] += diametro;

    //printf("%f %d %f %f %f %f\n",media,mediana,std,clustering,path,diametro);

    igraph_vector_int_destroy(&v);
}   

void generate_configuration_model(double p, int T){
    char file[200] = "./dados/degree.txt";
    int N = size_txt(file);
    
    double* resultados = (double *)calloc(6,sizeof(double*));
    
    //#pragma omp parallel for
    
    for (int i = 0; i < T; i++){

        //if(T!=1) printf("\e[1;1H\e[2J");
        igraph_t G = configuration_model(N,p,i);
        //igraph_configuration_model(N,p);
        get_params(resultados,&G);
        
        //printf("Rodou: %d\n",i+1);
        igraph_destroy(&G);
    }
    for (int i = 0; i < 6; i++)printf("%f\n",resultados[i]/T);
    
    //generate_resultados(resultados,T,"CM");

}