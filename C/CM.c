#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "mtwister.h"
#include "calc.h"
#include "rede.h"
#include <math.h>

struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
struct CM{
    int existir;
    int **n_existir;
    int **mat;
    struct Graph G;
};


int* get_degree(int N){
    FILE* file;
    file = fopen("./dados/degree.txt","r");
    int* degree = (int*) malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&degree[i]));
    fclose(file);
    return degree;
}

struct CM conf_model_p(struct CM MC,int* degree,int ego,double p){
    int N = 0;
    int *vizinhos = (int*) malloc(0*sizeof(int));

    for (int i = 0; i < MC.G.viz[ego][0]; i++){
        int vizinho = MC.G.viz[ego][i+1];
        if(degree[vizinho] !=0){
            vizinhos = (int*) realloc(vizinhos,(N+1)*sizeof(int));
            vizinhos[N] = vizinho;
            N++;
        }
    }
    if(N > 1){

        int *shuff = (int*) malloc(N*sizeof(int));

        for (int i = 0; i < N; i++){
            int site = vizinhos[i];

            for (int j = 0; j < N; j++) shuff[j] = vizinhos[j];
            shuff = randomize(shuff,N,i+site);
            shuff = ending(shuff,N,site,0);
            for (int j = 0; j < N-1; j++){
                if(degree[site] == 0) break;
                int vizinho = shuff[j];
                double rand = genrand64_real2();
                if(rand<=p){

                    int rep = check_existence(MC.mat,MC.G.edges,site,vizinho);
                    if((rep == 1) || (degree[vizinho] == 0)) continue;
                    rep = check_existence(MC.n_existir,MC.existir,site,vizinho);
                    if(rep == 1) continue;

                    MC.mat = (int**) realloc(MC.mat,(MC.G.edges+1)*sizeof(int*));
                    MC.mat[MC.G.edges] = (int*) malloc(2* sizeof(int));
                    MC.mat[MC.G.edges][0] = site;
                    MC.mat[MC.G.edges][1] = vizinho;
                    MC.G.edges += 1;

                    degree[site]--;
                    degree[vizinho]--;

                    MC.G.viz[site][0]++;
                    MC.G.viz[vizinho][0]++;
                    MC.G.viz[site] = (int*) realloc(MC.G.viz[site],(MC.G.viz[site][0]+1)*sizeof(int));
                    MC.G.viz[vizinho] = (int*) realloc(MC.G.viz[vizinho],(MC.G.viz[vizinho][0]+1)*sizeof(int));
                    MC.G.viz[site][MC.G.viz[site][0]] = vizinho;
                    MC.G.viz[vizinho][MC.G.viz[vizinho][0]] = site;
                }

                else{
                    int rep = check_existence(MC.n_existir,MC.existir,site,vizinho);
                    if(rep == 0){
                        MC.n_existir = (int**) realloc(MC.n_existir,(MC.existir+1)*sizeof(int*));
                        MC.n_existir[MC.existir] = (int*) malloc(2* sizeof(int));
                        MC.n_existir[MC.existir][0] = site;
                        MC.n_existir[MC.existir][1] = vizinho;
                        MC.existir += 1;
                    }
                }
            }
        }
        free(shuff);
    }
    free(vizinhos);
    return MC;
}

struct CM add_edge(struct CM MC,int* degree,int* shuff, int n,int site,double p){

    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(degree[site] == 0){
            if(p>0) MC = conf_model_p(MC,degree,site,p);
            break;
        }
        int rep = check_existence(MC.mat,MC.G.edges,site,vizinho);


        if((rep == 1) || (degree[vizinho] == 0)) continue;
        
        MC.mat = (int**) realloc(MC.mat,(MC.G.edges+1)*sizeof(int*));
        MC.mat[MC.G.edges] = (int*) malloc(2* sizeof(int));
        MC.mat[MC.G.edges][0] = site;
        MC.mat[MC.G.edges][1] = vizinho;
        MC.G.edges += 1;
        degree[site]--;
        degree[vizinho]--;
        
        MC.G.viz[site][0]++;
        MC.G.viz[vizinho][0]++;
        MC.G.viz[site] = (int*) realloc(MC.G.viz[site],(MC.G.viz[site][0]+1)*sizeof(int));
        MC.G.viz[vizinho] = (int*) realloc(MC.G.viz[vizinho],(MC.G.viz[vizinho][0]+1)*sizeof(int));
        MC.G.viz[site][MC.G.viz[site][0]] = vizinho;
        MC.G.viz[vizinho][MC.G.viz[vizinho][0]] = site;
    }
    return MC;
}

struct Graph configuration_model(int N, double p,int seed){

    struct CM MC;

    int* degree = get_degree(N);
    MC.G.Nodes = N;
    MC.G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){
        MC.G.viz[j] =(int*) malloc(1*sizeof(int));
        MC.G.viz[j][0] = 0;
    }
    MC.mat = (int **)malloc(0*sizeof(int*));
    MC.G.edges = 0;
    MC.existir = 0;
    MC.n_existir = (int **)malloc(0* sizeof(int*));

    init_genrand64(seed);
    quicksort(degree,0,N-1);
    int *shuff = (int*) malloc(N*sizeof(int));

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) shuff[j] = j;
        shuff = randomize(shuff,N,seed);
        shuff = ending(shuff,N, i,0);
        MC = add_edge(MC,degree,shuff,N-1,i,p);
    }
    //if(seed == 0) create_network(MC.G,p);
    for(int i = 0; i < MC.existir; i++)free(MC.n_existir[i]);
    for(int i = 0; i < MC.G.edges; i++) free(MC.mat[i]);
    free(MC.n_existir);
    free(MC.mat);
    free(shuff);
    free(degree);

    return MC.G;
}

void generate_configuration_model(double p, int T){

    int N = size_txt();
    
    double** resultados = (double **)malloc(T*sizeof(double*));
    //#pragma omp parallel for
    for (int i = 0; i < T; i++){

        if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        G = configuration_model(N,p,i);
        resultados[i] = (double*) malloc(7*sizeof(double));
        result(G,resultados[i]);
        printf("%d\n",i+1);
        
        for(int j = 0; j < N; j++)free(G.viz[j]);
        
        free(G.viz);
    }
    generate_resultados(resultados,T,"CM");

}