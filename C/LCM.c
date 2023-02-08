
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

struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
struct LCM{
    int existir;
    int **n_existir;
    int **degree;
    int **mat;
    int *faixa;
    struct Graph G;
};
int** get_array(int N){
    FILE* file;
    file = fopen("./dados/degrees.txt","r");
    int** degree = (int**) malloc(sizeof(int*)*N);

    for(int i = 0; i < N; i++) {
        degree[i] = (int*) malloc(sizeof(int)*5);
        if(fscanf(file,"%d\t%d\t%d\t%d\t%d\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

int somatorio(int** array,int elemento,int N){
    int soma = 0;
    for (int i = 0; i < N; i++) soma += array[elemento][i];
    return soma;
}

struct LCM local_conf_model_p(struct LCM Z,int ego,double p){
    int N = 0;
    int *vizinhos = (int*) malloc(0*sizeof(int));

    for (int i = 0; i < Z.G.viz[ego][0]; i++){
        int vizinho = Z.G.viz[ego][i+1];
        if(somatorio(Z.degree,vizinho,N) !=0){
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
                if(somatorio(Z.degree,site,N) == 0) break;

                int vizinho = shuff[j];
                double rand = genrand64_real2();
                if(rand<=p){

                    int rep = check_existence(Z.mat,Z.G.edges,site,vizinho);
                    if((rep == 1) || (Z.degree[vizinho][Z.faixa[site]] == 0) || (Z.degree[site][Z.faixa[vizinho]] == 0) ) continue;
                    rep = check_existence(Z.n_existir,Z.existir,site,vizinho);
                    if(rep == 1) continue;
                    int faixa1 = Z.faixa[site];
                    int faixa2 = Z.faixa[vizinho];

                    Z.mat = (int**) realloc(Z.mat,(Z.G.edges+1)*sizeof(int*));
                    Z.mat[Z.G.edges] = (int*) malloc(2* sizeof(int));
                    Z.mat[Z.G.edges][0] = site;
                    Z.mat[Z.G.edges][1] = vizinho;
                    Z.G.edges += 1;

                    Z.degree[site][faixa2]--;
                    Z.degree[vizinho][faixa1]--;

                    Z.G.viz[site][0]++;
                    Z.G.viz[vizinho][0]++;
                    Z.G.viz[site] = (int*) realloc(Z.G.viz[site],(Z.G.viz[site][0]+1)*sizeof(int));
                    Z.G.viz[vizinho] = (int*) realloc(Z.G.viz[vizinho],(Z.G.viz[vizinho][0]+1)*sizeof(int));
                    Z.G.viz[site][Z.G.viz[site][0]] = vizinho;
                    Z.G.viz[vizinho][Z.G.viz[vizinho][0]] = site;
                }

                else{
                    int rep = check_existence(Z.n_existir,Z.existir,site,vizinho);
                    if(rep == 0){
                        Z.n_existir = (int**) realloc(Z.n_existir,(Z.existir+1)*sizeof(int*));
                        Z.n_existir[Z.existir] = (int*) malloc(2* sizeof(int));
                        Z.n_existir[Z.existir][0] = site;
                        Z.n_existir[Z.existir][1] = vizinho;
                        Z.existir += 1;
                    }
                }
            }
        }
        free(shuff);
    }
    free(vizinhos);
    return Z;
}

struct LCM local_add_edge(struct LCM Z,int* shuff, int n,int site,double p){

    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(somatorio(Z.degree,site,5) == 0){
            if(p>0) Z = local_conf_model_p(Z,site,p);
            break;
        }
        int rep = check_existence(Z.mat,Z.G.edges,site,vizinho);


        if((rep == 1) || (Z.degree[vizinho][Z.faixa[site]] == 0) || (Z.degree[site][Z.faixa[vizinho]] == 0) ) continue;
        int faixa1 = Z.faixa[site];
        int faixa2 = Z.faixa[vizinho];
        
        Z.mat = (int**) realloc(Z.mat,(Z.G.edges+1)*sizeof(int*));
        Z.mat[Z.G.edges] = (int*) malloc(2* sizeof(int));
        Z.mat[Z.G.edges][0] = site;
        Z.mat[Z.G.edges][1] = vizinho;
        Z.G.edges += 1;

        Z.degree[site][faixa2]--;
        Z.degree[vizinho][faixa1]--;
        
        Z.G.viz[site][0]++;
        Z.G.viz[vizinho][0]++;
        Z.G.viz[site] = (int*) realloc(Z.G.viz[site],(Z.G.viz[site][0]+1)*sizeof(int));
        Z.G.viz[vizinho] = (int*) realloc(Z.G.viz[vizinho],(Z.G.viz[vizinho][0]+1)*sizeof(int));
        Z.G.viz[site][Z.G.viz[site][0]] = vizinho;
        Z.G.viz[vizinho][Z.G.viz[vizinho][0]] = site;
    }
    return Z;
}

struct Graph local_configuration_model(int N, double p,int seed){

    struct LCM Z;
    Z.G.Nodes = N;
    Z.G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){ 
        Z.G.viz[j] =(int*) malloc(1*sizeof(int));
        Z.G.viz[j][0] = 0;
    }
    Z.mat = (int **)malloc(0*sizeof(int*));
    Z.G.edges = 0;
    Z.existir = 0;
    Z.n_existir = (int **)malloc(0* sizeof(int*));
    Z.degree = get_array(N);
    Z.faixa = get_faixas(N);

    init_genrand64(seed);
    
    int *shuff = (int*) malloc(N*sizeof(int));

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) shuff[j] = j;
        shuff = randomize(shuff,N,seed);
        shuff = ending(shuff,N, i,0);
        Z = local_add_edge(Z,shuff,N-1,i,p);
    }
    for (int i = 0; i < Z.G.Nodes; i++) free(Z.degree[i]);
    for (int i = 0; i < Z.G.edges; i++){
        free(Z.mat[i]);
        free(Z.n_existir[i]);
    }
    
    free(Z.n_existir);
    free(Z.degree);
    free(Z.mat);
    free(shuff);
    return Z.G;
}

void generate_local_configuration_model(double p, int T){

    int N = size_txt();
    
    double** resultados = (double **)malloc(T*sizeof(double*));
    #pragma omp parallel for
    for (int i = 0; i < T; i++){

        if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        G = local_configuration_model(N,p,i);
        resultados[i] = (double*) malloc(7*sizeof(double));
        result(G,resultados[i]);
        printf("%d\n",i+1);
        //if(T == 1) create_network(G,p);
        
        for(int j = 0; j < N; j++)free(G.viz[j]);
        
        free(G.viz);
    }
    generate_resultados(resultados,T,"LCM");

}