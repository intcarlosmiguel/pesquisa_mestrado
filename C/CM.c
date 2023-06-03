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

/*struct Graph conf_model_p(struct Graph G,int* degree,int ego,double p){
    int N = 0;
    int *vizinhos = (int*) malloc(0*sizeof(int));

    for (int i = 0; i < G.viz[ego][0]; i++){
        int vizinho = G.viz[ego][i+1];
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

                    int rep = check_existence(mat,G.edges,site,vizinho);
                    if((rep == 1) || (degree[vizinho] == 0)) continue;
                    rep = check_existence(n_existir,existir,site,vizinho);
                    if(rep == 1) continue;

                    mat = (int**) realloc(mat,(G.edges+1)*sizeof(int*));
                    mat[G.edges] = (int*) malloc(2* sizeof(int));
                    mat[G.edges][0] = site;
                    mat[G.edges][1] = vizinho;
                    G.edges += 1;

                    degree[site]--;
                    degree[vizinho]--;

                    G.viz[site][0]++;
                    G.viz[vizinho][0]++;
                    G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
                    G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
                    G.viz[site][G.viz[site][0]] = vizinho;
                    G.viz[vizinho][G.viz[vizinho][0]] = site;
                }

                else{
                    int rep = check_existence(n_existir,existir,site,vizinho);
                    if(rep == 0){
                        n_existir = (int**) realloc(n_existir,(existir+1)*sizeof(int*));
                        n_existir[existir] = (int*) malloc(2* sizeof(int));
                        n_existir[existir][0] = site;
                        n_existir[existir][1] = vizinho;
                        existir += 1;
                    }
                }
            }
        }
        free(shuff);
    }
    free(vizinhos);
    return G;
}*/

int* reverse(int* arr,int N){
    int* new = (int*) malloc(N*sizeof(int));
    for (int i = 0; i < N; i++) new[i] = arr[N-(i+1)];
    return new;
}

struct Graph conf_model_p(struct Graph G,int* degree,int ego,double p){
    int N = 0;
    int *vizinhos = (int*) malloc(0*sizeof(int));
    //int *degrees = (int*) malloc(0*sizeof(int));
    for (int i = 0; i < G.viz[ego][0]; i++){
        int vizinho = G.viz[ego][i+1];
        if(degree[vizinho] !=0){
            vizinhos = (int*) realloc(vizinhos,(N+1)*sizeof(int));
            //degrees = (int*) realloc(vizinhos,(N+1)*sizeof(int));
            vizinhos[N] = vizinho;
            //degrees[N] = degree[vizinho];
            N++;
        }
    }
    if(N > 1){
        //sortIntByRef(vizinhos,degrees,N,sizeof(degrees[0]));
        //vizinhos = reverse(vizinhos,N);
        
        //free(degrees);
        int *shuff = (int*) malloc(N*sizeof(int));

        for (int i = 0; i < N; i++){

            int site = vizinhos[i];

            for (int j = 0; j < N; j++) shuff[j] = vizinhos[j];
            swap(&shuff[i],&shuff[N-1]);
            //shuff = realloc(shuff,(N-1)*sizeof(int));
            //shuff = randomize(shuff,N-1,site+i);
            shuff = randomize(shuff,N-1,i+site);
            //shuff = ending(shuff,N,site,0);

            for (int j = 0; j < N-1; j++){

                if(degree[site] == 0) break;

                int vizinho = shuff[j];

                if(genrand64_real2()<=p){
                    int rep = 0;
                    for (int k = 0; k < G.viz[site][0]; k++) if(vizinho == G.viz[site][k+1]) {rep = 1; break;}
                    if((rep == 1) || (degree[vizinho] == 0)) continue;
                    if(rep == 1) continue;

                    G.edges += 1;

                    degree[site]--;
                    degree[vizinho]--;

                    G.viz[site][0]++;
                    G.viz[vizinho][0]++;
                    G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
                    G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
                    G.viz[site][G.viz[site][0]] = vizinho;
                    G.viz[vizinho][G.viz[vizinho][0]] = site;
                }
            }
        }
        free(shuff);
    }
    free(vizinhos);
    return G;
}

struct Graph add_edge(struct Graph G,int* degree,int* shuff, int n,int site,double p){
    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(degree[site] == 0){
            //printf("%d\n",site);
            if(p>0) G = conf_model_p(G,degree,site,p);
            break;
        }
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
    }
    return G;
}

struct Graph configuration_model(int N, double p,int seed){

    struct Graph G;
    int* degree = get_degree(N);
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){ 
        G.viz[j] =(int*) malloc(1*sizeof(int));
        G.viz[j][0] = 0;
    }
    G.edges = 0;
    //existir = 0;
    //n_existir = (int **)malloc(0* sizeof(int*));

    init_genrand64(seed);
    degree = bubble_sort(degree,N);
    int *shuff = (int*) malloc(N*sizeof(int));
    //print_vetor(degree,G.Nodes);
    //printf("Erro 1\n");
    for (int i = 0; i < N; i++){
        //shuff = realloc(shuff,N*sizeof(int));
        for (int j = 0; j < N; j++) shuff[j] = j;
        swap(&shuff[i],&shuff[N-1]);
        //shuff = realloc(shuff,(N-1)*sizeof(int));
        shuff = randomize(shuff,N-1,seed+i);
        G = add_edge(G,degree,shuff,N-1,i,p);
    }

    //create_network(G,0.0);
    //for(int i = 0; i < existir; i++)free(n_existir[i]);
    //free(n_existir);
    free(shuff);
    free(degree);
    //printf("Erro\n");
    return G;
}

void generate_configuration_model(double p, int T){
    char file[200] = "./dados/degree.txt";
    int N = size_txt(file);
    
    double** resultados = (double **)malloc(T*sizeof(double*));
    //#pragma omp parallel for
    for (int i = 0; i < T; i++){

        //if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        G = configuration_model(N,p,i);
        resultados[i] = (double*) malloc(7*sizeof(double));
        result(G,resultados[i]);
        
        for(int j = 0; j < N; j++)free(G.viz[j]);
        
        free(G.viz);
        printf("Rodou: %d\n",i+1);
    }
    generate_resultados(resultados,T,"CM");

}