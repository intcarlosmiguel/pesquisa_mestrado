#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"


struct Graph{
    int **viz;
    int Nodes;
    int edges;
};


int** get_lig(int N){
    FILE* file;
    file = fopen("./dados/ligacoes.txt","r");
    int **degree = (int**) malloc(sizeof(int*)*N);
    for(int i = 0; i < N; i++){
        degree[i] = (int*) malloc(sizeof(int)*N);
        if(fscanf(file,"%d %d %d %d %d\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

int* get_vetor(int N,char arquivo[]){
    FILE* file;
    file = fopen(arquivo,"r");
    int* degree = (int*) malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&degree[i]));
    fclose(file);
    return degree;
}

double** get_probability(int N){
    FILE* file;
    file = fopen("./dados/probability.txt","r");
    double **degree = (double**) malloc(sizeof(double*)*N);
    for(int i = 0; i < N; i++){
        degree[i] = (double*) malloc(sizeof(double)*N);
        if(fscanf(file,"%lf %lf %lf %lf %lf\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

struct Graph SBM_p(int N){

    struct Graph G;

    int* faixas = get_vetor(N,"./dados/faixas.txt");
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){
        G.viz[j] =(int*) malloc(1*sizeof(int));
        G.viz[j][0] = 0;
    }
    G.edges = 0;

    double **probability = get_probability(5);
    for (int i = 0; i < N; i++){
        int faixa1 = faixas[i];
        for (int j = i+1; j < N; j++){
            double rand = genrand64_real2();
            int faixa2 = faixas[j];
            if(rand <= probability[faixa1][faixa2]){
                G.viz[i][0]++;
                G.viz[j][0]++;
                G.viz[i] = (int*) realloc(G.viz[i],(G.viz[i][0]+1)*sizeof(int));
                G.viz[j] = (int*) realloc(G.viz[j],(G.viz[j][0]+1)*sizeof(int));
                G.viz[i][G.viz[i][0]] = j;
                G.viz[j][G.viz[j][0]] = i;
                G.edges += 1;
            }
        }
        
    }
    free(faixas);
    for (int i = 0; i < 5; i++) free(probability[i]);
    free(probability);
    
    return G;
}

struct SBM{
    int** mat;
    struct Graph G;
};

int maior_entre(int u, int v){
    if(u > v) return u;
    else return v;
}
int menor_entre(int u, int v){
    if(u > v) return v;
    else return u;
}

int* find_faixa(int** array,int* faixa,int N,int faixa1,int faixa2,int* first){
    int* resultados = (int*) malloc(2*sizeof(int));
    for (int i = *first; i < N; i++){
        if((faixa[array[i][0]] == faixa1) && (faixa[array[i][1]] == faixa2)){
            resultados[0] = array[i][0];
            resultados[1] = array[i][1];
            *first = i;
            return resultados;
        }
        if((faixa[array[i][0]] == faixa2) && (faixa[array[i][1]] == faixa1)){
            resultados[1] = array[i][0];
            resultados[0] = array[i][1];
            *first = i;
            return resultados;
        }
    }
    resultados[0] = -1;
    resultados[1] = -1;
    *first = -1;
    return resultados;
}

int check_clustering(int u, int v,int** viz){
    for (int i = 1; i <= viz[u][0]; i++){
        int first = viz[u][i];
        for (int j = 1; j <= viz[first][0]; j++) if(viz[first][j] == v) return 1;
    }
    for (int i = 1; i <= viz[v][0]; i++){
        int first = viz[v][i];
        //if(v == 419)printf("Erro: %d,%d\n",first,v);
        for (int j = 1; j <= viz[first][0]; j++) if(viz[first][j] == u) return 1;
    }
    //printf("Deu erro no segundo\n");
    return 0;
}

struct SBM increase_clustering(struct SBM SBM2,int* faixas){
    for (int i = 0; i < SBM2.G.Nodes; i++){
        for (int j = 1; j <= SBM2.G.viz[i][0]; j++){

            int first = SBM2.G.viz[i][j];

            for (int k = j+1; k <= SBM2.G.viz[i][0]; k++){
                int second = SBM2.G.viz[i][k];
                //printf("%d - %d,%d\n",SBM2.G.viz[419][0],first,second);
                int rep = check_existence(SBM2.mat,SBM2.G.edges,first,second);

                if(rep == 0){
                    int faixa1 = faixas[first];
                    int faixa2 = faixas[second];
                    int primeiro = 0;
                    //if((first == 251) && (second == 18)) printf("A: %d - %d,%d\n",i,first,second);
                    int* resultados = find_faixa(SBM2.mat,faixas,SBM2.G.edges,faixa1,faixa2,&primeiro);
                    while(check_clustering(resultados[0],resultados[1],SBM2.G.viz) == 1){
                        free(resultados);
                        primeiro++;
                        resultados = find_faixa(SBM2.mat,faixas,SBM2.G.edges,faixa1,faixa2,&primeiro);
                        if(primeiro == -1) break;
                    }
                    if(primeiro!=-1){
                        //if((first == 251) && (second == 18))printf("%d,%d\n",SBM2.mat[primeiro][0],SBM2.mat[primeiro][1]);
                        SBM2.mat[primeiro][0] = first;
                        SBM2.mat[primeiro][1] = second;
                        //if((first == 251) && (second == 18))printf("%d,%d\n",SBM2.mat[primeiro][0],SBM2.mat[primeiro][1]);
                        //if((j == 10) && (k == 18)) printf("%d,%d\n",SBM2.mat[primeiro][0],SBM2.mat[primeiro][1]);
                        SBM2.G.viz[first][0] += 1;
                        SBM2.G.viz[second][0] += 1;
                        SBM2.G.viz[first] = (int*) realloc(SBM2.G.viz[first],(SBM2.G.viz[first][0]+1)*sizeof(int));
                        SBM2.G.viz[second] = (int*) realloc(SBM2.G.viz[second],(SBM2.G.viz[second][0]+1)*sizeof(int));
                        SBM2.G.viz[first][SBM2.G.viz[first][0]] = second;
                        SBM2.G.viz[second][SBM2.G.viz[second][0]] = first;

                        //if((first == 251) && (second == 18)) print_vetor(SBM2.G.viz[resultados[0]],SBM2.G.viz[resultados[0]][0]+1);
                        SBM2.G.viz[resultados[0]] = ending(SBM2.G.viz[resultados[0]],SBM2.G.viz[resultados[0]][0]+1, resultados[1],1);
                        SBM2.G.viz[resultados[1]] = ending(SBM2.G.viz[resultados[1]],SBM2.G.viz[resultados[1]][0]+1, resultados[0],1);
                        
                        SBM2.G.viz[resultados[0]] = (int*) realloc(SBM2.G.viz[resultados[0]],(SBM2.G.viz[resultados[0]][0])*sizeof(int));
                        SBM2.G.viz[resultados[1]] = (int*) realloc(SBM2.G.viz[resultados[1]],(SBM2.G.viz[resultados[1]][0])*sizeof(int));
                        SBM2.G.viz[resultados[0]][0] -= 1;
                        SBM2.G.viz[resultados[1]][0] -= 1;
                        //if((resultados[0] == 419) || (resultados[1] == 419)) printf("Remove: %d,%d - %d,%d\n",resultados[0],resultados[1],SBM2.G.viz[resultados[0]][0],SBM2.G.viz[resultados[1]][0]);
                    }
                    free(resultados);
                }
            }
            
        }
        
    }
    return SBM2;   
}



struct Graph SBM_edges(int N,int seed){

    struct SBM SBM2;
    int *faixas = get_vetor(N,"./dados/faixas.txt");
    SBM2.G.Nodes = N;
    SBM2.G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){
        SBM2.G.viz[j] =(int*) malloc(1*sizeof(int));
        SBM2.G.viz[j][0] = 0;
    }
    SBM2.mat = (int **)malloc(0*sizeof(int*));
    SBM2.G.edges = 0;

    int **ligacoes = get_lig(5);
   
    int edges = 0;
    for(int i = 0;i < 5; i++) for(int j = 0; j < 5; j++)  edges += ligacoes[i][j];
    init_genrand64(seed);
    int repeat = 0;
    while(SBM2.G.edges != edges){

        int i = genrand64_int63() % SBM2.G.Nodes;
        int j = genrand64_int63() % SBM2.G.Nodes;
        int rep = check_existence(SBM2.mat,SBM2.G.edges,i,j);
        int faixa1 = faixas[i];
        int faixa2 = faixas[j];
        int maior = maior_entre(faixa1,faixa2);
        int menor = menor_entre(faixa1,faixa2);

        if((rep == 0) && (ligacoes[menor][maior]>= 1) && (i!=j)){
            SBM2.mat = (int**) realloc(SBM2.mat,(SBM2.G.edges+1)*sizeof(int*));
            SBM2.mat[SBM2.G.edges] = (int*) malloc(2* sizeof(int));
            SBM2.mat[SBM2.G.edges][0] = i;
            SBM2.mat[SBM2.G.edges][1] = j;

            SBM2.G.viz[i][0]++;
            SBM2.G.viz[j][0]++;
            SBM2.G.viz[i] = (int*) realloc(SBM2.G.viz[i],(SBM2.G.viz[i][0]+1)*sizeof(int));
            SBM2.G.viz[j] = (int*) realloc(SBM2.G.viz[j],(SBM2.G.viz[j][0]+1)*sizeof(int));
            SBM2.G.viz[i][SBM2.G.viz[i][0]] = j;
            SBM2.G.viz[j][SBM2.G.viz[j][0]] = i;
            SBM2.G.edges += 1;
            if(ligacoes[menor][maior] >=1)ligacoes[menor][maior] -= 1;
            repeat = 0;
            //if((faixa1 < faixa2) && (ligacoes[faixa1][faixa2]>=2)) ligacoes[faixa1][faixa2]-=2;
            //else ligacoes[faixa2][faixa1]-=2;
        }
        else{
            repeat += 1;
        }
        if(repeat == 500) break;
    }
    
    SBM2 = increase_clustering(SBM2,faixas);
    
    free(faixas);
    for (int i = 0; i < 5; i++) free(ligacoes[i]);
    for (int i = 0; i < SBM2.G.edges; i++) free(SBM2.mat[i]);
    free(ligacoes);
    free(SBM2.mat);
    return SBM2.G;
}

void generate_SBM_p_model(int T,int model){

    int N = size_txt();
    
    double** resultados = (double **)malloc(T*sizeof(double*));
    int sum = 1;
    #pragma omp parallel for
    for (int i = 0; i < T; i++){

        //if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        if(model == 1) G = SBM_p(N);
        if(model == 2) G = SBM_edges(N,i);
        resultados[i] = (double*) malloc(7*sizeof(double));

        result(G,resultados[i]);
        printf("%d\n",i+1);
        
        for(int j = 0; j < N; j++)free(G.viz[j]);
        printf("%d\n",sum);
        sum++;
        free(G.viz);
    }
    //if(model == 1)generate_resultados(resultados,T,"SBM_p");
    if(model == 2)generate_resultados(resultados,T,"SBM_edges_clustering");

}
