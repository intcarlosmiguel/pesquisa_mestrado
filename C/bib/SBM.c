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

int* get_vetor(char* arquivo){
    FILE* file;
    int N = size_txt(arquivo);
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
    char file[200] = "./dados/faixas.txt";
    int* faixas = get_vetor("./dados/faixas.txt");

    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    int* sitio = (int *)malloc(N*sizeof(int));
    int* degree = (int *)malloc(N*sizeof(int));
    for (int j = 0; j < N; j++) G.viz[j] =(int*) calloc(1,sizeof(int));
    G.edges = 0;

    double **probability = get_probability(5);

    for (int i = 0; i < N; i++){
        int faixa1 = faixas[i];
        for (int j = i+1; j < N; j++){
            int faixa2 = faixas[j];
            if(genrand64_real2() <= probability[faixa1][faixa2]){
                G.viz[i][0]++;
                G.viz[j][0]++;
                G.viz[i] = (int*) realloc(G.viz[i],(G.viz[i][0]+1)*sizeof(int));
                G.viz[j] = (int*) realloc(G.viz[j],(G.viz[j][0]+1)*sizeof(int));
                G.viz[i][G.viz[i][0]] = j;
                G.viz[j][G.viz[j][0]] = i;
                G.edges += 1;
            }
        }
        sitio[i] = i;
        
    }
    for (int i = 0; i < G.Nodes; i++) degree[i] = G.viz[i][0];
    FILE *ARQUIVO,*arquivo;
    ARQUIVO = fopen("./output/degrees.txt","w");
    arquivo = fopen("./output/faixas.txt","w");
    sortIntByRef(sitio,degree,G.Nodes,sizeof(degree[0]));
    for (int i = N-1; i >=0 ; i--){
        int site = sitio[i];
        int* faixa = (int *) calloc(5,sizeof(int));
        if(G.viz[site][0] != 0) {
            for (int j = 0; j < G.viz[site][0]; j++) faixa[faixas[G.viz[site][j+1]]]++;
            fprintf(ARQUIVO,"%d\t%d\t%d\t%d\t%d\n",faixa[0],faixa[1],faixa[2],faixa[3],faixa[4]);
            fprintf(arquivo,"%d\n",faixas[i]);
        }
        free(faixa);
    }
    fclose(ARQUIVO);
    fclose(arquivo);
    
    free(faixas);
    for (int i = 0; i < 5; i++) free(probability[i]);
    free(probability);
    
    return G;
}

struct SBM{
    int** mat;
    int** ligacoes;
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
        for (int j = 1; j <= viz[first][0]; j++) if(viz[first][j] == u) return 1;
    }
    return 0;
}


/* int length() {
  int *arr = (int*) malloc(5*sizeof(int));
  arr[0] = 5;
  int *p = arr;
  int count = 0;
  
  while (*p++) {
    count++;
  }
  
  printf("The length of the array is: %d\n", count);
  return 0;
} */

struct SBM increase_clustering(struct SBM SBM2,int* faixas,double p){
    for (int i = 0; i < SBM2.G.Nodes; i++){
        for (int j = 1; j <= SBM2.G.viz[i][0]; j++){

            int first = SBM2.G.viz[i][j];

            for (int k = j+1; k <= SBM2.G.viz[i][0]; k++){
                int second = SBM2.G.viz[i][k];
                int rep = check_existence(SBM2.mat,SBM2.G.edges,first,second);

                if(rep == 0){
                    int faixa1 = faixas[first];
                    int faixa2 = faixas[second];
                    int primeiro = 0;
                    int* resultados = find_faixa(SBM2.mat,faixas,SBM2.G.edges,faixa1,faixa2,&primeiro);
                    double p_ = genrand64_real1();
                    while((check_clustering(resultados[0],resultados[1],SBM2.G.viz) == 1) || (p_ > p)){
                        free(resultados);
                        primeiro++;
                        resultados = find_faixa(SBM2.mat,faixas,SBM2.G.edges,faixa1,faixa2,&primeiro);
                        if(primeiro == -1) break;
                        p_ = genrand64_real1();
                        
                    }
                    if(primeiro!=-1){
                        SBM2.mat[primeiro][0] = first;
                        SBM2.mat[primeiro][1] = second;
                        SBM2.G.viz[first][0] += 1;
                        SBM2.G.viz[second][0] += 1;
                        SBM2.G.viz[first] = (int*) realloc(SBM2.G.viz[first],(SBM2.G.viz[first][0]+1)*sizeof(int));
                        SBM2.G.viz[second] = (int*) realloc(SBM2.G.viz[second],(SBM2.G.viz[second][0]+1)*sizeof(int));
                        SBM2.G.viz[first][SBM2.G.viz[first][0]] = second;
                        SBM2.G.viz[second][SBM2.G.viz[second][0]] = first;

                        SBM2.G.viz[resultados[0]] = ending(SBM2.G.viz[resultados[0]],SBM2.G.viz[resultados[0]][0]+1, resultados[1],1);
                        SBM2.G.viz[resultados[1]] = ending(SBM2.G.viz[resultados[1]],SBM2.G.viz[resultados[1]][0]+1, resultados[0],1);
                        
                        SBM2.G.viz[resultados[0]] = (int*) realloc(SBM2.G.viz[resultados[0]],(SBM2.G.viz[resultados[0]][0])*sizeof(int));
                        SBM2.G.viz[resultados[1]] = (int*) realloc(SBM2.G.viz[resultados[1]],(SBM2.G.viz[resultados[1]][0])*sizeof(int));
                        SBM2.G.viz[resultados[0]][0] -= 1;
                        SBM2.G.viz[resultados[1]][0] -= 1;
                    }
                    free(resultados);
                }
            }
            
        }
        
    }
    return SBM2;   
}

struct Graph append_neighbors(struct Graph G,int site,int vizinho){
    G.viz[site][0]++;
    G.viz[vizinho][0]++;
    G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
    G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
    G.viz[site][G.viz[site][0]] = vizinho;
    G.viz[vizinho][G.viz[vizinho][0]] = site;
    return G;
}

struct Graph miguel_clustering(struct Graph G,int* faixas,int** ligacoes,int site, int vizinho,double p){
    if(G.viz[site][0] > 1){
        for (int i = 1; i <= G.viz[site][0]; i++){
            int outro = G.viz[site][i];
            int faixa1 = faixas[outro];
            int faixa2 = faixas[vizinho];
            int maior = maior_entre(faixa1,faixa2);
            int menor = menor_entre(faixa1,faixa2);
            int rep = 0;
            for(int k = 0; k < G.viz[site][0];k++) if(G.viz[site][k+1] == outro) {rep = 1;break;}
            if((rep == 0) && (ligacoes[menor][maior]>= 1) && (outro!=vizinho) && (genrand64_real1() <= p)){
                G = append_neighbors(G,outro,vizinho);
                G.edges += 1;
                if(ligacoes[menor][maior] >=1)ligacoes[menor][maior] -= 1;
            }

        }
        
    }
    return G;
}

/*struct Graph SBM_edges(int N,int seed,double p){

    char file[200] = "./dados/faixas.txt";
    int *faixas = get_vetor("./dados/faixas.txt");
    struct Graph G;

    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++) G.viz[j] =(int*) calloc(1,sizeof(int));
    G.edges = 0;

    int** ligacoes = get_lig(5);
   
    int edges = 0;
    for(int i = 0;i < 5; i++) for(int j = 0; j < 5; j++)  edges += ligacoes[i][j];
    init_genrand64(seed);
    int repeat = 0;
    while(G.edges != edges){

        int i = genrand64_real1()*G.Nodes;
        int j = genrand64_real1()*G.Nodes;
        int rep = 0;
        for(int k = 0; k < G.viz[i][0];k++) if(G.viz[i][k+1] == j) {rep = 1;break;}
        int faixa1 = faixas[i];
        int faixa2 = faixas[j];
        int maior = maior_entre(faixa1,faixa2);
        int menor = menor_entre(faixa1,faixa2);

        if((rep == 0) && (ligacoes[menor][maior]>= 1) && (i!=j)){

            G = append_neighbors(G,i,j);
            G.edges += 1;
            if(ligacoes[menor][maior] >=1)ligacoes[menor][maior] -= 1;
            repeat = 0;

            //if(p > 0) SBM2 = miguel_clustering(SBM2,faixas,i,j,p);
            //if(p > 0) SBM2 = miguel_clustering(SBM2,faixas,j,i,p);

        }
        else{
            repeat += 1;
        }
        if(repeat == 500) break;
    }
    //if(p > 0) SBM2 = increase_clustering(SBM2,faixas,p);
    if(seed == 0){
        FILE *file;
        file = fopen("./output/SBM/faixas_SBM.txt","w");
        int* Nodes_faixa = (int*) malloc(5*sizeof(int));
        for(int j = 0; j < 5; j++) Nodes_faixa[j] = 0;
        for(int j = 0; j < G.Nodes; j++){
            for (int k = 1; k <= G.viz[j][0];k++) Nodes_faixa[faixas[G.viz[j][k]]]++;
            
            fprintf(file,"%d\t%d\t%d\t%d\t%d\n",Nodes_faixa[0],Nodes_faixa[1],Nodes_faixa[2],Nodes_faixa[3],Nodes_faixa[4]);
            for(int k = 0; k < 5; k++) Nodes_faixa[k] = 0;
        }
        fclose(file);
    }
    free(faixas);
    for (int i = 0; i < 5; i++) free(ligacoes[i]);
    free(ligacoes);
    return G;
}*/

void generate_SBM_p_model(int T,int model,double p){
    char file[200] = "./dados/degree.txt";
    int N = size_txt(file);
    double** resultados = (double **)malloc(T*sizeof(double*));
    //#pragma omp parallel for
    for (int i = 0; i < T; i++){

        //if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        if(model == 1) G = SBM_p(N);
        //if(model == 2) G = SBM_edges(N,i,p);
        resultados[i] = (double*) malloc(7*sizeof(double));

        result(G,resultados[i]);
        for(int j = 0; j < N; j++)free(G.viz[j]);
        free(G.viz);
    }
    generate_resultados(resultados,T,"SBM_p");
    //if(model == 2)generate_resultados(resultados,T,"SBM_edges_miguel_clustering");

}