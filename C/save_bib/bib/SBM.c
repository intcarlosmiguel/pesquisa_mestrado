#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <igraph.h>
#include <unistd.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"

struct Graph{
    int **viz;
    int Nodes;
    double edges;
};

double** _distribution_faixas(){
    double** distribution = (double**) malloc(5*sizeof(double*));
    char filename[200];
    FILE *file;
    int N;
    for (int i = 0; i < 5; i++){
        sprintf(filename,"./dados/distribution_%d.txt",i);
        file = fopen(filename,"r");
        N = size_txt(filename);
        distribution[i] = (double*) calloc(N,sizeof(double));
        for (int j = 0; j < N; j++) if(fscanf(file, "%lf\n", &distribution[i][j]));
        fclose(file);
    }
    return distribution;
}

double** _get_array(int N,char* ptr){
    FILE* file;
    file = fopen(ptr,"r");
    double** degree = (double**) malloc(sizeof(double*)*N);
    if (access(ptr, F_OK) == -1) printf("Arquivo não acessado!\n");
    for(int i = 0; i < N; i++) {
        degree[i] = (double*) malloc(sizeof(double)*5);
        if(fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

void cria_ligacao(struct Graph *G,int* degree,int* faixas,double** prob,double** mean,double** std,int* shuff,int site,double p,bool weight,igraph_vector_int_t* edges,igraph_vector_t* pesos){
    int i = 0,j = 0;
    bool rep;
    int faixa1;
    int faixa2;
    double duracao = -1;
    for (int vizinho = 0; vizinho < G->Nodes - (site); vizinho++){
        i = shuff[vizinho];
        if(degree[site] == 0){
            //if(p == 0) aumenta_clustering(G,degree,faixas,mean,std,weight,site,p,edges,pesos);
            break;
        }
        rep = 0;
        faixa1 = faixas[site];
        faixa2 = faixas[i];
        if(degree[i] == 0) continue;
        for (j = 0; j < G->viz[site][0]; j++) if(i == G->viz[site][j+1]) {rep = true; break;}
        
        if(rep) continue;
        
        if(prob[faixa1][faixa2] < genrand64_real1()) continue;

        degree[site]--;
        degree[i]--;
        
        G->viz[site][0]++;
        G->viz[i][0]++;
        G->viz[site] = (int*) realloc(G->viz[site],(G->viz[site][0]+1)*sizeof(int));
        G->viz[i] = (int*) realloc(G->viz[i],(G->viz[i][0]+1)*sizeof(int));
        G->viz[site][G->viz[site][0]] = i;
        G->viz[i][G->viz[i][0]] = site;
        igraph_vector_int_push_back(edges, site);
        igraph_vector_int_push_back(edges, i);
        if(weight){
            
            while(duracao < 0) duracao = normalRand(mean[faixa1][faixa2],std[faixa1][faixa2]);
            igraph_vector_push_back(pesos, duracao/1440);
            G->edges += duracao/1440;
        }
        else G->edges += 1;
        duracao = -1;
    }
}

int* _criarVetor(int pontoInicial, int pontoFinal, int tamanhoPasso) {
    int tamanhoVetor = (pontoFinal - pontoInicial) / tamanhoPasso + 1;
    int* vetor = (int*)malloc(tamanhoVetor * sizeof(int));

    if (vetor == NULL) {
        printf("Erro: Não foi possível alocar memória para o vetor.\n");
        exit(1);
    }

    for (int i = 0; i < tamanhoVetor; i++) vetor[i] = pontoInicial + i * tamanhoPasso;

    return vetor;
}

igraph_t SBM_p(int N,double p,bool weight,int seed,double* avg){

    struct Graph G;
    int i,j;
    double r;
    init_genrand64(seed);
    char file_names[200] = "./dados/multi_probability.txt";
    double** prob = _get_array(5,file_names);
    double** distribution = _distribution_faixas();
    double* p_faixas = malloc(5*sizeof(double));

    strcpy(file_names, "./dados/media.txt");
    double** mean = _get_array(5,file_names);
    strcpy(file_names, "./dados/std.txt");
    double** std = _get_array(5,file_names);

    // ============================ Criando Faixas Etárias ============================
    int* faixas = calloc(N,sizeof(int));
    p_faixas[0] = 0.268392;
    p_faixas[1] = 0.152334;
    p_faixas[2] = 0.302137;
    p_faixas[3] = 0.206757;
    p_faixas[4] = 0.070380;

    for (i = 1; i < 5; i++) p_faixas[i] += p_faixas[i-1];
    for (i = 0; i < N; i++){
        j = 0;
        r = genrand64_real1();
        while(r > p_faixas[j]) j++;
        faixas[i] = j;
    }
    free(p_faixas);
    // ============================ Criando Distribuição de Graus ============================ //

    int* graus = calloc(N,sizeof(int));
    for (i = 0; i < N; i++){
        j = 0;
        while(j == 0)j = empiric_distribution(distribution[faixas[i]]);
        graus[i] = j;
    }
    
    // ============================ Criando Rede ============================ //

    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (i = 0; i < N; i++){ 
        G.viz[i] =(int*) malloc(1*sizeof(int));
        G.viz[i][0] = 0;
    }

    bubbleSort_by(faixas,graus,G.Nodes);

    G.edges = 0;
    igraph_vector_int_t edges;
    igraph_vector_t pesos;
    igraph_vector_int_init(&edges, 0);
    igraph_vector_init(&pesos, 0);
    for (i = 0; i < N; i++){
        int* shuff = _criarVetor(i+1, N, 1);
        randomize(shuff,N-i,seed+i);
        cria_ligacao(&G,graus,faixas,prob,mean,std,shuff,i,p,weight,&edges,&pesos);
        free(shuff);
    }
    //int lig = 0;
    //for ( i = 0; i < N; i++)lig += graus[i];
    //print_vetor(graus,G.Nodes,sizeof(int));
    //printf("%d %d\n",lig,G.edges);
    
    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_empty(&Grafo, G.Nodes, IGRAPH_UNDIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    igraph_vector_t faixa_;
    igraph_vector_init(&faixa_, G.Nodes);
    igraph_cattribute_VAN_setv(&Grafo,"faixa",&faixa_);
    if(weight) SETEANV(&Grafo, "duracao", &pesos);
    
    *avg = (double)G.edges*2/(G.Nodes);

    //printf("%f\n",avg_degree);

    /* FILE *file;
    file = fopen("./output/modelo/teste.txt","a");
    for (i = 0; i < G.Nodes; i++){
        int* M = calloc(5,sizeof(int));
        for ( j = 0; j < G.viz[i][0]; j++) M[faixas[G.viz[i][j+1]]]++;
        fprintf(file,"%d %d %d %d %d %d\n",faixas[i],M[0],M[1],M[2],M[3],M[4]);
        free(M);
    }
    fclose(file); */
    for (i = 0; i < 5; i++){
        free(prob[i]);
        free(distribution[i]);
        free(mean[i]);
        free(std[i]);
    }
    for (i = 0; i < G.Nodes; i++) free(G.viz[i]);
    //printf("%f\n",avg_degree);
    free(graus);
    free(prob);
    free(distribution);
    free(faixas);
    free(mean);
    free(std);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&pesos);
    igraph_vector_destroy(&faixa_);
    return Grafo;
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

