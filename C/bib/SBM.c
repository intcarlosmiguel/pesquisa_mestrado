#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <igraph.h>

#include "mtwister.h"
#include "calc.h"


struct Graph{
    int **viz;
    int Nodes;
    double edges;
};

void append_neighbors(struct Graph *G,int site,int vizinho){
    G->viz[site][0]++;
    G->viz[vizinho][0]++;
    G->viz[site] = (int*) realloc(G->viz[site],(G->viz[site][0]+1)*sizeof(int));
    G->viz[vizinho] = (int*) realloc(G->viz[vizinho],(G->viz[vizinho][0]+1)*sizeof(int));
    G->viz[site][G->viz[site][0]] = vizinho;
    G->viz[vizinho][G->viz[vizinho][0]] = site;
}

/* void local_add_edge_p(struct Graph *G,int* degree,int* faixas,int site,double p,igraph_vector_int_t* edges){
    int* vizinhos = (int*) malloc(0*sizeof(int));
    int n = 0;
    int vizinho;
    for (int i = 0; i < G->viz[site][0]; i++){
        
        vizinho = G->viz[site][i+1];
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
            //ocal_conf_model_p(G,faixa,degree,p,shuff,n-(i+1),vizinho1,edges);
            
            free(shuff);
        }
    }
    free(vizinhos);
    return G;
} */

void SBM_add_edge(struct Graph* G,int* degree,int* faixa,int* shuff,double **constant,double** mean,double** std, int n,int site,double p,bool weight, igraph_vector_int_t* edges,igraph_vector_t* pesos){
    bool rep;
    int i,vizinho;
    double duracao = -1;
    for (i = 0; i < n; i++){
        vizinho = shuff[i];
        if(degree[site] == 0) break;
        rep = false;
        for (int j = 0; j < G->viz[site][0]; j++) if(vizinho == G->viz[site][j+1]) {rep = true; break;}

        int faixa1 = faixa[site];
        int faixa2 = faixa[vizinho];

        if((rep) || (degree[vizinho] == 0) ) continue;
        //if(0.2< genrand64_real1()) continue;
        //if(constant[faixa1][faixa2] < genrand64_real1()) continue;
        G->edges += 1;

        degree[site]--;
        degree[vizinho]--;
        
        append_neighbors(G,site,vizinho);

        igraph_vector_int_push_back(edges, site);
        igraph_vector_int_push_back(edges, vizinho);
        if(weight){
            
            while(duracao < 0) duracao = normalRand(mean[faixa1][faixa2],std[faixa1][faixa2]);
            igraph_vector_push_back(pesos, duracao/1440);
            G->edges += duracao/1440;
        }
        else G->edges += 1;
        duracao = -1;
    }
    //if((p > 0) && (somatorio(degree,site,5) == 0)) G = SBM_add_edge_p(G,degree,faixa,site,p,edges);
    //return G;
}

double** load_distribution_faixas(){
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

igraph_t SBM_p(int N,int seed,double p,bool weight,double* avg){

    init_genrand64(seed);
    struct Graph G;
    int* degree = (int *) calloc(N, sizeof(int));
    int* faixas = (int *) calloc(N, sizeof(int));
    int** site_per_faixa = (int **) malloc(5*sizeof(int*));
    int* n_faixas = (int *) calloc(5, sizeof(int));
    double* k_faixas = (double *) calloc(5, sizeof(double));
    char file_name[200] = "./dados/multi_probability_density.txt";
    double* p_faixas = (double *) calloc(5, sizeof(double));
    double** distribution = load_distribution_faixas();
    double** constant = double_get_array(file_name);

    strcpy(file_name, "./dados/media.txt");
    double** mean = double_get_array(file_name);
    strcpy(file_name, "./dados/std.txt");
    double** std = double_get_array(file_name);

    int i,j,k;
    double r;
    
    // POLYMOD//

    /*p_faixas[0] = 0.34393182;
    p_faixas[1] =  0.14065873;
    p_faixas[2] = 0.31957894;
    p_faixas[3] = 0.15781574;
    p_faixas[4] = 0.03801476;*/

    // BRASILEIRA //
    p_faixas[0] = 0.268392;
    p_faixas[1] = 0.152334;
    p_faixas[2] = 0.302137;
    p_faixas[3] = 0.206757;
    p_faixas[4] = 0.070380;

     // TESTE //
    p_faixas[0] = 0.05;
    p_faixas[1] = 0.1;
    p_faixas[2] = 0.3;
    p_faixas[3] = 0.2;
    p_faixas[4] = 0.35;
    
    // =================================== Gera faixas ================================== //
    for ( i = 0; i < 5; i++) site_per_faixa[i] = malloc(0*sizeof(int));
    
    for (j = 0; j < 5; j++)for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
    for (i = 1; i < 5; i++)p_faixas[i] += p_faixas[i-1];
    for (i = 0; i < N; i++){
        j = 0;
        r = genrand64_real1();
        while(r > p_faixas[j]) j++;
        faixas[i] = j;
        n_faixas[j]++;
        site_per_faixa[j] = (int*) realloc(site_per_faixa[j],n_faixas[j]*sizeof(int));
        site_per_faixa[j][n_faixas[j] - 1] = i;
    }

    // =================================== GERA GRAUS ================================== //
    for (i = 0; i < N; i++) while(degree[i] == 0) degree[i] = 100;//empiric_distribution(distribution[faixas[i]]);
    for (i = 0; i < 5; i++)free(distribution[i]);
    
    // =================================== GERA REDE ================================== //

    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (i = 0; i < N; i++){ 
        G.viz[i] = (int*) malloc(1*sizeof(int));
        G.viz[i][0] = 0;
    }
    G.edges = 0;
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    igraph_vector_t pesos;
    igraph_vector_init(&pesos, 0);
    bubbleSort_by(faixas, degree, N);
    //for ( j = 0; j < 5; j++) for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
    
    for (i = 0; i < N; i++){
        if(degree[i] == 0) continue;
        bool* check =  (bool *) calloc(5,sizeof(int));
        while(degree[i] != 0){
            seed++;
            r = genrand64_real1();
            j = 0;
            
            while(constant[faixas[i]][j] <= r)j++;
            int soma = 0;
            for(int l =0;l<5;l++) soma +=check[l];
            if(soma == 5){
                break;
            }
            if(check[j]) continue;
            site_per_faixa[j] = randomize(site_per_faixa[j],n_faixas[j],seed+i);
            for ( k = 0; k < n_faixas[j]; k++){
                int vizinho = site_per_faixa[j][k];
                
                double duracao = -1;
                bool rep = false;
                for (int l = 0; l < G.viz[i][0]; l++) if(vizinho == G.viz[i][l+1]) {rep = true; break;}
                if((rep) || (degree[vizinho] == 0) || (vizinho == i) ) continue;
                int faixa1 = faixas[i];
                int faixa2 = faixas[vizinho];

                degree[i]--;
                degree[vizinho]--;
                
                append_neighbors(&G,i,vizinho);

                igraph_vector_int_push_back(&edges, i);
                igraph_vector_int_push_back(&edges, vizinho);
                if(weight){
                    
                    while(duracao < 0) duracao = normalRand(mean[faixa1][faixa2],std[faixa1][faixa2]);
                    igraph_vector_push_back(&pesos, duracao/1440);
                    G.edges += duracao/1440;
                }
                else G.edges += 1;
                duracao = -1;
                break;
            }
            if(k  == n_faixas[j]) check[j] = true;
        }
        
        free(check);

        /* int* shuff = arange(i+1, N-1, 1);
        
        shuff = randomize(shuff,N-i-1,seed+i);
        //print_vetor(grau,5,sizeof(int));
        //printf("%d\n",degree[i]);

        //SBM_add_edge(&G,degree,faixas,shuff,constant,mean,std,N-i-1,i,p,weight,&edges,&pesos,grau);

        
        free(shuff); */
    }
    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_empty(&Grafo, G.Nodes, IGRAPH_UNDIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    igraph_vector_t faixa_;
    igraph_vector_init(&faixa_, G.Nodes);
    
    if(weight) SETEANV(&Grafo, "duracao", &pesos);
    
    *avg = (double)G.edges*2/(G.Nodes);

    k = 0;
    FILE *arquivo;
    arquivo = fopen("./testdist.txt","a");
    for (i = 0; i < N; i++){
        int* faixa = calloc(5,sizeof(int));
        k += degree[i];
        k_faixas[faixas[i]] += G.viz[i][0];
        for (j = 0; j < G.viz[i][0]; j++){
            faixa[faixas[G.viz[i][j+1]]]++;
        }
        fprintf(arquivo,"%d %d %d %d %d %d\n",faixas[i],faixa[0],faixa[1],faixa[2],faixa[3],faixa[4]);
        VECTOR(faixa_)[i] = faixas[i];
        free(G.viz[i]);
        free(faixa);
    }
    igraph_cattribute_VAN_setv(&Grafo,"faixa",&faixa_);
    fclose(arquivo);
    for (i = 0; i < 5; i++){
        free(constant[i]);
        free(mean[i]);
        free(std[i]);
    }
    /* arquivo = fopen("./resultados2.txt","a");
    printf("%f\n",G.edges/G.Nodes*2);
    fprintf(arquivo,"%f %f %f %f %f\n",k_faixas[0]/n_faixas[0],k_faixas[1]/n_faixas[1],k_faixas[2]/n_faixas[2],k_faixas[3]/n_faixas[3],k_faixas[4]/n_faixas[4]);
    fclose(arquivo); */
    //igraph_cattribute_VAN_setv(&Grafo,"faixa",&faixa_);
    free(n_faixas);
    free(faixas);
    free(degree);
    free(p_faixas);
    free(distribution);
    free(constant);
    free(G.viz);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&pesos);
    igraph_vector_destroy(&faixa_);
    return Grafo;
}




/* void generate_SBM_p_model(int T,int model,double p){
    char file[200] = "./dados/degree.txt";
    int N = size_txt(file);
    double** resultados = (double **)malloc(T*sizeof(double*));
    for (int i = 0; i < T; i++){

        //if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        if(model == 1) G = SBM_p(N,42,0);
        //if(model == 2) G = SBM_edges(N,i,p);
        resultados[i] = (double*) malloc(7*sizeof(double));

        result(G,resultados[i]);
        for(int j = 0; j < N; j++)free(G.viz[j]);
        free(G.viz);
    }
    generate_resultados(resultados,T,"SBM_p");
    //if(model == 2)generate_resultados(resultados,T,"SBM_edges_miguel_clustering");

} */
