
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include <igraph.h>
#include "mtwister.h"
#include "calc.h"
#include "vacinas.h"
struct Graph{
    int **viz;
    int Nodes;
    double edges;
};

float lll = 0;

int** get_array(int N,char* ptr){
    FILE* file;
    file = fopen(ptr,"r");
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

double** distribution_faixas(){
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

void append_vizinhos(struct Graph *G,int site,int vizinho,int plus){
    G->viz[site][0]++;
    G->viz[vizinho][0]++;
    G->viz[site] = (int*) realloc(G->viz[site],(G->viz[site][0]+1)*sizeof(int));
    G->viz[vizinho] = (int*) realloc(G->viz[vizinho],(G->viz[vizinho][0]+1)*sizeof(int));
    G->viz[site][G->viz[site][0]] = plus*vizinho;
    G->viz[vizinho][G->viz[vizinho][0]] = plus*site;
}

bool check_repetidos(struct Graph *G,int i,int vizinho){
    for (int j = 0; j < G->viz[i][0]; j++) if(vizinho == abs(G->viz[i][j+1])) return true;
    return false;
}

double generate_conections(struct Graph *G,int** degree, igraph_vector_t* faixas,double p,igraph_t* Grafo){
    int ligacoes_total = 0;
    int i,j;
    int** matriz = (int**) malloc(5*sizeof(int*));
    for (i = 0; i < 5; i++) matriz[i] = (int*) calloc(5,sizeof(int));
    for (i = 0; i < G->Nodes; i++){
        //print_vetor(degree[i],5,sizeof(int));
        for (j = 0; j < 5; j++) matriz[(int) VECTOR(*faixas)[i]][j] += degree[i][j];
        ligacoes_total += somatorio(degree,i,5);
    }
    double agrupamento_medio;
    igraph_transitivity_avglocal_undirected(Grafo,&agrupamento_medio,IGRAPH_TRANSITIVITY_ZERO);
    printf("%f\n",agrupamento_medio);
    for (i = 0; i < 5; i++) free(matriz[i]);
    free(matriz);
   
    /*if((p == 0) || (p ==0.25) || (p ==0.5) || p == (1.0)){
        FILE *arquivo;
        char file[800];
        sprintf(file,"./output/modelo/resultados_%.2f.txt",p);
        arquivo = fopen(file,"a");
        igraph_vector_t clustering;
        igraph_vector_init(&clustering, G->Nodes);
        igraph_transitivity_local_undirected(Grafo,&clustering,igraph_vss_all(),IGRAPH_TRANSITIVITY_ZERO);
        for (i = 0; i < G->Nodes; i++){
            int* grau_ = calloc(5,sizeof(int));
            for (j = 0; j < G->viz[i][0]; j++){
                int vizinho = G->viz[i][j+1];
                if(vizinho <0) continue;
                grau_[(int) VECTOR(*faixas)[vizinho]]++;
            }
            fprintf(arquivo,"%d %d %d %d %d %d %f\n", (int) VECTOR(*faixas)[i],grau_[0],grau_[1],grau_[2],grau_[3],grau_[4],VECTOR(clustering)[i]);
            free(grau_);
        }
        igraph_vector_destroy(&clustering);
        fclose(arquivo);
    }*/
    return ligacoes_total/G->edges;
}

struct Graph weighted_add_edge(struct Graph *G,double p,int** degree,int** site_per_faixas,double** mean,double** std,int* n_faixas,int site,int faixa,bool weight,igraph_vector_int_t* edges,igraph_vector_t* faixas,igraph_vector_t* pesos,bool** liga,igraph_matrix_t *W){
    int i,vizinho,faixa1;
    double duracao = -1;
    faixa1 = (int) VECTOR(*faixas)[site];
    for ( i = 0; i < n_faixas[faixa]; i++){

        if(degree[site][faixa] == 0) break;

        vizinho = site_per_faixas[faixa][i];

        if(site == vizinho) continue;

        if(genrand64_real1() <= p){
            
            
            if(degree[vizinho][faixa1] == 0) continue;

            if(liga[site][vizinho]) continue;
            liga[site][vizinho] = true;
            liga[vizinho][site] = true;
            degree[site][faixa]--;
            degree[vizinho][faixa1]--;
            
            //append_vizinhos(G,site,vizinho,1);
            igraph_vector_int_push_back(edges, site);
            igraph_vector_int_push_back(edges, vizinho);

            G->viz[site][0]++;
            G->viz[vizinho][0]++;
            G->viz[site][G->viz[site][0]] = vizinho;
            G->viz[vizinho][G->viz[vizinho][0]] = site;
            if(weight){
                
                while(duracao < 0) duracao = normalRand(mean[faixa1][faixa],std[faixa1][faixa]);
                igraph_vector_push_back(pesos, duracao/1440);
                G->edges += duracao/1440;
                MATRIX(*W,site,vizinho) = duracao/1440;
                MATRIX(*W,vizinho,site) = duracao/1440;
            }
            else G->edges += 1;
            duracao = -1;
        }
        else{
            
            if(liga[site][vizinho]) continue;
            liga[site][vizinho] = true;
            liga[vizinho][site] = true;
        }
    }
    return *G;
}



igraph_t local_configuration_model(int N, double p,int seed,const bool weight,double *avg,igraph_vector_int_t* centralidade,uint8_t estrategy,bool calcula,igraph_matrix_t *W){

    struct Graph G;
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    int i = 0,j = 0,k = 0;

    if(weight) igraph_matrix_init(W, N, N); // Inicializa a matriz sem tamanho específico
    else igraph_matrix_init(W, 0, 0); // Inicializa a matriz sem tamanho específico
    
    G.edges = 0;

    double r;
    igraph_vector_t faixas;
    igraph_vector_init(&faixas, G.Nodes);

    int** degree = (int**) malloc(N*sizeof(int*));
    int* grau = (int*) calloc(N,sizeof(int));
    int* n_faixas = (int *) calloc(N, sizeof(int));
    int** site_per_faixas = (int**) malloc(5*sizeof(int*));
    bool** liga = (bool**) malloc(N*sizeof(bool*));

    char file_name[200] = "./dados/multi_probability.txt";
    double* p_faixas = (double *) calloc(5, sizeof(double));
    double** distribution = distribution_faixas();
    double** constant = double_get_array(file_name);

    strcpy(file_name, "./dados/media.txt");
    double** mean = double_get_array(file_name);
    strcpy(file_name, "./dados/std.txt");
    double** std = double_get_array(file_name);
    
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

    double* sintomatico = (double*) malloc(5*sizeof(double));
    double* hospitalizacao = (double*) malloc(5*sizeof(double));
    double* morte = (double*) malloc(5*sizeof(double));

    morte[0] = (double)0.07710809281267686;
    morte[1] = (double)0.09561476547321676;
    morte[2] = (double)0.13779231777947765;
    morte[3] = (double)0.30044562859568047;
    morte[4] = (double)0.530315172817809;

    sintomatico[0] = 29.1/100;
    sintomatico[1] = 37.4/100;
    sintomatico[2] = 41.68/100;
    sintomatico[3] = 39.4/100;
    sintomatico[4] = 31.3/100; 

    hospitalizacao[0] = (double)1.04/100;
    hospitalizacao[1] = (double)1.33/100;
    hospitalizacao[2] = (double)1.38/100;
    hospitalizacao[3] = (double)7.6/100;
    hospitalizacao[4] = (double)24/100;


    for ( j = 0; j < 5; j++){
        site_per_faixas[j] = (int*) malloc(0*sizeof(int));
        for (i = 1; i < 5; i++){
            constant[j][i] += constant[j][i-1];
        }
    }
    for (i = 1; i < 5; i++) p_faixas[i] += p_faixas[i-1];
    for (i = 0; i < N; i++){

        liga[i] = (bool*) calloc(N,sizeof(bool));

        j = 0;
        r = genrand64_real1();
        while(r > p_faixas[j]) j++;
        VECTOR(faixas)[i] = j;
        n_faixas[j]++;
        site_per_faixas[j] = realloc(site_per_faixas[j],n_faixas[j]*sizeof(int));
        site_per_faixas[j][n_faixas[j] - 1] = i;
        k = 0;

        degree[i] = (int*) calloc(5,sizeof(int));

        while(k == 0)k = empiric_distribution(distribution[ (int) VECTOR(faixas)[i]]);
        grau[i] = k;
    }
    mergeSort(grau,&faixas,0,N-1);
    for (i = 0; i < N; i++){
        generate_multinomial(grau[i], 5, constant[(int) VECTOR(faixas)[i]], degree[i]);
        G.viz[i] =(int*) malloc((grau[i]+1)*sizeof(int));
        G.viz[i][0] = 0;
    }
    for (i = 0; i < 5; i++) free(distribution[i]);
    

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    igraph_vector_t pesos;
    igraph_vector_init(&pesos, 0);
    int site;
    uint8_t faixa;
    int vizinho;
    uint8_t faixa1;
    for (site = 0; site < N; site++){
        //printf("%d\n",site);
        for (faixa = 0; faixa < 5; faixa++){
            seed++;
            if(degree[site][faixa] == 0) continue;
            site_per_faixas[faixa] = randomize(site_per_faixas[faixa], n_faixas[faixa],seed);
            G = weighted_add_edge(&G,1,degree,site_per_faixas,mean,std,n_faixas,site,faixa,weight,&edges,&faixas,&pesos,liga,W);
        }
        if(p > 0){
            int **vizinhos = (int**) malloc(5*sizeof(int*));
            int *q_vizinhos = (int*) calloc(5,sizeof(int));
            for (i = 0; i < 5; i++) vizinhos[i] = (int*) malloc(0*sizeof(int));
            int n = 0;
            for (i = 0; i < G.viz[site][0]; i++){
                
                vizinho = G.viz[site][i+1];
                if(vizinho < 0) continue;
                if(somatorio(degree,vizinho,5) !=0){
                    faixa1 = (int)VECTOR(faixas)[vizinho];
                    q_vizinhos[faixa1]++;
                    vizinhos[faixa1] = realloc(vizinhos[faixa1],q_vizinhos[faixa1]*sizeof(int));
                    vizinhos[faixa1][q_vizinhos[faixa1] - 1] = vizinho;
                    n++;
                }
            }
            if(n > 1){
                
                for (int viz = 0; viz < G.viz[site][0]; viz++){
                    vizinho = G.viz[site][viz+1];
                    
                    if(vizinho < 0) continue;
                    for (faixa = 0; faixa < 5; faixa++){
                        seed++;
                        if(degree[vizinho][faixa] == 0) continue;
                        vizinhos[faixa] = randomize(vizinhos[faixa], q_vizinhos[faixa],seed);
                        G = weighted_add_edge(&G,p,degree,vizinhos,mean,std,q_vizinhos,vizinho,faixa,weight,&edges,&faixas,&pesos,liga,W);
                    }
                }
            }
            for(i = 0; i < 5; i++) free(vizinhos[i]);
            free(vizinhos);
            free(q_vizinhos);
        }

    }
    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_empty(&Grafo, G.Nodes, IGRAPH_UNDIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    igraph_cattribute_VAN_setv(&Grafo,"faixa",&faixas);
    *avg = (double)G.edges*2/(G.Nodes);
    //lll += generate_conections(&G,degree,&faixas,p,&Grafo);
    if(calcula){
        if(estrategy <11) *centralidade = traditional_centralities(&Grafo,estrategy);
        else{
            if(estrategy<17) *centralidade = propose_centralities(&Grafo,estrategy,morte,hospitalizacao,sintomatico);
            else{
                if(estrategy < 21)*centralidade = weighted_traditional_centralities(&Grafo,estrategy,&pesos);
                else *centralidade = new_centralities(&Grafo,estrategy,&pesos,&edges,&faixas);
            }
        }
    }
    for (i = 0; i < G.Nodes; i++){
        free(G.viz[i]);
        free(degree[i]);
        free(liga[i]);
    }
    free(liga);
    free(G.viz);
    free(degree);
    free(p_faixas);
    free(grau);
    for(i = 0;i < 5;i++){
        free(constant[i]);
        free(site_per_faixas[i]);
        free(mean[i]);
        free(std[i]);
    }
    free(site_per_faixas);
    free(constant);
    free(mean);
    free(std);
    free(n_faixas);
    free(distribution);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&faixas);

    
    return Grafo;
}

double calcular_correlacao(igraph_vector_t* vetor1, igraph_vector_t* vetor2) {
    double soma1 = 0, soma2 = 0, soma1Q = 0, soma2Q = 0, prodSoma = 0;
    int n = igraph_vector_size(vetor1);

    for (int i = 0; i < n; i++) {
        soma1 += VECTOR(*vetor1)[i];
        soma2 += VECTOR(*vetor2)[i];
        soma1Q += pow(VECTOR(*vetor1)[i], 2);
        soma2Q += pow(VECTOR(*vetor2)[i], 2);
        prodSoma += VECTOR(*vetor1)[i] * VECTOR(*vetor2)[i];
    }

    double numerador = prodSoma - (soma1 * soma2 / n);
    double denominador = sqrt((soma1Q - pow(soma1, 2) / n) * (soma2Q - pow(soma2, 2) / n));

    if (denominador == 0) {
        return 0;  // Evita divisão por zero
    } else {
        return numerador / denominador;
    }
}


void calcula_propriedades(igraph_t *Grafo,double p,int N, double *resultados,double* std_vector) {

    double caminho_medio;
    double diametro;
    double agrupamento_medio;
    double agrupamento_total;
    double correlation;

    igraph_vector_int_t degrees;
    igraph_vector_int_init(&degrees, N);
    igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    igraph_vector_t clustering;
    igraph_vector_init(&clustering, N);
    igraph_transitivity_local_undirected(Grafo,&clustering,igraph_vss_all(),IGRAPH_TRANSITIVITY_ZERO);

    igraph_vector_t graus;
    igraph_vector_init(&graus, N);
    for (int i = 0; i < N; i++) VECTOR(graus)[i] = (double) VECTOR(degrees)[i];
    
    correlation = calcular_correlacao(&graus,&clustering);

    double k_medio = (double) igraph_vector_sum(&graus)/N;

    igraph_vector_add_constant(&graus,-k_medio);
    igraph_vector_mul(&graus,&graus);
    igraph_vector_scale(&graus,1/N);

    double std = pow(igraph_vector_sum(&graus),0.5);

    // Mediana dos graus
    igraph_vector_sort(&graus);
    double mediana = (N%2 == 0)? VECTOR(graus)[(N-1)/2] : (VECTOR(graus)[N/2] + VECTOR(graus)[N/2+1])*0.5;

    // Calcula o menor caminho médio da Rede
    igraph_average_path_length(Grafo, &caminho_medio, NULL, IGRAPH_UNDIRECTED, 1);

    // Calcula o diâmetro da rede
    igraph_diameter(Grafo, &diametro, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);

    // Calcula o agrupamento médio
    igraph_transitivity_avglocal_undirected(Grafo,&agrupamento_medio,IGRAPH_TRANSITIVITY_ZERO);

    // Calcula o agrupamento total
    igraph_transitivity_undirected(Grafo,&agrupamento_total,IGRAPH_TRANSITIVITY_NAN);
    resultados[0] += k_medio;
    resultados[1] += mediana;
    resultados[2] += std;
    resultados[3] += agrupamento_medio;
    resultados[4] += agrupamento_total;
    resultados[5] += correlation;
    resultados[6] += caminho_medio;
    resultados[7] += diametro;

    std_vector[0] += k_medio*k_medio;
    std_vector[1] += mediana*mediana;
    std_vector[2] += std*std;
    std_vector[3] += agrupamento_medio*agrupamento_medio;
    std_vector[4] += agrupamento_total*agrupamento_total;
    std_vector[5] += correlation*correlation;
    std_vector[6] += caminho_medio*caminho_medio;
    std_vector[7] += diametro*diametro;

    igraph_vector_destroy(&graus);
    igraph_vector_destroy(&clustering);
    igraph_vector_int_destroy(&degrees);

}


void generate_local_configuration_model(double p, int redes,int seed){

    int N = 10000;
    int i;
    double a;

    double* resultados = (double *)calloc(8,sizeof(double));
    double* std = (double *)calloc(8,sizeof(double));
    igraph_vector_int_t centralidade;
    igraph_matrix_t W;
    for (i = 0; i < redes; i++){

        //if(redes!=1) printf("\e[1;1H\e[2J");

        igraph_t G;
        G = local_configuration_model(N,p,seed+i,false,&a,&centralidade,0,false,&W);
        if(redes!=1) calcula_propriedades(&G,p,N,resultados,std);
        igraph_destroy(&G);
       
    }
    for (i = 0; i < 8; i++){
        resultados[i] /= redes;
        std[i] /= redes;
        std[i] = sqrt(std[i] - pow(resultados[i],2));
        std[i] /= sqrt(redes);
    }
    if(redes > 1){
        FILE *file;
        char filecheck[800];
        sprintf(filecheck,"./output/modelo/resultados.txt");
        file = fopen(filecheck,"a");
        fprintf(file,"%.2f %f %f %f %f %f %f %f %f %f\n",p,resultados[0],resultados[1],resultados[2],resultados[3],resultados[4],resultados[5],resultados[6],resultados[7],lll/redes);
        fclose(file);

        sprintf(filecheck,"./output/modelo/std.txt");
        file = fopen(filecheck,"a");
        fprintf(file,"%.2f %f %f %f %f %f %f %f %f\n",p,std[0],std[1],std[2],std[3],std[4],std[5],std[6],std[7]);
        fclose(file);

        /* FILE *arquivo;
        arquivo = fopen("./output/modelo/clustering.txt","a");
        fprintf(arquivo,"%f %f %f\n",p,resultados[3],std[3]);
        fclose(arquivo); */
    }
    free(resultados);
    free(std);
    //igraph_vector_int_destroy(&centralidade);
    printf("Terminou\n");
}