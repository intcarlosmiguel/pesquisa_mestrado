
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"
#include <igraph.h>
struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
struct Retorno{
    int* faixa;
    int** degree;
};

double distance(double* point1,double* point2,int N){
    double distance2 = 0;
    for (int i = 0; i < N; i++) distance2 += pow(point1[i] - point2[i],2);
    distance2 = sqrt(distance2);
    return distance2;
}

double** double_get_array(int N,char* ptr){
    FILE* file;
    file = fopen(ptr,"r");
    double** degree = (double**) malloc(sizeof(double*)*N);

    for(int i = 0; i < N; i++) {
        degree[i] = (double*) malloc(sizeof(double)*5);
        if(fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

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

int vetor_soma(int* array, int L){
    int soma = 0;
    for (int i = 0; i < L; i++) soma += array[i];
    return soma;
}

int somatorio(int** array,int elemento,int N){
    int soma = 0;
    for (int i = 0; i < N; i++) soma += array[elemento][i];
    return soma;
}

int** M_bubbleSort(int **matriz, int linhas, int colunas) {
    int i, j;
    for (i = 0; i < linhas - 1; i++) {
        for (j = 0; j < linhas - i - 1; j++) {
            int somaAtual = 0, somaProxima = 0;
            int *linhaAtual = matriz[j];
            int *linhaProxima = matriz[j + 1];

            for (int k = 0; k < colunas; k++) {
                somaAtual += linhaAtual[k];
                somaProxima += linhaProxima[k];
            }

            if (somaAtual < somaProxima) {
                int *temp = matriz[j];
                matriz[j] = matriz[j + 1];
                matriz[j + 1] = temp;
            }
        }
    }
    return matriz;
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

void calc_metrics(struct Graph G,int* faixa){

    int i;
    int* n_faixa = calloc(5,sizeof(int));
    double* media_grau = calloc(5,sizeof(double));
    double k_medio = 0;

    for (i = 0; i < G.Nodes; i++){
        k_medio += G.viz[i][0];
        media_grau[faixa[i]] += G.viz[i][0];
        n_faixa[faixa[i]]++;
    }
    for ( i = 0; i < 5; i++) media_grau[i] /= n_faixa[i];
    printf("%f\n",k_medio/G.Nodes);
    print_vetor(media_grau,5,sizeof(double));
    free(n_faixa);
    free(media_grau);
}

void generate_file2(struct Graph G,int* faixa){

    int i,j;
    FILE* arquivo;
    FILE* arquivo2;
    char test[200] = "./test.txt";
    char test2[200] = "./test2.txt";

    arquivo = fopen(test,"a");
    arquivo2 = fopen(test2,"a");

    int** final = malloc(G.Nodes*sizeof(int*));
    for (i = 0; i < G.Nodes; i++){
        final[i] = calloc(5,sizeof(int));
        for (j = 0; j < G.viz[i][0]; j++) final[i][faixa[G.viz[i][j+1]]]++;
        for (j = 0; j < 5; j++){
            if(j!= 4) fprintf(arquivo,"%d\t",final[i][j]);
            else fprintf(arquivo,"%d\n",final[i][j]);
        }
        fprintf(arquivo2,"%d\n",faixa[i]);

    }
    fclose(arquivo);
    fclose(arquivo2);
}



struct Graph local_conf_model_p(struct Graph G,int* faixa,int** degree,double p,int* shuff,int n,int site,igraph_vector_int_t* edges){    

    for (int j = 0; j < n; j++){

        if(somatorio(degree,site,5) == 0) break;

        int vizinho = shuff[j];

        if(genrand64_real2()<=p){
            
            int rep = 0;
            for (int j = 0; j < G.viz[site][0]; j++) if(vizinho == abs(G.viz[site][j+1])) {rep = 1; break;}
            
            int faixa1 = faixa[site];
            int faixa2 = faixa[vizinho];
            if((rep == 1) || (degree[vizinho][faixa1] == 0) || (degree[site][faixa2] == 0) ) continue;
            G.edges += 1;

            degree[site][faixa2]--;
            degree[vizinho][faixa1]--;
            
            G.viz[site][0]++;
            G.viz[vizinho][0]++;
            G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
            G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
            G.viz[site][G.viz[site][0]] = vizinho;
            G.viz[vizinho][G.viz[vizinho][0]] = site;

            igraph_vector_int_push_back(edges, site);
            igraph_vector_int_push_back(edges, vizinho);
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

struct Graph local_add_edge_p(struct Graph G,int** degree,int* faixa,int site,double p,igraph_vector_int_t* edges){
    int *vizinhos = (int*) malloc(0*sizeof(int));
    int n = 0;
    for (int i = 0; i < G.viz[site][0]; i++){
        
        int vizinho = G.viz[site][i+1];
        if(vizinho < 0) continue;
        if(somatorio(degree,vizinho,5) !=0){
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
            G = local_conf_model_p(G,faixa,degree,p,shuff,n-(i+1),vizinho1,edges);
            
            free(shuff);
        }
    }
    free(vizinhos);
    return G;
}

struct Graph local_add_edge(struct Graph G,int** degree,int* faixa,int* shuff, int n,int site,double p,igraph_vector_int_t* edges){
    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(somatorio(degree,site,5) == 0) break;
        int rep = 0;
        for (int j = 0; j < G.viz[site][0]; j++) if(vizinho == G.viz[site][j+1]) {rep = 1; break;}
        

        int faixa1 = faixa[site];
        int faixa2 = faixa[vizinho];
        if((rep == 1) || (degree[vizinho][faixa1] == 0) || (degree[site][faixa2] == 0) ) continue;
        G.edges += 1;

        degree[site][faixa2]--;
        degree[vizinho][faixa1]--;
        
        G.viz[site][0]++;
        G.viz[vizinho][0]++;
        G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
        G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
        G.viz[site][G.viz[site][0]] = vizinho;
        G.viz[vizinho][G.viz[vizinho][0]] = site;

        igraph_vector_int_push_back(edges, site);
        igraph_vector_int_push_back(edges, vizinho);
    }
    if((p > 0) && (somatorio(degree,site,5) == 0)){
        G = local_add_edge_p(G,degree,faixa,site,p,edges);
    }
    return G;
}

struct Retorno generate_degree(int N,int seed){
    struct Retorno Return;

    init_genrand64(seed);

    char file_name[200] = "./dados/multi_probability.txt";
    double** constant = double_get_array(5,file_name);
    double** distribution = distribution_faixas();

    int* n_faixa = calloc(5,sizeof(int));

    Return.degree = (int**) malloc(N*sizeof(int*));
    int* grau = (int*) malloc(N*sizeof(int));
    double* p_faixas = malloc(5*sizeof(double));
    int i = 0,j = 0;
    double r;

    p_faixas[0] = 0.34393182;
    p_faixas[1] =  0.14065873;
    p_faixas[2] = 0.31957894;
    p_faixas[3] = 0.15781574;
    p_faixas[4] = 0.03801476;
    for ( j = 0; j < 5; j++) for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
    for (i = 1; i < 5; i++)p_faixas[i] += p_faixas[i-1];
    

    Return.faixa = calloc(N,sizeof(int));

    for (i = 0; i < N; i++){
        j = 0;
        r = genrand64_real1();
        while(r > p_faixas[j]) j++;
        n_faixa[j]++;
    }
    //for (i = 0; i < 4; i++) n_faixa[i+1] += n_faixa[i]; 

    j = 0;
    int s = n_faixa[0];

    for (i = 0; i < N; i++){
        if(i==s) {j++;s+= n_faixa[j];}
        Return.faixa[i] = j;
    }

    int k = 0;
    double x = 0;

    for (i = 0; i < N; i++){

        k = 0;

        Return.degree[i] = (int*) calloc(5,sizeof(int));

        while(k == 0)k = empiric_distribution(distribution[Return.faixa[i]]);

        grau[i] = k;

        generate_multinomial(k, 5, constant[Return.faixa[i]], Return.degree[i]);

    }

    //Return.degree = M_bubbleSort(Return.degree,N,5);
    //bubbleSort_by(Return.faixa, grau, N);

    for(i = 0;i < 5;i++){
        //media_grau[i] /= n_faixa[i];
        free(constant[i]);
        free(distribution[i]);
    }

    //print_vetor(media_grau,5,sizeof(double));

    free(constant);
    free(distribution);
    free(p_faixas);
    //free(media_grau);
    free(n_faixa);
    free(grau);
    return Return;
}

igraph_t local_configuration_model(int N, double p,int seed){

    struct Graph G;
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    int i = 0,j = 0;
    for (i = 0; i < N; i++){ 
        G.viz[i] =(int*) malloc(1*sizeof(int));
        G.viz[i][0] = 0;
    }
    
    G.edges = 0;

    int** degree;
    int* faixa;
    struct Retorno Return = generate_degree(N,seed);
    faixa = Return.faixa;
    degree = Return.degree;
    
    int *shuff = (int*) malloc(N*sizeof(int));
    

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++) shuff[j] = j;
        swap(&shuff[i],&shuff[N-1]);
        shuff = randomize(shuff,N-1,seed+i);

        G = local_add_edge(G,degree,faixa,shuff,N-1,i,p,&edges);
    }
    

    int** ligacoes_restantes = (int **)malloc(5*sizeof(int*));
    int ligacoes_total = 0;
    for (i = 0; i < 5; i++) ligacoes_restantes[i] = (int *)calloc(5,sizeof(int));
    int n = 0;
    int* ids = (int*) malloc(0*sizeof(int));

    for (i = 0; i < G.Nodes; i++){
        if(somatorio(degree,i,5) != 0){
            ligacoes_total += somatorio(degree,i,5);
            for (j = 0; j < 5; j++) ligacoes_restantes[j][faixa[i]] += degree[i][j];
            ids = realloc(ids,sizeof(int)*(n+1));
            ids[n] = i;
            n++;
        }
    }
    
    char file_names[200] = "./dados/multi_probability.txt";
    double** constant = double_get_array(5,file_names);
    double** distribution = distribution_faixas();
    double* p_faixas = malloc(5*sizeof(double));
    double r;

    p_faixas[0] = 0.34393182;
    p_faixas[1] =  0.14065873;
    p_faixas[2] = 0.31957894;
    p_faixas[3] = 0.15781574;
    p_faixas[4] = 0.03801476;

    for ( j = 0; j < 5; j++) for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
    for (i = 1; i < 5; i++) p_faixas[i] += p_faixas[i-1];
    int faixas = 0;
    while(ligacoes_total != 0){

        // Cria ligações
        int k = 0;
        faixas = 0;
        int* vetor = calloc(5,sizeof(int));
        init_genrand64(seed);
        seed++;
        r = genrand64_real1();
        while(r > p_faixas[faixas]) faixas++;
        while(k == 0)k = empiric_distribution(distribution[faixas]);
        
        if(ligacoes_total - k < 0){
            free(vetor);
            continue;
        }
        generate_multinomial(k, 5, constant[faixas], vetor);
        int menor = k;
        for ( i = 0; i < 5; i++){
            vetor[i] = vetor[i] < ligacoes_restantes[faixas][i] ? vetor[i] : ligacoes_restantes[faixas][i];
            if(ligacoes_restantes[faixas][i] < menor) if(ligacoes_restantes[faixas][i] != 0) menor = ligacoes_restantes[faixas][i];
        }
        if(vetor_soma(vetor,5) == 0){
            free(vetor);
            continue;
        }
        while(menor> vetor_soma(vetor,5)){
            double r = genrand64_real1();
            for(i = 0; i < 5; i++) if(r < constant[faixas][i]) break;
            if(vetor[i] != 0) vetor[i]++;
        }
        k = vetor_soma(vetor,5);
    
        // Adicionando o novo nó
        G.viz =(int**) realloc(G.viz,(G.Nodes+1)*sizeof(int*));
        G.viz[G.Nodes] = (int*) calloc(1,sizeof(int));
        G.viz[G.Nodes][0] = 0;
        degree = realloc(degree,(G.Nodes+1)*sizeof(int*));
        degree[G.Nodes] = calloc(5,sizeof(int));
        for (i = 0; i < 5; i++) degree[G.Nodes][i] = vetor[i];

        faixa = realloc(faixa,(G.Nodes+1)*sizeof(int));
        faixa[G.Nodes] = faixas;

        ids = randomize(ids,n,seed);
        
        // Adicionar ligações
        G = local_add_edge(G,degree,faixa,ids,n,G.Nodes,p,&edges);
        
        for (i = 0; i < 5; i++) ligacoes_restantes[faixas][i] -= vetor[i] - degree[G.Nodes][i];
        ligacoes_total = 0;
        for (i = 0; i < 5; i++) ligacoes_total += somatorio(ligacoes_restantes,i,5);
        G.Nodes++;

        free(vetor);
    }

    igraph_t Grafo;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_empty(&Grafo, G.Nodes, IGRAPH_UNDIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    igraph_vector_t faixa_;
    igraph_vector_init(&faixa_, G.Nodes);

    for (i = 0; i < G.Nodes; i++){
        free(G.viz[i]);
        VECTOR(faixa_)[i] = faixa[i];
        
    }
    char att[200] = "faixas";
    igraph_cattribute_VAN_setv(&Grafo,"faixa",&faixa_);
    for (i = 0; i < G.Nodes; i++) free(degree[i]);
    free(G.viz);
    free(degree);
    free(faixa);
    free(shuff);
    free(p_faixas);
    for(i = 0;i < 5;i++){
        free(constant[i]);
        free(distribution[i]);
    }
    free(constant);
    free(distribution);
    free(ids);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&faixa_);
    return Grafo;
}


void calcula_propriedades(igraph_t *Grafo,double p, double *resultados,double* prob,double* std_vector,double* k_means) {
    igraph_vector_int_t v;
    double clustering;
    double l;
    double d;
    double* idades = (double *)calloc(5,sizeof(double));
    double* k_mean = (double *)calloc(5,sizeof(double));
    
    int i;

    int n = igraph_vcount(Grafo);
    igraph_vector_t faixas;
    igraph_vector_init(&faixas, n);
    igraph_vector_int_init(&v, n);

    igraph_degree(Grafo, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);
    igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);
    FILE *file;
    char filecheck[800];
    sprintf(filecheck,"./output/modelo/k_idades_%.2f.txt",p);
    file = fopen(filecheck,"a");
    FILE *arquivo;
    sprintf(filecheck,"./output/modelo/matrix_%.2f.txt",p);
    arquivo = fopen(filecheck,"a");
    for (i = 0; i < n; i++){

        igraph_vector_int_t vizinhos;
        igraph_vector_int_init(&vizinhos, 0);
        igraph_neighbors(Grafo, &vizinhos, i,IGRAPH_ALL);
        int* matrix = (int *)calloc(5,sizeof(int));
        for (int j = 0; j < igraph_vector_int_size(&vizinhos); j++){
            int vizinho = VECTOR(vizinhos)[j];
            matrix[(int)VECTOR(faixas)[vizinho]]++;
        }
        fprintf(arquivo,"%d %d %d %d %d %d\n",(int)VECTOR(faixas)[i],matrix[0],matrix[1],matrix[2],matrix[3],matrix[4]);
        free(matrix);

        fprintf(file,"%ld %d\n",igraph_vector_int_size(&vizinhos), (int)VECTOR(faixas)[i]);
        idades[(int)VECTOR(faixas)[i]]++;
        k_mean[(int)VECTOR(faixas)[i]] += VECTOR(v)[i];

        igraph_vector_int_destroy(&vizinhos);
    }
    fclose(file);
    igraph_vector_int_sort(&v);
    fclose(arquivo);
    double N0 = 0;
    if(n%2 !=0) N0 = VECTOR(v)[(n-1)/2];
    else N0 = (VECTOR(v)[n/2] + VECTOR(v)[n/2+1])*0.5;
    double media = (double)igraph_vector_int_sum(&v)/n;
    igraph_vector_int_mul(&v,&v);
    double media2 = (double)igraph_vector_int_sum(&v)/n;
    double std = sqrt(media2 - pow(media,2));
    igraph_transitivity_avglocal_undirected(Grafo,&clustering,IGRAPH_TRANSITIVITY_ZERO);
    igraph_vector_int_destroy(&v);

    for ( i = 0; i < 5; i++){
        k_mean[i] /= idades[i];
        k_means[i] += k_mean[i];
        prob[i] += idades[i]/n;
    }
    
    igraph_average_path_length(Grafo, &l, NULL, IGRAPH_UNDIRECTED, 1);
    igraph_diameter(Grafo, &d, 0, 0, 0, 0, IGRAPH_UNDIRECTED, 1);

    resultados[0] += media;
    resultados[1] += N0;
    resultados[2] += std;
    resultados[3] += clustering;
    resultados[4] += l;
    resultados[5] += d;
    std_vector[0] += pow(media,2);
    std_vector[1] += pow(N0,2);
    std_vector[2] += pow(std,2);
    std_vector[3] += pow(clustering,2);
    std_vector[4] += pow(l,2);
    std_vector[5] += pow(d,2);
    free(idades);
    igraph_vector_destroy(&faixas);
}


void generate_local_configuration_model(double p, int redes,int seed){

    int N = 2029;
    int i;

    double* resultados = (double *)calloc(6,sizeof(double));
    double* std = (double *)calloc(6,sizeof(double));
    double* prob = (double *)calloc(5,sizeof(double));
    double* k_mean = (double *)calloc(5,sizeof(double));
    double clustering;
    for (i = 0; i < redes; i++){

        if(redes!=1) printf("\e[1;1H\e[2J");

        igraph_t G;
        G = local_configuration_model(N,p,seed+i);
        //igraph_transitivity_avglocal_undirected(&G,&clustering,IGRAPH_TRANSITIVITY_ZERO);
        //resultados[3] += clustering;
        //std[3] += pow(clustering,2);
        if(redes!=1) calcula_propriedades(&G,p,resultados,prob,std,k_mean);
        printf("%d\n",i+1);
        igraph_destroy(&G);
        
    }
   
    for (i = 0; i < 6; i++){
        resultados[i] /= redes;
        std[i] /= redes;
        std[i] = sqrt(std[i] - pow(resultados[i],2));
    }
    for (i = 0; i < 5; i++){
        prob[i] /= redes;
        k_mean[i] /= redes;
    }
    
    if(redes > 1){
        FILE *file;
        char filecheck[800];
        sprintf(filecheck,"./output/modelo/resultados.txt");
        file = fopen(filecheck,"a");
        fprintf(file,"================%.3f================\n",p);
        fprintf(file,"%f %f %f %f %f\n",k_mean[0],k_mean[1],k_mean[2],k_mean[3],k_mean[4]);
        fprintf(file,"%f %f %f %f %f\n",prob[0],prob[1],prob[2],prob[3],prob[4]);
        fprintf(file,"%f(%.2f) %f(%.2f) %f(%.2f) %f(%.2f) %f(%.2f) %f(%.2f) \n",resultados[0],std[0],resultados[1],std[1],resultados[2],std[2],resultados[3],std[3],resultados[4],std[4],resultados[5],std[5]);
        fclose(file);

        /* FILE *arquivo;
        arquivo = fopen("./output/modelo/clustering.txt","a");
        fprintf(arquivo,"%f %f %f\n",p,resultados[3],std[3]);
        fclose(arquivo); */
    }
    free(resultados);
    free(std);
    free(prob);
    printf("Terminou\n");
}