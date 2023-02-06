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
    int **mat;
    int **viz;
    int **n_existir;
    int Nodes;
    int edges;
    int existir;
    int* infect;
    int* infect2;
    int** degree;
    int* faixa;
};




double global_clustering(int** viz,int N){
    int c = 0;
    int n = 0;
    for (int i = 0; i < N; i++){
        c += combination(viz[i][0],2);
        if(combination(viz[i][0],2)!=0) for (int first = 1; first < viz[i][0]-1; first++) for (int second = first+1; second < viz[i][0]; second++) for(int vizinhos = 1; vizinhos<viz[second][0]; second++) if(viz[i][first] == viz[second][vizinhos]) {n++;break;}
    }
    return(double) n/c;
}

double avarage_clustering(int** viz,int N){
    int c = 0;
    double n = 0;
    double local = 0;
    for (int i = 0; i < N; i++){
        c = combination(viz[i][0],2);
        if(c!=0)
            for (int first = 1; first <=viz[i][0]; first++)
                for (int second = 1; second <=viz[i][0]; second++)
                    if(second!= first)
                        for(int vizinhos = 1; vizinhos<=viz[viz[i][second]][0]; vizinhos++) 
                            if(viz[i][first] == viz[viz[i][second]][vizinhos]) 
                                {n++;break;}
        if(c!=0)local += (double) n/(2*c);
        n = 0;
    }
    return(double) local/N;
}

double* list_clustering(int** viz,int N){
    double c = 0;
    double n = 0;
    double* clustering = (double*) malloc(N*sizeof(double));
    for (int i = 0; i < N; i++){
        c = combination(viz[i][0],2);
        if(c!=0) 
            for (int first = 1; first <=viz[i][0]; first++) 
                for (int second = 1; second <=viz[i][0]; second++) 
                    if(second!= first)
                        for(int vizinhos = 1; vizinhos<=viz[viz[i][second]][0]; vizinhos++) 
                            if(viz[i][first] == viz[viz[i][second]][vizinhos]) 
                                {n++;break;}
        clustering[i] = (c!=0)? n/(2*c): 0;
        n = 0;
    }
    return clustering;
}

double shortest_length(int** viz,int N,int site,double* diametro){

    double l = 0;
    int tam = 1;
    int infinity = N+2;
    int* lista = (int*) malloc(1*sizeof(int));
    int* distance = (int*) malloc(N*sizeof(int));
    int current = 0;

    for (int i = 0; i < N; i++) distance[i] = infinity;
    lista[0] = site;
    distance[site] = 0;
    //if(s == 76)printf("%d\n",site);
    while ((current<tam) && (tam != N)){

        int sitio = lista[current];
        int vizinhos = viz[sitio][0];
        for(int j = 1; j<=vizinhos;j++){
            int vizinho = viz[sitio][j];
            //if(s == 76)printf("%d,%d - %d/%d\n",sitio,vizinho,current,tam);
            if((distance[vizinho]>distance[sitio]+1)){
                distance[vizinho] = distance[sitio]+1;
                tam++;
                lista = (int*) realloc(lista,tam*sizeof(int));
                lista[tam-1] = vizinho;
                if(distance[sitio]+1 > *diametro) *diametro = (double) distance[sitio]+1;
            }

        }
        current++;
    }
    for (int i = 0; i < N; i++) if(distance[i] != infinity) l += distance[i];
    free(distance);
    free(lista);
    return (double) l/(N-1);
}

double av_path_length(int** viz,int N,double* diametro){
    double l = 0;
    for(int i = 0; i < N; i++) l += shortest_length(viz,N,i,diametro);
    return l/N;
}

double* degree_list(int** viz,int N){
    double* degree = (double*) malloc(sizeof(double*)*N);
    for (int i = 0; i < N; i++) degree[i] = (double) viz[i][0];
    return degree;
}

void degree_distribution(struct Graph G,double* media1,double* median1,double* std1){

    double media = 0;
    double media2 = 0;
    int* degree_ = (int*) malloc(sizeof(int*)*G.Nodes);
    int median = (G.Nodes+1)/2 - 1;
    media = (double)2*(G.edges)/G.Nodes;

    for (int i = 0; i < G.Nodes; i++){
        degree_[i] = G.viz[i][0];
        media2 += pow(degree_[i]-media,2);
    }
    
    degree_ = bubble_sort(degree_,G.Nodes);
    *media1 = media;
    *std1 = pow(media2/(G.Nodes-1),0.5);
    *median1 = degree_[median];

    free(degree_);
}

int check_existence(int** outro, int N,int site, int vizinho){
    for (int i = 0; i < N; i++){
        if((outro[i][0] == site) && (outro[i][1] == vizinho)) return 1;

        if((outro[i][0] == vizinho) && (outro[i][1] == site)) return 1;
    }
    return 0;
}

void infecao_si(struct Graph G,int* infect,int* infect2,double beta,int* tempo, int time){
    int i;
    for (i = 0; i < G.Nodes; i++){
        if(infect[i] == 0){
            int num_infect = 0;
            for (int j = 1; j < G.viz[i][0]+1; j++) num_infect += (infect[G.viz[i][j]] == 1)? 1: 0;
            double probability = (1 - pow(1 - beta,num_infect));
            infect2[i] = (genrand64_real2() < probability)? 1: 0;
        }
        else  infect2[i] = 1;
    }
    for (i = 0; i < G.Nodes; i++) tempo[time] += (G.infect2[i] == 1)? 1: 0;
}

struct Graph init(struct Graph G, double probability,int time,double beta){

    int i;
    G.infect = (int*) malloc(G.Nodes*sizeof(int));
    int* tempo = (int*) malloc(time*sizeof(int));
    for ( i = 0; i < time; i++) tempo[i] = 0;

    for (i = 0; i < G.Nodes; i++){
        G.infect[i] = (genrand64_real2() < probability)? 1: 0;
        G.infect2[i] = G.infect[i];
    }
    for (int i = 0; i < time; i++) ((i+1)%2 != 0 )? infecao_si(G,G.infect,G.infect2,beta,tempo,i) : infecao_si(G,G.infect2,G.infect,beta,tempo,i);
    
    return G;
}

int somatorio(int** array,int elemento,int N){
    int soma = 0;
    for (int i = 0; i < N; i++) soma += array[elemento][i];
    return soma;
}

struct Graph conf_model_p(struct Graph G,int ego,double p){
    int N = 0;
    int *vizinhos = (int*) malloc(0*sizeof(int));

    for (int i = 0; i < G.viz[ego][0]; i++){
        int vizinho = G.viz[ego][i+1];
        if(somatorio(G.degree,vizinho,N) !=0){
            vizinhos = (int*) realloc(vizinhos,(N+1)*sizeof(int));
            vizinhos[N] = vizinho;
            N++;
        }
    }
    //for (int i = 0; i < N; i++) printf("%d,",vizinhos[i]);
    //printf("Tamanho: %d\n",N);
    if(N > 1){

        int *shuff = (int*) malloc(N*sizeof(int));

        for (int i = 0; i < N; i++){
            int site = vizinhos[i];

            for (int j = 0; j < N; j++) shuff[j] = vizinhos[j];
            shuff = randomize(shuff,N,i+site);
            shuff = ending(shuff,N,site);
            //printf("Sitio: %d\n",site);
            for (int j = 0; j < N-1; j++){
                if(somatorio(G.degree,site,N) == 0) break;

                int vizinho = shuff[j];
                double rand = genrand64_real2();
                //printf("Vizinho: %d,%f\n",vizinho,rand);                
                if(rand<=p){

                    int rep = check_existence(G.mat,G.edges,site,vizinho);
                    if((rep == 1) || (G.degree[vizinho][G.faixa[site]] == 0) || (G.degree[site][G.faixa[vizinho]] == 0) ) continue;
                    rep = check_existence(G.n_existir,G.existir,site,vizinho);
                    if(rep == 1) continue;
                    int faixa1 = G.faixa[site];
                    int faixa2 = G.faixa[vizinho];

                    G.mat = (int**) realloc(G.mat,(G.edges+1)*sizeof(int*));
                    G.mat[G.edges] = (int*) malloc(2* sizeof(int));
                    G.mat[G.edges][0] = site;
                    G.mat[G.edges][1] = vizinho;
                    G.edges += 1;

                    G.degree[site][faixa2]--;
                    G.degree[vizinho][faixa1]--;

                    G.viz[site][0]++;
                    G.viz[vizinho][0]++;
                    G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
                    G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
                    G.viz[site][G.viz[site][0]] = vizinho;
                    G.viz[vizinho][G.viz[vizinho][0]] = site;
                }

                else{
                    //printf("%d\n",*n);
                    //if(*n > 0)printf("Teste: %p,%p\n",n_existe[0],n_existe);
                    int rep = check_existence(G.n_existir,G.existir,site,vizinho);
                    if(rep == 0){
                        G.n_existir = (int**) realloc(G.n_existir,(G.existir+1)*sizeof(int*));
                        G.n_existir[G.existir] = (int*) malloc(2* sizeof(int));
                        G.n_existir[G.existir][0] = site;
                        G.n_existir[G.existir][1] = vizinho;
                        G.existir += 1;
                    }
                }
            }
        }
        free(shuff);
    }
    free(vizinhos);
    return G;
}

struct Graph add_edge(struct Graph G,int* shuff, int n,int site,double p){

    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(somatorio(G.degree,site,5) == 0){
            //if(p>0) G = conf_model_p(G,site,p);
            break;
        }
        int rep = check_existence(G.mat,G.edges,site,vizinho);


        if((rep == 1) || (G.degree[vizinho][G.faixa[site]] == 0) || (G.degree[site][G.faixa[vizinho]] == 0) ) continue;
        int faixa1 = G.faixa[site];
        int faixa2 = G.faixa[vizinho];
        
        G.mat = (int**) realloc(G.mat,(G.edges+1)*sizeof(int*));
        G.mat[G.edges] = (int*) malloc(2* sizeof(int));
        G.mat[G.edges][0] = site;
        G.mat[G.edges][1] = vizinho;
        G.edges += 1;

        G.degree[site][faixa2]--;
        G.degree[vizinho][faixa1]--;
        
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
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){ 
        G.viz[j] =(int*) malloc(1*sizeof(int));
        G.viz[j][0] = 0;
    }
    G.mat = (int **)malloc(0*sizeof(int*));
    G.edges = 0;
    G.existir = 0;
    G.n_existir = (int **)malloc(0* sizeof(int*));
    G.degree = get_degree(N);
    G.faixa = get_faixas(N);

    init_genrand64(seed);
    //degree = bubble_sort(degree,N);
    
    int *shuff = (int*) malloc(N*sizeof(int));

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) shuff[j] = j;
        shuff = randomize(shuff,N,seed);
        shuff = ending(shuff,N, i);
        G = add_edge(G,shuff,N-1,i,p);
    }
    free(G.n_existir);
    free(shuff);
    return G;
}

void create_network(struct Graph G,double p){
    char arquivo[20];
    sprintf(arquivo, "dados_%.2f.txt", p);
    FILE* file;
    file = fopen(arquivo,"w");
    for (int i = 0; i < G.edges; i++) fprintf(file,"%d\t%d\n",G.mat[i][0],G.mat[i][1]);
    
    fclose(file);
}

void result(struct Graph G,double p,double* resultados){

    resultados[3] = avarage_clustering(G.viz,G.Nodes);
    resultados[5] = av_path_length(G.viz, G.Nodes,&resultados[6]);

    double *deg = degree_list(G.viz,G.Nodes);
    double *clustering = list_clustering(G.viz,G.Nodes);

    resultados[4] = correlation(deg,clustering,G.Nodes);
    degree_distribution(G,&resultados[0],&resultados[1],&resultados[2]);

    free(deg);
    free(clustering);
}

int generate_configuration_model(double p, int T){

    int N = size_txt();
    double** resultados = (double **)malloc(T*sizeof(double*));
    #pragma omp parallel for
    for (int i = 0; i < T; i++){

        if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        G = configuration_model(N,p,i);
        resultados[i] = (double*) malloc(7*sizeof(double));
        result(G,p,resultados[i]);
        printf("%d\n",i+1);
        if(T == 1) create_network(G,p);
        
        for(int j = 0; j < G.edges; j++) free(G.mat[j]);
        for(int j = 0; j < N; j++){free(G.viz[j]);free(G.degree[j]);}
        
        free(G.mat);
        free(G.viz);
        free(G.degree);
        free(G.faixa);
    }
    double media = 0;
    double median = 0;
    double std = 0;
    double as = 0;
    double l = 0;
    double r2 = 0;
    double diametro = 0;

    double media2 = 0;
    double median2 = 0;
    double std2 = 0;
    double as2 = 0;
    double l2 = 0;
    double r22 = 0;
    double diametro2 = 0;

    for (int i = 0; i < T; i++){

        media += resultados[i][0];
        media2 += pow(resultados[i][0],2);

        median += resultados[i][1];
        median2 += pow(resultados[i][1],2);

        std += resultados[i][2];
        std2 += pow(resultados[i][2],2);

        as += resultados[i][3];
        as2 += pow(resultados[i][3],2);

        r2 += resultados[i][4];
        r22 += pow(resultados[i][4],2);

        l += resultados[i][5];
        l2 += pow(resultados[i][5],2);

        diametro += resultados[i][6];
        diametro2 += pow(resultados[i][6],2);

    }

    media2 = pow(media2/T - pow(media/T,2),0.5);
    median2 = pow(median2/T - pow(median/T,2),0.5);;
    std2 = pow(std2/T - pow(std/T,2),0.5);
    as2 = pow(as2/T - pow(as/T,2),0.5);
    l2 = pow(l2/T - pow(l/T,2),0.5);
    r22 = pow(r22/T - pow(r2/T,2),0.5);
    diametro2 = pow(diametro2/T - pow(diametro/T,2),0.5);

    //printf("%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",p,media/T,median/T,std/T,as/T,r2/T,l/T,diametro/T);
    FILE* file;
    file = fopen("./resultados.txt","a");
    fprintf(file,"%f\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\n",p,media/T,media2,median/T,median2,std/T,std2,as/T,as2,r2/T,r22,l/T,l2,diametro/T,diametro2);
    fclose(file);
}

int main(int argc,char *argv[ ]){
    int T = atoi(argv[1]);
    //generate_configuration_model((double) 0,T);
    generate_configuration_model((double) 0.5,T);
    //generate_configuration_model((double) 1.0,T);
}
