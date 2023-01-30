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

struct Graph conf_model_p(struct Graph G,int* degree,int ego,double p){
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
                if(degree[site] == 0) break;
                int vizinho = shuff[j];
                double rand = genrand64_real2();
                //printf("Vizinho: %d,%f\n",vizinho,rand);                
                if(rand<=p){

                    int rep = check_existence(G.mat,G.edges,site,vizinho);
                    if((rep == 1) || (degree[vizinho] == 0)) continue;
                    rep = check_existence(G.n_existir,G.existir,site,vizinho);
                    if(rep == 1) continue;

                    G.mat = (int**) realloc(G.mat,(G.edges+1)*sizeof(int*));
                    G.mat[G.edges] = (int*) malloc(2* sizeof(int));
                    G.mat[G.edges][0] = site;
                    G.mat[G.edges][1] = vizinho;
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

struct Graph add_edge(struct Graph G,int* degree,int* shuff, int n,int site,double p){

    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        //printf("Vizinho: %d,%d\n",vizinho,degree[site]);
        if(degree[site] == 0){
            //for (int j = 0; j < viz[site][0]; j++) printf("%d,",viz[site][j+1]);
            if(p>0) G = conf_model_p(G,degree,site,p);
            //for (int j = 0; j < viz[site][0]; j++) printf("%d,",degree[viz[site][j+1]]);
            break;
        }
        //printf("Vizinho: %d\n",vizinho);
        int rep = check_existence(G.mat,G.edges,site,vizinho);


        if((rep == 1) || (degree[vizinho] == 0)) continue;
        
        G.mat = (int**) realloc(G.mat,(G.edges+1)*sizeof(int*));
        G.mat[G.edges] = (int*) malloc(2* sizeof(int));
        G.mat[G.edges][0] = site;
        G.mat[G.edges][1] = vizinho;
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

struct Graph configuration_model(int* degree,int N, double p,int seed){

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

    init_genrand64(seed);
    degree = bubble_sort(degree,N);
    
    int *shuff = (int*) malloc(N*sizeof(int));

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++) shuff[j] = j;
        shuff = randomize(shuff,N,seed);
        shuff = ending(shuff,N, i);
        //printf("Ego: %d\n",i);
        G = add_edge(G,degree,shuff,N-1,i,p);
    }
    free(G.n_existir);
    free(shuff);
    return G;
}

int* get_degree(int N){
    FILE* file;
    file = fopen("./degree.txt","r");
    int* degree = (int*) malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&degree[i]));
    fclose(file);
    return degree;
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

void generate_configuration_model(double p, int T){

    int N = size_txt();
    
    double** resultados = (double **)malloc(T*sizeof(double*));
    #pragma omp parallel for
    for (int i = 0; i < T; i++){

        if(T!=1) printf("\e[1;1H\e[2J");
        struct Graph G;
        int* degree = get_degree(N);
        G = configuration_model(degree,N,p,i);
        resultados[i] = (double*) malloc(7*sizeof(double));
        result(G,p,resultados[i]);
        printf("%d\n",i+1);
        if(T == 1) create_network(G,p);
        
        for(int j = 0; j < G.edges; j++) free(G.mat[j]);
        for(int j = 0; j < N; j++)free(G.viz[j]);
        
        free(G.mat);
        free(G.viz);
        free(degree);
    }

    double media = 0;
    double median = 0;
    double std = 0;
    double as = 0;
    double l = 0;
    double r2 = 0;
    double diametro = 0;

    for (int i = 0; i < T; i++){
        media += resultados[i][0];
        median += resultados[i][1];
        std += resultados[i][2];
        as += resultados[i][3];
        r2 += resultados[i][4];
        l += resultados[i][5];
        diametro += resultados[i][6];
    }
    FILE* file;
    file = fopen("./resultado.txt","a");
    fprintf(file,"%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",p,media/T,median/T,std/T,as/T,r2/T,l/T,diametro/T);
    fclose(file);
}

int main(int argc,char *argv[ ]){
    int T = atoi(argv[1]);
    //generate_configuration_model((double) i/100,T);
}
