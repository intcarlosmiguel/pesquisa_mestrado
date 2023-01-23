#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "mtwister.h"
#include "calc.h"
#include <math.h>

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
    int n = 0;
    double local = 0;
    for (int i = 0; i < N; i++){
        c = combination(viz[i][0],2);
        if(c!=0) for (int first = 1; first < viz[i][0]-1; first++) for (int second = first+1; second < viz[i][0]; second++) for(int vizinhos = 1; vizinhos<viz[second][0]; vizinhos++) if(viz[i][first] == viz[second][vizinhos]) {n++;break;}
        if(c!=0)local += (double) n/c;
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
        if(c!=0) for (int first = 1; first < viz[i][0]-1; first++) for (int second = first+1; second < viz[i][0]; second++) for(int vizinhos = 1; vizinhos<viz[second][0]; vizinhos++) if(viz[i][first] == viz[second][vizinhos]) {n++;break;}
        clustering[i] = (c!=0)? n/c: 0;
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
    while ((current<=tam) && (tam != N)){

        int sitio = lista[current];
        int vizinhos = viz[sitio][0];
        for(int j = 0; j<vizinhos;j++){
            int vizinho = viz[sitio][j+1];
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

void degree_distribution(int N,int** viz,int edges,double* media1,double* median1,double* std1){

    double media = 0;
    double media2 = 0;
    int* degree_ = (int*) malloc(sizeof(int*)*N);
    int median = (N+1)/2 - 1;
    media = (double)2*(edges)/N;

    for (int i = 0; i < N; i++){
        degree_[i] = viz[i][0];
        media2 += pow(degree_[i]-media,2);
    }
    
    degree_ = bubble_sort(degree_,N);

    *media1 += media;
    *std1 += pow(media2/(N-1),0.5);
    *median1 += degree_[median];
    
}

int** conf_model_p(int** viz,int* degree,int** mat,int ego,double p,int* edges){
    if(viz[ego][0] > 1){

        int N = viz[ego][0];
        int *vizinhos = (int*) malloc(N*sizeof(int));
        int *shuff = (int*) malloc(N*sizeof(int));

        for (int i = 0; i < N; i++) vizinhos[i] = viz[ego][i+1];
        //vizinhos = bubble_sort_by(vizinhos,deg,N);
        //print_vetor(vizinhos,N-1);
        for (int i = 0; i < N; i++){

            int site = vizinhos[i];
            //printf("%d\n",site);
            for (int j = 0; j < N; j++) shuff[j] = vizinhos[j];
            shuff = randomize(shuff,N,i+site);
            shuff = ending(shuff,N,site);
            //printf("Erro é antes do for!,%d,%d\n",degree[site],site);
            for (int j = 0; j < N-1; j++){

                if(degree[site] == 0) break;
                //printf("Erro no vizinho:%d\n",shuff[j]);
                int vizinho = shuff[j];
                //printf("Erro no vizinho:%d,%d\n",vizinho,degree[vizinho]);
                double rand = genrand64_real2();
                if(rand<=p){
                    int rep = 0;
                    for (int k = 0; k < *edges; k++){

                        if((mat[k][0] == site) && (mat[k][1] == vizinho)){
                            rep = 1;
                            break;
                        }

                        if((mat[k][0] == vizinho) && (mat[k][1] == site)){
                            rep = 1;
                            break;
                        }

                    }


                    if((rep == 1) || (degree[vizinho] == 0)) continue;

                    mat = (int**) realloc(mat,(*edges+1)*sizeof(int*));
                    mat[*edges] = (int*) malloc(2* sizeof(int));
                    mat[*edges][0] = site;
                    mat[*edges][1] = vizinho;
                    *edges += 1;

                    degree[site]--;
                    degree[vizinho]--;

                    viz[site][0]++;
                    viz[vizinho][0]++;
                    viz[site] = (int*) realloc(viz[site],(viz[site][0]+1)*sizeof(int));
                    viz[vizinho] = (int*) realloc(viz[vizinho],(viz[vizinho][0]+1)*sizeof(int));
                    viz[site][viz[site][0]] = vizinho;
                    viz[vizinho][viz[vizinho][0]] = site;
                }


            }
            
        }
        free(vizinhos);
        free(shuff);
    }
    return mat;
}

int** add_edge(int** viz,int* degree,int** mat,int* shuff,int* edges, int n,int site,double p){

    for (int i = 0; i < n; i++){

        int rep = 0;
        int vizinho = shuff[i];

        if(degree[site] == 0){
            if(site == 0) for (int j = 0; j < viz[site][0]; j++) printf("%d,",degree[viz[site][j+1]]);
            if(site == 0)printf("\n");
            mat = conf_model_p(viz,degree,mat,site,p,edges);
            if(site == 0) for (int j = 0; j < viz[site][0]; j++) printf("%d,",degree[viz[site][j+1]]);
            break;
        }

        for (int j = 0; j < *edges; j++){
            if((mat[j][0] == site) && (mat[j][1] == vizinho)){
                rep = 1;
                break;
            }

            if((mat[j][0] == vizinho) && (mat[j][1] == site)){
                rep = 1;
                break;
            }

        }


        if((rep == 1) || (degree[vizinho] == 0)) continue;
        
        mat = (int**) realloc(mat,(*edges+1)*sizeof(int*));
        mat[*edges] = (int*) malloc(2* sizeof(int));
        mat[*edges][0] = site;
        mat[*edges][1] = vizinho;
        *edges += 1;
        
        degree[site]--;
        degree[vizinho]--;
        
        viz[site][0]++;
        viz[vizinho][0]++;
        viz[site] = (int*) realloc(viz[site],(viz[site][0]+1)*sizeof(int));
        viz[vizinho] = (int*) realloc(viz[vizinho],(viz[vizinho][0]+1)*sizeof(int));
        viz[site][viz[site][0]] = vizinho;
        viz[vizinho][viz[vizinho][0]] = site;
    }
    return mat;
}

int** configuration_model(int** viz,int* degree,int N, double p,int* size,int seed){

    int edges = 0;
    init_genrand64(seed);
    degree = bubble_sort(degree,N);
    int **mat = (int **)malloc(0* sizeof(int*));
    int *shuff = (int*) malloc(N*sizeof(int));

    for (int i = 0; i < N; i++){

        for (int j = 0; j < N; j++) shuff[j] = j;
        shuff = randomize(shuff,N,seed);
        shuff = ending(shuff,N, i);
        //printf("%d\n",i);
        mat = add_edge(viz,degree,mat,shuff,&edges,N-1,i,p);
    }
    free(shuff);
    *size = edges;
    return mat;
}

int* get_degree(int N){
    FILE* file;
    file = fopen("./degree.txt","r");
    int* degree = (int*) malloc(sizeof(int)*N);
    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&degree[i]));
    fclose(file);
    return degree;
}

void generate_configuration_model(double p, int T){

    int N = size_txt();
    int edges = 0;
    double media = 0;
    double median = 0;
    double std = 0;
    double as = 0;
    double l = 0;
    double r2 = 0;
    double diametro = 0;
    double d = 0;

    for (int i = 0; i < T; i++){
        if(T!=1) printf("\e[1;1H\e[2J");
        int* degree = get_degree(N);
        edges = 0;
        int** viz = (int **)malloc(N*sizeof(int*));

        for (int i = 0; i < N; i++){
            viz[i] =(int*) malloc(1*sizeof(int));
            viz[i][0] = 0;
        }

        int** mat = configuration_model(viz,degree,N,p,&edges,i);
        //print_matrix(viz,N);
        as += avarage_clustering(viz,N);
        l += av_path_length(viz, N,&d);
        diametro += d;

        double *deg = degree_list(viz,N);
        double *clustering = list_clustering(viz,N);
        //r2 += correlation(deg,clustering,N);

        printf("%d\n",i+1);
        degree_distribution(N,viz,edges,&media,&median,&std);
        for(int j = 0; j < edges; j++) free(mat[j]);
        for(int j = 0; j < N; j++) free(viz[j]);
        free(mat);
        free(viz);
        free(degree);
        free(deg);
        free(clustering);
    }
    //if(T!=1) printf("\e[1;1H\e[2J");
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",media/T,median/T,std/T,as/T,l/T,diametro/T);
}

int main(int argc,char *argv[ ]){
    int T = atoi(argv[1]);
    generate_configuration_model(1.0,T);
}
