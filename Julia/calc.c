#include <stdlib.h>
#include <stdio.h>
#include "mtwister.h"
#include <string.h>
#include <math.h>
#define M_PI 3.14159265358979323846

double media(double* x,int N){
    double media_ = 0;
    for (int i = 0; i < N; i++) media_ += x[i];
    return media_/N;
}

double variance(double* x,int N){
    double var = 0;
    double x_ = media(x,N);
    for (int i = 0; i < N; i++) var += pow(x[i]-x_,2);
    return var/N;
}

double correlation(double* x,double* y,int N){
    double x_ = media(x,N);
    double y_ = media(y,N);
    //printf("%f,%f\n",x_,y_);
    double xy = 0;
    double var_x = 0;
    double var_y = 0;

    for (int i = 0; i < N; i++){
        xy += x[i]*y[i];
        var_x += pow(x[i] - x_,2);
        var_y += pow(y[i] - y_,2);
    }
    double rho = (xy/N - x_*y_)/pow(var_x/N*var_y/N,0.5);
    //printf("%f\n",rho);
    return rho;
}

int fatorial(int n,int p){
    if(n!=p) return n*fatorial(n-1,p);
    else return 1;
}

int combination(int n, int p){
    if(n<=1) return 0;
    return fatorial(n,n-p)/fatorial(p,0);
}

void swap (int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

int* ending(int* array, int n, int site,int inicio){
    for (int i = inicio; i < n; i++){
        if(array[i] == site){
            swap(&array[i], &array[n-1]);
            break;
        }
    }
    return array;
}

int* randomize (int* array, int n,int seed){
    init_genrand64(seed);
    for (int i = n-1; i > 0; i--){
        int j = genrand64_int64() % (i+1);
        swap(&array[i], &array[j]);
    }
    return array;
}

void print_vetor(int* array,int N){
    for (int i = 0; i < N; i++){
        if(i!=N-1) printf("%d\t",array[i]);
        else printf("%d\n",array[i]);
    }
}

int* bubble_sort_by(int *array,int *valor,int N){
    
    for (int i = 0 ; i < ( N - 1 ); i++){
        for (int j= 0 ; j < N - i - 1; j++){
            if(array[j] < array[j+1]){
                int temp = array[j];
                array[j]   = array[j+1];
                array[j+1] = temp;
                int v = valor[j];
                valor[j] = valor[j+1];
                valor[j+1] = v;
            }
        }
    }
    
    return array;
}

int* bubble_sort(int* array,int N){
    
    for (int i = 0 ; i < ( N - 1 ); i++){
        for (int j= 0 ; j < N - i - 1; j++){
            if(array[j] < array[j+1]){
                int temp=array[j];
                array[j]   = array[j+1];
                array[j+1] = temp;
            }
        }
    }
    
    return array;
}

int* remove_(int *vetor,int elemento,int N){
    int i;
    int* copia = (int*) malloc(sizeof(int)*(N-1));
    for (i = 0; i < N-1; i++){
        if(i<elemento) copia[i] = vetor[i];
        else copia[i] = vetor[i+1];
    }
    free(vetor);
    return copia;
}

int** get_degree2(int N){
    FILE* file;
    file = fopen("./degrees.txt","r");
    int** degree = (int**) malloc(sizeof(int*)*N);

    for(int i = 0; i < N; i++) {
        degree[i] = (int*) malloc(sizeof(int)*5);
        if(fscanf(file,"%d\t%d\t%d\t%d\t%d\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

int* get_faixas(int N){
    FILE* file;
    file = fopen("./output/SBM/site_faixas_SBM.txt","r");
    int* faixas = (int*) malloc(sizeof(int)*N);

    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&faixas[i]));
    fclose(file);
    return faixas;
}

int size_txt(){
    FILE* f;
    int L = 0;
    char c;
    f = fopen("./dados/degree.txt","r");
    for (c = getc(f); c != EOF; c = getc(f)) if (c == '\n') L = L + 1;
    return L+1;
}

void print_matrix(int** mat,int N,int n){
    for (int i = 0; i < N; i++){
        for (int j= 0; j < n; j++) printf("%d,",mat[i][j]);
        printf("\n");
    }
    printf("========================================================\n");
}

int check_existence(int** outro, int N, int site, int vizinho){
    for (int i = 0; i < N; i++){
        if(outro[i][0] == site){
            if(outro[i][1] == vizinho) return 1;
        }
        else if(outro[i][0] == vizinho){
            if(outro[i][1] == site) return 1;
        }
    }
    return 0;
}

void generate_resultados(double** resultados, int T,char arquivo[]){
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
    /*FILE* file;

    char destination[100] = "./output/resultados_";
    strcat(destination,arquivo);

    char ext[100] = ".txt";
    strcat(destination,ext);

    file = fopen(destination,"a");
    fprintf(file,"%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\t%.2f(%.2f)\n",media/T,media2,median/T,median2,std/T,std2,as/T,as2,r2/T,r22,l/T,l2,diametro/T,diametro2);
    fclose(file); */
}

double exponentialRand(double lambda) {
    if(lambda == 0) return 0;
    return -log(1 - genrand64_real1()) / lambda;
}

double normalRand(double mean, double stdDev) {
    double z = sqrt(-2 * log(genrand64_real1())) * cos(2 * M_PI * genrand64_real1());
    return mean + stdDev * z;
}
int partition(int *arr, int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] >= pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quicksort(int *arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quicksort(arr, low, pi - 1);
        quicksort(arr, pi + 1, high);
    }
}