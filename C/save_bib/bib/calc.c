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

void randomize (int* arr, int n,int seed){
    init_genrand64(seed);
    for (int i = n - 1; i > 0; i--) {
        int j = genrand64_int64() % (i + 1);
        
        // Troca os elementos usando ponteiros
        int temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}
int geometry(double lambda){
    double n = -1*(1/lambda)*log(genrand64_real1()/(exp(lambda) - 1));
    if(n <0) return geometry(lambda);
    return (int) n;
}

int potentia(double lambda,double A){
    double n = pow(A/genrand64_real1(),1/lambda);
    if((int)n <0) return potentia(lambda,A);
    return (int) n;
}

void print_vetor(void* array,int N,int check){
    if(check == sizeof(int)){
        int* intArray = (int*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%d ",intArray[i]);
            else printf("%d\n",intArray[i]);
        }
    }
    if(check == sizeof(double)){
        double* doubleArray = (double*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%f ",doubleArray[i]);
            else printf("%f\n",doubleArray[i]);
        }
    }
}

void bubbleSort(int *arr, int n, void *ref,int check) {
    switch (check){
    case sizeof(int):
        int temp;
        int* ref2 = (int*) malloc(n*sizeof(int));
        for (int i = 0; i < n; i++) ref2[i] = *(int *)(ref + i * sizeof(int)); 
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - i - 1; j++) {
                if (ref2[j] > ref2[j + 1]) {
                    temp = ref2[j];
                    ref2[j] = ref2[j + 1];
                    ref2[j + 1] = temp;
                    temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
        free(ref2);
        break;
    case sizeof(double):
        double temp_;
        double* ref2_ = (double*) malloc(n*sizeof(double));
        for (int i = 0; i < n; i++) ref2_[i] = *(double *)(ref + i * sizeof(double)); 
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - i - 1; j++) {
                if (ref2_[j] > ref2_[j + 1]) {
                    temp_ = ref2_[j];
                    ref2_[j] = ref2_[j + 1];
                    ref2_[j + 1] = temp_;
                    temp_ = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp_;
                }
            }
        }
        free(ref2_);
        break;
    default:
        break;
    }
}

void sortIntByRef(int *arr, void *ref, int n,int check) {
    // cria vetor de índices
    int idx[n];
    for (int i = 0; i < n; i++) {
        idx[i] = i;
    }

    // ordena o vetor de índices com base no vetor de referência
    bubbleSort(arr, n, ref,check);

    // reorganiza o vetor original com base no vetor de índices ordenado
    int sortedArr[n];
    for (int i = 0; i < n; i++) {
        sortedArr[i] = arr[idx[i]];
    }
    for (int i = 0; i < n; i++) {
        arr[i] = sortedArr[i];
    }
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

int* get_faixas(int N,char* str){
    FILE* file;
    file = fopen(str,"r");
    int* faixas = (int*) malloc(sizeof(int)*N);

    for(int i = 0; i < N; i++) if(fscanf(file,"%d\n",&faixas[i]));
    fclose(file);
    return faixas;
}

int size_txt(char *str){
    FILE* f;
    int L = 0;
    char c;
    f = fopen(str,"r");
    for (c = getc(f); c != EOF; c = getc(f)) if (c == '\n') L = L + 1;
    return L;
}

void print_matrix(int** mat,int N,int n){
    printf("========================================================\n");
    for (int i = 0; i < N; i++) print_vetor(mat[i],n,sizeof(mat[i][0]));
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

    printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",media/T,median/T,std/T,as/T,r2/T,l/T,diametro/T);
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
        if (arr[j] > pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quickSort(int *arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);

        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void generate_file(char* filename,void* array,int linhas,int colunas,int check){

    FILE *file;
    file = fopen(filename,"w");
    for (int i = 0; i < linhas; i++){
        char print[100] = "";
        char resultado[100] = "";
        for(int j = 0;j<colunas;j++){
            switch (check){
                case sizeof(int)/* constant-expression */:
                    if(j != colunas-1) sprintf(print,"%d ",((int**)array)[i][j]);
                    else  sprintf(print,"%d\n",((int**)array)[i][j]);
                    strcat(resultado, print);
                    break;
                case sizeof(double):
                    if(j != colunas-1) sprintf(print,"%f ",((double**)array)[i][j]);
                    else  sprintf(print,"%f\n",((double**)array)[i][j]);
                    strcat(resultado, print);
                    break;
                default:
                    break;
            }
        }
        fprintf(file,"%s",resultado);
    }
    fclose(file);
}

void load_file(char* filename,void* array,int check){
    int linhas = size_txt(filename);
    FILE *file;
    file = fopen(filename,"r");
    if(check == sizeof(int)){
        int* intArray = (int*)array;
        for (int i = 0; i < linhas; i++) if(fscanf(file,"%d\n",&intArray[i]));
    }
    if(check == sizeof(double)){
        double* doubleArray = (double*)array;
        for (int i = 0; i < linhas; i++) if(fscanf(file,"%lf\n",&doubleArray[i]));
    }
    fclose(file);
}

void generate_multinomial(int n, int k, double *probabilities, int *outcomes) {
    int i,j;

    // Gera os k-1 primeiros valores
    for (i = n; i > 0; i-- ) {
        double r = genrand64_real1();
        for(j = 0; j < k; j++) if(r < probabilities[j]) break;
        outcomes[j]++;
    }

    // O último valor é determinado pelo resto
}

int generalized_geometry(double lambda,double A){
    double n = -1*(1/lambda)*log(genrand64_real1()/A);
    if(n <0) return generalized_geometry(lambda,A);
    return (int) n;
}

int empiric_distribution(double* distribution){
    double r = genrand64_real1();
    int n = 0;
    while(r > distribution[n]) n++;
    return n;
}

void bubbleSort_by(int* v, int* v2, int n) {
    int i, j, temp;

    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            // Comparamos os valores em v2
            if (v2[j] < v2[j+1]) {
                // Troca os elementos em v
                temp = v[j];
                v[j] = v[j+1];
                v[j+1] = temp;

                // Troca os elementos em v2
                temp = v2[j];
                v2[j] = v2[j+1];
                v2[j+1] = temp;
            }
        }
    }
}