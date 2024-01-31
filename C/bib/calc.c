#include <stdlib.h>
#include <stdio.h>
#include "mtwister.h"
#include <string.h>
#include <math.h>
#include <igraph.h>

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
        char print[400] = "";
        char resultado[400] = "";
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

int* arange(int inicio, int fim, int passo){
    int N = (int)(fim - inicio)/passo+1;
    int* array = calloc(N,sizeof(int));
    for (int i = 0; i < N; i++) array[i] = inicio +passo*i;
    return array;
}

double** double_get_array(char* ptr){
    int linhas = size_txt(ptr);
    FILE* file;
    file = fopen(ptr,"r");
    double** array = (double**) malloc(sizeof(double*)*linhas);

    for(int i = 0; i < linhas; i++) {
        array[i] = (double*) malloc(sizeof(double)*5);
        if(fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\n",&array[i][0],&array[i][1],&array[i][2],&array[i][3],&array[i][4]));
    }
    fclose(file);
    return array;
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

// Função para mesclar dois subvetores de arr[] e reordenar o vetor brr[] de acordo
void merge(int *arr, igraph_vector_t *brr, int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    // Cria vetores temporários
    int *L = (int *)malloc(n1 * sizeof(int));
    int *R = (int *)malloc(n2 * sizeof(int));
    int *bL = (int *)malloc(n1 * sizeof(int));
    int *bR = (int *)malloc(n2 * sizeof(int));

    // Copia os dados para os vetores temporários L[], R[], bL[] e bR[]
    for (i = 0; i < n1; i++) {
        L[i] = arr[l + i];
        bL[i] = VECTOR(*brr)[l + i];
    }
    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1 + j];
        bR[j] = VECTOR(*brr)[m + 1 + j];
    }

    // Mescla os vetores temporários de volta para arr[l..r] e brr[l..r]
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (L[i] >= R[j]) {
            arr[k] = L[i];
            VECTOR(*brr)[k] = bL[i];
            i++;
        } else {
            arr[k] = R[j];
            VECTOR(*brr)[k] = bR[j];
            j++;
        }
        k++;
    }

    // Copia os elementos restantes de L[] e bL[], se houver
    while (i < n1) {
        arr[k] = L[i];
        VECTOR(*brr)[k] = bL[i];
        i++;
        k++;
    }

    // Copia os elementos restantes de R[] e bR[], se houver
    while (j < n2) {
        arr[k] = R[j];
        VECTOR(*brr)[k] = bR[j];
        j++;
        k++;
    }

    free(L);
    free(R);
    free(bL);
    free(bR);
}

// Função principal que implementa o MergeSort
void mergeSort(int *arr, igraph_vector_t* brr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;

        mergeSort(arr, brr, l, m);
        mergeSort(arr, brr, m + 1, r);

        merge(arr, brr, l, m, r);
    }
}
void print_side_by_side(int* a,int*b,int N){
    for (int  i = 0; i < N; i++) printf("%d %d %d\n",a[i],b[i],i);
    
}