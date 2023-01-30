#include <stdlib.h>
#include <stdio.h>
#include "mtwister.h"
#include <string.h>
#include <math.h>

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
    double cov = 0;

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
int* ending(int* array, int n, int site){
    for (int i = 0; i < n; i++){
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
        int j = rand() % (i+1);
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

int size_txt(){
    FILE* f;
    int L = 0;
    char c;
    f = fopen("./degree.txt","r");
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