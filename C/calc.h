#ifndef CALC_H
#define CALC_H

int fatorial(int n,int p);
int* get_faixas(int N);
int check_existence(int** outro, int N,int site, int vizinho);
int** get_degree2(int N);
double correlation(double* x,double* y,int N);
int combination(int n, int p);
void swap (int *a, int *b);
int* randomize (int* array, int n,int seed);
int* ending(int* array, int n, int site,int inicio);
int* bubble_sort_by(int *array,int *valor,int N);
void print_vetor(int* array,int N);
int* bubble_sort(int*vetor,int N);
int* remove_(int *vetor,int elemento,int N);
int size_txt();
void print_matrix(int** mat,int N,int n);
void generate_resultados(double** resultados, int T,char arquivo[]);

#endif