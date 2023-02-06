#ifndef CALC_H
#define CALC_H

int fatorial(int n,int p);
int* get_faixas(int N);
int** get_degree(int N);
double correlation(double* x,double* y,int N);
int combination(int n, int p);
void swap (int *a, int *b);
int* randomize (int* array, int n,int seed);
int* ending(int* array, int n, int site);
int* bubble_sort_by(int *array,int *valor,int N);
void print_vetor(int* array,int N);
int* bubble_sort(int*vetor,int N);
int* remove_(int *vetor,int elemento,int N);
int size_txt();
void print_matrix(int** mat,int N,int n);

#endif