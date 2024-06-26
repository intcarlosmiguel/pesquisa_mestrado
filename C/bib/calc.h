#ifndef CALC_H
#define CALC_H
#include <igraph.h>
int fatorial(int n,int p);
int* get_faixas(int N,char* str);
int check_existence(int** outro, int N,int site, int vizinho);
int** get_degree2(int N);
double correlation(double* x,double* y,int N);
int combination(int n, int p);
void swap (int *a, int *b);
int* randomize (int* array, int n,int seed);
int* ending(int* array, int n, int site,int inicio);
int* bubble_sort_by(int *array,int *valor,int N);
void print_vetor(void* array,int N,int check);
int* bubble_sort(int*vetor,int N);
int* remove_(int *vetor,int elemento,int N);
int size_txt(char *str);
void print_matrix(void** mat,int N,int n,int check);
void generate_resultados(double** resultados, int T,char arquivo[]);
double normalRand(double mean, double stdDev);
double exponentialRand(double lambda);
void quicksort(int *arr, int low, int high);
void sortIntByRef(int *arr, void *ref, int n,int check);
void quickSort(int *arr, int low, int high);
int geometry(double lambda);
int potentia(double lambda,double A);
void generate_file(char* filename,void* array,int linhas,int colunas,int check,int inicio,char* mode);
void load_file(char* filename,void* array,int check);
void generate_multinomial(int n, int k, double *probabilities, int *outcomes);
int generalized_geometry(double lambda,double A);
int empiric_distribution(double* distribution);
void bubbleSort_by(int* v, int* v2, int n);
double** double_get_array(char* ptr);
int* arange(int inicio, int fim, int passo);
void mergeSort(int *arr, igraph_vector_t* brr, int l, int r);

void print_side_by_side(int* a,int*b,int N);
void create_folder(int N);

#endif