
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"
#include "CM.h"
#include "SBM.h"

struct Graph{
    int **viz;
    int Nodes;
    int edges;
};
struct LCM{
    int existir;
    int **n_existir;
    int **degree;
    int *faixa;
    struct Graph G;
};

double** double_get_array(int N,char* ptr){
    FILE* file;
    file = fopen(ptr,"r");
    double** degree = (double**) malloc(sizeof(double*)*N);

    for(int i = 0; i < N; i++) {
        degree[i] = (double*) malloc(sizeof(double)*5);
        if(fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

int** get_array(int N,char* ptr){
    FILE* file;
    file = fopen(ptr,"r");
    int** degree = (int**) malloc(sizeof(int*)*N);

    for(int i = 0; i < N; i++) {
        degree[i] = (int*) malloc(sizeof(int)*5);
        if(fscanf(file,"%d\t%d\t%d\t%d\t%d\n",&degree[i][0],&degree[i][1],&degree[i][2],&degree[i][3],&degree[i][4]));
    }
    fclose(file);
    return degree;
}

int somatorio(int** array,int elemento,int N){
    int soma = 0;
    for (int i = 0; i < N; i++) soma += array[elemento][i];
    return soma;
}

int** M_bubbleSort(int **matriz, int linhas, int colunas) {
    int i, j;
    for (i = 0; i < linhas - 1; i++) {
        for (j = 0; j < linhas - i - 1; j++) {
            int somaAtual = 0, somaProxima = 0;
            int *linhaAtual = matriz[j];
            int *linhaProxima = matriz[j + 1];

            for (int k = 0; k < colunas; k++) {
                somaAtual += linhaAtual[k];
                somaProxima += linhaProxima[k];
            }

            if (somaAtual < somaProxima) {
                int *temp = matriz[j];
                matriz[j] = matriz[j + 1];
                matriz[j + 1] = temp;
            }
        }
    }
    return matriz;
}


/* struct LCM local_conf_model_p(struct LCM Z,int ego,double p){
    int N = 0;
    int *vizinhos = (int*) malloc(0*sizeof(int));

    for (int i = 0; i < G.viz[ego][0]; i++){
        int vizinho = Z.G.viz[ego][i+1];
        if(somatorio(Z.degree,vizinho,5) !=0){
            vizinhos = (int*) realloc(vizinhos,(N+1)*sizeof(int));
            vizinhos[N] = vizinho;
            N++;
        }
    }
    if(N > 1){

        int *shuff = (int*) malloc(N*sizeof(int));

        for (int i = 0; i < N; i++){
            int site = vizinhos[i];

            for (int j = 0; j < N; j++) shuff[j] = vizinhos[j];
            shuff = randomize(shuff,N,i+site);
            shuff = ending(shuff,N,site,0);
            for (int j = 0; j < N-1; j++){
                if(somatorio(Z.degree,site,N) == 0) break;

                int vizinho = shuff[j];
                double rand = genrand64_real2();
                if(rand<=p){

                    int rep = check_existence(Z.mat,Z.G.edges,site,vizinho);
                    if((rep == 1) || (Z.degree[vizinho][Z.faixa[site]] == 0) || (Z.degree[site][Z.faixa[vizinho]] == 0) ) continue;
                    rep = check_existence(Z.n_existir,Z.existir,site,vizinho);
                    if(rep == 1) continue;
                    int faixa1 = Z.faixa[site];
                    int faixa2 = Z.faixa[vizinho];

                    Z.mat = (int**) realloc(Z.mat,(Z.G.edges+1)*sizeof(int*));
                    Z.mat[Z.G.edges] = (int*) malloc(2* sizeof(int));
                    Z.mat[Z.G.edges][0] = site;
                    Z.mat[Z.G.edges][1] = vizinho;
                    Z.G.edges += 1;

                    Z.degree[site][faixa2]--;
                    Z.degree[vizinho][faixa1]--;

                    Z.G.viz[site][0]++;
                    Z.G.viz[vizinho][0]++;
                    Z.G.viz[site] = (int*) realloc(Z.G.viz[site],(Z.G.viz[site][0]+1)*sizeof(int));
                    Z.G.viz[vizinho] = (int*) realloc(Z.G.viz[vizinho],(Z.G.viz[vizinho][0]+1)*sizeof(int));
                    Z.G.viz[site][Z.G.viz[site][0]] = vizinho;
                    Z.G.viz[vizinho][Z.G.viz[vizinho][0]] = site;
                }

                else{
                    int rep = check_existence(Z.n_existir,Z.existir,site,vizinho);
                    if(rep == 0){
                        Z.n_existir = (int**) realloc(Z.n_existir,(Z.existir+1)*sizeof(int*));
                        Z.n_existir[Z.existir] = (int*) malloc(2* sizeof(int));
                        Z.n_existir[Z.existir][0] = site;
                        Z.n_existir[Z.existir][1] = vizinho;
                        Z.existir += 1;
                    }
                }
            }
        }
        free(shuff);
    }
    free(vizinhos);
    return Z;
} */

struct Graph local_add_edge(struct Graph G,int** degree,int* faixa,int* shuff, int n,int site,double p){
    for (int i = 0; i < n; i++){
        int vizinho = shuff[i];
        if(somatorio(degree,site,5) == 0){
            //if(p>0) Z = local_conf_model_p(Z,site,p);
            break;
        }
        int rep = 0;
        for (int j = 0; j < G.viz[site][0]; j++) if(vizinho == G.viz[site][j+1]) {rep = 1; break;}
        

        int faixa1 = faixa[site];
        int faixa2 = faixa[vizinho];
        if((rep == 1) || (degree[vizinho][faixa1] == 0) || (degree[site][faixa2] == 0) ) continue;
        G.edges += 1;

        degree[site][faixa2]--;
        degree[vizinho][faixa1]--;
        
        G.viz[site][0]++;
        G.viz[vizinho][0]++;
        G.viz[site] = (int*) realloc(G.viz[site],(G.viz[site][0]+1)*sizeof(int));
        G.viz[vizinho] = (int*) realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
        G.viz[site][G.viz[site][0]] = vizinho;
        G.viz[vizinho][G.viz[vizinho][0]] = site;
    }
    return G;
}

struct Graph local_configuration_model(int N, double p,int seed,int rede){
    struct Graph G;
    G.Nodes = N;
    G.viz = (int **)malloc(N*sizeof(int*));
    for (int j = 0; j < N; j++){ 
        G.viz[j] =(int*) malloc(1*sizeof(int));
        G.viz[j][0] = 0;
    }
    G.edges = 0;
    //existir = 0;
    //if(p > 0)n_existir = (int **)malloc(0* sizeof(int*));
    int** degree;
    int* faixa;
    if(rede == 0){
        char file_name[200] = "./output/teste.txt";
        degree = get_array(N,file_name);
        char sfile_name[200] = "./output/teste0.txt";
        faixa = get_faixas(N,sfile_name);
    }
    if(rede == 1){
        char file_name[200] = "./dados/geometry.txt";
        double** probability = double_get_array(5,file_name);
        char file_names[200] = "./dados/constant.txt";
        double** constant = double_get_array(5,file_names);
        degree = (int**) malloc(N*sizeof(int*));
        double* p_faixas = malloc(5*sizeof(double));
        int i = 0,j = 0;
        double r;

        p_faixas[0] = 0.28215076;
        p_faixas[1] = 0.11551187;
        p_faixas[2] = 0.32406636;
        p_faixas[3] = 0.2107681;
        p_faixas[4] = 0.06750292;

        for (i = 1; i < 5; i++) p_faixas[i] += p_faixas[i-1];

        faixa = calloc(N,sizeof(int));

        for (i = 0; i < N; i++){
            r = genrand64_real1();
            for (j = 0; j < 5; j++) if(r < p_faixas[j]) break;
            faixa[i] = j;
        }

        free(p_faixas);
        double x = 0;
        for (i = 0; i < N; i++){
            degree[i] = (int*) calloc(5,sizeof(int));
            //while(somatorio(degree,i,5) == 0) for (j = 0; j < 5; j++) degree[i][j] = geometry(probability[faixa[i]][j]);
            while(somatorio(degree,i,5) == 0) for (j = 0; j < 5; j++) degree[i][j] = potentia(probability[faixa[i]][j],constant[faixa[i]][j]);
            x += somatorio(degree,i,5);
        }
        degree = M_bubbleSort(degree,N,5);
        printf("%f\n",x/N);
    }
    else{
        char file_name[200] = "./dados/lambda.txt";
        double* probability = load_file(file_name,5);
        char file_names[200] = "./dados/multi_probability.txt";
        double** constant = double_get_array(5,file_names);
        //char file_namess[200] = "./dados/multi_constant.txt";
        //double* A = load_file(file_namess,5);

        degree = (int**) malloc(N*sizeof(int*));
        double* p_faixas = malloc(5*sizeof(double));
        int i = 0,j = 0;
        double r;

        p_faixas[0] = 0.28215076;
        p_faixas[1] = 0.11551187;
        p_faixas[2] = 0.32406636;
        p_faixas[3] = 0.2107681;
        p_faixas[4] = 0.06750292;
        printf("%f\n",constant[0][0]);
        for ( j = 0; j < 5; j++) for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
        for (i = 1; i < 5; i++)p_faixas[i] += p_faixas[i-1];
        

        faixa = calloc(N,sizeof(int));

        for (i = 0; i < N; i++){
            r = genrand64_real1();
            for (j = 0; j < 5; j++) if(r < p_faixas[j]) break;
            faixa[i] = j;
        }
        free(p_faixas);
        int k =0;
        double x = 0;
        for (i = 0; i < N; i++){
            k = 0;
            init_genrand64(seed+i);
            degree[i] = (int*) calloc(5,sizeof(int));
            while(k == 0)k = geometry(probability[faixa[i]]);
            //while(k == 0)k = generalized_geometry(probability[faixa[i]],A[faixa[i]]);
            generate_multinomial(k, 5, constant[faixa[i]], degree[i]);
            x += k;
        }
        printf("%f\n",x/N);
        degree = M_bubbleSort(degree,N,5);
    }
    //print_vetor(faixa,G.Nodes);
    int *shuff = (int*) malloc(N*sizeof(int));
    int i = 0,j = 0;
    int total = 0;
    for (i = 0; i < G.Nodes; i++)  total += somatorio(degree,i,5);
    for (i = 0; i < N; i++){
        
        for (j = 0; j < N; j++) shuff[j] = j;
        swap(&shuff[i],&shuff[N-1]);
        shuff = randomize(shuff,N-1,seed+i);

        G = local_add_edge(G,degree,faixa,shuff,N-1,i,p);
    }
    //print_matrix(degree,N,5);
    if(seed != 0){
        FILE *arquivo;
        arquivo = fopen("./dados/faixas_finais.txt","w");
        for (int i = 0; i < G.Nodes; i++) fprintf(arquivo,"%d\n",faixa[i]);
        fclose(arquivo);
    }
    double a = 0;
    int x = 0;
    for (i = 0; i < G.Nodes; i++)  x += somatorio(degree,i,5);
    for (i = 0; i < G.Nodes; i++) free(degree[i]);
    for (i = 0; i < G.Nodes; i++) {
        degree[i] = (int*) calloc(5,sizeof(int));
        for (j = 0; j < G.viz[i][0]; j++) degree[i][faixa[G.viz[i][j+1]]] += 1;
        a += G.viz[i][0];
    }
    printf("%f %f\n",a/N,(double)x/total);
    char file[200] = "./dados/graus_finais.txt";
    generate_file(file,degree,G.Nodes,5,sizeof(degree[0][0]));
    

    //if(p > 0) for (int i = 0; i < existir; i++) free(n_existir[i]);
    
    //if(p > 0)free(n_existir);
    free(degree);
    free(faixa);
    free(shuff);
    return G;
}

void generate_local_configuration_model(double p, int T){

    //char file[200] = "./output/teste.txt";
    int N = 2029;//size_txt(file)-1;
    
    double** resultados = (double **)malloc(T*sizeof(double*));
    for (int i = 0; i < T; i++){

        if(T!=1) printf("\e[1;1H\e[2J");

        struct Graph G;
        //G = local_configuration_model(N,p,i,1,0);

        resultados[i] = (double*) malloc(7*sizeof(double));
        result(G,resultados[i]);
        
        printf("%d\n",i+1);
        //if(T == 1) create_network(G,p);
        
        for(int j = 0; j < N; j++)free(G.viz[j]);
        
        free(G.viz);
    }
    generate_resultados(resultados,T,"LCM");

}