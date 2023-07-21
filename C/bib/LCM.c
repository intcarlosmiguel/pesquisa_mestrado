
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
struct Retorno{
    int* faixa;
    int** degree;
};

double distance(double* point1,double* point2,int N){
    double distance2 = 0;
    for (int i = 0; i < N; i++) distance2 += pow(point1[i] - point2[i],2);
    distance2 = sqrt(distance2);
    return distance2;
}

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

int vetor_soma(int* array, int L){
    int soma = 0;
    for (int i = 0; i < L; i++) soma += array[i];
    return soma;
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

double** distribution_faixas(){
    double** distribution = (double**) malloc(5*sizeof(double*));
    char filename[200];
    FILE *file;
    int N;
    for (int i = 0; i < 5; i++){
        sprintf(filename,"./dados/distribution_%d.txt",i);
        file = fopen(filename,"r");
        N = size_txt(filename);
        distribution[i] = (double*) calloc(N,sizeof(double));
        for (int j = 0; j < N; j++) if(fscanf(file, "%lf\n", &distribution[i][j]));
        fclose(file);
    }
    return distribution;
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

struct Retorno generate_degree(int N,int seed,int check){
    struct Retorno Return;/*
    if(check == 0){
        char file_name[200] = "./output/teste.txt";
        degree = get_array(N,file_name);
        char sfile_name[200] = "./output/teste0.txt";
        faixa = get_faixas(N,sfile_name);
    }
    if(check == 1){

        char file_name[200] = "./dados/geometry.txt";
        double** probability = double_get_array(5,file_name);
        char file_names[200] = "./dados/constant.txt";
        double** constant = double_get_array(5,file_names);

        degree = (int**) malloc(N*sizeof(int*));
        double* p_faixas = malloc(5*sizeof(double));
        int* grau = (int*) malloc(N*sizeof(int));
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
            //while(somatorio(degree,i,5) == 0) for (j = 0; j < 5; j++) degree[i][j] = potentia(probability[faixa[i]][j],constant[faixa[i]][j]);
            grau[i] = somatorio(degree,i,5);
            x += grau[i];
        }
        degree = M_bubbleSort(degree,N,5);
        printf("%f\n",x/N);
        bubbleSort_by(faixa, grau, N);
    }
    */if(check == 2){

        char file_name[200] = "./dados/multi_probability.txt";
        double** constant = double_get_array(5,file_name);
        double** distribution = distribution_faixas();

        int* n_faixa = calloc(5,sizeof(int));

        Return.degree = (int**) malloc(N*sizeof(int*));
        int* grau = (int*) malloc(N*sizeof(int));
        double* p_faixas = malloc(5*sizeof(double));
        int i = 0,j = 0;
        double r;

        p_faixas[0] = 0.34393182;
        p_faixas[1] =  0.14065873;
        p_faixas[2] = 0.31957894;
        p_faixas[3] = 0.15781574;
        p_faixas[4] = 0.03801476;
        for ( j = 0; j < 5; j++) for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
        for (i = 1; i < 5; i++)p_faixas[i] += p_faixas[i-1];
        

        Return.faixa = calloc(N,sizeof(int));

        for (i = 0; i < N; i++){
            j = 0;
            r = genrand64_real1();
            while(r > p_faixas[j]) j++;
            n_faixa[j]++;
        }
        //for (i = 0; i < 4; i++) n_faixa[i+1] += n_faixa[i]; 

        j = 0;
        int s = n_faixa[0];

        for (i = 0; i < N; i++){
            if(i==s) {j++;s+= n_faixa[j];}
            Return.faixa[i] = j;
        }

        int k = 0;
        double x = 0;

        for (i = 0; i < N; i++){

            k = 0;
            init_genrand64(seed+i);

            Return.degree[i] = (int*) calloc(5,sizeof(int));

            while(k == 0)k = empiric_distribution(distribution[Return.faixa[i]]);

            grau[i] = k;

            generate_multinomial(k, 5, constant[Return.faixa[i]], Return.degree[i]);

        }
        Return.degree = M_bubbleSort(Return.degree,N,5);
        bubbleSort_by(Return.faixa, grau, N);

        for(i = 0;i < 5;i++){
            //media_grau[i] /= n_faixa[i];
            free(constant[i]);
            free(distribution[i]);
        }

        //print_vetor(media_grau,5,sizeof(double));

        free(constant);
        free(distribution);
        free(p_faixas);
        //free(media_grau);
        free(n_faixa);
        free(grau);
    }
    if(check == 3){
        int i,j;
        char file_name[200] = "./dados/distribution.txt";
        int tamanho = size_txt(file_name);
        double* probability = (double*) calloc(tamanho,sizeof(double));
        load_file(file_name,probability,sizeof(double));
        
        double* p_faixas = malloc(5*sizeof(double));
        Return.degree = (int**) malloc(N*sizeof(int*));

        char file_names[200] = "./dados/multi_probability.txt";
        double** constant = double_get_array(5,file_names);
        p_faixas[0] = 0.34393182;
        p_faixas[1] =  0.14065873;
        p_faixas[2] = 0.31957894;
        p_faixas[3] = 0.15781574;
        p_faixas[4] = 0.03801476;

        for (i = 1; i < 5; i++)p_faixas[i] += p_faixas[i-1];
        Return.faixa = calloc(N,sizeof(int));
        for (i = 0; i < N; i++){
            init_genrand64(seed);
            double r = genrand64_real1();
            for (j = 0; j < 5; j++) if(r < p_faixas[j]) break;
            Return.faixa[i] = j;
            seed++;
        }
        free(p_faixas);
        double x = 0;
        for (i = 0; i < N; i++){
            int k = 0;
            init_genrand64(seed);
            Return.degree[i] = (int*) calloc(5,sizeof(int));
            //if(rede == 2) while(k == 0)k = geometry(probability[faixa[i]]);
            //while(k == 0)k = empiric_distribution(probability);
            generate_multinomial(k, 5, constant[Return.faixa[i]], Return.degree[i]);
            x += k;
            seed++;
        }
        printf("Grau médio: %f\n",x/N);
    }
    return Return;
}

void calc_metrics(struct Graph G,int* faixa){

    int i;
    int* n_faixa = calloc(5,sizeof(int));
    double* media_grau = calloc(5,sizeof(double));
    double k_medio = 0;

    for (i = 0; i < G.Nodes; i++){
        k_medio += G.viz[i][0];
        media_grau[faixa[i]] += G.viz[i][0];
        n_faixa[faixa[i]]++;
    }
    for ( i = 0; i < 5; i++) media_grau[i] /= n_faixa[i];
    printf("\n%f\n",k_medio/G.Nodes);
    print_vetor(media_grau,5,sizeof(double));
    free(n_faixa);
    free(media_grau);
}

void generate_file2(struct Graph G,int* faixa){

    int i,j;
    FILE* arquivo;
    FILE* arquivo2;
    char test[200] = "./test.txt";
    char test2[200] = "./test2.txt";

    arquivo = fopen(test,"a");
    arquivo2 = fopen(test2,"a");

    int** final = malloc(G.Nodes*sizeof(int*));
    for (i = 0; i < G.Nodes; i++){
        final[i] = calloc(5,sizeof(int));
        for (j = 0; j < G.viz[i][0]; j++) final[i][faixa[G.viz[i][j+1]]]++;
        for (j = 0; j < 5; j++){
            if(j!= 4) fprintf(arquivo,"%d\t",final[i][j]);
            else fprintf(arquivo,"%d\n",final[i][j]);
        }
        fprintf(arquivo2,"%d\n",faixa[i]);

    }
    fclose(arquivo);
    fclose(arquivo2);
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

    int** degree;
    int* faixa;
    struct Retorno Return = generate_degree(N,seed,2);
    faixa = Return.faixa;
    degree = Return.degree;
    
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
    
    //calc_metrics(G,faixa);

    int ligacoes_restantes = 0;
    int n = 0;
    int* ids = (int*) malloc(0*sizeof(int));
    for (i = 0; i < G.Nodes; i++){
        ligacoes_restantes += somatorio(degree,i,5);
        if(somatorio(degree,i,5)!= 0 ){
            ids = realloc(ids,sizeof(int)*(n+1));
            ids[n] = i;
            n++;
        }
    }
    
    char file_names[200] = "./dados/multi_probability.txt";
    double** constant = double_get_array(5,file_names);
    double** distribution = distribution_faixas();
    double* p_faixas = malloc(5*sizeof(double));
    double r;
    int cont = 0;

    p_faixas[0] = 0.34393182;
    p_faixas[1] =  0.14065873;
    p_faixas[2] = 0.31957894;
    p_faixas[3] = 0.15781574;
    p_faixas[4] = 0.03801476;

    for ( j = 0; j < 5; j++) for (i = 1; i < 5; i++) constant[j][i] += constant[j][i-1];
    for (i = 1; i < 5; i++) p_faixas[i] += p_faixas[i-1];
    
    while(ligacoes_restantes != 0){

        int k = 0;
        int faixas = 0;
        int* vetor = calloc(5,sizeof(int));

        r = genrand64_real1();
        while(r > p_faixas[faixas]) faixas++;
        while(k == 0)k = empiric_distribution(distribution[faixas]);
        
        if(ligacoes_restantes - k < 0) continue;

        int* ligacoes = (int*) calloc(k,sizeof(int));

        ids = randomize(ids,n,seed);
        seed++;
        cont = k;

        for (i = 0; i < n; i++){
            int site = ids[i];
            if(degree[site][faixas] != 0){
                vetor[faixa[site]]++;
                degree[site][faixas]--;
                ligacoes[(k - cont)] = site;
                cont--;
            }
            if(cont == 0) break;
        }
        ligacoes_restantes -= (k - cont);
        if(vetor_soma(vetor,5) != 0){
            cont = vetor_soma(vetor,5);

            G.viz = realloc(G.viz,(G.Nodes+1)*sizeof(int*));
            G.viz[G.Nodes] = calloc(cont+1,sizeof(int));
            G.viz[G.Nodes][0] = cont;

            for (i = 0; i < cont; i++) G.viz[G.Nodes][i+1] = ligacoes[i];

            for (i = 0; i < cont; i++){
                int vizinho = ligacoes[i];
                G.viz[vizinho][0]++;
                G.viz[vizinho] = realloc(G.viz[vizinho],(G.viz[vizinho][0]+1)*sizeof(int));
                G.viz[vizinho][G.viz[vizinho][0]] = G.Nodes;
            }

            faixa = realloc(faixa,(G.Nodes+1)*sizeof(int));
            faixa[G.Nodes] = faixas;
            G.edges += (k - cont);
            G.Nodes++;
        }
        free(ligacoes);
        free(vetor);
    }

    //generate_file(G,faixa);
    //calc_metrics(G,faixa);
    
    //if(p > 0) for (int i = 0; i < existir; i++) free(n_existir[i]);
    for (i = 0; i < N; i++) free(degree[i]);
    //if(p > 0)free(n_existir);
    free(degree);
    free(faixa);
    free(shuff);
    free(p_faixas);
    for(i = 0;i < 5;i++){
        free(constant[i]);
        free(distribution[i]);
    }
    free(constant);
    free(distribution);
    free(ids);
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