#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <stdbool.h>
#include <igraph.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "LCM.h"
#include "SBM.h"
#include "define.h"

double* sintomatico;
double* hospitalizacao;
double* morte;
int E0 = 10;

void constant_init(){
    sintomatico = (double*) malloc(5*sizeof(double));
    hospitalizacao = (double*) malloc(5*sizeof(double));
    morte = (double*) malloc(5*sizeof(double));
    sintomatico[0] = 29.1/100;
    sintomatico[1] = 37.4/100;
    sintomatico[2] = 41.68/100;
    sintomatico[3] = 39.4/100;
    sintomatico[4] = 31.3/100; 

    hospitalizacao[0] = (double)1.04/100;
    hospitalizacao[1] = (double)1.33/100;
    hospitalizacao[2] = (double)1.38/100;
    hospitalizacao[3] = (double)7.6/100;
    hospitalizacao[4] = (double)24/100;

    morte[0] = (double)0.07710809281267686;
    morte[1] = (double)0.09561476547321676;
    morte[2] = (double)0.13779231777947765;
    morte[3] = (double)0.30044562859568047;
    morte[4] = (double)0.530315172817809;
}

bool fileExists(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

int find_root(int root,int*cluster,int N){
    if(cluster[root] >= 0){
        cluster[root] = find_root(cluster[root],cluster,N);
        root = cluster[root];
    }
    return root;
}

int maior_cluster_infectados(igraph_t* Grafo,int N,char* estagio){
    int* cluster = malloc(N*sizeof(int));
    int maior_cluster = 0;
    int i,j,root_site,vizinho,root_vizinho;
    for (i = 0; i < N; i++)cluster[i] = -1;
    for (i = 0; i < N; i++){
        igraph_vector_int_t vizinhos;
        igraph_vector_int_init(&vizinhos, 0);
        igraph_neighbors(Grafo, &vizinhos, i,IGRAPH_ALL);
        if((estagio[i] == 'E') || (estagio[i] == 'I')|| (estagio[i] == 'A')){ 
            
            root_site = find_root(i,cluster,N);
            for (j = 0; j < igraph_vector_int_size(&vizinhos); j++){
                vizinho = VECTOR(vizinhos)[j];
                if((estagio[vizinho] == 'E') || (estagio[vizinho] == 'I')|| (estagio[vizinho] == 'A')){
                    
                    root_vizinho = find_root(vizinho,cluster,N);
                    if(root_site != root_vizinho){
                        cluster[vizinho] = root_site;
                        cluster[root_site] += -1;
                        if(-cluster[root_site] > maior_cluster) maior_cluster = -cluster[root_site];
                    }
                }
            }
        }
        igraph_vector_int_destroy(&vizinhos);
    }
    free(cluster);
    return maior_cluster;
}

double calc_estagio(int site,char* estagio,double* prob_estagio, igraph_t* Grafo,bool* vacinado,double* avg,bool weight){

    switch (estagio[site]){
        case 'S':// Suscetível
            double beta = 0;

            igraph_vector_int_t vizinhos;
            igraph_vector_int_init(&vizinhos, 0);
            igraph_neighbors(Grafo, &vizinhos, site,IGRAPH_ALL);
            igraph_integer_t eid;
            double peso;
            for (int j = 0; j < igraph_vector_int_size(&vizinhos); j++){
                
                int vizinho = VECTOR(vizinhos)[j];
                peso = 1;
                if (weight){
                    igraph_get_eid(Grafo, &eid, site, vizinho, IGRAPH_UNDIRECTED, 0);
                    peso = (double) igraph_cattribute_EAN(Grafo, "duracao", eid);
                    
                }
                
                if(estagio[vizinho] == 'E'){
                    if(!vacinado[vizinho]) beta += beta1*peso;
                    else beta += 0.058*beta1*peso; // 0.058
                }
                if(estagio[vizinho] == 'A'){
                    if(!vacinado[vizinho]) beta += beta2*peso;
                    else beta += 0.058*beta2*peso;// 0.058
                }
            }
            igraph_vector_int_destroy(&vizinhos);
            if(!vacinado[site]) return beta/ *avg;
            else return 0.058*beta/ *avg;// 0.058
            break;
        case 'E': // Exposto
            return sigma;
            break;
        case 'A': //Assintomático
            return gamma_A;
            break;
        case 'I': // Sintomático
            if(prob_estagio[site]< hospitalizacao[(int) VAN(Grafo, "faixa", site)]*(vacinado[site] ? 0.034 : 1)) return phi;
            else return gamma_I;
            break;
        case 'H': // Hospitalziado
            if(prob_estagio[site]<morte[(int) VAN(Grafo, "faixa", site)]*(vacinado[site] ? 0.034 : 1)) return delta;
            else return gamma_H;
            break;
        case 'R':
            return recupera;
            break;
        default:
            return 0;
            break;
    }
}

int find_lista(int site,int* lista){
    if(lista[site] != site) return find_lista(lista[site],lista);
    return site;
}
double* atualiza_estagio(double* rates,int site,double* rate,char* estagio,double* prob_estagio, igraph_t* Grafo,bool* vacinado,double* avg,bool weight){
    int i;
    *rate -= rates[site];
    
    rates[site] = calc_estagio(site, estagio,prob_estagio,Grafo,vacinado,avg,weight);
    *rate += rates[site];
    
    igraph_vector_int_t vizinhos;
    igraph_vector_int_init(&vizinhos, 0);
    igraph_neighbors(Grafo, &vizinhos, site,IGRAPH_ALL);
    for (i = 0; i < igraph_vector_int_size(&vizinhos); i++){
        int vizinho = VECTOR(vizinhos)[i];
        *rate -= rates[vizinho];
        rates[vizinho] = calc_estagio(vizinho, estagio,prob_estagio,Grafo,vacinado,avg,weight);
        *rate += rates[vizinho];
    }
    igraph_vector_int_destroy(&vizinhos);
    return rates;
}

double infect_init(char* estagio, double* rates,double* prob_estagio,bool* vacinado,igraph_t* Grafo,bool weight,double* avg_degree){

    int N = igraph_vcount(Grafo);
    
    int i = 0,site;
    double rate = 0,r;
    double EXPOSTOS = E0;
    for (i = 0; i < N; i++){
        prob_estagio[i] = genrand64_real1();
        estagio[i] = 'S';
    }
    while(EXPOSTOS!=0){
        r = genrand64_real1();
        site = r*N;
        if(estagio[site] == 'S'){
            estagio[site] = 'E';
            EXPOSTOS--;
            rates = atualiza_estagio(rates,site,&rate,estagio,prob_estagio,Grafo,vacinado,avg_degree,weight);
        }
    }
    return rate;
}

void infect(igraph_t* Grafo,double* rate, char* estagio,double* rates,double* prob_estagio,bool* vacinado,double** infect_time,const bool weight,const int primeiro_dia,const int ultimo_dia,int N_vacinas,double* estado,double* avg_degree){
    int N = igraph_vcount(Grafo);

    double Suscetiveis = estado[0];
    double Expostos = estado[1];
    double Assintomaticos = estado[2];
    double Sintomaticos = estado[3];
    double Hospitalizados = estado[4];
    double Recuperados = estado[5];
    double Mortos = estado[6];
    double Infectados = estado[7];
    double t = primeiro_dia;
    int i,site,s = 2*primeiro_dia;

    double tempo = 0;
    double ano = primeiro_dia;
    double Delta;
    double new_rate;
    while (ano < ultimo_dia){
        if(rate==0) break;
        tempo = exponentialRand(*rate);
        if(tempo ==0) break;
        Delta = *rate*genrand64_real1();
        new_rate = 0;
        ano += tempo;
        if(ano >= ultimo_dia) break;
        for (i = 0; i < N; i++){
            new_rate += rates[i];
            if(new_rate>=Delta) break;
        }
        site = i;
        switch (estagio[i]){
            case 'S':// Sucetível - 0
                estagio[i] = 'E';
                Expostos++;
                Suscetiveis--;
                break;
            case 'E': // Exposto - 1
                if(prob_estagio[i]<sintomatico[(int) VAN(Grafo, "faixa", i)]*(vacinado[i] ? 0.34693877551 : 1)){
                    estagio[i] = 'I';
                    Sintomaticos++;
                    Infectados++;
                }
                else{
                    estagio[i] = 'A';
                    Infectados++;
                    Assintomaticos++;
                    if((N_vacinas!=0) && (ano >= dia_infecao+1 ) && (!vacinado[i])){
                        vacinado[i] = true;
                        N_vacinas--;
                    }
                    
                }
                Expostos--;
                break;
            case 'A': //Assintomático
                estagio[i] = 'R';
                Assintomaticos--;
                Recuperados++;
                if((N_vacinas!=0) && (ano >= dia_infecao+1 ) && (!vacinado[i])){
                    vacinado[i] = true;
                    N_vacinas--;
                }
                break;
            case 'I': // Sintomático
                if(prob_estagio[i]<hospitalizacao[(int) VAN(Grafo, "faixa", site)]){
                    estagio[i] = 'H';
                    Hospitalizados++;
                }
                else{
                    estagio[i] = 'R';
                    Recuperados++;
                    if((N_vacinas!=0) && (ano >= dia_infecao+1 ) && (!vacinado[i])){
                        vacinado[i] = true;
                        N_vacinas--;
                    }
                }
                Sintomaticos--;
                break;
            case 'H': // Hospitalziado - 4
                if(prob_estagio[i]<morte[(int) VAN(Grafo, "faixa", site)]){
                    estagio[i] = 'D';
                    Mortos++;
                }
                else{
                    estagio[i] = 'R';
                    Recuperados++;
                    if((N_vacinas!=0) && (ano >= dia_infecao+1 ) && (!vacinado[i])){
                        vacinado[i] = true;
                        N_vacinas--;
                    }
                }
                Hospitalizados--;
                break;
            case 'R': // Recuperados 
                estagio[i] = 'S';
                Suscetiveis++;
                Recuperados--;
                if((N_vacinas!=0) && (ano >= dia_infecao+1) && (!vacinado[i])){
                    vacinado[i] = true;
                    N_vacinas--;
                }
                break;
            default:
                break;
        }
        prob_estagio[i] = genrand64_real1();
        if(ano > t){
            s++;
            t += 0.5;
            infect_time[s - 1][0] += Suscetiveis/N;
            infect_time[s - 1][1] += Expostos/N;
            infect_time[s - 1][2] += Assintomaticos/N;
            infect_time[s - 1][3] += Sintomaticos/N;
            infect_time[s - 1][4] += Hospitalizados/N;
            infect_time[s - 1][5] += Recuperados/N;
            infect_time[s - 1][6] += Mortos/N;
            infect_time[s - 1][7] += Infectados/N;
            //infect_time[s - 1][9] += (double)maior_cluster_infectados(Grafo,N,estagio)/N;
            //printf("%f %f %f %f %f %f %f %f %d %f %f\n",ano,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos,s,*rate,t);
            
        }
        //printf("%f %f %f %f %f %f %f %f %d %f\n",ano,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos,s,rate);
        if(t >= dias+1) break;
        if((Expostos == Assintomaticos) && (Assintomaticos == Sintomaticos) && (Sintomaticos ==Hospitalizados) && (Hospitalizados == 0)) break;

        rates = atualiza_estagio(rates,site,rate,estagio,prob_estagio,Grafo,vacinado,avg_degree,weight);
    }
    estado[0] = Suscetiveis;
    estado[1] = Expostos;
    estado[2] = Assintomaticos;
    estado[3] = Sintomaticos;
    estado[4] = Hospitalizados;
    estado[5] = Recuperados;
    estado[6] = Mortos;
    estado[7] = Infectados;
    if(ano <= dias){
        for (i = s; i < dias*2; i++){
            infect_time[i][0] += (double) Suscetiveis/N;
            infect_time[i][1] += (double) Expostos/N;
            infect_time[i][2] += (double) Assintomaticos/N;
            infect_time[i][3] += (double) Sintomaticos/N;
            infect_time[i][4] += (double) Hospitalizados/N;
            infect_time[i][5] += (double) Recuperados/N;
            infect_time[i][6] += (double) Mortos/N;
            infect_time[i][7] += (double) Infectados/N;
            //infect_time[i][8] += (double)maior_cluster_infectados(Grafo,N,estagio)/N;
        }
    }    
}

void save_file(double*** infect_time,const int* redes,const bool* weight,const double* N,const double* p,const int* estrategy){
    int i,j,k;
    char arquivo[800];
    char arquivo2[800];
    char *file_vacina[] = {"idade", "grau", "close", "harmonic","betwenness","eigenvector","eccentricity","clustering","kshell","random","pagerank","graumorte","probhosp","probmorte","probhospassin","probmortepassin","wclose","wharmonic","wbetwenness","weigenvector","wpagerank"};
    for (i = 0; i < infecao_total; i++)
        for (j = 0; j < q_resultados; j++) 
            infect_time[0][i][j] /= *redes;
    
    if(*weight)sprintf(arquivo,"./output/time/%d/ponderado/p/infect_%.2f.txt",(int)*N,*p);
    else sprintf(arquivo,"./output/time/%d/nponderado/p/infect_%.2f.txt",(int)*N,*p);

    if(!fileExists(arquivo)) generate_file(arquivo,infect_time[0],infecao_total,q_resultados,sizeof(double),0);
    double f;
    if(*weight) sprintf(arquivo2,"./output/vacina/%d/ponderado/%s_%.2f.txt",(int)*N,file_vacina[*estrategy],*p);
    else sprintf(arquivo2,"./output/vacina/%d/nponderado/%s_%.2f.txt",(int)*N,file_vacina[*estrategy],*p);
    FILE* file  = fopen(arquivo2,"w");;

    for (k = 0; k < 100; k++){

        f = (double) (k+1)/100;

        for (i = infecao_total+1; i < tempo_total; i++) for (j = 0; j < q_resultados; j++) infect_time[k][i][j] /= *redes;

        if((f == 0.75) ||(f == 0.25) ||(f == 0.50) || (f == 1.0)){
            if(*redes > cut_rede){
                if(*weight) sprintf(arquivo,"./output/time/%d/ponderado/%s_%.2f_%.2f.txt",(int)*N,file_vacina[*estrategy],*p,f);
                else sprintf(arquivo,"./output/time/%d/nponderado/%s_%.2f_%.2f.txt",(int)*N,file_vacina[*estrategy],*p,f);
                generate_file(arquivo,infect_time[k],tempo_total,q_resultados,sizeof(double),infecao_total+1);
            }
        }
        
        fprintf(file,"%f %f %f %f\n",f, infect_time[k][tempo_total-1][6], infect_time[k][tempo_total-1][4], infect_time[k][tempo_total-1][7]);
        
    }
    fclose(file);
    for (k = 0; k < 100; k++){
        for (i = 0; i < tempo_total; i++) free(infect_time[k][i]);
        free(infect_time[k]);
    }
    free(infect_time);
    printf("Terminou: %s %f",file_vacina[*estrategy],*p);
}

void generate_infect(double N,double p,int seed,int redes,int estrategy,bool weight){
    
    uint32_t i,j;
    create_folder((int)N);
    
    constant_init();
    double*** infect_time = (double***) malloc(100*sizeof(double**));
    for (i = 0; i < 100; i++){
        infect_time[i] = (double**) malloc(tempo_total*sizeof(double*));
        for (j = 0; j < tempo_total; j++)infect_time[i][j] = (double*) calloc(q_resultados,sizeof(double));        
    }
    int rede;
    omp_set_num_threads(THREADS);
    #pragma omp parallel for
    for (rede = 0; rede < redes; rede++){

        double avg_degree;
        igraph_vector_int_t centralidade;
        init_genrand64(rede*estrategy*(weight+1));

        igraph_t Grafo = local_configuration_model( N, p,rede,weight,&avg_degree,&centralidade,estrategy,true);
        
        double* prob_estagio = (double*) malloc(N*sizeof(double));
        double* rates = (double*) calloc(N,sizeof(double));
        double* estado = (double*) calloc(8,sizeof(double));
        char* estagio = malloc(N*sizeof(char));
        bool* vacinado =(bool*) calloc(N,sizeof(bool));

        double* copia_prob_estagio = (double*) malloc(N*sizeof(double));
        double* copia_estado = (double*) calloc(8,sizeof(double));
        double* copia_rates = (double*) calloc(N,sizeof(double));
        char* copia_estagio = malloc(N*sizeof(char));
        
        double rate = infect_init(estagio, rates,prob_estagio,vacinado,&Grafo,weight,&avg_degree);
        double copia_rate;
        estado[0] = N-10;
        estado[1] = 10;
        infect(&Grafo,&rate,estagio, rates,prob_estagio,vacinado,infect_time[0],weight,0,dia_infecao,0,estado,&avg_degree);
        int f,N_vacinas,sitio_vacinado;
        
        for (f = 1; f <= 100; f++){
            copia_rate = rate;
            N_vacinas = (int)f*N/100;
            for (i = 0; i < 8; i++)copia_estado[i] = estado[i];
            for (i = 0; i < N; i++){

                copia_prob_estagio[i] = prob_estagio[i];
                copia_rates[i] = rates[i];
                copia_estagio[i] = estagio[i];
                
                sitio_vacinado = VECTOR(centralidade)[i];

                if((estagio[sitio_vacinado] == 'S') || (estagio[sitio_vacinado] == 'R') || (estagio[sitio_vacinado] == 'A')){
                    if(N_vacinas != 0){
                        vacinado[sitio_vacinado] = true;
                        N_vacinas--;
                    }
                }
            }
            //print_vetor(copia_estado,8,sizeof(double));
            infect(&Grafo,&copia_rate,copia_estagio, copia_rates,copia_prob_estagio,vacinado,infect_time[f-1],weight,dia_infecao,dias,N_vacinas,copia_estado,&avg_degree);
            //if(f%10 == 0)printf("%d %.2f\n",rede+1,(double) f/100);
            free(vacinado);
            vacinado =(bool*) calloc(N,sizeof(bool));
        }
        

        
        igraph_vector_int_destroy(&centralidade);
        free(rates);
        free(prob_estagio);
        free(estagio);
        free(vacinado);
        free(estado);
        free(copia_rates);
        free(copia_prob_estagio);
        free(copia_estagio);
        free(copia_estado);
        igraph_destroy(&Grafo);
    }
    if(redes > cut_rede) save_file(infect_time,&redes,&weight,&N,&p,&estrategy);
    else{
        int k;
        for (i = 0; i < infecao_total; i++)
            for (j = 0; j < q_resultados; j++) 
                infect_time[0][i][j] /= redes;
        print_matrix((void**)infect_time[0],infecao_total,sizeof(double),sizeof(double));
        for (k = 0; k < 100; k++){

            for (i = infecao_total+1; i < tempo_total; i++) for (j = 0; j < q_resultados; j++) infect_time[k][i][j] /= redes;
            print_matrix((void**)infect_time[k],tempo_total,sizeof(double),sizeof(double));
        }
    }

    
}