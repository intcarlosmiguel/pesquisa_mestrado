#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <igraph.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
//#include "rede.h"
#include "LCM.h"
#include "SBM.h"

const double beta1 = (double)0.5;
const double beta2 = (double)0.41;
const double sigma = (double)1/5.1;
//double sigma;
const double gamma_A = (double)1/7;
const double phi = (double)1/6.937854956991619;
const double gamma_I = (double)1/7;
const double gamma_H = (double)1/12.095497827980246;
const double delta = (double)1/13.681339751546528;
const double recupera = (double)1/40;

double** infect_time;
int foi = 0;
const uint16_t dias = 190;
const uint16_t dia_infecao = 100;
const uint16_t tempo_total = dias*2;
int q_resultados = 12;
int certo = 0;

bool fileExists(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}




double calc_estagio(int site,char* estagio,double* prob_estagio, igraph_t* Grafo,bool* vacinado,double* hospitalizacao,double* morte,double* avg,bool weight){


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


double* atualiza_estagio(double* rates,int site,double* rate,char* estagio,double* prob_estagio, igraph_t* Grafo,bool* vacinado,double* hospitalizacao,double* morte,double* avg,bool weight){
    int i;
    *rate -= rates[site];
    rates[site] = calc_estagio(site, estagio,prob_estagio,Grafo,vacinado,hospitalizacao,morte,avg,weight);
    *rate += rates[site];

    igraph_vector_int_t vizinhos;
    igraph_vector_int_init(&vizinhos, 0);
    igraph_neighbors(Grafo, &vizinhos, site,IGRAPH_ALL);

    for (i = 0; i < igraph_vector_int_size(&vizinhos); i++){
        int vizinho = VECTOR(vizinhos)[i];
        *rate -= rates[vizinho];
        rates[vizinho] = calc_estagio(vizinho, estagio,prob_estagio,Grafo,vacinado,hospitalizacao,morte,avg,weight);
        *rate += rates[vizinho];
    }
    igraph_vector_int_destroy(&vizinhos);
    return rates;
}

void infect(int E0,double N,double p,int* seed,double** infect_time,double f,int vacina,int redes,bool weight, int cut_rede){
    double avg_degree;
    igraph_vector_int_t centralidade;
    int N_vacinas = f*N;
    igraph_t Grafo = local_configuration_model( N, p,*seed,weight,&avg_degree,&centralidade,vacina,(N_vacinas > 0)? true: false);
    
    

    double Suscetiveis = N - E0;
    double Expostos = E0;
    double Assintomaticos = 0;
    double Sintomaticos = 0;
    double Hospitalizados = 0;
    double Recuperados = 0;
    double Mortos = 0;
    double Infectados = 0;
    bool Vac = false;
    int site;

    double* prob_estagio = (double*) malloc(N*sizeof(double));
    
    double* rates = (double*) calloc(N,sizeof(double));
    char* estagio = malloc(N*sizeof(char));
    bool* vacinado = calloc(N,sizeof(bool));
    bool* foi_vacinado = calloc(N,sizeof(bool));

    init_genrand64(*seed);

    double* sintomatico = (double*) malloc(5*sizeof(double));
    double* hospitalizacao = (double*) malloc(5*sizeof(double));
    double* morte = (double*) malloc(5*sizeof(double));

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

    int i = 0;
    double rate = 0;
    for (i = 0; i < N; i++){
        prob_estagio[i] = genrand64_real1();
        estagio[i] = 'S';
    }


    while(E0!=0){
        double r = genrand64_real1();
        site = r*N;
        if(estagio[site] == 'S'){
            estagio[site] = 'E';
            E0--;
            rates = atualiza_estagio(rates,site,&rate,estagio,prob_estagio,&Grafo,vacinado,hospitalizacao,morte,&avg_degree,weight);
        }
    }
    double t = 0.;
    int j;
    uint16_t s = 0;

    uint16_t Nhospitalizados = 0;
    double tempo = 0;
    double ano = 0;
    double Delta;
    double new_rate;
    site = -1;
    while (ano < dias){
        if((ano>= dia_infecao) && (!Vac) && (f!= 0.0)){            
            for (j = 0; j < N; j++){
                int sitio_vacinado = VECTOR(centralidade)[j];
                if((estagio[sitio_vacinado] == 'S') || (estagio[sitio_vacinado] == 'R') || (estagio[sitio_vacinado] == 'A')){
                    if(N_vacinas!=0){
                        vacinado[sitio_vacinado] = true;
                        foi_vacinado[sitio_vacinado] = true;
                        N_vacinas--;
                    }
                    if(N_vacinas == 0) break;
                }
            }
            igraph_vector_int_destroy(&centralidade);
            Vac = true;
        }
        if(site >=0){
            rates = atualiza_estagio(rates,site,&rate,estagio,prob_estagio,&Grafo,vacinado,hospitalizacao,morte,&avg_degree,weight);
        }

        if(rate==0) break;
        tempo = exponentialRand(rate);
        if(tempo ==0) break;
        Delta = rate*genrand64_real1();
        new_rate = 0;
        ano += tempo;
        if(ano >= dias) break;
        if(rate<=0) break;
        for (i = 0; i < N; i++){
            new_rate += rates[i];
            if(new_rate>=Delta) break;
        }
        site = i;
        //if(i == N) site = i-1;
        switch (estagio[i]){
            case 'S':// Sucetível - 0
                estagio[i] = 'E';
                Expostos++;
                Suscetiveis--;
                break;
            case 'E': // Exposto - 1
                if(prob_estagio[i]<sintomatico[(int) VAN(&Grafo, "faixa", i)]*(vacinado[i] ? 0.34693877551 : 1)){
                    estagio[i] = 'I';
                    Sintomaticos++;
                    Infectados++;
                }
                else{
                    estagio[i] = 'A';
                    Infectados++;
                    Assintomaticos++;
                    if((N_vacinas!=0) && (ano >= dia_infecao+1 )){
                        if(!foi_vacinado[i]){
                            vacinado[i] = true;
                            foi_vacinado[i] = true;
                            N_vacinas--;
                        }
                    }
                    
                }
                Expostos--;
                break;
            case 'A': //Assintomático
                estagio[i] = 'R';
                Assintomaticos--;
                Recuperados++;
                if((N_vacinas!=0) && (ano >= dia_infecao+1 )){
                    if(!foi_vacinado[i]){
                        vacinado[i] = true;
                        foi_vacinado[i] = true;
                        N_vacinas--;
                    }
                }
                break;
            case 'I': // Sintomático
                if(prob_estagio[i]<hospitalizacao[(int) VAN(&Grafo, "faixa", site)]){
                    estagio[i] = 'H';
                    Hospitalizados++;
                    Nhospitalizados++;
                }
                else{
                    estagio[i] = 'R';
                    Recuperados++;
                    if((N_vacinas!=0) && (ano >= dia_infecao+1 )){
                        if(!foi_vacinado[i]){
                            vacinado[i] = true;
                            foi_vacinado[i] = true;
                            N_vacinas--;
                        }
                    }
                }
                Sintomaticos--;
                break;
            case 'H': // Hospitalziado - 4
                if(prob_estagio[i]<morte[(int) VAN(&Grafo, "faixa", site)]){
                    estagio[i] = 'D';
                    Mortos++;
                }
                else{
                    estagio[i] = 'R';
                    Recuperados++;
                    if((N_vacinas!=0) && (ano >= dia_infecao+1 )){
                        if(!foi_vacinado[i]){
                            vacinado[i] = true;
                            foi_vacinado[i] = true;
                            N_vacinas--;
                        }
                    }
                }
                Hospitalizados--;
                break;
            case 'R': // Recuperados 
                estagio[i] = 'S';
                Suscetiveis++;
                Recuperados--;
                if((N_vacinas!=0) && (ano >= dia_infecao+1 )){
                    if(!foi_vacinado[i]){
                        vacinado[i] = true;
                        foi_vacinado[i] = true;
                        N_vacinas--;
                    }
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
            infect_time[s - 1][7] += Nhospitalizados/N;
            infect_time[s - 1][8] += Mortos*Mortos/(N*N);
            infect_time[s - 1][9] += Nhospitalizados*Nhospitalizados/(N*N);
            infect_time[s - 1][10] += Infectados/N;
            infect_time[s - 1][11] += Infectados*Infectados/(N*N);
            //if(redes <= cut_rede)printf("%f %f %f %f %f %f %f %f %d\n",ano,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos,s);
        }
        if(t >= dias+1) break;
        if((Expostos == Assintomaticos) && (Assintomaticos == Sintomaticos) && (Sintomaticos ==Hospitalizados) && (Hospitalizados == 0)) break;
        

    }
    free(prob_estagio);
    free(morte);
    free(hospitalizacao);
    free(sintomatico);
    free(estagio);
    free(vacinado);
    
    //if(redes <= cut_rede)printf("S = %d\n",s);
    if(ano <= dias){
        for (i = s-1; i < dias*2; i++){
            infect_time[i][0] += (double) Suscetiveis/N;
            infect_time[i][1] += (double) Expostos/N;
            infect_time[i][2] += (double) Assintomaticos/N;
            infect_time[i][3] += (double) Sintomaticos/N;
            infect_time[i][4] += (double) Hospitalizados/N;
            infect_time[i][5] += (double) Recuperados/N;
            infect_time[i][6] += (double) Mortos/N;
            infect_time[i][7] += (double) Nhospitalizados/N;
            infect_time[i][8] += (double) Mortos*Mortos/(N*N);
            infect_time[i][9] += (double) Nhospitalizados*Nhospitalizados/(N*N);
            infect_time[i][10] += (double) Infectados/N;
            infect_time[i][11] += (double) Infectados*Infectados/(N*N);
        }
    }

    free(rates);
    igraph_destroy(&Grafo);
    free(foi_vacinado);
    *seed += 1;
}



void generate_infect(double N,double p,int seed, int redes,double f,int vacina,bool weight){
    
    
    int cut_rede = 5;
    
    uint16_t i;
    create_folder((int)N);
    
    infect_time = (double**) malloc(tempo_total*sizeof(double*));
    for (i = 0; i < tempo_total; i++) infect_time[i] = (double*) calloc(q_resultados,sizeof(double));        
    for (i = 0; i < redes; i++) infect(10,N,p,&seed,infect_time,f,vacina,redes,weight,cut_rede);
    
    for (i = 0; i < tempo_total; i++){
        if(redes <= cut_rede){
            //print_vetor(infect_time[i],q_resultados,sizeof(double));
            //printf("%f\n",quant[i]);
        }
        for (int j = 0; j < q_resultados; j++) infect_time[i][j] /= redes;
        //if(redes <= cut_rede) print_vetor(infect_time[i],q_resultados,sizeof(double));
        //if(redes <= cut_rede) printf(" =============================================== \n");
    }
    
    if((f == 0)){
        if(redes > cut_rede){
            char filecheck[800];
            if(weight) sprintf(filecheck,"./output/time/%d/ponderado/p/infect_%.2f.txt",(int)N,p);
            else sprintf(filecheck,"./output/time/%d/nponderado/p/infect_%.2f.txt",(int)N,p);
            generate_file(filecheck,infect_time,tempo_total,q_resultados,sizeof(infect_time[0][0]));
        }
    }
    char *file_vacina[] = {"idade", "grau", "close", "harmonic","betwenness","eigenvector","eccentricity","clustering","kshell","graumorte","probhosp","probmorte","probhospassin","probmortepassin","random"};
    if((f == 0.75) ||(f == 0.25) ||(f == 0.50) || (f == 1.0)){
        if(redes > cut_rede){
            char filename[800];
            if(weight) sprintf(filename,"./output/time/%d/ponderado/%s_%.2f_%.2f.txt",(int)N,file_vacina[vacina],p,f);
            else sprintf(filename,"./output/time/%d/nponderado/%s_%.2f_%.2f.txt",(int)N,file_vacina[vacina],p,f);
            generate_file(filename,infect_time,tempo_total,q_resultados,sizeof(double));
        }
    }
    int s = f*100;
    if ((redes > cut_rede) && (f !=0)){
        FILE *file;
        char filename[800];
        if(weight) sprintf(filename,"./output/vacina/%d/ponderado/%s_%.2f.txt",(int)N,file_vacina[vacina],p);
        else sprintf(filename,"./output/vacina/%d/nponderado/%s_%.2f.txt",(int)N,file_vacina[vacina],p);
        if(s%10 == 0)printf("%s ",file_vacina[vacina]);
        file = fopen(filename,"a");
        fprintf(file,"%f %f %f %f %f %f %f\n",f, infect_time[tempo_total-1][7], infect_time[tempo_total-1][6], infect_time[tempo_total-1][10],pow(infect_time[tempo_total-1][9] - pow(infect_time[tempo_total-1][7],2),0.5), pow(infect_time[tempo_total-1][8] - pow(infect_time[tempo_total-1][6],2),0.5),pow(infect_time[tempo_total-1][11] - pow(infect_time[tempo_total-1][10],2),0.5));
        fclose(file);
    }
    else printf("%.2f %f %f %f %f\n",f, infect_time[tempo_total-1][7], infect_time[tempo_total-1][6],pow(infect_time[tempo_total-1][9] - pow(infect_time[tempo_total-1][7],2),0.5), pow(infect_time[tempo_total-1][8] - pow(infect_time[tempo_total-1][6],2),0.5));
    //else printf("%.2f %f %f %f %f\n",f, infect_time[tempo_total-1][7], infect_time[tempo_total-1][6],pow(final[2] - final[0]*final[0],0.5), pow(final[3] - final[1]*final[1],0.5));

    for (i = 0; i < tempo_total; i++) free(infect_time[i]);
    free(infect_time);
    //if (redes != 1) printf("\e[1;1H\e[2J");
    
    if(s%10 == 0)printf("Terminou: %f %f\n",f,p);
}