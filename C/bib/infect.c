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
#include "rede.h"
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
double* quant;
int foi = 0;
const uint16_t dias = 150;
bool fileExists(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

bool* vacinacao_probability(uint8_t* estagio,double* prob_estagio,int* faixas, struct Graph G,bool* vacinado,double* hospitalizacao,double* morte,double* sintomatico, int Suscetiveis, int Recuperados, double f){

    uint16_t Vac = f*(Suscetiveis + Recuperados);

    double* probability = (double*) malloc(G.Nodes*sizeof(double));
    int* sitio = (int*) malloc(G.Nodes*sizeof(int));

    uint16_t i;
    double soma = 0;

    for (i = 0; i < G.Nodes; i++){
        probability[i] = 0;
        for (int j = 0; j < G.viz[i][0]; j++) probability[i] += sintomatico[faixas[G.viz[i][j+1]]]*hospitalizacao[faixas[G.viz[i][j+1]]]/**morte[faixas[G.viz[i][j+1]]]*/;
        probability[i] *= (1 - sintomatico[faixas[i]]);
        //probability[i] += sintomatico[faixas[i]]*morte[faixas[i]]*hospitalizacao[faixas[i]];
        sitio[i] = i;
    }

    sortIntByRef(sitio,probability,G.Nodes,sizeof(probability[0]));
    for (i = 0; i < G.Nodes; i++)
    for (i = G.Nodes-1; i >=0 ; i--){
        if(Vac == 0) break;
        if((estagio[sitio[i]] == 0) || (estagio[sitio[i]] == 5)){
            vacinado[sitio[i]] = true;
            Vac--;
        }
    }

    free(probability);
    free(sitio);
    return vacinado;
}

igraph_vector_int_t centrality(igraph_t* Grafo, uint8_t* estagio,int check){
    uint16_t i;
    int N = igraph_vcount(Grafo);
    igraph_vector_int_t centralidade;
    igraph_vector_int_init(&centralidade, N);
    switch (check){
        
        case 0: // Vacinação por idade
            
            igraph_vector_t faixas;
            igraph_vector_init(&faixas, 0);
            igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);
            igraph_vector_qsort_ind(&faixas,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&faixas);
            break;
        case 1: // Vacinação por grau
            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);
            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
            igraph_vector_int_qsort_ind(&degrees,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_int_destroy(&degrees);
            break;
        case 2:
            
            igraph_vector_t closeness;
            igraph_vector_int_t reachable_count;
            igraph_bool_t all_reachable;

            igraph_vector_int_init(&reachable_count, 0);
            igraph_vector_init(&closeness, N);

            igraph_closeness(Grafo, &closeness, &reachable_count, &all_reachable, igraph_vss_all(), IGRAPH_ALL, NULL, 1);
            igraph_vector_qsort_ind(&closeness,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&closeness);
            igraph_vector_int_destroy(&reachable_count);
            
            break;
        case 3: // Vacinação por Harmo
            igraph_vector_t Harmonic;
            igraph_vector_init(&Harmonic, N);
            igraph_harmonic_centrality(Grafo,&Harmonic,igraph_vss_all(),IGRAPH_ALL,NULL,0);
            igraph_vector_qsort_ind(&Harmonic,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Harmonic);
            break;
        case 4:
            igraph_vector_t Betweenness;
            igraph_vector_init(&Betweenness, N);
            igraph_betweenness(Grafo, &Betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL);
            igraph_vector_qsort_ind(&Betweenness,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Betweenness);
            break;
        case 5:
            igraph_vector_t eigenvector;
            igraph_real_t autovalor;
            igraph_vector_init(&eigenvector, 0);
            igraph_eigenvector_centrality(Grafo,&eigenvector,0,IGRAPH_UNDIRECTED,1,NULL,NULL);
            igraph_vector_qsort_ind(&eigenvector,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&eigenvector);
            break;
        case 6:
            igraph_vector_t eccentricity;
            igraph_vector_init(&eccentricity, N);
            igraph_eccentricity(Grafo,&eccentricity,igraph_vss_all(),IGRAPH_ALL);
            igraph_vector_qsort_ind(&eccentricity,&centralidade, IGRAPH_ASCENDING);
            igraph_vector_destroy(&eccentricity);
            break;
        case 7:
            igraph_vector_t clustering;
            igraph_vector_init(&clustering, N);
            igraph_transitivity_local_undirected(Grafo,&clustering,igraph_vss_all(),IGRAPH_TRANSITIVITY_ZERO);
            igraph_vector_qsort_ind(&clustering,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&clustering);
            break;
        case 8:
            igraph_vector_int_t k_shell;
            igraph_vector_int_init(&k_shell, N);
            igraph_coreness(Grafo, &k_shell, IGRAPH_ALL);
            igraph_vector_int_qsort_ind(&k_shell,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_int_destroy(&k_shell);
            break;
        default:// Aleatório
            igraph_vector_int_init_range(&centralidade,0,N);
            igraph_vector_int_shuffle(&centralidade);
            break;
    }
    return centralidade;
    
    //return vacinado;
}

double calc_estagio(int site,uint8_t* estagio,double* prob_estagio, igraph_t* Grafo,bool* vacinado,double* hospitalizacao,double* morte){


    switch (estagio[site]){
        case 0:// Suscetível
            double beta = 0;

            igraph_vector_int_t vizinhos;
            igraph_vector_int_init(&vizinhos, 0);
            igraph_neighbors(Grafo, &vizinhos, site,IGRAPH_ALL);

            for (int j = 0; j < igraph_vector_int_size(&vizinhos); j++){
                int vizinho = VECTOR(vizinhos)[j];
                if(estagio[vizinho] == 1){
                    if(!vacinado[vizinho]) beta += beta1;
                    else beta += 0.058*beta1; // 0.058
                }
                if(estagio[vizinho] == 2){
                    if(!vacinado[vizinho]) beta += beta2;
                    else beta += 0.058*beta2;// 0.058
                }
            }
            igraph_vector_int_destroy(&vizinhos);
            if(!vacinado[site]) return beta;
            else return 0.058*beta;// 0.058
            break;
        case 1: // Exposto
            return sigma;
            break;
        case 2: //Assintomático
            return gamma_A;
            break;
        case 3: // Sintomático
            if(prob_estagio[site]< hospitalizacao[(int) VAN(Grafo, "faixa", site)]*(vacinado[site] ? 0.034 : 1)) return phi;
            else return gamma_I;
            break;
        case 4: // Hospitalziado
            if(prob_estagio[site]<morte[(int) VAN(Grafo, "faixa", site)]*(vacinado[site] ? 0.034 : 1)) return delta;
            else return gamma_H;
            break;
        case 5:
            return recupera;
            break;
        default:
            return 0;
            break;
    }
}

void infect(int E0,int N0,double p,int seed,double** infect_time,double* quant,double* final,double f,int vacina,int redes){

    igraph_t Grafo = local_configuration_model(N0,p,seed);
    int N = igraph_vcount(&Grafo);

    int N_vacinas = f*N;

    uint16_t Suscetiveis = N - E0;
    uint16_t Expostos = E0;
    uint16_t Assintomaticos = 0;
    uint16_t Sintomaticos = 0;
    uint16_t Hospitalizados = 0;
    uint16_t Recuperados = 0;
    uint16_t Mortos = 0;
    bool Vac = false;
    int site;

    double* prob_estagio = (double*) malloc(N*sizeof(double));
    double* rates = (double*) calloc(N,sizeof(double));
    uint8_t* estagio = calloc(N,sizeof(uint8_t));
    bool* vacinado = calloc(N,sizeof(bool));
    bool* foi_vacinado = calloc(N,sizeof(bool));
    init_genrand64(seed);

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

    morte[0] = (double)0.408/100;
    morte[1] = (double)1.04/100;
    morte[2] = (double)3.89/100;
    morte[3] = (double)9.98/100;
    morte[4] = (double)17.5/100;

    int i = 0;
    double rate = 0;
    for (i = 0; i < N; i++) prob_estagio[i] = genrand64_real1();
    while(E0!=0){
        double r = genrand64_real1();
        site = r*N;
        if(estagio[site] == 0){
            estagio[site] = 1;
            E0--;
            rates[site] = calc_estagio(site, estagio,prob_estagio,&Grafo,vacinado,hospitalizacao,morte);
            rate += rates[site];
        }
    }
    double t = 0.5;
    int j;
    uint16_t s = 1;

    uint16_t Nhospitalizados = 0;
    double tempo = 0;
    double ano = 0;
    double Delta;
    double new_rate;
    int original = site;
    while (s != 0){
        //rate = 0;
        /* if(((int)ano == 60) || (ano == 0)){
            FILE *file;
            char filename[200];
            sprintf(filename,"./output/infect/grau_distribution_%d.txt",(int)ano);
            file = fopen(filename,"w");
            for (int k = 0; k < G.Nodes; k++) if((estagio[k] == 0) || (estagio[k] == 5)) fprintf(file,"%d %d\n",G.viz[k][0],faixas[k]);
            fclose(file);
        } */
        if((ano>= 61) && (!Vac) && (f!= 0.0)){
            igraph_vector_int_t centralidade = centrality(&Grafo,estagio,vacina);
            for (j = 0; j < N; j++){
                int sitio_vacinado = VECTOR(centralidade)[j];
                if((estagio[sitio_vacinado] == 0) || (estagio[sitio_vacinado] == 5) || (estagio[sitio_vacinado] == 2)){
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
        
        rate -= rates[site];
        rates[site] = calc_estagio(site, estagio,prob_estagio,&Grafo,vacinado,hospitalizacao,morte);
        rate += rates[site];

        igraph_vector_int_t vizinhos;
        igraph_vector_int_init(&vizinhos, 0);
        igraph_neighbors(&Grafo, &vizinhos, site,IGRAPH_ALL);

        for (j = 0; j < igraph_vector_int_size(&vizinhos); j++){
            int vizinho = VECTOR(vizinhos)[j];
            rate -= rates[vizinho];
            rates[vizinho] = calc_estagio(vizinho, estagio,prob_estagio,&Grafo,vacinado,hospitalizacao,morte);
            rate += rates[vizinho];
        }
        igraph_vector_int_destroy(&vizinhos);
        
        /*rate = 0;
        for (i = 0; i < N; i++){
            rates[i] = calc_estagio(i, estagio,prob_estagio,&Grafo,vacinado,hospitalizacao,morte);
            rate += rates[i];
        }*/
        if(rate==0) break;
        tempo = exponentialRand(rate);
        if(tempo ==0) break;
        Delta = rate*genrand64_real1();
        new_rate = 0;
        ano += tempo;
        if(ano > dias) break;
        if(rate<=0) break;
        for (i = 0; i < N; i++){
            new_rate += rates[i];
            if(new_rate>=Delta) break;
        }
        site = i;
        switch (estagio[i]){
            case 0:// Sucetível
                estagio[i] = 1;
                Expostos++;
                Suscetiveis--;
                break;
            case 1: // Exposto
                if(prob_estagio[i]<sintomatico[(int) VAN(&Grafo, "faixa", i)]*(vacinado[i] ? 0.34693877551 : 1)){
                    estagio[i] = 3;
                    Sintomaticos++;
                }
                else{
                    estagio[i] = 2;
                    Assintomaticos++;
                    if((N_vacinas!=0) && (ano >= 61 )){
                        if(!foi_vacinado[i]){
                            vacinado[i] = true;
                            foi_vacinado[i] = true;
                            N_vacinas--;
                        }
                    }
                    
                }
                Expostos--;
                break;
            case 2: //Assintomático
                estagio[i] = 5;
                Assintomaticos--;
                Recuperados++;
                if((N_vacinas!=0) && (ano >= 61 )){
                    if(!foi_vacinado[i]){
                        vacinado[i] = true;
                        foi_vacinado[i] = true;
                        N_vacinas--;
                    }
                }
                break;
            case 3: // Sintomático
                if(prob_estagio[i]<hospitalizacao[(int) VAN(&Grafo, "faixa", site)]){
                    estagio[i] = 4;
                    Hospitalizados++;
                    Nhospitalizados++;
                }
                else{
                    estagio[i] = 5;
                    Recuperados++;
                    if((N_vacinas!=0) && (ano >= 61 )){
                        if(!foi_vacinado[i]){
                            vacinado[i] = true;
                            foi_vacinado[i] = true;
                            N_vacinas--;
                        }
                    }
                }
                Sintomaticos--;
                break;
            case 4: // Hospitalziado
                if(prob_estagio[i]<morte[(int) VAN(&Grafo, "faixa", site)]){
                    estagio[i] = 6;
                    Mortos++;
                }
                else{
                    estagio[i] = 5;
                    Recuperados++;
                    if((N_vacinas!=0) && (ano >= 61 )){
                        if(!foi_vacinado[i]){
                            vacinado[i] = true;
                            foi_vacinado[i] = true;
                            N_vacinas--;
                        }
                    }
                }
                Hospitalizados--;
                break;
            case 5:
                estagio[i] = 0;
                Suscetiveis++;
                Recuperados--;
                if((N_vacinas!=0) && (ano >= 61 )){
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
        if(ano > t) {s++;t += 0.02;}
        if(redes == 1) printf("%f %d %d %d %d %d %d %d %f %f\n",ano,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos,quant[s-1], tempo);
        if((f == 0.5) || (f == 1.0)|| (f == 0.0)){
            infect_time[s - 1][0] += (double)Suscetiveis/N;
            infect_time[s - 1][1] += (double)Expostos/N;
            infect_time[s - 1][2] += (double)Assintomaticos/N;
            infect_time[s - 1][3] += (double)Sintomaticos/N;
            infect_time[s - 1][4] += (double)Hospitalizados/N;
            infect_time[s - 1][5] += (double)Recuperados/N;
            infect_time[s - 1][6] += (double)Mortos/N;
            quant[s -1]++;
        }
        if((Expostos == Assintomaticos) && (Assintomaticos == Sintomaticos) && (Sintomaticos ==Hospitalizados) && (Hospitalizados == 0)) break;
        

    }
    
    final[0] += (double) Nhospitalizados/N;
    final[1] += (double) Mortos/N;
    final[2] += (double) Nhospitalizados*Nhospitalizados/(N*N);
    final[3] += (double) Mortos*Mortos/(N*N);
    free(prob_estagio);
    free(morte);
    free(hospitalizacao);
    free(sintomatico);
    free(estagio);
    free(vacinado);
    
    free(rates);
    igraph_destroy(&Grafo);
    free(foi_vacinado);
    //printf("\e[1;1H\e[2J");
    //printf("%d %f\n",foi,ano);
    //foi++;
}

void generate_infect(uint16_t N,double p,int seed, int redes,double f,int vacina){
    //char file[200] = "./output/teste.txt";
    //uint16_t N = size_txt(file)-1;
    
    int tempo = dias*50;
    uint16_t i,j;

    //double* resultados = (double*) malloc(2*sizeof(double));
    infect_time = (double**) malloc(tempo*sizeof(double*));
    quant = (double*) calloc(tempo,sizeof(double));
    for (int i = 0; i < tempo; i++) infect_time[i] = (double*) calloc(7,sizeof(double));
    double* final = (double*) malloc(4*sizeof(double));
    final[0] = 0;
    final[1] = 0;
    final[2] = 0;
    final[3] = 0;

    for (j = 0; j < redes; j++) infect(1,N,p,seed+j ,infect_time,quant,final,f,vacina,redes);
    if((f == 0.5) || (f == 1.0) || (f == 0.0)){
        for (i = 0; i < tempo; i++){
            if(quant[i] == 0){
                tempo = i+1;
                infect_time = (double**) realloc(infect_time,tempo*sizeof(double*));
                break;
            }
            for (int j = 0; j < 7; j++) infect_time[i][j] /= quant[i];
        }
    }
    //for (i = 0; i < tempo; i++) infect_time[i][8] = pow(infect_time[i][8] - pow(infect_time[i][7],2),0.5);
    for(i = 0;i < 4;i++) final[i] /= redes;
    char *filecheck = "./time/infect.txt";
    if((f == 0)) generate_file(filecheck,infect_time,tempo,7,sizeof(infect_time[0][0]));

    if((f == 0.5) || (f == 1.0)){
        if(redes != 1){
            char filename[800];
            switch (vacina){
                case 0:
                    sprintf(filename,"./time/%d/infect_idade_%.2f_%.2f.txt",N,p,f);
                    break;
                case 1:
                    sprintf(filename,"./time/%d/infect_grau_%.2f_%.2f.txt",N,p,f);
                    break;
                case 2:
                    sprintf(filename,"./time/%d/infect_close_%.2f_%.2f.txt",N,p,f);
                    break;
                case 3:
                    sprintf(filename,"./time/%d/infect_harmonic_%.2f_%.2f.txt",N,p,f);
                    break;
                case 4:
                    sprintf(filename,"./time/%d/infect_betwenness_%.2f_%.2f.txt",N,p,f);
                    break;
                case 5:
                    sprintf(filename,"./time/%d/infect_eigenvector_%.2f_%.2f.txt",N,p,f);
                    break;
                case 6:
                    sprintf(filename,"./time/%d/infect_eccentricity_%.2f_%.2f.txt",N,p,f);
                    break;
                case 7:
                    sprintf(filename,"./time/%d/infect_clustering_%.2f_%.2f.txt",N,p,f);
                    break;
                case 8:
                    sprintf(filename,"./time/%d/infect_kshell_%.2f_%.2f.txt",N,p,f);
                    break;
                default:
                    sprintf(filename,"./time/%d/infect_random_%.2f_%.2f.txt",N,p,f);
                    break;
            }
            generate_file(filename,infect_time,tempo,7,sizeof(double));
        }
    }

    if (redes != 1){
        FILE *file;
        char filename[800];
        switch (vacina){
            case 0:
                sprintf(filename,"./vacina/%d/infect_vacina_idade_%.2f.txt",N,p);
                break;
            case 1:
                sprintf(filename,"./vacina/%d/infect_vacina_grau_%.2f.txt",N,p);
                break;
            case 2:
                sprintf(filename,"./vacina/%d/infect_vacina_close_%.2f.txt",N,p);
                break;
            case 3:
                sprintf(filename,"./vacina/%d/infect_vacina_harmonic_%.2f.txt",N,p);
                break;
            case 4:
                sprintf(filename,"./vacina/%d/infect_vacina_betweenness_%.2f.txt",N,p);
                break;
            case 5:
                sprintf(filename,"./vacina/%d/infect_vacina_eigenvector_%.2f.txt",N,p);
                break;
            case 6:
                sprintf(filename,"./vacina/%d/infect_vacina_eccentricity_%.2f.txt",N,p);
                break;
            case 7:
                sprintf(filename,"./vacina/%d/infect_vacina_clustering_%.2f.txt",N,p);
                break;
            case 8:
                sprintf(filename,"./vacina/%d/infect_vacina_kshell_%.2f.txt",N,p);
                break;
            default:
                sprintf(filename,"./vacina/%d/infect_vacina_random_%.2f.txt",N,p);
                break;
        }
        file = fopen(filename,"a");
        fprintf(file,"%f %f %f %f %f\n",f, final[0], final[1],pow(final[2] - final[0]*final[0],0.5), pow(final[3] - final[1]*final[1],0.5));
        fclose(file);
    }
    else printf("%.2f %f %f %f %f\n",f, final[0], final[1],pow(final[2] - final[0]*final[0],0.5), pow(final[3] - final[1]*final[1],0.5));

    for (i = 0; i < tempo; i++) free(infect_time[i]);
    free(infect_time);
    free(final);
    free(quant);
    //if (redes != 1) printf("\e[1;1H\e[2J");
    int s = f*100;
    if(s%10 == 0)printf("Terminou: %f %f\n",f,p);
}