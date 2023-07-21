#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

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

bool* vacinacao_clustering(bool* vacinado,struct Graph G, int Suscetiveis, int Recuperados, double f, uint8_t* estagio){
    
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    double* clustering;
    clustering = list_clustering(G.viz,G.Nodes);
    int* sitio = (int*) malloc(G.Nodes*sizeof(int));
    uint16_t i;

    for (i = 0; i < G.Nodes; i++) sitio[i] = i;
    sortIntByRef(sitio,clustering,G.Nodes,sizeof(clustering[0]));

    for (i = G.Nodes-1; i >=0 ; i--){
        if(Vac == 0) break;
        if((estagio[sitio[i]] == 0) || (estagio[sitio[i]] == 5)){
            vacinado[sitio[i]] = true;
            Vac--;
        }
    }
    free(sitio);
    free(clustering);
    return vacinado;
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

bool* vacinacao_proxima(bool* vacinado,struct Graph G, int Suscetiveis, int Recuperados, double f, uint8_t* estagio){
    
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    double* proxima = (double*) malloc(G.Nodes*sizeof(double));
    int* sitio = (int*) malloc(G.Nodes*sizeof(int));
    uint16_t i;
    int d;

    for (i = 0; i < G.Nodes; i++){
        proxima[i] = shortest_length(G.viz,G.Nodes,i,&d);
        if(proxima[i] != 0) proxima[i] = 1/proxima[i];
        sitio[i] = i;
    }
    sortIntByRef(sitio,proxima,G.Nodes,sizeof(proxima[0]));
    i = 0;
    for (i = G.Nodes-1; i >=0 ; i--){
        if(Vac == 0) break;
        if((estagio[sitio[i]] == 0) || (estagio[sitio[i]] == 5)){
            vacinado[sitio[i]] = true;
            Vac--;
        }
    }
    free(sitio);
    free(proxima);
    return vacinado;
}

bool* vacinacao_excentrica(bool* vacinado,struct Graph G, int Suscetiveis, int Recuperados, double f, uint8_t* estagio){
    
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    double* excentrica = (double*) malloc(G.Nodes*sizeof(double));
    int* sitio = (int*) malloc(G.Nodes*sizeof(int));
    uint16_t i;
    int d;

    for (i = 0; i < G.Nodes; i++){
        double _ = shortest_length(G.viz,G.Nodes,i,&d);
        if(d != 0) excentrica[i] = (double)1/d;
        sitio[i] = i;
    }
    sortIntByRef(sitio,excentrica,G.Nodes,sizeof(excentrica[0]));
    //for (i = 0; i < G.Nodes; i++) printf("%d %f %f\n",sitio[i],excentrica[i],excentrica[sitio[i]]);
    i = 0;
    for (i = G.Nodes-1; i >=0 ; i--){
        if(Vac == 0) break;
        if((estagio[sitio[i]] == 0) || (estagio[sitio[i]] == 5)){
            vacinado[sitio[i]] = true;
            Vac--;
        }
    }
    free(sitio);
    free(excentrica);
    return vacinado;
}

bool* vacinacao_aleatoria(bool* vacinado, int Suscetiveis, int Recuperados, double f, uint8_t* estagio,int N){
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    while(Vac!=0){
        double r = genrand64_real1();
        int N0 = r*N;
        if((estagio[N0] == 0) || (estagio[N0] == 5)){
            if(!vacinado[N0]){
                vacinado[N0] = true;
                Vac--;
            }
        }
    }
    return vacinado;
}

bool* vacinacao_grau(bool* vacinado,struct Graph G, int Suscetiveis, int Recuperados, double f, uint8_t* estagio){
    
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    int* grau = (int*) malloc(G.Nodes*sizeof(int));
    int* sitio = (int*) malloc(G.Nodes*sizeof(int));
    uint16_t i;

    for (i = 0; i < G.Nodes; i++) {grau[i] = G.viz[i][0];sitio[i] = i;}
    sortIntByRef(sitio,grau,G.Nodes,sizeof(grau[0]));
    free(grau);

    for (i = G.Nodes-1; i >=0 ; i--){
        if(Vac == 0) break;
        if((estagio[sitio[i]] == 0) || (estagio[sitio[i]] == 5)){
            vacinado[sitio[i]] = true;
            Vac--;
        }
    }
    free(sitio);
    return vacinado;
}

bool* vacinacao_idade(bool* vacinado,int* faixas, int Suscetiveis, int Recuperados, double f, uint8_t* estagio,int N){
    uint16_t i;
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    int* sitio = (int*) malloc(N*sizeof(int));

    for ( i = 0; i < N; i++) sitio[i] = i;

    sortIntByRef(sitio,faixas,N,sizeof(faixas[0]));

    for (i = N-1; i >=0 ; i--){
        if(Vac == 0) break;
        if((estagio[sitio[i]] == 0) || (estagio[sitio[i]] == 5)){
            vacinado[sitio[i]] = true;
            Vac--;
        }
    }
    free(sitio);
    return vacinado;
}

double calc_estagio(int site,uint8_t* estagio,double* prob_estagio,int* faixas, struct Graph G,bool* vacinado,double* hospitalizacao,double* morte){


    switch (estagio[site]){
        case 0:// Suscetível
            double beta = 0;

            for (int j = 0; j < G.viz[site][0]; j++){
                int vizinho = G.viz[site][j+1];
                if(estagio[vizinho] == 1){
                    if(!vacinado[vizinho]) beta += beta1;
                    else beta += 0.058*beta1; // 0.058
                }
                if(estagio[vizinho] == 2){
                    if(!vacinado[vizinho]) beta += beta2;
                    else beta += 0.058*beta2;// 0.058
                }
            }

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
            if(prob_estagio[site]< hospitalizacao[faixas[site]]*(vacinado[site] ? 0.034 : 1)) return phi;
            else return gamma_I;
            break;
        case 4: // Hospitalziado
            if(prob_estagio[site]<morte[faixas[site]]*(vacinado[site] ? 0.034 : 1)) return delta;
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

void infect(int S0,int E0,int N,int seed,double** infect_time,double* quant,double* final,double f,int vacina){
    uint16_t Suscetiveis = S0;
    uint16_t Expostos = E0;
    uint16_t Assintomaticos = 0;
    uint16_t Sintomaticos = 0;
    uint16_t Hospitalizados = 0;
    uint16_t Recuperados = 0;
    uint16_t Mortos = 0;
    bool Vac = false;
    int site;

    struct Graph G;
    G = local_configuration_model(N,0,seed,0);
    char file_name[200] = "./output/teste0.txt";
    int* faixas = get_faixas(G.Nodes,file_name);
    double* prob_estagio = (double*) malloc(G.Nodes*sizeof(double));
    double* rates = (double*) calloc(G.Nodes,sizeof(double));
    bool inicio = true;
    uint8_t* estagio = calloc(G.Nodes,sizeof(uint8_t));
    bool* vacinado = calloc(G.Nodes,sizeof(bool));
    init_genrand64(seed);

    int i = 0;
    for (i = 0; i < G.Nodes; i++){
        prob_estagio[i] = genrand64_real1();
        vacinado[i] = false;
    }
    while(E0!=0){
        double r = genrand64_real1();
        int N0 = r*G.Nodes;
        if(estagio[N0] == 0){
            estagio[N0] = 1;
            E0--;
        }
    }
    double t = 0.5;
    uint16_t s = 1;

    
    double* sintomatico = (double*) malloc(5*sizeof(double));
    double* hospitalizacao = (double*) malloc(5*sizeof(double));
    double* morte = (double*) malloc(5*sizeof(double));
    sintomatico[0] = 29.1/100;
    sintomatico[1] = 37.4/100;
    sintomatico[2] = 41.68/100;
    sintomatico[3] = 39.4/100;
    sintomatico[4] = 31.3/100; 
    /* sintomatico[0] = (double)47.3/100; // 0 - 20
    sintomatico[1] = (double)45.8/100; // 20 -  30
    sintomatico[2] = (double)0.4673015873; // 30 - 50
    sintomatico[3] = (double)0.4948096886; // 50 -70
    sintomatico[4] = (double)0.4; // 70+ */

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


    uint16_t Nhospitalizados = 0;
    double tempo = 0;
    double ano = 0;
    double rate = 0,Delta;
    double new_rate;
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
            switch (vacina){
            case 0:
                vacinado = vacinacao_aleatoria(vacinado,Suscetiveis, Recuperados,f, estagio,G.Nodes);
                break;
            case 1:
                vacinado = vacinacao_grau(vacinado,G,Suscetiveis, Recuperados,f, estagio);
                break;
            case 2:
                vacinado = vacinacao_idade(vacinado,faixas,Suscetiveis, Recuperados,f, estagio,G.Nodes);
                break;
            case 3:
                vacinado = vacinacao_excentrica(vacinado,G,Suscetiveis,Recuperados, f,estagio);
                break;
            case 4:
                vacinado = vacinacao_proxima(vacinado,G,Suscetiveis,Recuperados, f,estagio);
                break;
            case 5:
                vacinado = vacinacao_probability( estagio, prob_estagio,faixas, G, vacinado, hospitalizacao, morte, sintomatico, Suscetiveis,Recuperados, f);
                break;
            case 6:
                vacinado = vacinacao_clustering(vacinado,G,Suscetiveis,Recuperados, f,estagio);
                break;
            default:
                break;
            }
            
            Vac = true;
        }
        //rate = 0;
        //for (i = 0; i < G.Nodes; i++)  rate += calc_estagio(i, estagio,prob_estagio,faixas,G,vacinado,hospitalizacao,morte);
        if(inicio){
            for (i = 0; i < G.Nodes; i++) {
                rates[i] = calc_estagio(i, estagio,prob_estagio,faixas,G,vacinado,hospitalizacao,morte);
                rate += rates[i];
            }
            inicio = false;
        }
        else{
            rate -= rates[site];
            rates[site] = calc_estagio(site, estagio,prob_estagio,faixas,G,vacinado,hospitalizacao,morte);
            rate += rates[site];
            for (int j = 0; j < G.viz[site][0]; j++){
                int vizinho = G.viz[site][j+1];
                rate -= rates[vizinho];
                rates[vizinho] = calc_estagio(vizinho, estagio,prob_estagio,faixas,G,vacinado,hospitalizacao,morte);
                rate += rates[vizinho];
            }
            
        }
        if(rate==0) break;
        tempo = exponentialRand(rate);
        if(tempo ==0) break;
        Delta = rate*genrand64_real1();
        new_rate = 0;
        ano += tempo;
        if(ano > 3*365) break;

        for (i = 0; i < G.Nodes; i++){
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
                if(prob_estagio[i]<sintomatico[faixas[i]]*(vacinado[i] ? 0.34693877551 : 1)){
                    estagio[i] = 3;
                    Sintomaticos++;
                }
                else{
                    estagio[i] = 2;
                    Assintomaticos++;
                }
                Expostos--;
                break;
            case 2: //Assintomático
                estagio[i] = 5;
                Assintomaticos--;
                Recuperados++;
                break;
            case 3: // Sintomático
                if(prob_estagio[i]<hospitalizacao[faixas[i]]){
                    estagio[i] = 4;
                    Hospitalizados++;
                    Nhospitalizados++;
                }
                else{
                    estagio[i] = 5;
                    Recuperados++;
                }
                Sintomaticos--;
                break;
            case 4: // Hospitalziado
                if(prob_estagio[i]<morte[faixas[i]]){
                    estagio[i] = 6;
                    Mortos++;
                }
                else{
                    estagio[i] = 5;
                    Recuperados++;
                }
                Hospitalizados--;
                break;
            case 5:
                estagio[i] = 0;
                Suscetiveis++;
                Recuperados--;
                break;
            default:
                break;
        }
        prob_estagio[i] = genrand64_real1();
        if(ano > t) {s++;t += 0.5;}

        infect_time[s - 1][0] += Suscetiveis;
        infect_time[s - 1][1] += Expostos;
        infect_time[s - 1][2] += Assintomaticos;
        infect_time[s - 1][3] += Sintomaticos;
        infect_time[s - 1][4] += Hospitalizados;
        infect_time[s - 1][5] += Recuperados;
        infect_time[s - 1][6] += Mortos;
        quant[s -1]++;
        
        if((Expostos == Assintomaticos) && (Assintomaticos == Sintomaticos) && (Sintomaticos ==Hospitalizados) && (Hospitalizados == 0)) break;

    }
    final[0] += (double) Nhospitalizados;
    final[1] += (double) Mortos;
    final[2] += (double) Nhospitalizados*Nhospitalizados;
    final[3] += (double) Mortos*Mortos;
    free(prob_estagio);
    free(morte);
    free(faixas);
    free(hospitalizacao);
    free(sintomatico);
    free(estagio);
    free(vacinado);
    free(rates);
    for(int j = 0; j < N; j++)free(G.viz[j]);
    free(G.viz);
    //printf("\e[1;1H\e[2J");
}

void generate_infect(int seed, int redes,double f,int vacina){
    char file[200] = "./output/teste.txt";
    uint16_t N = size_txt(file)-1;
    uint16_t tempo = 365*3*2;
    uint16_t i,j;

    //double* resultados = (double*) malloc(2*sizeof(double));
    infect_time = (double**) malloc(tempo*sizeof(double*));
    quant = (double*) calloc(tempo,sizeof(double));
    for (int i = 0; i < tempo; i++){
        infect_time[i] = (double*) malloc(7*sizeof(double));
        for (int j = 0; j < 7; j++) infect_time[i][j] = 0;
        quant[i] = 0;
    }
    double* final = (double*) malloc(4*sizeof(double));
    final[0] = 0;
    final[1] = 0;
    final[2] = 0;
    final[3] = 0;

    #pragma omp parallel for
    for (j = 0; j < redes; j++) infect(N-1,1,N,seed+j ,infect_time,quant,final,f,vacina);
    for (i = 0; i < tempo; i++) for (int j = 0; j < 7; j++) infect_time[i][j] /= quant[i];
    //for (i = 0; i < tempo; i++) infect_time[i][8] = pow(infect_time[i][8] - pow(infect_time[i][7],2),0.5);
    for(i = 0;i < 4;i++) final[i] /= redes;
    
    if((f==0.5) && (redes != 1)){
        char filename[500];
        switch (vacina){
            case 0:
                sprintf(filename,"./time/infect_aleatorio.txt");
                break;
            case 1:
                sprintf(filename,"./time/infect_grau.txt");
                break;
            case 2:
                sprintf(filename,"./time/infect_idade.txt");
                break;
            case 3:
                sprintf(filename,"./time/infect_excentrico.txt");
                break;
            case 4:
                sprintf(filename,"./time/infect_proximo.txt");
                break;
            case 5:
                sprintf(filename,"./time/infect_probability_hospitalizado.txt");
                break;
            case 6:
                sprintf(filename,"./time/infect_clustering.txt");
                break;
            default:
                sprintf(filename,"./time/infect_without.txt");
                break;
        }
        generate_file(filename,infect_time,tempo,7,sizeof(infect_time[0][0]));
    }

    if (redes != 1){
        FILE *file;
        char filename[200];
        switch (vacina){
            case 0:
                sprintf(filename,"./vacina/infect_vacina_aleatorio.txt");
                break;
            case 1:
                sprintf(filename,"./vacina/infect_vacina_grau.txt");
                break;
            case 2:
                sprintf(filename,"./vacina/infect_vacina_idade.txt");
                break;
            case 3:
                sprintf(filename,"./vacina/infect_vacina_excentrico.txt");
                break;
            case 4:
                sprintf(filename,"./vacina/infect_vacina_proximo.txt");
                break;
            case 5:
                sprintf(filename,"./vacina/infect_vacina_probability_hospitalizado.txt");
                break;
            case 6:
                sprintf(filename,"./vacina/infect_vacina_clustering.txt");
                break;
            default:
                sprintf(filename,"./vacina/infect_vacina_without.txt");
                break;
        }
        file = fopen(filename,"a");
        fprintf(file,"%f %f %f %f %f\n",f, final[0], final[1],pow(final[2] - final[0]*final[0],0.5), pow(final[3] - final[1]*final[1],0.5));
        fclose(file);
    }
    else printf("%f %f %f %f %f\n",f, final[0], final[1],pow(final[2] - final[0]*final[0],0.5), pow(final[3] - final[1]*final[1],0.5));

    for (i = 0; i < tempo; i++) free(infect_time[i]);
    free(infect_time);
    free(final);
    free(quant);

    //if (redes != 1) printf("\e[1;1H\e[2J");
    printf("Terminou: %f\n",f);
}