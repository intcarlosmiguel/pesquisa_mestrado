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

int rodou;
const double beta1 = (double)0.5;
const double beta2 = (double)0.41;
const double sigma = (double)1/5.1;
const double gamma_A = (double)1/7;
const double phi = (double)0.025;
const double gamma_I = (double)1/7;
const double gamma_H = (double)1/14;
const double delta = (double)1/14;
const double recupera = (double)1/40;


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
    //quickSort(grau,sitio,0, G.Nodes-1);
    sitio = bubble_sort_by(sitio,grau,G.Nodes);
    i = 0;
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
    
    uint16_t Vac = f*(Suscetiveis + Recuperados);
    int* sitio = (int*) malloc(N*sizeof(int));

    for (uint16_t i = 0; i < N; i++) sitio[i] = i;
    sitio = bubble_sort_by(sitio,faixas,N);

    for (uint16_t i = N-1; i >=0 ; i--){
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
            double hs = hospitalizacao[faixas[site]];
            if(vacinado[site]) hs *= 0.034;

            if(prob_estagio[site]<hs) return phi;
            else return gamma_I;
            break;
        case 4: // Hospitalziado
            double mortalidade = morte[faixas[site]];
            if(vacinado[site]) mortalidade *= 0.034;

            if(prob_estagio[site]<mortalidade) return delta;
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

void infect(int S0,int E0,int N,int seed,double** infect_time,double* quant,double* final,double f){

    uint16_t Suscetiveis = S0;
    uint16_t Expostos = E0;
    uint16_t Assintomaticos = 0;
    uint16_t Sintomaticos = 0;
    uint16_t Hospitalizados = 0;
    uint16_t Recuperados = 0;
    uint16_t Mortos = 0;
    bool Vac = false;
    int site;
    
    //FILE *file;
    //char filename[40];

    //sprintf(filename,"./output/infect/%d.txt",E0);
    //file = fopen(filename,"a");

    struct Graph G;
    G = local_configuration_model(N,0,seed);
    int* faixas = get_faixas(G.Nodes);
    double* prob_estagio = (double*) malloc(G.Nodes*sizeof(double));
    double* rates = (double*) malloc(G.Nodes*sizeof(double));
    bool inicio = true;
    
    uint8_t* estagio = calloc(G.Nodes,sizeof(uint8_t));
    bool* vacinado = calloc(G.Nodes,sizeof(bool));
    //double* time = (double*) malloc(G.Nodes*sizeof(double));
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
    long double tempo = 0;
    long double ano = 0;
    double rate,Delta;
    double sintoma;
    double grau;
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
        if((ano>= 61) && (!Vac)){
            //vacinado = vacinacao_aleatoria(vacinado,Suscetiveis, Recuperados,f, estagio,G.Nodes);
            //vacinado = vacinacao_grau(vacinado,G,Suscetiveis, Recuperados,f, estagio);
            ///vacinado = vacinacao_idade(vacinado,faixas,Suscetiveis, Recuperados,f, estagio,G.Nodes);
            Vac = true;
        }

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
            new_rate += calc_estagio(i, estagio,prob_estagio,faixas,G,vacinado,hospitalizacao,morte);
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
                sintoma = sintomatico[faixas[i]];
                if(vacinado[i]) sintoma *= 0.34693877551;
                if(prob_estagio[i]<sintoma){
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

        grau = 0;
        for (i = 0; i < G.Nodes; i++) if((estagio[i] == 0) || (estagio[i] == 5))grau+= G.viz[i][0];
        
        infect_time[s - 1][7] += grau/(Suscetiveis+Recuperados);
        infect_time[s - 1][8] += infect_time[s - 1][7]*infect_time[s - 1][7];
        quant[s - 1]++;
        
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
    rodou++;
    //printf("\e[1;1H\e[2J");
    printf("%d\n",rodou);
}

void generate_file(char* filename,void* array,int linhas,int colunas,int check,double f){

    FILE *file;
    sprintf(filename,"%s_%.2f.txt",filename,f);
    file = fopen(filename,"w");
    for (uint16_t i = 0; i < linhas; i++){
        char print[100] = "";
        for(uint16_t j = 0;j<colunas;j++){
            switch (check){
                case sizeof(int)/* constant-expression */:
                    if(j != colunas-1) sprintf(print,"%s%d ",print,((int**)array)[i][j]);
                    else  sprintf(print,"%s%d\n",print,((int**)array)[i][j]);
                    break;
                case sizeof(double):
                    if(j != colunas-1) sprintf(print,"%s%f ",print,((double**)array)[i][j]);
                    else  sprintf(print,"%s%f\n",print,((double**)array)[i][j]);
                    break;
                default:
                    break;
                }
        }
        fprintf(file,"%s",print);
    }
    fclose(file);
}

void generate_infect(int seed, int redes,double f){
    rodou = 0;
    uint16_t N = size_txt();
    uint16_t tempo = 365*3*2;
    uint16_t i,j;

    //double* resultados = (double*) malloc(2*sizeof(double));
    double** infect_time = (double**) malloc(tempo*sizeof(double*));
    double* quant = (double*) malloc(tempo*sizeof(double));

    for (int i = 0; i < tempo; i++){
        infect_time[i] = (double*) malloc(9*sizeof(double));
        for (int j = 0; j < 9; j++) infect_time[i][j] = 0;
        quant[i] = 0;
    }

    double* final = (double*) malloc(4*sizeof(double));
    final[0] = 0;
    final[1] = 0;
    final[2] = 0;
    final[3] = 0;

    //#pragma omp parallel for
    for (j = 0; j < redes; j++)infect(N-1,1,N,seed+j ,infect_time,quant,final,f);
    
    for (i = 0; i < tempo; i++) for (int j = 0; j < 9; j++) infect_time[i][j] /= quant[i];
    for (i = 0; i < tempo; i++) infect_time[i][8] = pow(infect_time[i][8] - pow(infect_time[i][7],2),0.5);
    for(i = 0;i < 4;i++) final[i] /= redes;

    if(f==0.5){
        char filename[200] = "./output/infect/infect_vacinado_idade";
        generate_file(filename,infect_time,tempo,9,sizeof(infect_time[0][0]),f);
    }

    if (redes != 1){
        FILE *file;
        char filename[200];
        sprintf(filename,"./output/infect/infect_vacinado_idade.txt");
        file = fopen(filename,"a");
        fprintf(file,"%f %f %f %f %f\n",f, final[0], final[1],pow(final[2] - final[0]*final[0],0.5), pow(final[3] - final[1]*final[1],0.5));
        fclose(file);
    }


    for (i = 0; i < tempo; i++) free(infect_time[i]);
    free(infect_time);
    free(final);
    //free(resultados);
    free(quant);

    if (redes != 1) printf("\e[1;1H\e[2J");
    printf("Rodou: %f\n",f);
}