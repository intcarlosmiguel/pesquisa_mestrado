#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include "mtwister.h"
#include "calc.h"
#include "rede.h"
#include "LCM.h"
#include "SBM.h"

void infect(int S0,int E0,int N,int seed,double* s1,double* s2){

    int Suscetiveis = S0;
    int Expostos = E0;
    int Assintomaticos = 0;
    int Sintomaticos = 0;
    int Hospitalizados = 0;
    int Recuperados = 0;
    int Mortos = 0;
    
    FILE *file;
    char filename[40];

    sprintf(filename,"./output/infect/%d.txt",E0);
    file = fopen(filename,"w");

    struct Graph G;
    G = local_configuration_model(N,0,seed);
    int* faixas = get_faixas(G.Nodes);
    double* random = (double*) malloc(G.Nodes*sizeof(double));
    int* estagio = (int*) malloc(G.Nodes*sizeof(int));
    //double* time = (double*) malloc(G.Nodes*sizeof(double));
    init_genrand64(seed);


    int i = 0;
    for (i = 0; i < G.Nodes; i++){
        estagio[i] = 0;
        random[i] = genrand64_real1();
        //time[i] = 0;
    }

    while(E0!=0){
        double r = genrand64_real1();
        int N0 = r*G.Nodes;
        if(estagio[N0] == 0){
            estagio[N0] = 1;
            E0--;
        }
    }
    int s = 1;

    double beta1 = (double)0.5;
    double beta2 = (double)0.41;
    double sigma = (double)1/5.1;
    double gamma_A = (double)1/7;
    double phi = (double)0.025;
    double gamma_I = (double)1/7;
    double gamma_H = (double)1/14;
    double delta = (double)1/14;
    double recupera = (double)1/40;

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


    int Nhospitalizados = 0;
    double tempo = 0;
    //printf("%f %d %d %d %d %d %d %d\n",0.0,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos);
    while (s != 0){
        double rate = 0;

        for (i = 0; i < G.Nodes; i++){
            switch (estagio[i]){
                case 0:// Suscetível
                    double beta = 0;
                    for (int j = 0; j < G.viz[i][0]; j++){
                        int vizinho = G.viz[i][j+1];
                        if(estagio[vizinho] == 1) beta += beta1;
                        if(estagio[vizinho] == 2) beta += beta2;
                    }
                    rate += beta;
                    break;
                case 1: // Exposto
                    rate += sigma;
                    break;
                case 2: //Assintomático
                    rate += gamma_A;
                    break;
                case 3: // Sintomático
                    if(random[i]<hospitalizacao[faixas[i]]) rate += phi;
                    else rate += gamma_I;
                    break;
                case 4: // Hospitalziado
                    if(random[i]<morte[faixas[i]]) rate += delta;
                    else rate += gamma_H;
                    break;
                case 5:
                    rate += recupera;
                    break;
                default:
                    break;
            }
        }

        if(rate==0) break;
        tempo = exponentialRand(rate);
        if(tempo ==0) break;
        double Delta = rate*genrand64_real1();
        rate = 0;

        for (i = 0; i < G.Nodes; i++){
            
            switch (estagio[i]){
            case 0:// Sucetível
                double beta = 0;
                for (int j = 0; j < G.viz[i][0]; j++){
                    int vizinho = G.viz[i][j+1];
                    if(estagio[vizinho] == 1) beta += beta1;
                    if(estagio[vizinho] == 2) beta += beta2;
                }
                rate += beta;
                break;
            case 1: // Exposto
                rate += sigma;
                break;
            case 2: //Assintomático
                rate += gamma_A;
                break;
            case 3: // Sintomático
                if(random[i]<hospitalizacao[faixas[i]]) rate += phi;
                else rate += gamma_I;
                break;
            case 4: // Hospitalziado
                if(random[i]<morte[faixas[i]]) rate += delta;
                else rate += gamma_H;
                break;
            case 5:
                rate += recupera;
                break;
            default:
                break;
            }
            if(rate>=Delta) break;
        }
        //printf("%f,%f,%d,%d\n",rate,Delta,site,estagio[site]);

        
        

        switch (estagio[i]){
            case 0:// Sucetível
                estagio[i] = 1;
                Expostos++;
                Suscetiveis--;
                break;
            case 1: // Exposto
                if(random[i]<sintomatico[faixas[i]]){
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
                if(random[i]<hospitalizacao[faixas[i]]){
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
                if(random[i]<morte[faixas[i]]){
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
        /*for (int j = 0; j < G.Nodes; j++){
            if(j!=i) time[j] += tempo;
            else{
                FILE *mean_tempo;
                char filename[40];
                sprintf(filename,"./output/time_%d.txt",estagio[i]);
                mean_tempo = fopen(filename,"a");
                fprintf(mean_tempo,"%f\n",time[i]);
                fclose(mean_tempo);
                time[i] = 0;
            }
        }*/
        random[i] = genrand64_real1();
        printf("%f %d %d %d %d %d %d %d\n",tempo,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos);
        fprintf(file,"%f %d %d %d %d %d %d %d\n",tempo,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos);
        //s++;
        //if(Expostos == 0) break;
    }
    fclose(file);
    *s1 += Nhospitalizados;
    *s2 += Mortos;
    free(morte);
    free(random);
    free(faixas);
    free(hospitalizacao);
    free(sintomatico);
    free(estagio);
    for(int j = 0; j < N; j++)free(G.viz[j]);
    free(G.viz);
}

void generate_infect(int seed, int redes){

    int N = size_txt();
    FILE *file;
    file = fopen("./output/infect.txt","w");
    double* resultados = (double*) malloc(2*sizeof(double));
    //#pragma omp parallel for
    /* for (int i = 1; i < N; i++){
        printf("\e[1;1H\e[2J");
        //resultados[0] = 0;
        //resultados[1] = 0;
        //#pragma omp parallel for
        for (int j = 0; j < redes; j++) infect(N-i,i,N,i*j+ seed,&resultados[0],&resultados[1]);
        //resultados[0] /= redes;
        //resultados[1] /= redes;
        //fprintf(file,"%f %f\n",resultados[0],resultados[1]);
        printf("%d\n",i);
    } */
    //#pragma omp parallel for
    //for (int j = 0; j < redes; j++) infect(N-1,1,N,j+ seed,&resultados[0],&resultados[1]);
    //free(resultados);
    infect(N-800,800,N,seed,&resultados[0],&resultados[1]);
    fclose(file);
}