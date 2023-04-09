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

int rodou = 0;

void infect(int S0,int E0,int N,int seed,double** infect_time,double* quant,double* final,double f){

    int Suscetiveis = S0;
    int Expostos = E0;
    int Assintomaticos = 0;
    int Sintomaticos = 0;
    int Hospitalizados = 0;
    int Recuperados = 0;
    int Mortos = 0;
    int Vac = f*N;
    
    //FILE *file;
    //char filename[40];

    //sprintf(filename,"./output/infect/%d.txt",E0);
    //file = fopen(filename,"a");

    struct Graph G;
    G = local_configuration_model(N,0,seed);
    int* faixas = get_faixas(G.Nodes);
    double* random = (double*) malloc(G.Nodes*sizeof(double));
    int* estagio = (int*) malloc(G.Nodes*sizeof(int));
    int* vacinado = (int*) malloc(G.Nodes*sizeof(int));
    //double* time = (double*) malloc(G.Nodes*sizeof(double));
    init_genrand64(seed);


    int i = 0;
    for (i = 0; i < G.Nodes; i++){
        estagio[i] = 0;
        random[i] = genrand64_real1();
        vacinado[i] = 0;
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
    double t = 0.5;
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
    double ano = 0;
    //printf("%f %d %d %d %d %d %d %d\n",0.0,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos);
    while (s != 0){

        double rate = 0;
        if((ano>= 61) && (Vac !=0)){
            if(Vac < Suscetiveis + Recuperados){
                while(Vac!=0){
                    double r = genrand64_real1();
                    int N0 = r*G.Nodes;
                    if((estagio[N0] == 0) || (estagio[N0] == 5)){
                        vacinado[N0] = 1;
                        Vac--;
                    }
                }
            }
            else for(int k = 0;k < G.Nodes;k++)  if((estagio[k] == 0) || (estagio[k] == 5)) vacinado[k] = 1;
            Vac = 0;
        }


        for (i = 0; i < G.Nodes; i++){
            switch (estagio[i]){
                case 0:// Suscetível
                    double beta = 0;

                    for (int j = 0; j < G.viz[i][0]; j++){
                        int vizinho = G.viz[i][j+1];
                        if(estagio[vizinho] == 1){
                            if(vacinado[vizinho] == 0) beta += beta1;
                            else beta += 0.058*beta1;
                        }
                        if(estagio[vizinho] == 2){
                            if(vacinado[vizinho] == 0) beta += beta2;
                            else beta += 0.058*beta2;
                        }
                    }

                    if(vacinado[i] == 0) rate += beta;
                    else rate += 0.058*beta;
                    break;
                case 1: // Exposto
                    rate += sigma;
                    break;
                case 2: //Assintomático
                    rate += gamma_A;
                    break;
                case 3: // Sintomático
                    double hs = hospitalizacao[faixas[i]];
                    if(vacinado[i] == 1) hs *= 0.034;

                    if(random[i]<hs) rate += phi;
                    else rate += gamma_I;
                    break;
                case 4: // Hospitalziado
                    double mortalidade = morte[faixas[i]];
                    if(vacinado[i] == 1) mortalidade *= 0.034;

                    if(random[i]<mortalidade) rate += delta;
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
        ano += tempo;
        if(ano > 3*365) break;

        for (i = 0; i < G.Nodes; i++){
            
            switch (estagio[i]){
            case 0:// Suscetível
                double beta = 0;
                for (int j = 0; j < G.viz[i][0]; j++){
                    int vizinho = G.viz[i][j+1];
                    if(estagio[vizinho] == 1){
                        if(vacinado[vizinho] == 0) beta += beta1;
                        else beta += 0.058*beta1;
                    }
                    if(estagio[vizinho] == 2){
                        if(vacinado[vizinho] == 0) beta += beta2;
                        else beta += 0.058*beta2;
                    }
                }
                if(vacinado[i] == 0) rate += beta;
                else rate += 0.058*beta;
                break;
            case 1: // Exposto
                rate += sigma;
                break;
            case 2: //Assintomático
                rate += gamma_A;
                break;
            case 3: // Sintomático
                double hs = hospitalizacao[faixas[i]];
                if(vacinado[i] == 1) hs *= 0.034;

                if(random[i]<hs) rate += phi;
                else rate += gamma_I;
                break;
            case 4: // Hospitalziado
                double mortalidade = morte[faixas[i]];
                if(vacinado[i] == 1) mortalidade *= 0.034;

                if(random[i]<mortalidade) rate += delta;
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

        switch (estagio[i]){
            case 0:// Sucetível
                estagio[i] = 1;
                Expostos++;
                Suscetiveis--;
                break;
            case 1: // Exposto
                double sintoma = sintomatico[faixas[i]];
                if(vacinado[i] == 1) sintoma *= 0.34693877551;
                if(random[i]<sintoma){
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
        random[i] = genrand64_real1();
        //printf("%f %d %d %d %d %d %d %d\n",ano,Suscetiveis,Expostos,Assintomaticos,Sintomaticos,Hospitalizados,Recuperados,Mortos);
        if(ano > t) {s++;t += 0.5;}

        infect_time[s - 1][0] += Suscetiveis;
        infect_time[s - 1][1] += Expostos;
        infect_time[s - 1][2] += Assintomaticos;
        infect_time[s - 1][3] += Sintomaticos;
        infect_time[s - 1][4] += Hospitalizados;
        infect_time[s - 1][5] += Recuperados;
        infect_time[s - 1][6] += Mortos;
        quant[s - 1]++;
        if((Expostos == Assintomaticos) && (Assintomaticos == Sintomaticos) && (Sintomaticos ==Hospitalizados) && (Hospitalizados == 0)) break;

    }
    *final += Nhospitalizados;
    *(final+1) += Mortos;
    //fclose(file);
    free(morte);
    free(random);
    free(faixas);
    free(hospitalizacao);
    free(sintomatico);
    free(estagio);
    for(int j = 0; j < N; j++)free(G.viz[j]);
    free(G.viz);
    rodou++;
    printf("\e[1;1H\e[2J");
    printf("%d\n",rodou);
}

void generate_file(char* filename,void* array,int linhas,int colunas,int check,int x,double f){

    FILE *file;
    sprintf(filename,"%s_%d_%.2f.txt",filename,x,f);
    file = fopen(filename,"w");
    for (int i = 0; i < linhas; i++){
        char print[100] = "";
        for(int j = 0;j<colunas;j++){
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

    int N = size_txt();
    int tempo = 365*3*2;

    double* resultados = (double*) malloc(2*sizeof(double));
    double** infect_time = (double**) malloc(tempo*sizeof(double*));
    double* quant = (double*) malloc(tempo*sizeof(double));
    for (int i = 0; i < tempo; i++){
        infect_time[i] = (double*) malloc(7*sizeof(double));
        for (int j = 0; j < 7; j++) infect_time[i][j] = 0;
        quant[i] = 0;
    }

    double* final = (double*) malloc(2*sizeof(double));
    final[0] = 0;
    final[1] = 0;
    #pragma omp parallel for
    for (int j = 0; j < redes; j++) infect(N-1,1,N,seed+j,infect_time,quant,final,f);
    
    for (int i = 0; i < tempo; i++) for (int j = 0; j < 7; j++) infect_time[i][j] /= quant[i];
    
    char filename[40] = "./output/infect/infect";
    generate_file(filename,infect_time,tempo,7,sizeof(infect_time[0][0]),1,f);

    /* FILE *file;
    char filename[40];

    sprintf(filename,"./output/infect/infect_vacinado.txt");
    file = fopen(filename,"a");
    fprintf(file,"%f %f %f\n",f,(double) final[0]/redes,(double) final[1]/redes);
    fclose(file); */
}