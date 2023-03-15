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

struct iGraph{
    int* mortalidade;
    int* sintomatico;
    int* estagio;
    int** W;
    int* idade;
    double* random;
    struct Graph G;
};

int infect(struct Graph Grafo, int S0,int E0){

    struct iGraph rede;
    rede.G.Nodes = Grafo.Nodes;
    rede.G.edges = Grafo.edges;
    rede.G.viz = Grafo.viz;

    rede.random = (double*) malloc(rede.G.Nodes*sizeof(double));


    int i = 0;
    for (i = 0; i < rede.G.Nodes; i++){
        rede.estagio[i] = 0;
        rede.random[i] = genrand64_real1();
    }
    while(E0!=0){
        double rand = genrand64_real1();
        int N = rand*rede.G.Nodes;
        if(rede.estagio[N] == 0){
            rede.estagio[N] = 1;
            E0--;
        }
    }
    int s = 0;
    double sigma = 0;
    double gamma_A = 0;
    double hospitalizacao = 0;
    double phi = 0;
    double gamma_I = 0;
    double gamma_H = 0;
    double morte = 0;
    double delta = 0;
    double sintomatico = 0;
    while (s == 0){
        double rate = 0;
        for (i = 0; i < rede.G.Nodes; i++){
            switch (rede.estagio[i]){
            case 0:// Sucetível
                /* code */
                break;
            case 1: // Exposto
                rate += sigma;
                break;
            case 2: //Assintomático
                rate += gamma_A;
                break;
            case 3: // Sintomático
                double rand = genrand64_real1();
                if(rede.random[i]<hospitalizacao) rate += phi;
                else rate += gamma_I;
                break;
            case 4: // Hospitalziado
                double rand = genrand64_real1();
                if(rede.random[i]<morte) rate += delta;
                else rate += gamma_H;
                break;
            default:
                break;
            }
        }
        double tempo = exponentialRand(1/rate);
        double Delta = rate*genrand64_real1();
        rate = 0;
        for (i = 0; i < rede.G.Nodes; i++){
            switch (rede.estagio[i]){
            case 0:// Sucetível
                /* code */
                break;
            case 1: // Exposto
                rate += sigma;
                break;
            case 2: //Assintomático
                rate += gamma_A;
                break;
            case 3: // Sintomático
                double rand = genrand64_real1();
                if(rand<hospitalizacao) rate += phi;
                else rate += gamma_I;
                break;
            case 4: // Hospitalziado
                double rand = genrand64_real1();
                if(rand<morte) rate += delta;
                else rate += gamma_H;
                break;
            default:
                break;
            }
            if(rate>=Delta) break;
        }

        switch (rede.estagio[i]){
            case 0:// Sucetível
                rede.estagio[i] = 1;
                break;
            case 1: // Exposto
                if(rede.random[i]<sintomatico)rede.estagio[i] = 2;
                else rede.estagio[i] = 1;
                break;
            case 2: //Assintomático
                rede.estagio[i] = 5;
                break;
            case 3: // Sintomático
                if(rede.random[i]<hospitalizacao) rede.estagio[i] = 4;
                else rede.estagio[i] = 5;
                break;
            case 4: // Hospitalziado
                if(rede.random[i]<morte) rede.estagio[i] = 6;
                else rede.estagio[i] = 5;
                break;
            default:
                break;
        }
        rede.random[i] = genrand64_real1();

    }
    
    
}