#pragma once
#include "mtwister.h"
#include "calc.h"
#include <igraph.h>
#include <math.h>

igraph_vector_int_t traditional_centralities(igraph_t* Grafo,int estrategy){
    int N = igraph_vcount(Grafo);
    igraph_vector_int_t centralidade;
    igraph_vector_int_init(&centralidade, N);
    switch (estrategy){
        case 0:{ // Vacinação por idade
            
            igraph_vector_t faixas;
            igraph_vector_init(&faixas, 0);
            igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);
            igraph_vector_qsort_ind(&faixas,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&faixas);
            break;
        }
        case 1:{ // Vacinação por grau
            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);
            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
            igraph_vector_int_qsort_ind(&degrees,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_int_destroy(&degrees);
            break;
        }
        case 2:{
            
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
        }
        case 3:{ // Vacinação por Harmo
            igraph_vector_t Harmonic;
            igraph_vector_init(&Harmonic, N);
            igraph_harmonic_centrality(Grafo,&Harmonic,igraph_vss_all(),IGRAPH_ALL,NULL,0);
            igraph_vector_qsort_ind(&Harmonic,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Harmonic);
            break;
        }
        case 4:{
            igraph_vector_t Betweenness;
            igraph_vector_init(&Betweenness, N);
            igraph_betweenness(Grafo, &Betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL);
            igraph_vector_qsort_ind(&Betweenness,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Betweenness);
            break;
        }
        case 5:{
            igraph_vector_t eigenvector;
            igraph_vector_init(&eigenvector, 0);
            igraph_eigenvector_centrality(Grafo,&eigenvector,0,IGRAPH_UNDIRECTED,1,NULL,NULL);
            igraph_vector_qsort_ind(&eigenvector,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&eigenvector);
            break;
        }
        case 6:{
            igraph_vector_t eccentricity;
            igraph_vector_init(&eccentricity, N);
            igraph_eccentricity(Grafo,&eccentricity,igraph_vss_all(),IGRAPH_ALL);
            igraph_vector_qsort_ind(&eccentricity,&centralidade, IGRAPH_ASCENDING);
            igraph_vector_destroy(&eccentricity);
            break;
        }
        case 7:{ // Clsutering
            igraph_vector_t clustering;
            igraph_vector_init(&clustering, N);
            igraph_transitivity_local_undirected(Grafo,&clustering,igraph_vss_all(),IGRAPH_TRANSITIVITY_ZERO);
            igraph_vector_qsort_ind(&clustering,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&clustering);
            break;
        }
        case 8:{ // Kshell
            igraph_vector_int_t k_shell;
            igraph_vector_int_init(&k_shell, N);
            igraph_coreness(Grafo, &k_shell, IGRAPH_ALL);
            igraph_vector_int_qsort_ind(&k_shell,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_int_destroy(&k_shell);
            break;
        }
        case 9:{
            igraph_vector_int_init_range(&centralidade,0,N);
            igraph_vector_int_shuffle(&centralidade);
            break;
        }
        case 10:{
            igraph_vector_t pagerank;
            igraph_real_t value;
            igraph_real_t d = 0.85;  // Fator de amortecimento, tipicamente 0.85
            igraph_pagerank_algo_t algo = IGRAPH_PAGERANK_ALGO_PRPACK;  // Algoritmo PageRank a ser usado
            
            igraph_vector_init(&pagerank, 0);

            // Calcula o PageRank
            igraph_pagerank(Grafo, algo,&pagerank, &value,igraph_vss_all(), 0,d,NULL,NULL );
            igraph_vector_qsort_ind(&pagerank,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&pagerank);
            break;
        }
        default:{
            printf("Métrica tradicional não encontrada!\n");
        }
    }
    return centralidade;
}

igraph_vector_int_t propose_centralities(igraph_t* Grafo,int estrategy){
    int i,j,vizinho;
    int N = igraph_vcount(Grafo);
    igraph_vector_int_t centralidade;
    igraph_vector_int_init(&centralidade, N);

    igraph_vector_t faixas;
    igraph_vector_init(&faixas, 0);
    igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);

    igraph_vector_t prob;
    igraph_vector_init(&prob, N);

    igraph_vector_int_t degrees;
    igraph_vector_int_init(&degrees, 0);
    igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    
    double* sintomatico = (double*) malloc(5*sizeof(double));
    double* hospitalizacao = (double*) malloc(5*sizeof(double));
    double* morte = (double*) malloc(5*sizeof(double));

    morte[0] = (double)0.07710809281267686;
    morte[1] = (double)0.09561476547321676;
    morte[2] = (double)0.13779231777947765;
    morte[3] = (double)0.30044562859568047;
    morte[4] = (double)0.530315172817809;

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

    double w;
    for ( i= 0; i < N; i++){

        igraph_vector_int_t vizinhos;
        igraph_vector_int_init(&vizinhos, 0);
        igraph_neighbors(Grafo, &vizinhos, i,IGRAPH_ALL);
        VECTOR(prob)[i] = (double) VECTOR(degrees)[i]*morte[ (int) VECTOR(faixas)[i]];
        for (j = 0; j < igraph_vector_int_size(&vizinhos); j++){
            // Grau morte
            if(estrategy == 11) break;
            else{
                
                vizinho = (int) VECTOR(vizinhos)[j];
                w = sintomatico[(int) VECTOR(faixas)[vizinho]];
                //VECTOR(prob)[i] += ;

                if(estrategy >=12) // Probabilidade dos vizinhos serem hospitalizados;
                    w *= hospitalizacao[(int) VECTOR(faixas)[vizinho]];
                if((estrategy == 13)||(estrategy == 15)) // Probabilidade do vizinho morrer & Probabilidade dos vizinhos morrerem e o vértice ser assintomático;
                    w *= morte[(int) VECTOR(faixas)[vizinho]];
                if(estrategy >= 14) //Probabilidade dos vizinhos serem hospitalizados e o vértice ser assintomático;
                    w *= (1 - sintomatico[(int) VECTOR(faixas)[i]]);
                VECTOR(prob)[i] += w;
            }
        }
        igraph_vector_int_destroy(&vizinhos);
    }

    igraph_vector_qsort_ind(&prob,&centralidade, IGRAPH_DESCENDING);

    igraph_vector_destroy(&faixas);
    igraph_vector_destroy(&prob);
    free(morte);
    free(hospitalizacao);
    free(sintomatico);
    return centralidade;
}

igraph_vector_int_t weighted_traditional_centralities(igraph_t* Grafo,int estrategy,igraph_vector_t* pesos){
    int N = igraph_vcount(Grafo);
    igraph_vector_int_t centralidade;
    igraph_vector_int_init(&centralidade, N);
    switch (estrategy){
        case 16:{
            igraph_vector_t closeness;
            igraph_vector_int_t reachable_count;
            igraph_bool_t all_reachable;

            igraph_vector_int_init(&reachable_count, 0);
            igraph_vector_init(&closeness, N);

            igraph_closeness(Grafo, &closeness, &reachable_count, &all_reachable, igraph_vss_all(), IGRAPH_ALL, pesos, 1);
            igraph_vector_qsort_ind(&closeness,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&closeness);
            igraph_vector_int_destroy(&reachable_count);
            
            break;
        }
        case 17:{
            igraph_vector_t Harmonic;
            igraph_vector_init(&Harmonic, N);
            igraph_harmonic_centrality(Grafo,&Harmonic,igraph_vss_all(),IGRAPH_ALL,pesos,0);
            igraph_vector_qsort_ind(&Harmonic,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Harmonic);
            break;
        }
        case 18:{
            igraph_vector_t Betweenness;
            igraph_vector_init(&Betweenness, N);
            igraph_betweenness(Grafo, &Betweenness, igraph_vss_all(), IGRAPH_UNDIRECTED, pesos);
            igraph_vector_qsort_ind(&Betweenness,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Betweenness);
            break;
        }
        case 19:{
            igraph_vector_t eigenvector;
            igraph_vector_init(&eigenvector, 0);
            igraph_eigenvector_centrality(Grafo,&eigenvector,0,IGRAPH_UNDIRECTED,1,pesos,NULL);
            igraph_vector_qsort_ind(&eigenvector,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&eigenvector);
            break;
        }
        case 20:{
            igraph_vector_t pagerank;
            igraph_real_t value;
            igraph_real_t d = 0.85;  // Fator de amortecimento, tipicamente 0.85
            igraph_pagerank_algo_t algo = IGRAPH_PAGERANK_ALGO_PRPACK;  // Algoritmo PageRank a ser usado

            igraph_vector_init(&pagerank, 0);

            // Calcula o PageRank
            igraph_pagerank(Grafo, algo,&pagerank, &value,igraph_vss_all(), IGRAPH_DIRECTED,d,pesos,NULL );
            igraph_vector_qsort_ind(&pagerank,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&pagerank);
            break;
        }
        
        default:{
            printf("Métrica tradicional ponderada não encontrada!\n");
        }
    }
    return centralidade;
}

igraph_vector_int_t new_centralities(igraph_t* Grafo,int estrategy,igraph_vector_t* pesos,igraph_vector_int_t* edges,igraph_vector_t* faixas){

    uint32_t i;
    int N = igraph_vcount(Grafo);
    igraph_vector_int_t centralidade;
    igraph_vector_int_init(&centralidade, N);
    double* sintomatico = (double*) malloc(5*sizeof(double));
    double* morte = (double*) malloc(5*sizeof(double));;
    sintomatico[0] = 29.1/100;
    sintomatico[1] = 37.4/100;
    sintomatico[2] = 41.68/100;
    sintomatico[3] = 39.4/100;
    sintomatico[4] = 31.3/100;
    morte[0] = (double)0.07710809281267686;
    morte[1] = (double)0.09561476547321676;
    morte[2] = (double)0.13779231777947765;
    morte[3] = (double)0.30044562859568047;
    morte[4] = (double)0.530315172817809;
    switch (estrategy){
        
        case 21:{

            igraph_vector_t utility;
            igraph_vector_init(&utility, N);

            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);

            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

            int site1,site2;
            double termo_duplo,peso,w_i,w_j;
            int L = igraph_vector_int_size(edges);
            for (i = 0; i < L; i+=2 ){
                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];

                peso = 1;
                //VECTOR(utility)[site1] += (1 - sintomatico[(int) VECTOR(*faixas)[site2]])*peso/VECTOR(degrees)[site2] + morte[(int) VECTOR(*faixas)[site1]]*(1 - sintomatico[(int) VECTOR(*faixas)[site2]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];

                VECTOR(utility)[site1] += morte[(int) VECTOR(*faixas)[site2]]*peso/VECTOR(degrees)[site2] + morte[(int) VECTOR(*faixas)[site2]]*(1 - sintomatico[(int) VECTOR(*faixas)[site1]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];

                //VECTOR(utility)[site2] += (1 - sintomatico[(int) VECTOR(*faixas)[site1]])*peso/VECTOR(degrees)[site1] + morte[(int) VECTOR(*faixas)[site2]]*(1 - sintomatico[(int) VECTOR(*faixas)[site1]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];

                VECTOR(utility)[site2] += morte[(int) VECTOR(*faixas)[site1]]*peso/VECTOR(degrees)[site1] + morte[(int) VECTOR(*faixas)[site1]]*(1 - sintomatico[(int) VECTOR(*faixas)[site2]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];


            }
            for (i = 0; i < N; i++){
                w_i = 1 - sintomatico[(int) VECTOR(*faixas)[i]];
                //w_i = morte[(int) VECTOR(*faixas)[i]];
                VECTOR(utility)[i] += w_i;
            }
            igraph_vector_int_destroy(&degrees);
            igraph_vector_qsort_ind(&utility,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&utility);
            break;
        }
        case 22:{
            igraph_vector_t utility;
            igraph_vector_init(&utility, N);   

            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);

            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

            int site1,site2;
            double termo_duplo,peso,w_i,w_j;
            int L = igraph_vector_int_size(edges);
            for (i = 0; i < L; i+=2 ){
                peso = VECTOR(*pesos)[(int)i/2];

                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];

                //VECTOR(utility)[site1] += (1 - sintomatico[(int) VECTOR(*faixas)[site2]])*peso/VECTOR(degrees)[site2] + morte[(int) VECTOR(*faixas)[site1]]*(1 - sintomatico[(int) VECTOR(*faixas)[site2]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];

                VECTOR(utility)[site1] += morte[(int) VECTOR(*faixas)[site2]]*peso/VECTOR(degrees)[site2] + morte[(int) VECTOR(*faixas)[site2]]*(1 - sintomatico[(int) VECTOR(*faixas)[site1]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];

                //VECTOR(utility)[site2] += (1 - sintomatico[(int) VECTOR(*faixas)[site1]])*peso/VECTOR(degrees)[site1] + morte[(int) VECTOR(*faixas)[site2]]*(1 - sintomatico[(int) VECTOR(*faixas)[site1]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];

                VECTOR(utility)[site2] += morte[(int) VECTOR(*faixas)[site1]]*peso/VECTOR(degrees)[site1] + morte[(int) VECTOR(*faixas)[site1]]*(1 - sintomatico[(int) VECTOR(*faixas)[site2]])*pow(peso,2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];
            }
           for (i = 0; i < N; i++){
                w_i = 1 - sintomatico[(int) VECTOR(*faixas)[i]];
                //w_i = morte[(int) VECTOR(*faixas)[i]];
                VECTOR(utility)[i] += w_i;
            }
            igraph_vector_int_destroy(&degrees);
            igraph_vector_qsort_ind(&utility,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&utility);
            break;
        }
        case 23:{
            igraph_vector_t Laplacian_energy;
            igraph_vector_init(&Laplacian_energy, N);

            int site1,site2;
            uint32_t j;
            double peso,peso_ij,peso_ji,Normalized = 0;
            double* x_i = (double*) calloc(N,sizeof(double));
            double** X = (double**) malloc(N*sizeof(double*));

            for (i = 0; i < N; i++ ) X[i] =  (double*) calloc(N,sizeof(double));
            for (i = 0; i < igraph_vector_int_size(edges); i+=2 ){
                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];
                //peso_ij = 0.5*(1 - sintomatico[(int) VECTOR(*faixas)[site1]] +  morte[(int) VECTOR(*faixas)[site2]]);
                //peso_ij = 0.5*(1 - sintomatico[(int) VECTOR(*faixas)[site2]] +  morte[(int) VECTOR(*faixas)[site1]]);
                peso_ij = 0.5*(1 - sintomatico[(int) VECTOR(*faixas)[site2]] +  morte[(int) VECTOR(*faixas)[site1]]);
                peso_ij = 0.5*(1 - sintomatico[(int) VECTOR(*faixas)[site1]] +  morte[(int) VECTOR(*faixas)[site2]]);

                x_i[site1] += peso_ij;
                x_i[site2] += peso_ji;

                X[site1][site2] = 1;
                X[site2][site1] = 1;
            }
            for (i = 0; i < N; i++ ) {
                for (j = 0; j < N; j++ ) if(j!=i) VECTOR(Laplacian_energy)[i] += 2*X[j][i]*x_i[j]+ pow(X[j][i],2);
                free(X[i]);
            }
            free(x_i);
            free(X);

            igraph_vector_qsort_ind(&Laplacian_energy,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Laplacian_energy);
            break;
        }
        case 24:{
            igraph_vector_t Laplacian_energy;
            igraph_vector_init(&Laplacian_energy, N);

            int site1,site2;
            uint32_t j;
            double peso,peso_ij,peso_ji;
            double* x_i = (double*) calloc(N,sizeof(double));
            double** X = (double**) malloc(N*sizeof(double*));

            for (i = 0; i < N; i++ ) X[i] =  (double*) calloc(N,sizeof(double));

            for (i = 0; i < igraph_vector_int_size(edges); i+=2 ){

                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];

                //peso_ij = 0.5*VECTOR(*pesos)[(int)i/2]*(1 - sintomatico[(int) VECTOR(*faixas)[site2]] +  morte[(int) VECTOR(*faixas)[site1]]);
                //peso_ij = 0.5*VECTOR(*pesos)[(int)i/2]*(1 - sintomatico[(int) VECTOR(*faixas)[site1]] +  morte[(int) VECTOR(*faixas)[site2]]);
                peso_ij = 0.5*VECTOR(*pesos)[(int)i/2]*(1 - sintomatico[(int) VECTOR(*faixas)[site2]] +  morte[(int) VECTOR(*faixas)[site1]]);
                peso_ij = 0.5*VECTOR(*pesos)[(int)i/2]*(1 - sintomatico[(int) VECTOR(*faixas)[site1]] +  morte[(int) VECTOR(*faixas)[site2]]);
                
                x_i[site1] += peso_ij;
                x_i[site2] += peso_ji;

                X[site1][site2] = VECTOR(*pesos)[(int)i/2];
                X[site2][site1] = VECTOR(*pesos)[(int)i/2];
            }

            for (i = 0; i < N; i++ ) {
                for (j = 0; j < N; j++ ) VECTOR(Laplacian_energy)[i] += 2*X[j][i]*x_i[j]+ pow(X[j][i],2);
                free(X[i]);
            }
            free(x_i);
            free(X);
            igraph_vector_qsort_ind(&Laplacian_energy,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Laplacian_energy);   
            break;
        }
        case 25:{
            
            igraph_vector_t gravity;
            
            double d = 1/2;
            igraph_vector_init(&gravity, N);
            igraph_matrix_t distance;
            igraph_matrix_init(&distance, 0, 0);
            igraph_distances(Grafo, &distance, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);
            double w_i,w_j;
            for ( i = 0; i < N; i++){
               // w_i = morte[(int) VECTOR(*faixas)[i]];
                w_i = (1 - sintomatico[(int) VECTOR(*faixas)[i]]);
                for (int j = i+1; j < N; j++){
                    double increment = 1.0/ pow(MATRIX(distance,i,j),d);
                    //VECTOR(gravity)[i] += (1 - sintomatico[(int) VECTOR(*faixas)[j]])*increment;
                    //VECTOR(gravity)[j] += (1 - sintomatico[(int) VECTOR(*faixas)[i]])*increment;
                    VECTOR(gravity)[i] += morte[(int) VECTOR(*faixas)[j]]*increment;
                    VECTOR(gravity)[j] += morte[(int) VECTOR(*faixas)[i]]*increment;
                }
                VECTOR(gravity)[i] *= w_i;
            }
            igraph_matrix_destroy(&distance);
            igraph_vector_qsort_ind(&gravity,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&gravity);

            break;
        }
        case 26:{   

            
            igraph_vector_t gravity;
            
            double d = 1/2;
            igraph_vector_init(&gravity, N);
            igraph_matrix_t distance;
            igraph_matrix_init(&distance, 0, 0);
            igraph_distances_dijkstra(Grafo, &distance, igraph_vss_all(), igraph_vss_all(),pesos, IGRAPH_OUT);
            double w_i,w_j;
            for ( i = 0; i < N; i++){
                //w_i = morte[(int) VECTOR(*faixas)[i]];
                w_i = (1 - sintomatico[(int) VECTOR(*faixas)[i]]);
                for (int j = i+1; j < N; j++){
                    double increment = 1.0/ pow(MATRIX(distance,i,j),d);
                    //VECTOR(gravity)[i] += (1 - sintomatico[(int) VECTOR(*faixas)[j]])*increment;
                    //VECTOR(gravity)[j] += (1 - sintomatico[(int) VECTOR(*faixas)[i]])*increment;
                    VECTOR(gravity)[i] += morte[(int) VECTOR(*faixas)[j]]*increment;
                    VECTOR(gravity)[j] += morte[(int) VECTOR(*faixas)[i]]*increment;
                }
                VECTOR(gravity)[i] *= w_i;
            }
            igraph_matrix_destroy(&distance);
            igraph_vector_qsort_ind(&gravity,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&gravity);

            break;
        }
        
        default:{
            printf("Métrica tradicional ponderada não encontrada!\n");
        }
    }
    free(morte);
    free(sintomatico);
    return centralidade;
}
