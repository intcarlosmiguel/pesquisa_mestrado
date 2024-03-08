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

igraph_vector_int_t propose_centralities(igraph_t* Grafo,int estrategy,double* morte,double* hospitalizacao,double* sintomatico){
    uint16_t i;
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
    

    for ( i= 0; i < N; i++){

        igraph_vector_int_t vizinhos;
        igraph_vector_int_init(&vizinhos, 0);
        igraph_neighbors(Grafo, &vizinhos, i,IGRAPH_ALL);

        for ( int j = 0; j < igraph_vector_int_size(&vizinhos); j++){

            if(estrategy == 11) 
                VECTOR(prob)[i] = (double) VECTOR(degrees)[i]*morte[ (int) VECTOR(faixas)[i]];

            else{
                
                VECTOR(prob)[i] += sintomatico[(int) VECTOR(faixas)[VECTOR(vizinhos)[j]]];

                if(estrategy >=12) 
                    VECTOR(prob)[i] *= hospitalizacao[(int) VECTOR(faixas)[VECTOR(vizinhos)[j]]];
                if((estrategy == 13)||(estrategy == 15)) 
                    VECTOR(prob)[i] *= morte[(int) VECTOR(faixas)[VECTOR(vizinhos)[j]]];
                if(estrategy >= 14) 
                    VECTOR(prob)[i] *= (1 - sintomatico[(int) VECTOR(faixas)[i]]);
            }
        }
        igraph_vector_int_destroy(&vizinhos);
    }

    igraph_vector_qsort_ind(&prob,&centralidade, IGRAPH_DESCENDING);

    igraph_vector_destroy(&faixas);
    igraph_vector_destroy(&prob);

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
        }
        
        default:{
            printf("Métrica tradicional ponderada não encontrada!\n");
        }
    }
    return centralidade;
}

igraph_vector_int_t new_centralities(igraph_t* Grafo,int estrategy,igraph_vector_t* pesos,igraph_vector_int_t* edges){

    uint16_t i;
    int N = igraph_vcount(Grafo);
    igraph_vector_int_t centralidade;
    igraph_vector_int_init(&centralidade, N);
    
    switch (estrategy){
        
        case 21:{
            igraph_vector_t co_autorship;
            igraph_vector_init(&co_autorship, N);
            igraph_vector_add_constant(&co_autorship,1);

            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);

            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

            int site1,site2;
            double termo_duplo;

            for (i = 0; i < igraph_vector_int_size(edges); i+=2 ){
                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];
                termo_duplo = (double)pow(VECTOR(*pesos)[(int) i/2],2)/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];
                VECTOR(co_autorship)[site1] += (double)VECTOR(*pesos)[(int) i/2]/VECTOR(degrees)[site2] + termo_duplo;
                VECTOR(co_autorship)[site2] += (double)VECTOR(*pesos)[(int) i/2]/VECTOR(degrees)[site1] + termo_duplo;
            }
            igraph_vector_int_destroy(&degrees);
            igraph_vector_qsort_ind(&co_autorship,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&co_autorship);
            break;
        }
        case 22:{
            igraph_vector_t termo1;
            igraph_vector_init(&termo1, N);
            igraph_vector_t termo2;
            igraph_vector_init(&termo2, N);
            igraph_vector_t Laplacian_energy;
            igraph_vector_init(&termo1, N);
            igraph_vector_t termo4;
            igraph_vector_init(&termo2, N);

            int site1,site2;

            for (i = 0; i < igraph_vector_int_size(edges); i+=2 ){

                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];

                VECTOR(termo2)[site1] += pow(VECTOR(*pesos)[(int)i/2],2);
                VECTOR(termo1)[site1] += VECTOR(*pesos)[(int)i/2];
                VECTOR(termo2)[site2] += pow(VECTOR(*pesos)[(int)i/2],2);
                VECTOR(termo1)[site2] += VECTOR(*pesos)[(int)i/2];

                VECTOR(Laplacian_energy)[site1] -= pow(VECTOR(*pesos)[(int)i/2],2);
                VECTOR(termo4)[site1] -= VECTOR(*pesos)[(int)i/2];
                VECTOR(Laplacian_energy)[site2] -= pow(VECTOR(*pesos)[(int)i/2],2);
                VECTOR(termo4)[site2] -= VECTOR(*pesos)[(int)i/2];

            }
            double Normalized = 0;
            for (i = 0; i < N; i++ ){
                VECTOR(Laplacian_energy)[i] += VECTOR(termo2)[i];
                VECTOR(termo4)[i] += VECTOR(termo1)[i];
                VECTOR(termo1)[i] *= VECTOR(termo1)[i];
                VECTOR(termo4)[i] *= VECTOR(termo4)[i];
                Normalized += VECTOR(termo1)[i] + VECTOR(termo2)[i];
                VECTOR(Laplacian_energy)[i] += VECTOR(termo4)[i];
            }
            for (i = 0; i < N; i++ ) VECTOR(Laplacian_energy)[i] = Normalized - VECTOR(Laplacian_energy)[i];

            igraph_vector_destroy(&termo1);
            igraph_vector_destroy(&termo2);
            igraph_vector_destroy(&termo4);
            igraph_vector_qsort_ind(&Laplacian_energy,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&Laplacian_energy);
            break;
        }
        
        case 23:{
            
            igraph_vector_t gravity;
            
            int d = 1;
            igraph_vector_init(&gravity, N);
            igraph_vector_t faixas;
            igraph_vector_init(&faixas, 0);
            igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);

            igraph_matrix_t distance;
            igraph_matrix_init(&distance, 0, 0);

            igraph_distances(Grafo, &distance, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);

            for ( i = 0; i < N; i++){

                for (int j = i+1; i < N; i++){
                    VECTOR(gravity)[i] += VECTOR(faixas)[j]/pow(MATRIX(distance,i,j),d);
                    VECTOR(gravity)[j] += VECTOR(faixas)[i]/pow(MATRIX(distance,i,j),d);
                }
                VECTOR(gravity)[i] *= VECTOR(faixas)[i];
                
            }
            igraph_matrix_destroy(&distance);
            igraph_vector_qsort_ind(&gravity,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&gravity);

            break;
        }
        case 24:{
            
            igraph_vector_t gravity;
            
            int d = 1;
            igraph_vector_init(&gravity, N);
            igraph_vector_t faixas;
            igraph_vector_init(&faixas, 0);
            igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);

            igraph_matrix_t distance;
            igraph_matrix_init(&distance, 0, 0);

            igraph_distances_dijkstra(Grafo, &distance, igraph_vss_all(), igraph_vss_all(),pesos, IGRAPH_OUT);

            for ( i = 0; i < N; i++){

                for (int j = i+1; i < N; i++){
                    VECTOR(gravity)[i] += VECTOR(faixas)[j]/pow(MATRIX(distance,i,j),d);
                    VECTOR(gravity)[j] += VECTOR(faixas)[i]/pow(MATRIX(distance,i,j),d);
                }
                VECTOR(gravity)[i] *= VECTOR(faixas)[i];
                
            }
            igraph_matrix_destroy(&distance);
            igraph_vector_qsort_ind(&gravity,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&gravity);

            break;
        }
        case 25:{
            igraph_vector_t co_autorship;
            igraph_vector_init(&co_autorship, N);
            igraph_vector_add_constant(&co_autorship,1);

            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);

            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

            int site1,site2;
            double termo_duplo;

            for (i = 0; i < igraph_vector_int_size(edges); i+=2 ){
                site1 = VECTOR(*edges)[i+1];
                site2 = VECTOR(*edges)[i];
                termo_duplo = (double)1.0/(VECTOR(degrees)[site2])/VECTOR(degrees)[site1];
                VECTOR(co_autorship)[site1] += (double)1.0/VECTOR(degrees)[site2] + termo_duplo;
                VECTOR(co_autorship)[site2] += (double)1.0/VECTOR(degrees)[site1] + termo_duplo;
            }
            igraph_vector_int_destroy(&degrees);
            igraph_vector_qsort_ind(&co_autorship,&centralidade, IGRAPH_DESCENDING);
            igraph_vector_destroy(&co_autorship);
            break;
        }
        case 26:{

            igraph_vector_t gravity;
            igraph_vector_init(&gravity, N);
            double Normalized = 0;

            igraph_vector_t Energy_per_site;
            igraph_vector_init(&Energy_per_site, N);
            
            igraph_vector_t faixas;
            igraph_vector_init(&faixas, 0);
            igraph_cattribute_VANV(Grafo,"faixa",igraph_vss_all(),&faixas);

            igraph_matrix_t distances;
            igraph_matrix_init(&distances, 0, 0);
            igraph_distances(Grafo, &distances, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);

            igraph_vector_int_t degrees;
            igraph_vector_int_init(&degrees, 0);
            igraph_degree(Grafo, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

            int d = 1;

            for ( i = 0; i < N; i++){

                for (int j = i+1; i < N; i++){
                    VECTOR(gravity)[i] += VECTOR(faixas)[j]/pow(MATRIX(distances,i,j),d);
                    VECTOR(gravity)[j] += VECTOR(faixas)[i]/pow(MATRIX(distances,i,j),d);
                    VECTOR(Energy_per_site)[i] -= VECTOR(faixas)[j]/pow(MATRIX(distances,i,j),d);
                    VECTOR(Energy_per_site)[j] -= VECTOR(faixas)[i]/pow(MATRIX(distances,i,j),d);
                }
                VECTOR(gravity)[i] *= VECTOR(faixas)[i];
                Normalized += VECTOR(degrees)[i] + VECTOR(gravity)[i];
                VECTOR(Energy_per_site)[i] *= VECTOR(faixas)[i];
                VECTOR(Energy_per_site)[i] -= VECTOR(degrees)[i];
            }
            igraph_vector_add_constant(&Energy_per_site,Normalized);
            igraph_vector_qsort_ind(&Energy_per_site,&centralidade, IGRAPH_DESCENDING);

            igraph_matrix_destroy(&distances);
            igraph_vector_destroy(&gravity);
            igraph_vector_destroy(&Energy_per_site);
            igraph_vector_int_destroy(&degrees);
            break;
        }
        case 27:{
            
            igraph_vector_t efficiency;
            igraph_vector_init(&efficiency, N);

            igraph_matrix_t distance;
            igraph_matrix_init(&distance, 0, 0);
            igraph_distances(Grafo, &distance, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT);

            double Normalized = 0;
            for ( i = 0; i < N; i++){
                for ( int j = i+1; j < N; j++){
                    Normalized += 2.0/(MATRIX(distance,i,j));
                    VECTOR(efficiency)[i] -= 1.0/(MATRIX(distance,i,j));
                    VECTOR(efficiency)[j] -= 1.0/(MATRIX(distance,i,j));
                }
            }
            igraph_vector_add_constant(&efficiency,Normalized);
            igraph_vector_scale(&efficiency,1/Normalized);

            igraph_vector_qsort_ind(&efficiency,&centralidade, IGRAPH_ASCENDING);
            igraph_vector_destroy(&efficiency);
            break;
        }
        case 28:{
            
            igraph_vector_t efficiency;
            igraph_vector_init(&efficiency, N);

            igraph_matrix_t distance;
            igraph_matrix_init(&distance, 0, 0);
            igraph_distances_dijkstra(Grafo, &distance, igraph_vss_all(), igraph_vss_all(),pesos, IGRAPH_OUT);

            double Normalized = 0;
            for ( i = 0; i < N; i++){
                for ( int j = i+1; j < N; j++){
                    Normalized += 2.0/(MATRIX(distance,i,j));
                    VECTOR(efficiency)[i] -= 1.0/(MATRIX(distance,i,j));
                    VECTOR(efficiency)[j] -= 1.0/(MATRIX(distance,i,j));
                }
            }
            Normalized /= N*(N-1);
            igraph_vector_add_constant(&efficiency,Normalized);
            igraph_vector_scale(&efficiency,1/Normalized);

            igraph_vector_qsort_ind(&efficiency,&centralidade, IGRAPH_ASCENDING);
            igraph_vector_destroy(&efficiency);
            break;
        }
        default:{
            printf("Métrica tradicional ponderada não encontrada!\n");
        }
    }
    return centralidade;
}
