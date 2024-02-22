#pragma once
#include <stdint.h>
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
const uint16_t dias = 190;
const uint16_t dia_infecao = 100;
const uint16_t tempo_total = dias*2;
int q_resultados = 12;