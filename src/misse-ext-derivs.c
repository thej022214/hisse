/*
 *  misse-ext-derivs_c
 *
 *
 *  Created by Jeremy Beaulieu 9/25/2018
 *  Copyright 2018 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 54

static double params_misse[NUMELEMENTS];


void initmod_misse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_misse);
}


void misse_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E0A = y[0];
    double E0B = y[1];
    double E0C = y[2];
    double E0D = y[3];
    double E0E = y[4];
    double E0F = y[5];
    double E0G = y[6];
    double E0H = y[7];
    double E0I = y[8];
    double E0J = y[9];
    double E0K = y[10];
    double E0L = y[11];
    double E0M = y[12];
    double E0N = y[13];
    double E0O = y[14];
    double E0P = y[15];
    double E0Q = y[16];
    double E0R = y[17];
    double E0S = y[18];
    double E0T = y[19];
    double E0U = y[20];
    double E0V = y[21];
    double E0W = y[22];
    double E0X = y[23];
    double E0Y = y[24];
    double E0Z = y[25];
    
    double D0A = y[26];
    double D0B = y[27];
    double D0C = y[28];
    double D0D = y[29];
    double D0E = y[30];
    double D0F = y[31];
    double D0G = y[32];
    double D0H = y[33];
    double D0I = y[34];
    double D0J = y[35];
    double D0K = y[36];
    double D0L = y[37];
    double D0M = y[38];
    double D0N = y[39];
    double D0O = y[40];
    double D0P = y[41];
    double D0Q = y[42];
    double D0R = y[43];
    double D0S = y[44];
    double D0T = y[45];
    double D0U = y[46];
    double D0V = y[47];
    double D0W = y[48];
    double D0X = y[49];
    double D0Y = y[50];
    double D0Z = y[51];
    
    double
    lambda0A = params_misse[0],
    mu0A = params_misse[1],
    
    lambda0B = params_misse[2],
    mu0B = params_misse[3],
    
    lambda0C = params_misse[4],
    mu0C = params_misse[5],
    
    lambda0D = params_misse[6],
    mu0D = params_misse[7],
    
    lambda0E = params_misse[8],
    mu0E = params_misse[9],
    
    lambda0F = params_misse[10],
    mu0F = params_misse[11],
    
    lambda0G = params_misse[12],
    mu0G = params_misse[13],
    
    lambda0H = params_misse[14],
    mu0H = params_misse[15],
    
    lambda0I = params_misse[16],
    mu0I = params_misse[17],
    
    lambda0J = params_misse[18],
    mu0J = params_misse[19],
    
    lambda0K = params_misse[20],
    mu0K = params_misse[21],
    
    lambda0L = params_misse[22],
    mu0L = params_misse[23],
    
    lambda0M = params_misse[24],
    mu0M = params_misse[25],
    
    lambda0N = params_misse[26],
    mu0N = params_misse[27],
    
    lambda0O = params_misse[28],
    mu0O = params_misse[29],
    
    lambda0P = params_misse[30],
    mu0P = params_misse[31],
    
    lambda0Q = params_misse[32],
    mu0Q = params_misse[33],
    
    lambda0R = params_misse[34],
    mu0R = params_misse[35],
    
    lambda0S = params_misse[36],
    mu0S = params_misse[37],
    
    lambda0T = params_misse[38],
    mu0T = params_misse[39],
    
    lambda0U = params_misse[40],
    mu0U = params_misse[41],
    
    lambda0V = params_misse[42],
    mu0V = params_misse[43],
    
    lambda0W = params_misse[44],
    mu0W = params_misse[45],
    
    lambda0X = params_misse[46],
    mu0X = params_misse[47],
    
    lambda0Y = params_misse[48],
    mu0Y = params_misse[49],
    
    lambda0Z = params_misse[50],
    mu0Z = params_misse[51];
    
    int hidden_states = params_misse[53];
    double q0[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int i;
    for ( i = 0; i < 26; i++ ) {
        if(i < hidden_states){
            q0[i] = params_misse[52];
        }
    }

    
    ydot[0] = mu0A - (lambda0A + q0[0] + mu0A) * E0A + lambda0A*E0A*E0A + (q0[1]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[1] = mu0B - (lambda0B + q0[1] + mu0B) * E0B + lambda0B*E0B*E0B + (q0[0]*E0A + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[2] = mu0C - (lambda0C + q0[2] + mu0C) * E0C + lambda0C*E0C*E0C + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[3] = mu0D - (lambda0D + q0[2] + mu0D) * E0D + lambda0D*E0D*E0D + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[4] = mu0E - (lambda0E + q0[2] + mu0E) * E0E + lambda0E*E0E*E0E + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[5] = mu0F - (lambda0F + q0[2] + mu0F) * E0F + lambda0F*E0F*E0F + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[6] = mu0G - (lambda0G + q0[2] + mu0G) * E0G + lambda0G*E0G*E0G + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[7] = mu0H - (lambda0H + q0[2] + mu0H) * E0H + lambda0H*E0H*E0H + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[8] = mu0I - (lambda0I + q0[2] + mu0I) * E0I + lambda0I*E0I*E0I + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[9] = mu0J - (lambda0J + q0[2] + mu0J) * E0J + lambda0J*E0J*E0J + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[10] = mu0K - (lambda0K + q0[2] + mu0K) * E0K + lambda0K*E0K*E0K + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[11] = mu0L - (lambda0L + q0[2] + mu0L) * E0L + lambda0L*E0L*E0L + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[12] = mu0M - (lambda0M + q0[2] + mu0M) * E0M + lambda0M*E0M*E0M + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[13] = mu0N - (lambda0N + q0[2] + mu0N) * E0N + lambda0N*E0N*E0N + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[14] = mu0O - (lambda0O + q0[2] + mu0O) * E0O + lambda0O*E0O*E0O + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[15] = mu0P - (lambda0P + q0[2] + mu0P) * E0P + lambda0P*E0P*E0P + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[16] = mu0Q - (lambda0Q + q0[2] + mu0Q) * E0Q + lambda0Q*E0Q*E0Q + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[17] = mu0R - (lambda0R + q0[2] + mu0R) * E0R + lambda0R*E0R*E0R + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[18] = mu0S - (lambda0S + q0[2] + mu0S) * E0S + lambda0S*E0S*E0S + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[19] = mu0T - (lambda0T + q0[2] + mu0T) * E0T + lambda0T*E0T*E0T + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[20] = mu0U - (lambda0U + q0[2] + mu0U) * E0U + lambda0U*E0U*E0U + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[21] = mu0V - (lambda0V + q0[2] + mu0V) * E0V + lambda0V*E0V*E0V + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[22] = mu0W - (lambda0W + q0[2] + mu0W) * E0W + lambda0W*E0W*E0W + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0X + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[23] = mu0X - (lambda0X + q0[2] + mu0X) * E0X + lambda0X*E0X*E0X + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0Y + q0[2]*E0Z);
    
    ydot[24] = mu0Y - (lambda0Y + q0[2] + mu0Y) * E0Y + lambda0Y*E0Y*E0Y + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Z);
    
    ydot[25] = mu0Z - (lambda0Z + q0[2] + mu0Z) * E0Z + lambda0Z*E0Z*E0Z + (q0[2]*E0A + q0[2]*E0B + q0[2]*E0C + q0[2]*E0D + q0[2]*E0E + q0[2]*E0F + q0[2]*E0G + q0[2]*E0H + q0[2]*E0I + q0[2]*E0J + q0[2]*E0K + q0[2]*E0L + q0[2]*E0M + q0[2]*E0N + q0[2]*E0O + q0[2]*E0P + q0[2]*E0Q + q0[2]*E0R + q0[2]*E0S + q0[2]*E0T + q0[2]*E0U + q0[2]*E0V + q0[2]*E0W + q0[2]*E0X + q0[2]*E0Y);
    
    
    ydot[26] =  - (lambda0A + q0[0] + mu0A) * D0A + 2*lambda0A*E0A*D0A + (q0[1]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[27] =  - (lambda0B + q0[1] + mu0B) * D0B + 2*lambda0B*E0B*D0B + (q0[0]*D0A + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[28] =  - (lambda0C + q0[2] + mu0C) * D0C + 2*lambda0C*E0C*D0C + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[29] =  - (lambda0D + q0[2] + mu0D) * D0D + 2*lambda0D*E0D*D0D + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[30] =  - (lambda0E + q0[2] + mu0E) * D0E + 2*lambda0E*E0E*D0E + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[31] =  - (lambda0F + q0[2] + mu0F) * D0F + 2*lambda0F*E0F*D0F + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[32] =  - (lambda0G + q0[2] + mu0G) * D0G + 2*lambda0G*E0G*D0G + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[33] =  - (lambda0H + q0[2] + mu0H) * D0H + 2*lambda0H*E0H*D0H + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[34] =  - (lambda0I + q0[2] + mu0I) * D0I + 2*lambda0I*E0I*D0I + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[35] =  - (lambda0J + q0[2] + mu0J) * D0J + 2*lambda0J*E0J*D0J + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[36] =  - (lambda0K + q0[2] + mu0K) * D0K + 2*lambda0K*E0K*D0K + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[37] =  - (lambda0L + q0[2] + mu0L) * D0L + 2*lambda0L*E0L*D0L + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[38] =  - (lambda0M + q0[2] + mu0M) * D0M + 2*lambda0M*E0M*D0M + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[39] =  - (lambda0N + q0[2] + mu0N) * D0N + 2*lambda0N*E0N*D0N + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[40] =  - (lambda0O + q0[2] + mu0O) * D0O + 2*lambda0O*E0O*D0O + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[41] =  - (lambda0P + q0[2] + mu0P) * D0P + 2*lambda0P*E0P*D0P + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[42] =  - (lambda0Q + q0[2] + mu0Q) * D0Q + 2*lambda0Q*E0Q*D0Q + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[43] =  - (lambda0R + q0[2] + mu0R) * D0R + 2*lambda0R*E0R*D0R + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[44] =  - (lambda0S + q0[2] + mu0S) * D0S + 2*lambda0S*E0S*D0S + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[45] =  - (lambda0T + q0[2] + mu0T) * D0T + 2*lambda0T*E0T*D0T + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[46] =  - (lambda0U + q0[2] + mu0U) * D0U + 2*lambda0U*E0U*D0U + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[47] =  - (lambda0V + q0[2] + mu0V) * D0V + 2*lambda0V*E0V*D0V + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[48] =  - (lambda0W + q0[2] + mu0W) * D0W + 2*lambda0W*E0W*D0W + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0X + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[49] =  - (lambda0X + q0[2] + mu0X) * D0X + 2*lambda0X*E0X*D0X + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0Y + q0[2]*D0Z);
    
    ydot[50] =  - (lambda0Y + q0[2] + mu0Y) * D0Y + 2*lambda0Y*E0Y*D0Y + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Z);
    
    ydot[51] =  - (lambda0Z + q0[2] + mu0Z) * D0Z + 2*lambda0Z*E0Z*D0Z + (q0[2]*D0A + q0[2]*D0B + q0[2]*D0C + q0[2]*D0D + q0[2]*D0E + q0[2]*D0F + q0[2]*D0G + q0[2]*D0H + q0[2]*D0I + q0[2]*D0J + q0[2]*D0K + q0[2]*D0L + q0[2]*D0M + q0[2]*D0N + q0[2]*D0O + q0[2]*D0P + q0[2]*D0Q + q0[2]*D0R + q0[2]*D0S + q0[2]*D0T + q0[2]*D0U + q0[2]*D0V + q0[2]*D0W + q0[2]*D0X + q0[2]*D0Y);
    
}






