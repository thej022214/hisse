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
#define NUMELEMENTS 55

static double params_misse[NUMELEMENTS];


void initmod_misse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_misse);
}


void change_rate(double *q0, int length, int hidden_states, int index, double rate){
    
    if(index < hidden_states){
        int i;
        for ( i = 0; i < length; i++ ) {
            if(i < hidden_states){
                q0[i] = rate;
            }
        }
    }
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
    mu0Z = params_misse[51],
    
    rate = params_misse[52];
    
    int hidden_states = params_misse[53];
    
    double qA[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qA, 26, hidden_states, 0, rate);
    double qB[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qB, 26, hidden_states, 1, rate);
    double qC[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qC, 26, hidden_states, 2, rate);
    double qD[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qD, 26, hidden_states, 3, rate);
    double qE[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qE, 26, hidden_states, 4, rate);
    double qF[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qF, 26, hidden_states, 5, rate);
    double qG[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qG, 26, hidden_states, 6, rate);
    double qH[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qH, 26, hidden_states, 7, rate);
    double qI[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qI, 26, hidden_states, 8, rate);
    double qJ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qJ, 26, hidden_states, 9, rate);
    double qK[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qK, 26, hidden_states, 10, rate);
    double qL[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qL, 26, hidden_states, 11, rate);
    double qM[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qM, 26, hidden_states, 12, rate);
    double qN[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qN, 26, hidden_states, 13, rate);
    double qO[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qO, 26, hidden_states, 14, rate);
    double qP[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qP, 26, hidden_states, 15, rate);
    double qQ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qQ, 26, hidden_states, 16, rate);
    double qR[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qR, 26, hidden_states, 17, rate);
    double qS[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qS, 26, hidden_states, 18, rate);
    double qT[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qT, 26, hidden_states, 19, rate);
    double qU[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qU, 26, hidden_states, 20, rate);
    double qV[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qV, 26, hidden_states, 21, rate);
    double qW[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qW, 26, hidden_states, 22, rate);
    double qX[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qX, 26, hidden_states, 23, rate);
    double qY[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qY, 26, hidden_states, 24, rate);
    double qZ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qZ, 26, hidden_states, 25, rate);
    
    double psi = params_misse[54];
    
    
    /* The E's */
    ydot[0] = mu0A - (lambda0A + psi + (qA[1] + qA[2] + qA[3] + qA[4] + qA[5] + qA[6] + qA[7] + qA[8] + qA[9] + qA[10] + qA[11] + qA[12] + qA[13] + qA[14] + qA[15] + qA[16] + qA[17] + qA[18] + qA[19] + qA[20] + qA[21] + qA[22] + qA[23] + qA[24] + qA[25]) + mu0A) * E0A + lambda0A*E0A*E0A + (qA[1]*E0B + qA[2]*E0C + qA[3]*E0D + qA[4]*E0E + qA[5]*E0F + qA[6]*E0G + qA[7]*E0H + qA[8]*E0I + qA[9]*E0J + qA[10]*E0K + qA[11]*E0L + qA[12]*E0M + qA[13]*E0N + qA[14]*E0O + qA[15]*E0P + qA[16]*E0Q + qA[17]*E0R + qA[18]*E0S + qA[19]*E0T + qA[20]*E0U + qA[21]*E0V + qA[22]*E0W + qA[23]*E0X + qA[24]*E0Y + qA[25]*E0Z);
    
    ydot[1] = mu0B - (lambda0B + psi + (qB[0] + qB[2] + qB[3] + qB[4] + qB[5] + qB[6] + qB[7] + qB[8] + qB[9] + qB[10] + qB[11] + qB[12] + qB[13] + qB[14] + qB[15] + qB[16] + qB[17] + qB[18] + qB[19] + qB[20] + qB[21] + qB[22] + qB[23] + qB[24] + qB[25]) + mu0B) * E0B + lambda0B*E0B*E0B + (qB[0]*E0A + qB[2]*E0C + qB[3]*E0D + qB[4]*E0E + qB[5]*E0F + qB[6]*E0G + qB[7]*E0H + qB[8]*E0I + qB[9]*E0J + qB[10]*E0K + qB[11]*E0L + qB[12]*E0M + qB[13]*E0N + qB[14]*E0O + qB[15]*E0P + qB[16]*E0Q + qB[17]*E0R + qB[18]*E0S + qB[19]*E0T + qB[20]*E0U + qB[21]*E0V + qB[22]*E0W + qB[23]*E0X + qB[24]*E0Y + qB[25]*E0Z);
    
    ydot[2] = mu0C - (lambda0C + psi + (qC[0] + qC[1] + qC[3] + qC[4] + qC[5] + qC[6] + qC[7] + qC[8] + qC[9] + qC[10] + qC[11] + qC[12] + qC[13] + qC[14] + qC[15] + qC[16] + qC[17] + qC[18] + qC[19] + qC[20] + qC[21] + qC[22] + qC[23] + qC[24] + qC[25]) + mu0C) * E0C + lambda0C*E0C*E0C + (qC[0]*E0A + qC[1]*E0B + qC[3]*E0D + qC[4]*E0E + qC[5]*E0F + qC[6]*E0G + qC[7]*E0H + qC[8]*E0I + qC[9]*E0J + qC[10]*E0K + qC[11]*E0L + qC[12]*E0M + qC[13]*E0N + qC[14]*E0O + qC[15]*E0P + qC[16]*E0Q + qC[17]*E0R + qC[18]*E0S + qC[19]*E0T + qC[20]*E0U + qC[21]*E0V + qC[22]*E0W + qC[23]*E0X + qC[24]*E0Y + qC[25]*E0Z);
    
    ydot[3] = mu0D - (lambda0D + psi + (qD[0] + qD[1] + qD[2] + qD[4] + qD[5] + qD[6] + qD[7] + qD[8] + qD[9] + qD[10] + qD[11] + qD[12] + qD[13] + qD[14] + qD[15] + qD[16] + qD[17] + qD[18] + qD[19] + qD[20] + qD[21] + qD[22] + qD[23] + qD[24] + qD[25]) + mu0D) * E0D + lambda0D*E0D*E0D + (qD[0]*E0A + qD[1]*E0B + qD[2]*E0C + qD[4]*E0E + qD[5]*E0F + qD[6]*E0G + qD[7]*E0H + qD[8]*E0I + qD[9]*E0J + qD[10]*E0K + qD[11]*E0L + qD[12]*E0M + qD[13]*E0N + qD[14]*E0O + qD[15]*E0P + qD[16]*E0Q + qD[17]*E0R + qD[18]*E0S + qD[19]*E0T + qD[20]*E0U + qD[21]*E0V + qD[22]*E0W + qD[23]*E0X + qD[24]*E0Y + qD[25]*E0Z);
    
    ydot[4] = mu0E - (lambda0E + psi + (qE[0] + qE[1] + qE[2] + qE[3] + qE[5] + qE[6] + qE[7] + qE[8] + qE[9] + qE[10] + qE[11] + qE[12] + qE[13] + qE[14] + qE[15] + qE[16] + qE[17] + qE[18] + qE[19] + qE[20] + qE[21] + qE[22] + qE[23] + qE[24] + qE[25]) + mu0E) * E0E + lambda0E*E0E*E0E + (qE[0]*E0A + qE[1]*E0B + qE[2]*E0C + qE[3]*E0D + qE[5]*E0F + qE[6]*E0G + qE[7]*E0H + qE[8]*E0I + qE[9]*E0J + qE[10]*E0K + qE[11]*E0L + qE[12]*E0M + qE[13]*E0N + qE[14]*E0O + qE[15]*E0P + qE[16]*E0Q + qE[17]*E0R + qE[18]*E0S + qE[19]*E0T + qE[20]*E0U + qE[21]*E0V + qE[22]*E0W + qE[23]*E0X + qE[24]*E0Y + qE[25]*E0Z);
    
    ydot[5] = mu0F - (lambda0F + psi + (qF[0] + qF[1] + qF[2] + qF[3] + qF[4] + qF[6] + qF[7] + qF[8] + qF[9] + qF[10] + qF[11] + qF[12] + qF[13] + qF[14] + qF[15] + qF[16] + qF[17] + qF[18] + qF[19] + qF[20] + qF[21] + qF[22] + qF[23] + qF[24] + qF[25]) + mu0F) * E0F + lambda0F*E0F*E0F + (qF[0]*E0A + qF[1]*E0B + qF[2]*E0C + qF[3]*E0D + qF[4]*E0E + qF[6]*E0G + qF[7]*E0H + qF[8]*E0I + qF[9]*E0J + qF[10]*E0K + qF[11]*E0L + qF[12]*E0M + qF[13]*E0N + qF[14]*E0O + qF[15]*E0P + qF[16]*E0Q + qF[17]*E0R + qF[18]*E0S + qF[19]*E0T + qF[20]*E0U + qF[21]*E0V + qF[22]*E0W + qF[23]*E0X + qF[24]*E0Y + qF[25]*E0Z);
    
    ydot[6] = mu0G - (lambda0G + psi + (qG[0] + qG[1] + qG[2] + qG[3] + qG[4] + qG[5] + qG[7] + qG[8] + qG[9] + qG[10] + qG[11] + qG[12] + qG[13] + qG[14] + qG[15] + qG[16] + qG[17] + qG[18] + qG[19] + qG[20] + qG[21] + qG[22] + qG[23] + qG[24] + qG[25]) + mu0G) * E0G + lambda0G*E0G*E0G + (qG[0]*E0A + qG[1]*E0B + qG[2]*E0C + qG[3]*E0D + qG[4]*E0E + qG[5]*E0F + qG[7]*E0H + qG[8]*E0I + qG[9]*E0J + qG[10]*E0K + qG[11]*E0L + qG[12]*E0M + qG[13]*E0N + qG[14]*E0O + qG[15]*E0P + qG[16]*E0Q + qG[17]*E0R + qG[18]*E0S + qG[19]*E0T + qG[20]*E0U + qG[21]*E0V + qG[22]*E0W + qG[23]*E0X + qG[24]*E0Y + qG[25]*E0Z);
    
    ydot[7] = mu0H - (lambda0H + psi + (qH[0] + qH[1] + qH[2] + qH[3] + qH[4] + qH[5] + qH[6] + qH[8] + qH[9] + qH[10] + qH[11] + qH[12] + qH[13] + qH[14] + qH[15] + qH[16] + qH[17] + qH[18] + qH[19] + qH[20] + qH[21] + qH[22] + qH[23] + qH[24] + qH[25]) + mu0H) * E0H + lambda0H*E0H*E0H + (qH[0]*E0A + qH[1]*E0B + qH[2]*E0C + qH[3]*E0D + qH[4]*E0E + qH[5]*E0F + qH[6]*E0G + qH[8]*E0I + qH[9]*E0J + qH[10]*E0K + qH[11]*E0L + qH[12]*E0M + qH[13]*E0N + qH[14]*E0O + qH[15]*E0P + qH[16]*E0Q + qH[17]*E0R + qH[18]*E0S + qH[19]*E0T + qH[20]*E0U + qH[21]*E0V + qH[22]*E0W + qH[23]*E0X + qH[24]*E0Y + qH[25]*E0Z);
    
    ydot[8] = mu0I - (lambda0I + psi + (qI[0] + qI[1] + qI[2] + qI[3] + qI[4] + qI[5] + qI[6] + qI[7] + qI[9] + qI[10] + qI[11] + qI[12] + qI[13] + qI[14] + qI[15] + qI[16] + qI[17] + qI[18] + qI[19] + qI[20] + qI[21] + qI[22] + qI[23] + qI[24] + qI[25]) + mu0I) * E0I + lambda0I*E0I*E0I + (qI[0]*E0A + qI[1]*E0B + qI[2]*E0C + qI[3]*E0D + qI[4]*E0E + qI[5]*E0F + qI[6]*E0G + qI[7]*E0H + qI[9]*E0J + qI[10]*E0K + qI[11]*E0L + qI[12]*E0M + qI[13]*E0N + qI[14]*E0O + qI[15]*E0P + qI[16]*E0Q + qI[17]*E0R + qI[18]*E0S + qI[19]*E0T + qI[20]*E0U + qI[21]*E0V + qI[22]*E0W + qI[23]*E0X + qI[24]*E0Y + qI[25]*E0Z);
    
    ydot[9] = mu0J - (lambda0J + psi + (qJ[0] + qJ[1] + qJ[2] + qJ[3] + qJ[4] + qJ[5] + qJ[6] + qJ[7] + qJ[8] + qJ[10] + qJ[11] + qJ[12] + qJ[13] + qJ[14] + qJ[15] + qJ[16] + qJ[17] + qJ[18] + qJ[19] + qJ[20] + qJ[21] + qJ[22] + qJ[23] + qJ[24] + qJ[25]) + mu0J) * E0J + lambda0J*E0J*E0J + (qJ[0]*E0A + qJ[1]*E0B + qJ[2]*E0C + qJ[3]*E0D + qJ[4]*E0E + qJ[5]*E0F + qJ[6]*E0G + qJ[7]*E0H + qJ[8]*E0I + qJ[10]*E0K + qJ[11]*E0L + qJ[12]*E0M + qJ[13]*E0N + qJ[14]*E0O + qJ[15]*E0P + qJ[16]*E0Q + qJ[17]*E0R + qJ[18]*E0S + qJ[19]*E0T + qJ[20]*E0U + qJ[21]*E0V + qJ[22]*E0W + qJ[23]*E0X + qJ[24]*E0Y + qJ[25]*E0Z);
    
    ydot[10] = mu0K - (lambda0K + psi + (qK[0] + qK[1] + qK[2] + qK[3] + qK[4] + qK[5] + qK[6] + qK[7] + qK[8] + qK[9] + qK[11] + qK[12] + qK[13] + qK[14] + qK[15] + qK[16] + qK[17] + qK[18] + qK[19] + qK[20] + qK[21] + qK[22] + qK[23] + qK[24] + qK[25]) + mu0K) * E0K + lambda0K*E0K*E0K + (qK[0]*E0A + qK[1]*E0B + qK[2]*E0C + qK[3]*E0D + qK[4]*E0E + qK[5]*E0F + qK[6]*E0G + qK[7]*E0H + qK[8]*E0I + qK[9]*E0J + qK[11]*E0L + qK[12]*E0M + qK[13]*E0N + qK[14]*E0O + qK[15]*E0P + qK[16]*E0Q + qK[17]*E0R + qK[18]*E0S + qK[19]*E0T + qK[20]*E0U + qK[21]*E0V + qK[22]*E0W + qK[23]*E0X + qK[24]*E0Y + qK[25]*E0Z);
    
    ydot[11] = mu0L - (lambda0L + psi + (qL[0] + qL[1] + qL[2] + qL[3] + qL[4] + qL[5] + qL[6] + qL[7] + qL[8] + qL[9] + qL[10] + qL[12] + qL[13] + qL[14] + qL[15] + qL[16] + qL[17] + qL[18] + qL[19] + qL[20] + qL[21] + qL[22] + qL[23] + qL[24] + qL[25]) + mu0L) * E0L + lambda0L*E0L*E0L + (qL[0]*E0A + qL[1]*E0B + qL[2]*E0C + qL[3]*E0D + qL[4]*E0E + qL[5]*E0F + qL[6]*E0G + qL[7]*E0H + qL[8]*E0I + qL[9]*E0J + qL[10]*E0K + qL[12]*E0M + qL[13]*E0N + qL[14]*E0O + qL[15]*E0P + qL[16]*E0Q + qL[17]*E0R + qL[18]*E0S + qL[19]*E0T + qL[20]*E0U + qL[21]*E0V + qL[22]*E0W + qL[23]*E0X + qL[24]*E0Y + qL[25]*E0Z);
    
    ydot[12] = mu0M - (lambda0M + psi + (qM[0] + qM[1] + qM[2] + qM[3] + qM[4] + qM[5] + qM[6] + qM[7] + qM[8] + qM[9] + qM[10] + qM[11] + qM[13] + qM[14] + qM[15] + qM[16] + qM[17] + qM[18] + qM[19] + qM[20] + qM[21] + qM[22] + qM[23] + qM[24] + qM[25]) + mu0M) * E0M + lambda0M*E0M*E0M + (qM[0]*E0A + qM[1]*E0B + qM[2]*E0C + qM[3]*E0D + qM[4]*E0E + qM[5]*E0F + qM[6]*E0G + qM[7]*E0H + qM[8]*E0I + qM[9]*E0J + qM[10]*E0K + qM[11]*E0L + qM[13]*E0N + qM[14]*E0O + qM[15]*E0P + qM[16]*E0Q + qM[17]*E0R + qM[18]*E0S + qM[19]*E0T + qM[20]*E0U + qM[21]*E0V + qM[22]*E0W + qM[23]*E0X + qM[24]*E0Y + qM[25]*E0Z);
    
    ydot[13] = mu0N - (lambda0N + psi + (qN[0] + qN[1] + qN[2] + qN[3] + qN[4] + qN[5] + qN[6] + qN[7] + qN[8] + qN[9] + qN[10] + qN[11] + qN[12] + qN[14] + qN[15] + qN[16] + qN[17] + qN[18] + qN[19] + qN[20] + qN[21] + qN[22] + qN[23] + qN[24] + qN[25]) + mu0N) * E0N + lambda0N*E0N*E0N + (qN[0]*E0A + qN[1]*E0B + qN[2]*E0C + qN[3]*E0D + qN[4]*E0E + qN[5]*E0F + qN[6]*E0G + qN[7]*E0H + qN[8]*E0I + qN[9]*E0J + qN[10]*E0K + qN[11]*E0L + qN[12]*E0M + qN[14]*E0O + qN[15]*E0P + qN[16]*E0Q + qN[17]*E0R + qN[18]*E0S + qN[19]*E0T + qN[20]*E0U + qN[21]*E0V + qN[22]*E0W + qN[23]*E0X + qN[24]*E0Y + qN[25]*E0Z);
    
    ydot[14] = mu0O - (lambda0O + psi + (qO[0] + qO[1] + qO[2] + qO[3] + qO[4] + qO[5] + qO[6] + qO[7] + qO[8] + qO[9] + qO[10] + qO[11] + qO[12] + qO[13] + qO[15] + qO[16] + qO[17] + qO[18] + qO[19] + qO[20] + qO[21] + qO[22] + qO[23] + qO[24] + qO[25]) + mu0O) * E0O + lambda0O*E0O*E0O + (qO[0]*E0A + qO[1]*E0B + qO[2]*E0C + qO[3]*E0D + qO[4]*E0E + qO[5]*E0F + qO[6]*E0G + qO[7]*E0H + qO[8]*E0I + qO[9]*E0J + qO[10]*E0K + qO[11]*E0L + qO[12]*E0M + qO[13]*E0N + qO[15]*E0P + qO[16]*E0Q + qO[17]*E0R + qO[18]*E0S + qO[19]*E0T + qO[20]*E0U + qO[21]*E0V + qO[22]*E0W + qO[23]*E0X + qO[24]*E0Y + qO[25]*E0Z);
    
    ydot[15] = mu0P - (lambda0P + psi + (qP[0] + qP[1] + qP[2] + qP[3] + qP[4] + qP[5] + qP[6] + qP[7] + qP[8] + qP[9] + qP[10] + qP[11] + qP[12] + qP[13] + qP[14] + qP[16] + qP[17] + qP[18] + qP[19] + qP[20] + qP[21] + qP[22] + qP[23] + qP[24] + qP[25]) + mu0P) * E0P + lambda0P*E0P*E0P + (qP[0]*E0A + qP[1]*E0B + qP[2]*E0C + qP[3]*E0D + qP[4]*E0E + qP[5]*E0F + qP[6]*E0G + qP[7]*E0H + qP[8]*E0I + qP[9]*E0J + qP[10]*E0K + qP[11]*E0L + qP[12]*E0M + qP[13]*E0N + qP[14]*E0O + qP[16]*E0Q + qP[17]*E0R + qP[18]*E0S + qP[19]*E0T + qP[20]*E0U + qP[21]*E0V + qP[22]*E0W + qP[23]*E0X + qP[24]*E0Y + qP[25]*E0Z);
    
    ydot[16] = mu0Q - (lambda0Q + psi + (qQ[0] + qQ[1] + qQ[2] + qQ[3] + qQ[4] + qQ[5] + qQ[6] + qQ[7] + qQ[8] + qQ[9] + qQ[10] + qQ[11] + qQ[12] + qQ[13] + qQ[14] + qQ[15] + qQ[17] + qQ[18] + qQ[19] + qQ[20] + qQ[21] + qQ[22] + qQ[23] + qQ[24] + qQ[25]) + mu0Q) * E0Q + lambda0Q*E0Q*E0Q + (qQ[0]*E0A + qQ[1]*E0B + qQ[2]*E0C + qQ[3]*E0D + qQ[4]*E0E + qQ[5]*E0F + qQ[6]*E0G + qQ[7]*E0H + qQ[8]*E0I + qQ[9]*E0J + qQ[10]*E0K + qQ[11]*E0L + qQ[12]*E0M + qQ[13]*E0N + qQ[14]*E0O + qQ[15]*E0P + qQ[17]*E0R + qQ[18]*E0S + qQ[19]*E0T + qQ[20]*E0U + qQ[21]*E0V + qQ[22]*E0W + qQ[23]*E0X + qQ[24]*E0Y + qQ[25]*E0Z);
    
    ydot[17] = mu0R - (lambda0R + psi + (qR[0] + qR[1] + qR[2] + qR[3] + qR[4] + qR[5] + qR[6] + qR[7] + qR[8] + qR[9] + qR[10] + qR[11] + qR[12] + qR[13] + qR[14] + qR[15] + qR[16] + qR[18] + qR[19] + qR[20] + qR[21] + qR[22] + qR[23] + qR[24] + qR[25]) + mu0R) * E0R + lambda0R*E0R*E0R + (qR[0]*E0A + qR[1]*E0B + qR[2]*E0C + qR[3]*E0D + qR[4]*E0E + qR[5]*E0F + qR[6]*E0G + qR[7]*E0H + qR[8]*E0I + qR[9]*E0J + qR[10]*E0K + qR[11]*E0L + qR[12]*E0M + qR[13]*E0N + qR[14]*E0O + qR[15]*E0P + qR[16]*E0Q + qR[18]*E0S + qR[19]*E0T + qR[20]*E0U + qR[21]*E0V + qR[22]*E0W + qR[23]*E0X + qR[24]*E0Y + qR[25]*E0Z);
    
    ydot[18] = mu0S - (lambda0S + psi + (qS[0] + qS[1] + qS[2] + qS[3] + qS[4] + qS[5] + qS[6] + qS[7] + qS[8] + qS[9] + qS[10] + qS[11] + qS[12] + qS[13] + qS[14] + qS[15] + qS[16] + qS[17] + qS[19] + qS[20] + qS[21] + qS[22] + qS[23] + qS[24] + qS[25]) + mu0S) * E0S + lambda0S*E0S*E0S + (qS[0]*E0A + qS[1]*E0B + qS[2]*E0C + qS[3]*E0D + qS[4]*E0E + qS[5]*E0F + qS[6]*E0G + qS[7]*E0H + qS[8]*E0I + qS[9]*E0J + qS[10]*E0K + qS[11]*E0L + qS[12]*E0M + qS[13]*E0N + qS[14]*E0O + qS[15]*E0P + qS[16]*E0Q + qS[17]*E0R + qS[19]*E0T + qS[20]*E0U + qS[21]*E0V + qS[22]*E0W + qS[23]*E0X + qS[24]*E0Y + qS[25]*E0Z);
    
    ydot[19] = mu0T - (lambda0T + psi + (qT[0] + qT[1] + qT[2] + qT[3] + qT[4] + qT[5] + qT[6] + qT[7] + qT[8] + qT[9] + qT[10] + qT[11] + qT[12] + qT[13] + qT[14] + qT[15] + qT[16] + qT[17] + qT[18] + qT[20] + qT[21] + qT[22] + qT[23] + qT[24] + qT[25]) + mu0T) * E0T + lambda0T*E0T*E0T + (qT[0]*E0A + qT[1]*E0B + qT[2]*E0C + qT[3]*E0D + qT[4]*E0E + qT[5]*E0F + qT[6]*E0G + qT[7]*E0H + qT[8]*E0I + qT[9]*E0J + qT[10]*E0K + qT[11]*E0L + qT[12]*E0M + qT[13]*E0N + qT[14]*E0O + qT[15]*E0P + qT[16]*E0Q + qT[17]*E0R + qT[18]*E0S + qT[20]*E0U + qT[21]*E0V + qT[22]*E0W + qT[23]*E0X + qT[24]*E0Y + qT[25]*E0Z);
    
    ydot[20] = mu0U - (lambda0U + psi + (qU[0] + qU[1] + qU[2] + qU[3] + qU[4] + qU[5] + qU[6] + qU[7] + qU[8] + qU[9] + qU[10] + qU[11] + qU[12] + qU[13] + qU[14] + qU[15] + qU[16] + qU[17] + qU[18] + qU[19] + qU[21] + qU[22] + qU[23] + qU[24] + qU[25]) + mu0U) * E0U + lambda0U*E0U*E0U + (qU[0]*E0A + qU[1]*E0B + qU[2]*E0C + qU[3]*E0D + qU[4]*E0E + qU[5]*E0F + qU[6]*E0G + qU[7]*E0H + qU[8]*E0I + qU[9]*E0J + qU[10]*E0K + qU[11]*E0L + qU[12]*E0M + qU[13]*E0N + qU[14]*E0O + qU[15]*E0P + qU[16]*E0Q + qU[17]*E0R + qU[18]*E0S + qU[19]*E0T + qU[21]*E0V + qU[22]*E0W + qU[23]*E0X + qU[24]*E0Y + qU[25]*E0Z);
    
    ydot[21] = mu0V - (lambda0V + psi + (qV[0] + qV[1] + qV[2] + qV[3] + qV[4] + qV[5] + qV[6] + qV[7] + qV[8] + qV[9] + qV[10] + qV[11] + qV[12] + qV[13] + qV[14] + qV[15] + qV[16] + qV[17] + qV[18] + qV[19] + qV[20] + qV[22] + qV[23] + qV[24] + qV[25]) + mu0V) * E0V + lambda0V*E0V*E0V + (qV[0]*E0A + qV[1]*E0B + qV[2]*E0C + qV[3]*E0D + qV[4]*E0E + qV[5]*E0F + qV[6]*E0G + qV[7]*E0H + qV[8]*E0I + qV[9]*E0J + qV[10]*E0K + qV[11]*E0L + qV[12]*E0M + qV[13]*E0N + qV[14]*E0O + qV[15]*E0P + qV[16]*E0Q + qV[17]*E0R + qV[18]*E0S + qV[19]*E0T + qV[20]*E0U + qV[22]*E0W + qV[23]*E0X + qV[24]*E0Y + qV[25]*E0Z);
    
    ydot[22] = mu0W - (lambda0W + psi + (qW[0] + qW[1] + qW[2] + qW[3] + qW[4] + qW[5] + qW[6] + qW[7] + qW[8] + qW[9] + qW[10] + qW[11] + qW[12] + qW[13] + qW[14] + qW[15] + qW[16] + qW[17] + qW[18] + qW[19] + qW[20] + qW[21] + qW[23] + qW[24] + qW[25]) + mu0W) * E0W + lambda0W*E0W*E0W + (qW[0]*E0A + qW[1]*E0B + qW[2]*E0C + qW[3]*E0D + qW[4]*E0E + qW[5]*E0F + qW[6]*E0G + qW[7]*E0H + qW[8]*E0I + qW[9]*E0J + qW[10]*E0K + qW[11]*E0L + qW[12]*E0M + qW[13]*E0N + qW[14]*E0O + qW[15]*E0P + qW[16]*E0Q + qW[17]*E0R + qW[18]*E0S + qW[19]*E0T + qW[20]*E0U + qW[21]*E0V + qW[23]*E0X + qW[24]*E0Y + qW[25]*E0Z);
    
    ydot[23] = mu0X - (lambda0X + psi + (qX[0] + qX[1] + qX[2] + qX[3] + qX[4] + qX[5] + qX[6] + qX[7] + qX[8] + qX[9] + qX[10] + qX[11] + qX[12] + qX[13] + qX[14] + qX[15] + qX[16] + qX[17] + qX[18] + qX[19] + qX[20] + qX[21] + qX[22] + qX[24] + qX[25]) + mu0X) * E0X + lambda0X*E0X*E0X + (qX[0]*E0A + qX[1]*E0B + qX[2]*E0C + qX[3]*E0D + qX[4]*E0E + qX[5]*E0F + qX[6]*E0G + qX[7]*E0H + qX[8]*E0I + qX[9]*E0J + qX[10]*E0K + qX[11]*E0L + qX[12]*E0M + qX[13]*E0N + qX[14]*E0O + qX[15]*E0P + qX[16]*E0Q + qX[17]*E0R + qX[18]*E0S + qX[19]*E0T + qX[20]*E0U + qX[21]*E0V + qX[22]*E0W + qX[24]*E0Y + qX[25]*E0Z);
    
    ydot[24] = mu0Y - (lambda0Y + psi + (qY[0] + qY[1] + qY[2] + qY[3] + qY[4] + qY[5] + qY[6] + qY[7] + qY[8] + qY[9] + qY[10] + qY[11] + qY[12] + qY[13] + qY[14] + qY[15] + qY[16] + qY[17] + qY[18] + qY[19] + qY[20] + qY[21] + qY[22] + qY[23] + qY[25]) + mu0Y) * E0Y + lambda0Y*E0Y*E0Y + (qY[0]*E0A + qY[1]*E0B + qY[2]*E0C + qY[3]*E0D + qY[4]*E0E + qY[5]*E0F + qY[6]*E0G + qY[7]*E0H + qY[8]*E0I + qY[9]*E0J + qY[10]*E0K + qY[11]*E0L + qY[12]*E0M + qY[13]*E0N + qY[14]*E0O + qY[15]*E0P + qY[16]*E0Q + qY[17]*E0R + qY[18]*E0S + qY[19]*E0T + qY[20]*E0U + qY[21]*E0V + qY[22]*E0W + qY[23]*E0X + qY[25]*E0Z);
    
    ydot[25] = mu0Z - (lambda0Z + psi + (qZ[0] + qZ[1] + qZ[2] + qZ[3] + qZ[4] + qZ[5] + qZ[6] + qZ[7] + qZ[8] + qZ[9] + qZ[10] + qZ[11] + qZ[12] + qZ[13] + qZ[14] + qZ[15] + qZ[16] + qZ[17] + qZ[18] + qZ[19] + qZ[20] + qZ[21] + qZ[22] + qZ[23] + qZ[24]) + mu0Z) * E0Z + lambda0Z*E0Z*E0Z + (qZ[0]*E0A + qZ[1]*E0B + qZ[2]*E0C + qZ[3]*E0D + qZ[4]*E0E + qZ[5]*E0F + qZ[6]*E0G + qZ[7]*E0H + qZ[8]*E0I + qZ[9]*E0J + qZ[10]*E0K + qZ[11]*E0L + qZ[12]*E0M + qZ[13]*E0N + qZ[14]*E0O + qZ[15]*E0P + qZ[16]*E0Q + qZ[17]*E0R + qZ[18]*E0S + qZ[19]*E0T + qZ[20]*E0U + qZ[21]*E0V + qZ[22]*E0W + qZ[23]*E0X + qZ[24]*E0Y);
    
    
    /* The D's */
    ydot[26] =  - (lambda0A + psi + (qA[1] + qA[2] + qA[3] + qA[4] + qA[5] + qA[6] + qA[7] + qA[8] + qA[9] + qA[10] + qA[11] + qA[12] + qA[13] + qA[14] + qA[15] + qA[16] + qA[17] + qA[18] + qA[19] + qA[20] + qA[21] + qA[22] + qA[23] + qA[24] + qA[25]) + mu0A) * D0A + 2*lambda0A*E0A*D0A + (qA[1]*D0B + qA[2]*D0C + qA[3]*D0D + qA[4]*D0E + qA[5]*D0F + qA[6]*D0G + qA[7]*D0H + qA[8]*D0I + qA[9]*D0J + qA[10]*D0K + qA[11]*D0L + qA[12]*D0M + qA[13]*D0N + qA[14]*D0O + qA[15]*D0P + qA[16]*D0Q + qA[17]*D0R + qA[18]*D0S + qA[19]*D0T + qA[20]*D0U + qA[21]*D0V + qA[22]*D0W + qA[23]*D0X + qA[24]*D0Y + qA[25]*D0Z);
    
    ydot[27] =  - (lambda0B + psi + (qB[0] + qB[2] + qB[3] + qB[4] + qB[5] + qB[6] + qB[7] + qB[8] + qB[9] + qB[10] + qB[11] + qB[12] + qB[13] + qB[14] + qB[15] + qB[16] + qB[17] + qB[18] + qB[19] + qB[20] + qB[21] + qB[22] + qB[23] + qB[24] + qB[25]) + mu0B) * D0B + 2*lambda0B*E0B*D0B + (qB[0]*D0A + qB[2]*D0C + qB[3]*D0D + qB[4]*D0E + qB[5]*D0F + qB[6]*D0G + qB[7]*D0H + qB[8]*D0I + qB[9]*D0J + qB[10]*D0K + qB[11]*D0L + qB[12]*D0M + qB[13]*D0N + qB[14]*D0O + qB[15]*D0P + qB[16]*D0Q + qB[17]*D0R + qB[18]*D0S + qB[19]*D0T + qB[20]*D0U + qB[21]*D0V + qB[22]*D0W + qB[23]*D0X + qB[24]*D0Y + qB[25]*D0Z);
    
    ydot[28] =  - (lambda0C + psi + (qC[0] + qC[1] + qC[3] + qC[4] + qC[5] + qC[6] + qC[7] + qC[8] + qC[9] + qC[10] + qC[11] + qC[12] + qC[13] + qC[14] + qC[15] + qC[16] + qC[17] + qC[18] + qC[19] + qC[20] + qC[21] + qC[22] + qC[23] + qC[24] + qC[25]) + mu0C) * D0C + 2*lambda0C*E0C*D0C + (qC[0]*D0A + qC[1]*D0B + qC[3]*D0D + qC[4]*D0E + qC[5]*D0F + qC[6]*D0G + qC[7]*D0H + qC[8]*D0I + qC[9]*D0J + qC[10]*D0K + qC[11]*D0L + qC[12]*D0M + qC[13]*D0N + qC[14]*D0O + qC[15]*D0P + qC[16]*D0Q + qC[17]*D0R + qC[18]*D0S + qC[19]*D0T + qC[20]*D0U + qC[21]*D0V + qC[22]*D0W + qC[23]*D0X + qC[24]*D0Y + qC[25]*D0Z);
    
    ydot[29] =  - (lambda0D + psi + (qD[0] + qD[1] + qD[2] + qD[4] + qD[5] + qD[6] + qD[7] + qD[8] + qD[9] + qD[10] + qD[11] + qD[12] + qD[13] + qD[14] + qD[15] + qD[16] + qD[17] + qD[18] + qD[19] + qD[20] + qD[21] + qD[22] + qD[23] + qD[24] + qD[25])+ mu0D) * D0D + 2*lambda0D*E0D*D0D + (qD[0]*D0A + qD[1]*D0B + qD[2]*D0C + qD[4]*D0E + qD[5]*D0F + qD[6]*D0G + qD[7]*D0H + qD[8]*D0I + qD[9]*D0J + qD[10]*D0K + qD[11]*D0L + qD[12]*D0M + qD[13]*D0N + qD[14]*D0O + qD[15]*D0P + qD[16]*D0Q + qD[17]*D0R + qD[18]*D0S + qD[19]*D0T + qD[20]*D0U + qD[21]*D0V + qD[22]*D0W + qD[23]*D0X + qD[24]*D0Y + qD[25]*D0Z);
    
    ydot[30] =  - (lambda0E + psi + (qE[0] + qE[1] + qE[2] + qE[3] + qE[5] + qE[6] + qE[7] + qE[8] + qE[9] + qE[10] + qE[11] + qE[12] + qE[13] + qE[14] + qE[15] + qE[16] + qE[17] + qE[18] + qE[19] + qE[20] + qE[21] + qE[22] + qE[23] + qE[24] + qE[25])+ mu0E) * D0E + 2*lambda0E*E0E*D0E + (qE[0]*D0A + qE[1]*D0B + qE[2]*D0C + qE[3]*D0D + qE[5]*D0F + qE[6]*D0G + qE[7]*D0H + qE[8]*D0I + qE[9]*D0J + qE[10]*D0K + qE[11]*D0L + qE[12]*D0M + qE[13]*D0N + qE[14]*D0O + qE[15]*D0P + qE[16]*D0Q + qE[17]*D0R + qE[18]*D0S + qE[19]*D0T + qE[20]*D0U + qE[21]*D0V + qE[22]*D0W + qE[23]*D0X + qE[24]*D0Y + qE[25]*D0Z);
    
    ydot[31] =  - (lambda0F + psi + (qF[0] + qF[1] + qF[2] + qF[3] + qF[4] + qF[6] + qF[7] + qF[8] + qF[9] + qF[10] + qF[11] + qF[12] + qF[13] + qF[14] + qF[15] + qF[16] + qF[17] + qF[18] + qF[19] + qF[20] + qF[21] + qF[22] + qF[23] + qF[24] + qF[25]) + mu0F) * D0F + 2*lambda0F*E0F*D0F + (qF[0]*D0A + qF[1]*D0B + qF[2]*D0C + qF[3]*D0D + qF[4]*D0E + qF[6]*D0G + qF[7]*D0H + qF[8]*D0I + qF[9]*D0J + qF[10]*D0K + qF[11]*D0L + qF[12]*D0M + qF[13]*D0N + qF[14]*D0O + qF[15]*D0P + qF[16]*D0Q + qF[17]*D0R + qF[18]*D0S + qF[19]*D0T + qF[20]*D0U + qF[21]*D0V + qF[22]*D0W + qF[23]*D0X + qF[24]*D0Y + qF[25]*D0Z);
    
    ydot[32] =  - (lambda0G + psi + (qG[0] + qG[1] + qG[2] + qG[3] + qG[4] + qG[5] + qG[7] + qG[8] + qG[9] + qG[10] + qG[11] + qG[12] + qG[13] + qG[14] + qG[15] + qG[16] + qG[17] + qG[18] + qG[19] + qG[20] + qG[21] + qG[22] + qG[23] + qG[24] + qG[25]) + mu0G) * D0G + 2*lambda0G*E0G*D0G + (qG[0]*D0A + qG[1]*D0B + qG[2]*D0C + qG[3]*D0D + qG[4]*D0E + qG[5]*D0F + qG[7]*D0H + qG[8]*D0I + qG[9]*D0J + qG[10]*D0K + qG[11]*D0L + qG[12]*D0M + qG[13]*D0N + qG[14]*D0O + qG[15]*D0P + qG[16]*D0Q + qG[17]*D0R + qG[18]*D0S + qG[19]*D0T + qG[20]*D0U + qG[21]*D0V + qG[22]*D0W + qG[23]*D0X + qG[24]*D0Y + qG[25]*D0Z);
    
    ydot[33] =  - (lambda0H + psi + (qH[0] + qH[1] + qH[2] + qH[3] + qH[4] + qH[5] + qH[6] + qH[8] + qH[9] + qH[10] + qH[11] + qH[12] + qH[13] + qH[14] + qH[15] + qH[16] + qH[17] + qH[18] + qH[19] + qH[20] + qH[21] + qH[22] + qH[23] + qH[24] + qH[25]) + mu0H) * D0H + 2*lambda0H*E0H*D0H + (qH[0]*D0A + qH[1]*D0B + qH[2]*D0C + qH[3]*D0D + qH[4]*D0E + qH[5]*D0F + qH[6]*D0G + qH[8]*D0I + qH[9]*D0J + qH[10]*D0K + qH[11]*D0L + qH[12]*D0M + qH[13]*D0N + qH[14]*D0O + qH[15]*D0P + qH[16]*D0Q + qH[17]*D0R + qH[18]*D0S + qH[19]*D0T + qH[20]*D0U + qH[21]*D0V + qH[22]*D0W + qH[23]*D0X + qH[24]*D0Y + qH[25]*D0Z);
    
    ydot[34] =  - (lambda0I + psi + (qI[0] + qI[1] + qI[2] + qI[3] + qI[4] + qI[5] + qI[6] + qI[7] + qI[9] + qI[10] + qI[11] + qI[12] + qI[13] + qI[14] + qI[15] + qI[16] + qI[17] + qI[18] + qI[19] + qI[20] + qI[21] + qI[22] + qI[23] + qI[24] + qI[25]) + mu0I) * D0I + 2*lambda0I*E0I*D0I + (qI[0]*D0A + qI[1]*D0B + qI[2]*D0C + qI[3]*D0D + qI[4]*D0E + qI[5]*D0F + qI[6]*D0G + qI[7]*D0H + qI[9]*D0J + qI[10]*D0K + qI[11]*D0L + qI[12]*D0M + qI[13]*D0N + qI[14]*D0O + qI[15]*D0P + qI[16]*D0Q + qI[17]*D0R + qI[18]*D0S + qI[19]*D0T + qI[20]*D0U + qI[21]*D0V + qI[22]*D0W + qI[23]*D0X + qI[24]*D0Y + qI[25]*D0Z);
    
    ydot[35] =  - (lambda0J + psi + (qJ[0] + qJ[1] + qJ[2] + qJ[3] + qJ[4] + qJ[5] + qJ[6] + qJ[7] + qJ[8] + qJ[10] + qJ[11] + qJ[12] + qJ[13] + qJ[14] + qJ[15] + qJ[16] + qJ[17] + qJ[18] + qJ[19] + qJ[20] + qJ[21] + qJ[22] + qJ[23] + qJ[24] + qJ[25]) + mu0J) * D0J + 2*lambda0J*E0J*D0J + (qJ[0]*D0A + qJ[1]*D0B + qJ[2]*D0C + qJ[3]*D0D + qJ[4]*D0E + qJ[5]*D0F + qJ[6]*D0G + qJ[7]*D0H + qJ[8]*D0I + qJ[10]*D0K + qJ[11]*D0L + qJ[12]*D0M + qJ[13]*D0N + qJ[14]*D0O + qJ[15]*D0P + qJ[16]*D0Q + qJ[17]*D0R + qJ[18]*D0S + qJ[19]*D0T + qJ[20]*D0U + qJ[21]*D0V + qJ[22]*D0W + qJ[23]*D0X + qJ[24]*D0Y + qJ[25]*D0Z);
    
    ydot[36] =  - (lambda0K + psi + (qK[0] + qK[1] + qK[2] + qK[3] + qK[4] + qK[5] + qK[6] + qK[7] + qK[8] + qK[9] + qK[11] + qK[12] + qK[13] + qK[14] + qK[15] + qK[16] + qK[17] + qK[18] + qK[19] + qK[20] + qK[21] + qK[22] + qK[23] + qK[24] + qK[25]) + mu0K) * D0K + 2*lambda0K*E0K*D0K + (qK[0]*D0A + qK[1]*D0B + qK[2]*D0C + qK[3]*D0D + qK[4]*D0E + qK[5]*D0F + qK[6]*D0G + qK[7]*D0H + qK[8]*D0I + qK[9]*D0J + qK[11]*D0L + qK[12]*D0M + qK[13]*D0N + qK[14]*D0O + qK[15]*D0P + qK[16]*D0Q + qK[17]*D0R + qK[18]*D0S + qK[19]*D0T + qK[20]*D0U + qK[21]*D0V + qK[22]*D0W + qK[23]*D0X + qK[24]*D0Y + qK[25]*D0Z);
    
    ydot[37] =  - (lambda0L + psi + (qL[0] + qL[1] + qL[2] + qL[3] + qL[4] + qL[5] + qL[6] + qL[7] + qL[8] + qL[9] + qL[10] + qL[12] + qL[13] + qL[14] + qL[15] + qL[16] + qL[17] + qL[18] + qL[19] + qL[20] + qL[21] + qL[22] + qL[23] + qL[24] + qL[25]) + mu0L) * D0L + 2*lambda0L*E0L*D0L + (qL[0]*D0A + qL[1]*D0B + qL[2]*D0C + qL[3]*D0D + qL[4]*D0E + qL[5]*D0F + qL[6]*D0G + qL[7]*D0H + qL[8]*D0I + qL[9]*D0J + qL[10]*D0K + qL[12]*D0M + qL[13]*D0N + qL[14]*D0O + qL[15]*D0P + qL[16]*D0Q + qL[17]*D0R + qL[18]*D0S + qL[19]*D0T + qL[20]*D0U + qL[21]*D0V + qL[22]*D0W + qL[23]*D0X + qL[24]*D0Y + qL[25]*D0Z);
    
    ydot[38] =  - (lambda0M + psi + (qM[0] + qM[1] + qM[2] + qM[3] + qM[4] + qM[5] + qM[6] + qM[7] + qM[8] + qM[9] + qM[10] + qM[11] + qM[13] + qM[14] + qM[15] + qM[16] + qM[17] + qM[18] + qM[19] + qM[20] + qM[21] + qM[22] + qM[23] + qM[24] + qM[25]) + mu0M) * D0M + 2*lambda0M*E0M*D0M + (qM[0]*D0A + qM[1]*D0B + qM[2]*D0C + qM[3]*D0D + qM[4]*D0E + qM[5]*D0F + qM[6]*D0G + qM[7]*D0H + qM[8]*D0I + qM[9]*D0J + qM[10]*D0K + qM[11]*D0L + qM[13]*D0N + qM[14]*D0O + qM[15]*D0P + qM[16]*D0Q + qM[17]*D0R + qM[18]*D0S + qM[19]*D0T + qM[20]*D0U + qM[21]*D0V + qM[22]*D0W + qM[23]*D0X + qM[24]*D0Y + qM[25]*D0Z);
    
    ydot[39] =  - (lambda0N + psi + (qN[0] + qN[1] + qN[2] + qN[3] + qN[4] + qN[5] + qN[6] + qN[7] + qN[8] + qN[9] + qN[10] + qN[11] + qN[12] + qN[14] + qN[15] + qN[16] + qN[17] + qN[18] + qN[19] + qN[20] + qN[21] + qN[22] + qN[23] + qN[24] + qN[25]) + mu0N) * D0N + 2*lambda0N*E0N*D0N + (qN[0]*D0A + qN[1]*D0B + qN[2]*D0C + qN[3]*D0D + qN[4]*D0E + qN[5]*D0F + qN[6]*D0G + qN[7]*D0H + qN[8]*D0I + qN[9]*D0J + qN[10]*D0K + qN[11]*D0L + qN[12]*D0M + qN[14]*D0O + qN[15]*D0P + qN[16]*D0Q + qN[17]*D0R + qN[18]*D0S + qN[19]*D0T + qN[20]*D0U + qN[21]*D0V + qN[22]*D0W + qN[23]*D0X + qN[24]*D0Y + qN[25]*D0Z);
    
    ydot[40] =  - (lambda0O + psi + (qO[0] + qO[1] + qO[2] + qO[3] + qO[4] + qO[5] + qO[6] + qO[7] + qO[8] + qO[9] + qO[10] + qO[11] + qO[12] + qO[13] + qO[15] + qO[16] + qO[17] + qO[18] + qO[19] + qO[20] + qO[21] + qO[22] + qO[23] + qO[24] + qO[25]) + mu0O) * D0O + 2*lambda0O*E0O*D0O + (qO[0]*D0A + qO[1]*D0B + qO[2]*D0C + qO[3]*D0D + qO[4]*D0E + qO[5]*D0F + qO[6]*D0G + qO[7]*D0H + qO[8]*D0I + qO[9]*D0J + qO[10]*D0K + qO[11]*D0L + qO[12]*D0M + qO[13]*D0N + qO[15]*D0P + qO[16]*D0Q + qO[17]*D0R + qO[18]*D0S + qO[19]*D0T + qO[20]*D0U + qO[21]*D0V + qO[22]*D0W + qO[23]*D0X + qO[24]*D0Y + qO[25]*D0Z);
    
    ydot[41] =  - (lambda0P + psi + (qP[0] + qP[1] + qP[2] + qP[3] + qP[4] + qP[5] + qP[6] + qP[7] + qP[8] + qP[9] + qP[10] + qP[11] + qP[12] + qP[13] + qP[14] + qP[16] + qP[17] + qP[18] + qP[19] + qP[20] + qP[21] + qP[22] + qP[23] + qP[24] + qP[25]) + mu0P) * D0P + 2*lambda0P*E0P*D0P + (qP[0]*D0A + qP[1]*D0B + qP[2]*D0C + qP[3]*D0D + qP[4]*D0E + qP[5]*D0F + qP[6]*D0G + qP[7]*D0H + qP[8]*D0I + qP[9]*D0J + qP[10]*D0K + qP[11]*D0L + qP[12]*D0M + qP[13]*D0N + qP[14]*D0O + qP[16]*D0Q + qP[17]*D0R + qP[18]*D0S + qP[19]*D0T + qP[20]*D0U + qP[21]*D0V + qP[22]*D0W + qP[23]*D0X + qP[24]*D0Y + qP[25]*D0Z);
    
    ydot[42] =  - (lambda0Q + psi + (qQ[0] + qQ[1] + qQ[2] + qQ[3] + qQ[4] + qQ[5] + qQ[6] + qQ[7] + qQ[8] + qQ[9] + qQ[10] + qQ[11] + qQ[12] + qQ[13] + qQ[14] + qQ[15] + qQ[17] + qQ[18] + qQ[19] + qQ[20] + qQ[21] + qQ[22] + qQ[23] + qQ[24] + qQ[25]) + mu0Q) * D0Q + 2*lambda0Q*E0Q*D0Q + (qQ[0]*D0A + qQ[1]*D0B + qQ[2]*D0C + qQ[3]*D0D + qQ[4]*D0E + qQ[5]*D0F + qQ[6]*D0G + qQ[7]*D0H + qQ[8]*D0I + qQ[9]*D0J + qQ[10]*D0K + qQ[11]*D0L + qQ[12]*D0M + qQ[13]*D0N + qQ[14]*D0O + qQ[15]*D0P + qQ[17]*D0R + qQ[18]*D0S + qQ[19]*D0T + qQ[20]*D0U + qQ[21]*D0V + qQ[22]*D0W + qQ[23]*D0X + qQ[24]*D0Y + qQ[25]*D0Z);
    
    ydot[43] =  - (lambda0R + psi + (qR[0] + qR[1] + qR[2] + qR[3] + qR[4] + qR[5] + qR[6] + qR[7] + qR[8] + qR[9] + qR[10] + qR[11] + qR[12] + qR[13] + qR[14] + qR[15] + qR[16] + qR[18] + qR[19] + qR[20] + qR[21] + qR[22] + qR[23] + qR[24] + qR[25]) + mu0R) * D0R + 2*lambda0R*E0R*D0R + (qR[0]*D0A + qR[1]*D0B + qR[2]*D0C + qR[3]*D0D + qR[4]*D0E + qR[5]*D0F + qR[6]*D0G + qR[7]*D0H + qR[8]*D0I + qR[9]*D0J + qR[10]*D0K + qR[11]*D0L + qR[12]*D0M + qR[13]*D0N + qR[14]*D0O + qR[15]*D0P + qR[16]*D0Q + qR[18]*D0S + qR[19]*D0T + qR[20]*D0U + qR[21]*D0V + qR[22]*D0W + qR[23]*D0X + qR[24]*D0Y + qR[25]*D0Z);
    
    ydot[44] =  - (lambda0S + psi + (qS[0] + qS[1] + qS[2] + qS[3] + qS[4] + qS[5] + qS[6] + qS[7] + qS[8] + qS[9] + qS[10] + qS[11] + qS[12] + qS[13] + qS[14] + qS[15] + qS[16] + qS[17] + qS[19] + qS[20] + qS[21] + qS[22] + qS[23] + qS[24] + qS[25]) + mu0S) * D0S + 2*lambda0S*E0S*D0S + (qS[0]*D0A + qS[1]*D0B + qS[2]*D0C + qS[3]*D0D + qS[4]*D0E + qS[5]*D0F + qS[6]*D0G + qS[7]*D0H + qS[8]*D0I + qS[9]*D0J + qS[10]*D0K + qS[11]*D0L + qS[12]*D0M + qS[13]*D0N + qS[14]*D0O + qS[15]*D0P + qS[16]*D0Q + qS[17]*D0R + qS[19]*D0T + qS[20]*D0U + qS[21]*D0V + qS[22]*D0W + qS[23]*D0X + qS[24]*D0Y + qS[25]*D0Z);
    
    ydot[45] =  - (lambda0T + psi + (qT[0] + qT[1] + qT[2] + qT[3] + qT[4] + qT[5] + qT[6] + qT[7] + qT[8] + qT[9] + qT[10] + qT[11] + qT[12] + qT[13] + qT[14] + qT[15] + qT[16] + qT[17] + qT[18] + qT[20] + qT[21] + qT[22] + qT[23] + qT[24] + qT[25]) + mu0T) * D0T + 2*lambda0T*E0T*D0T + (qT[0]*D0A + qT[1]*D0B + qT[2]*D0C + qT[3]*D0D + qT[4]*D0E + qT[5]*D0F + qT[6]*D0G + qT[7]*D0H + qT[8]*D0I + qT[9]*D0J + qT[10]*D0K + qT[11]*D0L + qT[12]*D0M + qT[13]*D0N + qT[14]*D0O + qT[15]*D0P + qT[16]*D0Q + qT[17]*D0R + qT[18]*D0S + qT[20]*D0U + qT[21]*D0V + qT[22]*D0W + qT[23]*D0X + qT[24]*D0Y + qT[25]*D0Z);
    
    ydot[46] =  - (lambda0U + psi + (qU[0] + qU[1] + qU[2] + qU[3] + qU[4] + qU[5] + qU[6] + qU[7] + qU[8] + qU[9] + qU[10] + qU[11] + qU[12] + qU[13] + qU[14] + qU[15] + qU[16] + qU[17] + qU[18] + qU[19] + qU[21] + qU[22] + qU[23] + qU[24] + qU[25]) + mu0U) * D0U + 2*lambda0U*E0U*D0U + (qU[0]*D0A + qU[1]*D0B + qU[2]*D0C + qU[3]*D0D + qU[4]*D0E + qU[5]*D0F + qU[6]*D0G + qU[7]*D0H + qU[8]*D0I + qU[9]*D0J + qU[10]*D0K + qU[11]*D0L + qU[12]*D0M + qU[13]*D0N + qU[14]*D0O + qU[15]*D0P + qU[16]*D0Q + qU[17]*D0R + qU[18]*D0S + qU[19]*D0T + qU[21]*D0V + qU[22]*D0W + qU[23]*D0X + qU[24]*D0Y + qU[25]*D0Z);
    
    ydot[47] =  - (lambda0V + psi + (qV[0] + qV[1] + qV[2] + qV[3] + qV[4] + qV[5] + qV[6] + qV[7] + qV[8] + qV[9] + qV[10] + qV[11] + qV[12] + qV[13] + qV[14] + qV[15] + qV[16] + qV[17] + qV[18] + qV[19] + qV[20] + qV[22] + qV[23] + qV[24] + qV[25]) + mu0V) * D0V + 2*lambda0V*E0V*D0V + (qV[0]*D0A + qV[1]*D0B + qV[2]*D0C + qV[3]*D0D + qV[4]*D0E + qV[5]*D0F + qV[6]*D0G + qV[7]*D0H + qV[8]*D0I + qV[9]*D0J + qV[10]*D0K + qV[11]*D0L + qV[12]*D0M + qV[13]*D0N + qV[14]*D0O + qV[15]*D0P + qV[16]*D0Q + qV[17]*D0R + qV[18]*D0S + qV[19]*D0T + qV[20]*D0U + qV[22]*D0W + qV[23]*D0X + qV[24]*D0Y + qV[25]*D0Z);
    
    ydot[48] =  - (lambda0W + psi + (qW[0] + qW[1] + qW[2] + qW[3] + qW[4] + qW[5] + qW[6] + qW[7] + qW[8] + qW[9] + qW[10] + qW[11] + qW[12] + qW[13] + qW[14] + qW[15] + qW[16] + qW[17] + qW[18] + qW[19] + qW[20] + qW[21] + qW[23] + qW[24] + qW[25]) + mu0W) * D0W + 2*lambda0W*E0W*D0W + (qW[0]*D0A + qW[1]*D0B + qW[2]*D0C + qW[3]*D0D + qW[4]*D0E + qW[5]*D0F + qW[6]*D0G + qW[7]*D0H + qW[8]*D0I + qW[9]*D0J + qW[10]*D0K + qW[11]*D0L + qW[12]*D0M + qW[13]*D0N + qW[14]*D0O + qW[15]*D0P + qW[16]*D0Q + qW[17]*D0R + qW[18]*D0S + qW[19]*D0T + qW[20]*D0U + qW[21]*D0V + qW[23]*D0X + qW[24]*D0Y + qW[25]*D0Z);
    
    ydot[49] =  - (lambda0X + psi + (qX[0] + qX[1] + qX[2] + qX[3] + qX[4] + qX[5] + qX[6] + qX[7] + qX[8] + qX[9] + qX[10] + qX[11] + qX[12] + qX[13] + qX[14] + qX[15] + qX[16] + qX[17] + qX[18] + qX[19] + qX[20] + qX[21] + qX[22] + qX[24] + qX[25]) + mu0X) * D0X + 2*lambda0X*E0X*D0X + (qX[0]*D0A + qX[1]*D0B + qX[2]*D0C + qX[3]*D0D + qX[4]*D0E + qX[5]*D0F + qX[6]*D0G + qX[7]*D0H + qX[8]*D0I + qX[9]*D0J + qX[10]*D0K + qX[11]*D0L + qX[12]*D0M + qX[13]*D0N + qX[14]*D0O + qX[15]*D0P + qX[16]*D0Q + qX[17]*D0R + qX[18]*D0S + qX[19]*D0T + qX[20]*D0U + qX[21]*D0V + qX[22]*D0W + qX[24]*D0Y + qX[25]*D0Z);
    
    ydot[50] =  - (lambda0Y + psi + (qY[0] + qY[1] + qY[2] + qY[3] + qY[4] + qY[5] + qY[6] + qY[7] + qY[8] + qY[9] + qY[10] + qY[11] + qY[12] + qY[13] + qY[14] + qY[15] + qY[16] + qY[17] + qY[18] + qY[19] + qY[20] + qY[21] + qY[22] + qY[23] + qY[25]) + mu0Y) * D0Y + 2*lambda0Y*E0Y*D0Y + (qY[0]*D0A + qY[1]*D0B + qY[2]*D0C + qY[3]*D0D + qY[4]*D0E + qY[5]*D0F + qY[6]*D0G + qY[7]*D0H + qY[8]*D0I + qY[9]*D0J + qY[10]*D0K + qY[11]*D0L + qY[12]*D0M + qY[13]*D0N + qY[14]*D0O + qY[15]*D0P + qY[16]*D0Q + qY[17]*D0R + qY[18]*D0S + qY[19]*D0T + qY[20]*D0U + qY[21]*D0V + qY[22]*D0W + qY[23]*D0X + qY[25]*D0Z);
    
    ydot[51] =  - (lambda0Z + psi + (qZ[0] + qZ[1] + qZ[2] + qZ[3] + qZ[4] + qZ[5] + qZ[6] + qZ[7] + qZ[8] + qZ[9] + qZ[10] + qZ[11] + qZ[12] + qZ[13] + qZ[14] + qZ[15] + qZ[16] + qZ[17] + qZ[18] + qZ[19] + qZ[20] + qZ[21] + qZ[22] + qZ[23] + qZ[24]) + mu0Z) * D0Z + 2*lambda0Z*E0Z*D0Z + (qZ[0]*D0A + qZ[1]*D0B + qZ[2]*D0C + qZ[3]*D0D + qZ[4]*D0E + qZ[5]*D0F + qZ[6]*D0G + qZ[7]*D0H + qZ[8]*D0I + qZ[9]*D0J + qZ[10]*D0K + qZ[11]*D0L + qZ[12]*D0M + qZ[13]*D0N + qZ[14]*D0O + qZ[15]*D0P + qZ[16]*D0Q + qZ[17]*D0R + qZ[18]*D0S + qZ[19]*D0T + qZ[20]*D0U + qZ[21]*D0V + qZ[22]*D0W + qZ[23]*D0X + qZ[24]*D0Y);
    
}


void misse_strat_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
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
    mu0Z = params_misse[51],
    
    rate = params_misse[52];
    
    int hidden_states = params_misse[53];
    
    double qA[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qA, 26, hidden_states, 0, rate);
    double qB[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qB, 26, hidden_states, 1, rate);
    double qC[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qC, 26, hidden_states, 2, rate);
    double qD[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qD, 26, hidden_states, 3, rate);
    double qE[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qE, 26, hidden_states, 4, rate);
    double qF[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qF, 26, hidden_states, 5, rate);
    double qG[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qG, 26, hidden_states, 6, rate);
    double qH[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qH, 26, hidden_states, 7, rate);
    double qI[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qI, 26, hidden_states, 8, rate);
    double qJ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qJ, 26, hidden_states, 9, rate);
    double qK[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qK, 26, hidden_states, 10, rate);
    double qL[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qL, 26, hidden_states, 11, rate);
    double qM[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qM, 26, hidden_states, 12, rate);
    double qN[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qN, 26, hidden_states, 13, rate);
    double qO[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qO, 26, hidden_states, 14, rate);
    double qP[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qP, 26, hidden_states, 15, rate);
    double qQ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qQ, 26, hidden_states, 16, rate);
    double qR[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qR, 26, hidden_states, 17, rate);
    double qS[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qS, 26, hidden_states, 18, rate);
    double qT[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qT, 26, hidden_states, 19, rate);
    double qU[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qU, 26, hidden_states, 20, rate);
    double qV[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qV, 26, hidden_states, 21, rate);
    double qW[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qW, 26, hidden_states, 22, rate);
    double qX[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qX, 26, hidden_states, 23, rate);
    double qY[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qY, 26, hidden_states, 24, rate);
    double qZ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    change_rate(qZ, 26, hidden_states, 25, rate);
    
    double psi = params_misse[54];
    
    
    /* The E's */
    ydot[0] = mu0A - (lambda0A + psi + (qA[1] + qA[2] + qA[3] + qA[4] + qA[5] + qA[6] + qA[7] + qA[8] + qA[9] + qA[10] + qA[11] + qA[12] + qA[13] + qA[14] + qA[15] + qA[16] + qA[17] + qA[18] + qA[19] + qA[20] + qA[21] + qA[22] + qA[23] + qA[24] + qA[25]) + mu0A) * E0A + lambda0A*E0A*E0A + (qA[1]*E0B + qA[2]*E0C + qA[3]*E0D + qA[4]*E0E + qA[5]*E0F + qA[6]*E0G + qA[7]*E0H + qA[8]*E0I + qA[9]*E0J + qA[10]*E0K + qA[11]*E0L + qA[12]*E0M + qA[13]*E0N + qA[14]*E0O + qA[15]*E0P + qA[16]*E0Q + qA[17]*E0R + qA[18]*E0S + qA[19]*E0T + qA[20]*E0U + qA[21]*E0V + qA[22]*E0W + qA[23]*E0X + qA[24]*E0Y + qA[25]*E0Z);
    
    ydot[1] = mu0B - (lambda0B + psi + (qB[0] + qB[2] + qB[3] + qB[4] + qB[5] + qB[6] + qB[7] + qB[8] + qB[9] + qB[10] + qB[11] + qB[12] + qB[13] + qB[14] + qB[15] + qB[16] + qB[17] + qB[18] + qB[19] + qB[20] + qB[21] + qB[22] + qB[23] + qB[24] + qB[25]) + mu0B) * E0B + lambda0B*E0B*E0B + (qB[0]*E0A + qB[2]*E0C + qB[3]*E0D + qB[4]*E0E + qB[5]*E0F + qB[6]*E0G + qB[7]*E0H + qB[8]*E0I + qB[9]*E0J + qB[10]*E0K + qB[11]*E0L + qB[12]*E0M + qB[13]*E0N + qB[14]*E0O + qB[15]*E0P + qB[16]*E0Q + qB[17]*E0R + qB[18]*E0S + qB[19]*E0T + qB[20]*E0U + qB[21]*E0V + qB[22]*E0W + qB[23]*E0X + qB[24]*E0Y + qB[25]*E0Z);
    
    ydot[2] = mu0C - (lambda0C + psi + (qC[0] + qC[1] + qC[3] + qC[4] + qC[5] + qC[6] + qC[7] + qC[8] + qC[9] + qC[10] + qC[11] + qC[12] + qC[13] + qC[14] + qC[15] + qC[16] + qC[17] + qC[18] + qC[19] + qC[20] + qC[21] + qC[22] + qC[23] + qC[24] + qC[25]) + mu0C) * E0C + lambda0C*E0C*E0C + (qC[0]*E0A + qC[1]*E0B + qC[3]*E0D + qC[4]*E0E + qC[5]*E0F + qC[6]*E0G + qC[7]*E0H + qC[8]*E0I + qC[9]*E0J + qC[10]*E0K + qC[11]*E0L + qC[12]*E0M + qC[13]*E0N + qC[14]*E0O + qC[15]*E0P + qC[16]*E0Q + qC[17]*E0R + qC[18]*E0S + qC[19]*E0T + qC[20]*E0U + qC[21]*E0V + qC[22]*E0W + qC[23]*E0X + qC[24]*E0Y + qC[25]*E0Z);
    
    ydot[3] = mu0D - (lambda0D + psi + (qD[0] + qD[1] + qD[2] + qD[4] + qD[5] + qD[6] + qD[7] + qD[8] + qD[9] + qD[10] + qD[11] + qD[12] + qD[13] + qD[14] + qD[15] + qD[16] + qD[17] + qD[18] + qD[19] + qD[20] + qD[21] + qD[22] + qD[23] + qD[24] + qD[25]) + mu0D) * E0D + lambda0D*E0D*E0D + (qD[0]*E0A + qD[1]*E0B + qD[2]*E0C + qD[4]*E0E + qD[5]*E0F + qD[6]*E0G + qD[7]*E0H + qD[8]*E0I + qD[9]*E0J + qD[10]*E0K + qD[11]*E0L + qD[12]*E0M + qD[13]*E0N + qD[14]*E0O + qD[15]*E0P + qD[16]*E0Q + qD[17]*E0R + qD[18]*E0S + qD[19]*E0T + qD[20]*E0U + qD[21]*E0V + qD[22]*E0W + qD[23]*E0X + qD[24]*E0Y + qD[25]*E0Z);
    
    ydot[4] = mu0E - (lambda0E + psi + (qE[0] + qE[1] + qE[2] + qE[3] + qE[5] + qE[6] + qE[7] + qE[8] + qE[9] + qE[10] + qE[11] + qE[12] + qE[13] + qE[14] + qE[15] + qE[16] + qE[17] + qE[18] + qE[19] + qE[20] + qE[21] + qE[22] + qE[23] + qE[24] + qE[25]) + mu0E) * E0E + lambda0E*E0E*E0E + (qE[0]*E0A + qE[1]*E0B + qE[2]*E0C + qE[3]*E0D + qE[5]*E0F + qE[6]*E0G + qE[7]*E0H + qE[8]*E0I + qE[9]*E0J + qE[10]*E0K + qE[11]*E0L + qE[12]*E0M + qE[13]*E0N + qE[14]*E0O + qE[15]*E0P + qE[16]*E0Q + qE[17]*E0R + qE[18]*E0S + qE[19]*E0T + qE[20]*E0U + qE[21]*E0V + qE[22]*E0W + qE[23]*E0X + qE[24]*E0Y + qE[25]*E0Z);
    
    ydot[5] = mu0F - (lambda0F + psi + (qF[0] + qF[1] + qF[2] + qF[3] + qF[4] + qF[6] + qF[7] + qF[8] + qF[9] + qF[10] + qF[11] + qF[12] + qF[13] + qF[14] + qF[15] + qF[16] + qF[17] + qF[18] + qF[19] + qF[20] + qF[21] + qF[22] + qF[23] + qF[24] + qF[25]) + mu0F) * E0F + lambda0F*E0F*E0F + (qF[0]*E0A + qF[1]*E0B + qF[2]*E0C + qF[3]*E0D + qF[4]*E0E + qF[6]*E0G + qF[7]*E0H + qF[8]*E0I + qF[9]*E0J + qF[10]*E0K + qF[11]*E0L + qF[12]*E0M + qF[13]*E0N + qF[14]*E0O + qF[15]*E0P + qF[16]*E0Q + qF[17]*E0R + qF[18]*E0S + qF[19]*E0T + qF[20]*E0U + qF[21]*E0V + qF[22]*E0W + qF[23]*E0X + qF[24]*E0Y + qF[25]*E0Z);
    
    ydot[6] = mu0G - (lambda0G + psi + (qG[0] + qG[1] + qG[2] + qG[3] + qG[4] + qG[5] + qG[7] + qG[8] + qG[9] + qG[10] + qG[11] + qG[12] + qG[13] + qG[14] + qG[15] + qG[16] + qG[17] + qG[18] + qG[19] + qG[20] + qG[21] + qG[22] + qG[23] + qG[24] + qG[25]) + mu0G) * E0G + lambda0G*E0G*E0G + (qG[0]*E0A + qG[1]*E0B + qG[2]*E0C + qG[3]*E0D + qG[4]*E0E + qG[5]*E0F + qG[7]*E0H + qG[8]*E0I + qG[9]*E0J + qG[10]*E0K + qG[11]*E0L + qG[12]*E0M + qG[13]*E0N + qG[14]*E0O + qG[15]*E0P + qG[16]*E0Q + qG[17]*E0R + qG[18]*E0S + qG[19]*E0T + qG[20]*E0U + qG[21]*E0V + qG[22]*E0W + qG[23]*E0X + qG[24]*E0Y + qG[25]*E0Z);
    
    ydot[7] = mu0H - (lambda0H + psi + (qH[0] + qH[1] + qH[2] + qH[3] + qH[4] + qH[5] + qH[6] + qH[8] + qH[9] + qH[10] + qH[11] + qH[12] + qH[13] + qH[14] + qH[15] + qH[16] + qH[17] + qH[18] + qH[19] + qH[20] + qH[21] + qH[22] + qH[23] + qH[24] + qH[25]) + mu0H) * E0H + lambda0H*E0H*E0H + (qH[0]*E0A + qH[1]*E0B + qH[2]*E0C + qH[3]*E0D + qH[4]*E0E + qH[5]*E0F + qH[6]*E0G + qH[8]*E0I + qH[9]*E0J + qH[10]*E0K + qH[11]*E0L + qH[12]*E0M + qH[13]*E0N + qH[14]*E0O + qH[15]*E0P + qH[16]*E0Q + qH[17]*E0R + qH[18]*E0S + qH[19]*E0T + qH[20]*E0U + qH[21]*E0V + qH[22]*E0W + qH[23]*E0X + qH[24]*E0Y + qH[25]*E0Z);
    
    ydot[8] = mu0I - (lambda0I + psi + (qI[0] + qI[1] + qI[2] + qI[3] + qI[4] + qI[5] + qI[6] + qI[7] + qI[9] + qI[10] + qI[11] + qI[12] + qI[13] + qI[14] + qI[15] + qI[16] + qI[17] + qI[18] + qI[19] + qI[20] + qI[21] + qI[22] + qI[23] + qI[24] + qI[25]) + mu0I) * E0I + lambda0I*E0I*E0I + (qI[0]*E0A + qI[1]*E0B + qI[2]*E0C + qI[3]*E0D + qI[4]*E0E + qI[5]*E0F + qI[6]*E0G + qI[7]*E0H + qI[9]*E0J + qI[10]*E0K + qI[11]*E0L + qI[12]*E0M + qI[13]*E0N + qI[14]*E0O + qI[15]*E0P + qI[16]*E0Q + qI[17]*E0R + qI[18]*E0S + qI[19]*E0T + qI[20]*E0U + qI[21]*E0V + qI[22]*E0W + qI[23]*E0X + qI[24]*E0Y + qI[25]*E0Z);
    
    ydot[9] = mu0J - (lambda0J + psi + (qJ[0] + qJ[1] + qJ[2] + qJ[3] + qJ[4] + qJ[5] + qJ[6] + qJ[7] + qJ[8] + qJ[10] + qJ[11] + qJ[12] + qJ[13] + qJ[14] + qJ[15] + qJ[16] + qJ[17] + qJ[18] + qJ[19] + qJ[20] + qJ[21] + qJ[22] + qJ[23] + qJ[24] + qJ[25]) + mu0J) * E0J + lambda0J*E0J*E0J + (qJ[0]*E0A + qJ[1]*E0B + qJ[2]*E0C + qJ[3]*E0D + qJ[4]*E0E + qJ[5]*E0F + qJ[6]*E0G + qJ[7]*E0H + qJ[8]*E0I + qJ[10]*E0K + qJ[11]*E0L + qJ[12]*E0M + qJ[13]*E0N + qJ[14]*E0O + qJ[15]*E0P + qJ[16]*E0Q + qJ[17]*E0R + qJ[18]*E0S + qJ[19]*E0T + qJ[20]*E0U + qJ[21]*E0V + qJ[22]*E0W + qJ[23]*E0X + qJ[24]*E0Y + qJ[25]*E0Z);
    
    ydot[10] = mu0K - (lambda0K + psi + (qK[0] + qK[1] + qK[2] + qK[3] + qK[4] + qK[5] + qK[6] + qK[7] + qK[8] + qK[9] + qK[11] + qK[12] + qK[13] + qK[14] + qK[15] + qK[16] + qK[17] + qK[18] + qK[19] + qK[20] + qK[21] + qK[22] + qK[23] + qK[24] + qK[25]) + mu0K) * E0K + lambda0K*E0K*E0K + (qK[0]*E0A + qK[1]*E0B + qK[2]*E0C + qK[3]*E0D + qK[4]*E0E + qK[5]*E0F + qK[6]*E0G + qK[7]*E0H + qK[8]*E0I + qK[9]*E0J + qK[11]*E0L + qK[12]*E0M + qK[13]*E0N + qK[14]*E0O + qK[15]*E0P + qK[16]*E0Q + qK[17]*E0R + qK[18]*E0S + qK[19]*E0T + qK[20]*E0U + qK[21]*E0V + qK[22]*E0W + qK[23]*E0X + qK[24]*E0Y + qK[25]*E0Z);
    
    ydot[11] = mu0L - (lambda0L + psi + (qL[0] + qL[1] + qL[2] + qL[3] + qL[4] + qL[5] + qL[6] + qL[7] + qL[8] + qL[9] + qL[10] + qL[12] + qL[13] + qL[14] + qL[15] + qL[16] + qL[17] + qL[18] + qL[19] + qL[20] + qL[21] + qL[22] + qL[23] + qL[24] + qL[25]) + mu0L) * E0L + lambda0L*E0L*E0L + (qL[0]*E0A + qL[1]*E0B + qL[2]*E0C + qL[3]*E0D + qL[4]*E0E + qL[5]*E0F + qL[6]*E0G + qL[7]*E0H + qL[8]*E0I + qL[9]*E0J + qL[10]*E0K + qL[12]*E0M + qL[13]*E0N + qL[14]*E0O + qL[15]*E0P + qL[16]*E0Q + qL[17]*E0R + qL[18]*E0S + qL[19]*E0T + qL[20]*E0U + qL[21]*E0V + qL[22]*E0W + qL[23]*E0X + qL[24]*E0Y + qL[25]*E0Z);
    
    ydot[12] = mu0M - (lambda0M + psi + (qM[0] + qM[1] + qM[2] + qM[3] + qM[4] + qM[5] + qM[6] + qM[7] + qM[8] + qM[9] + qM[10] + qM[11] + qM[13] + qM[14] + qM[15] + qM[16] + qM[17] + qM[18] + qM[19] + qM[20] + qM[21] + qM[22] + qM[23] + qM[24] + qM[25]) + mu0M) * E0M + lambda0M*E0M*E0M + (qM[0]*E0A + qM[1]*E0B + qM[2]*E0C + qM[3]*E0D + qM[4]*E0E + qM[5]*E0F + qM[6]*E0G + qM[7]*E0H + qM[8]*E0I + qM[9]*E0J + qM[10]*E0K + qM[11]*E0L + qM[13]*E0N + qM[14]*E0O + qM[15]*E0P + qM[16]*E0Q + qM[17]*E0R + qM[18]*E0S + qM[19]*E0T + qM[20]*E0U + qM[21]*E0V + qM[22]*E0W + qM[23]*E0X + qM[24]*E0Y + qM[25]*E0Z);
    
    ydot[13] = mu0N - (lambda0N + psi + (qN[0] + qN[1] + qN[2] + qN[3] + qN[4] + qN[5] + qN[6] + qN[7] + qN[8] + qN[9] + qN[10] + qN[11] + qN[12] + qN[14] + qN[15] + qN[16] + qN[17] + qN[18] + qN[19] + qN[20] + qN[21] + qN[22] + qN[23] + qN[24] + qN[25]) + mu0N) * E0N + lambda0N*E0N*E0N + (qN[0]*E0A + qN[1]*E0B + qN[2]*E0C + qN[3]*E0D + qN[4]*E0E + qN[5]*E0F + qN[6]*E0G + qN[7]*E0H + qN[8]*E0I + qN[9]*E0J + qN[10]*E0K + qN[11]*E0L + qN[12]*E0M + qN[14]*E0O + qN[15]*E0P + qN[16]*E0Q + qN[17]*E0R + qN[18]*E0S + qN[19]*E0T + qN[20]*E0U + qN[21]*E0V + qN[22]*E0W + qN[23]*E0X + qN[24]*E0Y + qN[25]*E0Z);
    
    ydot[14] = mu0O - (lambda0O + psi + (qO[0] + qO[1] + qO[2] + qO[3] + qO[4] + qO[5] + qO[6] + qO[7] + qO[8] + qO[9] + qO[10] + qO[11] + qO[12] + qO[13] + qO[15] + qO[16] + qO[17] + qO[18] + qO[19] + qO[20] + qO[21] + qO[22] + qO[23] + qO[24] + qO[25]) + mu0O) * E0O + lambda0O*E0O*E0O + (qO[0]*E0A + qO[1]*E0B + qO[2]*E0C + qO[3]*E0D + qO[4]*E0E + qO[5]*E0F + qO[6]*E0G + qO[7]*E0H + qO[8]*E0I + qO[9]*E0J + qO[10]*E0K + qO[11]*E0L + qO[12]*E0M + qO[13]*E0N + qO[15]*E0P + qO[16]*E0Q + qO[17]*E0R + qO[18]*E0S + qO[19]*E0T + qO[20]*E0U + qO[21]*E0V + qO[22]*E0W + qO[23]*E0X + qO[24]*E0Y + qO[25]*E0Z);
    
    ydot[15] = mu0P - (lambda0P + psi + (qP[0] + qP[1] + qP[2] + qP[3] + qP[4] + qP[5] + qP[6] + qP[7] + qP[8] + qP[9] + qP[10] + qP[11] + qP[12] + qP[13] + qP[14] + qP[16] + qP[17] + qP[18] + qP[19] + qP[20] + qP[21] + qP[22] + qP[23] + qP[24] + qP[25]) + mu0P) * E0P + lambda0P*E0P*E0P + (qP[0]*E0A + qP[1]*E0B + qP[2]*E0C + qP[3]*E0D + qP[4]*E0E + qP[5]*E0F + qP[6]*E0G + qP[7]*E0H + qP[8]*E0I + qP[9]*E0J + qP[10]*E0K + qP[11]*E0L + qP[12]*E0M + qP[13]*E0N + qP[14]*E0O + qP[16]*E0Q + qP[17]*E0R + qP[18]*E0S + qP[19]*E0T + qP[20]*E0U + qP[21]*E0V + qP[22]*E0W + qP[23]*E0X + qP[24]*E0Y + qP[25]*E0Z);
    
    ydot[16] = mu0Q - (lambda0Q + psi + (qQ[0] + qQ[1] + qQ[2] + qQ[3] + qQ[4] + qQ[5] + qQ[6] + qQ[7] + qQ[8] + qQ[9] + qQ[10] + qQ[11] + qQ[12] + qQ[13] + qQ[14] + qQ[15] + qQ[17] + qQ[18] + qQ[19] + qQ[20] + qQ[21] + qQ[22] + qQ[23] + qQ[24] + qQ[25]) + mu0Q) * E0Q + lambda0Q*E0Q*E0Q + (qQ[0]*E0A + qQ[1]*E0B + qQ[2]*E0C + qQ[3]*E0D + qQ[4]*E0E + qQ[5]*E0F + qQ[6]*E0G + qQ[7]*E0H + qQ[8]*E0I + qQ[9]*E0J + qQ[10]*E0K + qQ[11]*E0L + qQ[12]*E0M + qQ[13]*E0N + qQ[14]*E0O + qQ[15]*E0P + qQ[17]*E0R + qQ[18]*E0S + qQ[19]*E0T + qQ[20]*E0U + qQ[21]*E0V + qQ[22]*E0W + qQ[23]*E0X + qQ[24]*E0Y + qQ[25]*E0Z);
    
    ydot[17] = mu0R - (lambda0R + psi + (qR[0] + qR[1] + qR[2] + qR[3] + qR[4] + qR[5] + qR[6] + qR[7] + qR[8] + qR[9] + qR[10] + qR[11] + qR[12] + qR[13] + qR[14] + qR[15] + qR[16] + qR[18] + qR[19] + qR[20] + qR[21] + qR[22] + qR[23] + qR[24] + qR[25]) + mu0R) * E0R + lambda0R*E0R*E0R + (qR[0]*E0A + qR[1]*E0B + qR[2]*E0C + qR[3]*E0D + qR[4]*E0E + qR[5]*E0F + qR[6]*E0G + qR[7]*E0H + qR[8]*E0I + qR[9]*E0J + qR[10]*E0K + qR[11]*E0L + qR[12]*E0M + qR[13]*E0N + qR[14]*E0O + qR[15]*E0P + qR[16]*E0Q + qR[18]*E0S + qR[19]*E0T + qR[20]*E0U + qR[21]*E0V + qR[22]*E0W + qR[23]*E0X + qR[24]*E0Y + qR[25]*E0Z);
    
    ydot[18] = mu0S - (lambda0S + psi + (qS[0] + qS[1] + qS[2] + qS[3] + qS[4] + qS[5] + qS[6] + qS[7] + qS[8] + qS[9] + qS[10] + qS[11] + qS[12] + qS[13] + qS[14] + qS[15] + qS[16] + qS[17] + qS[19] + qS[20] + qS[21] + qS[22] + qS[23] + qS[24] + qS[25]) + mu0S) * E0S + lambda0S*E0S*E0S + (qS[0]*E0A + qS[1]*E0B + qS[2]*E0C + qS[3]*E0D + qS[4]*E0E + qS[5]*E0F + qS[6]*E0G + qS[7]*E0H + qS[8]*E0I + qS[9]*E0J + qS[10]*E0K + qS[11]*E0L + qS[12]*E0M + qS[13]*E0N + qS[14]*E0O + qS[15]*E0P + qS[16]*E0Q + qS[17]*E0R + qS[19]*E0T + qS[20]*E0U + qS[21]*E0V + qS[22]*E0W + qS[23]*E0X + qS[24]*E0Y + qS[25]*E0Z);
    
    ydot[19] = mu0T - (lambda0T + psi + (qT[0] + qT[1] + qT[2] + qT[3] + qT[4] + qT[5] + qT[6] + qT[7] + qT[8] + qT[9] + qT[10] + qT[11] + qT[12] + qT[13] + qT[14] + qT[15] + qT[16] + qT[17] + qT[18] + qT[20] + qT[21] + qT[22] + qT[23] + qT[24] + qT[25]) + mu0T) * E0T + lambda0T*E0T*E0T + (qT[0]*E0A + qT[1]*E0B + qT[2]*E0C + qT[3]*E0D + qT[4]*E0E + qT[5]*E0F + qT[6]*E0G + qT[7]*E0H + qT[8]*E0I + qT[9]*E0J + qT[10]*E0K + qT[11]*E0L + qT[12]*E0M + qT[13]*E0N + qT[14]*E0O + qT[15]*E0P + qT[16]*E0Q + qT[17]*E0R + qT[18]*E0S + qT[20]*E0U + qT[21]*E0V + qT[22]*E0W + qT[23]*E0X + qT[24]*E0Y + qT[25]*E0Z);
    
    ydot[20] = mu0U - (lambda0U + psi + (qU[0] + qU[1] + qU[2] + qU[3] + qU[4] + qU[5] + qU[6] + qU[7] + qU[8] + qU[9] + qU[10] + qU[11] + qU[12] + qU[13] + qU[14] + qU[15] + qU[16] + qU[17] + qU[18] + qU[19] + qU[21] + qU[22] + qU[23] + qU[24] + qU[25]) + mu0U) * E0U + lambda0U*E0U*E0U + (qU[0]*E0A + qU[1]*E0B + qU[2]*E0C + qU[3]*E0D + qU[4]*E0E + qU[5]*E0F + qU[6]*E0G + qU[7]*E0H + qU[8]*E0I + qU[9]*E0J + qU[10]*E0K + qU[11]*E0L + qU[12]*E0M + qU[13]*E0N + qU[14]*E0O + qU[15]*E0P + qU[16]*E0Q + qU[17]*E0R + qU[18]*E0S + qU[19]*E0T + qU[21]*E0V + qU[22]*E0W + qU[23]*E0X + qU[24]*E0Y + qU[25]*E0Z);
    
    ydot[21] = mu0V - (lambda0V + psi + (qV[0] + qV[1] + qV[2] + qV[3] + qV[4] + qV[5] + qV[6] + qV[7] + qV[8] + qV[9] + qV[10] + qV[11] + qV[12] + qV[13] + qV[14] + qV[15] + qV[16] + qV[17] + qV[18] + qV[19] + qV[20] + qV[22] + qV[23] + qV[24] + qV[25]) + mu0V) * E0V + lambda0V*E0V*E0V + (qV[0]*E0A + qV[1]*E0B + qV[2]*E0C + qV[3]*E0D + qV[4]*E0E + qV[5]*E0F + qV[6]*E0G + qV[7]*E0H + qV[8]*E0I + qV[9]*E0J + qV[10]*E0K + qV[11]*E0L + qV[12]*E0M + qV[13]*E0N + qV[14]*E0O + qV[15]*E0P + qV[16]*E0Q + qV[17]*E0R + qV[18]*E0S + qV[19]*E0T + qV[20]*E0U + qV[22]*E0W + qV[23]*E0X + qV[24]*E0Y + qV[25]*E0Z);
    
    ydot[22] = mu0W - (lambda0W + psi + (qW[0] + qW[1] + qW[2] + qW[3] + qW[4] + qW[5] + qW[6] + qW[7] + qW[8] + qW[9] + qW[10] + qW[11] + qW[12] + qW[13] + qW[14] + qW[15] + qW[16] + qW[17] + qW[18] + qW[19] + qW[20] + qW[21] + qW[23] + qW[24] + qW[25]) + mu0W) * E0W + lambda0W*E0W*E0W + (qW[0]*E0A + qW[1]*E0B + qW[2]*E0C + qW[3]*E0D + qW[4]*E0E + qW[5]*E0F + qW[6]*E0G + qW[7]*E0H + qW[8]*E0I + qW[9]*E0J + qW[10]*E0K + qW[11]*E0L + qW[12]*E0M + qW[13]*E0N + qW[14]*E0O + qW[15]*E0P + qW[16]*E0Q + qW[17]*E0R + qW[18]*E0S + qW[19]*E0T + qW[20]*E0U + qW[21]*E0V + qW[23]*E0X + qW[24]*E0Y + qW[25]*E0Z);
    
    ydot[23] = mu0X - (lambda0X + psi + (qX[0] + qX[1] + qX[2] + qX[3] + qX[4] + qX[5] + qX[6] + qX[7] + qX[8] + qX[9] + qX[10] + qX[11] + qX[12] + qX[13] + qX[14] + qX[15] + qX[16] + qX[17] + qX[18] + qX[19] + qX[20] + qX[21] + qX[22] + qX[24] + qX[25]) + mu0X) * E0X + lambda0X*E0X*E0X + (qX[0]*E0A + qX[1]*E0B + qX[2]*E0C + qX[3]*E0D + qX[4]*E0E + qX[5]*E0F + qX[6]*E0G + qX[7]*E0H + qX[8]*E0I + qX[9]*E0J + qX[10]*E0K + qX[11]*E0L + qX[12]*E0M + qX[13]*E0N + qX[14]*E0O + qX[15]*E0P + qX[16]*E0Q + qX[17]*E0R + qX[18]*E0S + qX[19]*E0T + qX[20]*E0U + qX[21]*E0V + qX[22]*E0W + qX[24]*E0Y + qX[25]*E0Z);
    
    ydot[24] = mu0Y - (lambda0Y + psi + (qY[0] + qY[1] + qY[2] + qY[3] + qY[4] + qY[5] + qY[6] + qY[7] + qY[8] + qY[9] + qY[10] + qY[11] + qY[12] + qY[13] + qY[14] + qY[15] + qY[16] + qY[17] + qY[18] + qY[19] + qY[20] + qY[21] + qY[22] + qY[23] + qY[25]) + mu0Y) * E0Y + lambda0Y*E0Y*E0Y + (qY[0]*E0A + qY[1]*E0B + qY[2]*E0C + qY[3]*E0D + qY[4]*E0E + qY[5]*E0F + qY[6]*E0G + qY[7]*E0H + qY[8]*E0I + qY[9]*E0J + qY[10]*E0K + qY[11]*E0L + qY[12]*E0M + qY[13]*E0N + qY[14]*E0O + qY[15]*E0P + qY[16]*E0Q + qY[17]*E0R + qY[18]*E0S + qY[19]*E0T + qY[20]*E0U + qY[21]*E0V + qY[22]*E0W + qY[23]*E0X + qY[25]*E0Z);
    
    ydot[25] = mu0Z - (lambda0Z + psi + (qZ[0] + qZ[1] + qZ[2] + qZ[3] + qZ[4] + qZ[5] + qZ[6] + qZ[7] + qZ[8] + qZ[9] + qZ[10] + qZ[11] + qZ[12] + qZ[13] + qZ[14] + qZ[15] + qZ[16] + qZ[17] + qZ[18] + qZ[19] + qZ[20] + qZ[21] + qZ[22] + qZ[23] + qZ[24]) + mu0Z) * E0Z + lambda0Z*E0Z*E0Z + (qZ[0]*E0A + qZ[1]*E0B + qZ[2]*E0C + qZ[3]*E0D + qZ[4]*E0E + qZ[5]*E0F + qZ[6]*E0G + qZ[7]*E0H + qZ[8]*E0I + qZ[9]*E0J + qZ[10]*E0K + qZ[11]*E0L + qZ[12]*E0M + qZ[13]*E0N + qZ[14]*E0O + qZ[15]*E0P + qZ[16]*E0Q + qZ[17]*E0R + qZ[18]*E0S + qZ[19]*E0T + qZ[20]*E0U + qZ[21]*E0V + qZ[22]*E0W + qZ[23]*E0X + qZ[24]*E0Y);
    
    
    /* The D's */
    ydot[26] =  - (lambda0A + psi + (qA[1] + qA[2] + qA[3] + qA[4] + qA[5] + qA[6] + qA[7] + qA[8] + qA[9] + qA[10] + qA[11] + qA[12] + qA[13] + qA[14] + qA[15] + qA[16] + qA[17] + qA[18] + qA[19] + qA[20] + qA[21] + qA[22] + qA[23] + qA[24] + qA[25]) + mu0A) * D0A + (qA[1]*D0B + qA[2]*D0C + qA[3]*D0D + qA[4]*D0E + qA[5]*D0F + qA[6]*D0G + qA[7]*D0H + qA[8]*D0I + qA[9]*D0J + qA[10]*D0K + qA[11]*D0L + qA[12]*D0M + qA[13]*D0N + qA[14]*D0O + qA[15]*D0P + qA[16]*D0Q + qA[17]*D0R + qA[18]*D0S + qA[19]*D0T + qA[20]*D0U + qA[21]*D0V + qA[22]*D0W + qA[23]*D0X + qA[24]*D0Y + qA[25]*D0Z);
    
    ydot[27] =  - (lambda0B + psi + (qB[0] + qB[2] + qB[3] + qB[4] + qB[5] + qB[6] + qB[7] + qB[8] + qB[9] + qB[10] + qB[11] + qB[12] + qB[13] + qB[14] + qB[15] + qB[16] + qB[17] + qB[18] + qB[19] + qB[20] + qB[21] + qB[22] + qB[23] + qB[24] + qB[25]) + mu0B) * D0B + (qB[0]*D0A + qB[2]*D0C + qB[3]*D0D + qB[4]*D0E + qB[5]*D0F + qB[6]*D0G + qB[7]*D0H + qB[8]*D0I + qB[9]*D0J + qB[10]*D0K + qB[11]*D0L + qB[12]*D0M + qB[13]*D0N + qB[14]*D0O + qB[15]*D0P + qB[16]*D0Q + qB[17]*D0R + qB[18]*D0S + qB[19]*D0T + qB[20]*D0U + qB[21]*D0V + qB[22]*D0W + qB[23]*D0X + qB[24]*D0Y + qB[25]*D0Z);
    
    ydot[28] =  - (lambda0C + psi + (qC[0] + qC[1] + qC[3] + qC[4] + qC[5] + qC[6] + qC[7] + qC[8] + qC[9] + qC[10] + qC[11] + qC[12] + qC[13] + qC[14] + qC[15] + qC[16] + qC[17] + qC[18] + qC[19] + qC[20] + qC[21] + qC[22] + qC[23] + qC[24] + qC[25]) + mu0C) * D0C + (qC[0]*D0A + qC[1]*D0B + qC[3]*D0D + qC[4]*D0E + qC[5]*D0F + qC[6]*D0G + qC[7]*D0H + qC[8]*D0I + qC[9]*D0J + qC[10]*D0K + qC[11]*D0L + qC[12]*D0M + qC[13]*D0N + qC[14]*D0O + qC[15]*D0P + qC[16]*D0Q + qC[17]*D0R + qC[18]*D0S + qC[19]*D0T + qC[20]*D0U + qC[21]*D0V + qC[22]*D0W + qC[23]*D0X + qC[24]*D0Y + qC[25]*D0Z);
    
    ydot[29] =  - (lambda0D + psi + (qD[0] + qD[1] + qD[2] + qD[4] + qD[5] + qD[6] + qD[7] + qD[8] + qD[9] + qD[10] + qD[11] + qD[12] + qD[13] + qD[14] + qD[15] + qD[16] + qD[17] + qD[18] + qD[19] + qD[20] + qD[21] + qD[22] + qD[23] + qD[24] + qD[25])+ mu0D) * D0D + (qD[0]*D0A + qD[1]*D0B + qD[2]*D0C + qD[4]*D0E + qD[5]*D0F + qD[6]*D0G + qD[7]*D0H + qD[8]*D0I + qD[9]*D0J + qD[10]*D0K + qD[11]*D0L + qD[12]*D0M + qD[13]*D0N + qD[14]*D0O + qD[15]*D0P + qD[16]*D0Q + qD[17]*D0R + qD[18]*D0S + qD[19]*D0T + qD[20]*D0U + qD[21]*D0V + qD[22]*D0W + qD[23]*D0X + qD[24]*D0Y + qD[25]*D0Z);
    
    ydot[30] =  - (lambda0E + psi + (qE[0] + qE[1] + qE[2] + qE[3] + qE[5] + qE[6] + qE[7] + qE[8] + qE[9] + qE[10] + qE[11] + qE[12] + qE[13] + qE[14] + qE[15] + qE[16] + qE[17] + qE[18] + qE[19] + qE[20] + qE[21] + qE[22] + qE[23] + qE[24] + qE[25])+ mu0E) * D0E + (qE[0]*D0A + qE[1]*D0B + qE[2]*D0C + qE[3]*D0D + qE[5]*D0F + qE[6]*D0G + qE[7]*D0H + qE[8]*D0I + qE[9]*D0J + qE[10]*D0K + qE[11]*D0L + qE[12]*D0M + qE[13]*D0N + qE[14]*D0O + qE[15]*D0P + qE[16]*D0Q + qE[17]*D0R + qE[18]*D0S + qE[19]*D0T + qE[20]*D0U + qE[21]*D0V + qE[22]*D0W + qE[23]*D0X + qE[24]*D0Y + qE[25]*D0Z);
    
    ydot[31] =  - (lambda0F + psi + (qF[0] + qF[1] + qF[2] + qF[3] + qF[4] + qF[6] + qF[7] + qF[8] + qF[9] + qF[10] + qF[11] + qF[12] + qF[13] + qF[14] + qF[15] + qF[16] + qF[17] + qF[18] + qF[19] + qF[20] + qF[21] + qF[22] + qF[23] + qF[24] + qF[25]) + mu0F) * D0F + (qF[0]*D0A + qF[1]*D0B + qF[2]*D0C + qF[3]*D0D + qF[4]*D0E + qF[6]*D0G + qF[7]*D0H + qF[8]*D0I + qF[9]*D0J + qF[10]*D0K + qF[11]*D0L + qF[12]*D0M + qF[13]*D0N + qF[14]*D0O + qF[15]*D0P + qF[16]*D0Q + qF[17]*D0R + qF[18]*D0S + qF[19]*D0T + qF[20]*D0U + qF[21]*D0V + qF[22]*D0W + qF[23]*D0X + qF[24]*D0Y + qF[25]*D0Z);
    
    ydot[32] =  - (lambda0G + psi + (qG[0] + qG[1] + qG[2] + qG[3] + qG[4] + qG[5] + qG[7] + qG[8] + qG[9] + qG[10] + qG[11] + qG[12] + qG[13] + qG[14] + qG[15] + qG[16] + qG[17] + qG[18] + qG[19] + qG[20] + qG[21] + qG[22] + qG[23] + qG[24] + qG[25]) + mu0G) * D0G + (qG[0]*D0A + qG[1]*D0B + qG[2]*D0C + qG[3]*D0D + qG[4]*D0E + qG[5]*D0F + qG[7]*D0H + qG[8]*D0I + qG[9]*D0J + qG[10]*D0K + qG[11]*D0L + qG[12]*D0M + qG[13]*D0N + qG[14]*D0O + qG[15]*D0P + qG[16]*D0Q + qG[17]*D0R + qG[18]*D0S + qG[19]*D0T + qG[20]*D0U + qG[21]*D0V + qG[22]*D0W + qG[23]*D0X + qG[24]*D0Y + qG[25]*D0Z);
    
    ydot[33] =  - (lambda0H + psi + (qH[0] + qH[1] + qH[2] + qH[3] + qH[4] + qH[5] + qH[6] + qH[8] + qH[9] + qH[10] + qH[11] + qH[12] + qH[13] + qH[14] + qH[15] + qH[16] + qH[17] + qH[18] + qH[19] + qH[20] + qH[21] + qH[22] + qH[23] + qH[24] + qH[25]) + mu0H) * D0H + (qH[0]*D0A + qH[1]*D0B + qH[2]*D0C + qH[3]*D0D + qH[4]*D0E + qH[5]*D0F + qH[6]*D0G + qH[8]*D0I + qH[9]*D0J + qH[10]*D0K + qH[11]*D0L + qH[12]*D0M + qH[13]*D0N + qH[14]*D0O + qH[15]*D0P + qH[16]*D0Q + qH[17]*D0R + qH[18]*D0S + qH[19]*D0T + qH[20]*D0U + qH[21]*D0V + qH[22]*D0W + qH[23]*D0X + qH[24]*D0Y + qH[25]*D0Z);
    
    ydot[34] =  - (lambda0I + psi + (qI[0] + qI[1] + qI[2] + qI[3] + qI[4] + qI[5] + qI[6] + qI[7] + qI[9] + qI[10] + qI[11] + qI[12] + qI[13] + qI[14] + qI[15] + qI[16] + qI[17] + qI[18] + qI[19] + qI[20] + qI[21] + qI[22] + qI[23] + qI[24] + qI[25]) + mu0I) * D0I + (qI[0]*D0A + qI[1]*D0B + qI[2]*D0C + qI[3]*D0D + qI[4]*D0E + qI[5]*D0F + qI[6]*D0G + qI[7]*D0H + qI[9]*D0J + qI[10]*D0K + qI[11]*D0L + qI[12]*D0M + qI[13]*D0N + qI[14]*D0O + qI[15]*D0P + qI[16]*D0Q + qI[17]*D0R + qI[18]*D0S + qI[19]*D0T + qI[20]*D0U + qI[21]*D0V + qI[22]*D0W + qI[23]*D0X + qI[24]*D0Y + qI[25]*D0Z);
    
    ydot[35] =  - (lambda0J + psi + (qJ[0] + qJ[1] + qJ[2] + qJ[3] + qJ[4] + qJ[5] + qJ[6] + qJ[7] + qJ[8] + qJ[10] + qJ[11] + qJ[12] + qJ[13] + qJ[14] + qJ[15] + qJ[16] + qJ[17] + qJ[18] + qJ[19] + qJ[20] + qJ[21] + qJ[22] + qJ[23] + qJ[24] + qJ[25]) + mu0J) * D0J + (qJ[0]*D0A + qJ[1]*D0B + qJ[2]*D0C + qJ[3]*D0D + qJ[4]*D0E + qJ[5]*D0F + qJ[6]*D0G + qJ[7]*D0H + qJ[8]*D0I + qJ[10]*D0K + qJ[11]*D0L + qJ[12]*D0M + qJ[13]*D0N + qJ[14]*D0O + qJ[15]*D0P + qJ[16]*D0Q + qJ[17]*D0R + qJ[18]*D0S + qJ[19]*D0T + qJ[20]*D0U + qJ[21]*D0V + qJ[22]*D0W + qJ[23]*D0X + qJ[24]*D0Y + qJ[25]*D0Z);
    
    ydot[36] =  - (lambda0K + psi + (qK[0] + qK[1] + qK[2] + qK[3] + qK[4] + qK[5] + qK[6] + qK[7] + qK[8] + qK[9] + qK[11] + qK[12] + qK[13] + qK[14] + qK[15] + qK[16] + qK[17] + qK[18] + qK[19] + qK[20] + qK[21] + qK[22] + qK[23] + qK[24] + qK[25]) + mu0K) * D0K + (qK[0]*D0A + qK[1]*D0B + qK[2]*D0C + qK[3]*D0D + qK[4]*D0E + qK[5]*D0F + qK[6]*D0G + qK[7]*D0H + qK[8]*D0I + qK[9]*D0J + qK[11]*D0L + qK[12]*D0M + qK[13]*D0N + qK[14]*D0O + qK[15]*D0P + qK[16]*D0Q + qK[17]*D0R + qK[18]*D0S + qK[19]*D0T + qK[20]*D0U + qK[21]*D0V + qK[22]*D0W + qK[23]*D0X + qK[24]*D0Y + qK[25]*D0Z);
    
    ydot[37] =  - (lambda0L + psi + (qL[0] + qL[1] + qL[2] + qL[3] + qL[4] + qL[5] + qL[6] + qL[7] + qL[8] + qL[9] + qL[10] + qL[12] + qL[13] + qL[14] + qL[15] + qL[16] + qL[17] + qL[18] + qL[19] + qL[20] + qL[21] + qL[22] + qL[23] + qL[24] + qL[25]) + mu0L) * D0L + (qL[0]*D0A + qL[1]*D0B + qL[2]*D0C + qL[3]*D0D + qL[4]*D0E + qL[5]*D0F + qL[6]*D0G + qL[7]*D0H + qL[8]*D0I + qL[9]*D0J + qL[10]*D0K + qL[12]*D0M + qL[13]*D0N + qL[14]*D0O + qL[15]*D0P + qL[16]*D0Q + qL[17]*D0R + qL[18]*D0S + qL[19]*D0T + qL[20]*D0U + qL[21]*D0V + qL[22]*D0W + qL[23]*D0X + qL[24]*D0Y + qL[25]*D0Z);
    
    ydot[38] =  - (lambda0M + psi + (qM[0] + qM[1] + qM[2] + qM[3] + qM[4] + qM[5] + qM[6] + qM[7] + qM[8] + qM[9] + qM[10] + qM[11] + qM[13] + qM[14] + qM[15] + qM[16] + qM[17] + qM[18] + qM[19] + qM[20] + qM[21] + qM[22] + qM[23] + qM[24] + qM[25]) + mu0M) * D0M + (qM[0]*D0A + qM[1]*D0B + qM[2]*D0C + qM[3]*D0D + qM[4]*D0E + qM[5]*D0F + qM[6]*D0G + qM[7]*D0H + qM[8]*D0I + qM[9]*D0J + qM[10]*D0K + qM[11]*D0L + qM[13]*D0N + qM[14]*D0O + qM[15]*D0P + qM[16]*D0Q + qM[17]*D0R + qM[18]*D0S + qM[19]*D0T + qM[20]*D0U + qM[21]*D0V + qM[22]*D0W + qM[23]*D0X + qM[24]*D0Y + qM[25]*D0Z);
    
    ydot[39] =  - (lambda0N + psi + (qN[0] + qN[1] + qN[2] + qN[3] + qN[4] + qN[5] + qN[6] + qN[7] + qN[8] + qN[9] + qN[10] + qN[11] + qN[12] + qN[14] + qN[15] + qN[16] + qN[17] + qN[18] + qN[19] + qN[20] + qN[21] + qN[22] + qN[23] + qN[24] + qN[25]) + mu0N) * D0N + (qN[0]*D0A + qN[1]*D0B + qN[2]*D0C + qN[3]*D0D + qN[4]*D0E + qN[5]*D0F + qN[6]*D0G + qN[7]*D0H + qN[8]*D0I + qN[9]*D0J + qN[10]*D0K + qN[11]*D0L + qN[12]*D0M + qN[14]*D0O + qN[15]*D0P + qN[16]*D0Q + qN[17]*D0R + qN[18]*D0S + qN[19]*D0T + qN[20]*D0U + qN[21]*D0V + qN[22]*D0W + qN[23]*D0X + qN[24]*D0Y + qN[25]*D0Z);
    
    ydot[40] =  - (lambda0O + psi + (qO[0] + qO[1] + qO[2] + qO[3] + qO[4] + qO[5] + qO[6] + qO[7] + qO[8] + qO[9] + qO[10] + qO[11] + qO[12] + qO[13] + qO[15] + qO[16] + qO[17] + qO[18] + qO[19] + qO[20] + qO[21] + qO[22] + qO[23] + qO[24] + qO[25]) + mu0O) * D0O + (qO[0]*D0A + qO[1]*D0B + qO[2]*D0C + qO[3]*D0D + qO[4]*D0E + qO[5]*D0F + qO[6]*D0G + qO[7]*D0H + qO[8]*D0I + qO[9]*D0J + qO[10]*D0K + qO[11]*D0L + qO[12]*D0M + qO[13]*D0N + qO[15]*D0P + qO[16]*D0Q + qO[17]*D0R + qO[18]*D0S + qO[19]*D0T + qO[20]*D0U + qO[21]*D0V + qO[22]*D0W + qO[23]*D0X + qO[24]*D0Y + qO[25]*D0Z);
    
    ydot[41] =  - (lambda0P + psi + (qP[0] + qP[1] + qP[2] + qP[3] + qP[4] + qP[5] + qP[6] + qP[7] + qP[8] + qP[9] + qP[10] + qP[11] + qP[12] + qP[13] + qP[14] + qP[16] + qP[17] + qP[18] + qP[19] + qP[20] + qP[21] + qP[22] + qP[23] + qP[24] + qP[25]) + mu0P) * D0P + (qP[0]*D0A + qP[1]*D0B + qP[2]*D0C + qP[3]*D0D + qP[4]*D0E + qP[5]*D0F + qP[6]*D0G + qP[7]*D0H + qP[8]*D0I + qP[9]*D0J + qP[10]*D0K + qP[11]*D0L + qP[12]*D0M + qP[13]*D0N + qP[14]*D0O + qP[16]*D0Q + qP[17]*D0R + qP[18]*D0S + qP[19]*D0T + qP[20]*D0U + qP[21]*D0V + qP[22]*D0W + qP[23]*D0X + qP[24]*D0Y + qP[25]*D0Z);
    
    ydot[42] =  - (lambda0Q + psi + (qQ[0] + qQ[1] + qQ[2] + qQ[3] + qQ[4] + qQ[5] + qQ[6] + qQ[7] + qQ[8] + qQ[9] + qQ[10] + qQ[11] + qQ[12] + qQ[13] + qQ[14] + qQ[15] + qQ[17] + qQ[18] + qQ[19] + qQ[20] + qQ[21] + qQ[22] + qQ[23] + qQ[24] + qQ[25]) + mu0Q) * D0Q + (qQ[0]*D0A + qQ[1]*D0B + qQ[2]*D0C + qQ[3]*D0D + qQ[4]*D0E + qQ[5]*D0F + qQ[6]*D0G + qQ[7]*D0H + qQ[8]*D0I + qQ[9]*D0J + qQ[10]*D0K + qQ[11]*D0L + qQ[12]*D0M + qQ[13]*D0N + qQ[14]*D0O + qQ[15]*D0P + qQ[17]*D0R + qQ[18]*D0S + qQ[19]*D0T + qQ[20]*D0U + qQ[21]*D0V + qQ[22]*D0W + qQ[23]*D0X + qQ[24]*D0Y + qQ[25]*D0Z);
    
    ydot[43] =  - (lambda0R + psi + (qR[0] + qR[1] + qR[2] + qR[3] + qR[4] + qR[5] + qR[6] + qR[7] + qR[8] + qR[9] + qR[10] + qR[11] + qR[12] + qR[13] + qR[14] + qR[15] + qR[16] + qR[18] + qR[19] + qR[20] + qR[21] + qR[22] + qR[23] + qR[24] + qR[25]) + mu0R) * D0R + (qR[0]*D0A + qR[1]*D0B + qR[2]*D0C + qR[3]*D0D + qR[4]*D0E + qR[5]*D0F + qR[6]*D0G + qR[7]*D0H + qR[8]*D0I + qR[9]*D0J + qR[10]*D0K + qR[11]*D0L + qR[12]*D0M + qR[13]*D0N + qR[14]*D0O + qR[15]*D0P + qR[16]*D0Q + qR[18]*D0S + qR[19]*D0T + qR[20]*D0U + qR[21]*D0V + qR[22]*D0W + qR[23]*D0X + qR[24]*D0Y + qR[25]*D0Z);
    
    ydot[44] =  - (lambda0S + psi + (qS[0] + qS[1] + qS[2] + qS[3] + qS[4] + qS[5] + qS[6] + qS[7] + qS[8] + qS[9] + qS[10] + qS[11] + qS[12] + qS[13] + qS[14] + qS[15] + qS[16] + qS[17] + qS[19] + qS[20] + qS[21] + qS[22] + qS[23] + qS[24] + qS[25]) + mu0S) * D0S + (qS[0]*D0A + qS[1]*D0B + qS[2]*D0C + qS[3]*D0D + qS[4]*D0E + qS[5]*D0F + qS[6]*D0G + qS[7]*D0H + qS[8]*D0I + qS[9]*D0J + qS[10]*D0K + qS[11]*D0L + qS[12]*D0M + qS[13]*D0N + qS[14]*D0O + qS[15]*D0P + qS[16]*D0Q + qS[17]*D0R + qS[19]*D0T + qS[20]*D0U + qS[21]*D0V + qS[22]*D0W + qS[23]*D0X + qS[24]*D0Y + qS[25]*D0Z);
    
    ydot[45] =  - (lambda0T + psi + (qT[0] + qT[1] + qT[2] + qT[3] + qT[4] + qT[5] + qT[6] + qT[7] + qT[8] + qT[9] + qT[10] + qT[11] + qT[12] + qT[13] + qT[14] + qT[15] + qT[16] + qT[17] + qT[18] + qT[20] + qT[21] + qT[22] + qT[23] + qT[24] + qT[25]) + mu0T) * D0T + (qT[0]*D0A + qT[1]*D0B + qT[2]*D0C + qT[3]*D0D + qT[4]*D0E + qT[5]*D0F + qT[6]*D0G + qT[7]*D0H + qT[8]*D0I + qT[9]*D0J + qT[10]*D0K + qT[11]*D0L + qT[12]*D0M + qT[13]*D0N + qT[14]*D0O + qT[15]*D0P + qT[16]*D0Q + qT[17]*D0R + qT[18]*D0S + qT[20]*D0U + qT[21]*D0V + qT[22]*D0W + qT[23]*D0X + qT[24]*D0Y + qT[25]*D0Z);
    
    ydot[46] =  - (lambda0U + psi + (qU[0] + qU[1] + qU[2] + qU[3] + qU[4] + qU[5] + qU[6] + qU[7] + qU[8] + qU[9] + qU[10] + qU[11] + qU[12] + qU[13] + qU[14] + qU[15] + qU[16] + qU[17] + qU[18] + qU[19] + qU[21] + qU[22] + qU[23] + qU[24] + qU[25]) + mu0U) * D0U + (qU[0]*D0A + qU[1]*D0B + qU[2]*D0C + qU[3]*D0D + qU[4]*D0E + qU[5]*D0F + qU[6]*D0G + qU[7]*D0H + qU[8]*D0I + qU[9]*D0J + qU[10]*D0K + qU[11]*D0L + qU[12]*D0M + qU[13]*D0N + qU[14]*D0O + qU[15]*D0P + qU[16]*D0Q + qU[17]*D0R + qU[18]*D0S + qU[19]*D0T + qU[21]*D0V + qU[22]*D0W + qU[23]*D0X + qU[24]*D0Y + qU[25]*D0Z);
    
    ydot[47] =  - (lambda0V + psi + (qV[0] + qV[1] + qV[2] + qV[3] + qV[4] + qV[5] + qV[6] + qV[7] + qV[8] + qV[9] + qV[10] + qV[11] + qV[12] + qV[13] + qV[14] + qV[15] + qV[16] + qV[17] + qV[18] + qV[19] + qV[20] + qV[22] + qV[23] + qV[24] + qV[25]) + mu0V) * D0V + (qV[0]*D0A + qV[1]*D0B + qV[2]*D0C + qV[3]*D0D + qV[4]*D0E + qV[5]*D0F + qV[6]*D0G + qV[7]*D0H + qV[8]*D0I + qV[9]*D0J + qV[10]*D0K + qV[11]*D0L + qV[12]*D0M + qV[13]*D0N + qV[14]*D0O + qV[15]*D0P + qV[16]*D0Q + qV[17]*D0R + qV[18]*D0S + qV[19]*D0T + qV[20]*D0U + qV[22]*D0W + qV[23]*D0X + qV[24]*D0Y + qV[25]*D0Z);
    
    ydot[48] =  - (lambda0W + psi + (qW[0] + qW[1] + qW[2] + qW[3] + qW[4] + qW[5] + qW[6] + qW[7] + qW[8] + qW[9] + qW[10] + qW[11] + qW[12] + qW[13] + qW[14] + qW[15] + qW[16] + qW[17] + qW[18] + qW[19] + qW[20] + qW[21] + qW[23] + qW[24] + qW[25]) + mu0W) * D0W + (qW[0]*D0A + qW[1]*D0B + qW[2]*D0C + qW[3]*D0D + qW[4]*D0E + qW[5]*D0F + qW[6]*D0G + qW[7]*D0H + qW[8]*D0I + qW[9]*D0J + qW[10]*D0K + qW[11]*D0L + qW[12]*D0M + qW[13]*D0N + qW[14]*D0O + qW[15]*D0P + qW[16]*D0Q + qW[17]*D0R + qW[18]*D0S + qW[19]*D0T + qW[20]*D0U + qW[21]*D0V + qW[23]*D0X + qW[24]*D0Y + qW[25]*D0Z);
    
    ydot[49] =  - (lambda0X + psi + (qX[0] + qX[1] + qX[2] + qX[3] + qX[4] + qX[5] + qX[6] + qX[7] + qX[8] + qX[9] + qX[10] + qX[11] + qX[12] + qX[13] + qX[14] + qX[15] + qX[16] + qX[17] + qX[18] + qX[19] + qX[20] + qX[21] + qX[22] + qX[24] + qX[25]) + mu0X) * D0X + (qX[0]*D0A + qX[1]*D0B + qX[2]*D0C + qX[3]*D0D + qX[4]*D0E + qX[5]*D0F + qX[6]*D0G + qX[7]*D0H + qX[8]*D0I + qX[9]*D0J + qX[10]*D0K + qX[11]*D0L + qX[12]*D0M + qX[13]*D0N + qX[14]*D0O + qX[15]*D0P + qX[16]*D0Q + qX[17]*D0R + qX[18]*D0S + qX[19]*D0T + qX[20]*D0U + qX[21]*D0V + qX[22]*D0W + qX[24]*D0Y + qX[25]*D0Z);
    
    ydot[50] =  - (lambda0Y + psi + (qY[0] + qY[1] + qY[2] + qY[3] + qY[4] + qY[5] + qY[6] + qY[7] + qY[8] + qY[9] + qY[10] + qY[11] + qY[12] + qY[13] + qY[14] + qY[15] + qY[16] + qY[17] + qY[18] + qY[19] + qY[20] + qY[21] + qY[22] + qY[23] + qY[25]) + mu0Y) * D0Y + (qY[0]*D0A + qY[1]*D0B + qY[2]*D0C + qY[3]*D0D + qY[4]*D0E + qY[5]*D0F + qY[6]*D0G + qY[7]*D0H + qY[8]*D0I + qY[9]*D0J + qY[10]*D0K + qY[11]*D0L + qY[12]*D0M + qY[13]*D0N + qY[14]*D0O + qY[15]*D0P + qY[16]*D0Q + qY[17]*D0R + qY[18]*D0S + qY[19]*D0T + qY[20]*D0U + qY[21]*D0V + qY[22]*D0W + qY[23]*D0X + qY[25]*D0Z);
    
    ydot[51] =  - (lambda0Z + psi + (qZ[0] + qZ[1] + qZ[2] + qZ[3] + qZ[4] + qZ[5] + qZ[6] + qZ[7] + qZ[8] + qZ[9] + qZ[10] + qZ[11] + qZ[12] + qZ[13] + qZ[14] + qZ[15] + qZ[16] + qZ[17] + qZ[18] + qZ[19] + qZ[20] + qZ[21] + qZ[22] + qZ[23] + qZ[24]) + mu0Z) * D0Z + (qZ[0]*D0A + qZ[1]*D0B + qZ[2]*D0C + qZ[3]*D0D + qZ[4]*D0E + qZ[5]*D0F + qZ[6]*D0G + qZ[7]*D0H + qZ[8]*D0I + qZ[9]*D0J + qZ[10]*D0K + qZ[11]*D0L + qZ[12]*D0M + qZ[13]*D0N + qZ[14]*D0O + qZ[15]*D0P + qZ[16]*D0Q + qZ[17]*D0R + qZ[18]*D0S + qZ[19]*D0T + qZ[20]*D0U + qZ[21]*D0V + qZ[22]*D0W + qZ[23]*D0X + qZ[24]*D0Y);
    
}


