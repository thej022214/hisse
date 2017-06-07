/*
 *  higeosse-ext-derivs_c
 *
 *
 *  Created by Jeremy Beaulieu 6/06/2017
 *  Copyright 2017 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 95

static double params_higeosse[NUMELEMENTS];


void initmod_higeosse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_higeosse);
}



void higeosse_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E_0A = y[0];
    double E_1A = y[1];
    double E_01A = y[2];
    double E_0B = y[3];
    double E_1B = y[4];
    double E_01B = y[5];
    double E_0C = y[6];
    double E_1C = y[7];
    double E_01C = y[8];
    double E_0D = y[9];
    double E_1D = y[10];
    double E_01D = y[11];
    double E_0E = y[12];
    double E_1E = y[13];
    double E_01E = y[14];
    
    double D_N0A = y[15];
    double D_N1A = y[16];
    double D_N01A = y[17];
    double D_N0B = y[18];
    double D_N1B = y[19];
    double D_N01B = y[20];
    double D_N0C = y[21];
    double D_N1C = y[22];
    double D_N01C = y[23];
    double D_N0D = y[24];
    double D_N1D = y[25];
    double D_N01D = y[26];
    double D_N0E = y[27];
    double D_N1E = y[28];
    double D_N01E = y[29];
    
    double
    s0A = params_higeosse[0],
    s1A = params_higeosse[1],
    s01A = params_higeosse[2],
    x0A = params_higeosse[3],
    x1A = params_higeosse[4],
    d0A_01A = params_higeosse[5],
    d1A_01A = params_higeosse[6],
    d01A_0A = x1A,
    d01A_1A = x0A,
    
    d0A_0B = params_higeosse[7],
    d0A_0C = params_higeosse[8],
    d0A_0D = params_higeosse[9],
    d0A_0E = params_higeosse[10],
    d1A_1B = params_higeosse[11],
    d1A_1C = params_higeosse[12],
    d1A_1D = params_higeosse[13],
    d1A_1E = params_higeosse[14],
    d01A_01B = params_higeosse[15],
    d01A_01C = params_higeosse[16],
    d01A_01D = params_higeosse[17],
    d01A_01E = params_higeosse[18],
    
    s0B = params_higeosse[19],
    s1B = params_higeosse[20],
    s01B = params_higeosse[21],
    x0B = params_higeosse[22],
    x1B = params_higeosse[23],
    d0B_01B = params_higeosse[24],
    d1B_01B = params_higeosse[25],
    d01B_0B = x1B,
    d01B_1B = x0B,

    d0B_0A = params_higeosse[26],
    d0B_0C = params_higeosse[27],
    d0B_0D = params_higeosse[28],
    d0B_0E = params_higeosse[29],
    d1B_1A = params_higeosse[30],
    d1B_1C = params_higeosse[31],
    d1B_1D = params_higeosse[32],
    d1B_1E = params_higeosse[33],
    d01B_01A = params_higeosse[34],
    d01B_01C = params_higeosse[35],
    d01B_01D = params_higeosse[36],
    d01B_01E = params_higeosse[37],
    
    s0C = params_higeosse[38],
    s1C = params_higeosse[39],
    s01C = params_higeosse[40],
    x0C = params_higeosse[41],
    x1C = params_higeosse[42],
    d0C_01C = params_higeosse[43],
    d1C_01C = params_higeosse[44],
    d01C_0C = x1C,
    d01C_1C = x0C,

    d0C_0A = params_higeosse[45],
    d0C_0B = params_higeosse[46],
    d0C_0D = params_higeosse[47],
    d0C_0E = params_higeosse[48],
    d1C_1A = params_higeosse[49],
    d1C_1B = params_higeosse[50],
    d1C_1D = params_higeosse[51],
    d1C_1E = params_higeosse[52],
    d01C_01A = params_higeosse[53],
    d01C_01B = params_higeosse[54],
    d01C_01D = params_higeosse[55],
    d01C_01E = params_higeosse[56],
    
    s0D = params_higeosse[57],
    s1D = params_higeosse[58],
    s01D = params_higeosse[59],
    x0D = params_higeosse[60],
    x1D = params_higeosse[61],
    d0D_01D = params_higeosse[62],
    d1D_01D = params_higeosse[63],
    d01D_0D = x1D,
    d01D_1D = x0D,

    d0D_0A = params_higeosse[64],
    d0D_0B = params_higeosse[65],
    d0D_0C = params_higeosse[66],
    d0D_0E = params_higeosse[67],
    d1D_1A = params_higeosse[68],
    d1D_1B = params_higeosse[69],
    d1D_1C = params_higeosse[70],
    d1D_1E = params_higeosse[71],
    d01D_01A = params_higeosse[72],
    d01D_01B = params_higeosse[73],
    d01D_01C = params_higeosse[74],
    d01D_01E = params_higeosse[75],
    
    s0E = params_higeosse[76],
    s1E = params_higeosse[77],
    s01E = params_higeosse[78],
    x0E = params_higeosse[79],
    x1E = params_higeosse[80],
    d0E_01E = params_higeosse[81],
    d1E_01E = params_higeosse[82],
    d01E_0E = x1E,
    d01E_1E = x0E,

    d0E_0A = params_higeosse[83],
    d0E_0B = params_higeosse[84],
    d0E_0C = params_higeosse[85],
    d0E_0D = params_higeosse[86],
    d1E_1A = params_higeosse[87],
    d1E_1B = params_higeosse[88],
    d1E_1C = params_higeosse[89],
    d1E_1D = params_higeosse[90],
    d01E_01A = params_higeosse[91],
    d01E_01B = params_higeosse[92],
    d01E_01C = params_higeosse[93],
    d01E_01D = params_higeosse[94];

    
    /* The E's */
    /*  dE_0A / dt  */
    ydot[0]  = -(s0A + d0A_01A + d0A_0B + d0A_0C + d0A_0D + d0A_0E + x0A) * E_0A + (d0A_01A * E_01A + d0A_0B * E_0B + d0A_0C * E_0C + d0A_0D * E_0D + d0A_0E * E_0E) + x0A + (s0A * E_0A * E_0A);
    
    /*  E_1A / dt  */
    ydot[1]  = -(s1A + d1A_01A + d1A_1B + d1A_1C + d1A_1D + d1A_1E + x1A) * E_1A + (d1A_01A * E_01A + d1A_1B * E_1B + d1A_1C * E_1C + d1A_1D * E_1D + d1A_1E * E_1E) + x1A + (s1A * E_1A * E_1A);
    
    /*  dE_01A / dt  */
    ydot[2] = -(s01A + s0A + s1A + d01A_0A + d01A_1A + d01A_01B + d01A_01C + d01A_01D + d01A_01E) * E_01A + (d01A_0A * E_0A + d01A_1A * E_1A + d01A_01B * E_01B + d01A_01C * E_01C + d01A_01D * E_01D + d01A_01E * E_01E) + s0A * E_0A * E_01A + s1A * E_01A * E_1A + s01A * E_0A * E_1A;
    
    
    /*  dE_0B / dt  */
    ydot[3]  = -(s0B + d0B_01B + d0B_0A + d0B_0C + d0B_0D + d0B_0E + x0B) * E_0B + (d0B_01B * E_01B + d0B_0A * E_0A + d0B_0C * E_0C + d0B_0D * E_0D + d0B_0E * E_0E) + x0B + (s0B * E_0B * E_0B);
    
    /*  E_1B / dt  */
    ydot[4]  = -(s1B + d1B_01B + d1B_1A + d1B_1C + d1B_1D + d1B_1E + x1B) * E_1B + (d1B_01B * E_01B + d1B_1A * E_1A + d1B_1C * E_1C + d1B_1D * E_1D + d1B_1E * E_1E) + x1B + (s1B * E_1B * E_1B);
    
    /*  dE_01B / dt  */
    ydot[5] = -(s01B + s0B + s1B + d01B_0B + d01B_1B + d01B_01A + d01B_01C + d01B_01D + d01B_01E) * E_01B + (d01B_0B * E_0B + d01B_1B * E_1B + d01B_01A * E_01A + d01B_01C * E_01C + d01B_01D * E_01D + d01B_01E * E_01E) + s0B * E_0B * E_01B + s1B * E_01B * E_1B + s01B * E_0B * E_1B;


    /*  dE_0C / dt  */
    ydot[6]  = -(s0C + d0C_01C + d0C_0B + d0C_0A + d0C_0D + d0C_0E + x0C) * E_0C + (d0C_01C * E_01C + d0C_0B * E_0B + d0C_0A * E_0A + d0C_0D * E_0D + d0C_0E * E_0E) + x0C + (s0C * E_0C * E_0C);
    
    /*  E_1C / dt  */
    ydot[7]  = -(s1C + d1C_01C + d1C_1B + d1C_1A + d1C_1D + d1C_1E + x1C) * E_1C + (d1C_01C * E_01C + d1C_1B * E_1B + d1C_1A * E_1A + d1C_1D * E_1D + d1C_1E * E_1E) + x1C + (s1C * E_1C * E_1C);
    
    /*  dE_01C/ dt  */
    ydot[8] = -(s01C + s0C + s1C + d01C_0C + d01C_1C + d01C_01B + d01C_01A + d01C_01D + d01C_01E) * E_01C + (d01C_0C * E_0C + d01C_1C * E_1C + d01C_01B * E_01B + d01C_01A * E_01A + d01C_01D * E_01D + d01C_01E * E_01E) + s0C * E_0C * E_01C + s1C * E_01C * E_1C + s01C * E_0C * E_1C;

    
    /*  dE_0D / dt  */
    ydot[9]  = -(s0D + d0D_01D + d0D_0B + d0D_0A + d0D_0C + d0D_0E + x0D) * E_0D + (d0D_01D * E_01D + d0D_0B * E_0B + d0D_0A * E_0A + d0D_0C * E_0C + d0D_0E * E_0E) + x0D + (s0D * E_0D * E_0D);
    
    /*  E_1D / dt  */
    ydot[10] = -(s1D + d1D_01D + d1D_1B + d1D_1A + d1D_1C + d1D_1E + x1D) * E_1D + (d1D_01D * E_01D + d1D_1B * E_1B + d1D_1A * E_1A + d1D_1C * E_1C + d1D_1E * E_1E) + x1D + (s1D * E_1D * E_1D);
    
    /*  dE_01D/ dt  */
    ydot[11] = -(s01D + s0D + s1D + d01D_0D + d01D_1D + d01D_01B + d01D_01A + d01D_01C + d01D_01E) * E_01D + (d01D_0D * E_0D + d01D_1D * E_1D + d01D_01B * E_01B + d01D_01A * E_01A + d01D_01C * E_01C + d01D_01E * E_01E) + s0D * E_0D * E_01D + s1D * E_01D * E_1D + s01D * E_0D * E_1D;
    

    /*  dE_0E / dt  */
    ydot[12] = -(s0E + d0E_01E + d0E_0B + d0E_0A + d0E_0C + d0E_0D + x0E) * E_0E + (d0E_01E * E_01E + d0E_0B * E_0B + d0E_0A * E_0A + d0E_0C * E_0C + d0E_0D * E_0D) + x0E + (s0E * E_0E * E_0E);
    
    /*  E_1E / dt  */
    ydot[13] = -(s1E + d1E_01E + d1E_1B + d1E_1A + d1E_1C + d1E_1D + x1E) * E_1E + (d1E_01E * E_01E + d1E_1B * E_1B + d1E_1A * E_1A + d1E_1C * E_1C + d1E_1D * E_1D) + x1E + (s1E * E_1E * E_1E);
    
    /*  dE_01E/ dt  */
    ydot[14] = -(s01E + s0E + s1E + d01E_0E + d01E_1E + d01E_01B + d01E_01A + d01E_01C + d01E_01D) * E_01E + (d01E_0E * E_0E + d01E_1E * E_1E + d01E_01B * E_01B + d01E_01A * E_01A + d01E_01C * E_01C + d01E_01D * E_01D) + s0E * E_0E * E_01E + s1E * E_01E * E_1E + s01E * E_0E * E_1E;

    
    
    
    /* The D's */
    
    /*  dD_N0A / dt  */
    ydot[15] = -(s0A + d0A_01A + d0A_0B + d0A_0C + d0A_0D + d0A_0E + x0A) * D_N0A + (d0A_01A * D_N01A + d0A_0B * D_N0B + d0A_0C * D_N0C + d0A_0D * D_N0D + d0A_0E * D_N0E) + s0A * (D_N0A * E_0A + D_N0A * E_0A);
    
    /*  dD_N1A / dt  */
    ydot[16] = -(s1A + d1A_01A + d1A_1B + d1A_1C + d1A_1D + d1A_1E + x1A) * D_N1A + (d1A_01A * D_N01A + d1A_1B * D_N1B + d1A_1C * D_N1C + d1A_1D * D_N1D + d1A_1E * D_N1E) + s1A * (D_N1A * E_1A + D_N1A * E_1A);
    
    /*  dD_N01A / dt  */
    ydot[17] = -(s01A + s0A + s1A + d01A_0A + d01A_1A + d01A_01B + d01A_01C + d01A_01D + d01A_01E) * D_N01A + (d01A_0A * D_N0A + d01A_1A * D_N1A + d01A_01B * D_N01B + d01A_01C * D_N01C + d01A_01D * D_N01D + d01A_01E * D_N01E) + s01A * (D_N0A * E_1A + D_N1A * E_0A) + s0A * (E_0A * D_N01A + E_01A * D_N0A) + s1A * (E_1A * D_N01A + E_01A * D_N1A);

    
    /*  dD_N0B / dt  */
    ydot[18] = -(s0B + d0B_01B + d0B_0A + d0B_0C + d0B_0D + d0B_0E + x0B) * D_N0B + (d0B_01B * D_N01B + d0B_0A * D_N0A + d0B_0C * D_N0C + d0B_0D * D_N0D + d0B_0E * D_N0E) + s0B * (D_N0B * E_0B + D_N0B * E_0B);
    
    /*  dD_N1B / dt  */
    ydot[19] = -(s1B + d1B_01B + d1B_1A + d1B_1C + d1B_1D + d1B_1E + x1B) * D_N1B + (d1B_01B * D_N01B + d1B_1A * D_N1A + d1B_1C * D_N1C + d1B_1D * D_N1D + d1B_1E * D_N1E) + s1B * (D_N1B * E_1B + D_N1B * E_1B);
    
    /*  dD_N01B / dt  */
    ydot[20] = -(s01B + s0B + s1B + d01B_0B + d01B_1B + d01B_01A + d01B_01C + d01B_01D + d01B_01E) * D_N01B + (d01B_0B * D_N0B + d01B_1B * D_N1B + d01B_01A * D_N01A + d01B_01C * D_N01C + d01B_01D * D_N01D + d01B_01E * D_N01E) + s01B * (D_N0B * E_1B + D_N1B * E_0B) + s0B * (E_0B * D_N01B + E_01B * D_N0B) + s1B * (E_1B * D_N01B + E_01B * D_N1B);

    
    /*  dD_N0C / dt  */
    ydot[21] = -(s0C + d0C_01C + d0C_0B + d0C_0A + d0C_0D + d0C_0E + x0C) * D_N0C + (d0C_01C * D_N01C + d0C_0B * D_N0B + d0C_0A * D_N0A + d0C_0D * D_N0D + d0C_0E * D_N0E) + s0C * (D_N0C * E_0C + D_N0C * E_0C);
    
    /*  dD_N1C / dt  */
    ydot[22] = -(s1C + d1C_01C + d1C_1B + d1C_1A + d1C_1D + d1C_1E + x1C) * D_N1C + (d1C_01C * D_N01C + d1C_1B * D_N1B + d1C_1A * D_N1A + d1C_1D * D_N1D + d1C_1E * D_N1E) + s1C * (D_N1C * E_1C + D_N1C * E_1C);
    
    /*  dD_N01C / dt  */
    ydot[23] = -(s01C + s0C + s1C + d01C_0C + d01C_1C + d01C_01A + d01C_01B + d01C_01D + d01C_01E) * D_N01C + (d01C_0C * D_N0C + d01C_1C * D_N1C + d01C_01B * D_N01B + d01C_01A * D_N01A + d01C_01D * D_N01D + d01C_01E * D_N01E) + s01C * (D_N0C * E_1C + D_N1C * E_0C) + s0C * (E_0C * D_N01C + E_01C * D_N0C) + s1C * (E_1C * D_N01C + E_01C * D_N1C);

    
    /*  dD_N0D / dt  */
    ydot[24] = -(s0D + d0D_01D + d0D_0B + d0D_0A + d0D_0C + d0D_0E + x0D) * D_N0D + (d0D_01D * D_N01D + d0D_0B * D_N0B + d0D_0A * D_N0A + d0D_0C * D_N0C + d0D_0E * D_N0E) + s0D * (D_N0D * E_0D + D_N0D * E_0D);
    
    /*  dD_N1D / dt  */
    ydot[25] = -(s1D + d1D_01D + d1D_1B + d1D_1A + d1D_1C + d1D_1E + x1D) * D_N1D + (d1D_01D * D_N01D + d1D_1B * D_N1B + d1D_1A * D_N1A + d1D_1C * D_N1C + d1D_1E * D_N1E) + s1D * (D_N1D * E_1D + D_N1D * E_1D);
    
    /*  dD_N01D / dt  */
    ydot[26] = -(s01D + s0D + s1D + d01D_0D + d01D_1D + d01D_01A + d01D_01B + d01D_01C + d01D_01E) * D_N01D + (d01D_0D * D_N0D + d01D_1D * D_N1D + d01D_01B * D_N01B + d01D_01A * D_N01A + d01D_01C * D_N01C + d01D_01E * D_N01E) + s01D * (D_N0D * E_1D + D_N1D * E_0D) + s0D * (E_0D * D_N01D + E_01D * D_N0D) + s1D * (E_1D * D_N01D + E_01D * D_N1D);

    
    /*  dD_N0E / dt  */
    ydot[27] = -(s0E + d0E_01E + d0E_0B + d0E_0A + d0E_0C + d0E_0D + x0E) * D_N0E + (d0E_01E * D_N01E + d0E_0B * D_N0B + d0E_0A * D_N0A + d0E_0C * D_N0C + d0E_0D * D_N0D) + s0E * (D_N0E * E_0E + D_N0E * E_0E);
    
    /*  dD_N1E / dt  */
    ydot[28] = -(s1E + d1E_01E + d1E_1B + d1E_1A + d1E_1C + d1E_1D + x1E) * D_N1E + (d1E_01E * D_N01E + d1E_1B * D_N1B + d1E_1A * D_N1A + d1E_1C * D_N1C + d1E_1D * D_N1D) + s1E * (D_N1E * E_1E + D_N1E * E_1E);
    
    /*  dD_N01E / dt  */
    ydot[29] = -(s01E + s0E + s1E + d01E_0E + d01E_1E + d01E_01A + d01E_01B + d01E_01C + d01E_01D) * D_N01E + (d01E_0E * D_N0E + d01E_1E * D_N1E + d01E_01B * D_N01B + d01E_01A * D_N01A + d01E_01C * D_N01C + d01E_01D * D_N01D) + s01E * (D_N0E * E_1E + D_N1E * E_0E) + s0E * (E_0E * D_N01E + E_01E * D_N0E) + s1E * (E_1E * D_N01E + E_01E * D_N1E);

}













