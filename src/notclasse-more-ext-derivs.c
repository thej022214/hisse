/*
 *  notclasse-more-ext-derivs_c
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
#define NUMELEMENTS 120

static double params_hinoclass[NUMELEMENTS];


void initmod_hinoclass(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_hinoclass);
}



void notclasse_more_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
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
    s0A = params_hinoclass[0],
    s1A = params_hinoclass[1],
    s01A = params_hinoclass[2],
    x0A = params_hinoclass[3],
    x1A = params_hinoclass[4],
    x01A = params_hinoclass[5],
    d0A_1A  = params_hinoclass[6],
    d0A_01A = params_hinoclass[7],
    d1A_0A = params_hinoclass[8],
    d1A_01A = params_hinoclass[9],
    d01A_0A = params_hinoclass[10],
    d01A_1A = params_hinoclass[11],
    
    d0A_0B = params_hinoclass[12],
    d0A_0C = params_hinoclass[13],
    d0A_0D = params_hinoclass[14],
    d0A_0E = params_hinoclass[15],
    d1A_1B = params_hinoclass[16],
    d1A_1C = params_hinoclass[17],
    d1A_1D = params_hinoclass[18],
    d1A_1E = params_hinoclass[19],
    d01A_01B = params_hinoclass[20],
    d01A_01C = params_hinoclass[21],
    d01A_01D = params_hinoclass[22],
    d01A_01E = params_hinoclass[23],
    
    s0B = params_hinoclass[24],
    s1B = params_hinoclass[25],
    s01B = params_hinoclass[26],
    x0B = params_hinoclass[27],
    x1B = params_hinoclass[28],
    x01B = params_hinoclass[29],
    d0B_1B  = params_hinoclass[30],
    d0B_01B = params_hinoclass[31],
    d1B_0B = params_hinoclass[32],
    d1B_01B = params_hinoclass[33],
    d01B_0B = params_hinoclass[34],
    d01B_1B = params_hinoclass[35],
    d0B_0A = params_hinoclass[36],
    d0B_0C = params_hinoclass[37],
    d0B_0D = params_hinoclass[38],
    d0B_0E = params_hinoclass[39],
    d1B_1A = params_hinoclass[40],
    d1B_1C = params_hinoclass[41],
    d1B_1D = params_hinoclass[42],
    d1B_1E = params_hinoclass[43],
    d01B_01A = params_hinoclass[44],
    d01B_01C = params_hinoclass[45],
    d01B_01D = params_hinoclass[46],
    d01B_01E = params_hinoclass[47],
    
    s0C = params_hinoclass[48],
    s1C = params_hinoclass[49],
    s01C = params_hinoclass[50],
    x0C = params_hinoclass[51],
    x1C = params_hinoclass[52],
    x01C = params_hinoclass[53],
    d0C_1C  = params_hinoclass[54],
    d0C_01C = params_hinoclass[55],
    d1C_0C = params_hinoclass[56],
    d1C_01C = params_hinoclass[57],
    d01C_0C = params_hinoclass[58],
    d01C_1C = params_hinoclass[59],
    d0C_0A = params_hinoclass[60],
    d0C_0B = params_hinoclass[61],
    d0C_0D = params_hinoclass[62],
    d0C_0E = params_hinoclass[63],
    d1C_1A = params_hinoclass[64],
    d1C_1B = params_hinoclass[65],
    d1C_1D = params_hinoclass[66],
    d1C_1E = params_hinoclass[67],
    d01C_01A = params_hinoclass[68],
    d01C_01B = params_hinoclass[69],
    d01C_01D = params_hinoclass[70],
    d01C_01E = params_hinoclass[71],
    
    s0D = params_hinoclass[72],
    s1D = params_hinoclass[73],
    s01D = params_hinoclass[74],
    x0D = params_hinoclass[75],
    x1D = params_hinoclass[76],
    x01D = params_hinoclass[77],
    d0D_1D  = params_hinoclass[78],
    d0D_01D = params_hinoclass[79],
    d1D_0D = params_hinoclass[80],
    d1D_01D = params_hinoclass[81],
    d01D_0D = params_hinoclass[82],
    d01D_1D = params_hinoclass[83],
    d0D_0A = params_hinoclass[84],
    d0D_0B = params_hinoclass[85],
    d0D_0C = params_hinoclass[86],
    d0D_0E = params_hinoclass[87],
    d1D_1A = params_hinoclass[88],
    d1D_1B = params_hinoclass[89],
    d1D_1C = params_hinoclass[90],
    d1D_1E = params_hinoclass[91],
    d01D_01A = params_hinoclass[92],
    d01D_01B = params_hinoclass[93],
    d01D_01C = params_hinoclass[94],
    d01D_01E = params_hinoclass[95],
    
    s0E = params_hinoclass[96],
    s1E = params_hinoclass[97],
    s01E = params_hinoclass[98],
    x0E = params_hinoclass[99],
    x1E = params_hinoclass[100],
    x01E = params_hinoclass[101],
    d0E_1E  = params_hinoclass[102],
    d0E_01E = params_hinoclass[103],
    d1E_0E = params_hinoclass[104],
    d1E_01E = params_hinoclass[105],
    d01E_0E = params_hinoclass[106],
    d01E_1E = params_hinoclass[107],
    d0E_0A = params_hinoclass[108],
    d0E_0B = params_hinoclass[109],
    d0E_0C = params_hinoclass[110],
    d0E_0D = params_hinoclass[111],
    d1E_1A = params_hinoclass[112],
    d1E_1B = params_hinoclass[113],
    d1E_1C = params_hinoclass[114],
    d1E_1D = params_hinoclass[115],
    d01E_01A = params_hinoclass[116],
    d01E_01B = params_hinoclass[117],
    d01E_01C = params_hinoclass[118],
    d01E_01D = params_hinoclass[119];

    
    /* The E's */
    /*  dE_0A / dt  */
    ydot[0]  = -(s0A + d0A_1A + d0A_01A + d0A_0B + d0A_0C + d0A_0D + d0A_0E + x0A) * E_0A + (d0A_1A * E_1A + d0A_01A * E_01A + d0A_0B * E_0B + d0A_0C * E_0C + d0A_0D * E_0D + d0A_0E * E_0E) + x0A + (s0A * E_0A * E_0A);
    
    /*  dE_1A / dt  */
    ydot[1]  = -(s1A + d1A_0A + d1A_01A + d1A_1B + d1A_1C + d1A_1D + d1A_1E + x1A) * E_1A + (d1A_0A * E_0A + d1A_01A * E_01A + d1A_1B * E_1B + d1A_1C * E_1C + d1A_1D * E_1D + d1A_1E * E_1E) + x1A + (s1A * E_1A * E_1A);
    
    /*  dE_01A / dt  */
    ydot[2] = -(s01A + d01A_0A + d01A_1A + d01A_01B + d01A_01C + d01A_01D + d01A_01E + x01A) * E_01A + (d01A_0A * E_0A + d01A_1A * E_1A + d01A_01B * E_01B + d01A_01C * E_01C + d01A_01D * E_01D + d01A_01E * E_01E) + x01A + (s01A * E_01A * E_01A);
    
    
    /*  dE_0B / dt  */
    ydot[3]  = -(s0B + d0B_1B + d0B_01B + d0B_0A + d0B_0C + d0B_0D + d0B_0E + x0B) * E_0B + (d0B_1B * E_1B + d0B_01B * E_01B + d0B_0A * E_0A + d0B_0C * E_0C + d0B_0D * E_0D + d0B_0E * E_0E) + x0B + (s0B * E_0B * E_0B);
    
    /*  dE_1B / dt  */
    ydot[4]  = -(s1B + d1B_0B + d1B_01B + d1B_1A + d1B_1C + d1B_1D + d1B_1E + x1B) * E_1B + (d1B_0B * E_0B + d1B_01B * E_01B + d1B_1A * E_1A + d1B_1C * E_1C + d1B_1D * E_1D + d1B_1E * E_1E) + x1B + (s1B * E_1B * E_1B);
    
    /*  dE_01B / dt  */
    ydot[5] = -(s01B + d01B_0B + d01B_1B + d01B_01A + d01B_01C + d01B_01D + d01B_01E + x01B) * E_01B + (d01B_0B * E_0B + d01B_1B * E_1B + d01B_01A * E_01A + d01B_01C * E_01C + d01B_01D * E_01D + d01B_01E * E_01E) + x01B + (s01B * E_01B * E_01B);


    /*  dE_0C / dt  */
    ydot[6]  = -(s0C + d0C_1C + d0C_01C + d0C_0B + d0C_0A + d0C_0D + d0C_0E + x0C) * E_0C + (d0C_1C * E_1C + d0C_01C * E_01C + d0C_0B * E_0B + d0C_0A * E_0A + d0C_0D * E_0D + d0C_0E * E_0E) + x0C + (s0C * E_0C * E_0C);
    
    /*  dE_1C / dt  */
    ydot[7]  = -(s1C + d1C_0C + d1C_01C + d1C_1B + d1C_1A + d1C_1D + d1C_1E + x1C) * E_1C + (d1C_0C * E_0C + d1C_01C * E_01C + d1C_1B * E_1B + d1C_1A * E_1A + d1C_1D * E_1D + d1C_1E * E_1E) + x1C + (s1C * E_1C * E_1C);
    
    /*  dE_01C/ dt  */
    ydot[8] = -(s01C + d01C_0C + d01C_1C + d01C_01B + d01C_01A + d01C_01D + d01C_01E + x01C) * E_01C + (d01C_0C * E_0C + d01C_1C * E_1C + d01C_01B * E_01B + d01C_01A * E_01A + d01C_01D * E_01D + d01C_01E * E_01E) + x01C + (s01C * E_01C * E_01C);

    
    /*  dE_0D / dt  */
    ydot[9]  = -(s0D + d0D_1D + d0D_01D + d0D_0B + d0D_0A + d0D_0C + d0D_0E + x0D) * E_0D + (d0D_1D * E_1D + d0D_01D * E_01D + d0D_0B * E_0B + d0D_0A * E_0A + d0D_0C * E_0C + d0D_0E * E_0E) + x0D + (s0D * E_0D * E_0D);
    
    /*  dE_1D / dt  */
    ydot[10] = -(s1D + d1D_0D + d1D_01D + d1D_1B + d1D_1A + d1D_1C + d1D_1E + x1D) * E_1D + (d1D_0D * E_0D + d1D_01D * E_01D + d1D_1B * E_1B + d1D_1A * E_1A + d1D_1C * E_1C + d1D_1E * E_1E) + x1D + (s1D * E_1D * E_1D);
    
    /*  dE_01D/ dt  */
    ydot[11] = -(s01D + d01D_0D + d01D_1D + d01D_01B + d01D_01A + d01D_01C + d01D_01E + x01D) * E_01D + (d01D_0D * E_0D + d01D_1D * E_1D + d01D_01B * E_01B + d01D_01A * E_01A + d01D_01C * E_01C + d01D_01E * E_01E) + x01D + (s01D * E_01D * E_01D);
    

    /*  dE_0E / dt  */
    ydot[12] = -(s0E + d0E_1E + d0E_01E + d0E_0B + d0E_0A + d0E_0C + d0E_0D + x0E) * E_0E + (d0E_1E * E_1E + d0E_01E * E_01E + d0E_0B * E_0B + d0E_0A * E_0A + d0E_0C * E_0C + d0E_0D * E_0D) + x0E + (s0E * E_0E * E_0E);
    
    /*  dE_1E / dt  */
    ydot[13] = -(s1E + d1E_0E + d1E_01E + d1E_1B + d1E_1A + d1E_1C + d1E_1D + x1E) * E_1E + (d1E_0E * E_0E + d1E_01E * E_01E + d1E_1B * E_1B + d1E_1A * E_1A + d1E_1C * E_1C + d1E_1D * E_1D) + x1E + (s1E * E_1E * E_1E);
    
    /*  dE_01E/ dt  */
    ydot[14] = -(s01E + d01E_0E + d01E_1E + d01E_01B + d01E_01A + d01E_01C + d01E_01D + x01E) * E_01E + (d01E_0E * E_0E + d01E_1E * E_1E + d01E_01B * E_01B + d01E_01A * E_01A + d01E_01C * E_01C + d01E_01D * E_01D) + x01E + (s01E * E_01E * E_01E);

    
    
    
    /* The D's */
    
    /*  dD_N0A / dt  */
    ydot[15] = -(s0A + d0A_1A + d0A_01A + d0A_0B + d0A_0C + d0A_0D + d0A_0E + x0A) * D_N0A + (d0A_1A * D_N1A + d0A_01A * D_N01A + d0A_0B * D_N0B + d0A_0C * D_N0C + d0A_0D * D_N0D + d0A_0E * D_N0E) + s0A * (D_N0A * E_0A + D_N0A * E_0A);
    
    /*  dD_N1A / dt  */
    ydot[16] = -(s1A + d1A_0A + d1A_01A + d1A_1B + d1A_1C + d1A_1D + d1A_1E + x1A) * D_N1A + (d1A_0A * D_N0A + d1A_01A * D_N01A + d1A_1B * D_N1B + d1A_1C * D_N1C + d1A_1D * D_N1D + d1A_1E * D_N1E) + s1A * (D_N1A * E_1A + D_N1A * E_1A);
    
    /*  dD_N01A / dt  */
    ydot[17] = -(s01A + d01A_0A + d01A_1A + d01A_01B + d01A_01C + d01A_01D + d01A_01E + x01A) * D_N01A + (d01A_0A * D_N0A + d01A_1A * D_N1A + d01A_01B * D_N01B + d01A_01C * D_N01C + d01A_01D * D_N01D + d01A_01E * D_N01E) + s01A * (D_N01A * E_01A + D_N01A * E_01A);

    
    /*  dD_N0B / dt  */
    ydot[18] = -(s0B + d0B_1B + d0B_01B + d0B_0A + d0B_0C + d0B_0D + d0B_0E + x0B) * D_N0B + (d0B_1B * D_N1B + d0B_01B * D_N01B + d0B_0A * D_N0A + d0B_0C * D_N0C + d0B_0D * D_N0D + d0B_0E * D_N0E) + s0B * (D_N0B * E_0B + D_N0B * E_0B);
    
    /*  dD_N1B / dt  */
    ydot[19] = -(s1B + d1B_0B + d1B_01B + d1B_1A + d1B_1C + d1B_1D + d1B_1E + x1B) * D_N1B + (d1B_0B * D_N0B + d1B_01B * D_N01B + d1B_1A * D_N1A + d1B_1C * D_N1C + d1B_1D * D_N1D + d1B_1E * D_N1E) + s1B * (D_N1B * E_1B + D_N1B * E_1B);
    
    /*  dD_N01B / dt  */
    ydot[20] = -(s01B + d01B_0B + d01B_1B + d01B_01A + d01B_01C + d01B_01D + d01B_01E + x01B) * D_N01B + (d01B_0B * D_N0B + d01B_1B * D_N1B + d01B_01A * D_N01A + d01B_01C * D_N01C + d01B_01D * D_N01D + d01B_01E * D_N01E) + s01B * (D_N01B * E_01B + D_N01B * E_01B);

    
    /*  dD_N0C / dt  */
    ydot[21] = -(s0C + d0C_1C + d0C_01C + d0C_0B + d0C_0A + d0C_0D + d0C_0E + x0C) * D_N0C + (d0C_1C * D_N1C + d0C_01C * D_N01C + d0C_0B * D_N0B + d0C_0A * D_N0A + d0C_0D * D_N0D + d0C_0E * D_N0E) + s0C * (D_N0C * E_0C + D_N0C * E_0C);
    
    /*  dD_N1C / dt  */
    ydot[22] = -(s1C + d1C_0C + d1C_01C + d1C_1B + d1C_1A + d1C_1D + d1C_1E + x1C) * D_N1C + (d1C_0C * D_N0C + d1C_01C * D_N01C + d1C_1B * D_N1B + d1C_1A * D_N1A + d1C_1D * D_N1D + d1C_1E * D_N1E) + s1C * (D_N1C * E_1C + D_N1C * E_1C);
    
    /*  dD_N01C / dt  */
    ydot[23] = -(s01C + d01C_0C + d01C_1C + d01C_01A + d01C_01B + d01C_01D + d01C_01E + x01C) * D_N01C + (d01C_0C * D_N0C + d01C_1C * D_N1C + d01C_01B * D_N01B + d01C_01A * D_N01A + d01C_01D * D_N01D + d01C_01E * D_N01E) + s01C * (D_N01C * E_01C + D_N01C * E_01C);

    
    /*  dD_N0D / dt  */
    ydot[24] = -(s0D + d0D_1D + d0D_01D + d0D_0B + d0D_0A + d0D_0C + d0D_0E + x0D) * D_N0D + (d0D_1D * D_N1D + d0D_01D * D_N01D + d0D_0B * D_N0B + d0D_0A * D_N0A + d0D_0C * D_N0C + d0D_0E * D_N0E) + s0D * (D_N0D * E_0D + D_N0D * E_0D);
    
    /*  dD_N1D / dt  */
    ydot[25] = -(s1D + d1D_0D + d1D_01D + d1D_1B + d1D_1A + d1D_1C + d1D_1E + x1D) * D_N1D + (d1D_0D * D_N0D + d1D_01D * D_N01D + d1D_1B * D_N1B + d1D_1A * D_N1A + d1D_1C * D_N1C + d1D_1E * D_N1E) + s1D * (D_N1D * E_1D + D_N1D * E_1D);
    
    /*  dD_N01D / dt  */
    ydot[26] = -(s01D + d01D_0D + d01D_1D + d01D_01A + d01D_01B + d01D_01C + d01D_01E + x01D) * D_N01D + (d01D_0D * D_N0D + d01D_1D * D_N1D + d01D_01B * D_N01B + d01D_01A * D_N01A + d01D_01C * D_N01C + d01D_01E * D_N01E) + s01D * (D_N01D * E_01D + D_N01D * E_01D);

    
    /*  dD_N0E / dt  */
    ydot[27] = -(s0E + d0E_1E + d0E_01E + d0E_0B + d0E_0A + d0E_0C + d0E_0D + x0E) * D_N0E + (d0E_1E * D_N1E + d0E_01E * D_N01E + d0E_0B * D_N0B + d0E_0A * D_N0A + d0E_0C * D_N0C + d0E_0D * D_N0D) + s0E * (D_N0E * E_0E + D_N0E * E_0E);
    
    /*  dD_N1E / dt  */
    ydot[28] = -(s1E + d1E_0E + d1E_01E + d1E_1B + d1E_1A + d1E_1C + d1E_1D + x1E) * D_N1E + (d1E_0E * D_N0E + d1E_01E * D_N01E + d1E_1B * D_N1B + d1E_1A * D_N1A + d1E_1C * D_N1C + d1E_1D * D_N1D) + s1E * (D_N1E * E_1E + D_N1E * E_1E);
    
    /*  dD_N01E / dt  */
    ydot[29] = -(s01E + d01E_0E + d01E_1E + d01E_01A + d01E_01B + d01E_01C + d01E_01D + x01E) * D_N01E + (d01E_0E * D_N0E + d01E_1E * D_N1E + d01E_01B * D_N01B + d01E_01A * D_N01A + d01E_01C * D_N01C + d01E_01D * D_N01D) + s01E * (D_N01E * E_01E + D_N01E * E_01E);

}













