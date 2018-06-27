/*
 *  canonical_geosse-ext-derivs_c
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
#define NUMELEMENTS 11

static double params_geosse[NUMELEMENTS];


void initmod_geosse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_geosse);
}


//void geosse_canonical_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
//    double E_2 = y[0];
//    double E_3 = y[1];
//    double E_1 = y[2];
    
//    double D_N2 = y[3];
//    double D_N3 = y[4];
//    double D_N1 = y[5];
    
//    double
//    sA  = params_geohisse[0],     /* speciation within region A */
//    sB  = params_geohisse[1],     /* speciation within region B */
//    sAB = params_geohisse[2],     /* between-region speciation  */
//    xA  = params_geohisse[3],     /* extinction from region A   */
//    xB  = params_geohisse[4],     /* extinction from region B   */
//    dA  = params_geohisse[5],     /* dispersal from A to B      */
//    dB  = params_geohisse[6];     /* dispersal from B to A      */

    /*  dE_2 / dt  */
//    ydot[0] = -(sA + dA + xA) * E_2
//    + xA + dA * E_1 + sA * E_2 * E_2;
    
    /*  E_3 / dt  */
//    ydot[1] = -(sB + dB + xB) * E_3
//    + xB + dB * E_1 + sB * E_3 * E_3;
   
    /*  dE_1 / dt  */
//    ydot[2] = -(sA + sB + xA + xB + sAB) * E_1
//    + xA * E_3 + xB * E_2
//    + sA * E_1 * E_2 + sB * E_1 * E_3 + sAB * E_2 * E_3;

    /*  dD_N2 / dt  */
//    ydot[3] = -(sA + dA + xA) * D_N2
//    + dA * D_N1 + 2 * sA * D_N2 * E_2;
    
    /*  dD_N3 / dt  */
//    ydot[4] = -(sB + dB + xB) * D_N3
//    + dB * D_N1 + 2 * sB * D_N3 * E_3;

    /*  dD_N1 / dt  */
//    ydot[5] = -(sA + sB + sAB + xA + xB) * D_N1
//    + xA * D_N3 + xB * D_N2
//    + sA * (E_2 * D_N1 + E_1 * D_N2)
//    + sB * (E_3 * D_N1 + E_1 * D_N3)
//   + sAB * (E_2 * D_N3 + E_3 * D_N2);
    
//}


void classe_geosse_equivalent_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E_0 = y[0];
    double E_1 = y[1];
    double E_2 = y[2];
    
    double D_N0 = y[3];
    double D_N1 = y[4];
    double D_N2 = y[5];
    
    double
    s00A  = params_geosse[0],     /* speciation within region A */
    s11A  = params_geosse[1],     /* speciation within region B */
    s01A = params_geosse[2],     /* between-region speciation  */
    x00A  = params_geosse[3],     /* extinction from region A   */
    x11A  = params_geosse[4],     /* extinction from region B   */
    d00A_11A  = params_geosse[5],   /* jumps from 0 to 1          */
    d00A_01A  = params_geosse[6],  /* dispersal from A to AB     */
    d11A_00A = params_geosse[7],    /* jumps from 1 to 0          */
    d11A_01A = params_geosse[8],   /* dispersal from B to AB     */
    d01A_00A = params_geosse[9],   /* true extirpation rate      */
    d01A_11A = params_geosse[10];  /* true extirpation rate      */
    
    ydot[0] =  -(s00A + d00A_11A + d00A_01A +x00A) * E00A + (d00A_11A*E11A + d00A_01A*E01A) + x00A+s00A*E00A*E00A;
    
    ydot[1] =  -(s11A + (d11A_00A + d11A_01A)+x11A) * E11A + (d11A_00A*E00A + d11A_01A*E01A) + x11A+s11A*E11A*E11A;
    
    ydot[2] =  -(s01A + s00A + s11A + (d01A_00A + d01A_11A)) * E01A + (d01A_00A*E00A + d01A_11A*E11A) + s00A*E01A*E00A + s11A*E01A*E11A + s01A*E11A*E00A;
    
    ydot[3] =  -(s00A+(d00A_11A + d00A_01A)+x00A) * D00A + (d00A_11A*D11A + d00A_01A*D01A) + 2*s00A*E00A*D00A;
    
    ydot[4] =  -(s11A+(d11A_00A + d11A_01A)+x11A) * D11A + (d11A_00A*D00A + d11A_01A*D01A) + 2*s11A*E11A*D11A;
    
    ydot[5] =  -(s01A + s00A + s11A + (d01A_00A + d01A_11A)) * D01A + (d01A_00A*D00A + d01A_11A*D11A) + s00A*(E00A*D01A+E01A*D00A) + s11A*(E11A*D01A + E01A*D11A) + s01A*(E00A*D11A + E11A*D00A);
    
}

