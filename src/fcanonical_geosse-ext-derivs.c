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

static double params_fgeosse[NUMELEMENTS];


void initmod_fgeosse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_fgeosse);
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


void fclasse_geosse_equivalent_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E_0 = y[0];
    double E_1 = y[1];
    double E_2 = y[2];
    
    double D_N0 = y[3];
    double D_N1 = y[4];
    double D_N2 = y[5];
    
    double
    sA  = params_fgeosse[0],     /* speciation within region A */
    sB  = params_fgeosse[1],     /* speciation within region B */
    sAB = params_fgeosse[2],     /* between-region speciation  */
    xA  = params_fgeosse[3],     /* extinction from region A   */
    xB  = params_fgeosse[4],     /* extinction from region B   */
    d0_1  = params_fgeosse[5],   /* jumps from 0 to 1          */
    d0_01  = params_fgeosse[6],  /* dispersal from A to AB     */
    d1_0 = params_fgeosse[7],    /* jumps from 1 to 0          */
    d1_01 = params_fgeosse[8],   /* dispersal from B to AB     */
    d01_0 = params_fgeosse[9],   /* true extirpation rate      */
    d01_1 = params_fgeosse[10];  /* true extirpation rate      */
    
    /*  dE_2 / dt  */
    ydot[0] = -(sA + d0_1 + d0_01 + xA) * E_0 + (d0_1 * E_1 + d0_01 * E_2) + xA + (sA * E_0 * E_0);
    
    /*  dE_3 / dt  */
    ydot[1] = -(sB + d1_0 + d1_01 + xB) * E_1 + (d1_0 * E_0 + d1_01 * E_2) + xB + (sB * E_1 * E_1);
    
    /*  dE_1 / dt  */
    ydot[2] = -(sAB + sA + sB + d01_0 + d01_1) * E_2 + (d01_0 * E_0 + d01_1 * E_1) + sA * E_0 * E_2 + sB * E_2 * E_1 + sAB * E_0 * E_1;
    
    /*  dD_N2 / dt  */
    ydot[3] = -(sA + d0_1 + d0_01 + xA) * D_N0 + (d0_1 * D_N1 + d0_01 * D_N2) + sA * (D_N0 * E_0 + D_N0 * E_0);
    
    /*  dD_N3 / dt  */
    ydot[4] = -(sB + d1_0 + d1_01 + xB) * D_N1 + (d1_0 * D_N0 + d1_01 * D_N2) + sB * (D_N1 * E_1 + D_N1 * E_1);
    
    /*  dD_N1 / dt  */
    ydot[5] = -(sAB + sA + sB + d01_0 + d01_1) * D_N2 + (d01_0 * D_N0 + d01_1 * D_N1) + sAB * (D_N0 * E_1 + D_N1 * E_0) + sA * (E_0 * D_N2 + E_2 * D_N0) + sB * (E_1 * D_N2 + E_2 * D_N1);
    
}



