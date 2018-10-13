/*
 *  notclasse-ext-derivs_c
 *
 *
 *  Created by Jeremy Beaulieu 10/15/2017
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

static double params_fnoclass[NUMELEMENTS];


void initmod_fnoclass(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_fnoclass);
}



void fnotclasse_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E_0 = y[0];
    double E_1 = y[1];
    double E_2 = y[2];
    
    double D_N0 = y[3];
    double D_N1 = y[4];
    double D_N2 = y[5];

    double
    s0  = params_fnoclass[0],     /* speciation within region 0 */
    s1  = params_fnoclass[1],     /* speciation within region 1 */
    s01 = params_fnoclass[2],     /* between-region speciation  */
    x0  = params_fnoclass[3],     /* extinction from region 0   */
    x1  = params_fnoclass[4],     /* extinction from 1          */
    x01 = 0,                      /* extinction from 01         */

    d0_1  = params_fnoclass[5],   /* jumps from 0 to 1          */
    d0_01 = params_fnoclass[6],   /* dispersal from A to AB     */
    d1_0  = params_fnoclass[7],   /* jumps from 1 to 0          */
    d1_01 = params_fnoclass[8],   /* dispersal from B to AB     */
    d01_0 = params_fnoclass[9],   /* true extirpation rate      */
    d01_1 = params_fnoclass[10];  /* true extirpation rate      */

    ydot[0] = -(x0 + d0_1 + d0_01 + s0) * E_0 + (s0 * E_0 * E_0) + x0 + (d0_1 * E_1 + d0_01 * E_2);
    ydot[1] = -(x1 + d1_0 + d1_01 + s1) * E_1 + (s1 * E_1 * E_1) + x1 + (d1_0 * E_0 + d1_01 * E_2);
    ydot[2] = -(x01 + d01_0 + d01_1 + s01) * E_2 + (s01 * E_2 * E_2) + x01 + (d01_0 * E_0 + d01_1 * E_1);
    
    ydot[3] = -(x0 + d0_1 + d0_01 + s0) * D_N0 + 2 * s0 * E_0 * D_N0 + (d0_1 * D_N1 + d0_01 * D_N2);
    ydot[4] = -(x1 + d1_0 + d1_01 + s1) * D_N1 + 2 * s1 * E_1 * D_N1 + (d1_0 * D_N0 + d1_01 * D_N2);
    ydot[5] = -(x01 + d01_0 + d01_1 + s01) * D_N2 + 2 * s01 * E_2 * D_N2 + (d01_0 * D_N0 + d01_1 * D_N1);

}
