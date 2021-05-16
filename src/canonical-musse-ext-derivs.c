/*
 *  canonical-musse-ext-derivs_c
 *
 *
 *  Created by Jeremy Beaulieu 4/25/2018
 *  Copyright 2018 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 21

static double params_musse[NUMELEMENTS];


void initmod_musse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_musse);
}


void musse_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E00A = y[0];
    double E01A = y[1];
    double E10A = y[2];
    double E11A = y[3];
    
    double D00A = y[4];
    double D01A = y[5];
    double D10A = y[6];
    double D11A = y[7];
    
    double
    lambda00A = params_musse[0],
    lambda01A = params_musse[1],
    lambda10A = params_musse[2],
    lambda11A = params_musse[3],
    mu00A = params_musse[4],
    mu01A = params_musse[5],
    mu10A = params_musse[6],
    mu11A = params_musse[7],
    q00A_01A = params_musse[8],
    q00A_10A = params_musse[9],
    q00A_11A = params_musse[10],
    q01A_00A = params_musse[11],
    q01A_10A = params_musse[12],
    q01A_11A = params_musse[13],
    q10A_00A = params_musse[14],
    q10A_01A = params_musse[15],
    q10A_11A = params_musse[16],
    q11A_00A = params_musse[17],
    q11A_01A = params_musse[18],
    q11A_10A = params_musse[19],
    psi = params_musse[20];
    
    /* The E's */
    ydot[0] = mu00A - (lambda00A + psi + (q00A_01A + q00A_10A + q00A_11A) + mu00A) * E00A + lambda00A*E00A*E00A + (q00A_01A*E01A + q00A_10A*E10A + q00A_11A*E11A);
    
    ydot[1] = mu01A - (lambda01A + psi + (q01A_00A + q01A_10A + q01A_11A) + mu01A) * E01A + lambda01A*E01A*E01A + (q01A_00A*E00A + q01A_10A*E10A + q01A_11A*E11A);
    
    ydot[2] = mu10A - (lambda10A + psi + (q10A_00A + q10A_01A + q10A_11A) + mu10A) * E10A + lambda10A*E10A*E10A + (q10A_00A*E00A + q10A_01A*E01A + q10A_11A*E11A);
    
    ydot[3] = mu11A - (lambda11A + psi + (q11A_00A + q11A_01A + q11A_10A) + mu11A) * E11A + lambda11A*E11A*E11A + (q11A_00A*E00A + q11A_01A*E01A + q11A_10A*E10A);
    
    
    /* The D's */
    ydot[4] =  - (lambda00A + psi + (q00A_01A + q00A_10A + q00A_11A) + mu00A) * D00A + 2*lambda00A*E00A*D00A + (q00A_01A*D01A + q00A_10A*D10A + q00A_11A*D11A);
    
    ydot[5] =  - (lambda01A + psi + (q01A_00A + q01A_10A + q01A_11A) + mu01A) * D01A + 2*lambda01A*E01A*D01A + (q01A_00A*D00A + q01A_10A*D10A + q01A_11A*D11A);
    
    ydot[6] =  - (lambda10A + psi + (q10A_00A + q10A_01A + q10A_11A) + mu10A) * D10A + 2*lambda10A*E10A*D10A + (q10A_00A*D00A + q10A_01A*D01A + q10A_11A*D11A);
    
    ydot[7] =  - (lambda11A + psi + (q11A_00A + q11A_01A + q11A_10A) + mu11A) * D11A + 2*lambda11A*E11A*D11A + (q11A_00A*D00A + q11A_01A*D01A + q11A_10A*D10A);
}



