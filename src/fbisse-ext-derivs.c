/*
 *  fbisse-ext-derivs_c
 *
 *
 *  Created by Jeremy Beaulieu 8/2/2015
 *  Copyright 2019 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 7

static double params_fbisse[NUMELEMENTS];

void initmod_fbisse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_fbisse);
}


void maddison_DE_fbisse(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E0 = y[0], E1 = y[1];
    double D0 = y[2], D1 = y[3];
    
    double
    lambda0 = params_fbisse[0],
    lambda1 = params_fbisse[1],
    mu0 = params_fbisse[2],
    mu1 = params_fbisse[3],
    q01 = params_fbisse[4],
    q10 = params_fbisse[5],
    psi = params_fbisse[6];
    
    ydot[0] = -(mu0 + psi + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1;
    ydot[1] = -(mu1 + psi + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0;
    ydot[2] = -(mu0 + psi + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1;
    ydot[3] = -(mu1 + psi + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0;
}


void maddison_DE_strat_fbisse(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E0 = y[0], E1 = y[1];
    double D0 = y[2], D1 = y[3];
    
    double
    lambda0 = params_fbisse[0],
    lambda1 = params_fbisse[1],
    mu0 = params_fbisse[2],
    mu1 = params_fbisse[3],
    q01 = params_fbisse[4],
    q10 = params_fbisse[5],
    psi = params_fbisse[6];
    
    ydot[0] = -(mu0 + psi + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1;
    ydot[1] = -(mu1 + psi + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0;
    ydot[2] = -(mu0 + psi + q01 + lambda0) * D0 + q01 * D1;
    ydot[3] = -(mu1 + psi + q10 + lambda1) * D1 + q10 * D0;
}



