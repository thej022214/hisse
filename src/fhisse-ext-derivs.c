/*
 *  fhisse-ext-derivs_c
 *  
 *
 *  Created by Jeremy Beaulieu 8/3/2019
 *  Copyright 2019 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 48

static double params_fhisse[NUMELEMENTS];


void initmod_fhisse(void (* odeparms)(int *, double *)){
	int N = NUMELEMENTS;
	odeparms(&N, params_fhisse);
}


void maddison_DE_fhisse(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
	
    double E0A = y[0], E1A = y[1], E0B = y[2], E1B = y[3], E0C = y[4], E1C = y[5], E0D = y[6], E1D = y[7];
	double D0A = y[8], D1A = y[9], D0B = y[10], D1B = y[11], D0C = y[12], D1C = y[13], D0D = y[14], D1D = y[15];

    double 
	lambda0A = params_fhisse[0],
    lambda1A = params_fhisse[1],
    mu0A = params_fhisse[2],
    mu1A = params_fhisse[3],
    q0A1A = params_fhisse[4],
    q1A0A = params_fhisse[5],
    q0A0B = params_fhisse[6],
    q0A0C = params_fhisse[7],
    q0A0D = params_fhisse[8],
    q1A1B = params_fhisse[9],
    q1A1C = params_fhisse[10],
    q1A1D = params_fhisse[11],

    lambda0B = params_fhisse[12],
    lambda1B = params_fhisse[13],
    mu0B = params_fhisse[14],
    mu1B = params_fhisse[15],
    q0B1B = params_fhisse[16],
    q1B0B = params_fhisse[17],
    q0B0A = params_fhisse[18],
    q0B0C = params_fhisse[19],
    q0B0D = params_fhisse[20],
    q1B1A = params_fhisse[21],
    q1B1C = params_fhisse[22],
    q1B1D = params_fhisse[23],

    lambda0C = params_fhisse[24],
    lambda1C = params_fhisse[25],
    mu0C = params_fhisse[26],
    mu1C = params_fhisse[27],
    q0C1C = params_fhisse[28],
    q1C0C = params_fhisse[29],
    q0C0A = params_fhisse[30],
    q0C0B = params_fhisse[31],
    q0C0D = params_fhisse[32],
    q1C1A = params_fhisse[33],
    q1C1B = params_fhisse[34],
    q1C1D = params_fhisse[35],

    lambda0D = params_fhisse[36],
    lambda1D = params_fhisse[37],
    mu0D = params_fhisse[38],
    mu1D = params_fhisse[39],
    q0D1D = params_fhisse[40],
    q1D0D = params_fhisse[41],
	q0D0A = params_fhisse[42],
	q0D0B = params_fhisse[43],
	q0D0C = params_fhisse[44],
	q1D1A = params_fhisse[45],
	q1D1B = params_fhisse[46],
    q1D1C = params_fhisse[47];
	
	ydot[0] = -(mu0A + (q0A0B + q0A0C + q0A0D + q0A1A) + lambda0A) * E0A + lambda0A * E0A * E0A + mu0A + (q0A0B*E0B + q0A0C*E0C + q0A0D*E0D + q0A1A*E1A);
    ydot[1] = -(mu1A + (q1A0A + q1A1B + q1A1C + q1A1D) + lambda1A) * E1A + lambda1A * E1A * E1A + mu1A + (q1A0A*E0A + q1A1B*E1B + q1A1C*E1C + q1A1D*E1D);
    ydot[2] = -(mu0B + (q0B0A + q0B0C + q0B0D + q0B1B) + lambda0B) * E0B + lambda0B * E0B * E0B + mu0B + (q0B0A*E0A + q0B0C*E0C + q0B0D*E0D + q0B1B*E1B);
    ydot[3] = -(mu1B + (q1B0B + q1B1A + q1B1C + q1B1D) + lambda1B) * E1B + lambda1B * E1B * E1B + mu1B + (q1B0B*E0B + q1B1A*E1A + q1B1C*E1C + q1B1D*E1D);
    ydot[4] = -(mu0C + (q0C0A + q0C0B + q0C0D + q0C1C) + lambda0C) * E0C + lambda0C * E0C * E0C + mu0C + (q0C0A*E0A + q0C0B*E0B + q0C0D*E0D + q0C1C*E1C);
    ydot[5] = -(mu1C + (q1C0C + q1C1A + q1C1B + q1C1D) + lambda1C) * E1C + lambda1C * E1C * E1C + mu1C + (q1C0C*E0C + q1C1A*E1A + q1C1B*E1B + q1C1D*E1D);
    ydot[6] = -(mu0D + (q0D0A + q0D0B + q0D0C + q0D1D) + lambda0D) * E0D + lambda0D * E0D * E0D + mu0D + (q0D0A*E0A + q0D0B*E0B + q0D0C*E0C + q0D1D*E1D);
	ydot[7] = -(mu1D + (q1D0D + q1D1A + q1D1B + q1D1C) + lambda1D) * E1D + lambda1D * E1D * E1D + mu1D + (q1D0D*E0D + q1D1A*E1A + q1D1B*E1B + q1D1C*E1C);

	ydot[8]  = -(mu0A + (q0A0B + q0A0C + q0A0D + q0A1A) + lambda0A) * D0A + 2 * lambda0A * E0A * D0A + (q0A0B*D0B + q0A0C*D0C + q0A0D*D0D + q0A1A*D1A);
    ydot[9] = -(mu1A + (q1A0A + q1A1B + q1A1C + q1A1D) + lambda1A) * D1A + 2 * lambda1A * E1A * D1A + (q1A0A*D0A + q1A1B*D1B + q1A1C*D1C + q1A1D*D1D);
    ydot[10]  = -(mu0B + (q0B0A + q0B0C + q0B0D + q0B1B) + lambda0B) * D0B + 2 * lambda0B * E0B * D0B + (q0B0A*D0A + q0B0C*D0C + q0B0D*D0D + q0B1B*D1B);
    ydot[11] = -(mu1B + (q1B0B + q1B1A + q1B1C + q1B1D) + lambda1B) * D1B + 2 * lambda1B * E1B * D1B + (q1B0B*D0B + q1B1A*D1A + q1B1C*D1C + q1B1D*D1D);
    ydot[12] = -(mu0C + (q0C0A + q0C0B + q0C0D + q0C1C) + lambda0C) * D0C + 2 * lambda0C * E0C * D0C + (q0C0A*D0A + q0C0B*D0B + q0C0D*D0D + q0C1C*D1C);
    ydot[13] = -(mu1C + (q1C0C + q1C1A + q1C1B + q1C1D) + lambda1C) * D1C + 2 * lambda1C * E1C * D1C + (q1C0C*D0C + q1C1A*D1A + q1C1B*D1B + q1C1D*D1D);
    ydot[14] = -(mu0D + (q0D0A + q0D0B + q0D0C + q0D1D) + lambda0D) * D0D + 2 * lambda0D * E0D * D0D + (q0D0A*D0A + q0D0B*D0B + q0D0C*D0C + q0D1D*D1D);
	ydot[15] = -(mu1D + (q1D0D + q1D1A + q1D1B + q1D1C) + lambda1D) * D1D + 2 * lambda1D * E1D * D1D + (q1D0D*D0D + q1D1A*D1A + q1D1B*D1B + q1D1C*D1C);
}


