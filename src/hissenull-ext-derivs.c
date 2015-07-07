/*
 *  hissenull-ext-derivs_c
 *  
 *
 *  Created by Jeremy Beaulieu 5/26/2015
 *  Copyright 2015 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 51

static double params_hisse_null[NUMELEMENTS]; 


void initmod_hisse_null(void (* odeparms)(int *, double *)){
	int N = NUMELEMENTS;
	odeparms(&N, params_hisse_null);
}


//If we are calling this function then we are at a node:
void set_birth_hisse_null_void(double *birth, double *point_time, double *tot_time, double *x_turnover0A, double *x_eps0A, double *x_turnover0B, double *x_eps0B,double *x_turnover0C, double *x_eps0C, double *x_turnover0D, double *x_eps0D, double *x_turnover1A, double *x_eps1A, double *x_turnover1B, double *x_eps1B,double *x_turnover1C, double *x_eps1C, double *x_turnover1D, double *x_eps1D, double *q0B0A, double *q0C0A, double *q0D0A, double *q1A0A, double *q0A0B, double *q0C0B, double *q0D0B, double *q1B0B, double *q0A0C, double *q0B0C, double *q0D0C, double *q1C0C, double *q0A0D, double *q0B0D, double *q0C0D, double *q1D0D, double *q0A1A, double *q1B1A, double *q1C1A, double *q1D1A, double *q0B1B, double *q1A1B, double *q1C1B, double *q1D1B, double *q0C1C, double *q1A1C, double *q1B1C, double *q1D1C, double *q0D1D, double *q1A1D, double *q1B1D, double *q1C1D, double *focal_edge_length, double *tipward_age, int *state){

	if(*state == 0){
		*birth = *x_turnover0A / (1 + *x_eps0A);
	}
	if(*state == 1){
		*birth = *x_turnover0B / (1 + *x_eps0B);
	}
	if(*state == 2){
		*birth = *x_turnover0C / (1 + *x_eps0C);
	}
	if(*state == 3){
		*birth = *x_turnover0D / (1 + *x_eps0D);
	}
	if(*state == 4){
		*birth = *x_turnover1A / (1 + *x_eps1A);
	}
	if(*state == 5){
		*birth = *x_turnover1B / (1 + *x_eps1B);
	}
	if(*state == 6){
		*birth = *x_turnover1C / (1 + *x_eps1C);
	}
	if(*state == 7){
		*birth = *x_turnover1D / (1 + *x_eps1D);
	}
}


//If we are calling this function then we are somewhere along a branch:
double set_birth_hisse_null(double *point_time, double *tot_time, double *x_turnover0A, double *x_eps0A, double *x_turnover0B, double *x_eps0B,double *x_turnover0C, double *x_eps0C, double *x_turnover0D, double *x_eps0D, double *x_turnover1A, double *x_eps1A, double *x_turnover1B, double *x_eps1B,double *x_turnover1C, double *x_eps1C, double *x_turnover1D, double *x_eps1D, double *q0B0A, double *q0C0A, double *q0D0A, double *q1A0A, double *q0A0B, double *q0C0B, double *q0D0B, double *q1B0B, double *q0A0C, double *q0B0C, double *q0D0C, double *q1C0C, double *q0A0D, double *q0B0D, double *q0C0D, double *q1D0D, double *q0A1A, double *q1B1A, double *q1C1A, double *q1D1A, double *q0B1B, double *q1A1B, double *q1C1B, double *q1D1B, double *q0C1C, double *q1A1C, double *q1B1C, double *q1D1C, double *q0D1D, double *q1A1D, double *q1B1D, double *q1C1D, double *focal_edge_length, double *tipward_age, int *state){

	if(*state == 0){
		double birth_rate = *x_turnover0A / (1 + *x_eps0A);
		return(birth_rate);
	}
	if(*state == 1){
		double birth_rate = *x_turnover0B / (1 + *x_eps0B);
		return(birth_rate);
	}
	if(*state == 2){
		double birth_rate = *x_turnover0C / (1 + *x_eps0C);
		return(birth_rate);
	}
	if(*state == 3){
		double birth_rate = *x_turnover0D / (1 + *x_eps0D);
		return(birth_rate);
	}
	if(*state == 4){
		double birth_rate = *x_turnover1A / (1 + *x_eps1A);
		return(birth_rate);
	}
	if(*state == 5){
		double birth_rate = *x_turnover1B / (1 + *x_eps1B);
		return(birth_rate);
	}
	if(*state == 6){
		double birth_rate = *x_turnover1C / (1 + *x_eps1C);
		return(birth_rate);
	}
	if(*state == 7){
		double birth_rate = *x_turnover1D / (1 + *x_eps1D);
		return(birth_rate);
	}	
	return 0;
}


//If we are calling this function then we are at a node:
void set_death_hisse_null_void(double *death, double *point_time, double *tot_time, double *x_turnover0A, double *x_eps0A, double *x_turnover0B, double *x_eps0B,double *x_turnover0C, double *x_eps0C, double *x_turnover0D, double *x_eps0D, double *x_turnover1A, double *x_eps1A, double *x_turnover1B, double *x_eps1B,double *x_turnover1C, double *x_eps1C, double *x_turnover1D, double *x_eps1D, double *q0B0A, double *q0C0A, double *q0D0A, double *q1A0A, double *q0A0B, double *q0C0B, double *q0D0B, double *q1B0B, double *q0A0C, double *q0B0C, double *q0D0C, double *q1C0C, double *q0A0D, double *q0B0D, double *q0C0D, double *q1D0D, double *q0A1A, double *q1B1A, double *q1C1A, double *q1D1A, double *q0B1B, double *q1A1B, double *q1C1B, double *q1D1B, double *q0C1C, double *q1A1C, double *q1B1C, double *q1D1C, double *q0D1D, double *q1A1D, double *q1B1D, double *q1C1D, double *focal_edge_length, double *tipward_age, int *state){

	if(*state == 0){
		*death = (*x_turnover0A * *x_eps0A) / (1 + *x_eps0A);
	}
	if(*state == 1){
		*death = (*x_turnover0B * *x_eps0B) / (1 + *x_eps0B);
	}
	if(*state == 2){
		*death = (*x_turnover0C * *x_eps0C) / (1 + *x_eps0C);
	}
	if(*state == 3){
		*death = (*x_turnover0D * *x_eps0D) / (1 + *x_eps0D);
	}
	if(*state == 4){
		*death = (*x_turnover1A * *x_eps1A) / (1 + *x_eps1A);
	}
	if(*state == 5){
		*death = (*x_turnover1B * *x_eps1B) / (1 + *x_eps1B);
	}
	if(*state == 6){
		*death = (*x_turnover1C * *x_eps1C) / (1 + *x_eps1C);
	}
	if(*state == 7){
		*death = (*x_turnover1D * *x_eps1D) / (1 + *x_eps1D);
	}
}


//If we are calling this function then we are somewhere along a branch:
double set_death_hisse_null(double *point_time, double *tot_time, double *x_turnover0A, double *x_eps0A, double *x_turnover0B, double *x_eps0B,double *x_turnover0C, double *x_eps0C, double *x_turnover0D, double *x_eps0D, double *x_turnover1A, double *x_eps1A, double *x_turnover1B, double *x_eps1B,double *x_turnover1C, double *x_eps1C, double *x_turnover1D, double *x_eps1D, double *q0B0A, double *q0C0A, double *q0D0A, double *q1A0A, double *q0A0B, double *q0C0B, double *q0D0B, double *q1B0B, double *q0A0C, double *q0B0C, double *q0D0C, double *q1C0C, double *q0A0D, double *q0B0D, double *q0C0D, double *q1D0D, double *q0A1A, double *q1B1A, double *q1C1A, double *q1D1A, double *q0B1B, double *q1A1B, double *q1C1B, double *q1D1B, double *q0C1C, double *q1A1C, double *q1B1C, double *q1D1C, double *q0D1D, double *q1A1D, double *q1B1D, double *q1C1D, double *focal_edge_length, double *tipward_age, int *state){
	
	if(*state == 0){
		double death_rate = (*x_turnover0A * *x_eps0A) / (1 + *x_eps0A);
		return(death_rate);
	}
	if(*state == 1){
		double death_rate = (*x_turnover0B * *x_eps0B) / (1 + *x_eps0B);
		return(death_rate);
	}
	if(*state == 2){
		double death_rate = (*x_turnover0C * *x_eps0C) / (1 + *x_eps0C);
		return(death_rate);
	}
	if(*state == 3){
		double death_rate = (*x_turnover0D * *x_eps0D) / (1 + *x_eps0D);
		return(death_rate);
	}
	if(*state == 4){
		double death_rate = (*x_turnover1A * *x_eps1A) / (1 + *x_eps1A);
		return(death_rate);
	}
	if(*state == 5){
		double death_rate = (*x_turnover1B * *x_eps1B) / (1 + *x_eps1B);
		return(death_rate);
	}
	if(*state == 6){
		double death_rate = (*x_turnover1C * *x_eps1C) / (1 + *x_eps1C);
		return(death_rate);
	}
	if(*state == 7){
		double death_rate = (*x_turnover1D * *x_eps1D) / (1 + *x_eps1D);
		return(death_rate);
	}	
	return 0;
}


void maddison_DE_hisse_null(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
	
	double E0A = y[0], E0B = y[1], E0C = y[2], E0D = y[3], E1A = y[4], E1B = y[5], E1C = y[6], E1D = y[7];
	double D0A = y[8], D0B = y[9], D0C = y[10], D0D = y[11], D1A = y[12], D1B = y[13], D1C = y[14], D1D = y[15];

	double point_time = *t;
	double tot_time = params_hisse_null[0],

	x_turnover0A = params_hisse_null[1],
	x_turnover0B = params_hisse_null[2],
	x_turnover0C = params_hisse_null[3],
	x_turnover0D = params_hisse_null[4],
	x_turnover1A = params_hisse_null[5],
	x_turnover1B = params_hisse_null[6],
	x_turnover1C = params_hisse_null[7],
	x_turnover1D = params_hisse_null[8],
	
	x_eps0A = params_hisse_null[9],
	x_eps0B = params_hisse_null[10],
	x_eps0C = params_hisse_null[11],
	x_eps0D = params_hisse_null[12],
	x_eps1A = params_hisse_null[13],
	x_eps1B = params_hisse_null[14],
	x_eps1C = params_hisse_null[15],
	x_eps1D = params_hisse_null[16],
	
	q0B0A = params_hisse_null[17],
	q0C0A = params_hisse_null[18],
	q0D0A = params_hisse_null[19],
	q1A0A = params_hisse_null[20],
	q0A0B = params_hisse_null[21],
	q0C0B = params_hisse_null[22],
	q0D0B = params_hisse_null[23],
	q1B0B = params_hisse_null[24],
	q0A0C = params_hisse_null[25],
	q0B0C = params_hisse_null[26],
	q0D0C = params_hisse_null[27],
	q1C0C = params_hisse_null[28],
	q0A0D = params_hisse_null[29],
	q0B0D = params_hisse_null[30],
	q0C0D = params_hisse_null[31],
	q1D0D = params_hisse_null[32],
	q0A1A = params_hisse_null[33],
	q1B1A = params_hisse_null[34], 
	q1C1A = params_hisse_null[35],
	q1D1A = params_hisse_null[36],
	q0B1B = params_hisse_null[37],
	q1A1B = params_hisse_null[38], 
	q1C1B = params_hisse_null[39],
	q1D1B = params_hisse_null[40],
	q0C1C = params_hisse_null[41],
	q1A1C = params_hisse_null[42],
	q1B1C = params_hisse_null[43],
	q1D1C = params_hisse_null[44],
	q0D1D = params_hisse_null[45],
	q1A1D = params_hisse_null[46],
	q1B1D = params_hisse_null[47],
	q1C1D = params_hisse_null[48],
	
	focal_edge_length = params_hisse_null[49],
	tipward_age = params_hisse_null[50];

	int state = 0;
	double birth_rate0A = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate0A = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);

	state = 1;
	double birth_rate0B = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate0B = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	
	state = 2;
	double birth_rate0C = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate0C = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	
	state = 3;
	double birth_rate0D = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate0D = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);

	state = 4;
	double birth_rate1A = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate1A = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	
	state = 5;
	double birth_rate1B = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate1B = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	
	state = 6;
	double birth_rate1C = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate1C = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	
	state = 7;
	double birth_rate1D = set_birth_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);
	double death_rate1D = set_death_hisse_null(&point_time, &tot_time, &x_turnover0A, &x_eps0A, &x_turnover0B, &x_eps0B,&x_turnover0C, &x_eps0C, &x_turnover0D, &x_eps0D, &x_turnover1A, &x_eps1A, &x_turnover1B, &x_eps1B,&x_turnover1C, &x_eps1C, &x_turnover1D, &x_eps1D, &q0B0A, &q0C0A, &q0D0A, &q1A0A, &q0A0B, &q0C0B, &q0D0B, &q1B0B, &q0A0C, &q0B0C, &q0D0C, &q1C0C, &q0A0D, &q0B0D, &q0C0D, &q1D0D, &q0A1A, &q1B1A, &q1C1A, &q1D1A, &q0B1B, &q1A1B, &q1C1B, &q1D1B, &q0C1C, &q1A1C, &q1B1C, &q1D1C, &q0D1D, &q1A1D, &q1B1D, &q1C1D, &focal_edge_length, &tipward_age, &state);

	ydot[0] = -(death_rate0A + (q0A0B + q0A0C + q0A0D + q0A1A) + birth_rate0A) * E0A + birth_rate0A * E0A * E0A + death_rate0A + (q0A0B*E0B + q0A0C*E0C + q0A0D*E0D + q0A1A*E1A);
	ydot[1] = -(death_rate0B + (q0B0A + q0B0C + q0B0D + q0B1B) + birth_rate0B) * E0B + birth_rate0B * E0B * E0B + death_rate0B + (q0B0A*E0A + q0B0C*E0C + q0B0D*E0D + q0B1B*E1B);
	ydot[2] = -(death_rate0C + (q0C0A + q0C0B + q0C0D + q0C1C) + birth_rate0C) * E0C + birth_rate0C * E0C * E0C + death_rate0C + (q0C0A*E0A + q0C0B*E0B + q0C0D*E0D + q0C1C*E1C);
	ydot[3] = -(death_rate0D + (q0D0A + q0D0B + q0D0C + q0D1D) + birth_rate0D) * E0D + birth_rate0D * E0D * E0D + death_rate0D + (q0D0A*E0A + q0D0B*E0B + q0D0C*E0C + q0D1D*E1D);
	ydot[4] = -(death_rate1A + (q1A0A + q1A1B + q1A1C + q1A1D) + birth_rate1A) * E1A + birth_rate1A * E1A * E1A + death_rate1A + (q1A0A*E0A + q1A1B*E1B + q1A1C*E1C + q1A1D*E1D);
	ydot[5] = -(death_rate1B + (q1B0B + q1B1A + q1B1C + q1B1D) + birth_rate1B) * E1B + birth_rate1B * E1B * E1B + death_rate1B + (q1B0B*E0B + q1B1A*E1A + q1B1C*E1C + q1B1D*E1D);
	ydot[6] = -(death_rate1C + (q1C0C + q1C1A + q1C1B + q1C1D) + birth_rate1C) * E1C + birth_rate1C * E1C * E1C + death_rate1C + (q1C0C*E0C + q1C1A*E1A + q1C1B*E1B + q1C1D*E1D);
	ydot[7] = -(death_rate1D + (q1D0D + q1D1A + q1D1B + q1D1C) + birth_rate1D) * E1D + birth_rate1D * E1D * E1D + death_rate1D + (q1D0D*E0D + q1D1A*E1A + q1D1B*E1B + q1D1C*E1C);

	ydot[8]  = -(death_rate0A + (q0A0B + q0A0C + q0A0D + q0A1A) + birth_rate0A) * D0A + 2 * birth_rate0A * E0A * D0A + (q0A0B*D0B + q0A0C*D0C + q0A0D*D0D + q0A1A*D1A);
	ydot[9]  = -(death_rate0B + (q0B0A + q0B0C + q0B0D + q0B1B) + birth_rate0B) * D0B + 2 * birth_rate0B * E0B * D0B + (q0B0A*D0A + q0B0C*D0C + q0B0D*D0D + q0B1B*D1B);
	ydot[10] = -(death_rate0C + (q0C0A + q0C0B + q0C0D + q0C1C) + birth_rate0C) * D0C + 2 * birth_rate0C * E0C * D0C + (q0C0A*D0A + q0C0B*D0B + q0C0D*D0D + q0C1C*D1C);
	ydot[11] = -(death_rate0D + (q0D0A + q0D0B + q0D0C + q0D1D) + birth_rate0D) * D0D + 2 * birth_rate0D * E0D * D0D + (q0D0A*D0A + q0D0B*D0B + q0D0C*D0C + q0D1D*D1D);
	ydot[12] = -(death_rate1A + (q1A0A + q1A1B + q1A1C + q1A1D) + birth_rate1A) * D1A + 2 * birth_rate1A * E1A * D1A + (q1A0A*D0A + q1A1B*D1B + q1A1C*D1C + q1A1D*D1D);
	ydot[13] = -(death_rate1B + (q1B0B + q1B1A + q1B1C + q1B1D) + birth_rate1B) * D1B + 2 * birth_rate1B * E1B * D1B + (q1B0B*D0B + q1B1A*D1A + q1B1C*D1C + q1B1D*D1D);
	ydot[14] = -(death_rate1C + (q1C0C + q1C1A + q1C1B + q1C1D) + birth_rate1C) * D1C + 2 * birth_rate1C * E1C * D1C + (q1C0C*D0C + q1C1A*D1A + q1C1B*D1B + q1C1D*D1D);
	ydot[15] = -(death_rate1D + (q1D0D + q1D1A + q1D1B + q1D1C) + birth_rate1D) * D1D + 2 * birth_rate1D * E1D * D1D + (q1D0D*D0D + q1D1A*D1A + q1D1B*D1B + q1D1C*D1C);

}


