/*
 *  bisse-ext-derivs_c
 *  
 *
 *  Created by Jeremy Beaulieu 11/5/2014
 *  Copyright 2014 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 28

static double params_bisse[NUMELEMENTS];

void initmod_bisse(void (* odeparms)(int *, double *)){
	int N = NUMELEMENTS;
	odeparms(&N, params_bisse);
}


double set_transitions_bisse(double *point_time, double *timeslice, double *rate, double *slice_factor){
	if(*timeslice < *point_time){
		double factor_slice = *slice_factor;
		double new_rate = *rate * factor_slice;
		return(new_rate);
	}else{
		double new_rate = *rate;
		return(new_rate);
	}
}


double set_parameter_bisse(int turnover, int at_speciation, double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha, double *turnover_trend_beta, double *turnover_beta_factor, double *turnover_slice_factor, double *eps_trend_alpha, double *eps_trend_beta, double *eps_beta_factor, double *eps_slice_factor, double *x_turnover, double *x_eps, double *q01, double *q10, double *focal_edge_length, double *tipward_age){

	int lg=0;
	double param_mean = 0;
	if(turnover==1){
		if(*timeslice < *point_time){
			*x_turnover = *x_turnover * *turnover_slice_factor;
		}
		double x_a = *x_turnover;
		//printf("%f\nturnover", *x_turnover);
		double proportion_up_from_root = (*tot_time - *point_time ) / *tot_time;
		double rescaled_up_from_root = 0.1 + proportion_up_from_root * (0.9 - 0.1);
		double param_beta_scaling = *turnover_beta_factor * dbeta(rescaled_up_from_root, *turnover_trend_alpha, *turnover_trend_beta, lg);
		//printf("%f\nparam_beta_scaling", param_beta_scaling);
		param_mean = param_beta_scaling * x_a;
		//printf("%f\nparam_mean", param_mean);
	}else{
		//The following rescales the rate by the slice factor. As of now, we only allow a single timeslice:
		if(*timeslice < *point_time){
			*x_eps = *x_eps * *eps_slice_factor;
		}
		double x_a = *x_eps;
		double proportion_up_from_root = (*tot_time - *point_time ) / *tot_time;
		double rescaled_up_from_root = 0.1 + proportion_up_from_root * (0.9 - 0.1);
		double param_beta_scaling = *eps_beta_factor * dbeta(rescaled_up_from_root, *eps_trend_alpha, *eps_trend_beta, lg);
		param_mean = param_beta_scaling * x_a;
		//printf("%f\nextinction", param_mean);
	}
	return param_mean;
}

//If we are calling this function then we are at a node:
void set_birth_bisse_void(double *birth, double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1, double *q01, double *q10, double *focal_edge_length, double *tipward_age, int *state){

	int turnover = 1;
	int at_speciation = 1;
	if(*state == 0){
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		*birth = turnover_rate / (1 + extinction_fraction);
	}else{
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		*birth = turnover_rate / (1 + extinction_fraction);
	}
}


//If we are calling this function then we are somewhere along a branch:
double set_birth_bisse(double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1, double *q01, double *q10, double *focal_edge_length, double *tipward_age, int *state){

	int turnover = 1;
	int at_speciation = 0;
	//printf("%d\n", state);
	if(*state == 0){
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		//printf("%f\nturnover_rate", turnover_rate);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		double birth_rate = turnover_rate / (1 + extinction_fraction);
		//printf("%f\nbirthrate", birth_rate);
		return birth_rate;
	}else{
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		double birth_rate = turnover_rate / (1 + extinction_fraction);
		//printf("%f\nbirthrate", birth_rate);
		return birth_rate;
	}
}


//If we are calling this function then we are at a node:
void set_death_bisse_void(double *death, double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1, double *q01, double *q10, double *focal_edge_length, double *tipward_age, int state){

	int turnover = 1;
	int at_speciation = 1;
	if(state == 0){
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		*death = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
	}else{
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		*death = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
	}
}


//If we are calling this function then we are somewhere along a branch:
double set_death_bisse(double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1, double *q01, double *q10, double *focal_edge_length, double *tipward_age,int *state){
	
	int turnover = 1;
	int at_speciation = 0;
	if(*state == 0){
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, focal_edge_length, tipward_age);
		double death_rate = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
		//printf("%f\ndeathrate", death_rate);
		return death_rate;
	}else{
		double turnover_rate = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_bisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, focal_edge_length, tipward_age);
		double death_rate = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
		//printf("%f\ndeathrate", death_rate);
		return death_rate;
	}
}


void maddison_DE_bisse(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
	
	double E0 = y[0], E1 = y[1];
	double D0 = y[2], D1 = y[3];

	double point_time = *t;
	double tot_time = params_bisse[0],
	timeslice = params_bisse[1],

	//Why not include all this stuff?
	turnover_trend_alpha0 = params_bisse[2],
	turnover_trend_beta0 = params_bisse[3],
	turnover_beta_factor0 = params_bisse[4],
	turnover_slice_factor0 = params_bisse[5],
	eps_trend_alpha0 = params_bisse[6],
	eps_trend_beta0 = params_bisse[7],
	eps_beta_factor0 = params_bisse[8],
	eps_slice_factor0 = params_bisse[9],

	turnover_trend_alpha1 = params_bisse[10],
	turnover_trend_beta1 = params_bisse[11],
	turnover_beta_factor1 = params_bisse[12],
	turnover_slice_factor1 = params_bisse[13],
	eps_trend_alpha1 = params_bisse[14],
	eps_trend_beta1 = params_bisse[15],
	eps_beta_factor1 = params_bisse[16],
	eps_slice_factor1 = params_bisse[17],
	
	x_turnover0 = params_bisse[18],
	x_eps0 = params_bisse[19],
	x_turnover1 = params_bisse[20],
	x_eps1 = params_bisse[21],
	q01 = params_bisse[22],
	q10 = params_bisse[23],
	q01_slice_factor = params_bisse[24],
	q10_slice_factor = params_bisse[25],
	focal_edge_length = params_bisse[26],
	tipward_age = params_bisse[27];
	//State 0:
	int state = 0;
	double birth_rate0 = set_birth_bisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1, &q01, &q10, &focal_edge_length, &tipward_age, &state);
	//printf("%fbirthrate0 \n", birth_rate0);
	double death_rate0 = set_death_bisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1, &q01, &q10, &focal_edge_length, &tipward_age,&state);
	//printf("%fdeathrate0 \n", death_rate0);
	//State 1:
	state = 1;
	double birth_rate1 = set_birth_bisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1, &q01, &q10, &focal_edge_length, &tipward_age, &state);
	//printf("%fbirthrate1 \n", birth_rate1);
	double death_rate1 = set_death_bisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1, &q01, &q10, &focal_edge_length, &tipward_age,&state);
	//printf("%fdeathrate1 \n", death_rate1);

	double nq01 = set_transitions_bisse(&point_time, &timeslice, &q01, &q01_slice_factor);
	double nq10 = set_transitions_bisse(&point_time, &timeslice, &q10, &q10_slice_factor);
	
	ydot[0] = -(death_rate0 + nq01 + birth_rate0) * E0 + birth_rate0 * E0 * E0 + death_rate0 + nq01 * E1;
	ydot[1] = -(death_rate1 + nq10 + birth_rate1) * E1 + birth_rate1 * E1 * E1 + death_rate1 + nq10 * E0;
	ydot[2] = -(death_rate0 + nq01 + birth_rate0) * D0 + 2 * birth_rate0 * E0 * D0 + nq01 * D1;
	ydot[3] = -(death_rate1 + nq10 + birth_rate1) * D1 + 2 * birth_rate1 * E1 * D1 + nq10 * D0;
}


