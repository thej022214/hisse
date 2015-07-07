/*
 *  hisse-ext-derivs_c
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
#define NUMELEMENTS 68

static double params_hisse[NUMELEMENTS]; 


void initmod_hisse(void (* odeparms)(int *, double *)){
	int N = NUMELEMENTS;
	odeparms(&N, params_hisse);
}


double set_transitions_hisse(double *point_time, double *timeslice, double *rate, double *slice_factor){
	if(*timeslice < *point_time){
		double new_rate = *rate * *slice_factor;
		return(new_rate);
	}else{
		return(*rate);
	}
}
							

void set_parameter_hisse_sim(double *param, int *turnover, int *n_taxa, double *tot_time, double *point_time, double *turnover_param_anc, double *turnover_param_indep, double *turnover_sigma_indep, double *turnover_scaling, double *turnover_weight_logistic, double *turnover_trend_alpha, double *turnover_trend_beta, double *turnover_scale_factor, double *turnover_sigma_anc, double *eps_param_anc, double *eps_param_indep, double *eps_sigma_indep, double *eps_scaling, double *eps_weight_logistic, double *eps_trend_alpha, double *eps_trend_beta, double *eps_scale_factor, double *eps_sigma_anc, double *turn_k, double *eps_k, double *turnover_kick_rate, double *turnover_kick_stdev, double *eps_kick_rate, double *eps_kick_stdev, double *split_times, double *turnover_splits, double *eps_splits, double *x_at_turnover, double *x_at_eps, int *length_split_times){
	//For some reason if we are at the present this returns n_taxa = 0. In these cases just fill it to be equal to the total diversity of the tree:
	if(*n_taxa == 0){
		*n_taxa = *length_split_times + 1;
	}
	int lg=0;
	//printf("%d\n", *turnover);
	if(*turnover==1){
		//The following line may be used eventually for the time slice model. But keep an eye on position -- any discrepancies are likely to be due to this:
		//double param_indep = turnover_splits[position];
		//
		double x_a = *x_at_turnover * *turnover_scaling;
		double proportion_up_from_root = (*tot_time - *point_time ) / *tot_time;
		double rescaled_up_from_root = 0.1 + proportion_up_from_root * (0.9 - 0.1);
		double param_scaling = *turnover_scale_factor * dbeta(rescaled_up_from_root, *turnover_trend_alpha, *turnover_trend_beta, lg);
		//double logistic_scaling = 1.0 - (*turnover_weight_logistic * (*n_taxa / *turn_k));
		*param = param_scaling * x_a;
	}else{
		//The following line may be used eventually for the time slice model. But keep an eye on position -- any discrepancies are likely to be due to this:
		//double param_indep = eps_splits[position];
		//
		double x_a = *x_at_eps * *eps_scaling;
		double proportion_up_from_root = (*tot_time - *point_time ) / *tot_time;
		double rescaled_up_from_root = 0.1 + proportion_up_from_root * (0.9 - 0.1);
		double param_scaling = *eps_scale_factor * dbeta(rescaled_up_from_root, *eps_trend_alpha, *eps_trend_beta, lg);
		//double logistic_scaling = 1.0 - (*eps_weight_logistic * (*n_taxa / *eps_k));
		*param = param_scaling * x_a;
	}
}


double set_parameter_hisse(int turnover, int at_speciation, double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha, double *turnover_trend_beta, double *turnover_beta_factor, double *turnover_slice_factor, double *eps_trend_alpha, double *eps_trend_beta, double *eps_beta_factor, double *eps_slice_factor, double *x_turnover, double *x_eps, double *q01, double *q10, double *q0A, double *qA0, double *q1B, double *qB1, double *q0B, double *qB0, double *q1A, double *qA1, double *qBA, double *qAB, double *focal_edge_length, double *tipward_age){

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
void set_birth_hisse_void(double *birth, double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *turnover_trend_alphaA, double *turnover_trend_betaA, double *turnover_beta_factorA, double *turnover_slice_factorA, double *eps_trend_alphaA, double *eps_trend_betaA, double *eps_beta_factorA, double *eps_slice_factorA, double *turnover_trend_alphaB, double *turnover_trend_betaB, double *turnover_beta_factorB, double *turnover_slice_factorB, double *eps_trend_alphaB, double *eps_trend_betaB, double *eps_beta_factorB, double *eps_slice_factorB, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1,double *x_turnoverA, double *x_epsA, double *x_turnoverB, double *x_epsB, double *q01, double *q10, double *q0A, double *qA0, double *q1B, double *qB1, double *q0B, double *qB0, double *q1A, double *qA1, double *qBA, double *qAB, double *focal_edge_length, double *tipward_age, int *state){

	int turnover = 1;
	int at_speciation = 1;
	if(*state == 0){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*birth = turnover_rate / (1 + extinction_fraction);
	}
	if(*state == 1){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*birth = turnover_rate / (1 + extinction_fraction);
	}
	if(*state == 2){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*birth = turnover_rate / (1 + extinction_fraction);
	}
	if(*state == 3){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*birth = turnover_rate / (1 + extinction_fraction);
	}
	
}


//If we are calling this function then we are somewhere along a branch:
double set_birth_hisse(double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *turnover_trend_alphaA, double *turnover_trend_betaA, double *turnover_beta_factorA, double *turnover_slice_factorA, double *eps_trend_alphaA, double *eps_trend_betaA, double *eps_beta_factorA, double *eps_slice_factorA, double *turnover_trend_alphaB, double *turnover_trend_betaB, double *turnover_beta_factorB, double *turnover_slice_factorB, double *eps_trend_alphaB, double *eps_trend_betaB, double *eps_beta_factorB, double *eps_slice_factorB, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1,double *x_turnoverA, double *x_epsA, double *x_turnoverB, double *x_epsB, double *q01, double *q10, double *q0A, double *qA0, double *q1B, double *qB1, double *q0B, double *qB0, double *q1A, double *qA1, double *qBA, double *qAB, double *focal_edge_length, double *tipward_age, int *state){

	int turnover = 1;
	int at_speciation = 0;
	//printf("%d\n", state);
	if(*state == 0){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double birth_rate = turnover_rate / (1 + extinction_fraction);
		//printf("%f\nbirthrate", birth_rate);
		return birth_rate;
	}
	if(*state == 1){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double birth_rate = turnover_rate / (1 + extinction_fraction);
		return birth_rate;
	}
	if(*state == 2){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double birth_rate = turnover_rate / (1 + extinction_fraction);
		return birth_rate;
	}
	if(*state == 3){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double birth_rate = turnover_rate / (1 + extinction_fraction);
		return birth_rate;
	}
	return 0;
}


//If we are calling this function then we are at a node:
void set_death_hisse_void(double *death, double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *turnover_trend_alphaA, double *turnover_trend_betaA, double *turnover_beta_factorA, double *turnover_slice_factorA, double *eps_trend_alphaA, double *eps_trend_betaA, double *eps_beta_factorA, double *eps_slice_factorA, double *turnover_trend_alphaB, double *turnover_trend_betaB, double *turnover_beta_factorB, double *turnover_slice_factorB, double *eps_trend_alphaB, double *eps_trend_betaB, double *eps_beta_factorB, double *eps_slice_factorB, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1,double *x_turnoverA, double *x_epsA, double *x_turnoverB, double *x_epsB, double *q01, double *q10, double *q0A, double *qA0, double *q1B, double *qB1, double *q0B, double *qB0, double *q1A, double *qA1, double *qBA, double *qAB, double *focal_edge_length, double *tipward_age, int state){

	int turnover = 1;
	int at_speciation = 1;
	if(state == 0){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB,focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*death = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
	}
	if(state == 1){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*death = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
	}
	if(state == 2){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*death = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);		
	}
	if(state == 3){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		*death = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
	}
	
}


//If we are calling this function then we are somewhere along a branch:
double set_death_hisse(double *point_time, double *tot_time, double *timeslice, double *turnover_trend_alpha0, double *turnover_trend_beta0, double *turnover_beta_factor0, double *turnover_slice_factor0, double *eps_trend_alpha0, double *eps_trend_beta0, double *eps_beta_factor0, double *eps_slice_factor0, double *turnover_trend_alpha1, double *turnover_trend_beta1, double *turnover_beta_factor1, double *turnover_slice_factor1, double *eps_trend_alpha1, double *eps_trend_beta1, double *eps_beta_factor1, double *eps_slice_factor1, double *turnover_trend_alphaA, double *turnover_trend_betaA, double *turnover_beta_factorA, double *turnover_slice_factorA, double *eps_trend_alphaA, double *eps_trend_betaA, double *eps_beta_factorA, double *eps_slice_factorA, double *turnover_trend_alphaB, double *turnover_trend_betaB, double *turnover_beta_factorB, double *turnover_slice_factorB, double *eps_trend_alphaB, double *eps_trend_betaB, double *eps_beta_factorB, double *eps_slice_factorB, double *x_turnover0, double *x_eps0, double *x_turnover1, double *x_eps1,double *x_turnoverA, double *x_epsA, double *x_turnoverB, double *x_epsB, double *q01, double *q10, double *q0A, double *qA0, double *q1B, double *qB1, double *q0B, double *qB0, double *q1A, double *qA1, double *qBA, double *qAB, double *focal_edge_length, double *tipward_age, int *state){
	
	int turnover = 1;
	int at_speciation = 0;
	if(*state == 0){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha0, turnover_trend_beta0, turnover_beta_factor0, turnover_slice_factor0, eps_trend_alpha0, eps_trend_beta0, eps_beta_factor0, eps_slice_factor0, x_turnover0, x_eps0, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double death_rate = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
		return death_rate;
	}
	if(*state == 1){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice,  turnover_trend_alpha1, turnover_trend_beta1, turnover_beta_factor1, turnover_slice_factor1, eps_trend_alpha1, eps_trend_beta1, eps_beta_factor1, eps_slice_factor1, x_turnover1, x_eps1, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double death_rate = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
		return death_rate;
	}
	if(*state == 2){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaA, turnover_trend_betaA, turnover_beta_factorA, turnover_slice_factorA, eps_trend_alphaA, eps_trend_betaA, eps_beta_factorA, eps_slice_factorA, x_turnoverA, x_epsA, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double death_rate = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
		return death_rate;
	}
	if(*state == 3){
		double turnover_rate = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		turnover = 0;
		double extinction_fraction = set_parameter_hisse(turnover, at_speciation, point_time, tot_time, timeslice, turnover_trend_alphaB, turnover_trend_betaB, turnover_beta_factorB, turnover_slice_factorB, eps_trend_alphaB, eps_trend_betaB, eps_beta_factorB, eps_slice_factorB, x_turnoverB, x_epsB, q01, q10, q0A, qA0, q1B, qB1, q0B, qB0, q1A, qA1, qBA, qAB, focal_edge_length, tipward_age);
		double death_rate = (extinction_fraction * turnover_rate) / (1 + extinction_fraction);
		return death_rate;
	}
	return 0;
}


void maddison_DE_hisse(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
	
	double E0 = y[0], E1 = y[1], EA = y[2], EB = y[3];
	double D0 = y[4], D1 = y[5], DA = y[6], DB = y[7];

	double point_time = *t;
	double tot_time = params_hisse[0],
	timeslice = params_hisse[1],

	turnover_trend_alpha0 = params_hisse[2],
	turnover_trend_beta0 = params_hisse[3],
	turnover_beta_factor0 = params_hisse[4],
	turnover_slice_factor0 = params_hisse[5],
	eps_trend_alpha0 = params_hisse[6],
	eps_trend_beta0 = params_hisse[7],
	eps_beta_factor0 = params_hisse[8],
	eps_slice_factor0 = params_hisse[9],

	turnover_trend_alpha1 = params_hisse[10],
	turnover_trend_beta1 = params_hisse[11],
	turnover_beta_factor1 = params_hisse[12],
	turnover_slice_factor1 = params_hisse[13],
	eps_trend_alpha1 = params_hisse[14],
	eps_trend_beta1 = params_hisse[15],
	eps_beta_factor1 = params_hisse[16],
	eps_slice_factor1 = params_hisse[17],

	turnover_trend_alphaA = params_hisse[18],
	turnover_trend_betaA = params_hisse[19],
	turnover_beta_factorA = params_hisse[20],
	turnover_slice_factorA = params_hisse[21],
	eps_trend_alphaA = params_hisse[22],
	eps_trend_betaA = params_hisse[23],
	eps_beta_factorA = params_hisse[24],
	eps_slice_factorA = params_hisse[25],
	
	turnover_trend_alphaB = params_hisse[26],
	turnover_trend_betaB = params_hisse[27],
	turnover_beta_factorB = params_hisse[28],
	turnover_slice_factorB = params_hisse[29],
	eps_trend_alphaB = params_hisse[30],
	eps_trend_betaB = params_hisse[31],
	eps_beta_factorB = params_hisse[32],
	eps_slice_factorB = params_hisse[33],
		
	x_turnover0 = params_hisse[34],
	x_eps0 = params_hisse[35],
	x_turnover1 = params_hisse[36],
	x_eps1 = params_hisse[37],

	x_turnoverA = params_hisse[38],
	x_epsA = params_hisse[39],
	x_turnoverB = params_hisse[40],
	x_epsB = params_hisse[41],
		
	q01 = params_hisse[42],
	q10 = params_hisse[43],
	q0A = params_hisse[44], 
	qA0 = params_hisse[45],
	q1B = params_hisse[46],
	qB1 = params_hisse[47],
	q0B = params_hisse[48],
	qB0 = params_hisse[49],
	q1A = params_hisse[50],
	qA1 = params_hisse[51],
	qBA = params_hisse[52],
	qAB = params_hisse[53],

	q01_slice_factor = params_hisse[54],
	q10_slice_factor = params_hisse[55],
	q0A_slice_factor = params_hisse[56],
	qA0_slice_factor = params_hisse[57],
	q1B_slice_factor = params_hisse[58],
	qB1_slice_factor = params_hisse[59],
	q0B_slice_factor = params_hisse[60],
	qB0_slice_factor = params_hisse[61],
	q1A_slice_factor = params_hisse[62],
	qA1_slice_factor = params_hisse[63],
	qBA_slice_factor = params_hisse[64],
	qAB_slice_factor = params_hisse[65],

	focal_edge_length = params_hisse[66],
	tipward_age = params_hisse[67];
	//State 0:
	int state = 0;
	double birth_rate0 = set_birth_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	//printf("%fbirthrate0 \n", birth_rate0);
	double death_rate0 = set_death_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	//printf("%fdeathrate0 \n", death_rate0);
	//State 1:
	state = 1;
	double birth_rate1 = set_birth_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	//printf("%fbirthrate1 \n", birth_rate1);
	double death_rate1 = set_death_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	//printf("%fdeathrate1 \n", death_rate1);
	
	state = 2;
	double birth_rateA = set_birth_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	//printf("%fbirthrateA \n", birth_rateA);
	double death_rateA = set_death_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	
	state = 3;
	double birth_rateB = set_birth_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	//printf("%fbirthrateB \n", birth_rateB);
	double death_rateB = set_death_hisse(&point_time, &tot_time, &timeslice, &turnover_trend_alpha0, &turnover_trend_beta0, &turnover_beta_factor0, &turnover_slice_factor0, &eps_trend_alpha0, &eps_trend_beta0, &eps_beta_factor0, &eps_slice_factor0, &turnover_trend_alpha1, &turnover_trend_beta1, &turnover_beta_factor1, &turnover_slice_factor1, &eps_trend_alpha1, &eps_trend_beta1, &eps_beta_factor1, &eps_slice_factor1, &turnover_trend_alphaA, &turnover_trend_betaA, &turnover_beta_factorA, &turnover_slice_factorA, &eps_trend_alphaA, &eps_trend_betaA, &eps_beta_factorA, &eps_slice_factorA, &turnover_trend_alphaB, &turnover_trend_betaB, &turnover_beta_factorB, &turnover_slice_factorB, &eps_trend_alphaB, &eps_trend_betaB, &eps_beta_factorB, &eps_slice_factorB, &x_turnover0, &x_eps0, &x_turnover1, &x_eps1,&x_turnoverA, &x_epsA, &x_turnoverB, &x_epsB, &q01, &q10, &q0A, &qA0, &q1B, &qB1, &q0B, &qB0, &q1A, &qA1, &qBA, &qAB, &focal_edge_length, &tipward_age, &state);
	
	double nq01 = set_transitions_hisse(&point_time, &timeslice, &q01, &q01_slice_factor);
	double nq10 = set_transitions_hisse(&point_time, &timeslice, &q10, &q10_slice_factor);
	double nq0A = set_transitions_hisse(&point_time, &timeslice, &q0A, &q0A_slice_factor);
	double nqA0 = set_transitions_hisse(&point_time, &timeslice, &qA0, &qA0_slice_factor);
	double nq1B = set_transitions_hisse(&point_time, &timeslice, &q1B, &q1B_slice_factor);
	double nqB1 = set_transitions_hisse(&point_time, &timeslice, &qB1, &qB1_slice_factor);
	double nq0B = set_transitions_hisse(&point_time, &timeslice, &q0B, &q0B_slice_factor);
	double nqB0 = set_transitions_hisse(&point_time, &timeslice, &qB0, &qB0_slice_factor);
	double nq1A = set_transitions_hisse(&point_time, &timeslice, &q1A, &q1A_slice_factor);
	double nqA1 = set_transitions_hisse(&point_time, &timeslice, &qA1, &qA1_slice_factor);
	double nqBA = set_transitions_hisse(&point_time, &timeslice, &qBA, &qBA_slice_factor);
	double nqAB = set_transitions_hisse(&point_time, &timeslice, &qAB, &qAB_slice_factor);
	
	//Yes, yes, I wrote this out instead of being cool and making a sweet function to deal with the summations. Sue me. (plus I don't know C very well...)
	ydot[0] = -(death_rate0 + (nq01 + nq0A + nq0B) + birth_rate0) * E0 + birth_rate0 * E0 * E0 + death_rate0 + (nq01*E1 + nq0A*EA + nq0B*EB);
	ydot[1] = -(death_rate1 + (nq10 + nq1A + nq1B) + birth_rate1) * E1 + birth_rate1 * E1 * E1 + death_rate1 + (nq10*E0 + nq1A*EA + nq1B*EB);
	ydot[2] = -(death_rateA + (nqA0 + nqA1 + nqAB) + birth_rateA) * EA + birth_rateA * EA * EA + death_rateA + (nqA0*E0 + nqA1*E1 + nqAB*EB);
	ydot[3] = -(death_rateB + (nqB0 + nqB1 + nqBA) + birth_rateB) * EB + birth_rateB * EB * EB + death_rateB + (nqB0*E0 + nqB1*E1 + nqBA*EA);

	ydot[4] = -(death_rate0 + (nq01 + nq0A + nq0B) + birth_rate0) * D0 + 2 * birth_rate0 * E0 * D0 + (nq01*D1 + nq0A*DA + nq0B*DB);
	ydot[5] = -(death_rate1 + (nq10 + nq1A + nq1B) + birth_rate1) * D1 + 2 * birth_rate1 * E1 * D1 + (nq10*D0 + nq1A*DA + nq1B*DB);
	ydot[6] = -(death_rateA + (nqA0 + nqA1 + nqAB) + birth_rateA) * DA + 2 * birth_rateA * EA * DA + (nqA0*D0 + nqA1*D1 + nqAB*DB);
	ydot[7] = -(death_rateB + (nqB0 + nqB1 + nqBA) + birth_rateB) * DB + 2 * birth_rateB * EB * DB + (nqB0*D0 + nqB1*D1 + nqBA*DA);

}


