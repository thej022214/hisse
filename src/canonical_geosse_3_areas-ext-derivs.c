/*
 *  canonical_geosse_3-areas-ext-derivs_c
 *
 *  Modified by Daniel Caetano from Jeremy Beaulieu 6/06/2017
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 27
// This is the number of total parameters in the model.

static double params_geosse_3_areas[NUMELEMENTS];

void initmod_geosse_3_areas(void (* odeparms)(int *, double *)){
  int N = NUMELEMENTS;
  odeparms(&N, params_geosse_3_areas);
}

void geosse_3_areas_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){

  // Here we will use a letter notation. endemic areas are: A, B, C. Widespread areas are AB, BC, AC.
  // The three range widespread are ABC is not considered here to facilitate model implementation.
  // All parameters will be associated with these nomenclature.
  
  double
    E_A = y[0],      /* ODE for extinction */
    E_B = y[1],
    E_C = y[2],
    E_AB = y[3],
    E_BC = y[4],
    E_AC = y[5],            
    
    D_NA = y[6],     /* ODE for other processes */
    D_NB = y[7],
    D_NC = y[8],
    D_NAB = y[9],
    D_NBC = y[10],
    D_NAC = y[11],            
    
    sA  = params_geosse_3_areas[0],     /* speciation within region A */
    sB  = params_geosse_3_areas[1],     /* speciation within region B */
    sC  = params_geosse_3_areas[2],     /* speciation within region C */
    sAB = params_geosse_3_areas[3],     /* between-region speciation  */
    sBC = params_geosse_3_areas[4],     /* between-region speciation  */
    sAC = params_geosse_3_areas[5],     /* between-region speciation  */        

    // These are true extinction of endemic lineages.
    xA  = params_geosse_3_areas[6],     /* extinction from region A   */
    xB  = params_geosse_3_areas[7],     /* extinction from region B   */
    xC  = params_geosse_3_areas[8],     /* extinction from region C   */

    // These are transitions between endemic areas. These parameters are only used in the +jump models.
    // The original GeoSSE process assumes that these jumps have 0 rate.
    dA_B = params_geosse_3_areas[9],   /* jumps from A to B          */
    dB_A = params_geosse_3_areas[10],   /* jumps from B to A          */
    dB_C = params_geosse_3_areas[11],   /* jumps from B to C          */
    dC_B = params_geosse_3_areas[12],   /* jumps from C to B          */
    dA_C = params_geosse_3_areas[13],   /* jumps from A to C          */
    dC_A = params_geosse_3_areas[14],   /* jumps from C to A          */

    // These are dispersal parameters:
    dA_AB = params_geosse_3_areas[15],  /* dispersal from A to AB     */
    dA_AC = params_geosse_3_areas[16],  /* dispersal from A to AC     */
    dB_AB = params_geosse_3_areas[17],  /* dispersal from B to AB     */
    dB_BC = params_geosse_3_areas[18],  /* dispersal from B to BC     */        
    dC_AC = params_geosse_3_areas[19],  /* dispersal from C to AC     */
    dC_BC = params_geosse_3_areas[20],  /* dispersal from C to BC     */

    // These are the local extinction or extirpation parameters.
    dAB_A = params_geosse_3_areas[21],   /* true extirpation rate  */
    dAC_A = params_geosse_3_areas[22],   /* true extirpation rate  */
    dAB_B = params_geosse_3_areas[23],   /* true extirpation rate  */
    dBC_B = params_geosse_3_areas[24],   /* true extirpation rate  */
    dAC_C = params_geosse_3_areas[25],   /* true extirpation rate  */
    dBC_C = params_geosse_3_areas[26];   /* true extirpation rate  */
  // The extinction ODEs
  
  /*  dE_A / dt  */
  ydot[0] = -(sA + dA_B + dA_C + dA_AB + dA_AC + xA) * E_A + xA + (dA_B * E_B + dA_C * E_C + dA_AB * E_AB + dA_AC * E_AC) + (sA * E_A * E_A);

  /*  dE_B / dt  */
  ydot[1] = -(sB + dB_A + dB_C + dB_AB + dB_BC + xB) * E_B + xB + (dB_A * E_A + dB_C * E_C + dB_AB * E_AB + dB_BC * E_BC) + (sB * E_B * E_B);

  /*  dE_C / dt  */
  ydot[2] = -(sC + dC_A + dC_B + dC_AC + dC_BC + xC) * E_C + xC + (dC_A * E_A + dC_B * E_B + dC_AC * E_AC + dC_BC * E_BC) + (sC * E_C * E_C);
   
  /*  dE_AB / dt  */
  ydot[3] = -(sA + sB + sAB + dAB_A + dAB_B) * E_AB + (dAB_A * E_A + dAB_B * E_B) + sA * E_A * E_AB + sB * E_B * E_AB + sAB * E_A * E_B;

  /*  dE_BC / dt  */
  ydot[4] = -(sB + sC + sBC + dBC_B + dBC_C) * E_BC + (dBC_B * E_B + dBC_C * E_C) + sB * E_B * E_BC + sC * E_C * E_BC + sBC * E_B * E_C;

  /*  dE_AC / dt  */
  ydot[5] = -(sA + sC + sAC + dAC_A + dAC_C) * E_AC + (dAC_A * E_A + dAC_C * E_C) + sA * E_A * E_AC + sC * E_C * E_AC + sAC * E_A * E_C;

  // The events ODEs - dispersion and speciation
  
  /*  dD_NA / dt  */
  ydot[6] = -(sA + dA_B + dA_C + dA_AB + dA_AC + xA) * D_NA + (dA_C * D_NC + dA_B * D_NB + dA_AC * D_NAC + dA_AB * D_NAB) + sA * (D_NA * E_A + D_NA * E_A);

  /*  dD_NB / dt  */
  ydot[7] = -(sB + dB_A + dB_C + dB_AB + dB_BC + xB) * D_NB + (dB_C * D_NC + dB_A * D_NA + dB_BC * D_NBC + dB_BC * D_NBC) + sB * (D_NB * E_B + D_NB * E_B);

  /*  dD_NC / dt  */
  ydot[8] = -(sC + dC_A + dC_B + dC_BC + dC_AC + xC) * D_NC + (dC_A * D_NA + dC_B * D_NB + dC_BC * D_NBC + dC_AC * D_NAC) + sC * (D_NC * E_C + D_NC * E_C);
    
  /*  dD_NAB / dt  */
  ydot[9] = -(sA + sB + sAB + dAB_A + dAB_B) * D_NAB + (dAB_A * D_NA + dAB_B * D_NB) + sAB * (D_NA * E_B + D_NB * E_A) + sA * (E_A * D_NAB + E_AB * D_NA) + sB * (E_B * D_NAB + E_AB * D_NB);

  /*  dD_NBC / dt  */
  ydot[10] = -(sB + sC + sBC + dBC_B + dBC_C) * D_NBC + (dBC_B * D_NB + dBC_C * D_NC) + sBC * (D_NB * E_C + D_NC * E_B) + sB * (E_B * D_NBC + E_BC * D_NB) + sC * (E_C * D_NBC + E_BC * D_NC);

  /*  dD_NAC / dt  */
  ydot[11] = -(sA + sC + sAC + dAC_A + dAC_C) * D_NAC + (dAC_A * D_NA + dAC_C * D_NC) + sAC * (D_NA * E_C + D_NC * E_A) + sA * (E_A * D_NAC + E_AC * D_NA) + sC * (E_C * D_NAC + E_AC * D_NC);
    
}
