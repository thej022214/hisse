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

  // Here we will use a letter notation. endemic areas are: 1, 2, 3. Widespread areas are 12, 23, 13.
  // The three range widespread are 123 is not considered here to facilitate model implementation.
  // All parameters will be associated with these nomenclature.

  // Need to transform the nomenclature to numeric.
  
  double
    E_1 = y[0],      /* ODE for extinction */
    E_2 = y[1],
    E_3 = y[2],
    E_12 = y[3],
    E_23 = y[4],
    E_13 = y[5],            
    
    D_N1 = y[6],     /* ODE for other processes */
    D_N2 = y[7],
    D_N3 = y[8],
    D_N12 = y[9],
    D_N23 = y[10],
    D_N13 = y[11],            
    
    s1  = params_geosse_3_areas[0],     /* speciation within region 1 */
    s2  = params_geosse_3_areas[1],     /* speciation within region 2 */
    s3  = params_geosse_3_areas[2],     /* speciation within region 3 */
    s12 = params_geosse_3_areas[3],     /* between-region speciation  */
    s23 = params_geosse_3_areas[4],     /* between-region speciation  */
    s13 = params_geosse_3_areas[5],     /* between-region speciation  */        

    // These are true extinction of endemic lineages.
    x1  = params_geosse_3_areas[6],     /* extinction from region 1   */
    x2  = params_geosse_3_areas[7],     /* extinction from region 2   */
    x3  = params_geosse_3_areas[8],     /* extinction from region 3   */

    // These are transitions between endemic areas. These parameters are only used in the +jump models.
    // The original GeoSSE process assumes that these jumps have 0 rate.
    d1_2 = params_geosse_3_areas[9],   /* jumps from 1 to 2          */
    d2_1 = params_geosse_3_areas[10],   /* jumps from 2 to 1          */
    d2_3 = params_geosse_3_areas[11],   /* jumps from 2 to 3          */
    d3_2 = params_geosse_3_areas[12],   /* jumps from 3 to 2          */
    d1_3 = params_geosse_3_areas[13],   /* jumps from 1 to 3          */
    d3_1 = params_geosse_3_areas[14],   /* jumps from 3 to 1          */

    // These are dispersal parameters:
    d1_12 = params_geosse_3_areas[15],  /* dispersal from 1 to 12     */
    d1_13 = params_geosse_3_areas[16],  /* dispersal from 1 to 13     */
    d2_12 = params_geosse_3_areas[17],  /* dispersal from 2 to 12     */
    d2_23 = params_geosse_3_areas[18],  /* dispersal from 2 to 23     */        
    d3_13 = params_geosse_3_areas[19],  /* dispersal from 3 to 13     */
    d3_23 = params_geosse_3_areas[20],  /* dispersal from 3 to 23     */

    // These are the local extinction or extirpation parameters.
    d12_1 = params_geosse_3_areas[21],   /* true extirpation rate  */
    d13_1 = params_geosse_3_areas[22],   /* true extirpation rate  */
    d12_2 = params_geosse_3_areas[23],   /* true extirpation rate  */
    d23_2 = params_geosse_3_areas[24],   /* true extirpation rate  */
    d13_3 = params_geosse_3_areas[25],   /* true extirpation rate  */
    d23_3 = params_geosse_3_areas[26];   /* true extirpation rate  */
  
  // The extinction ODEs
  
  /*  dE_1 / dt  */
  ydot[0] = -(s1 + d1_2 + d1_3 + d1_12 + d1_13 + x1) * E_1 + x1 + d1_2 * E_2 + d1_3 * E_3 + d1_12 * E_12 + d1_13 * E_13 + s1 * E_1 * E_1;

  /*  dE_2 / dt  */
  ydot[1] = -(s2 + d2_1 + d2_3 + d2_12 + d2_23 + x2) * E_2 + x2 + d2_1 * E_1 + d2_3 * E_3 + d2_12 * E_12 + d2_23 * E_23 + s2 * E_2 * E_2;

  /*  dE_3 / dt  */
  ydot[2] = -(s3 + d3_1 + d3_2 + d3_13 + d3_23 + x3) * E_3 + x3 + d3_1 * E_1 + d3_2 * E_2 + d3_13 * E_13 + d3_23 * E_23 + (s3 * E_3 * E_3);
   
  /*  dE_12 / dt  */
  ydot[3] = -(s1 + s2 + s12 + d12_1 + d12_2) * E_12 + d12_1 * E_1 + d12_2 * E_2 + s1 * E_1 * E_12 + s2 * E_2 * E_12 + s12 * E_1 * E_2;

  /*  dE_23 / dt  */
  ydot[4] = -(s2 + s3 + s23 + d23_2 + d23_3) * E_23 + d23_2 * E_2 + d23_3 * E_3 + s2 * E_2 * E_23 + s3 * E_3 * E_23 + s23 * E_2 * E_3;

  /*  dE_13 / dt  */
  ydot[5] = -(s1 + s3 + s13 + d13_1 + d13_3) * E_13 + d13_1 * E_1 + d13_3 * E_3 + s1 * E_1 * E_13 + s3 * E_3 * E_13 + s13 * E_1 * E_3;

  // The events ODEs - dispersion and speciation
  
  /*  dD_N1 / dt  */
  /* ydot[6] = -(s1 + d1_2 + d1_3 + d1_12 + d1_13 + x1) * D_N1 + (d1_3 * D_N3 + d1_2 * D_N2 + d1_13 * D_N13 + d1_12 * D_N12) + s1 * (D_N1 * E_1 + D_N1 * E_1); */
  ydot[6] = -(s1 + d1_2 + d1_3 + d1_12 + d1_13 + x1) * D_N1 + d1_3 * D_N3 + d1_2 * D_N2 + d1_13 * D_N13 + d1_12 * D_N12 + 2.0 * s1 * D_N1 * E_1;

  /*  dD_N2 / dt  */
  ydot[7] = -(s2 + d2_1 + d2_3 + d2_12 + d2_23 + x2) * D_N2 + d2_3 * D_N3 + d2_1 * D_N1 + d2_23 * D_N23 + d2_12 * D_N12 + 2.0 * s2 * D_N2 * E_2;

  /*  dD_N3 / dt  */
  ydot[8] = -(s3 + d3_1 + d3_2 + d3_23 + d3_13 + x3) * D_N3 + d3_1 * D_N1 + d3_2 * D_N2 + d3_23 * D_N23 + d3_13 * D_N13 + 2.0 * s3 * D_N3 * E_3;
    
  /*  dD_N12 / dt  */
  ydot[9] = -(s1 + s2 + s12 + d12_1 + d12_2) * D_N12 + d12_1 * D_N1 + d12_2 * D_N2 + s12 * (D_N1 * E_2 + D_N2 * E_1) + s1 * (E_1 * D_N12 + E_12 * D_N1) + s2 * (E_2 * D_N12 + E_12 * D_N2);

  /*  dD_N23 / dt  */
  ydot[10] = -(s2 + s3 + s23 + d23_2 + d23_3) * D_N23 + d23_2 * D_N2 + d23_3 * D_N3 + s23 * (D_N2 * E_3 + D_N3 * E_2) + s2 * (E_2 * D_N23 + E_23 * D_N2) + s3 * (E_3 * D_N23 + E_23 * D_N3);

  /*  dD_N13 / dt  */
  ydot[11] = -(s1 + s3 + s13 + d13_1 + d13_3) * D_N13 + d13_1 * D_N1 + d13_3 * D_N3 + s13 * (D_N1 * E_3 + D_N3 * E_1) + s1 * (E_1 * D_N13 + E_13 * D_N1) + s3 * (E_3 * D_N13 + E_13 * D_N3);
  
}
