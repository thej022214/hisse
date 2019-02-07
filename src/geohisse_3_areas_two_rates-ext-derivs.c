/*
 *  geosse_3_areas_two_rates-ext-derivs.c
 *
 *  By Daniel Caetano Feb/2019
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 66
// This is the number of total parameters in the model.

static double params_3_areas_two_rates[NUMELEMENTS];

void initmod_geosse_3_areas_two_rates(void (* odeparms)(int *, double *)){
  int N = NUMELEMENTS;
  odeparms(&N, params_3_areas_two_rates);
}

void geosse_3_areas_two_rates_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){

  // Here we will use a letter notation. endemic areas are: 1, 2, 3. Widespread areas are 12, 23, 13.
  // The three range widespread are 123 is not considered here to facilitate model implementation.
  // All parameters will be associated with these nomenclature.

  // Need to transform the nomenclature to numeric.
  
  double
    E_1A = y[0],      /* ODE for extinction */
    E_2A = y[1],
    E_3A = y[2],
    E_12A = y[3],
    E_23A = y[4],
    E_13A = y[5],

    E_1B = y[6],
    E_2B = y[7],
    E_3B = y[8],
    E_12B = y[9],
    E_23B = y[10],
    E_13B = y[11],            
    
    D_N1A = y[12],     /* ODE for other processes */
    D_N2A = y[13],
    D_N3A = y[14],
    D_N12A = y[15],
    D_N23A = y[16],
    D_N13A = y[17],
    
    D_N1B = y[18],
    D_N2B = y[19],
    D_N3B = y[20],
    D_N12B = y[21],
    D_N23B = y[22],
    D_N13B = y[23],            
    
    s1A  = params_3_areas_two_rates[0],     /* speciation within region 1 */
    s2A  = params_3_areas_two_rates[1],     /* speciation within region 2 */
    s3A  = params_3_areas_two_rates[2],     /* speciation within region 3 */
    s12A = params_3_areas_two_rates[3],     /* between-region speciation  */
    s23A = params_3_areas_two_rates[4],     /* between-region speciation  */
    s13A = params_3_areas_two_rates[5],     /* between-region speciation  */        

    s1B  = params_3_areas_two_rates[6],     /* speciation within region 1 */
    s2B  = params_3_areas_two_rates[7],     /* speciation within region 2 */
    s3B  = params_3_areas_two_rates[8],     /* speciation within region 3 */
    s12B = params_3_areas_two_rates[9],     /* between-region speciation  */
    s23B = params_3_areas_two_rates[10],     /* between-region speciation  */
    s13B = params_3_areas_two_rates[11],     /* between-region speciation  */        
    
    // These are true extinction of endemic lineages.
    x1A  = params_3_areas_two_rates[12],     /* extinction from region 1   */
    x2A  = params_3_areas_two_rates[13],     /* extinction from region 2   */
    x3A  = params_3_areas_two_rates[14],     /* extinction from region 3   */

    x1B  = params_3_areas_two_rates[15],     /* extinction from region 1   */
    x2B  = params_3_areas_two_rates[16],     /* extinction from region 2   */
    x3B  = params_3_areas_two_rates[17],     /* extinction from region 3   */

    // Transitions between areas along the branches:
    // Transitions for layer A.
    d1_2A = params_3_areas_two_rates[18],   /* jumps from 1 to 2          */
    d2_1A = params_3_areas_two_rates[19],   /* jumps from 2 to 1          */
    d2_3A = params_3_areas_two_rates[20],   /* jumps from 2 to 3          */
    d3_2A = params_3_areas_two_rates[21],   /* jumps from 3 to 2          */
    d1_3A = params_3_areas_two_rates[22],   /* jumps from 1 to 3          */
    d3_1A = params_3_areas_two_rates[23],   /* jumps from 3 to 1          */
    
    d1_12A = params_3_areas_two_rates[24],  /* dispersal from 1 to 12     */
    d1_13A = params_3_areas_two_rates[25],  /* dispersal from 1 to 13     */
    d2_12A = params_3_areas_two_rates[26],  /* dispersal from 2 to 12     */
    d2_23A = params_3_areas_two_rates[27],  /* dispersal from 2 to 23     */        
    d3_13A = params_3_areas_two_rates[28],  /* dispersal from 3 to 13     */
    d3_23A = params_3_areas_two_rates[29],  /* dispersal from 3 to 23     */

    d12_1A = params_3_areas_two_rates[30],   /* true extirpation rate  */
    d13_1A = params_3_areas_two_rates[31],   /* true extirpation rate  */
    d12_2A = params_3_areas_two_rates[32],   /* true extirpation rate  */
    d23_2A = params_3_areas_two_rates[33],   /* true extirpation rate  */
    d13_3A = params_3_areas_two_rates[34],   /* true extirpation rate  */
    d23_3A = params_3_areas_two_rates[35],   /* true extirpation rate  */

    // Transitions for layer B.
    d1_2B = params_3_areas_two_rates[36],   /* jumps from 1 to 2          */
    d2_1B = params_3_areas_two_rates[37],   /* jumps from 2 to 1          */
    d2_3B = params_3_areas_two_rates[38],   /* jumps from 2 to 3          */
    d3_2B = params_3_areas_two_rates[39],   /* jumps from 3 to 2          */
    d1_3B = params_3_areas_two_rates[40],   /* jumps from 1 to 3          */
    d3_1B = params_3_areas_two_rates[41],   /* jumps from 3 to 1          */
    
    d1_12B = params_3_areas_two_rates[42],  /* dispersal from 1 to 12     */
    d1_13B = params_3_areas_two_rates[43],  /* dispersal from 1 to 13     */
    d2_12B = params_3_areas_two_rates[44],  /* dispersal from 2 to 12     */
    d2_23B = params_3_areas_two_rates[45],  /* dispersal from 2 to 23     */        
    d3_13B = params_3_areas_two_rates[46],  /* dispersal from 3 to 13     */
    d3_23B = params_3_areas_two_rates[47],  /* dispersal from 3 to 23     */

    d12_1B = params_3_areas_two_rates[48],   /* true extirpation rate  */
    d13_1B = params_3_areas_two_rates[49],   /* true extirpation rate  */
    d12_2B = params_3_areas_two_rates[50],   /* true extirpation rate  */
    d23_2B = params_3_areas_two_rates[51],   /* true extirpation rate  */
    d13_3B = params_3_areas_two_rates[52],   /* true extirpation rate  */
    d23_3B = params_3_areas_two_rates[53],   /* true extirpation rate  */

    // These are the transitions between the layers.
    d1A_1B = params_3_areas_two_rates[54],
    d2A_2B = params_3_areas_two_rates[55],
    d3A_3B = params_3_areas_two_rates[56],
    d1B_1A = params_3_areas_two_rates[57],
    d2B_2A = params_3_areas_two_rates[58],
    d3B_3A = params_3_areas_two_rates[59],

    d12A_12B = params_3_areas_two_rates[60],
    d23A_23B = params_3_areas_two_rates[61],
    d13A_13B = params_3_areas_two_rates[62],
    d12B_12A = params_3_areas_two_rates[63],
    d23B_23A = params_3_areas_two_rates[64],
    d13B_13A = params_3_areas_two_rates[65];
  
  // The extinction ODEs

  // Extinction of the endemic lineages:
  
  /*  dE_1A / dt  */
  ydot[0] = -(s1A + d1_2A + d1_3A + d1_12A + d1_13A + x1A + d1A_1B) * E_1A + x1A + d1A_1B * E_1B + d1_2A * E_2A + d1_3A * E_3A + d1_12A * E_12A + d1_13A * E_13A + s1A * E_1A * E_1A;

  /*  dE_2A / dt  */
  ydot[1] = -(s2A + d2_1A + d2_3A + d2_12A + d2_23A + x2A + d2A_2B) * E_2A + x2A + d2A_2B * E_2B + d2_1A * E_1A + d2_3A * E_3A + d2_12A * E_12A + d2_23A * E_23A + s2A * E_2A * E_2A;

  /*  dE_3A / dt  */
  ydot[2] = -(s3A + d3_1A + d3_2A + d3_13A + d3_23A + x3A + d3A_3B) * E_3A + x3A + d3A_3B * E_3B + d3_1A * E_1A + d3_2A * E_2A + d3_13A * E_13A + d3_23A * E_23A + s3A * E_3A * E_3A;

  /*  dE_12A / dt  */
  ydot[3] = -(s1A + s2A + s12A + d12_1A + d12_2A + d12A_12B) * E_12A + d12_1A * E_1A + d12_2A * E_2A + d12A_12B * E_12B + s1A * E_1A * E_12A + s2A * E_2A * E_12A + s12A * E_1A * E_2A;

  /*  dE_23A / dt  */
  ydot[4] = -(s2A + s3A + s23A + d23_2A + d23_3A + d23A_23B) * E_23A + d23_2A * E_2A + d23_3A * E_3A + d23A_23B * E_23B + s2A * E_2A * E_23A + s3A * E_3A * E_23A + s23A * E_2A * E_3A;

  /*  dE_13A / dt  */
  ydot[5] = -(s1A + s3A + s13A + d13_1A + d13_3A + d13A_13B) * E_13A + d13_1A * E_1A + d13_3A * E_3A + d13A_13B * E_13B  + s1A * E_1A * E_13A + s3A * E_3A * E_13A + s13A * E_1A * E_3A;
  
  /*  dE_1B / dt  */
  ydot[6] = -(s1B + d1_2B + d1_3B + d1_12B + d1_13B + x1B + d1B_1A) * E_1B + x1B + d1B_1A * E_1A + d1_2B * E_2B + d1_3B * E_3B + d1_12B * E_12B + d1_13B * E_13B + s1B * E_1B * E_1B;

  /*  dE_2B / dt  */
  ydot[7] = -(s2B + d2_1B + d2_3B + d2_12B + d2_23B + x2B + d2B_2A) * E_2B + x2B + d2B_2A * E_2A + d2_1B * E_1B + d2_3B * E_3B + d2_12B * E_12B + d2_23B * E_23B + s2B * E_2B * E_2B;

  /*  dE_3B / dt  */
  ydot[8] = -(s3B + d3_1B + d3_2B + d3_13B + d3_23B + x3B + d3B_3A) * E_3B + x3B + d3B_3A * E_3A + d3_1B * E_1B + d3_2B * E_2B + d3_13B * E_13B + d3_23B * E_23B + s3B * E_3B * E_3B;

  /*  dE_12B / dt  */
  ydot[9] = -(s1B + s2B + s12B + d12_1B + d12_2B + d12B_12A) * E_12B + d12_1B * E_1B + d12_2B * E_2B + d12B_12A * E_12A + s1B * E_1B * E_12B + s2B * E_2B * E_12B + s12B * E_1B * E_2B;

  /*  dE_23B / dt  */
  ydot[10] = -(s2B + s3B + s23B + d23_2B + d23_3B + d23B_23A) * E_23B + d23_2B * E_2B + d23_3B * E_3B + d23B_23A * E_23A + s2B * E_2B * E_23B + s3B * E_3B * E_23B + s23B * E_2B * E_3B;

  /*  dE_13B / dt  */
  ydot[11] = -(s1B + s3B + s13B + d13_1B + d13_3B + d13B_13A) * E_13B + d13_1B * E_1B + d13_3B * E_3B + d13B_13A * E_13A  + s1B * E_1B * E_13B + s3B * E_3B * E_13B + s13B * E_1B * E_3B;

  // The events ODEs - dispersion and speciation
  
  /*  dD_N1A / dt  */
  ydot[12] = -(s1A + d1_2A + d1_3A + d1_12A + d1_13A + x1A + d1A_1B) * D_N1A + d1_3A * D_N3A + d1_2A * D_N2A + d1A_1B * D_N1B + d1_13A * D_N13A + d1_12A * D_N12A + 2.0 * s1A * D_N1A * E_1A;

  /*  dD_N2A / dt  */
  ydot[13] = -(s2A + d2_1A + d2_3A + d2_12A + d2_23A + x2A + d2A_2B) * D_N2A + d2_3A * D_N3A + d2_1A * D_N1A + d2_23A * D_N23A + d2_12A * D_N12A + d2A_2B + D_N2B + 2.0 * s2A * D_N2A * E_2A;

  /*  dD_N3A / dt  */
  ydot[14] = -(s3A + d3_1A + d3_2A + d3_23A + d3_13A + x3A + d3A_3B) * D_N3A + d3_1A * D_N1A + d3_2A * D_N2A + d3_23A * D_N23A + d3_13A * D_N13A + d3A_3B * D_N3B + 2.0 * s3A * D_N3A * E_3A;

  /*  dD_N12A / dt  */
  ydot[15] = -(s1A + s2A + s12A + d12_1A + d12_2A + d12A_12B) * D_N12A + d12_1A * D_N1A + d12_2A * D_N2A + d12A_12B * D_N12B + s12A * (D_N1A * E_2A + D_N2A * E_1A) + s1A * (E_1A * D_N12A + E_12A * D_N1A) + s2A * (E_2A * D_N12A + E_12A * D_N2A);

  /*  dD_N23A / dt  */
  ydot[16] = -(s2A + s3A + s23A + d23_2A + d23_3A + d23A_23B) * D_N23A + d23_2A * D_N2A + d23_3A * D_N3A + d23A_23B * D_N23B + s23A * (D_N2A * E_3A + D_N3A * E_2A) + s2A * (E_2A * D_N23A + E_23A * D_N2A) + s3A * (E_3A * D_N23A + E_23A * D_N3A);

  /*  dD_N13A / dt  */
  ydot[17] = -(s1A + s3A + s13A + d13_1A + d13_3A + d13A_13B) * D_N13A + d13_1A * D_N1A + d13_3A * D_N3A + d13A_13B * D_N13B + s13A * (D_N1A * E_3A + D_N3A * E_1A) + s1A * (E_1A * D_N13A + E_13A * D_N1A) + s3A * (E_3A * D_N13A + E_13A * D_N3A);

  /*  dD_N1B / dt  */
  ydot[18] = -(s1B + d1_2B + d1_3B + d1_12B + d1_13B + x1B + d1B_1A) * D_N1B + d1_3B * D_N3B + d1_2B * D_N2B + d1B_1A * D_N1A + d1_13B * D_N13B + d1_12B * D_N12B + 2.0 * s1B * D_N1B * E_1B;

  /*  dD_N2B / dt  */
  ydot[19] = -(s2B + d2_1B + d2_3B + d2_12B + d2_23B + x2B + d2B_2A) * D_N2B + d2_3B * D_N3B + d2_1B * D_N1B + d2_23B * D_N23B + d2_12B * D_N12B + d2B_2A + D_N2A + 2.0 * s2B * D_N2B * E_2B;

  /*  dD_N3B / dt  */
  ydot[20] = -(s3B + d3_1B + d3_2B + d3_23B + d3_13B + x3B + d3B_3A) * D_N3B + d3_1B * D_N1B + d3_2B * D_N2B + d3_23B * D_N23B + d3_13B * D_N13B + d3B_3A * D_N3A + 2.0 * s3B * D_N3B * E_3B;

  /*  dD_N12B / dt  */
  ydot[21] = -(s1B + s2B + s12B + d12_1B + d12_2B + d12B_12A) * D_N12B + d12_1B * D_N1B + d12_2B * D_N2B + d12B_12A * D_N12A + s12B * (D_N1B * E_2B + D_N2B * E_1B) + s1B * (E_1B * D_N12B + E_12B * D_N1B) + s2B * (E_2B * D_N12B + E_12B * D_N2B);

  /*  dD_N23B / dt  */
  ydot[22] = -(s2B + s3B + s23B + d23_2B + d23_3B + d23B_23A) * D_N23B + d23_2B * D_N2B + d23_3B * D_N3B + d23B_23A * D_N23A + s23B * (D_N2B * E_3B + D_N3B * E_2B) + s2B * (E_2B * D_N23B + E_23B * D_N2B) + s3B * (E_3B * D_N23B + E_23B * D_N3B);

  /*  dD_N13B / dt  */
  ydot[23] = -(s1B + s3B + s13B + d13_1B + d13_3B + d13B_13A) * D_N13B + d13_1B * D_N1B + d13_3B * D_N3B + d13B_13A * D_N13A + s13B * (D_N1B * E_3B + D_N3B * E_1B) + s1B * (E_1B * D_N13B + E_13B * D_N1B) + s3B * (E_3B * D_N13B + E_13B * D_N3B);
}

