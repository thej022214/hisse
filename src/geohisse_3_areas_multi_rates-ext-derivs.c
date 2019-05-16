/*
 *  geosse_3_areas_multi_rates-ext-derivs.c
 *
 *  By Daniel Caetano Feb/2019
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 336
// This is the total number of parameters for the 12 hidden states model.

static double params_3_areas_multi_rates[NUMELEMENTS];

void initmod_geosse_3_areas_multi_rates(void (* odeparms)(int *, double *)){
  int N = NUMELEMENTS;
  odeparms(&N, params_3_areas_multi_rates);
}

void geosse_3_areas_multi_rates_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){

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

    E_1C = y[12],
    E_2C = y[13],
    E_3C = y[14],
    E_12C = y[15],
    E_23C = y[16],
    E_13C = y[17],

    E_1D = y[18],
    E_2D = y[19],
    E_3D = y[20],
    E_12D = y[21],
    E_23D = y[22],
    E_13D = y[23],            

    E_1E = y[24],
    E_2E = y[25],
    E_3E = y[26],
    E_12E = y[27],
    E_23E = y[28],
    E_13E = y[29],

    E_1F = y[30],
    E_2F = y[31],
    E_3F = y[32],
    E_12F = y[33],
    E_23F = y[34],
    E_13F = y[35],            

    E_1G = y[36],
    E_2G = y[37],
    E_3G = y[38],
    E_12G = y[39],
    E_23G = y[40],
    E_13G = y[41],

    E_1H = y[42],
    E_2H = y[43],
    E_3H = y[44],
    E_12H = y[45],
    E_23H = y[46],
    E_13H = y[47],

    E_1I = y[48],
    E_2I = y[49],
    E_3I = y[50],
    E_12I = y[51],
    E_23I = y[52],
    E_13I = y[53],

    E_1J = y[54],
    E_2J = y[55],
    E_3J = y[56],
    E_12J = y[57],
    E_23J = y[58],
    E_13J = y[59],

    E_1K = y[60],
    E_2K = y[61],
    E_3K = y[62],
    E_12K = y[63],
    E_23K = y[64],
    E_13K = y[65],

    E_1L = y[66],
    E_2L = y[67],
    E_3L = y[68],
    E_12L = y[69],
    E_23L = y[70],
    E_13L = y[71],

    /* The probability of the state on time t. */
    
    D_N1A = y[72],
    D_N2A = y[73],
    D_N3A = y[74],
    D_N12A = y[75],
    D_N23A = y[76],
    D_N13A = y[77],

    D_N1B = y[78],
    D_N2B = y[79],
    D_N3B = y[80],
    D_N12B = y[81],
    D_N23B = y[82],
    D_N13B = y[83],            

    D_N1C = y[84],
    D_N2C = y[85],
    D_N3C = y[86],
    D_N12C = y[87],
    D_N23C = y[88],
    D_N13C = y[89],

    D_N1D = y[90],
    D_N2D = y[91],
    D_N3D = y[92],
    D_N12D = y[93],
    D_N23D = y[94],
    D_N13D = y[95],            

    D_N1E = y[96],
    D_N2E = y[97],
    D_N3E = y[98],
    D_N12E = y[99],
    D_N23E = y[100],
    D_N13E = y[101],

    D_N1F = y[102],
    D_N2F = y[103],
    D_N3F = y[104],
    D_N12F = y[105],
    D_N23F = y[106],
    D_N13F = y[107],            

    D_N1G = y[108],
    D_N2G = y[109],
    D_N3G = y[110],
    D_N12G = y[111],
    D_N23G = y[112],
    D_N13G = y[113],

    D_N1H = y[114],
    D_N2H = y[115],
    D_N3H = y[116],
    D_N12H = y[117],
    D_N23H = y[118],
    D_N13H = y[119],

    D_N1I = y[120],
    D_N2I = y[121],
    D_N3I = y[122],
    D_N12I = y[123],
    D_N23I = y[124],
    D_N13I = y[125],

    D_N1J = y[126],
    D_N2J = y[127],
    D_N3J = y[128],
    D_N12J = y[129],
    D_N23J = y[130],
    D_N13J = y[131],

    D_N1K = y[132],
    D_N2K = y[133],
    D_N3K = y[134],
    D_N12K = y[135],
    D_N23K = y[136],
    D_N13K = y[137],

    D_N1L = y[138],
    D_N2L = y[139],
    D_N3L = y[140],
    D_N12L = y[141],
    D_N23L = y[142],
    D_N13L = y[143],
    
    /* The speciation parameters. */
   
    s1A  = params_3_areas_multi_rates[0],     /* speciation within region 1 */
    s2A  = params_3_areas_multi_rates[1],     /* speciation within region 2 */
    s3A  = params_3_areas_multi_rates[2],     /* speciation within region 3 */
    s12A = params_3_areas_multi_rates[3],     /* between-region speciation  */
    s23A = params_3_areas_multi_rates[4],     /* between-region speciation  */
    s13A = params_3_areas_multi_rates[5],     /* between-region speciation  */        

    s1B  = params_3_areas_multi_rates[6],
    s2B  = params_3_areas_multi_rates[7],
    s3B  = params_3_areas_multi_rates[8],
    s12B = params_3_areas_multi_rates[9],
    s23B = params_3_areas_multi_rates[10],
    s13B = params_3_areas_multi_rates[11],

    s1C  = params_3_areas_multi_rates[12],
    s2C  = params_3_areas_multi_rates[13],
    s3C  = params_3_areas_multi_rates[14],
    s12C = params_3_areas_multi_rates[15],
    s23C = params_3_areas_multi_rates[16],
    s13C = params_3_areas_multi_rates[17],

    s1D  = params_3_areas_multi_rates[18],
    s2D  = params_3_areas_multi_rates[19],
    s3D  = params_3_areas_multi_rates[20],
    s12D = params_3_areas_multi_rates[21],
    s23D = params_3_areas_multi_rates[22],
    s13D = params_3_areas_multi_rates[23],

    s1E  = params_3_areas_multi_rates[24],
    s2E  = params_3_areas_multi_rates[25],
    s3E  = params_3_areas_multi_rates[26],
    s12E = params_3_areas_multi_rates[27],
    s23E = params_3_areas_multi_rates[28],
    s13E = params_3_areas_multi_rates[29],
    
    s1F  = params_3_areas_multi_rates[30],
    s2F  = params_3_areas_multi_rates[31],
    s3F  = params_3_areas_multi_rates[32],
    s12F = params_3_areas_multi_rates[33],
    s23F = params_3_areas_multi_rates[34],
    s13F = params_3_areas_multi_rates[35],

    s1G  = params_3_areas_multi_rates[36],
    s2G  = params_3_areas_multi_rates[37],
    s3G  = params_3_areas_multi_rates[38],
    s12G = params_3_areas_multi_rates[39],
    s23G = params_3_areas_multi_rates[40],
    s13G = params_3_areas_multi_rates[41],

    s1H  = params_3_areas_multi_rates[42],
    s2H  = params_3_areas_multi_rates[43],
    s3H  = params_3_areas_multi_rates[44],
    s12H = params_3_areas_multi_rates[45],
    s23H = params_3_areas_multi_rates[46],
    s13H = params_3_areas_multi_rates[47],

    s1I  = params_3_areas_multi_rates[48],
    s2I  = params_3_areas_multi_rates[49],
    s3I  = params_3_areas_multi_rates[50],
    s12I = params_3_areas_multi_rates[51],
    s23I = params_3_areas_multi_rates[52],
    s13I = params_3_areas_multi_rates[53],
    
    s1J  = params_3_areas_multi_rates[54],
    s2J  = params_3_areas_multi_rates[55],
    s3J  = params_3_areas_multi_rates[56],
    s12J = params_3_areas_multi_rates[57],
    s23J = params_3_areas_multi_rates[58],
    s13J = params_3_areas_multi_rates[59],

    s1K  = params_3_areas_multi_rates[60],
    s2K  = params_3_areas_multi_rates[61],
    s3K  = params_3_areas_multi_rates[62],
    s12K = params_3_areas_multi_rates[63],
    s23K = params_3_areas_multi_rates[64],
    s13K = params_3_areas_multi_rates[65],

    s1L  = params_3_areas_multi_rates[66],
    s2L  = params_3_areas_multi_rates[67],
    s3L  = params_3_areas_multi_rates[68],
    s12L = params_3_areas_multi_rates[69],
    s23L = params_3_areas_multi_rates[70],
    s13L = params_3_areas_multi_rates[71],
    
    // These are true extinction of endemic lineages.
    x1A  = params_3_areas_multi_rates[72],     /* extinction from region 1   */
    x2A  = params_3_areas_multi_rates[73],     /* extinction from region 2   */
    x3A  = params_3_areas_multi_rates[74],     /* extinction from region 3   */

    x1B  = params_3_areas_multi_rates[75],
    x2B  = params_3_areas_multi_rates[76],
    x3B  = params_3_areas_multi_rates[77],

    x1C  = params_3_areas_multi_rates[78],
    x2C  = params_3_areas_multi_rates[79],
    x3C  = params_3_areas_multi_rates[80],

    x1D  = params_3_areas_multi_rates[81],
    x2D  = params_3_areas_multi_rates[82],
    x3D  = params_3_areas_multi_rates[83],
    
    x1E  = params_3_areas_multi_rates[84],
    x2E  = params_3_areas_multi_rates[85],
    x3E  = params_3_areas_multi_rates[86],

    x1F  = params_3_areas_multi_rates[87],
    x2F  = params_3_areas_multi_rates[88],
    x3F  = params_3_areas_multi_rates[89],

    x1G  = params_3_areas_multi_rates[90],
    x2G  = params_3_areas_multi_rates[91],
    x3G  = params_3_areas_multi_rates[92],

    x1H  = params_3_areas_multi_rates[93],
    x2H  = params_3_areas_multi_rates[94],
    x3H  = params_3_areas_multi_rates[95],

    x1I  = params_3_areas_multi_rates[96],
    x2I  = params_3_areas_multi_rates[97],
    x3I  = params_3_areas_multi_rates[98],

    x1J  = params_3_areas_multi_rates[99],
    x2J  = params_3_areas_multi_rates[100],
    x3J  = params_3_areas_multi_rates[101],

    x1K  = params_3_areas_multi_rates[102],
    x2K  = params_3_areas_multi_rates[103],
    x3K  = params_3_areas_multi_rates[104],

    x1L  = params_3_areas_multi_rates[105],
    x2L  = params_3_areas_multi_rates[106],
    x3L  = params_3_areas_multi_rates[107],

    // Transitions between areas along the branches:
    // Transitions for layer A.
    d1_2A = params_3_areas_multi_rates[108],   /* jumps from 1 to 2          */
    d2_1A = params_3_areas_multi_rates[109],   /* jumps from 2 to 1          */
    d2_3A = params_3_areas_multi_rates[110],   /* jumps from 2 to 3          */
    d3_2A = params_3_areas_multi_rates[111],   /* jumps from 3 to 2          */
    d1_3A = params_3_areas_multi_rates[112],   /* jumps from 1 to 3          */
    d3_1A = params_3_areas_multi_rates[113],   /* jumps from 3 to 1          */
    
    d1_12A = params_3_areas_multi_rates[114],  /* dispersal from 1 to 12     */
    d1_13A = params_3_areas_multi_rates[115],  /* dispersal from 1 to 13     */
    d2_12A = params_3_areas_multi_rates[116],  /* dispersal from 2 to 12     */
    d2_23A = params_3_areas_multi_rates[117],  /* dispersal from 2 to 23     */        
    d3_13A = params_3_areas_multi_rates[118],  /* dispersal from 3 to 13     */
    d3_23A = params_3_areas_multi_rates[119],  /* dispersal from 3 to 23     */

    d12_1A = params_3_areas_multi_rates[120],   /* true extirpation rate  */
    d13_1A = params_3_areas_multi_rates[121],   /* true extirpation rate  */
    d12_2A = params_3_areas_multi_rates[122],   /* true extirpation rate  */
    d23_2A = params_3_areas_multi_rates[123],   /* true extirpation rate  */
    d13_3A = params_3_areas_multi_rates[124],   /* true extirpation rate  */
    d23_3A = params_3_areas_multi_rates[125],   /* true extirpation rate  */

    d1_2B = params_3_areas_multi_rates[126],
    d2_1B = params_3_areas_multi_rates[127],
    d2_3B = params_3_areas_multi_rates[128],
    d3_2B = params_3_areas_multi_rates[129],
    d1_3B = params_3_areas_multi_rates[130],
    d3_1B = params_3_areas_multi_rates[131],
    
    d1_12B = params_3_areas_multi_rates[132],
    d1_13B = params_3_areas_multi_rates[133],
    d2_12B = params_3_areas_multi_rates[134],
    d2_23B = params_3_areas_multi_rates[135],
    d3_13B = params_3_areas_multi_rates[136],
    d3_23B = params_3_areas_multi_rates[137],

    d12_1B = params_3_areas_multi_rates[138],
    d13_1B = params_3_areas_multi_rates[139],
    d12_2B = params_3_areas_multi_rates[140],
    d23_2B = params_3_areas_multi_rates[141],
    d13_3B = params_3_areas_multi_rates[142],
    d23_3B = params_3_areas_multi_rates[143],

    d1_2C = params_3_areas_multi_rates[144],
    d2_1C = params_3_areas_multi_rates[145],
    d2_3C = params_3_areas_multi_rates[146],
    d3_2C = params_3_areas_multi_rates[147],
    d1_3C = params_3_areas_multi_rates[148],
    d3_1C = params_3_areas_multi_rates[149],
    
    d1_12C = params_3_areas_multi_rates[150],
    d1_13C = params_3_areas_multi_rates[151],
    d2_12C = params_3_areas_multi_rates[152],
    d2_23C = params_3_areas_multi_rates[153],
    d3_13C = params_3_areas_multi_rates[154],
    d3_23C = params_3_areas_multi_rates[155],

    d12_1C = params_3_areas_multi_rates[156],
    d13_1C = params_3_areas_multi_rates[157],
    d12_2C = params_3_areas_multi_rates[158],
    d23_2C = params_3_areas_multi_rates[159],
    d13_3C = params_3_areas_multi_rates[160],
    d23_3C = params_3_areas_multi_rates[161],

    d1_2D = params_3_areas_multi_rates[162],
    d2_1D = params_3_areas_multi_rates[163],
    d2_3D = params_3_areas_multi_rates[164],
    d3_2D = params_3_areas_multi_rates[165],
    d1_3D = params_3_areas_multi_rates[166],
    d3_1D = params_3_areas_multi_rates[167],
    
    d1_12D = params_3_areas_multi_rates[168],
    d1_13D = params_3_areas_multi_rates[169],
    d2_12D = params_3_areas_multi_rates[170],
    d2_23D = params_3_areas_multi_rates[171],
    d3_13D = params_3_areas_multi_rates[172],
    d3_23D = params_3_areas_multi_rates[173],

    d12_1D = params_3_areas_multi_rates[174],
    d13_1D = params_3_areas_multi_rates[175],
    d12_2D = params_3_areas_multi_rates[176],
    d23_2D = params_3_areas_multi_rates[177],
    d13_3D = params_3_areas_multi_rates[178],
    d23_3D = params_3_areas_multi_rates[179],

    d1_2E = params_3_areas_multi_rates[180],
    d2_1E = params_3_areas_multi_rates[181],
    d2_3E = params_3_areas_multi_rates[182],
    d3_2E = params_3_areas_multi_rates[183],
    d1_3E = params_3_areas_multi_rates[184],
    d3_1E = params_3_areas_multi_rates[185],
    
    d1_12E = params_3_areas_multi_rates[186],
    d1_13E = params_3_areas_multi_rates[187],
    d2_12E = params_3_areas_multi_rates[188],
    d2_23E = params_3_areas_multi_rates[189],
    d3_13E = params_3_areas_multi_rates[190],
    d3_23E = params_3_areas_multi_rates[191],

    d12_1E = params_3_areas_multi_rates[192],
    d13_1E = params_3_areas_multi_rates[193],
    d12_2E = params_3_areas_multi_rates[194],
    d23_2E = params_3_areas_multi_rates[195],
    d13_3E = params_3_areas_multi_rates[196],
    d23_3E = params_3_areas_multi_rates[197],

    d1_2F = params_3_areas_multi_rates[198],
    d2_1F = params_3_areas_multi_rates[199],
    d2_3F = params_3_areas_multi_rates[200],
    d3_2F = params_3_areas_multi_rates[201],
    d1_3F = params_3_areas_multi_rates[202],
    d3_1F = params_3_areas_multi_rates[203],
    
    d1_12F = params_3_areas_multi_rates[204],
    d1_13F = params_3_areas_multi_rates[205],
    d2_12F = params_3_areas_multi_rates[206],
    d2_23F = params_3_areas_multi_rates[207],
    d3_13F = params_3_areas_multi_rates[208],
    d3_23F = params_3_areas_multi_rates[209],

    d12_1F = params_3_areas_multi_rates[210],
    d13_1F = params_3_areas_multi_rates[211],
    d12_2F = params_3_areas_multi_rates[212],
    d23_2F = params_3_areas_multi_rates[213],
    d13_3F = params_3_areas_multi_rates[214],
    d23_3F = params_3_areas_multi_rates[215],

    d1_2G = params_3_areas_multi_rates[216],
    d2_1G = params_3_areas_multi_rates[217],
    d2_3G = params_3_areas_multi_rates[218],
    d3_2G = params_3_areas_multi_rates[219],
    d1_3G = params_3_areas_multi_rates[220],
    d3_1G = params_3_areas_multi_rates[221],
    
    d1_12G = params_3_areas_multi_rates[222],
    d1_13G = params_3_areas_multi_rates[223],
    d2_12G = params_3_areas_multi_rates[224],
    d2_23G = params_3_areas_multi_rates[225],
    d3_13G = params_3_areas_multi_rates[226],
    d3_23G = params_3_areas_multi_rates[227],

    d12_1G = params_3_areas_multi_rates[228],
    d13_1G = params_3_areas_multi_rates[229],
    d12_2G = params_3_areas_multi_rates[230],
    d23_2G = params_3_areas_multi_rates[231],
    d13_3G = params_3_areas_multi_rates[232],
    d23_3G = params_3_areas_multi_rates[233],

    d1_2H = params_3_areas_multi_rates[234],
    d2_1H = params_3_areas_multi_rates[235],
    d2_3H = params_3_areas_multi_rates[236],
    d3_2H = params_3_areas_multi_rates[237],
    d1_3H = params_3_areas_multi_rates[238],
    d3_1H = params_3_areas_multi_rates[239],
    
    d1_12H = params_3_areas_multi_rates[240],
    d1_13H = params_3_areas_multi_rates[241],
    d2_12H = params_3_areas_multi_rates[242],
    d2_23H = params_3_areas_multi_rates[243],
    d3_13H = params_3_areas_multi_rates[244],
    d3_23H = params_3_areas_multi_rates[245],

    d12_1H = params_3_areas_multi_rates[246],
    d13_1H = params_3_areas_multi_rates[247],
    d12_2H = params_3_areas_multi_rates[248],
    d23_2H = params_3_areas_multi_rates[249],
    d13_3H = params_3_areas_multi_rates[250],
    d23_3H = params_3_areas_multi_rates[251],

    d1_2I = params_3_areas_multi_rates[252],
    d2_1I = params_3_areas_multi_rates[253],
    d2_3I = params_3_areas_multi_rates[254],
    d3_2I = params_3_areas_multi_rates[255],
    d1_3I = params_3_areas_multi_rates[256],
    d3_1I = params_3_areas_multi_rates[257],
    
    d1_12I = params_3_areas_multi_rates[258],
    d1_13I = params_3_areas_multi_rates[259],
    d2_12I = params_3_areas_multi_rates[260],
    d2_23I = params_3_areas_multi_rates[261],
    d3_13I = params_3_areas_multi_rates[262],
    d3_23I = params_3_areas_multi_rates[263],

    d12_1I = params_3_areas_multi_rates[264],
    d13_1I = params_3_areas_multi_rates[265],
    d12_2I = params_3_areas_multi_rates[266],
    d23_2I = params_3_areas_multi_rates[267],
    d13_3I = params_3_areas_multi_rates[268],
    d23_3I = params_3_areas_multi_rates[269],

    d1_2J = params_3_areas_multi_rates[270],
    d2_1J = params_3_areas_multi_rates[271],
    d2_3J = params_3_areas_multi_rates[272],
    d3_2J = params_3_areas_multi_rates[273],
    d1_3J = params_3_areas_multi_rates[274],
    d3_1J = params_3_areas_multi_rates[275],
    
    d1_12J = params_3_areas_multi_rates[276],
    d1_13J = params_3_areas_multi_rates[277],
    d2_12J = params_3_areas_multi_rates[278],
    d2_23J = params_3_areas_multi_rates[279],
    d3_13J = params_3_areas_multi_rates[280],
    d3_23J = params_3_areas_multi_rates[281],

    d12_1J = params_3_areas_multi_rates[282],
    d13_1J = params_3_areas_multi_rates[283],
    d12_2J = params_3_areas_multi_rates[284],
    d23_2J = params_3_areas_multi_rates[285],
    d13_3J = params_3_areas_multi_rates[286],
    d23_3J = params_3_areas_multi_rates[287],

    d1_2K = params_3_areas_multi_rates[288],
    d2_1K = params_3_areas_multi_rates[289],
    d2_3K = params_3_areas_multi_rates[290],
    d3_2K = params_3_areas_multi_rates[291],
    d1_3K = params_3_areas_multi_rates[292],
    d3_1K = params_3_areas_multi_rates[293],
    
    d1_12K = params_3_areas_multi_rates[294],
    d1_13K = params_3_areas_multi_rates[295],
    d2_12K = params_3_areas_multi_rates[296],
    d2_23K = params_3_areas_multi_rates[297],
    d3_13K = params_3_areas_multi_rates[298],
    d3_23K = params_3_areas_multi_rates[299],

    d12_1K = params_3_areas_multi_rates[300],
    d13_1K = params_3_areas_multi_rates[301],
    d12_2K = params_3_areas_multi_rates[302],
    d23_2K = params_3_areas_multi_rates[303],
    d13_3K = params_3_areas_multi_rates[304],
    d23_3K = params_3_areas_multi_rates[305],

    d1_2L = params_3_areas_multi_rates[306],
    d2_1L = params_3_areas_multi_rates[307],
    d2_3L = params_3_areas_multi_rates[308],
    d3_2L = params_3_areas_multi_rates[309],
    d1_3L = params_3_areas_multi_rates[310],
    d3_1L = params_3_areas_multi_rates[311],
    
    d1_12L = params_3_areas_multi_rates[312],
    d1_13L = params_3_areas_multi_rates[313],
    d2_12L = params_3_areas_multi_rates[314],
    d2_23L = params_3_areas_multi_rates[315],
    d3_13L = params_3_areas_multi_rates[316],
    d3_23L = params_3_areas_multi_rates[317],

    d12_1L = params_3_areas_multi_rates[318],
    d13_1L = params_3_areas_multi_rates[319],
    d12_2L = params_3_areas_multi_rates[320],
    d23_2L = params_3_areas_multi_rates[321],
    d13_3L = params_3_areas_multi_rates[322],
    d23_3L = params_3_areas_multi_rates[323],
    
    /* These are the transitions between the layers. */
    /* Need to use one parameter for each layer to prevent the model of visiting layers that are not included in the model. This can be problematic for the inference. */
    dLAY_A = params_3_areas_multi_rates[324],
    dLAY_B = params_3_areas_multi_rates[325],
    dLAY_C = params_3_areas_multi_rates[326],
    dLAY_D = params_3_areas_multi_rates[327],
    dLAY_E = params_3_areas_multi_rates[328],
    dLAY_F = params_3_areas_multi_rates[329],
    dLAY_G = params_3_areas_multi_rates[330],
    dLAY_H = params_3_areas_multi_rates[331],
    dLAY_I = params_3_areas_multi_rates[332],
    dLAY_J = params_3_areas_multi_rates[333],
    dLAY_K = params_3_areas_multi_rates[334],
    dLAY_L = params_3_areas_multi_rates[335];

  /* Now need to include all the hidden rates as well as transitions between these layers. Note that the rate of transition between the layers is global. So we have a single parameter to model the rate that these layers jump between each other. */
  /* A desirable upgrade for this program would be a dynamic way to reduce the integration elements following the probabilities of the process, such that the integration would not take place including a hidden layer that has 0 probability. */
  
  // The extinction ODEs

  // Extinction of the endemic lineages:
  
  /*  dE_1A / dt  */
  ydot[0] = - (s1A + d1_2A + d1_3A + d1_12A + d1_13A + x1A) * E_1A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1A
    + x1A
    + d1_2A * E_2A
    + d1_3A * E_3A
    + d1_12A * E_12A
    + d1_13A * E_13A
    + s1A * E_1A * E_1A
    + dLAY_B * E_1B
    + dLAY_C * E_1C
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2A / dt  */
  ydot[1] = -(s2A + d2_1A + d2_3A + d2_12A + d2_23A + x2A) * E_2A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2A
    + x2A
    + d2_1A * E_1A
    + d2_3A * E_3A
    + d2_12A * E_12A
    + d2_23A * E_23A
    + s2A * E_2A * E_2A
    + dLAY_B * E_2B
    + dLAY_C * E_2C
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3A / dt  */
  ydot[2] = -(s3A + d3_1A + d3_2A + d3_13A + d3_23A + x3A) * E_3A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3A
    + x3A
    + d3_1A * E_1A
    + d3_2A * E_2A
    + d3_13A * E_13A
    + d3_23A * E_23A
    + s3A * E_3A * E_3A
    + dLAY_B * E_3B
    + dLAY_C * E_3C
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12A / dt  */
  ydot[3] = -(s1A + s2A + s12A + d12_1A + d12_2A) * E_12A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12A
    + d12_1A * E_1A
    + d12_2A * E_2A
    + s1A * E_1A * E_12A
    + s2A * E_2A * E_12A
    + s12A * E_1A * E_2A
    + dLAY_B * E_12B
    + dLAY_C * E_12C
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;
  
  /*  dE_23A / dt  */
  ydot[4] = -(s2A + s3A + s23A + d23_2A + d23_3A) * E_23A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23A
    + d23_2A * E_2A
    + d23_3A * E_3A
    + s2A * E_2A * E_23A
    + s3A * E_3A * E_23A
    + s23A * E_2A * E_3A
    + dLAY_B * E_23B
    + dLAY_C * E_23C
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13A / dt  */
  ydot[5] = -(s1A + s3A + s13A + d13_1A + d13_3A) * E_13A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13A
    + d13_1A * E_1A
    + d13_3A * E_3A
    + s1A * E_1A * E_13A
    + s3A * E_3A * E_13A
    + s13A * E_1A * E_3A
    + dLAY_B * E_13B
    + dLAY_C * E_13C
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;
  
  /*  dE_1B / dt  */
  ydot[6] = -(s1B + d1_2B + d1_3B + d1_12B + d1_13B + x1B) * E_1B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1B
    + x1B
    + d1_2B * E_2B
    + d1_3B * E_3B
    + d1_12B * E_12B
    + d1_13B * E_13B
    + s1B * E_1B * E_1B
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2B / dt  */
  ydot[7] = -(s2B + d2_1B + d2_3B + d2_12B + d2_23B + x2B) * E_2B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2B
    + x2B
    + d2_1B * E_1B
    + d2_3B * E_3B
    + d2_12B * E_12B
    + d2_23B * E_23B
    + s2B * E_2B * E_2B
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3B / dt  */
  ydot[8] = -(s3B + d3_1B + d3_2B + d3_13B + d3_23B + x3B) * E_3B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3B
    + x3B
    + d3_1B * E_1B
    + d3_2B * E_2B
    + d3_13B * E_13B
    + d3_23B * E_23B
    + s3B * E_3B * E_3B
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12B / dt  */
  ydot[9] = -(s1B + s2B + s12B + d12_1B + d12_2B) * E_12B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12B
    + d12_1B * E_1B
    + d12_2B * E_2B
    + s1B * E_1B * E_12B
    + s2B * E_2B * E_12B
    + s12B * E_1B * E_2B
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23B / dt  */
  ydot[10] = -(s2B + s3B + s23B + d23_2B + d23_3B) * E_23B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23B
    + d23_2B * E_2B
    + d23_3B * E_3B
    + s2B * E_2B * E_23B
    + s3B * E_3B * E_23B
    + s23B * E_2B * E_3B
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13B / dt  */
  ydot[11] = -(s1B + s3B + s13B + d13_1B + d13_3B) * E_13B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13B
    + d13_1B * E_1B
    + d13_3B * E_3B
    + s1B * E_1B * E_13B
    + s3B * E_3B * E_13B
    + s13B * E_1B * E_3B
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

    /*  dE_1C / dt  */
  ydot[12] = -(s1C + d1_2C + d1_3C + d1_12C + d1_13C + x1C) * E_1C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1C
    + x1C
    + d1_2C * E_2C
    + d1_3C * E_3C
    + d1_12C * E_12C
    + d1_13C * E_13C
    + s1C * E_1C * E_1C
    + dLAY_A * E_1A
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2C / dt  */
  ydot[13] = -(s2C + d2_1C + d2_3C + d2_12C + d2_23C + x2C) * E_2C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2C
    + x2C
    + d2_1C * E_1C
    + d2_3C * E_3C
    + d2_12C * E_12C
    + d2_23C * E_23C
    + s2C * E_2C * E_2C
    + dLAY_A * E_2A
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3C / dt  */
  ydot[14] = -(s3C + d3_1C + d3_2C + d3_13C + d3_23C + x3C) * E_3C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3C
    + x3C
    + d3_1C * E_1C
    + d3_2C * E_2C
    + d3_13C * E_13C
    + d3_23C * E_23C
    + s3C * E_3C * E_3C
    + dLAY_A * E_3A
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12C / dt  */
  ydot[15] = -(s1C + s2C + s12C + d12_1C + d12_2C) * E_12C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12C
    + d12_1C * E_1C
    + d12_2C * E_2C
    + s1C * E_1C * E_12C
    + s2C * E_2C * E_12C
    + s12C * E_1C * E_2C
    + dLAY_A * E_12A
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23C / dt  */
  ydot[16] = -(s2C + s3C + s23C + d23_2C + d23_3C) * E_23C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23C
    + d23_2C * E_2C
    + d23_3C * E_3C
    + s2C * E_2C * E_23C
    + s3C * E_3C * E_23C
    + s23C * E_2C * E_3C
    + dLAY_A * E_23A
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13C / dt  */
  ydot[17] = -(s1C + s3C + s13C + d13_1C + d13_3C) * E_13C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13C
    + d13_1C * E_1C
    + d13_3C * E_3C
    + s1C * E_1C * E_13C
    + s3C * E_3C * E_13C
    + s13C * E_1C * E_3C
    + dLAY_A * E_13A
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

    /*  dE_1D / dt  */
  ydot[18] = -(s1D + d1_2D + d1_3D + d1_12D + d1_13D + x1D) * E_1D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1D
    + x1D
    + d1_2D * E_2D
    + d1_3D * E_3D
    + d1_12D * E_12D
    + d1_13D * E_13D
    + s1D * E_1D * E_1D
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2D / dt  */
  ydot[19] = -(s2D + d2_1D + d2_3D + d2_12D + d2_23D + x2D) * E_2D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2D
    + x2D
    + d2_1D * E_1D
    + d2_3D * E_3D
    + d2_12D * E_12D
    + d2_23D * E_23D
    + s2D * E_2D * E_2D
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3D / dt  */
  ydot[20] = -(s3D + d3_1D + d3_2D + d3_13D + d3_23D + x3D) * E_3D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3D
    + x3D
    + d3_1D * E_1D
    + d3_2D * E_2D
    + d3_13D * E_13D
    + d3_23D * E_23D
    + s3D * E_3D * E_3D
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12D / dt  */
  ydot[21] = -(s1D + s2D + s12D + d12_1D + d12_2D) * E_12D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12D
    + d12_1D * E_1D
    + d12_2D * E_2D
    + s1D * E_1D * E_12D
    + s2D * E_2D * E_12D
    + s12D * E_1D * E_2D
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23D / dt  */
  ydot[22] = -(s2D + s3D + s23D + d23_2D + d23_3D) * E_23D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23D
    + d23_2D * E_2D
    + d23_3D * E_3D
    + s2D * E_2D * E_23D
    + s3D * E_3D * E_23D
    + s23D * E_2D * E_3D
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13D / dt  */
  ydot[23] = -(s1D + s3D + s13D + d13_1D + d13_3D) * E_13D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13D
    + d13_1D * E_1D
    + d13_3D * E_3D
    + s1D * E_1D * E_13D
    + s3D * E_3D * E_13D
    + s13D * E_1D * E_3D
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

  /*  dE_1E / dt  */
  ydot[24] = -(s1E + d1_2E + d1_3E + d1_12E + d1_13E + x1E) * E_1E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1E
    + x1E
    + d1_2E * E_2E
    + d1_3E * E_3E
    + d1_12E * E_12E
    + d1_13E * E_13E
    + s1E * E_1E * E_1E
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2E / dt  */
  ydot[25] = -(s2E + d2_1E + d2_3E + d2_12E + d2_23E + x2E) * E_2E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2E
    + x2E
    + d2_1E * E_1E
    + d2_3E * E_3E
    + d2_12E * E_12E
    + d2_23E * E_23E
    + s2E * E_2E * E_2E
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3E / dt  */
  ydot[26] = -(s3E + d3_1E + d3_2E + d3_13E + d3_23E + x3E) * E_3E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3E
    + x3E
    + d3_1E * E_1E
    + d3_2E * E_2E
    + d3_13E * E_13E
    + d3_23E * E_23E
    + s3E * E_3E * E_3E
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12E / dt  */
  ydot[27] = -(s1E + s2E + s12E + d12_1E + d12_2E) * E_12E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12E
    + d12_1E * E_1E
    + d12_2E * E_2E
    + s1E * E_1E * E_12E
    + s2E * E_2E * E_12E
    + s12E * E_1E * E_2E
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23E / dt  */
  ydot[28] = -(s2E + s3E + s23E + d23_2E + d23_3E) * E_23E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23E
    + d23_2E * E_2E
    + d23_3E * E_3E
    + s2E * E_2E * E_23E
    + s3E * E_3E * E_23E
    + s23E * E_2E * E_3E
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13E / dt  */
  ydot[29] = -(s1E + s3E + s13E + d13_1E + d13_3E) * E_13E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13E
    + d13_1E * E_1E
    + d13_3E * E_3E
    + s1E * E_1E * E_13E
    + s3E * E_3E * E_13E
    + s13E * E_1E * E_3E
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

  /*  dE_1F / dt  */
  ydot[30] = -(s1F + d1_2F + d1_3F + d1_12F + d1_13F + x1F) * E_1F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1F
    + x1F
    + d1_2F * E_2F
    + d1_3F * E_3F
    + d1_12F * E_12F
    + d1_13F * E_13F
    + s1F * E_1F * E_1F
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2F / dt  */
  ydot[31] = -(s2F + d2_1F + d2_3F + d2_12F + d2_23F + x2F) * E_2F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2F
    + x2F
    + d2_1F * E_1F
    + d2_3F * E_3F
    + d2_12F * E_12F
    + d2_23F * E_23F
    + s2F * E_2F * E_2F
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3F / dt  */
  ydot[32] = -(s3F + d3_1F + d3_2F + d3_13F + d3_23F + x3F) * E_3F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3F
    + x3F
    + d3_1F * E_1F
    + d3_2F * E_2F
    + d3_13F * E_13F
    + d3_23F * E_23F
    + s3F * E_3F * E_3F
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12F / dt  */
  ydot[33] = -(s1F + s2F + s12F + d12_1F + d12_2F) * E_12F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12F
    + d12_1F * E_1F
    + d12_2F * E_2F
    + s1F * E_1F * E_12F
    + s2F * E_2F * E_12F
    + s12F * E_1F * E_2F
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23F / dt  */
  ydot[34] = -(s2F + s3F + s23F + d23_2F + d23_3F) * E_23F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23F
    + d23_2F * E_2F
    + d23_3F * E_3F
    + s2F * E_2F * E_23F
    + s3F * E_3F * E_23F
    + s23F * E_2F * E_3F
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13F / dt  */
  ydot[35] = -(s1F + s3F + s13F + d13_1F + d13_3F) * E_13F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13F
    + d13_1F * E_1F
    + d13_3F * E_3F
    + s1F * E_1F * E_13F
    + s3F * E_3F * E_13F
    + s13F * E_1F * E_3F
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

  /*  dE_1G / dt  */
  ydot[36] = -(s1G + d1_2G + d1_3G + d1_12G + d1_13G + x1G) * E_1G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1G
    + x1G
    + d1_2G * E_2G
    + d1_3G * E_3G
    + d1_12G * E_12G
    + d1_13G * E_13G
    + s1G * E_1G * E_1G
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2G / dt  */
  ydot[37] = -(s2G + d2_1G + d2_3G + d2_12G + d2_23G + x2G) * E_2G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2G
    + x2G
    + d2_1G * E_1G
    + d2_3G * E_3G
    + d2_12G * E_12G
    + d2_23G * E_23G
    + s2G * E_2G * E_2G
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3G / dt  */
  ydot[38] = -(s3G + d3_1G + d3_2G + d3_13G + d3_23G + x3G) * E_3G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3G
    + x3G
    + d3_1G * E_1G
    + d3_2G * E_2G
    + d3_13G * E_13G
    + d3_23G * E_23G
    + s3G * E_3G * E_3G
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12G / dt  */
  ydot[39] = -(s1G + s2G + s12G + d12_1G + d12_2G) * E_12G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12G
    + d12_1G * E_1G
    + d12_2G * E_2G
    + s1G * E_1G * E_12G
    + s2G * E_2G * E_12G
    + s12G * E_1G * E_2G
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23G / dt  */
  ydot[40] = -(s2G + s3G + s23G + d23_2G + d23_3G) * E_23G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23G
    + d23_2G * E_2G
    + d23_3G * E_3G
    + s2G * E_2G * E_23G
    + s3G * E_3G * E_23G
    + s23G * E_2G * E_3G
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13G / dt  */
  ydot[41] = -(s1G + s3G + s13G + d13_1G + d13_3G) * E_13G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13G
    + d13_1G * E_1G
    + d13_3G * E_3G
    + s1G * E_1G * E_13G
    + s3G * E_3G * E_13G
    + s13G * E_1G * E_3G
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

  /*  dE_1H / dt  */
  ydot[42] = -(s1H + d1_2H + d1_3H + d1_12H + d1_13H + x1H) * E_1H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_1H
    + x1H
    + d1_2H * E_2H
    + d1_3H * E_3H
    + d1_12H * E_12H
    + d1_13H * E_13H
    + s1H * E_1H * E_1H
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2H / dt  */
  ydot[43] = -(s2H + d2_1H + d2_3H + d2_12H + d2_23H + x2H) * E_2H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_2H
    + x2H
    + d2_1H * E_1H
    + d2_3H * E_3H
    + d2_12H * E_12H
    + d2_23H * E_23H
    + s2H * E_2H * E_2H
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3H / dt  */
  ydot[44] = -(s3H + d3_1H + d3_2H + d3_13H + d3_23H + x3H) * E_3H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_3H
    + x3H
    + d3_1H * E_1H
    + d3_2H * E_2H
    + d3_13H * E_13H
    + d3_23H * E_23H
    + s3H * E_3H * E_3H
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12H / dt  */
  ydot[45] = -(s1H + s2H + s12H + d12_1H + d12_2H) * E_12H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_12H
    + d12_1H * E_1H
    + d12_2H * E_2H
    + s1H * E_1H * E_12H
    + s2H * E_2H * E_12H
    + s12H * E_1H * E_2H
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23H / dt  */
  ydot[46] = -(s2H + s3H + s23H + d23_2H + d23_3H) * E_23H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_23H
    + d23_2H * E_2H
    + d23_3H * E_3H
    + s2H * E_2H * E_23H
    + s3H * E_3H * E_23H
    + s23H * E_2H * E_3H
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13H / dt  */
  ydot[47] = -(s1H + s3H + s13H + d13_1H + d13_3H) * E_13H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * E_13H
    + d13_1H * E_1H
    + d13_3H * E_3H
    + s1H * E_1H * E_13H
    + s3H * E_3H * E_13H
    + s13H * E_1H * E_3H
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

    /*  dE_1I / dt  */
  ydot[48] = -(s1I + d1_2I + d1_3I + d1_12I + d1_13I + x1I) * E_1I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * E_1I
    + x1I
    + d1_2I * E_2I
    + d1_3I * E_3I
    + d1_12I * E_12I
    + d1_13I * E_13I
    + s1I * E_1I * E_1I
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_J * E_1J
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2I / dt  */
  ydot[49] = -(s2I + d2_1I + d2_3I + d2_12I + d2_23I + x2I) * E_2I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * E_2I
    + x2I
    + d2_1I * E_1I
    + d2_3I * E_3I
    + d2_12I * E_12I
    + d2_23I * E_23I
    + s2I * E_2I * E_2I
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_J * E_2J
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3I / dt  */
  ydot[50] = -(s3I + d3_1I + d3_2I + d3_13I + d3_23I + x3I) * E_3I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * E_3I
    + x3I
    + d3_1I * E_1I
    + d3_2I * E_2I
    + d3_13I * E_13I
    + d3_23I * E_23I
    + s3I * E_3I * E_3I
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_J * E_3J
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12I / dt  */
  ydot[51] = -(s1I + s2I + s12I + d12_1I + d12_2I) * E_12I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * E_12I
    + d12_1I * E_1I
    + d12_2I * E_2I
    + s1I * E_1I * E_12I
    + s2I * E_2I * E_12I
    + s12I * E_1I * E_2I
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_J * E_12J
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23I / dt  */
  ydot[52] = -(s2I + s3I + s23I + d23_2I + d23_3I) * E_23I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * E_23I
    + d23_2I * E_2I
    + d23_3I * E_3I
    + s2I * E_2I * E_23I
    + s3I * E_3I * E_23I
    + s23I * E_2I * E_3I
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_J * E_23J
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13I / dt  */
  ydot[53] = -(s1I + s3I + s13I + d13_1I + d13_3I) * E_13I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * E_13I
    + d13_1I * E_1I
    + d13_3I * E_3I
    + s1I * E_1I * E_13I
    + s3I * E_3I * E_13I
    + s13I * E_1I * E_3I
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_J * E_13J
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

  /*  dE_1J / dt  */
  ydot[54] = -(s1J + d1_2J + d1_3J + d1_12J + d1_13J + x1J) * E_1J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * E_1J
    + x1J
    + d1_2J * E_2J
    + d1_3J * E_3J
    + d1_12J * E_12J
    + d1_13J * E_13J
    + s1J * E_1J * E_1J
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_K * E_1K
    + dLAY_L * E_1L;

  /*  dE_2J / dt  */
  ydot[55] = -(s2J + d2_1J + d2_3J + d2_12J + d2_23J + x2J) * E_2J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * E_2J
    + x2J
    + d2_1J * E_1J
    + d2_3J * E_3J
    + d2_12J * E_12J
    + d2_23J * E_23J
    + s2J * E_2J * E_2J
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_K * E_2K
    + dLAY_L * E_2L;

  /*  dE_3J / dt  */
  ydot[56] = -(s3J + d3_1J + d3_2J + d3_13J + d3_23J + x3J) * E_3J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * E_3J
    + x3J
    + d3_1J * E_1J
    + d3_2J * E_2J
    + d3_13J * E_13J
    + d3_23J * E_23J
    + s3J * E_3J * E_3J
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_K * E_3K
    + dLAY_L * E_3L;

  /*  dE_12J / dt  */
  ydot[57] = -(s1J + s2J + s12J + d12_1J + d12_2J) * E_12J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * E_12J
    + d12_1J * E_1J
    + d12_2J * E_2J
    + s1J * E_1J * E_12J
    + s2J * E_2J * E_12J
    + s12J * E_1J * E_2J
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_K * E_12K
    + dLAY_L * E_12L;

  /*  dE_23J / dt  */
  ydot[58] = -(s2J + s3J + s23J + d23_2J + d23_3J) * E_23J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * E_23J
    + d23_2J * E_2J
    + d23_3J * E_3J
    + s2J * E_2J * E_23J
    + s3J * E_3J * E_23J
    + s23J * E_2J * E_3J
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_K * E_23K
    + dLAY_L * E_23L;

  /*  dE_13J / dt  */
  ydot[59] = -(s1J + s3J + s13J + d13_1J + d13_3J) * E_13J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * E_13J
    + d13_1J * E_1J
    + d13_3J * E_3J
    + s1J * E_1J * E_13J
    + s3J * E_3J * E_13J
    + s13J * E_1J * E_3J
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_K * E_13K
    + dLAY_L * E_13L;

  /*  dE_1K / dt  */
  ydot[60] = -(s1K + d1_2K + d1_3K + d1_12K + d1_13K + x1K) * E_1K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * E_1K
    + x1K
    + d1_2K * E_2K
    + d1_3K * E_3K
    + d1_12K * E_12K
    + d1_13K * E_13K
    + s1K * E_1K * E_1K
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_L * E_1L;

  /*  dE_2K / dt  */
  ydot[61] = -(s2K + d2_1K + d2_3K + d2_12K + d2_23K + x2K) * E_2K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * E_2K
    + x2K
    + d2_1K * E_1K
    + d2_3K * E_3K
    + d2_12K * E_12K
    + d2_23K * E_23K
    + s2K * E_2K * E_2K
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_L * E_2L;

  /*  dE_3K / dt  */
  ydot[62] = -(s3K + d3_1K + d3_2K + d3_13K + d3_23K + x3K) * E_3K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * E_3K
    + x3K
    + d3_1K * E_1K
    + d3_2K * E_2K
    + d3_13K * E_13K
    + d3_23K * E_23K
    + s3K * E_3K * E_3K
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_L * E_3L;

  /*  dE_12K / dt  */
  ydot[63] = -(s1K + s2K + s12K + d12_1K + d12_2K) * E_12K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * E_12K
    + d12_1K * E_1K
    + d12_2K * E_2K
    + s1K * E_1K * E_12K
    + s2K * E_2K * E_12K
    + s12K * E_1K * E_2K
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_L * E_12L;

  /*  dE_23K / dt  */
  ydot[64] = -(s2K + s3K + s23K + d23_2K + d23_3K) * E_23K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * E_23K
    + d23_2K * E_2K
    + d23_3K * E_3K
    + s2K * E_2K * E_23K
    + s3K * E_3K * E_23K
    + s23K * E_2K * E_3K
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_L * E_23L;

  /*  dE_13K / dt  */
  ydot[65] = -(s1K + s3K + s13K + d13_1K + d13_3K) * E_13K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * E_13K
    + d13_1K * E_1K
    + d13_3K * E_3K
    + s1K * E_1K * E_13K
    + s3K * E_3K * E_13K
    + s13K * E_1K * E_3K
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_L * E_13L;

  /*  dE_1L / dt  */
  ydot[66] = -(s1L + d1_2L + d1_3L + d1_12L + d1_13L + x1L) * E_1L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * E_1L
    + x1L
    + d1_2L * E_2L
    + d1_3L * E_3L
    + d1_12L * E_12L
    + d1_13L * E_13L
    + s1L * E_1L * E_1L
    + dLAY_A * E_1A
    + dLAY_C * E_1C
    + dLAY_B * E_1B
    + dLAY_D * E_1D
    + dLAY_E * E_1E
    + dLAY_F * E_1F
    + dLAY_G * E_1G
    + dLAY_H * E_1H
    + dLAY_I * E_1I
    + dLAY_J * E_1J
    + dLAY_L * E_1K;

  /*  dE_2L / dt  */
  ydot[67] = -(s2L + d2_1L + d2_3L + d2_12L + d2_23L + x2L) * E_2L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * E_2L
    + x2L
    + d2_1L * E_1L
    + d2_3L * E_3L
    + d2_12L * E_12L
    + d2_23L * E_23L
    + s2L * E_2L * E_2L
    + dLAY_A * E_2A
    + dLAY_C * E_2C
    + dLAY_B * E_2B
    + dLAY_D * E_2D
    + dLAY_E * E_2E
    + dLAY_F * E_2F
    + dLAY_G * E_2G
    + dLAY_H * E_2H
    + dLAY_I * E_2I
    + dLAY_J * E_2J
    + dLAY_K * E_2K;

  /*  dE_3L / dt  */
  ydot[68] = -(s3L + d3_1L + d3_2L + d3_13L + d3_23L + x3L) * E_3L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * E_3L
    + x3L
    + d3_1L * E_1L
    + d3_2L * E_2L
    + d3_13L * E_13L
    + d3_23L * E_23L
    + s3L * E_3L * E_3L
    + dLAY_A * E_3A
    + dLAY_C * E_3C
    + dLAY_B * E_3B
    + dLAY_D * E_3D
    + dLAY_E * E_3E
    + dLAY_F * E_3F
    + dLAY_G * E_3G
    + dLAY_H * E_3H
    + dLAY_I * E_3I
    + dLAY_J * E_3J
    + dLAY_K * E_3K;

  /*  dE_12L / dt  */
  ydot[69] = -(s1L + s2L + s12L + d12_1L + d12_2L) * E_12L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * E_12L
    + d12_1L * E_1L
    + d12_2L * E_2L
    + s1L * E_1L * E_12L
    + s2L * E_2L * E_12L
    + s12L * E_1L * E_2L
    + dLAY_A * E_12A
    + dLAY_C * E_12C
    + dLAY_B * E_12B
    + dLAY_D * E_12D
    + dLAY_E * E_12E
    + dLAY_F * E_12F
    + dLAY_G * E_12G
    + dLAY_H * E_12H
    + dLAY_I * E_12I
    + dLAY_J * E_12J
    + dLAY_K * E_12K;

  /*  dE_23L / dt  */
  ydot[70] = -(s2L + s3L + s23L + d23_2L + d23_3L) * E_23L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * E_23L
    + d23_2L * E_2L
    + d23_3L * E_3L
    + s2L * E_2L * E_23L
    + s3L * E_3L * E_23L
    + s23L * E_2L * E_3L
    + dLAY_A * E_23A
    + dLAY_C * E_23C
    + dLAY_B * E_23B
    + dLAY_D * E_23D
    + dLAY_E * E_23E
    + dLAY_F * E_23F
    + dLAY_G * E_23G
    + dLAY_H * E_23H
    + dLAY_I * E_23I
    + dLAY_J * E_23J
    + dLAY_K * E_23K;

  /*  dE_13L / dt  */
  ydot[71] = -(s1L + s3L + s13L + d13_1L + d13_3L) * E_13L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * E_13L
    + d13_1L * E_1L
    + d13_3L * E_3L
    + s1L * E_1L * E_13L
    + s3L * E_3L * E_13L
    + s13L * E_1L * E_3L
    + dLAY_A * E_13A
    + dLAY_C * E_13C
    + dLAY_B * E_13B
    + dLAY_D * E_13D
    + dLAY_E * E_13E
    + dLAY_F * E_13F
    + dLAY_G * E_13G
    + dLAY_H * E_13H
    + dLAY_I * E_13I
    + dLAY_J * E_13J
    + dLAY_K * E_13K;
  
  // The events ODEs - dispersion and speciation
  
  /*  dD_N1A / dt  */
  ydot[72] = -(s1A + d1_2A + d1_3A + d1_12A + d1_13A + x1A) * D_N1A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1A
    + d1_3A * D_N3A
    + d1_2A * D_N2A
    + d1_13A * D_N13A
    + d1_12A * D_N12A
    + 2.0 * s1A * D_N1A * E_1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2A / dt  */
  ydot[73] = -(s2A + d2_1A + d2_3A + d2_12A + d2_23A + x2A) * D_N2A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2A
    + d2_3A * D_N3A
    + d2_1A * D_N1A
    + d2_23A * D_N23A
    + d2_12A * D_N12A
    + 2.0 * s2A * D_N2A * E_2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;
    
  /*  dD_N3A / dt  */
  ydot[74] = -(s3A + d3_1A + d3_2A + d3_23A + d3_13A + x3A) * D_N3A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3A
    + d3_1A * D_N1A
    + d3_2A * D_N2A
    + d3_23A * D_N23A
    + d3_13A * D_N13A
    + 2.0 * s3A * D_N3A * E_3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12A / dt  */
  ydot[75] = -(s1A + s2A + s12A + d12_1A + d12_2A) * D_N12A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12A
    + d12_1A * D_N1A
    + d12_2A * D_N2A
    + s12A * (D_N1A * E_2A + D_N2A * E_1A)
    + s1A * (E_1A * D_N12A + E_12A * D_N1A)
    + s2A * (E_2A * D_N12A + E_12A * D_N2A)
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23A / dt  */
  ydot[76] = -(s2A + s3A + s23A + d23_2A + d23_3A) * D_N23A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23A
    + d23_2A * D_N2A
    + d23_3A * D_N3A
    + s23A * (D_N2A * E_3A + D_N3A * E_2A)
    + s2A * (E_2A * D_N23A + E_23A * D_N2A)
    + s3A * (E_3A * D_N23A + E_23A * D_N3A)
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13A / dt  */
  ydot[77] = -(s1A + s3A + s13A + d13_1A + d13_3A) * D_N13A
    - (dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13A
    + d13_1A * D_N1A
    + d13_3A * D_N3A
    + s13A * (D_N1A * E_3A + D_N3A * E_1A)
    + s1A * (E_1A * D_N13A + E_13A * D_N1A)
    + s3A * (E_3A * D_N13A + E_13A * D_N3A)
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

  /*  dD_N1B / dt  */
  ydot[78] = -(s1B + d1_2B + d1_3B + d1_12B + d1_13B + x1B) * D_N1B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1B
    + d1_3B * D_N3B
    + d1_2B * D_N2B
    + d1_13B * D_N13B
    + d1_12B * D_N12B
    + 2.0 * s1B * D_N1B * E_1B
    + dLAY_A * D_N1A
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2B / dt  */
  ydot[79] = -(s2B + d2_1B + d2_3B + d2_12B + d2_23B + x2B) * D_N2B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2B
    + d2_3B * D_N3B
    + d2_1B * D_N1B
    + d2_23B * D_N23B
    + d2_12B * D_N12B
    + 2.0 * s2B * D_N2B * E_2B
    + dLAY_A * D_N2A
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3B / dt  */
  ydot[80] = -(s3B + d3_1B + d3_2B + d3_23B + d3_13B + x3B) * D_N3B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3B
    + d3_1B * D_N1B
    + d3_2B * D_N2B
    + d3_23B * D_N23B
    + d3_13B * D_N13B
    + 2.0 * s3B * D_N3B * E_3B
    + dLAY_A * D_N3A
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12B / dt  */
  ydot[81] = -(s1B + s2B + s12B + d12_1B + d12_2B) * D_N12B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12B
    + d12_1B * D_N1B
    + d12_2B * D_N2B
    + s12B * (D_N1B * E_2B + D_N2B * E_1B)
    + s1B * (E_1B * D_N12B + E_12B * D_N1B)
    + s2B * (E_2B * D_N12B + E_12B * D_N2B)
    + dLAY_A * D_N12A
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23B / dt  */
  ydot[82] = -(s2B + s3B + s23B + d23_2B + d23_3B) * D_N23B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23B
    + d23_2B * D_N2B
    + d23_3B * D_N3B
    + s23B * (D_N2B * E_3B + D_N3B * E_2B)
    + s2B * (E_2B * D_N23B + E_23B * D_N2B)
    + s3B * (E_3B * D_N23B + E_23B * D_N3B)
    + dLAY_A * D_N23A
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13B / dt  */
  ydot[83] = -(s1B + s3B + s13B + d13_1B + d13_3B) * D_N13B
    - (dLAY_A + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13B
    + d13_1B * D_N1B
    + d13_3B * D_N3B
    + s13B * (D_N1B * E_3B + D_N3B * E_1B)
    + s1B * (E_1B * D_N13B + E_13B * D_N1B)
    + s3B * (E_3B * D_N13B + E_13B * D_N3B)
    + dLAY_A * D_N13A
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

    /*  dD_N1C / dt  */
  ydot[84] = -(s1C + d1_2C + d1_3C + d1_12C + d1_13C + x1C) * D_N1C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1C
    + d1_3C * D_N3C
    + d1_2C * D_N2C
    + d1_13C * D_N13C
    + d1_12C * D_N12C
    + 2.0 * s1C * D_N1C * E_1C
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2C / dt  */
  ydot[85] = -(s2C + d2_1C + d2_3C + d2_12C + d2_23C + x2C) * D_N2C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2C
    + d2_3C * D_N3C
    + d2_1C * D_N1C
    + d2_23C * D_N23C
    + d2_12C * D_N12C
    + 2.0 * s2C * D_N2C * E_2C
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3C / dt  */
  ydot[86] = -(s3C + d3_1C + d3_2C + d3_23C + d3_13C + x3C) * D_N3C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3C
    + d3_1C * D_N1C
    + d3_2C * D_N2C
    + d3_23C * D_N23C
    + d3_13C * D_N13C
    + 2.0 * s3C * D_N3C * E_3C
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12C / dt  */
  ydot[87] = -(s1C + s2C + s12C + d12_1C + d12_2C) * D_N12C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12C
    + d12_1C * D_N1C
    + d12_2C * D_N2C
    + s12C * (D_N1C * E_2C + D_N2C * E_1C)
    + s1C * (E_1C * D_N12C + E_12C * D_N1C)
    + s2C * (E_2C * D_N12C + E_12C * D_N2C)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23C / dt  */
  ydot[88] = -(s2C + s3C + s23C + d23_2C + d23_3C) * D_N23C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23C
    + d23_2C * D_N2C
    + d23_3C * D_N3C
    + s23C * (D_N2C * E_3C + D_N3C * E_2C)
    + s2C * (E_2C * D_N23C + E_23C * D_N2C)
    + s3C * (E_3C * D_N23C + E_23C * D_N3C)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13C / dt  */
  ydot[89] = -(s1C + s3C + s13C + d13_1C + d13_3C) * D_N13C
    - (dLAY_A + dLAY_B + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13C
    + d13_1C * D_N1C
    + d13_3C * D_N3C
    + s13C * (D_N1C * E_3C + D_N3C * E_1C)
    + s1C * (E_1C * D_N13C + E_13C * D_N1C)
    + s3C * (E_3C * D_N13C + E_13C * D_N3C)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

  /*  dD_N1D / dt  */
  ydot[90] = -(s1D + d1_2D + d1_3D + d1_12D + d1_13D + x1D) * D_N1D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1D
    + d1_3D * D_N3D
    + d1_2D * D_N2D
    + d1_13D * D_N13D
    + d1_12D * D_N12D
    + 2.0 * s1D * D_N1D * E_1D
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2D / dt  */
  ydot[91] = -(s2D + d2_1D + d2_3D + d2_12D + d2_23D + x2D) * D_N2D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2D
    + d2_3D * D_N3D
    + d2_1D * D_N1D
    + d2_23D * D_N23D
    + d2_12D * D_N12D
    + 2.0 * s2D * D_N2D * E_2D
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3D / dt  */
  ydot[92] = -(s3D + d3_1D + d3_2D + d3_23D + d3_13D + x3D) * D_N3D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3D
    + d3_1D * D_N1D
    + d3_2D * D_N2D
    + d3_23D * D_N23D
    + d3_13D * D_N13D
    + 2.0 * s3D * D_N3D * E_3D
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12D / dt  */
  ydot[93] = -(s1D + s2D + s12D + d12_1D + d12_2D) * D_N12D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12D
    + d12_1D * D_N1D
    + d12_2D * D_N2D
    + s12D * (D_N1D * E_2D + D_N2D * E_1D)
    + s1D * (E_1D * D_N12D + E_12D * D_N1D)
    + s2D * (E_2D * D_N12D + E_12D * D_N2D)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23D / dt  */
  ydot[94] = -(s2D + s3D + s23D + d23_2D + d23_3D) * D_N23D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23D
    + d23_2D * D_N2D
    + d23_3D * D_N3D
    + s23D * (D_N2D * E_3D + D_N3D * E_2D)
    + s2D * (E_2D * D_N23D + E_23D * D_N2D)
    + s3D * (E_3D * D_N23D + E_23D * D_N3D)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13D / dt  */
  ydot[95] = -(s1D + s3D + s13D + d13_1D + d13_3D) * D_N13D
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13D
    + d13_1D * D_N1D
    + d13_3D * D_N3D
    + s13D * (D_N1D * E_3D + D_N3D * E_1D)
    + s1D * (E_1D * D_N13D + E_13D * D_N1D)
    + s3D * (E_3D * D_N13D + E_13D * D_N3D)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

  /*  dD_N1E / dt  */
  ydot[96] = -(s1E + d1_2E + d1_3E + d1_12E + d1_13E + x1E) * D_N1E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1E
    + d1_3E * D_N3E
    + d1_2E * D_N2E
    + d1_13E * D_N13E
    + d1_12E * D_N12E
    + 2.0 * s1E * D_N1E * E_1E
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2E / dt  */
  ydot[97] = -(s2E + d2_1E + d2_3E + d2_12E + d2_23E + x2E) * D_N2E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2E
    + d2_3E * D_N3E
    + d2_1E * D_N1E
    + d2_23E * D_N23E
    + d2_12E * D_N12E
    + 2.0 * s2E * D_N2E * E_2E
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3E / dt  */
  ydot[98] = -(s3E + d3_1E + d3_2E + d3_23E + d3_13E + x3E) * D_N3E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3E
    + d3_1E * D_N1E
    + d3_2E * D_N2E
    + d3_23E * D_N23E
    + d3_13E * D_N13E
    + 2.0 * s3E * D_N3E * E_3E
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12E / dt  */
  ydot[99] = -(s1E + s2E + s12E + d12_1E + d12_2E) * D_N12E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12E
    + d12_1E * D_N1E
    + d12_2E * D_N2E
    + s12E * (D_N1E * E_2E + D_N2E * E_1E)
    + s1E * (E_1E * D_N12E + E_12E * D_N1E)
    + s2E * (E_2E * D_N12E + E_12E * D_N2E)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23E / dt  */
  ydot[100] = -(s2E + s3E + s23E + d23_2E + d23_3E) * D_N23E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23E
    + d23_2E * D_N2E
    + d23_3E * D_N3E
    + s23E * (D_N2E * E_3E + D_N3E * E_2E)
    + s2E * (E_2E * D_N23E + E_23E * D_N2E)
    + s3E * (E_3E * D_N23E + E_23E * D_N3E)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13E / dt  */
  ydot[101] = -(s1E + s3E + s13E + d13_1E + d13_3E) * D_N13E
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13E
    + d13_1E * D_N1E
    + d13_3E * D_N3E
    + s13E * (D_N1E * E_3E + D_N3E * E_1E)
    + s1E * (E_1E * D_N13E + E_13E * D_N1E)
    + s3E * (E_3E * D_N13E + E_13E * D_N3E)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

  /*  dD_N1F / dt  */
  ydot[102] = -(s1F + d1_2F + d1_3F + d1_12F + d1_13F + x1F) * D_N1F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1F
    + d1_3F * D_N3F
    + d1_2F * D_N2F
    + d1_13F * D_N13F
    + d1_12F * D_N12F
    + 2.0 * s1F * D_N1F * E_1F
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2F / dt  */
  ydot[103] = -(s2F + d2_1F + d2_3F + d2_12F + d2_23F + x2F) * D_N2F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2F
    + d2_3F * D_N3F
    + d2_1F * D_N1F
    + d2_23F * D_N23F
    + d2_12F * D_N12F
    + 2.0 * s2F * D_N2F * E_2F
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3F / dt  */
  ydot[104] = -(s3F + d3_1F + d3_2F + d3_23F + d3_13F + x3F) * D_N3F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3F
    + d3_1F * D_N1F
    + d3_2F * D_N2F
    + d3_23F * D_N23F
    + d3_13F * D_N13F
    + 2.0 * s3F * D_N3F * E_3F
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12F / dt  */
  ydot[105] = -(s1F + s2F + s12F + d12_1F + d12_2F) * D_N12F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12F
    + d12_1F * D_N1F
    + d12_2F * D_N2F
    + s12F * (D_N1F * E_2F + D_N2F * E_1F)
    + s1F * (E_1F * D_N12F + E_12F * D_N1F)
    + s2F * (E_2F * D_N12F + E_12F * D_N2F)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23F / dt  */
  ydot[106] = -(s2F + s3F + s23F + d23_2F + d23_3F) * D_N23F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23F
    + d23_2F * D_N2F
    + d23_3F * D_N3F
    + s23F * (D_N2F * E_3F + D_N3F * E_2F)
    + s2F * (E_2F * D_N23F + E_23F * D_N2F)
    + s3F * (E_3F * D_N23F + E_23F * D_N3F)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13F / dt  */
  ydot[107] = -(s1F + s3F + s13F + d13_1F + d13_3F) * D_N13F
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13F
    + d13_1F * D_N1F
    + d13_3F * D_N3F
    + s13F * (D_N1F * E_3F + D_N3F * E_1F)
    + s1F * (E_1F * D_N13F + E_13F * D_N1F)
    + s3F * (E_3F * D_N13F + E_13F * D_N3F)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

  /*  dD_N1G / dt  */
  ydot[108] = -(s1G + d1_2G + d1_3G + d1_12G + d1_13G + x1G) * D_N1G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1G
    + d1_3G * D_N3G
    + d1_2G * D_N2G
    + d1_13G * D_N13G
    + d1_12G * D_N12G
    + 2.0 * s1G * D_N1G * E_1G
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;

  /*  dD_N2G / dt  */
  ydot[109] = -(s2G + d2_1G + d2_3G + d2_12G + d2_23G + x2G) * D_N2G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2G
    + d2_3G * D_N3G
    + d2_1G * D_N1G
    + d2_23G * D_N23G
    + d2_12G * D_N12G
    + 2.0 * s2G * D_N2G * E_2G
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3G / dt  */
  ydot[110] = -(s3G + d3_1G + d3_2G + d3_23G + d3_13G + x3G) * D_N3G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3G
    + d3_1G * D_N1G
    + d3_2G * D_N2G
    + d3_23G * D_N23G
    + d3_13G * D_N13G
    + 2.0 * s3G * D_N3G * E_3G
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12G / dt  */
  ydot[111] = -(s1G + s2G + s12G + d12_1G + d12_2G) * D_N12G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12G
    + d12_1G * D_N1G
    + d12_2G * D_N2G
    + s12G * (D_N1G * E_2G + D_N2G * E_1G)
    + s1G * (E_1G * D_N12G + E_12G * D_N1G)
    + s2G * (E_2G * D_N12G + E_12G * D_N2G)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;

  /*  dD_N23G / dt  */
  ydot[112] = -(s2G + s3G + s23G + d23_2G + d23_3G) * D_N23G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23G
    + d23_2G * D_N2G
    + d23_3G * D_N3G
    + s23G * (D_N2G * E_3G + D_N3G * E_2G)
    + s2G * (E_2G * D_N23G + E_23G * D_N2G)
    + s3G * (E_3G * D_N23G + E_23G * D_N3G)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;

  /*  dD_N13G / dt  */
  ydot[113] = -(s1G + s3G + s13G + d13_1G + d13_3G) * D_N13G
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_H + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13G
    + d13_1G * D_N1G
    + d13_3G * D_N3G
    + s13G * (D_N1G * E_3G + D_N3G * E_1G)
    + s1G * (E_1G * D_N13G + E_13G * D_N1G)
    + s3G * (E_3G * D_N13G + E_13G * D_N3G)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;

  /*  dD_N1H / dt  */
  ydot[114] = -(s1H + d1_2H + d1_3H + d1_12H + d1_13H + x1H) * D_N1H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N1H
    + d1_3H * D_N3H
    + d1_2H * D_N2H
    + d1_13H * D_N13H
    + d1_12H * D_N12H
    + 2.0 * s1H * D_N1H * E_1H
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;    

  /*  dD_N2H / dt  */
  ydot[115] = -(s2H + d2_1H + d2_3H + d2_12H + d2_23H + x2H) * D_N2H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N2H
    + d2_3H * D_N3H
    + d2_1H * D_N1H
    + d2_23H * D_N23H
    + d2_12H * D_N12H
    + 2.0 * s2H * D_N2H * E_2H
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;

  /*  dD_N3H / dt  */
  ydot[116] = -(s3H + d3_1H + d3_2H + d3_23H + d3_13H + x3H) * D_N3H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N3H
    + d3_1H * D_N1H
    + d3_2H * D_N2H
    + d3_23H * D_N23H
    + d3_13H * D_N13H
    + 2.0 * s3H * D_N3H * E_3H
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;

  /*  dD_N12H / dt  */
  ydot[117] = -(s1H + s2H + s12H + d12_1H + d12_2H) * D_N12H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N12H
    + d12_1H * D_N1H
    + d12_2H * D_N2H
    + s12H * (D_N1H * E_2H + D_N2H * E_1H)
    + s1H * (E_1H * D_N12H + E_12H * D_N1H)
    + s2H * (E_2H * D_N12H + E_12H * D_N2H)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;    

  /*  dD_N23H / dt  */
  ydot[118] = -(s2H + s3H + s23H + d23_2H + d23_3H) * D_N23H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N23H
    + d23_2H * D_N2H
    + d23_3H * D_N3H
    + s23H * (D_N2H * E_3H + D_N3H * E_2H)
    + s2H * (E_2H * D_N23H + E_23H * D_N2H)
    + s3H * (E_3H * D_N23H + E_23H * D_N3H)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;    

  /*  dD_N13H / dt  */
  ydot[119] = -(s1H + s3H + s13H + d13_1H + d13_3H) * D_N13H
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_I + dLAY_J + dLAY_K + dLAY_L) * D_N13H
    + d13_1H * D_N1H
    + d13_3H * D_N3H
    + s13H * (D_N1H * E_3H + D_N3H * E_1H)
    + s1H * (E_1H * D_N13H + E_13H * D_N1H)
    + s3H * (E_3H * D_N13H + E_13H * D_N3H)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;    

  /*  dD_N1I / dt  */
  ydot[120] = -(s1I + d1_2I + d1_3I + d1_12I + d1_13I + x1I) * D_N1I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * D_N1I
    + d1_3I * D_N3I
    + d1_2I * D_N2I
    + d1_13I * D_N13I
    + d1_12I * D_N12I
    + 2.0 * s1I * D_N1I * E_1I
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;    

  /*  dD_N2I / dt  */
  ydot[121] = -(s2I + d2_1I + d2_3I + d2_12I + d2_23I + x2I) * D_N2I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * D_N2I
    + d2_3I * D_N3I
    + d2_1I * D_N1I
    + d2_23I * D_N23I
    + d2_12I * D_N12I
    + 2.0 * s2I * D_N2I * E_2I
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;    

  /*  dD_N3I / dt  */
  ydot[122] = -(s3I + d3_1I + d3_2I + d3_23I + d3_13I + x3I) * D_N3I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * D_N3I
    + d3_1I * D_N1I
    + d3_2I * D_N2I
    + d3_23I * D_N23I
    + d3_13I * D_N13I
    + 2.0 * s3I * D_N3I * E_3I
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;    

  /*  dD_N12I / dt  */
  ydot[123] = -(s1I + s2I + s12I + d12_1I + d12_2I) * D_N12I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * D_N12I
    + d12_1I * D_N1I
    + d12_2I * D_N2I
    + s12I * (D_N1I * E_2I + D_N2I * E_1I)
    + s1I * (E_1I * D_N12I + E_12I * D_N1I)
    + s2I * (E_2I * D_N12I + E_12I * D_N2I)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;    

  /*  dD_N23I / dt  */
  ydot[124] = -(s2I + s3I + s23I + d23_2I + d23_3I) * D_N23I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * D_N23I
    + d23_2I * D_N2I
    + d23_3I * D_N3I
    + s23I * (D_N2I * E_3I + D_N3I * E_2I)
    + s2I * (E_2I * D_N23I + E_23I * D_N2I)
    + s3I * (E_3I * D_N23I + E_23I * D_N3I)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;    

  /*  dD_N13I / dt  */
  ydot[125] = -(s1I + s3I + s13I + d13_1I + d13_3I) * D_N13I
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_J + dLAY_K + dLAY_L) * D_N13I
    + d13_1I * D_N1I
    + d13_3I * D_N3I
    + s13I * (D_N1I * E_3I + D_N3I * E_1I)
    + s1I * (E_1I * D_N13I + E_13I * D_N1I)
    + s3I * (E_3I * D_N13I + E_13I * D_N3I)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;    

  /*  dD_N1J / dt  */
  ydot[126] = -(s1J + d1_2J + d1_3J + d1_12J + d1_13J + x1J) * D_N1J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * D_N1J
    + d1_3J * D_N3J
    + d1_2J * D_N2J
    + d1_13J * D_N13J
    + d1_12J * D_N12J
    + 2.0 * s1J * D_N1J * E_1J
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_K * D_N1K
    + dLAY_L * D_N1L;    

  /*  dD_N2J / dt  */
  ydot[127] = -(s2J + d2_1J + d2_3J + d2_12J + d2_23J + x2J) * D_N2J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * D_N2J
    + d2_3J * D_N3J
    + d2_1J * D_N1J
    + d2_23J * D_N23J
    + d2_12J * D_N12J
    + 2.0 * s2J * D_N2J * E_2J
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_K * D_N2K
    + dLAY_L * D_N2L;    

  /*  dD_N3J / dt  */
  ydot[128] = -(s3J + d3_1J + d3_2J + d3_23J + d3_13J + x3J) * D_N3J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * D_N3J
    + d3_1J * D_N1J
    + d3_2J * D_N2J
    + d3_23J * D_N23J
    + d3_13J * D_N13J
    + 2.0 * s3J * D_N3J * E_3J
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_K * D_N3K
    + dLAY_L * D_N3L;    

  /*  dD_N12J / dt  */
  ydot[129] = -(s1J + s2J + s12J + d12_1J + d12_2J) * D_N12J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * D_N12J
    + d12_1J * D_N1J
    + d12_2J * D_N2J
    + s12J * (D_N1J * E_2J + D_N2J * E_1J)
    + s1J * (E_1J * D_N12J + E_12J * D_N1J)
    + s2J * (E_2J * D_N12J + E_12J * D_N2J)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_K * D_N12K
    + dLAY_L * D_N12L;    

  /*  dD_N23J / dt  */
  ydot[130] = -(s2J + s3J + s23J + d23_2J + d23_3J) * D_N23J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * D_N23J
    + d23_2J * D_N2J
    + d23_3J * D_N3J
    + s23J * (D_N2J * E_3J + D_N3J * E_2J)
    + s2J * (E_2J * D_N23J + E_23J * D_N2J)
    + s3J * (E_3J * D_N23J + E_23J * D_N3J)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_K * D_N23K
    + dLAY_L * D_N23L;    

  /*  dD_N13J / dt  */
  ydot[131] = -(s1J + s3J + s13J + d13_1J + d13_3J) * D_N13J
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_K + dLAY_L) * D_N13J
    + d13_1J * D_N1J
    + d13_3J * D_N3J
    + s13J * (D_N1J * E_3J + D_N3J * E_1J)
    + s1J * (E_1J * D_N13J + E_13J * D_N1J)
    + s3J * (E_3J * D_N13J + E_13J * D_N3J)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_K * D_N13K
    + dLAY_L * D_N13L;    

  /*  dD_N1K / dt  */
  ydot[132] = -(s1K + d1_2K + d1_3K + d1_12K + d1_13K + x1K) * D_N1K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * D_N1K
    + d1_3K * D_N3K
    + d1_2K * D_N2K
    + d1_13K * D_N13K
    + d1_12K * D_N12K
    + 2.0 * s1K * D_N1K * E_1K
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_L * D_N1L;    

  /*  dD_N2K / dt  */
  ydot[133] = -(s2K + d2_1K + d2_3K + d2_12K + d2_23K + x2K) * D_N2K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * D_N2K
    + d2_3K * D_N3K
    + d2_1K * D_N1K
    + d2_23K * D_N23K
    + d2_12K * D_N12K
    + 2.0 * s2K * D_N2K * E_2K
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_L * D_N2L;    

  /*  dD_N3K / dt  */
  ydot[134] = -(s3K + d3_1K + d3_2K + d3_23K + d3_13K + x3K) * D_N3K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * D_N3K
    + d3_1K * D_N1K
    + d3_2K * D_N2K
    + d3_23K * D_N23K
    + d3_13K * D_N13K
    + 2.0 * s3K * D_N3K * E_3K
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_L * D_N3L;    

  /*  dD_N12K / dt  */
  ydot[135] = -(s1K + s2K + s12K + d12_1K + d12_2K) * D_N12K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * D_N12K
    + d12_1K * D_N1K
    + d12_2K * D_N2K
    + s12K * (D_N1K * E_2K + D_N2K * E_1K)
    + s1K * (E_1K * D_N12K + E_12K * D_N1K)
    + s2K * (E_2K * D_N12K + E_12K * D_N2K)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_L * D_N12L;    

  /*  dD_N23K / dt  */
  ydot[136] = -(s2K + s3K + s23K + d23_2K + d23_3K) * D_N23K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * D_N23K
    + d23_2K * D_N2K
    + d23_3K * D_N3K
    + s23K * (D_N2K * E_3K + D_N3K * E_2K)
    + s2K * (E_2K * D_N23K + E_23K * D_N2K)
    + s3K * (E_3K * D_N23K + E_23K * D_N3K)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_L * D_N23L;    

  /*  dD_N13K / dt  */
  ydot[137] = -(s1K + s3K + s13K + d13_1K + d13_3K) * D_N13K
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_L) * D_N13K
    + d13_1K * D_N1K
    + d13_3K * D_N3K
    + s13K * (D_N1K * E_3K + D_N3K * E_1K)
    + s1K * (E_1K * D_N13K + E_13K * D_N1K)
    + s3K * (E_3K * D_N13K + E_13K * D_N3K)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_L * D_N13L;    

  /*  dD_N1L / dt  */
  ydot[138] = -(s1L + d1_2L + d1_3L + d1_12L + d1_13L + x1L) * D_N1L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * D_N1L
    + d1_3L * D_N3L
    + d1_2L * D_N2L
    + d1_13L * D_N13L
    + d1_12L * D_N12L
    + 2.0 * s1L * D_N1L * E_1L
    + dLAY_A * D_N1A
    + dLAY_B * D_N1B
    + dLAY_C * D_N1C
    + dLAY_D * D_N1D
    + dLAY_E * D_N1E
    + dLAY_F * D_N1F
    + dLAY_G * D_N1G
    + dLAY_H * D_N1H
    + dLAY_I * D_N1I
    + dLAY_J * D_N1J
    + dLAY_K * D_N1K;

  /*  dD_N2L / dt  */
  ydot[139] = -(s2L + d2_1L + d2_3L + d2_12L + d2_23L + x2L) * D_N2L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * D_N2L
    + d2_3L * D_N3L
    + d2_1L * D_N1L
    + d2_23L * D_N23L
    + d2_12L * D_N12L
    + 2.0 * s2L * D_N2L * E_2L
    + dLAY_A * D_N2A
    + dLAY_B * D_N2B
    + dLAY_C * D_N2C
    + dLAY_D * D_N2D
    + dLAY_E * D_N2E
    + dLAY_F * D_N2F
    + dLAY_G * D_N2G
    + dLAY_H * D_N2H
    + dLAY_I * D_N2I
    + dLAY_J * D_N2J
    + dLAY_K * D_N2K;

  /*  dD_N3L / dt  */
  ydot[140] = -(s3L + d3_1L + d3_2L + d3_23L + d3_13L + x3L) * D_N3L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * D_N3L
    + d3_1L * D_N1L
    + d3_2L * D_N2L
    + d3_23L * D_N23L
    + d3_13L * D_N13L
    + 2.0 * s3L * D_N3L * E_3L
    + dLAY_A * D_N3A
    + dLAY_B * D_N3B
    + dLAY_C * D_N3C
    + dLAY_D * D_N3D
    + dLAY_E * D_N3E
    + dLAY_F * D_N3F
    + dLAY_G * D_N3G
    + dLAY_H * D_N3H
    + dLAY_I * D_N3I
    + dLAY_J * D_N3J
    + dLAY_K * D_N3K;

  /*  dD_N12L / dt  */
  ydot[141] = -(s1L + s2L + s12L + d12_1L + d12_2L) * D_N12L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * D_N12L
    + d12_1L * D_N1L
    + d12_2L * D_N2L
    + s12L * (D_N1L * E_2L + D_N2L * E_1L)
    + s1L * (E_1L * D_N12L + E_12L * D_N1L)
    + s2L * (E_2L * D_N12L + E_12L * D_N2L)
    + dLAY_A * D_N12A
    + dLAY_B * D_N12B
    + dLAY_C * D_N12C
    + dLAY_D * D_N12D
    + dLAY_E * D_N12E
    + dLAY_F * D_N12F
    + dLAY_G * D_N12G
    + dLAY_H * D_N12H
    + dLAY_I * D_N12I
    + dLAY_J * D_N12J
    + dLAY_K * D_N12K;

  /*  dD_N23L / dt  */
  ydot[142] = -(s2L + s3L + s23L + d23_2L + d23_3L) * D_N23L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * D_N23L
    + d23_2L * D_N2L
    + d23_3L * D_N3L
    + s23L * (D_N2L * E_3L + D_N3L * E_2L)
    + s2L * (E_2L * D_N23L + E_23L * D_N2L)
    + s3L * (E_3L * D_N23L + E_23L * D_N3L)
    + dLAY_A * D_N23A
    + dLAY_B * D_N23B
    + dLAY_C * D_N23C
    + dLAY_D * D_N23D
    + dLAY_E * D_N23E
    + dLAY_F * D_N23F
    + dLAY_G * D_N23G
    + dLAY_H * D_N23H
    + dLAY_I * D_N23I
    + dLAY_J * D_N23J
    + dLAY_K * D_N23K;

  /*  dD_N13L / dt  */
  ydot[143] = -(s1L + s3L + s13L + d13_1L + d13_3L) * D_N13L
    - (dLAY_A + dLAY_B + dLAY_C + dLAY_D + dLAY_E + dLAY_F + dLAY_G + dLAY_H + dLAY_I + dLAY_J + dLAY_K) * D_N13L
    + d13_1L * D_N1L
    + d13_3L * D_N3L
    + s13L * (D_N1L * E_3L + D_N3L * E_1L)
    + s1L * (E_1L * D_N13L + E_13L * D_N1L)
    + s3L * (E_3L * D_N13L + E_13L * D_N3L)
    + dLAY_A * D_N13A
    + dLAY_B * D_N13B
    + dLAY_C * D_N13C
    + dLAY_D * D_N13D
    + dLAY_E * D_N13E
    + dLAY_F * D_N13F
    + dLAY_G * D_N13G
    + dLAY_H * D_N13H
    + dLAY_I * D_N13I
    + dLAY_J * D_N13J
    + dLAY_K * D_N13K;
}

