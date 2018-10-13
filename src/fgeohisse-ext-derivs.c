
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 380

static double params_fgeohisse[NUMELEMENTS];


void initmod_fgeohisse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_fgeohisse);
}


void fgeohisse_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E00A = y[0];
    double E11A = y[1];
    double E01A = y[2];
    double E00B = y[3];
    double E11B = y[4];
    double E01B = y[5];
    double E00C = y[6];
    double E11C = y[7];
    double E01C = y[8];
    double E00D = y[9];
    double E11D = y[10];
    double E01D = y[11];
    double E00E = y[12];
    double E11E = y[13];
    double E01E = y[14];
    double E00F = y[15];
    double E11F = y[16];
    double E01F = y[17];
    double E00G = y[18];
    double E11G = y[19];
    double E01G = y[20];
    double E00H = y[21];
    double E11H = y[22];
    double E01H = y[23];
    double E00I = y[24];
    double E11I = y[25];
    double E01I = y[26];
    double E00J = y[27];
    double E11J = y[28];
    double E01J = y[29];

    double D00A = y[30];
    double D11A = y[31];
    double D01A = y[32];
    double D00B = y[33];
    double D11B = y[34];
    double D01B = y[35];
    double D00C = y[36];
    double D11C = y[37];
    double D01C = y[38];
    double D00D = y[39];
    double D11D = y[40];
    double D01D = y[41];
    double D00E = y[42];
    double D11E = y[43];
    double D01E = y[44];
    double D00F = y[45];
    double D11F = y[46];
    double D01F = y[47];
    double D00G = y[48];
    double D11G = y[49];
    double D01G = y[50];
    double D00H = y[51];
    double D11H = y[52];
    double D01H = y[53];
    double D00I = y[54];
    double D11I = y[55];
    double D01I = y[56];
    double D00J = y[57];
    double D11J = y[58];
    double D01J = y[59];
    
    double
    s00A = params_fgeohisse[0],
    s11A = params_fgeohisse[1],
    s01A = params_fgeohisse[2],
    x00A = params_fgeohisse[3],
    x11A = params_fgeohisse[4],
    d00A_11A = params_fgeohisse[5],
    d00A_01A = params_fgeohisse[6],
    d11A_00A = params_fgeohisse[7],
    d11A_01A = params_fgeohisse[8],
    d01A_00A = params_fgeohisse[9],
    d01A_11A = params_fgeohisse[10],
    
    d00A_00B = params_fgeohisse[11],
    d00A_00C = params_fgeohisse[12],
    d00A_00D = params_fgeohisse[13],
    d00A_00E = params_fgeohisse[14],
    d00A_00F = params_fgeohisse[15],
    d00A_00G = params_fgeohisse[16],
    d00A_00H = params_fgeohisse[17],
    d00A_00I = params_fgeohisse[18],
    d00A_00J = params_fgeohisse[19],
    d11A_11B = params_fgeohisse[20],
    d11A_11C = params_fgeohisse[21],
    d11A_11D = params_fgeohisse[22],
    d11A_11E = params_fgeohisse[23],
    d11A_11F = params_fgeohisse[24],
    d11A_11G = params_fgeohisse[25],
    d11A_11H = params_fgeohisse[26],
    d11A_11I = params_fgeohisse[27],
    d11A_11J = params_fgeohisse[28],
    d01A_01B = params_fgeohisse[29],
    d01A_01C = params_fgeohisse[30],
    d01A_01D = params_fgeohisse[31],
    d01A_01E = params_fgeohisse[32],
    d01A_01F = params_fgeohisse[33],
    d01A_01G = params_fgeohisse[34],
    d01A_01H = params_fgeohisse[35],
    d01A_01I = params_fgeohisse[36],
    d01A_01J = params_fgeohisse[37],

    s00B = params_fgeohisse[38],
    s11B = params_fgeohisse[39],
    s01B = params_fgeohisse[40],
    x00B = params_fgeohisse[41],
    x11B = params_fgeohisse[42],
    d00B_11B = params_fgeohisse[43],
    d00B_01B = params_fgeohisse[44],
    d11B_00B = params_fgeohisse[45],
    d11B_01B = params_fgeohisse[46],
    d01B_00B = params_fgeohisse[47],
    d01B_11B = params_fgeohisse[48],
    
    d00B_00A = params_fgeohisse[49],
    d00B_00C = params_fgeohisse[50],
    d00B_00D = params_fgeohisse[51],
    d00B_00E = params_fgeohisse[52],
    d00B_00F = params_fgeohisse[53],
    d00B_00G = params_fgeohisse[54],
    d00B_00H = params_fgeohisse[55],
    d00B_00I = params_fgeohisse[56],
    d00B_00J = params_fgeohisse[57],
    d11B_11A = params_fgeohisse[58],
    d11B_11C = params_fgeohisse[59],
    d11B_11D = params_fgeohisse[60],
    d11B_11E = params_fgeohisse[61],
    d11B_11F = params_fgeohisse[62],
    d11B_11G = params_fgeohisse[63],
    d11B_11H = params_fgeohisse[64],
    d11B_11I = params_fgeohisse[65],
    d11B_11J = params_fgeohisse[66],
    d01B_01A = params_fgeohisse[67],
    d01B_01C = params_fgeohisse[68],
    d01B_01D = params_fgeohisse[69],
    d01B_01E = params_fgeohisse[70],
    d01B_01F = params_fgeohisse[71],
    d01B_01G = params_fgeohisse[72],
    d01B_01H = params_fgeohisse[73],
    d01B_01I = params_fgeohisse[74],
    d01B_01J = params_fgeohisse[75],

    s00C = params_fgeohisse[76],
    s11C = params_fgeohisse[77],
    s01C = params_fgeohisse[78],
    x00C = params_fgeohisse[79],
    x11C = params_fgeohisse[80],
    d00C_11C = params_fgeohisse[81],
    d00C_01C = params_fgeohisse[82],
    d11C_00C = params_fgeohisse[83],
    d11C_01C = params_fgeohisse[84],
    d01C_00C = params_fgeohisse[85],
    d01C_11C = params_fgeohisse[86],
    
    d00C_00A = params_fgeohisse[87],
    d00C_00B = params_fgeohisse[88],
    d00C_00D = params_fgeohisse[89],
    d00C_00E = params_fgeohisse[90],
    d00C_00F = params_fgeohisse[91],
    d00C_00G = params_fgeohisse[92],
    d00C_00H = params_fgeohisse[93],
    d00C_00I = params_fgeohisse[94],
    d00C_00J = params_fgeohisse[95],
    d11C_11A = params_fgeohisse[96],
    d11C_11B = params_fgeohisse[97],
    d11C_11D = params_fgeohisse[98],
    d11C_11E = params_fgeohisse[99],
    d11C_11F = params_fgeohisse[100],
    d11C_11G = params_fgeohisse[101],
    d11C_11H = params_fgeohisse[102],
    d11C_11I = params_fgeohisse[103],
    d11C_11J = params_fgeohisse[104],
    d01C_01A = params_fgeohisse[105],
    d01C_01B = params_fgeohisse[106],
    d01C_01D = params_fgeohisse[107],
    d01C_01E = params_fgeohisse[108],
    d01C_01F = params_fgeohisse[109],
    d01C_01G = params_fgeohisse[110],
    d01C_01H = params_fgeohisse[111],
    d01C_01I = params_fgeohisse[112],
    d01C_01J = params_fgeohisse[113],

    s00D = params_fgeohisse[114],
    s11D = params_fgeohisse[115],
    s01D = params_fgeohisse[116],
    x00D = params_fgeohisse[117],
    x11D = params_fgeohisse[118],
    d00D_11D = params_fgeohisse[119],
    d00D_01D = params_fgeohisse[120],
    d11D_00D = params_fgeohisse[121],
    d11D_01D = params_fgeohisse[122],
    d01D_00D = params_fgeohisse[123],
    d01D_11D = params_fgeohisse[124],
    
    d00D_00A = params_fgeohisse[125],
    d00D_00B = params_fgeohisse[126],
    d00D_00C = params_fgeohisse[127],
    d00D_00E = params_fgeohisse[128],
    d00D_00F = params_fgeohisse[129],
    d00D_00G = params_fgeohisse[130],
    d00D_00H = params_fgeohisse[131],
    d00D_00I = params_fgeohisse[132],
    d00D_00J = params_fgeohisse[133],
    d11D_11A = params_fgeohisse[134],
    d11D_11B = params_fgeohisse[135],
    d11D_11C = params_fgeohisse[136],
    d11D_11E = params_fgeohisse[137],
    d11D_11F = params_fgeohisse[138],
    d11D_11G = params_fgeohisse[139],
    d11D_11H = params_fgeohisse[140],
    d11D_11I = params_fgeohisse[141],
    d11D_11J = params_fgeohisse[142],
    d01D_01A = params_fgeohisse[143],
    d01D_01B = params_fgeohisse[144],
    d01D_01C = params_fgeohisse[145],
    d01D_01E = params_fgeohisse[146],
    d01D_01F = params_fgeohisse[147],
    d01D_01G = params_fgeohisse[148],
    d01D_01H = params_fgeohisse[149],
    d01D_01I = params_fgeohisse[150],
    d01D_01J = params_fgeohisse[151],

    s00E = params_fgeohisse[152],
    s11E = params_fgeohisse[153],
    s01E = params_fgeohisse[154],
    x00E = params_fgeohisse[155],
    x11E = params_fgeohisse[156],
    d00E_11E = params_fgeohisse[157],
    d00E_01E = params_fgeohisse[158],
    d11E_00E = params_fgeohisse[159],
    d11E_01E = params_fgeohisse[160],
    d01E_00E = params_fgeohisse[161],
    d01E_11E = params_fgeohisse[162],
    
    d00E_00A = params_fgeohisse[163],
    d00E_00B = params_fgeohisse[164],
    d00E_00C = params_fgeohisse[165],
    d00E_00D = params_fgeohisse[166],
    d00E_00F = params_fgeohisse[167],
    d00E_00G = params_fgeohisse[168],
    d00E_00H = params_fgeohisse[169],
    d00E_00I = params_fgeohisse[170],
    d00E_00J = params_fgeohisse[171],
    d11E_11A = params_fgeohisse[172],
    d11E_11B = params_fgeohisse[173],
    d11E_11C = params_fgeohisse[174],
    d11E_11D = params_fgeohisse[175],
    d11E_11F = params_fgeohisse[176],
    d11E_11G = params_fgeohisse[177],
    d11E_11H = params_fgeohisse[178],
    d11E_11I = params_fgeohisse[179],
    d11E_11J = params_fgeohisse[180],
    d01E_01A = params_fgeohisse[181],
    d01E_01B = params_fgeohisse[182],
    d01E_01C = params_fgeohisse[183],
    d01E_01D = params_fgeohisse[184],
    d01E_01F = params_fgeohisse[185],
    d01E_01G = params_fgeohisse[186],
    d01E_01H = params_fgeohisse[187],
    d01E_01I = params_fgeohisse[188],
    d01E_01J = params_fgeohisse[189],

    s00F = params_fgeohisse[190],
    s11F = params_fgeohisse[191],
    s01F = params_fgeohisse[192],
    x00F = params_fgeohisse[193],
    x11F = params_fgeohisse[194],
    d00F_11F = params_fgeohisse[195],
    d00F_01F = params_fgeohisse[196],
    d11F_00F = params_fgeohisse[197],
    d11F_01F = params_fgeohisse[198],
    d01F_00F = params_fgeohisse[199],
    d01F_11F = params_fgeohisse[200],
    
    d00F_00A = params_fgeohisse[201],
    d00F_00B = params_fgeohisse[202],
    d00F_00C = params_fgeohisse[203],
    d00F_00D = params_fgeohisse[204],
    d00F_00E = params_fgeohisse[205],
    d00F_00G = params_fgeohisse[206],
    d00F_00H = params_fgeohisse[207],
    d00F_00I = params_fgeohisse[208],
    d00F_00J = params_fgeohisse[209],
    d11F_11A = params_fgeohisse[210],
    d11F_11B = params_fgeohisse[211],
    d11F_11C = params_fgeohisse[212],
    d11F_11D = params_fgeohisse[213],
    d11F_11E = params_fgeohisse[214],
    d11F_11G = params_fgeohisse[215],
    d11F_11H = params_fgeohisse[216],
    d11F_11I = params_fgeohisse[217],
    d11F_11J = params_fgeohisse[218],
    d01F_01A = params_fgeohisse[219],
    d01F_01B = params_fgeohisse[220],
    d01F_01C = params_fgeohisse[221],
    d01F_01D = params_fgeohisse[222],
    d01F_01E = params_fgeohisse[223],
    d01F_01G = params_fgeohisse[224],
    d01F_01H = params_fgeohisse[225],
    d01F_01I = params_fgeohisse[226],
    d01F_01J = params_fgeohisse[227],

    s00G = params_fgeohisse[228],
    s11G = params_fgeohisse[229],
    s01G = params_fgeohisse[230],
    x00G = params_fgeohisse[231],
    x11G = params_fgeohisse[232],
    d00G_11G = params_fgeohisse[233],
    d00G_01G = params_fgeohisse[234],
    d11G_00G = params_fgeohisse[235],
    d11G_01G = params_fgeohisse[236],
    d01G_00G = params_fgeohisse[237],
    d01G_11G = params_fgeohisse[238],
    
    d00G_00A = params_fgeohisse[239],
    d00G_00B = params_fgeohisse[240],
    d00G_00C = params_fgeohisse[241],
    d00G_00D = params_fgeohisse[242],
    d00G_00E = params_fgeohisse[243],
    d00G_00F = params_fgeohisse[244],
    d00G_00H = params_fgeohisse[245],
    d00G_00I = params_fgeohisse[246],
    d00G_00J = params_fgeohisse[247],
    d11G_11A = params_fgeohisse[248],
    d11G_11B = params_fgeohisse[249],
    d11G_11C = params_fgeohisse[250],
    d11G_11D = params_fgeohisse[251],
    d11G_11E = params_fgeohisse[252],
    d11G_11F = params_fgeohisse[253],
    d11G_11H = params_fgeohisse[254],
    d11G_11I = params_fgeohisse[255],
    d11G_11J = params_fgeohisse[256],
    d01G_01A = params_fgeohisse[257],
    d01G_01B = params_fgeohisse[258],
    d01G_01C = params_fgeohisse[259],
    d01G_01D = params_fgeohisse[260],
    d01G_01E = params_fgeohisse[261],
    d01G_01F = params_fgeohisse[262],
    d01G_01H = params_fgeohisse[263],
    d01G_01I = params_fgeohisse[264],
    d01G_01J = params_fgeohisse[265],
    
    s00H = params_fgeohisse[266],
    s11H = params_fgeohisse[267],
    s01H = params_fgeohisse[268],
    x00H = params_fgeohisse[269],
    x11H = params_fgeohisse[270],
    d00H_11H = params_fgeohisse[271],
    d00H_01H = params_fgeohisse[272],
    d11H_00H = params_fgeohisse[273],
    d11H_01H = params_fgeohisse[274],
    d01H_00H = params_fgeohisse[275],
    d01H_11H = params_fgeohisse[276],
    
    d00H_00A = params_fgeohisse[277],
    d00H_00B = params_fgeohisse[278],
    d00H_00C = params_fgeohisse[279],
    d00H_00D = params_fgeohisse[280],
    d00H_00E = params_fgeohisse[281],
    d00H_00F = params_fgeohisse[282],
    d00H_00G = params_fgeohisse[283],
    d00H_00I = params_fgeohisse[284],
    d00H_00J = params_fgeohisse[285],
    d11H_11A = params_fgeohisse[286],
    d11H_11B = params_fgeohisse[287],
    d11H_11C = params_fgeohisse[288],
    d11H_11D = params_fgeohisse[289],
    d11H_11E = params_fgeohisse[290],
    d11H_11F = params_fgeohisse[291],
    d11H_11G = params_fgeohisse[292],
    d11H_11I = params_fgeohisse[293],
    d11H_11J = params_fgeohisse[294],
    d01H_01A = params_fgeohisse[295],
    d01H_01B = params_fgeohisse[296],
    d01H_01C = params_fgeohisse[297],
    d01H_01D = params_fgeohisse[298],
    d01H_01E = params_fgeohisse[299],
    d01H_01F = params_fgeohisse[300],
    d01H_01G = params_fgeohisse[301],
    d01H_01I = params_fgeohisse[302],
    d01H_01J = params_fgeohisse[303],

    s00I = params_fgeohisse[304],
    s11I = params_fgeohisse[305],
    s01I = params_fgeohisse[306],
    x00I = params_fgeohisse[307],
    x11I = params_fgeohisse[308],
    d00I_11I = params_fgeohisse[309],
    d00I_01I = params_fgeohisse[310],
    d11I_00I = params_fgeohisse[311],
    d11I_01I = params_fgeohisse[312],
    d01I_00I = params_fgeohisse[313],
    d01I_11I = params_fgeohisse[314],
    
    d00I_00A = params_fgeohisse[315],
    d00I_00B = params_fgeohisse[316],
    d00I_00C = params_fgeohisse[317],
    d00I_00D = params_fgeohisse[318],
    d00I_00E = params_fgeohisse[319],
    d00I_00F = params_fgeohisse[320],
    d00I_00G = params_fgeohisse[321],
    d00I_00H = params_fgeohisse[322],
    d00I_00J = params_fgeohisse[323],
    d11I_11A = params_fgeohisse[324],
    d11I_11B = params_fgeohisse[325],
    d11I_11C = params_fgeohisse[326],
    d11I_11D = params_fgeohisse[327],
    d11I_11E = params_fgeohisse[328],
    d11I_11F = params_fgeohisse[329],
    d11I_11G = params_fgeohisse[330],
    d11I_11H = params_fgeohisse[331],
    d11I_11J = params_fgeohisse[332],
    d01I_01A = params_fgeohisse[333],
    d01I_01B = params_fgeohisse[334],
    d01I_01C = params_fgeohisse[335],
    d01I_01D = params_fgeohisse[336],
    d01I_01E = params_fgeohisse[337],
    d01I_01F = params_fgeohisse[338],
    d01I_01G = params_fgeohisse[339],
    d01I_01H = params_fgeohisse[340],
    d01I_01J = params_fgeohisse[341],

    s00J = params_fgeohisse[342],
    s11J = params_fgeohisse[343],
    s01J = params_fgeohisse[344],
    x00J = params_fgeohisse[345],
    x11J = params_fgeohisse[346],
    d00J_11J = params_fgeohisse[347],
    d00J_01J = params_fgeohisse[348],
    d11J_00J = params_fgeohisse[349],
    d11J_01J = params_fgeohisse[350],
    d01J_00J = params_fgeohisse[351],
    d01J_11J = params_fgeohisse[352],
    
    d00J_00A = params_fgeohisse[353],
    d00J_00B = params_fgeohisse[354],
    d00J_00C = params_fgeohisse[355],
    d00J_00D = params_fgeohisse[356],
    d00J_00E = params_fgeohisse[357],
    d00J_00F = params_fgeohisse[358],
    d00J_00G = params_fgeohisse[359],
    d00J_00H = params_fgeohisse[360],
    d00J_00I = params_fgeohisse[361],
    d11J_11A = params_fgeohisse[362],
    d11J_11B = params_fgeohisse[363],
    d11J_11C = params_fgeohisse[364],
    d11J_11D = params_fgeohisse[365],
    d11J_11E = params_fgeohisse[366],
    d11J_11F = params_fgeohisse[367],
    d11J_11G = params_fgeohisse[368],
    d11J_11H = params_fgeohisse[369],
    d11J_11I = params_fgeohisse[370],
    d01J_01A = params_fgeohisse[371],
    d01J_01B = params_fgeohisse[372],
    d01J_01C = params_fgeohisse[373],
    d01J_01D = params_fgeohisse[374],
    d01J_01E = params_fgeohisse[375],
    d01J_01F = params_fgeohisse[376],
    d01J_01G = params_fgeohisse[377],
    d01J_01H = params_fgeohisse[378],
    d01J_01I = params_fgeohisse[379];


    /* The E's */
    
    ydot[0] =  -(s00A+(d00A_11A + d00A_01A + d00A_00B + d00A_00C + d00A_00D + d00A_00E + d00A_00F + d00A_00G + d00A_00H + d00A_00I + d00A_00J)+x00A) * E00A + (d00A_11A*E11A + d00A_01A*E01A + d00A_00B*E00B + d00A_00C*E00C + d00A_00D*E00D + d00A_00E*E00E + d00A_00F*E00F + d00A_00G*E00G + d00A_00H*E00H + d00A_00I*E00I + d00A_00J*E00J) + x00A+s00A*E00A*E00A;
    
    ydot[1] =  -(s11A+(d11A_00A + d11A_01A + d11A_11B + d11A_11C + d11A_11D + d11A_11E + d11A_11F + d11A_11G + d11A_11H + d11A_11I + d11A_11J)+x11A) * E11A + (d11A_00A*E00A + d11A_01A*E01A + d11A_11B*E11B + d11A_11C*E11C + d11A_11D*E11D + d11A_11E*E11E + d11A_11F*E11F + d11A_11G*E11G + d11A_11H*E11H + d11A_11I*E11I + d11A_11J*E11J) + x11A+s11A*E11A*E11A;
    
    ydot[2] =  -(s01A + s00A + s11A + (d01A_00A + d01A_11A + d01A_01B + d01A_01C + d01A_01D + d01A_01E + d01A_01F + d01A_01G + d01A_01H + d01A_01I + d01A_01J)) * E01A + (d01A_00A*E00A + d01A_11A*E11A + d01A_01B*E01B + d01A_01C*E01C + d01A_01D*E01D + d01A_01E*E01E + d01A_01F*E01F + d01A_01G*E01G + d01A_01H*E01H + d01A_01I*E01I + d01A_01J*E01J) + s00A*E01A*E00A + s11A*E01A*E11A + s01A*E11A*E00A;
    
    ydot[3] =  -(s00B+(d00B_00A + d00B_11B + d00B_01B + d00B_00C + d00B_00D + d00B_00E + d00B_00F + d00B_00G + d00B_00H + d00B_00I + d00B_00J)+x00B) * E00B + (d00B_00A*E00A + d00B_11B*E11B + d00B_01B*E01B + d00B_00C*E00C + d00B_00D*E00D + d00B_00E*E00E + d00B_00F*E00F + d00B_00G*E00G + d00B_00H*E00H + d00B_00I*E00I + d00B_00J*E00J) + x00B+s00B*E00B*E00B;
    
    ydot[4] =  -(s11B+(d11B_11A + d11B_00B + d11B_01B + d11B_11C + d11B_11D + d11B_11E + d11B_11F + d11B_11G + d11B_11H + d11B_11I + d11B_11J)+x11B) * E11B + (d11B_11A*E11A + d11B_00B*E00B + d11B_01B*E01B + d11B_11C*E11C + d11B_11D*E11D + d11B_11E*E11E + d11B_11F*E11F + d11B_11G*E11G + d11B_11H*E11H + d11B_11I*E11I + d11B_11J*E11J) + x11B+s11B*E11B*E11B;
    
    ydot[5] =  -(s01B + s00B + s11B + (d01B_01A + d01B_00B + d01B_11B + d01B_01C + d01B_01D + d01B_01E + d01B_01F + d01B_01G + d01B_01H + d01B_01I + d01B_01J)) * E01B + (d01B_01A*E01A + d01B_00B*E00B + d01B_11B*E11B + d01B_01C*E01C + d01B_01D*E01D + d01B_01E*E01E + d01B_01F*E01F + d01B_01G*E01G + d01B_01H*E01H + d01B_01I*E01I + d01B_01J*E01J) + s00B*E01B*E00B + s11B*E01B*E11B + s01B*E11B*E00B;
    
    ydot[6] =  -(s00C+(d00C_00A + d00C_00B + d00C_11C + d00C_01C + d00C_00D + d00C_00E + d00C_00F + d00C_00G + d00C_00H + d00C_00I + d00C_00J)+x00C) * E00C + (d00C_00A*E00A + d00C_00B*E00B + d00C_11C*E11C + d00C_01C*E01C + d00C_00D*E00D + d00C_00E*E00E + d00C_00F*E00F + d00C_00G*E00G + d00C_00H*E00H + d00C_00I*E00I + d00C_00J*E00J) + x00C+s00C*E00C*E00C;
    
    ydot[7] =  -(s11C+(d11C_11A + d11C_11B + d11C_00C + d11C_01C + d11C_11D + d11C_11E + d11C_11F + d11C_11G + d11C_11H + d11C_11I + d11C_11J)+x11C) * E11C + (d11C_11A*E11A + d11C_11B*E11B + d11C_00C*E00C + d11C_01C*E01C + d11C_11D*E11D + d11C_11E*E11E + d11C_11F*E11F + d11C_11G*E11G + d11C_11H*E11H + d11C_11I*E11I + d11C_11J*E11J) + x11C+s11C*E11C*E11C;
    
    ydot[8] =  -(s01C + s00C + s11C + (d01C_01A + d01C_01B + d01C_00C + d01C_11C + d01C_01D + d01C_01E + d01C_01F + d01C_01G + d01C_01H + d01C_01I + d01C_01J)) * E01C + (d01C_01A*E01A + d01C_01B*E01B + d01C_00C*E00C + d01C_11C*E11C + d01C_01D*E01D + d01C_01E*E01E + d01C_01F*E01F + d01C_01G*E01G + d01C_01H*E01H + d01C_01I*E01I + d01C_01J*E01J) + s00C*E01C*E00C + s11C*E01C*E11C + s01C*E11C*E00C;
    
    ydot[9] =  -(s00D+(d00D_00A + d00D_00B + d00D_00C + d00D_11D + d00D_01D + d00D_00E + d00D_00F + d00D_00G + d00D_00H + d00D_00I + d00D_00J)+x00D) * E00D + (d00D_00A*E00A + d00D_00B*E00B + d00D_00C*E00C + d00D_11D*E11D + d00D_01D*E01D + d00D_00E*E00E + d00D_00F*E00F + d00D_00G*E00G + d00D_00H*E00H + d00D_00I*E00I + d00D_00J*E00J) + x00D+s00D*E00D*E00D;
    
    ydot[10] =  -(s11D+(d11D_11A + d11D_11B + d11D_11C + d11D_00D + d11D_01D + d11D_11E + d11D_11F + d11D_11G + d11D_11H + d11D_11I + d11D_11J)+x11D) * E11D + (d11D_11A*E11A + d11D_11B*E11B + d11D_11C*E11C + d11D_00D*E00D + d11D_01D*E01D + d11D_11E*E11E + d11D_11F*E11F + d11D_11G*E11G + d11D_11H*E11H + d11D_11I*E11I + d11D_11J*E11J) + x11D+s11D*E11D*E11D;
    
    ydot[11] =  -(s01D + s00D + s11D + (d01D_01A + d01D_01B + d01D_01C + d01D_00D + d01D_11D + d01D_01E + d01D_01F + d01D_01G + d01D_01H + d01D_01I + d01D_01J)) * E01D + (d01D_01A*E01A + d01D_01B*E01B + d01D_01C*E01C + d01D_00D*E00D + d01D_11D*E11D + d01D_01E*E01E + d01D_01F*E01F + d01D_01G*E01G + d01D_01H*E01H + d01D_01I*E01I + d01D_01J*E01J) + s00D*E01D*E00D + s11D*E01D*E11D + s01D*E11D*E00D;
    
    ydot[12] =  -(s00E+(d00E_00A + d00E_00B + d00E_00C + d00E_00D + d00E_11E + d00E_01E + d00E_00F + d00E_00G + d00E_00H + d00E_00I + d00E_00J)+x00E) * E00E + (d00E_00A*E00A + d00E_00B*E00B + d00E_00C*E00C + d00E_00D*E00D + d00E_11E*E11E + d00E_01E*E01E + d00E_00F*E00F + d00E_00G*E00G + d00E_00H*E00H + d00E_00I*E00I + d00E_00J*E00J) + x00E+s00E*E00E*E00E;
    
    ydot[13] =  -(s11E+(d11E_11A + d11E_11B + d11E_11C + d11E_11D + d11E_00E + d11E_01E + d11E_11F + d11E_11G + d11E_11H + d11E_11I + d11E_11J)+x11E) * E11E + (d11E_11A*E11A + d11E_11B*E11B + d11E_11C*E11C + d11E_11D*E11D + d11E_00E*E00E + d11E_01E*E01E + d11E_11F*E11F + d11E_11G*E11G + d11E_11H*E11H + d11E_11I*E11I + d11E_11J*E11J) + x11E+s11E*E11E*E11E;
    
    ydot[14] =  -(s01E + s00E + s11E + (d01E_01A + d01E_01B + d01E_01C + d01E_01D + d01E_00E + d01E_11E + d01E_01F + d01E_01G + d01E_01H + d01E_01I + d01E_01J)) * E01E + (d01E_01A*E01A + d01E_01B*E01B + d01E_01C*E01C + d01E_01D*E01D + d01E_00E*E00E + d01E_11E*E11E + d01E_01F*E01F + d01E_01G*E01G + d01E_01H*E01H + d01E_01I*E01I + d01E_01J*E01J) + s00E*E01E*E00E + s11E*E01E*E11E + s01E*E11E*E00E;
    
    ydot[15] =  -(s00F+(d00F_00A + d00F_00B + d00F_00C + d00F_00D + d00F_00E + d00F_11F + d00F_01F + d00F_00G + d00F_00H + d00F_00I + d00F_00J)+x00F) * E00F + (d00F_00A*E00A + d00F_00B*E00B + d00F_00C*E00C + d00F_00D*E00D + d00F_00E*E00E + d00F_11F*E11F + d00F_01F*E01F + d00F_00G*E00G + d00F_00H*E00H + d00F_00I*E00I + d00F_00J*E00J) + x00F+s00F*E00F*E00F;
    
    ydot[16] =  -(s11F+(d11F_11A + d11F_11B + d11F_11C + d11F_11D + d11F_11E + d11F_00F + d11F_01F + d11F_11G + d11F_11H + d11F_11I + d11F_11J)+x11F) * E11F + (d11F_11A*E11A + d11F_11B*E11B + d11F_11C*E11C + d11F_11D*E11D + d11F_11E*E11E + d11F_00F*E00F + d11F_01F*E01F + d11F_11G*E11G + d11F_11H*E11H + d11F_11I*E11I + d11F_11J*E11J) + x11F+s11F*E11F*E11F;
    
    ydot[17] =  -(s01F + s00F + s11F + (d01F_01A + d01F_01B + d01F_01C + d01F_01D + d01F_01E + d01F_00F + d01F_11F + d01F_01G + d01F_01H + d01F_01I + d01F_01J)) * E01F + (d01F_01A*E01A + d01F_01B*E01B + d01F_01C*E01C + d01F_01D*E01D + d01F_01E*E01E + d01F_00F*E00F + d01F_11F*E11F + d01F_01G*E01G + d01F_01H*E01H + d01F_01I*E01I + d01F_01J*E01J) + s00F*E01F*E00F + s11F*E01F*E11F + s01F*E11F*E00F;
    
    ydot[18] =  -(s00G+(d00G_00A + d00G_00B + d00G_00C + d00G_00D + d00G_00E + d00G_00F + d00G_11G + d00G_01G + d00G_00H + d00G_00I + d00G_00J)+x00G) * E00G + (d00G_00A*E00A + d00G_00B*E00B + d00G_00C*E00C + d00G_00D*E00D + d00G_00E*E00E + d00G_00F*E00F + d00G_11G*E11G + d00G_01G*E01G + d00G_00H*E00H + d00G_00I*E00I + d00G_00J*E00J) + x00G+s00G*E00G*E00G;
    
    ydot[19] =  -(s11G+(d11G_11A + d11G_11B + d11G_11C + d11G_11D + d11G_11E + d11G_11F + d11G_00G + d11G_01G + d11G_11H + d11G_11I + d11G_11J)+x11G) * E11G + (d11G_11A*E11A + d11G_11B*E11B + d11G_11C*E11C + d11G_11D*E11D + d11G_11E*E11E + d11G_11F*E11F + d11G_00G*E00G + d11G_01G*E01G + d11G_11H*E11H + d11G_11I*E11I + d11G_11J*E11J) + x11G+s11G*E11G*E11G;
    
    ydot[20] =  -(s01G + s00G + s11G + (d01G_01A + d01G_01B + d01G_01C + d01G_01D + d01G_01E + d01G_01F + d01G_00G + d01G_11G + d01G_01H + d01G_01I + d01G_01J)) * E01G + (d01G_01A*E01A + d01G_01B*E01B + d01G_01C*E01C + d01G_01D*E01D + d01G_01E*E01E + d01G_01F*E01F + d01G_00G*E00G + d01G_11G*E11G + d01G_01H*E01H + d01G_01I*E01I + d01G_01J*E01J) + s00G*E01G*E00G + s11G*E01G*E11G + s01G*E11G*E00G;
    
    ydot[21] =  -(s00H+(d00H_00A + d00H_00B + d00H_00C + d00H_00D + d00H_00E + d00H_00F + d00H_00G + d00H_11H + d00H_01H + d00H_00I + d00H_00J)+x00H) * E00H + (d00H_00A*E00A + d00H_00B*E00B + d00H_00C*E00C + d00H_00D*E00D + d00H_00E*E00E + d00H_00F*E00F + d00H_00G*E00G + d00H_11H*E11H + d00H_01H*E01H + d00H_00I*E00I + d00H_00J*E00J) + x00H+s00H*E00H*E00H;
    
    ydot[22] =  -(s11H+(d11H_11A + d11H_11B + d11H_11C + d11H_11D + d11H_11E + d11H_11F + d11H_11G + d11H_00H + d11H_01H + d11H_11I + d11H_11J)+x11H) * E11H + (d11H_11A*E11A + d11H_11B*E11B + d11H_11C*E11C + d11H_11D*E11D + d11H_11E*E11E + d11H_11F*E11F + d11H_11G*E11G + d11H_00H*E00H + d11H_01H*E01H + d11H_11I*E11I + d11H_11J*E11J) + x11H+s11H*E11H*E11H;
    
    ydot[23] =  -(s01H + s00H + s11H + (d01H_01A + d01H_01B + d01H_01C + d01H_01D + d01H_01E + d01H_01F + d01H_01G + d01H_00H + d01H_11H + d01H_01I + d01H_01J)) * E01H + (d01H_01A*E01A + d01H_01B*E01B + d01H_01C*E01C + d01H_01D*E01D + d01H_01E*E01E + d01H_01F*E01F + d01H_01G*E01G + d01H_00H*E00H + d01H_11H*E11H + d01H_01I*E01I + d01H_01J*E01J) + s00H*E01H*E00H + s11H*E01H*E11H + s01H*E11H*E00H;
    
    ydot[24] =  -(s00I+(d00I_00A + d00I_00B + d00I_00C + d00I_00D + d00I_00E + d00I_00F + d00I_00G + d00I_00H + d00I_11I + d00I_01I + d00I_00J)+x00I) * E00I + (d00I_00A*E00A + d00I_00B*E00B + d00I_00C*E00C + d00I_00D*E00D + d00I_00E*E00E + d00I_00F*E00F + d00I_00G*E00G + d00I_00H*E00H + d00I_11I*E11I + d00I_01I*E01I + d00I_00J*E00J) + x00I+s00I*E00I*E00I;
    
    ydot[25] =  -(s11I+(d11I_11A + d11I_11B + d11I_11C + d11I_11D + d11I_11E + d11I_11F + d11I_11G + d11I_11H + d11I_00I + d11I_01I + d11I_11J)+x11I) * E11I + (d11I_11A*E11A + d11I_11B*E11B + d11I_11C*E11C + d11I_11D*E11D + d11I_11E*E11E + d11I_11F*E11F + d11I_11G*E11G + d11I_11H*E11H + d11I_00I*E00I + d11I_01I*E01I + d11I_11J*E11J) + x11I+s11I*E11I*E11I;
    
    ydot[26] =  -(s01I + s00I + s11I + (d01I_01A + d01I_01B + d01I_01C + d01I_01D + d01I_01E + d01I_01F + d01I_01G + d01I_01H + d01I_00I + d01I_11I + d01I_01J)) * E01I + (d01I_01A*E01A + d01I_01B*E01B + d01I_01C*E01C + d01I_01D*E01D + d01I_01E*E01E + d01I_01F*E01F + d01I_01G*E01G + d01I_01H*E01H + d01I_00I*E00I + d01I_11I*E11I + d01I_01J*E01J) + s00I*E01I*E00I + s11I*E01I*E11I + s01I*E11I*E00I;
    
    ydot[27] =  -(s00J+(d00J_00A + d00J_00B + d00J_00C + d00J_00D + d00J_00E + d00J_00F + d00J_00G + d00J_00H + d00J_00I + d00J_11J + d00J_01J)+x00J) * E00J + (d00J_00A*E00A + d00J_00B*E00B + d00J_00C*E00C + d00J_00D*E00D + d00J_00E*E00E + d00J_00F*E00F + d00J_00G*E00G + d00J_00H*E00H + d00J_00I*E00I + d00J_11J*E11J + d00J_01J*E01J) + x00J+s00J*E00J*E00J;
    
    ydot[28] =  -(s11J+(d11J_11A + d11J_11B + d11J_11C + d11J_11D + d11J_11E + d11J_11F + d11J_11G + d11J_11H + d11J_11I + d11J_00J + d11J_01J)+x11J) * E11J + (d11J_11A*E11A + d11J_11B*E11B + d11J_11C*E11C + d11J_11D*E11D + d11J_11E*E11E + d11J_11F*E11F + d11J_11G*E11G + d11J_11H*E11H + d11J_11I*E11I + d11J_00J*E00J + d11J_01J*E01J) + x11J+s11J*E11J*E11J;
    
    ydot[29] =  -(s01J + s00J + s11J + (d01J_01A + d01J_01B + d01J_01C + d01J_01D + d01J_01E + d01J_01F + d01J_01G + d01J_01H + d01J_01I + d01J_00J + d01J_11J)) * E01J + (d01J_01A*E01A + d01J_01B*E01B + d01J_01C*E01C + d01J_01D*E01D + d01J_01E*E01E + d01J_01F*E01F + d01J_01G*E01G + d01J_01H*E01H + d01J_01I*E01I + d01J_00J*E00J + d01J_11J*E11J) + s00J*E01J*E00J + s11J*E01J*E11J + s01J*E11J*E00J;
    
    
    
    
    /* The D's */
    ydot[30] =  -(s00A+(d00A_11A + d00A_01A + d00A_00B + d00A_00C + d00A_00D + d00A_00E + d00A_00F + d00A_00G + d00A_00H + d00A_00I + d00A_00J)+x00A) * D00A + (d00A_11A*D11A + d00A_01A*D01A + d00A_00B*D00B + d00A_00C*D00C + d00A_00D*D00D + d00A_00E*D00E + d00A_00F*D00F + d00A_00G*D00G + d00A_00H*D00H + d00A_00I*D00I + d00A_00J*D00J) + 2*s00A*E00A*D00A;
    
    ydot[31] =  -(s11A+(d11A_00A + d11A_01A + d11A_11B + d11A_11C + d11A_11D + d11A_11E + d11A_11F + d11A_11G + d11A_11H + d11A_11I + d11A_11J)+x11A) * D11A + (d11A_00A*D00A + d11A_01A*D01A + d11A_11B*D11B + d11A_11C*D11C + d11A_11D*D11D + d11A_11E*D11E + d11A_11F*D11F + d11A_11G*D11G + d11A_11H*D11H + d11A_11I*D11I + d11A_11J*D11J) + 2*s11A*E11A*D11A;
    
    ydot[32] =  -(s01A + s00A + s11A + (d01A_00A + d01A_11A + d01A_01B + d01A_01C + d01A_01D + d01A_01E + d01A_01F + d01A_01G + d01A_01H + d01A_01I + d01A_01J)) * D01A + (d01A_00A*D00A + d01A_11A*D11A + d01A_01B*D01B + d01A_01C*D01C + d01A_01D*D01D + d01A_01E*D01E + d01A_01F*D01F + d01A_01G*D01G + d01A_01H*D01H + d01A_01I*D01I + d01A_01J*D01J) + s00A*(E00A*D01A+E01A*D00A) + s11A*(E11A*D01A + E01A*D11A) + s01A*(E00A*D11A + E11A*D00A);
    
    ydot[33] =  -(s00B+(d00B_00A + d00B_11B + d00B_01B + d00B_00C + d00B_00D + d00B_00E + d00B_00F + d00B_00G + d00B_00H + d00B_00I + d00B_00J)+x00B) * D00B + (d00B_00A*D00A + d00B_11B*D11B + d00B_01B*D01B + d00B_00C*D00C + d00B_00D*D00D + d00B_00E*D00E + d00B_00F*D00F + d00B_00G*D00G + d00B_00H*D00H + d00B_00I*D00I + d00B_00J*D00J) + 2*s00B*E00B*D00B;
    
    ydot[34] =  -(s11B+(d11B_11A + d11B_00B + d11B_01B + d11B_11C + d11B_11D + d11B_11E + d11B_11F + d11B_11G + d11B_11H + d11B_11I + d11B_11J)+x11B) * D11B + (d11B_11A*D11A + d11B_00B*D00B + d11B_01B*D01B + d11B_11C*D11C + d11B_11D*D11D + d11B_11E*D11E + d11B_11F*D11F + d11B_11G*D11G + d11B_11H*D11H + d11B_11I*D11I + d11B_11J*D11J) + 2*s11B*E11B*D11B;
    
    ydot[35] =  -(s01B + s00B + s11B + (d01B_01A + d01B_00B + d01B_11B + d01B_01C + d01B_01D + d01B_01E + d01B_01F + d01B_01G + d01B_01H + d01B_01I + d01B_01J)) * D01B + (d01B_01A*D01A + d01B_00B*D00B + d01B_11B*D11B + d01B_01C*D01C + d01B_01D*D01D + d01B_01E*D01E + d01B_01F*D01F + d01B_01G*D01G + d01B_01H*D01H + d01B_01I*D01I + d01B_01J*D01J) + s00B*(E00B*D01B+E01B*D00B) + s11B*(E11B*D01B + E01B*D11B) + s01B*(E00B*D11B + E11B*D00B);
    
    ydot[36] =  -(s00C+(d00C_00A + d00C_00B + d00C_11C + d00C_01C + d00C_00D + d00C_00E + d00C_00F + d00C_00G + d00C_00H + d00C_00I + d00C_00J)+x00C) * D00C + (d00C_00A*D00A + d00C_00B*D00B + d00C_11C*D11C + d00C_01C*D01C + d00C_00D*D00D + d00C_00E*D00E + d00C_00F*D00F + d00C_00G*D00G + d00C_00H*D00H + d00C_00I*D00I + d00C_00J*D00J) + 2*s00C*E00C*D00C;
    
    ydot[37] =  -(s11C+(d11C_11A + d11C_11B + d11C_00C + d11C_01C + d11C_11D + d11C_11E + d11C_11F + d11C_11G + d11C_11H + d11C_11I + d11C_11J)+x11C) * D11C + (d11C_11A*D11A + d11C_11B*D11B + d11C_00C*D00C + d11C_01C*D01C + d11C_11D*D11D + d11C_11E*D11E + d11C_11F*D11F + d11C_11G*D11G + d11C_11H*D11H + d11C_11I*D11I + d11C_11J*D11J) + 2*s11C*E11C*D11C;
    
    ydot[38] =  -(s01C + s00C + s11C + (d01C_01A + d01C_01B + d01C_00C + d01C_11C + d01C_01D + d01C_01E + d01C_01F + d01C_01G + d01C_01H + d01C_01I + d01C_01J)) * D01C + (d01C_01A*D01A + d01C_01B*D01B + d01C_00C*D00C + d01C_11C*D11C + d01C_01D*D01D + d01C_01E*D01E + d01C_01F*D01F + d01C_01G*D01G + d01C_01H*D01H + d01C_01I*D01I + d01C_01J*D01J) + s00C*(E00C*D01C+E01C*D00C) + s11C*(E11C*D01C + E01C*D11C) + s01C*(E00C*D11C + E11C*D00C);
    
    ydot[39] =  -(s00D+(d00D_00A + d00D_00B + d00D_00C + d00D_11D + d00D_01D + d00D_00E + d00D_00F + d00D_00G + d00D_00H + d00D_00I + d00D_00J)+x00D) * D00D + (d00D_00A*D00A + d00D_00B*D00B + d00D_00C*D00C + d00D_11D*D11D + d00D_01D*D01D + d00D_00E*D00E + d00D_00F*D00F + d00D_00G*D00G + d00D_00H*D00H + d00D_00I*D00I + d00D_00J*D00J) + 2*s00D*E00D*D00D;
    
    ydot[40] =  -(s11D+(d11D_11A + d11D_11B + d11D_11C + d11D_00D + d11D_01D + d11D_11E + d11D_11F + d11D_11G + d11D_11H + d11D_11I + d11D_11J)+x11D) * D11D + (d11D_11A*D11A + d11D_11B*D11B + d11D_11C*D11C + d11D_00D*D00D + d11D_01D*D01D + d11D_11E*D11E + d11D_11F*D11F + d11D_11G*D11G + d11D_11H*D11H + d11D_11I*D11I + d11D_11J*D11J) + 2*s11D*E11D*D11D;
    
    ydot[41] =  -(s01D + s00D + s11D + (d01D_01A + d01D_01B + d01D_01C + d01D_00D + d01D_11D + d01D_01E + d01D_01F + d01D_01G + d01D_01H + d01D_01I + d01D_01J)) * D01D + (d01D_01A*D01A + d01D_01B*D01B + d01D_01C*D01C + d01D_00D*D00D + d01D_11D*D11D + d01D_01E*D01E + d01D_01F*D01F + d01D_01G*D01G + d01D_01H*D01H + d01D_01I*D01I + d01D_01J*D01J) + s00D*(E00D*D01D+E01D*D00D) + s11D*(E11D*D01D + E01D*D11D) + s01D*(E00D*D11D + E11D*D00D);
    
    ydot[42] =  -(s00E+(d00E_00A + d00E_00B + d00E_00C + d00E_00D + d00E_11E + d00E_01E + d00E_00F + d00E_00G + d00E_00H + d00E_00I + d00E_00J)+x00E) * D00E + (d00E_00A*D00A + d00E_00B*D00B + d00E_00C*D00C + d00E_00D*D00D + d00E_11E*D11E + d00E_01E*D01E + d00E_00F*D00F + d00E_00G*D00G + d00E_00H*D00H + d00E_00I*D00I + d00E_00J*D00J) + 2*s00E*E00E*D00E;
    
    ydot[43] =  -(s11E+(d11E_11A + d11E_11B + d11E_11C + d11E_11D + d11E_00E + d11E_01E + d11E_11F + d11E_11G + d11E_11H + d11E_11I + d11E_11J)+x11E) * D11E + (d11E_11A*D11A + d11E_11B*D11B + d11E_11C*D11C + d11E_11D*D11D + d11E_00E*D00E + d11E_01E*D01E + d11E_11F*D11F + d11E_11G*D11G + d11E_11H*D11H + d11E_11I*D11I + d11E_11J*D11J) + 2*s11E*E11E*D11E;
    
    ydot[44] =  -(s01E + s00E + s11E + (d01E_01A + d01E_01B + d01E_01C + d01E_01D + d01E_00E + d01E_11E + d01E_01F + d01E_01G + d01E_01H + d01E_01I + d01E_01J)) * D01E + (d01E_01A*D01A + d01E_01B*D01B + d01E_01C*D01C + d01E_01D*D01D + d01E_00E*D00E + d01E_11E*D11E + d01E_01F*D01F + d01E_01G*D01G + d01E_01H*D01H + d01E_01I*D01I + d01E_01J*D01J) + s00E*(E00E*D01E+E01E*D00E) + s11E*(E11E*D01E + E01E*D11E) + s01E*(E00E*D11E + E11E*D00E);
    
    ydot[45] =  -(s00F+(d00F_00A + d00F_00B + d00F_00C + d00F_00D + d00F_00E + d00F_11F + d00F_01F + d00F_00G + d00F_00H + d00F_00I + d00F_00J)+x00F) * D00F + (d00F_00A*D00A + d00F_00B*D00B + d00F_00C*D00C + d00F_00D*D00D + d00F_00E*D00E + d00F_11F*D11F + d00F_01F*D01F + d00F_00G*D00G + d00F_00H*D00H + d00F_00I*D00I + d00F_00J*D00J) + 2*s00F*E00F*D00F;
    
    ydot[46] =  -(s11F+(d11F_11A + d11F_11B + d11F_11C + d11F_11D + d11F_11E + d11F_00F + d11F_01F + d11F_11G + d11F_11H + d11F_11I + d11F_11J)+x11F) * D11F + (d11F_11A*D11A + d11F_11B*D11B + d11F_11C*D11C + d11F_11D*D11D + d11F_11E*D11E + d11F_00F*D00F + d11F_01F*D01F + d11F_11G*D11G + d11F_11H*D11H + d11F_11I*D11I + d11F_11J*D11J) + 2*s11F*E11F*D11F;
    
    ydot[47] =  -(s01F + s00F + s11F + (d01F_01A + d01F_01B + d01F_01C + d01F_01D + d01F_01E + d01F_00F + d01F_11F + d01F_01G + d01F_01H + d01F_01I + d01F_01J)) * D01F + (d01F_01A*D01A + d01F_01B*D01B + d01F_01C*D01C + d01F_01D*D01D + d01F_01E*D01E + d01F_00F*D00F + d01F_11F*D11F + d01F_01G*D01G + d01F_01H*D01H + d01F_01I*D01I + d01F_01J*D01J) + s00F*(E00F*D01F+E01F*D00F) + s11F*(E11F*D01F + E01F*D11F) + s01F*(E00F*D11F + E11F*D00F);
    
    ydot[48] =  -(s00G+(d00G_00A + d00G_00B + d00G_00C + d00G_00D + d00G_00E + d00G_00F + d00G_11G + d00G_01G + d00G_00H + d00G_00I + d00G_00J)+x00G) * D00G + (d00G_00A*D00A + d00G_00B*D00B + d00G_00C*D00C + d00G_00D*D00D + d00G_00E*D00E + d00G_00F*D00F + d00G_11G*D11G + d00G_01G*D01G + d00G_00H*D00H + d00G_00I*D00I + d00G_00J*D00J) + 2*s00G*E00G*D00G;
    
    ydot[49] =  -(s11G+(d11G_11A + d11G_11B + d11G_11C + d11G_11D + d11G_11E + d11G_11F + d11G_00G + d11G_01G + d11G_11H + d11G_11I + d11G_11J)+x11G) * D11G + (d11G_11A*D11A + d11G_11B*D11B + d11G_11C*D11C + d11G_11D*D11D + d11G_11E*D11E + d11G_11F*D11F + d11G_00G*D00G + d11G_01G*D01G + d11G_11H*D11H + d11G_11I*D11I + d11G_11J*D11J) + 2*s11G*E11G*D11G;
    
    ydot[50] =  -(s01G + s00G + s11G + (d01G_01A + d01G_01B + d01G_01C + d01G_01D + d01G_01E + d01G_01F + d01G_00G + d01G_11G + d01G_01H + d01G_01I + d01G_01J)) * D01G + (d01G_01A*D01A + d01G_01B*D01B + d01G_01C*D01C + d01G_01D*D01D + d01G_01E*D01E + d01G_01F*D01F + d01G_00G*D00G + d01G_11G*D11G + d01G_01H*D01H + d01G_01I*D01I + d01G_01J*D01J) + s00G*(E00G*D01G+E01G*D00G) + s11G*(E11G*D01G + E01G*D11G) + s01G*(E00G*D11G + E11G*D00G);
    
    ydot[51] =  -(s00H+(d00H_00A + d00H_00B + d00H_00C + d00H_00D + d00H_00E + d00H_00F + d00H_00G + d00H_11H + d00H_01H + d00H_00I + d00H_00J)+x00H) * D00H + (d00H_00A*D00A + d00H_00B*D00B + d00H_00C*D00C + d00H_00D*D00D + d00H_00E*D00E + d00H_00F*D00F + d00H_00G*D00G + d00H_11H*D11H + d00H_01H*D01H + d00H_00I*D00I + d00H_00J*D00J) + 2*s00H*E00H*D00H;
    
    ydot[52] =  -(s11H+(d11H_11A + d11H_11B + d11H_11C + d11H_11D + d11H_11E + d11H_11F + d11H_11G + d11H_00H + d11H_01H + d11H_11I + d11H_11J)+x11H) * D11H + (d11H_11A*D11A + d11H_11B*D11B + d11H_11C*D11C + d11H_11D*D11D + d11H_11E*D11E + d11H_11F*D11F + d11H_11G*D11G + d11H_00H*D00H + d11H_01H*D01H + d11H_11I*D11I + d11H_11J*D11J) + 2*s11H*E11H*D11H;
    
    ydot[53] =  -(s01H + s00H + s11H + (d01H_01A + d01H_01B + d01H_01C + d01H_01D + d01H_01E + d01H_01F + d01H_01G + d01H_00H + d01H_11H + d01H_01I + d01H_01J)) * D01H + (d01H_01A*D01A + d01H_01B*D01B + d01H_01C*D01C + d01H_01D*D01D + d01H_01E*D01E + d01H_01F*D01F + d01H_01G*D01G + d01H_00H*D00H + d01H_11H*D11H + d01H_01I*D01I + d01H_01J*D01J) + s00H*(E00H*D01H+E01H*D00H) + s11H*(E11H*D01H + E01H*D11H) + s01H*(E00H*D11H + E11H*D00H);
    
    ydot[54] =  -(s00I+(d00I_00A + d00I_00B + d00I_00C + d00I_00D + d00I_00E + d00I_00F + d00I_00G + d00I_00H + d00I_11I + d00I_01I + d00I_00J)+x00I) * D00I + (d00I_00A*D00A + d00I_00B*D00B + d00I_00C*D00C + d00I_00D*D00D + d00I_00E*D00E + d00I_00F*D00F + d00I_00G*D00G + d00I_00H*D00H + d00I_11I*D11I + d00I_01I*D01I + d00I_00J*D00J) + 2*s00I*E00I*D00I;
    
    ydot[55] =  -(s11I+(d11I_11A + d11I_11B + d11I_11C + d11I_11D + d11I_11E + d11I_11F + d11I_11G + d11I_11H + d11I_00I + d11I_01I + d11I_11J)+x11I) * D11I + (d11I_11A*D11A + d11I_11B*D11B + d11I_11C*D11C + d11I_11D*D11D + d11I_11E*D11E + d11I_11F*D11F + d11I_11G*D11G + d11I_11H*D11H + d11I_00I*D00I + d11I_01I*D01I + d11I_11J*D11J) + 2*s11I*E11I*D11I;
    
    ydot[56] =  -(s01I + s00I + s11I + (d01I_01A + d01I_01B + d01I_01C + d01I_01D + d01I_01E + d01I_01F + d01I_01G + d01I_01H + d01I_00I + d01I_11I + d01I_01J)) * D01I + (d01I_01A*D01A + d01I_01B*D01B + d01I_01C*D01C + d01I_01D*D01D + d01I_01E*D01E + d01I_01F*D01F + d01I_01G*D01G + d01I_01H*D01H + d01I_00I*D00I + d01I_11I*D11I + d01I_01J*D01J) + s00I*(E00I*D01I+E01I*D00I) + s11I*(E11I*D01I + E01I*D11I) + s01I*(E00I*D11I + E11I*D00I);
    
    ydot[57] =  -(s00J+(d00J_00A + d00J_00B + d00J_00C + d00J_00D + d00J_00E + d00J_00F + d00J_00G + d00J_00H + d00J_00I + d00J_11J + d00J_01J)+x00J) * D00J + (d00J_00A*D00A + d00J_00B*D00B + d00J_00C*D00C + d00J_00D*D00D + d00J_00E*D00E + d00J_00F*D00F + d00J_00G*D00G + d00J_00H*D00H + d00J_00I*D00I + d00J_11J*D11J + d00J_01J*D01J) + 2*s00J*E00J*D00J;
    
    ydot[58] =  -(s11J+(d11J_11A + d11J_11B + d11J_11C + d11J_11D + d11J_11E + d11J_11F + d11J_11G + d11J_11H + d11J_11I + d11J_00J + d11J_01J)+x11J) * D11J + (d11J_11A*D11A + d11J_11B*D11B + d11J_11C*D11C + d11J_11D*D11D + d11J_11E*D11E + d11J_11F*D11F + d11J_11G*D11G + d11J_11H*D11H + d11J_11I*D11I + d11J_00J*D00J + d11J_01J*D01J) + 2*s11J*E11J*D11J;
    
    ydot[59] =  -(s01J + s00J + s11J + (d01J_01A + d01J_01B + d01J_01C + d01J_01D + d01J_01E + d01J_01F + d01J_01G + d01J_01H + d01J_01I + d01J_00J + d01J_11J)) * D01J + (d01J_01A*D01A + d01J_01B*D01B + d01J_01C*D01C + d01J_01D*D01D + d01J_01E*D01E + d01J_01F*D01F + d01J_01G*D01G + d01J_01H*D01H + d01J_01I*D01I + d01J_00J*D00J + d01J_11J*D11J) + s00J*(E00J*D01J+E01J*D00J) + s11J*(E11J*D01J + E01J*D11J) + s01J*(E00J*D11J + E11J*D00J);
    
}
