/*
 *  muhisse-ext-derivs_c
 *
 *
 *  Created by Jeremy Beaulieu 4/25/2018
 *  Copyright 2018 Awesome Inc_ All rights reserved_
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 385

static double params_muhisse[NUMELEMENTS];


void initmod_muhisse(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_muhisse);
}




void muhisse_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E00A = y[0];
    double E01A = y[1];
    double E10A = y[2];
    double E11A = y[3];
    double E00B = y[4];
    double E01B = y[5];
    double E10B = y[6];
    double E11B = y[7];
    double E00C = y[8];
    double E01C = y[9];
    double E10C = y[10];
    double E11C = y[11];
    double E00D = y[12];
    double E01D = y[13];
    double E10D = y[14];
    double E11D = y[15];
    double E00E = y[16];
    double E01E = y[17];
    double E10E = y[18];
    double E11E = y[19];
    double E00F = y[20];
    double E01F = y[21];
    double E10F = y[22];
    double E11F = y[23];
    double E00G = y[24];
    double E01G = y[25];
    double E10G = y[26];
    double E11G = y[27];
    double E00H = y[28];
    double E01H = y[29];
    double E10H = y[30];
    double E11H = y[31];
    
    double D00A = y[32];
    double D01A = y[33];
    double D10A = y[34];
    double D11A = y[35];
    double D00B = y[36];
    double D01B = y[37];
    double D10B = y[38];
    double D11B = y[39];
    double D00C = y[40];
    double D01C = y[41];
    double D10C = y[42];
    double D11C = y[43];
    double D00D = y[44];
    double D01D = y[45];
    double D10D = y[46];
    double D11D = y[47];
    double D00E = y[48];
    double D01E = y[49];
    double D10E = y[50];
    double D11E = y[51];
    double D00F = y[52];
    double D01F = y[53];
    double D10F = y[54];
    double D11F = y[55];
    double D00G = y[56];
    double D01G = y[57];
    double D10G = y[58];
    double D11G = y[59];
    double D00H = y[60];
    double D01H = y[61];
    double D10H = y[62];
    double D11H = y[63];
    
    
    double
    lambda00A = params_muhisse[0],
    lambda01A = params_muhisse[1],
    lambda10A = params_muhisse[2],
    lambda11A = params_muhisse[3],
    mu00A = params_muhisse[4],
    mu01A = params_muhisse[5],
    mu10A = params_muhisse[6],
    mu11A = params_muhisse[7],
    q00A_01A = params_muhisse[8],
    q00A_10A = params_muhisse[9],
    q00A_11A = params_muhisse[10],
    q01A_00A = params_muhisse[11],
    q01A_10A = params_muhisse[12],
    q01A_11A = params_muhisse[13],
    q10A_00A = params_muhisse[14],
    q10A_01A = params_muhisse[15],
    q10A_11A = params_muhisse[16],
    q11A_00A = params_muhisse[17],
    q11A_01A = params_muhisse[18],
    q11A_10A = params_muhisse[19],
    q00A_00B = params_muhisse[20],
    q00A_00C = params_muhisse[21],
    q00A_00D = params_muhisse[22],
    q00A_00E = params_muhisse[23],
    q00A_00F = params_muhisse[24],
    q00A_00G = params_muhisse[25],
    q00A_00H = params_muhisse[26],
    q01A_01B = params_muhisse[27],
    q01A_01C = params_muhisse[28],
    q01A_01D = params_muhisse[29],
    q01A_01E = params_muhisse[30],
    q01A_01F = params_muhisse[31],
    q01A_01G = params_muhisse[32],
    q01A_01H = params_muhisse[33],
    q10A_10B = params_muhisse[34],
    q10A_10C = params_muhisse[35],
    q10A_10D = params_muhisse[36],
    q10A_10E = params_muhisse[37],
    q10A_10F = params_muhisse[38],
    q10A_10G = params_muhisse[39],
    q10A_10H = params_muhisse[40],
    q11A_11B = params_muhisse[41],
    q11A_11C = params_muhisse[42],
    q11A_11D = params_muhisse[43],
    q11A_11E = params_muhisse[44],
    q11A_11F = params_muhisse[45],
    q11A_11G = params_muhisse[46],
    q11A_11H = params_muhisse[47],
    lambda00B = params_muhisse[48],
    lambda01B = params_muhisse[49],
    lambda10B = params_muhisse[50],
    lambda11B = params_muhisse[51],
    mu00B = params_muhisse[52],
    mu01B = params_muhisse[53],
    mu10B = params_muhisse[54],
    mu11B = params_muhisse[55],
    q00B_01B = params_muhisse[56],
    q00B_10B = params_muhisse[57],
    q00B_11B = params_muhisse[58],
    q01B_00B = params_muhisse[59],
    q01B_10B = params_muhisse[60],
    q01B_11B = params_muhisse[61],
    q10B_00B = params_muhisse[62],
    q10B_01B = params_muhisse[63],
    q10B_11B = params_muhisse[64],
    q11B_00B = params_muhisse[65],
    q11B_01B = params_muhisse[66],
    q11B_10B = params_muhisse[67],
    q00B_00A = params_muhisse[68],
    q00B_00C = params_muhisse[69],
    q00B_00D = params_muhisse[70],
    q00B_00E = params_muhisse[71],
    q00B_00F = params_muhisse[72],
    q00B_00G = params_muhisse[73],
    q00B_00H = params_muhisse[74],
    q01B_01A = params_muhisse[75],
    q01B_01C = params_muhisse[76],
    q01B_01D = params_muhisse[77],
    q01B_01E = params_muhisse[78],
    q01B_01F = params_muhisse[79],
    q01B_01G = params_muhisse[80],
    q01B_01H = params_muhisse[81],
    q10B_10A = params_muhisse[82],
    q10B_10C = params_muhisse[83],
    q10B_10D = params_muhisse[84],
    q10B_10E = params_muhisse[85],
    q10B_10F = params_muhisse[86],
    q10B_10G = params_muhisse[87],
    q10B_10H = params_muhisse[88],
    q11B_11A = params_muhisse[89],
    q11B_11C = params_muhisse[90],
    q11B_11D = params_muhisse[91],
    q11B_11E = params_muhisse[92],
    q11B_11F = params_muhisse[93],
    q11B_11G = params_muhisse[94],
    q11B_11H = params_muhisse[95],
    lambda00C = params_muhisse[96],
    lambda01C = params_muhisse[97],
    lambda10C = params_muhisse[98],
    lambda11C = params_muhisse[99],
    mu00C = params_muhisse[100],
    mu01C = params_muhisse[101],
    mu10C = params_muhisse[102],
    mu11C = params_muhisse[103],
    q00C_01C = params_muhisse[104],
    q00C_10C = params_muhisse[105],
    q00C_11C = params_muhisse[106],
    q01C_00C = params_muhisse[107],
    q01C_10C = params_muhisse[108],
    q01C_11C = params_muhisse[109],
    q10C_00C = params_muhisse[110],
    q10C_01C = params_muhisse[111],
    q10C_11C = params_muhisse[112],
    q11C_00C = params_muhisse[113],
    q11C_01C = params_muhisse[114],
    q11C_10C = params_muhisse[115],
    q00C_00A = params_muhisse[116],
    q00C_00B = params_muhisse[117],
    q00C_00D = params_muhisse[118],
    q00C_00E = params_muhisse[119],
    q00C_00F = params_muhisse[120],
    q00C_00G = params_muhisse[121],
    q00C_00H = params_muhisse[122],
    q01C_01A = params_muhisse[123],
    q01C_01B = params_muhisse[124],
    q01C_01D = params_muhisse[125],
    q01C_01E = params_muhisse[126],
    q01C_01F = params_muhisse[127],
    q01C_01G = params_muhisse[128],
    q01C_01H = params_muhisse[129],
    q10C_10A = params_muhisse[130],
    q10C_10B = params_muhisse[131],
    q10C_10D = params_muhisse[132],
    q10C_10E = params_muhisse[133],
    q10C_10F = params_muhisse[134],
    q10C_10G = params_muhisse[135],
    q10C_10H = params_muhisse[136],
    q11C_11A = params_muhisse[137],
    q11C_11B = params_muhisse[138],
    q11C_11D = params_muhisse[139],
    q11C_11E = params_muhisse[140],
    q11C_11F = params_muhisse[141],
    q11C_11G = params_muhisse[142],
    q11C_11H = params_muhisse[143],
    lambda00D = params_muhisse[144],
    lambda01D = params_muhisse[145],
    lambda10D = params_muhisse[146],
    lambda11D = params_muhisse[147],
    mu00D = params_muhisse[148],
    mu01D = params_muhisse[149],
    mu10D = params_muhisse[150],
    mu11D = params_muhisse[151],
    q00D_01D = params_muhisse[152],
    q00D_10D = params_muhisse[153],
    q00D_11D = params_muhisse[154],
    q01D_00D = params_muhisse[155],
    q01D_10D = params_muhisse[156],
    q01D_11D = params_muhisse[157],
    q10D_00D = params_muhisse[158],
    q10D_01D = params_muhisse[159],
    q10D_11D = params_muhisse[160],
    q11D_00D = params_muhisse[161],
    q11D_01D = params_muhisse[162],
    q11D_10D = params_muhisse[163],
    q00D_00A = params_muhisse[164],
    q00D_00B = params_muhisse[165],
    q00D_00C = params_muhisse[166],
    q00D_00E = params_muhisse[167],
    q00D_00F = params_muhisse[168],
    q00D_00G = params_muhisse[169],
    q00D_00H = params_muhisse[170],
    q01D_01A = params_muhisse[171],
    q01D_01B = params_muhisse[172],
    q01D_01C = params_muhisse[173],
    q01D_01E = params_muhisse[174],
    q01D_01F = params_muhisse[175],
    q01D_01G = params_muhisse[176],
    q01D_01H = params_muhisse[177],
    q10D_10A = params_muhisse[178],
    q10D_10B = params_muhisse[179],
    q10D_10C = params_muhisse[180],
    q10D_10E = params_muhisse[181],
    q10D_10F = params_muhisse[182],
    q10D_10G = params_muhisse[183],
    q10D_10H = params_muhisse[184],
    q11D_11A = params_muhisse[185],
    q11D_11B = params_muhisse[186],
    q11D_11C = params_muhisse[187],
    q11D_11E = params_muhisse[188],
    q11D_11F = params_muhisse[189],
    q11D_11G = params_muhisse[190],
    q11D_11H = params_muhisse[191],
    lambda00E = params_muhisse[192],
    lambda01E = params_muhisse[193],
    lambda10E = params_muhisse[194],
    lambda11E = params_muhisse[195],
    mu00E = params_muhisse[196],
    mu01E = params_muhisse[197],
    mu10E = params_muhisse[198],
    mu11E = params_muhisse[199],
    q00E_01E = params_muhisse[200],
    q00E_10E = params_muhisse[201],
    q00E_11E = params_muhisse[202],
    q01E_00E = params_muhisse[203],
    q01E_10E = params_muhisse[204],
    q01E_11E = params_muhisse[205],
    q10E_00E = params_muhisse[206],
    q10E_01E = params_muhisse[207],
    q10E_11E = params_muhisse[208],
    q11E_00E = params_muhisse[209],
    q11E_01E = params_muhisse[210],
    q11E_10E = params_muhisse[211],
    q00E_00A = params_muhisse[212],
    q00E_00B = params_muhisse[213],
    q00E_00C = params_muhisse[214],
    q00E_00D = params_muhisse[215],
    q00E_00F = params_muhisse[216],
    q00E_00G = params_muhisse[217],
    q00E_00H = params_muhisse[218],
    q01E_01A = params_muhisse[219],
    q01E_01B = params_muhisse[220],
    q01E_01C = params_muhisse[221],
    q01E_01D = params_muhisse[222],
    q01E_01F = params_muhisse[223],
    q01E_01G = params_muhisse[224],
    q01E_01H = params_muhisse[225],
    q10E_10A = params_muhisse[226],
    q10E_10B = params_muhisse[227],
    q10E_10C = params_muhisse[228],
    q10E_10D = params_muhisse[229],
    q10E_10F = params_muhisse[230],
    q10E_10G = params_muhisse[231],
    q10E_10H = params_muhisse[232],
    q11E_11A = params_muhisse[233],
    q11E_11B = params_muhisse[234],
    q11E_11C = params_muhisse[235],
    q11E_11D = params_muhisse[236],
    q11E_11F = params_muhisse[237],
    q11E_11G = params_muhisse[238],
    q11E_11H = params_muhisse[239],
    lambda00F = params_muhisse[240],
    lambda01F = params_muhisse[241],
    lambda10F = params_muhisse[242],
    lambda11F = params_muhisse[243],
    mu00F = params_muhisse[244],
    mu01F = params_muhisse[245],
    mu10F = params_muhisse[246],
    mu11F = params_muhisse[247],
    q00F_01F = params_muhisse[248],
    q00F_10F = params_muhisse[249],
    q00F_11F = params_muhisse[250],
    q01F_00F = params_muhisse[251],
    q01F_10F = params_muhisse[252],
    q01F_11F = params_muhisse[253],
    q10F_00F = params_muhisse[254],
    q10F_01F = params_muhisse[255],
    q10F_11F = params_muhisse[256],
    q11F_00F = params_muhisse[257],
    q11F_01F = params_muhisse[258],
    q11F_10F = params_muhisse[259],
    q00F_00A = params_muhisse[260],
    q00F_00B = params_muhisse[261],
    q00F_00C = params_muhisse[262],
    q00F_00D = params_muhisse[263],
    q00F_00E = params_muhisse[264],
    q00F_00G = params_muhisse[265],
    q00F_00H = params_muhisse[266],
    q01F_01A = params_muhisse[267],
    q01F_01B = params_muhisse[268],
    q01F_01C = params_muhisse[269],
    q01F_01D = params_muhisse[270],
    q01F_01E = params_muhisse[271],
    q01F_01G = params_muhisse[272],
    q01F_01H = params_muhisse[273],
    q10F_10A = params_muhisse[274],
    q10F_10B = params_muhisse[275],
    q10F_10C = params_muhisse[276],
    q10F_10D = params_muhisse[277],
    q10F_10E = params_muhisse[278],
    q10F_10G = params_muhisse[279],
    q10F_10H = params_muhisse[280],
    q11F_11A = params_muhisse[281],
    q11F_11B = params_muhisse[282],
    q11F_11C = params_muhisse[283],
    q11F_11D = params_muhisse[284],
    q11F_11E = params_muhisse[285],
    q11F_11G = params_muhisse[286],
    q11F_11H = params_muhisse[287],
    lambda00G = params_muhisse[288],
    lambda01G = params_muhisse[289],
    lambda10G = params_muhisse[290],
    lambda11G = params_muhisse[291],
    mu00G = params_muhisse[292],
    mu01G = params_muhisse[293],
    mu10G = params_muhisse[294],
    mu11G = params_muhisse[295],
    q00G_01G = params_muhisse[296],
    q00G_10G = params_muhisse[297],
    q00G_11G = params_muhisse[298],
    q01G_00G = params_muhisse[299],
    q01G_10G = params_muhisse[300],
    q01G_11G = params_muhisse[301],
    q10G_00G = params_muhisse[302],
    q10G_01G = params_muhisse[303],
    q10G_11G = params_muhisse[304],
    q11G_00G = params_muhisse[305],
    q11G_01G = params_muhisse[306],
    q11G_10G = params_muhisse[307],
    q00G_00A = params_muhisse[308],
    q00G_00B = params_muhisse[309],
    q00G_00C = params_muhisse[310],
    q00G_00D = params_muhisse[311],
    q00G_00E = params_muhisse[312],
    q00G_00F = params_muhisse[313],
    q00G_00H = params_muhisse[314],
    q01G_01A = params_muhisse[315],
    q01G_01B = params_muhisse[316],
    q01G_01C = params_muhisse[317],
    q01G_01D = params_muhisse[318],
    q01G_01E = params_muhisse[319],
    q01G_01F = params_muhisse[320],
    q01G_01H = params_muhisse[321],
    q10G_10A = params_muhisse[322],
    q10G_10B = params_muhisse[323],
    q10G_10C = params_muhisse[324],
    q10G_10D = params_muhisse[325],
    q10G_10E = params_muhisse[326],
    q10G_10F = params_muhisse[327],
    q10G_10H = params_muhisse[328],
    q11G_11A = params_muhisse[329],
    q11G_11B = params_muhisse[330],
    q11G_11C = params_muhisse[331],
    q11G_11D = params_muhisse[332],
    q11G_11E = params_muhisse[333],
    q11G_11F = params_muhisse[334],
    q11G_11H = params_muhisse[335],
    lambda00H = params_muhisse[336],
    lambda01H = params_muhisse[337],
    lambda10H = params_muhisse[338],
    lambda11H = params_muhisse[339],
    mu00H = params_muhisse[340],
    mu01H = params_muhisse[341],
    mu10H = params_muhisse[342],
    mu11H = params_muhisse[343],
    q00H_01H = params_muhisse[344],
    q00H_10H = params_muhisse[345],
    q00H_11H = params_muhisse[346],
    q01H_00H = params_muhisse[347],
    q01H_10H = params_muhisse[348],
    q01H_11H = params_muhisse[349],
    q10H_00H = params_muhisse[350],
    q10H_01H = params_muhisse[351],
    q10H_11H = params_muhisse[352],
    q11H_00H = params_muhisse[353],
    q11H_01H = params_muhisse[354],
    q11H_10H = params_muhisse[355],
    q00H_00A = params_muhisse[356],
    q00H_00B = params_muhisse[357],
    q00H_00C = params_muhisse[358],
    q00H_00D = params_muhisse[359],
    q00H_00E = params_muhisse[360],
    q00H_00F = params_muhisse[361],
    q00H_00G = params_muhisse[362],
    q01H_01A = params_muhisse[363],
    q01H_01B = params_muhisse[364],
    q01H_01C = params_muhisse[365],
    q01H_01D = params_muhisse[366],
    q01H_01E = params_muhisse[367],
    q01H_01F = params_muhisse[368],
    q01H_01G = params_muhisse[369],
    q10H_10A = params_muhisse[370],
    q10H_10B = params_muhisse[371],
    q10H_10C = params_muhisse[372],
    q10H_10D = params_muhisse[373],
    q10H_10E = params_muhisse[374],
    q10H_10F = params_muhisse[375],
    q10H_10G = params_muhisse[376],
    q11H_11A = params_muhisse[377],
    q11H_11B = params_muhisse[378],
    q11H_11C = params_muhisse[379],
    q11H_11D = params_muhisse[380],
    q11H_11E = params_muhisse[381],
    q11H_11F = params_muhisse[382],
    q11H_11G = params_muhisse[383],
    psi = params_muhisse[384];
    
    /* The E's */
    ydot[0] = mu00A - (lambda00A + psi + (q00A_01A + q00A_10A + q00A_11A + q00A_00B + q00A_00C + q00A_00D + q00A_00E + q00A_00F + q00A_00G + q00A_00H) + mu00A) * E00A + lambda00A*E00A*E00A + (q00A_01A*E01A + q00A_10A*E10A + q00A_11A*E11A + q00A_00B*E00B + q00A_00C*E00C + q00A_00D*E00D + q00A_00E*E00E + q00A_00F*E00F + q00A_00G*E00G + q00A_00H*E00H);
    
    ydot[1] = mu01A - (lambda01A + psi + (q01A_00A + q01A_10A + q01A_11A + q01A_01B + q01A_01C + q01A_01D + q01A_01E + q01A_01F + q01A_01G + q01A_01H) + mu01A) * E01A + lambda01A*E01A*E01A + (q01A_00A*E00A + q01A_10A*E10A + q01A_11A*E11A + q01A_01B*E01B + q01A_01C*E01C + q01A_01D*E01D + q01A_01E*E01E + q01A_01F*E01F + q01A_01G*E01G + q01A_01H*E01H);
    
    ydot[2] = mu10A - (lambda10A + psi + (q10A_00A + q10A_01A + q10A_11A + q10A_10B + q10A_10C + q10A_10D + q10A_10E + q10A_10F + q10A_10G + q10A_10H) + mu10A) * E10A + lambda10A*E10A*E10A + (q10A_00A*E00A + q10A_01A*E01A + q10A_11A*E11A + q10A_10B*E10B + q10A_10C*E10C + q10A_10D*E10D + q10A_10E*E10E + q10A_10F*E10F + q10A_10G*E10G + q10A_10H*E10H);
    
    ydot[3] = mu11A - (lambda11A + psi + (q11A_00A + q11A_01A + q11A_10A + q11A_11B + q11A_11C + q11A_11D + q11A_11E + q11A_11F + q11A_11G + q11A_11H) + mu11A) * E11A + lambda11A*E11A*E11A + (q11A_00A*E00A + q11A_01A*E01A + q11A_10A*E10A + q11A_11B*E11B + q11A_11C*E11C + q11A_11D*E11D + q11A_11E*E11E + q11A_11F*E11F + q11A_11G*E11G + q11A_11H*E11H);
    
    ydot[4] = mu00B - (lambda00B + psi + (q00B_00A + q00B_01B + q00B_10B + q00B_11B + q00B_00C + q00B_00D + q00B_00E + q00B_00F + q00B_00G + q00B_00H) + mu00B) * E00B + lambda00B*E00B*E00B + (q00B_00A*E00A + q00B_01B*E01B + q00B_10B*E10B + q00B_11B*E11B + q00B_00C*E00C + q00B_00D*E00D + q00B_00E*E00E + q00B_00F*E00F + q00B_00G*E00G + q00B_00H*E00H);
    
    ydot[5] = mu01B - (lambda01B + psi + (q01B_01A + q01B_00B + q01B_10B + q01B_11B + q01B_01C + q01B_01D + q01B_01E + q01B_01F + q01B_01G + q01B_01H) + mu01B) * E01B + lambda01B*E01B*E01B + (q01B_01A*E01A + q01B_00B*E00B + q01B_10B*E10B + q01B_11B*E11B + q01B_01C*E01C + q01B_01D*E01D + q01B_01E*E01E + q01B_01F*E01F + q01B_01G*E01G + q01B_01H*E01H);
    
    ydot[6] = mu10B - (lambda10B + psi + (q10B_10A + q10B_00B + q10B_01B + q10B_11B + q10B_10C + q10B_10D + q10B_10E + q10B_10F + q10B_10G + q10B_10H) + mu10B) * E10B + lambda10B*E10B*E10B + (q10B_10A*E10A + q10B_00B*E00B + q10B_01B*E01B + q10B_11B*E11B + q10B_10C*E10C + q10B_10D*E10D + q10B_10E*E10E + q10B_10F*E10F + q10B_10G*E10G + q10B_10H*E10H);
    
    ydot[7] = mu11B - (lambda11B + psi + (q11B_11A + q11B_00B + q11B_01B + q11B_10B + q11B_11C + q11B_11D + q11B_11E + q11B_11F + q11B_11G + q11B_11H) + mu11B) * E11B + lambda11B*E11B*E11B + (q11B_11A*E11A + q11B_00B*E00B + q11B_01B*E01B + q11B_10B*E10B + q11B_11C*E11C + q11B_11D*E11D + q11B_11E*E11E + q11B_11F*E11F + q11B_11G*E11G + q11B_11H*E11H);
    
    ydot[8] = mu00C - (lambda00C + psi + (q00C_00A + q00C_00B + q00C_01C + q00C_10C + q00C_11C + q00C_00D + q00C_00E + q00C_00F + q00C_00G + q00C_00H) + mu00C) * E00C + lambda00C*E00C*E00C + (q00C_00A*E00A + q00C_00B*E00B + q00C_01C*E01C + q00C_10C*E10C + q00C_11C*E11C + q00C_00D*E00D + q00C_00E*E00E + q00C_00F*E00F + q00C_00G*E00G + q00C_00H*E00H);
    
    ydot[9] = mu01C - (lambda01C + psi + (q01C_01A + q01C_01B + q01C_00C + q01C_10C + q01C_11C + q01C_01D + q01C_01E + q01C_01F + q01C_01G + q01C_01H) + mu01C) * E01C + lambda01C*E01C*E01C + (q01C_01A*E01A + q01C_01B*E01B + q01C_00C*E00C + q01C_10C*E10C + q01C_11C*E11C + q01C_01D*E01D + q01C_01E*E01E + q01C_01F*E01F + q01C_01G*E01G + q01C_01H*E01H);
    
    ydot[10] = mu10C - (lambda10C + psi + (q10C_10A + q10C_10B + q10C_00C + q10C_01C + q10C_11C + q10C_10D + q10C_10E + q10C_10F + q10C_10G + q10C_10H) + mu10C) * E10C + lambda10C*E10C*E10C + (q10C_10A*E10A + q10C_10B*E10B + q10C_00C*E00C + q10C_01C*E01C + q10C_11C*E11C + q10C_10D*E10D + q10C_10E*E10E + q10C_10F*E10F + q10C_10G*E10G + q10C_10H*E10H);
    
    ydot[11] = mu11C - (lambda11C + psi + (q11C_11A + q11C_11B + q11C_00C + q11C_01C + q11C_10C + q11C_11D + q11C_11E + q11C_11F + q11C_11G + q11C_11H) + mu11C) * E11C + lambda11C*E11C*E11C + (q11C_11A*E11A + q11C_11B*E11B + q11C_00C*E00C + q11C_01C*E01C + q11C_10C*E10C + q11C_11D*E11D + q11C_11E*E11E + q11C_11F*E11F + q11C_11G*E11G + q11C_11H*E11H);
    
    ydot[12] = mu00D - (lambda00D + psi + (q00D_00A + q00D_00B + q00D_00C + q00D_01D + q00D_10D + q00D_11D + q00D_00E + q00D_00F + q00D_00G + q00D_00H) + mu00D) * E00D + lambda00D*E00D*E00D + (q00D_00A*E00A + q00D_00B*E00B + q00D_00C*E00C + q00D_01D*E01D + q00D_10D*E10D + q00D_11D*E11D + q00D_00E*E00E + q00D_00F*E00F + q00D_00G*E00G + q00D_00H*E00H);
    
    ydot[13] = mu01D - (lambda01D + psi + (q01D_01A + q01D_01B + q01D_01C + q01D_00D + q01D_10D + q01D_11D + q01D_01E + q01D_01F + q01D_01G + q01D_01H) + mu01D) * E01D + lambda01D*E01D*E01D + (q01D_01A*E01A + q01D_01B*E01B + q01D_01C*E01C + q01D_00D*E00D + q01D_10D*E10D + q01D_11D*E11D + q01D_01E*E01E + q01D_01F*E01F + q01D_01G*E01G + q01D_01H*E01H);
    
    ydot[14] = mu10D - (lambda10D + psi + (q10D_10A + q10D_10B + q10D_10C + q10D_00D + q10D_01D + q10D_11D + q10D_10E + q10D_10F + q10D_10G + q10D_10H) + mu10D) * E10D + lambda10D*E10D*E10D + (q10D_10A*E10A + q10D_10B*E10B + q10D_10C*E10C + q10D_00D*E00D + q10D_01D*E01D + q10D_11D*E11D + q10D_10E*E10E + q10D_10F*E10F + q10D_10G*E10G + q10D_10H*E10H);
    
    ydot[15] = mu11D - (lambda11D + psi + (q11D_11A + q11D_11B + q11D_11C + q11D_00D + q11D_01D + q11D_10D + q11D_11E + q11D_11F + q11D_11G + q11D_11H) + mu11D) * E11D + lambda11D*E11D*E11D + (q11D_11A*E11A + q11D_11B*E11B + q11D_11C*E11C + q11D_00D*E00D + q11D_01D*E01D + q11D_10D*E10D + q11D_11E*E11E + q11D_11F*E11F + q11D_11G*E11G + q11D_11H*E11H);
    
    ydot[16] = mu00E - (lambda00E + psi + (q00E_00A + q00E_00B + q00E_00C + q00E_00D + q00E_01E + q00E_10E + q00E_11E + q00E_00F + q00E_00G + q00E_00H) + mu00E) * E00E + lambda00E*E00E*E00E + (q00E_00A*E00A + q00E_00B*E00B + q00E_00C*E00C + q00E_00D*E00D + q00E_01E*E01E + q00E_10E*E10E + q00E_11E*E11E + q00E_00F*E00F + q00E_00G*E00G + q00E_00H*E00H);
    
    ydot[17] = mu01E - (lambda01E + psi + (q01E_01A + q01E_01B + q01E_01C + q01E_01D + q01E_00E + q01E_10E + q01E_11E + q01E_01F + q01E_01G + q01E_01H) + mu01E) * E01E + lambda01E*E01E*E01E + (q01E_01A*E01A + q01E_01B*E01B + q01E_01C*E01C + q01E_01D*E01D + q01E_00E*E00E + q01E_10E*E10E + q01E_11E*E11E + q01E_01F*E01F + q01E_01G*E01G + q01E_01H*E01H);
    
    ydot[18] = mu10E - (lambda10E + psi + (q10E_10A + q10E_10B + q10E_10C + q10E_10D + q10E_00E + q10E_01E + q10E_11E + q10E_10F + q10E_10G + q10E_10H) + mu10E) * E10E + lambda10E*E10E*E10E + (q10E_10A*E10A + q10E_10B*E10B + q10E_10C*E10C + q10E_10D*E10D + q10E_00E*E00E + q10E_01E*E01E + q10E_11E*E11E + q10E_10F*E10F + q10E_10G*E10G + q10E_10H*E10H);
    
    ydot[19] = mu11E - (lambda11E + psi + (q11E_11A + q11E_11B + q11E_11C + q11E_11D + q11E_00E + q11E_01E + q11E_10E + q11E_11F + q11E_11G + q11E_11H) + mu11E) * E11E + lambda11E*E11E*E11E + (q11E_11A*E11A + q11E_11B*E11B + q11E_11C*E11C + q11E_11D*E11D + q11E_00E*E00E + q11E_01E*E01E + q11E_10E*E10E + q11E_11F*E11F + q11E_11G*E11G + q11E_11H*E11H);
    
    ydot[20] = mu00F - (lambda00F + psi + (q00F_00A + q00F_00B + q00F_00C + q00F_00D + q00F_00E + q00F_01F + q00F_10F + q00F_11F + q00F_00G + q00F_00H) + mu00F) * E00F + lambda00F*E00F*E00F + (q00F_00A*E00A + q00F_00B*E00B + q00F_00C*E00C + q00F_00D*E00D + q00F_00E*E00E + q00F_01F*E01F + q00F_10F*E10F + q00F_11F*E11F + q00F_00G*E00G + q00F_00H*E00H);
    
    ydot[21] = mu01F - (lambda01F + psi + (q01F_01A + q01F_01B + q01F_01C + q01F_01D + q01F_01E + q01F_00F + q01F_10F + q01F_11F + q01F_01G + q01F_01H) + mu01F) * E01F + lambda01F*E01F*E01F + (q01F_01A*E01A + q01F_01B*E01B + q01F_01C*E01C + q01F_01D*E01D + q01F_01E*E01E + q01F_00F*E00F + q01F_10F*E10F + q01F_11F*E11F + q01F_01G*E01G + q01F_01H*E01H);
    
    ydot[22] = mu10F - (lambda10F + psi + (q10F_10A + q10F_10B + q10F_10C + q10F_10D + q10F_10E + q10F_00F + q10F_01F + q10F_11F + q10F_10G + q10F_10H) + mu10F) * E10F + lambda10F*E10F*E10F + (q10F_10A*E10A + q10F_10B*E10B + q10F_10C*E10C + q10F_10D*E10D + q10F_10E*E10E + q10F_00F*E00F + q10F_01F*E01F + q10F_11F*E11F + q10F_10G*E10G + q10F_10H*E10H);
    
    ydot[23] = mu11F - (lambda11F + psi + (q11F_11A + q11F_11B + q11F_11C + q11F_11D + q11F_11E + q11F_00F + q11F_01F + q11F_10F + q11F_11G + q11F_11H) + mu11F) * E11F + lambda11F*E11F*E11F + (q11F_11A*E11A + q11F_11B*E11B + q11F_11C*E11C + q11F_11D*E11D + q11F_11E*E11E + q11F_00F*E00F + q11F_01F*E01F + q11F_10F*E10F + q11F_11G*E11G + q11F_11H*E11H);
    
    ydot[24] = mu00G - (lambda00G + psi + (q00G_00A + q00G_00B + q00G_00C + q00G_00D + q00G_00E + q00G_00F + q00G_01G + q00G_10G + q00G_11G + q00G_00H) + mu00G) * E00G + lambda00G*E00G*E00G + (q00G_00A*E00A + q00G_00B*E00B + q00G_00C*E00C + q00G_00D*E00D + q00G_00E*E00E + q00G_00F*E00F + q00G_01G*E01G + q00G_10G*E10G + q00G_11G*E11G + q00G_00H*E00H);
    
    ydot[25] = mu01G - (lambda01G + psi + (q01G_01A + q01G_01B + q01G_01C + q01G_01D + q01G_01E + q01G_01F + q01G_00G + q01G_10G + q01G_11G + q01G_01H) + mu01G) * E01G + lambda01G*E01G*E01G + (q01G_01A*E01A + q01G_01B*E01B + q01G_01C*E01C + q01G_01D*E01D + q01G_01E*E01E + q01G_01F*E01F + q01G_00G*E00G + q01G_10G*E10G + q01G_11G*E11G + q01G_01H*E01H);
    
    ydot[26] = mu10G - (lambda10G + psi + (q10G_10A + q10G_10B + q10G_10C + q10G_10D + q10G_10E + q10G_10F + q10G_00G + q10G_01G + q10G_11G + q10G_10H) + mu10G) * E10G + lambda10G*E10G*E10G + (q10G_10A*E10A + q10G_10B*E10B + q10G_10C*E10C + q10G_10D*E10D + q10G_10E*E10E + q10G_10F*E10F + q10G_00G*E00G + q10G_01G*E01G + q10G_11G*E11G + q10G_10H*E10H);
    
    ydot[27] = mu11G - (lambda11G + psi + (q11G_11A + q11G_11B + q11G_11C + q11G_11D + q11G_11E + q11G_11F + q11G_00G + q11G_01G + q11G_10G + q11G_11H) + mu11G) * E11G + lambda11G*E11G*E11G + (q11G_11A*E11A + q11G_11B*E11B + q11G_11C*E11C + q11G_11D*E11D + q11G_11E*E11E + q11G_11F*E11F + q11G_00G*E00G + q11G_01G*E01G + q11G_10G*E10G + q11G_11H*E11H);
    
    ydot[28] = mu00H - (lambda00H + psi + (q00H_00A + q00H_00B + q00H_00C + q00H_00D + q00H_00E + q00H_00F + q00H_00G + q00H_01H + q00H_10H + q00H_11H) + mu00H) * E00H + lambda00H*E00H*E00H + (q00H_00A*E00A + q00H_00B*E00B + q00H_00C*E00C + q00H_00D*E00D + q00H_00E*E00E + q00H_00F*E00F + q00H_00G*E00G + q00H_01H*E01H + q00H_10H*E10H + q00H_11H*E11H);
    
    ydot[29] = mu01H - (lambda01H + psi + (q01H_01A + q01H_01B + q01H_01C + q01H_01D + q01H_01E + q01H_01F + q01H_01G + q01H_00H + q01H_10H + q01H_11H) + mu01H) * E01H + lambda01H*E01H*E01H + (q01H_01A*E01A + q01H_01B*E01B + q01H_01C*E01C + q01H_01D*E01D + q01H_01E*E01E + q01H_01F*E01F + q01H_01G*E01G + q01H_00H*E00H + q01H_10H*E10H + q01H_11H*E11H);
    
    ydot[30] = mu10H - (lambda10H + psi + (q10H_10A + q10H_10B + q10H_10C + q10H_10D + q10H_10E + q10H_10F + q10H_10G + q10H_00H + q10H_01H + q10H_11H) + mu10H) * E10H + lambda10H*E10H*E10H + (q10H_10A*E10A + q10H_10B*E10B + q10H_10C*E10C + q10H_10D*E10D + q10H_10E*E10E + q10H_10F*E10F + q10H_10G*E10G + q10H_00H*E00H + q10H_01H*E01H + q10H_11H*E11H);
    
    ydot[31] = mu11H - (lambda11H + psi + (q11H_11A + q11H_11B + q11H_11C + q11H_11D + q11H_11E + q11H_11F + q11H_11G + q11H_00H + q11H_01H + q11H_10H) + mu11H) * E11H + lambda11H*E11H*E11H + (q11H_11A*E11A + q11H_11B*E11B + q11H_11C*E11C + q11H_11D*E11D + q11H_11E*E11E + q11H_11F*E11F + q11H_11G*E11G + q11H_00H*E00H + q11H_01H*E01H + q11H_10H*E10H);
    
    
    /* The D's */
    ydot[32] =  - (lambda00A + psi + (q00A_01A + q00A_10A + q00A_11A + q00A_00B + q00A_00C + q00A_00D + q00A_00E + q00A_00F + q00A_00G + q00A_00H) + mu00A) * D00A + 2*lambda00A*E00A*D00A + (q00A_01A*D01A + q00A_10A*D10A + q00A_11A*D11A + q00A_00B*D00B + q00A_00C*D00C + q00A_00D*D00D + q00A_00E*D00E + q00A_00F*D00F + q00A_00G*D00G + q00A_00H*D00H);
    
    ydot[33] =  - (lambda01A + psi + (q01A_00A + q01A_10A + q01A_11A + q01A_01B + q01A_01C + q01A_01D + q01A_01E + q01A_01F + q01A_01G + q01A_01H) + mu01A) * D01A + 2*lambda01A*E01A*D01A + (q01A_00A*D00A + q01A_10A*D10A + q01A_11A*D11A + q01A_01B*D01B + q01A_01C*D01C + q01A_01D*D01D + q01A_01E*D01E + q01A_01F*D01F + q01A_01G*D01G + q01A_01H*D01H);
    
    ydot[34] =  - (lambda10A + psi + (q10A_00A + q10A_01A + q10A_11A + q10A_10B + q10A_10C + q10A_10D + q10A_10E + q10A_10F + q10A_10G + q10A_10H) + mu10A) * D10A + 2*lambda10A*E10A*D10A + (q10A_00A*D00A + q10A_01A*D01A + q10A_11A*D11A + q10A_10B*D10B + q10A_10C*D10C + q10A_10D*D10D + q10A_10E*D10E + q10A_10F*D10F + q10A_10G*D10G + q10A_10H*D10H);
    
    ydot[35] =  - (lambda11A + psi + (q11A_00A + q11A_01A + q11A_10A + q11A_11B + q11A_11C + q11A_11D + q11A_11E + q11A_11F + q11A_11G + q11A_11H) + mu11A) * D11A + 2*lambda11A*E11A*D11A + (q11A_00A*D00A + q11A_01A*D01A + q11A_10A*D10A + q11A_11B*D11B + q11A_11C*D11C + q11A_11D*D11D + q11A_11E*D11E + q11A_11F*D11F + q11A_11G*D11G + q11A_11H*D11H);
    
    ydot[36] =  - (lambda00B + psi + (q00B_00A + q00B_01B + q00B_10B + q00B_11B + q00B_00C + q00B_00D + q00B_00E + q00B_00F + q00B_00G + q00B_00H) + mu00B) * D00B + 2*lambda00B*E00B*D00B + (q00B_00A*D00A + q00B_01B*D01B + q00B_10B*D10B + q00B_11B*D11B + q00B_00C*D00C + q00B_00D*D00D + q00B_00E*D00E + q00B_00F*D00F + q00B_00G*D00G + q00B_00H*D00H);
    
    ydot[37] =  - (lambda01B + psi + (q01B_01A + q01B_00B + q01B_10B + q01B_11B + q01B_01C + q01B_01D + q01B_01E + q01B_01F + q01B_01G + q01B_01H) + mu01B) * D01B + 2*lambda01B*E01B*D01B + (q01B_01A*D01A + q01B_00B*D00B + q01B_10B*D10B + q01B_11B*D11B + q01B_01C*D01C + q01B_01D*D01D + q01B_01E*D01E + q01B_01F*D01F + q01B_01G*D01G + q01B_01H*D01H);
    
    ydot[38] =  - (lambda10B + psi + (q10B_10A + q10B_00B + q10B_01B + q10B_11B + q10B_10C + q10B_10D + q10B_10E + q10B_10F + q10B_10G + q10B_10H) + mu10B) * D10B + 2*lambda10B*E10B*D10B + (q10B_10A*D10A + q10B_00B*D00B + q10B_01B*D01B + q10B_11B*D11B + q10B_10C*D10C + q10B_10D*D10D + q10B_10E*D10E + q10B_10F*D10F + q10B_10G*D10G + q10B_10H*D10H);
    
    ydot[39] =  - (lambda11B + psi + (q11B_11A + q11B_00B + q11B_01B + q11B_10B + q11B_11C + q11B_11D + q11B_11E + q11B_11F + q11B_11G + q11B_11H) + mu11B) * D11B + 2*lambda11B*E11B*D11B + (q11B_11A*D11A + q11B_00B*D00B + q11B_01B*D01B + q11B_10B*D10B + q11B_11C*D11C + q11B_11D*D11D + q11B_11E*D11E + q11B_11F*D11F + q11B_11G*D11G + q11B_11H*D11H);
    
    ydot[40] =  - (lambda00C + psi + (q00C_00A + q00C_00B + q00C_01C + q00C_10C + q00C_11C + q00C_00D + q00C_00E + q00C_00F + q00C_00G + q00C_00H) + mu00C) * D00C + 2*lambda00C*E00C*D00C + (q00C_00A*D00A + q00C_00B*D00B + q00C_01C*D01C + q00C_10C*D10C + q00C_11C*D11C + q00C_00D*D00D + q00C_00E*D00E + q00C_00F*D00F + q00C_00G*D00G + q00C_00H*D00H);
    
    ydot[41] =  - (lambda01C + psi + (q01C_01A + q01C_01B + q01C_00C + q01C_10C + q01C_11C + q01C_01D + q01C_01E + q01C_01F + q01C_01G + q01C_01H) + mu01C) * D01C + 2*lambda01C*E01C*D01C + (q01C_01A*D01A + q01C_01B*D01B + q01C_00C*D00C + q01C_10C*D10C + q01C_11C*D11C + q01C_01D*D01D + q01C_01E*D01E + q01C_01F*D01F + q01C_01G*D01G + q01C_01H*D01H);
    
    ydot[42] =  - (lambda10C + psi + (q10C_10A + q10C_10B + q10C_00C + q10C_01C + q10C_11C + q10C_10D + q10C_10E + q10C_10F + q10C_10G + q10C_10H) + mu10C) * D10C + 2*lambda10C*E10C*D10C + (q10C_10A*D10A + q10C_10B*D10B + q10C_00C*D00C + q10C_01C*D01C + q10C_11C*D11C + q10C_10D*D10D + q10C_10E*D10E + q10C_10F*D10F + q10C_10G*D10G + q10C_10H*D10H);
    
    ydot[43] =  - (lambda11C + psi + (q11C_11A + q11C_11B + q11C_00C + q11C_01C + q11C_10C + q11C_11D + q11C_11E + q11C_11F + q11C_11G + q11C_11H) + mu11C) * D11C + 2*lambda11C*E11C*D11C + (q11C_11A*D11A + q11C_11B*D11B + q11C_00C*D00C + q11C_01C*D01C + q11C_10C*D10C + q11C_11D*D11D + q11C_11E*D11E + q11C_11F*D11F + q11C_11G*D11G + q11C_11H*D11H);
    
    ydot[44] =  - (lambda00D + psi + (q00D_00A + q00D_00B + q00D_00C + q00D_01D + q00D_10D + q00D_11D + q00D_00E + q00D_00F + q00D_00G + q00D_00H) + mu00D) * D00D + 2*lambda00D*E00D*D00D + (q00D_00A*D00A + q00D_00B*D00B + q00D_00C*D00C + q00D_01D*D01D + q00D_10D*D10D + q00D_11D*D11D + q00D_00E*D00E + q00D_00F*D00F + q00D_00G*D00G + q00D_00H*D00H);
    
    ydot[45] =  - (lambda01D + psi + (q01D_01A + q01D_01B + q01D_01C + q01D_00D + q01D_10D + q01D_11D + q01D_01E + q01D_01F + q01D_01G + q01D_01H) + mu01D) * D01D + 2*lambda01D*E01D*D01D + (q01D_01A*D01A + q01D_01B*D01B + q01D_01C*D01C + q01D_00D*D00D + q01D_10D*D10D + q01D_11D*D11D + q01D_01E*D01E + q01D_01F*D01F + q01D_01G*D01G + q01D_01H*D01H);
    
    ydot[46] =  - (lambda10D + psi + (q10D_10A + q10D_10B + q10D_10C + q10D_00D + q10D_01D + q10D_11D + q10D_10E + q10D_10F + q10D_10G + q10D_10H) + mu10D) * D10D + 2*lambda10D*E10D*D10D + (q10D_10A*D10A + q10D_10B*D10B + q10D_10C*D10C + q10D_00D*D00D + q10D_01D*D01D + q10D_11D*D11D + q10D_10E*D10E + q10D_10F*D10F + q10D_10G*D10G + q10D_10H*D10H);
    
    ydot[47] =  - (lambda11D + psi + (q11D_11A + q11D_11B + q11D_11C + q11D_00D + q11D_01D + q11D_10D + q11D_11E + q11D_11F + q11D_11G + q11D_11H) + mu11D) * D11D + 2*lambda11D*E11D*D11D + (q11D_11A*D11A + q11D_11B*D11B + q11D_11C*D11C + q11D_00D*D00D + q11D_01D*D01D + q11D_10D*D10D + q11D_11E*D11E + q11D_11F*D11F + q11D_11G*D11G + q11D_11H*D11H);
    
    ydot[48] =  - (lambda00E + psi + (q00E_00A + q00E_00B + q00E_00C + q00E_00D + q00E_01E + q00E_10E + q00E_11E + q00E_00F + q00E_00G + q00E_00H) + mu00E) * D00E + 2*lambda00E*E00E*D00E + (q00E_00A*D00A + q00E_00B*D00B + q00E_00C*D00C + q00E_00D*D00D + q00E_01E*D01E + q00E_10E*D10E + q00E_11E*D11E + q00E_00F*D00F + q00E_00G*D00G + q00E_00H*D00H);
    
    ydot[49] =  - (lambda01E + psi + (q01E_01A + q01E_01B + q01E_01C + q01E_01D + q01E_00E + q01E_10E + q01E_11E + q01E_01F + q01E_01G + q01E_01H) + mu01E) * D01E + 2*lambda01E*E01E*D01E + (q01E_01A*D01A + q01E_01B*D01B + q01E_01C*D01C + q01E_01D*D01D + q01E_00E*D00E + q01E_10E*D10E + q01E_11E*D11E + q01E_01F*D01F + q01E_01G*D01G + q01E_01H*D01H);
    
    ydot[50] =  - (lambda10E + psi + (q10E_10A + q10E_10B + q10E_10C + q10E_10D + q10E_00E + q10E_01E + q10E_11E + q10E_10F + q10E_10G + q10E_10H) + mu10E) * D10E + 2*lambda10E*E10E*D10E + (q10E_10A*D10A + q10E_10B*D10B + q10E_10C*D10C + q10E_10D*D10D + q10E_00E*D00E + q10E_01E*D01E + q10E_11E*D11E + q10E_10F*D10F + q10E_10G*D10G + q10E_10H*D10H);
    
    ydot[51] =  - (lambda11E + psi + (q11E_11A + q11E_11B + q11E_11C + q11E_11D + q11E_00E + q11E_01E + q11E_10E + q11E_11F + q11E_11G + q11E_11H) + mu11E) * D11E + 2*lambda11E*E11E*D11E + (q11E_11A*D11A + q11E_11B*D11B + q11E_11C*D11C + q11E_11D*D11D + q11E_00E*D00E + q11E_01E*D01E + q11E_10E*D10E + q11E_11F*D11F + q11E_11G*D11G + q11E_11H*D11H);
    
    ydot[52] =  - (lambda00F + psi + (q00F_00A + q00F_00B + q00F_00C + q00F_00D + q00F_00E + q00F_01F + q00F_10F + q00F_11F + q00F_00G + q00F_00H) + mu00F) * D00F + 2*lambda00F*E00F*D00F + (q00F_00A*D00A + q00F_00B*D00B + q00F_00C*D00C + q00F_00D*D00D + q00F_00E*D00E + q00F_01F*D01F + q00F_10F*D10F + q00F_11F*D11F + q00F_00G*D00G + q00F_00H*D00H);
    
    ydot[53] =  - (lambda01F + psi + (q01F_01A + q01F_01B + q01F_01C + q01F_01D + q01F_01E + q01F_00F + q01F_10F + q01F_11F + q01F_01G + q01F_01H) + mu01F) * D01F + 2*lambda01F*E01F*D01F + (q01F_01A*D01A + q01F_01B*D01B + q01F_01C*D01C + q01F_01D*D01D + q01F_01E*D01E + q01F_00F*D00F + q01F_10F*D10F + q01F_11F*D11F + q01F_01G*D01G + q01F_01H*D01H);
    
    ydot[54] =  - (lambda10F + psi + (q10F_10A + q10F_10B + q10F_10C + q10F_10D + q10F_10E + q10F_00F + q10F_01F + q10F_11F + q10F_10G + q10F_10H) + mu10F) * D10F + 2*lambda10F*E10F*D10F + (q10F_10A*D10A + q10F_10B*D10B + q10F_10C*D10C + q10F_10D*D10D + q10F_10E*D10E + q10F_00F*D00F + q10F_01F*D01F + q10F_11F*D11F + q10F_10G*D10G + q10F_10H*D10H);
    
    ydot[55] =  - (lambda11F + psi + (q11F_11A + q11F_11B + q11F_11C + q11F_11D + q11F_11E + q11F_00F + q11F_01F + q11F_10F + q11F_11G + q11F_11H) + mu11F) * D11F + 2*lambda11F*E11F*D11F + (q11F_11A*D11A + q11F_11B*D11B + q11F_11C*D11C + q11F_11D*D11D + q11F_11E*D11E + q11F_00F*D00F + q11F_01F*D01F + q11F_10F*D10F + q11F_11G*D11G + q11F_11H*D11H);
    
    ydot[56] =  - (lambda00G + psi + (q00G_00A + q00G_00B + q00G_00C + q00G_00D + q00G_00E + q00G_00F + q00G_01G + q00G_10G + q00G_11G + q00G_00H) + mu00G) * D00G + 2*lambda00G*E00G*D00G + (q00G_00A*D00A + q00G_00B*D00B + q00G_00C*D00C + q00G_00D*D00D + q00G_00E*D00E + q00G_00F*D00F + q00G_01G*D01G + q00G_10G*D10G + q00G_11G*D11G + q00G_00H*D00H);
    
    ydot[57] =  - (lambda01G + psi + (q01G_01A + q01G_01B + q01G_01C + q01G_01D + q01G_01E + q01G_01F + q01G_00G + q01G_10G + q01G_11G + q01G_01H) + mu01G) * D01G + 2*lambda01G*E01G*D01G + (q01G_01A*D01A + q01G_01B*D01B + q01G_01C*D01C + q01G_01D*D01D + q01G_01E*D01E + q01G_01F*D01F + q01G_00G*D00G + q01G_10G*D10G + q01G_11G*D11G + q01G_01H*D01H);
    
    ydot[58] =  - (lambda10G + psi + (q10G_10A + q10G_10B + q10G_10C + q10G_10D + q10G_10E + q10G_10F + q10G_00G + q10G_01G + q10G_11G + q10G_10H) + mu10G) * D10G + 2*lambda10G*E10G*D10G + (q10G_10A*D10A + q10G_10B*D10B + q10G_10C*D10C + q10G_10D*D10D + q10G_10E*D10E + q10G_10F*D10F + q10G_00G*D00G + q10G_01G*D01G + q10G_11G*D11G + q10G_10H*D10H);
    
    ydot[59] =  - (lambda11G + psi + (q11G_11A + q11G_11B + q11G_11C + q11G_11D + q11G_11E + q11G_11F + q11G_00G + q11G_01G + q11G_10G + q11G_11H) + mu11G) * D11G + 2*lambda11G*E11G*D11G + (q11G_11A*D11A + q11G_11B*D11B + q11G_11C*D11C + q11G_11D*D11D + q11G_11E*D11E + q11G_11F*D11F + q11G_00G*D00G + q11G_01G*D01G + q11G_10G*D10G + q11G_11H*D11H);
    
    ydot[60] =  - (lambda00H + psi + (q00H_00A + q00H_00B + q00H_00C + q00H_00D + q00H_00E + q00H_00F + q00H_00G + q00H_01H + q00H_10H + q00H_11H) + mu00H) * D00H + 2*lambda00H*E00H*D00H + (q00H_00A*D00A + q00H_00B*D00B + q00H_00C*D00C + q00H_00D*D00D + q00H_00E*D00E + q00H_00F*D00F + q00H_00G*D00G + q00H_01H*D01H + q00H_10H*D10H + q00H_11H*D11H);
    
    ydot[61] =  - (lambda01H + psi + (q01H_01A + q01H_01B + q01H_01C + q01H_01D + q01H_01E + q01H_01F + q01H_01G + q01H_00H + q01H_10H + q01H_11H) + mu01H) * D01H + 2*lambda01H*E01H*D01H + (q01H_01A*D01A + q01H_01B*D01B + q01H_01C*D01C + q01H_01D*D01D + q01H_01E*D01E + q01H_01F*D01F + q01H_01G*D01G + q01H_00H*D00H + q01H_10H*D10H + q01H_11H*D11H);
    
    ydot[62] =  - (lambda10H + psi + (q10H_10A + q10H_10B + q10H_10C + q10H_10D + q10H_10E + q10H_10F + q10H_10G + q10H_00H + q10H_01H + q10H_11H) + mu10H) * D10H + 2*lambda10H*E10H*D10H + (q10H_10A*D10A + q10H_10B*D10B + q10H_10C*D10C + q10H_10D*D10D + q10H_10E*D10E + q10H_10F*D10F + q10H_10G*D10G + q10H_00H*D00H + q10H_01H*D01H + q10H_11H*D11H);
    
    ydot[63] =  - (lambda11H + psi + (q11H_11A + q11H_11B + q11H_11C + q11H_11D + q11H_11E + q11H_11F + q11H_11G + q11H_00H + q11H_01H + q11H_10H) + mu11H) * D11H + 2*lambda11H*E11H*D11H + (q11H_11A*D11A + q11H_11B*D11B + q11H_11C*D11C + q11H_11D*D11D + q11H_11E*D11E + q11H_11F*D11F + q11H_11G*D11G + q11H_00H*D00H + q11H_01H*D01H + q11H_10H*D10H);
    
}





void muhisse_strat_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
    double E00A = y[0];
    double E01A = y[1];
    double E10A = y[2];
    double E11A = y[3];
    double E00B = y[4];
    double E01B = y[5];
    double E10B = y[6];
    double E11B = y[7];
    double E00C = y[8];
    double E01C = y[9];
    double E10C = y[10];
    double E11C = y[11];
    double E00D = y[12];
    double E01D = y[13];
    double E10D = y[14];
    double E11D = y[15];
    double E00E = y[16];
    double E01E = y[17];
    double E10E = y[18];
    double E11E = y[19];
    double E00F = y[20];
    double E01F = y[21];
    double E10F = y[22];
    double E11F = y[23];
    double E00G = y[24];
    double E01G = y[25];
    double E10G = y[26];
    double E11G = y[27];
    double E00H = y[28];
    double E01H = y[29];
    double E10H = y[30];
    double E11H = y[31];
    
    double D00A = y[32];
    double D01A = y[33];
    double D10A = y[34];
    double D11A = y[35];
    double D00B = y[36];
    double D01B = y[37];
    double D10B = y[38];
    double D11B = y[39];
    double D00C = y[40];
    double D01C = y[41];
    double D10C = y[42];
    double D11C = y[43];
    double D00D = y[44];
    double D01D = y[45];
    double D10D = y[46];
    double D11D = y[47];
    double D00E = y[48];
    double D01E = y[49];
    double D10E = y[50];
    double D11E = y[51];
    double D00F = y[52];
    double D01F = y[53];
    double D10F = y[54];
    double D11F = y[55];
    double D00G = y[56];
    double D01G = y[57];
    double D10G = y[58];
    double D11G = y[59];
    double D00H = y[60];
    double D01H = y[61];
    double D10H = y[62];
    double D11H = y[63];
    
    
    double
    lambda00A = params_muhisse[0],
    lambda01A = params_muhisse[1],
    lambda10A = params_muhisse[2],
    lambda11A = params_muhisse[3],
    mu00A = params_muhisse[4],
    mu01A = params_muhisse[5],
    mu10A = params_muhisse[6],
    mu11A = params_muhisse[7],
    q00A_01A = params_muhisse[8],
    q00A_10A = params_muhisse[9],
    q00A_11A = params_muhisse[10],
    q01A_00A = params_muhisse[11],
    q01A_10A = params_muhisse[12],
    q01A_11A = params_muhisse[13],
    q10A_00A = params_muhisse[14],
    q10A_01A = params_muhisse[15],
    q10A_11A = params_muhisse[16],
    q11A_00A = params_muhisse[17],
    q11A_01A = params_muhisse[18],
    q11A_10A = params_muhisse[19],
    q00A_00B = params_muhisse[20],
    q00A_00C = params_muhisse[21],
    q00A_00D = params_muhisse[22],
    q00A_00E = params_muhisse[23],
    q00A_00F = params_muhisse[24],
    q00A_00G = params_muhisse[25],
    q00A_00H = params_muhisse[26],
    q01A_01B = params_muhisse[27],
    q01A_01C = params_muhisse[28],
    q01A_01D = params_muhisse[29],
    q01A_01E = params_muhisse[30],
    q01A_01F = params_muhisse[31],
    q01A_01G = params_muhisse[32],
    q01A_01H = params_muhisse[33],
    q10A_10B = params_muhisse[34],
    q10A_10C = params_muhisse[35],
    q10A_10D = params_muhisse[36],
    q10A_10E = params_muhisse[37],
    q10A_10F = params_muhisse[38],
    q10A_10G = params_muhisse[39],
    q10A_10H = params_muhisse[40],
    q11A_11B = params_muhisse[41],
    q11A_11C = params_muhisse[42],
    q11A_11D = params_muhisse[43],
    q11A_11E = params_muhisse[44],
    q11A_11F = params_muhisse[45],
    q11A_11G = params_muhisse[46],
    q11A_11H = params_muhisse[47],
    lambda00B = params_muhisse[48],
    lambda01B = params_muhisse[49],
    lambda10B = params_muhisse[50],
    lambda11B = params_muhisse[51],
    mu00B = params_muhisse[52],
    mu01B = params_muhisse[53],
    mu10B = params_muhisse[54],
    mu11B = params_muhisse[55],
    q00B_01B = params_muhisse[56],
    q00B_10B = params_muhisse[57],
    q00B_11B = params_muhisse[58],
    q01B_00B = params_muhisse[59],
    q01B_10B = params_muhisse[60],
    q01B_11B = params_muhisse[61],
    q10B_00B = params_muhisse[62],
    q10B_01B = params_muhisse[63],
    q10B_11B = params_muhisse[64],
    q11B_00B = params_muhisse[65],
    q11B_01B = params_muhisse[66],
    q11B_10B = params_muhisse[67],
    q00B_00A = params_muhisse[68],
    q00B_00C = params_muhisse[69],
    q00B_00D = params_muhisse[70],
    q00B_00E = params_muhisse[71],
    q00B_00F = params_muhisse[72],
    q00B_00G = params_muhisse[73],
    q00B_00H = params_muhisse[74],
    q01B_01A = params_muhisse[75],
    q01B_01C = params_muhisse[76],
    q01B_01D = params_muhisse[77],
    q01B_01E = params_muhisse[78],
    q01B_01F = params_muhisse[79],
    q01B_01G = params_muhisse[80],
    q01B_01H = params_muhisse[81],
    q10B_10A = params_muhisse[82],
    q10B_10C = params_muhisse[83],
    q10B_10D = params_muhisse[84],
    q10B_10E = params_muhisse[85],
    q10B_10F = params_muhisse[86],
    q10B_10G = params_muhisse[87],
    q10B_10H = params_muhisse[88],
    q11B_11A = params_muhisse[89],
    q11B_11C = params_muhisse[90],
    q11B_11D = params_muhisse[91],
    q11B_11E = params_muhisse[92],
    q11B_11F = params_muhisse[93],
    q11B_11G = params_muhisse[94],
    q11B_11H = params_muhisse[95],
    lambda00C = params_muhisse[96],
    lambda01C = params_muhisse[97],
    lambda10C = params_muhisse[98],
    lambda11C = params_muhisse[99],
    mu00C = params_muhisse[100],
    mu01C = params_muhisse[101],
    mu10C = params_muhisse[102],
    mu11C = params_muhisse[103],
    q00C_01C = params_muhisse[104],
    q00C_10C = params_muhisse[105],
    q00C_11C = params_muhisse[106],
    q01C_00C = params_muhisse[107],
    q01C_10C = params_muhisse[108],
    q01C_11C = params_muhisse[109],
    q10C_00C = params_muhisse[110],
    q10C_01C = params_muhisse[111],
    q10C_11C = params_muhisse[112],
    q11C_00C = params_muhisse[113],
    q11C_01C = params_muhisse[114],
    q11C_10C = params_muhisse[115],
    q00C_00A = params_muhisse[116],
    q00C_00B = params_muhisse[117],
    q00C_00D = params_muhisse[118],
    q00C_00E = params_muhisse[119],
    q00C_00F = params_muhisse[120],
    q00C_00G = params_muhisse[121],
    q00C_00H = params_muhisse[122],
    q01C_01A = params_muhisse[123],
    q01C_01B = params_muhisse[124],
    q01C_01D = params_muhisse[125],
    q01C_01E = params_muhisse[126],
    q01C_01F = params_muhisse[127],
    q01C_01G = params_muhisse[128],
    q01C_01H = params_muhisse[129],
    q10C_10A = params_muhisse[130],
    q10C_10B = params_muhisse[131],
    q10C_10D = params_muhisse[132],
    q10C_10E = params_muhisse[133],
    q10C_10F = params_muhisse[134],
    q10C_10G = params_muhisse[135],
    q10C_10H = params_muhisse[136],
    q11C_11A = params_muhisse[137],
    q11C_11B = params_muhisse[138],
    q11C_11D = params_muhisse[139],
    q11C_11E = params_muhisse[140],
    q11C_11F = params_muhisse[141],
    q11C_11G = params_muhisse[142],
    q11C_11H = params_muhisse[143],
    lambda00D = params_muhisse[144],
    lambda01D = params_muhisse[145],
    lambda10D = params_muhisse[146],
    lambda11D = params_muhisse[147],
    mu00D = params_muhisse[148],
    mu01D = params_muhisse[149],
    mu10D = params_muhisse[150],
    mu11D = params_muhisse[151],
    q00D_01D = params_muhisse[152],
    q00D_10D = params_muhisse[153],
    q00D_11D = params_muhisse[154],
    q01D_00D = params_muhisse[155],
    q01D_10D = params_muhisse[156],
    q01D_11D = params_muhisse[157],
    q10D_00D = params_muhisse[158],
    q10D_01D = params_muhisse[159],
    q10D_11D = params_muhisse[160],
    q11D_00D = params_muhisse[161],
    q11D_01D = params_muhisse[162],
    q11D_10D = params_muhisse[163],
    q00D_00A = params_muhisse[164],
    q00D_00B = params_muhisse[165],
    q00D_00C = params_muhisse[166],
    q00D_00E = params_muhisse[167],
    q00D_00F = params_muhisse[168],
    q00D_00G = params_muhisse[169],
    q00D_00H = params_muhisse[170],
    q01D_01A = params_muhisse[171],
    q01D_01B = params_muhisse[172],
    q01D_01C = params_muhisse[173],
    q01D_01E = params_muhisse[174],
    q01D_01F = params_muhisse[175],
    q01D_01G = params_muhisse[176],
    q01D_01H = params_muhisse[177],
    q10D_10A = params_muhisse[178],
    q10D_10B = params_muhisse[179],
    q10D_10C = params_muhisse[180],
    q10D_10E = params_muhisse[181],
    q10D_10F = params_muhisse[182],
    q10D_10G = params_muhisse[183],
    q10D_10H = params_muhisse[184],
    q11D_11A = params_muhisse[185],
    q11D_11B = params_muhisse[186],
    q11D_11C = params_muhisse[187],
    q11D_11E = params_muhisse[188],
    q11D_11F = params_muhisse[189],
    q11D_11G = params_muhisse[190],
    q11D_11H = params_muhisse[191],
    lambda00E = params_muhisse[192],
    lambda01E = params_muhisse[193],
    lambda10E = params_muhisse[194],
    lambda11E = params_muhisse[195],
    mu00E = params_muhisse[196],
    mu01E = params_muhisse[197],
    mu10E = params_muhisse[198],
    mu11E = params_muhisse[199],
    q00E_01E = params_muhisse[200],
    q00E_10E = params_muhisse[201],
    q00E_11E = params_muhisse[202],
    q01E_00E = params_muhisse[203],
    q01E_10E = params_muhisse[204],
    q01E_11E = params_muhisse[205],
    q10E_00E = params_muhisse[206],
    q10E_01E = params_muhisse[207],
    q10E_11E = params_muhisse[208],
    q11E_00E = params_muhisse[209],
    q11E_01E = params_muhisse[210],
    q11E_10E = params_muhisse[211],
    q00E_00A = params_muhisse[212],
    q00E_00B = params_muhisse[213],
    q00E_00C = params_muhisse[214],
    q00E_00D = params_muhisse[215],
    q00E_00F = params_muhisse[216],
    q00E_00G = params_muhisse[217],
    q00E_00H = params_muhisse[218],
    q01E_01A = params_muhisse[219],
    q01E_01B = params_muhisse[220],
    q01E_01C = params_muhisse[221],
    q01E_01D = params_muhisse[222],
    q01E_01F = params_muhisse[223],
    q01E_01G = params_muhisse[224],
    q01E_01H = params_muhisse[225],
    q10E_10A = params_muhisse[226],
    q10E_10B = params_muhisse[227],
    q10E_10C = params_muhisse[228],
    q10E_10D = params_muhisse[229],
    q10E_10F = params_muhisse[230],
    q10E_10G = params_muhisse[231],
    q10E_10H = params_muhisse[232],
    q11E_11A = params_muhisse[233],
    q11E_11B = params_muhisse[234],
    q11E_11C = params_muhisse[235],
    q11E_11D = params_muhisse[236],
    q11E_11F = params_muhisse[237],
    q11E_11G = params_muhisse[238],
    q11E_11H = params_muhisse[239],
    lambda00F = params_muhisse[240],
    lambda01F = params_muhisse[241],
    lambda10F = params_muhisse[242],
    lambda11F = params_muhisse[243],
    mu00F = params_muhisse[244],
    mu01F = params_muhisse[245],
    mu10F = params_muhisse[246],
    mu11F = params_muhisse[247],
    q00F_01F = params_muhisse[248],
    q00F_10F = params_muhisse[249],
    q00F_11F = params_muhisse[250],
    q01F_00F = params_muhisse[251],
    q01F_10F = params_muhisse[252],
    q01F_11F = params_muhisse[253],
    q10F_00F = params_muhisse[254],
    q10F_01F = params_muhisse[255],
    q10F_11F = params_muhisse[256],
    q11F_00F = params_muhisse[257],
    q11F_01F = params_muhisse[258],
    q11F_10F = params_muhisse[259],
    q00F_00A = params_muhisse[260],
    q00F_00B = params_muhisse[261],
    q00F_00C = params_muhisse[262],
    q00F_00D = params_muhisse[263],
    q00F_00E = params_muhisse[264],
    q00F_00G = params_muhisse[265],
    q00F_00H = params_muhisse[266],
    q01F_01A = params_muhisse[267],
    q01F_01B = params_muhisse[268],
    q01F_01C = params_muhisse[269],
    q01F_01D = params_muhisse[270],
    q01F_01E = params_muhisse[271],
    q01F_01G = params_muhisse[272],
    q01F_01H = params_muhisse[273],
    q10F_10A = params_muhisse[274],
    q10F_10B = params_muhisse[275],
    q10F_10C = params_muhisse[276],
    q10F_10D = params_muhisse[277],
    q10F_10E = params_muhisse[278],
    q10F_10G = params_muhisse[279],
    q10F_10H = params_muhisse[280],
    q11F_11A = params_muhisse[281],
    q11F_11B = params_muhisse[282],
    q11F_11C = params_muhisse[283],
    q11F_11D = params_muhisse[284],
    q11F_11E = params_muhisse[285],
    q11F_11G = params_muhisse[286],
    q11F_11H = params_muhisse[287],
    lambda00G = params_muhisse[288],
    lambda01G = params_muhisse[289],
    lambda10G = params_muhisse[290],
    lambda11G = params_muhisse[291],
    mu00G = params_muhisse[292],
    mu01G = params_muhisse[293],
    mu10G = params_muhisse[294],
    mu11G = params_muhisse[295],
    q00G_01G = params_muhisse[296],
    q00G_10G = params_muhisse[297],
    q00G_11G = params_muhisse[298],
    q01G_00G = params_muhisse[299],
    q01G_10G = params_muhisse[300],
    q01G_11G = params_muhisse[301],
    q10G_00G = params_muhisse[302],
    q10G_01G = params_muhisse[303],
    q10G_11G = params_muhisse[304],
    q11G_00G = params_muhisse[305],
    q11G_01G = params_muhisse[306],
    q11G_10G = params_muhisse[307],
    q00G_00A = params_muhisse[308],
    q00G_00B = params_muhisse[309],
    q00G_00C = params_muhisse[310],
    q00G_00D = params_muhisse[311],
    q00G_00E = params_muhisse[312],
    q00G_00F = params_muhisse[313],
    q00G_00H = params_muhisse[314],
    q01G_01A = params_muhisse[315],
    q01G_01B = params_muhisse[316],
    q01G_01C = params_muhisse[317],
    q01G_01D = params_muhisse[318],
    q01G_01E = params_muhisse[319],
    q01G_01F = params_muhisse[320],
    q01G_01H = params_muhisse[321],
    q10G_10A = params_muhisse[322],
    q10G_10B = params_muhisse[323],
    q10G_10C = params_muhisse[324],
    q10G_10D = params_muhisse[325],
    q10G_10E = params_muhisse[326],
    q10G_10F = params_muhisse[327],
    q10G_10H = params_muhisse[328],
    q11G_11A = params_muhisse[329],
    q11G_11B = params_muhisse[330],
    q11G_11C = params_muhisse[331],
    q11G_11D = params_muhisse[332],
    q11G_11E = params_muhisse[333],
    q11G_11F = params_muhisse[334],
    q11G_11H = params_muhisse[335],
    lambda00H = params_muhisse[336],
    lambda01H = params_muhisse[337],
    lambda10H = params_muhisse[338],
    lambda11H = params_muhisse[339],
    mu00H = params_muhisse[340],
    mu01H = params_muhisse[341],
    mu10H = params_muhisse[342],
    mu11H = params_muhisse[343],
    q00H_01H = params_muhisse[344],
    q00H_10H = params_muhisse[345],
    q00H_11H = params_muhisse[346],
    q01H_00H = params_muhisse[347],
    q01H_10H = params_muhisse[348],
    q01H_11H = params_muhisse[349],
    q10H_00H = params_muhisse[350],
    q10H_01H = params_muhisse[351],
    q10H_11H = params_muhisse[352],
    q11H_00H = params_muhisse[353],
    q11H_01H = params_muhisse[354],
    q11H_10H = params_muhisse[355],
    q00H_00A = params_muhisse[356],
    q00H_00B = params_muhisse[357],
    q00H_00C = params_muhisse[358],
    q00H_00D = params_muhisse[359],
    q00H_00E = params_muhisse[360],
    q00H_00F = params_muhisse[361],
    q00H_00G = params_muhisse[362],
    q01H_01A = params_muhisse[363],
    q01H_01B = params_muhisse[364],
    q01H_01C = params_muhisse[365],
    q01H_01D = params_muhisse[366],
    q01H_01E = params_muhisse[367],
    q01H_01F = params_muhisse[368],
    q01H_01G = params_muhisse[369],
    q10H_10A = params_muhisse[370],
    q10H_10B = params_muhisse[371],
    q10H_10C = params_muhisse[372],
    q10H_10D = params_muhisse[373],
    q10H_10E = params_muhisse[374],
    q10H_10F = params_muhisse[375],
    q10H_10G = params_muhisse[376],
    q11H_11A = params_muhisse[377],
    q11H_11B = params_muhisse[378],
    q11H_11C = params_muhisse[379],
    q11H_11D = params_muhisse[380],
    q11H_11E = params_muhisse[381],
    q11H_11F = params_muhisse[382],
    q11H_11G = params_muhisse[383],
    psi = params_muhisse[384];
    
    /* The E's */
    ydot[0] = mu00A - (lambda00A + psi + (q00A_01A + q00A_10A + q00A_11A + q00A_00B + q00A_00C + q00A_00D + q00A_00E + q00A_00F + q00A_00G + q00A_00H) + mu00A) * E00A + lambda00A*E00A*E00A + (q00A_01A*E01A + q00A_10A*E10A + q00A_11A*E11A + q00A_00B*E00B + q00A_00C*E00C + q00A_00D*E00D + q00A_00E*E00E + q00A_00F*E00F + q00A_00G*E00G + q00A_00H*E00H);
    
    ydot[1] = mu01A - (lambda01A + psi + (q01A_00A + q01A_10A + q01A_11A + q01A_01B + q01A_01C + q01A_01D + q01A_01E + q01A_01F + q01A_01G + q01A_01H) + mu01A) * E01A + lambda01A*E01A*E01A + (q01A_00A*E00A + q01A_10A*E10A + q01A_11A*E11A + q01A_01B*E01B + q01A_01C*E01C + q01A_01D*E01D + q01A_01E*E01E + q01A_01F*E01F + q01A_01G*E01G + q01A_01H*E01H);
    
    ydot[2] = mu10A - (lambda10A + psi + (q10A_00A + q10A_01A + q10A_11A + q10A_10B + q10A_10C + q10A_10D + q10A_10E + q10A_10F + q10A_10G + q10A_10H) + mu10A) * E10A + lambda10A*E10A*E10A + (q10A_00A*E00A + q10A_01A*E01A + q10A_11A*E11A + q10A_10B*E10B + q10A_10C*E10C + q10A_10D*E10D + q10A_10E*E10E + q10A_10F*E10F + q10A_10G*E10G + q10A_10H*E10H);
    
    ydot[3] = mu11A - (lambda11A + psi + (q11A_00A + q11A_01A + q11A_10A + q11A_11B + q11A_11C + q11A_11D + q11A_11E + q11A_11F + q11A_11G + q11A_11H) + mu11A) * E11A + lambda11A*E11A*E11A + (q11A_00A*E00A + q11A_01A*E01A + q11A_10A*E10A + q11A_11B*E11B + q11A_11C*E11C + q11A_11D*E11D + q11A_11E*E11E + q11A_11F*E11F + q11A_11G*E11G + q11A_11H*E11H);
    
    ydot[4] = mu00B - (lambda00B + psi + (q00B_00A + q00B_01B + q00B_10B + q00B_11B + q00B_00C + q00B_00D + q00B_00E + q00B_00F + q00B_00G + q00B_00H) + mu00B) * E00B + lambda00B*E00B*E00B + (q00B_00A*E00A + q00B_01B*E01B + q00B_10B*E10B + q00B_11B*E11B + q00B_00C*E00C + q00B_00D*E00D + q00B_00E*E00E + q00B_00F*E00F + q00B_00G*E00G + q00B_00H*E00H);
    
    ydot[5] = mu01B - (lambda01B + psi + (q01B_01A + q01B_00B + q01B_10B + q01B_11B + q01B_01C + q01B_01D + q01B_01E + q01B_01F + q01B_01G + q01B_01H) + mu01B) * E01B + lambda01B*E01B*E01B + (q01B_01A*E01A + q01B_00B*E00B + q01B_10B*E10B + q01B_11B*E11B + q01B_01C*E01C + q01B_01D*E01D + q01B_01E*E01E + q01B_01F*E01F + q01B_01G*E01G + q01B_01H*E01H);
    
    ydot[6] = mu10B - (lambda10B + psi + (q10B_10A + q10B_00B + q10B_01B + q10B_11B + q10B_10C + q10B_10D + q10B_10E + q10B_10F + q10B_10G + q10B_10H) + mu10B) * E10B + lambda10B*E10B*E10B + (q10B_10A*E10A + q10B_00B*E00B + q10B_01B*E01B + q10B_11B*E11B + q10B_10C*E10C + q10B_10D*E10D + q10B_10E*E10E + q10B_10F*E10F + q10B_10G*E10G + q10B_10H*E10H);
    
    ydot[7] = mu11B - (lambda11B + psi + (q11B_11A + q11B_00B + q11B_01B + q11B_10B + q11B_11C + q11B_11D + q11B_11E + q11B_11F + q11B_11G + q11B_11H) + mu11B) * E11B + lambda11B*E11B*E11B + (q11B_11A*E11A + q11B_00B*E00B + q11B_01B*E01B + q11B_10B*E10B + q11B_11C*E11C + q11B_11D*E11D + q11B_11E*E11E + q11B_11F*E11F + q11B_11G*E11G + q11B_11H*E11H);
    
    ydot[8] = mu00C - (lambda00C + psi + (q00C_00A + q00C_00B + q00C_01C + q00C_10C + q00C_11C + q00C_00D + q00C_00E + q00C_00F + q00C_00G + q00C_00H) + mu00C) * E00C + lambda00C*E00C*E00C + (q00C_00A*E00A + q00C_00B*E00B + q00C_01C*E01C + q00C_10C*E10C + q00C_11C*E11C + q00C_00D*E00D + q00C_00E*E00E + q00C_00F*E00F + q00C_00G*E00G + q00C_00H*E00H);
    
    ydot[9] = mu01C - (lambda01C + psi + (q01C_01A + q01C_01B + q01C_00C + q01C_10C + q01C_11C + q01C_01D + q01C_01E + q01C_01F + q01C_01G + q01C_01H) + mu01C) * E01C + lambda01C*E01C*E01C + (q01C_01A*E01A + q01C_01B*E01B + q01C_00C*E00C + q01C_10C*E10C + q01C_11C*E11C + q01C_01D*E01D + q01C_01E*E01E + q01C_01F*E01F + q01C_01G*E01G + q01C_01H*E01H);
    
    ydot[10] = mu10C - (lambda10C + psi + (q10C_10A + q10C_10B + q10C_00C + q10C_01C + q10C_11C + q10C_10D + q10C_10E + q10C_10F + q10C_10G + q10C_10H) + mu10C) * E10C + lambda10C*E10C*E10C + (q10C_10A*E10A + q10C_10B*E10B + q10C_00C*E00C + q10C_01C*E01C + q10C_11C*E11C + q10C_10D*E10D + q10C_10E*E10E + q10C_10F*E10F + q10C_10G*E10G + q10C_10H*E10H);
    
    ydot[11] = mu11C - (lambda11C + psi + (q11C_11A + q11C_11B + q11C_00C + q11C_01C + q11C_10C + q11C_11D + q11C_11E + q11C_11F + q11C_11G + q11C_11H) + mu11C) * E11C + lambda11C*E11C*E11C + (q11C_11A*E11A + q11C_11B*E11B + q11C_00C*E00C + q11C_01C*E01C + q11C_10C*E10C + q11C_11D*E11D + q11C_11E*E11E + q11C_11F*E11F + q11C_11G*E11G + q11C_11H*E11H);
    
    ydot[12] = mu00D - (lambda00D + psi + (q00D_00A + q00D_00B + q00D_00C + q00D_01D + q00D_10D + q00D_11D + q00D_00E + q00D_00F + q00D_00G + q00D_00H) + mu00D) * E00D + lambda00D*E00D*E00D + (q00D_00A*E00A + q00D_00B*E00B + q00D_00C*E00C + q00D_01D*E01D + q00D_10D*E10D + q00D_11D*E11D + q00D_00E*E00E + q00D_00F*E00F + q00D_00G*E00G + q00D_00H*E00H);
    
    ydot[13] = mu01D - (lambda01D + psi + (q01D_01A + q01D_01B + q01D_01C + q01D_00D + q01D_10D + q01D_11D + q01D_01E + q01D_01F + q01D_01G + q01D_01H) + mu01D) * E01D + lambda01D*E01D*E01D + (q01D_01A*E01A + q01D_01B*E01B + q01D_01C*E01C + q01D_00D*E00D + q01D_10D*E10D + q01D_11D*E11D + q01D_01E*E01E + q01D_01F*E01F + q01D_01G*E01G + q01D_01H*E01H);
    
    ydot[14] = mu10D - (lambda10D + psi + (q10D_10A + q10D_10B + q10D_10C + q10D_00D + q10D_01D + q10D_11D + q10D_10E + q10D_10F + q10D_10G + q10D_10H) + mu10D) * E10D + lambda10D*E10D*E10D + (q10D_10A*E10A + q10D_10B*E10B + q10D_10C*E10C + q10D_00D*E00D + q10D_01D*E01D + q10D_11D*E11D + q10D_10E*E10E + q10D_10F*E10F + q10D_10G*E10G + q10D_10H*E10H);
    
    ydot[15] = mu11D - (lambda11D + psi + (q11D_11A + q11D_11B + q11D_11C + q11D_00D + q11D_01D + q11D_10D + q11D_11E + q11D_11F + q11D_11G + q11D_11H) + mu11D) * E11D + lambda11D*E11D*E11D + (q11D_11A*E11A + q11D_11B*E11B + q11D_11C*E11C + q11D_00D*E00D + q11D_01D*E01D + q11D_10D*E10D + q11D_11E*E11E + q11D_11F*E11F + q11D_11G*E11G + q11D_11H*E11H);
    
    ydot[16] = mu00E - (lambda00E + psi + (q00E_00A + q00E_00B + q00E_00C + q00E_00D + q00E_01E + q00E_10E + q00E_11E + q00E_00F + q00E_00G + q00E_00H) + mu00E) * E00E + lambda00E*E00E*E00E + (q00E_00A*E00A + q00E_00B*E00B + q00E_00C*E00C + q00E_00D*E00D + q00E_01E*E01E + q00E_10E*E10E + q00E_11E*E11E + q00E_00F*E00F + q00E_00G*E00G + q00E_00H*E00H);
    
    ydot[17] = mu01E - (lambda01E + psi + (q01E_01A + q01E_01B + q01E_01C + q01E_01D + q01E_00E + q01E_10E + q01E_11E + q01E_01F + q01E_01G + q01E_01H) + mu01E) * E01E + lambda01E*E01E*E01E + (q01E_01A*E01A + q01E_01B*E01B + q01E_01C*E01C + q01E_01D*E01D + q01E_00E*E00E + q01E_10E*E10E + q01E_11E*E11E + q01E_01F*E01F + q01E_01G*E01G + q01E_01H*E01H);
    
    ydot[18] = mu10E - (lambda10E + psi + (q10E_10A + q10E_10B + q10E_10C + q10E_10D + q10E_00E + q10E_01E + q10E_11E + q10E_10F + q10E_10G + q10E_10H) + mu10E) * E10E + lambda10E*E10E*E10E + (q10E_10A*E10A + q10E_10B*E10B + q10E_10C*E10C + q10E_10D*E10D + q10E_00E*E00E + q10E_01E*E01E + q10E_11E*E11E + q10E_10F*E10F + q10E_10G*E10G + q10E_10H*E10H);
    
    ydot[19] = mu11E - (lambda11E + psi + (q11E_11A + q11E_11B + q11E_11C + q11E_11D + q11E_00E + q11E_01E + q11E_10E + q11E_11F + q11E_11G + q11E_11H) + mu11E) * E11E + lambda11E*E11E*E11E + (q11E_11A*E11A + q11E_11B*E11B + q11E_11C*E11C + q11E_11D*E11D + q11E_00E*E00E + q11E_01E*E01E + q11E_10E*E10E + q11E_11F*E11F + q11E_11G*E11G + q11E_11H*E11H);
    
    ydot[20] = mu00F - (lambda00F + psi + (q00F_00A + q00F_00B + q00F_00C + q00F_00D + q00F_00E + q00F_01F + q00F_10F + q00F_11F + q00F_00G + q00F_00H) + mu00F) * E00F + lambda00F*E00F*E00F + (q00F_00A*E00A + q00F_00B*E00B + q00F_00C*E00C + q00F_00D*E00D + q00F_00E*E00E + q00F_01F*E01F + q00F_10F*E10F + q00F_11F*E11F + q00F_00G*E00G + q00F_00H*E00H);
    
    ydot[21] = mu01F - (lambda01F + psi + (q01F_01A + q01F_01B + q01F_01C + q01F_01D + q01F_01E + q01F_00F + q01F_10F + q01F_11F + q01F_01G + q01F_01H) + mu01F) * E01F + lambda01F*E01F*E01F + (q01F_01A*E01A + q01F_01B*E01B + q01F_01C*E01C + q01F_01D*E01D + q01F_01E*E01E + q01F_00F*E00F + q01F_10F*E10F + q01F_11F*E11F + q01F_01G*E01G + q01F_01H*E01H);
    
    ydot[22] = mu10F - (lambda10F + psi + (q10F_10A + q10F_10B + q10F_10C + q10F_10D + q10F_10E + q10F_00F + q10F_01F + q10F_11F + q10F_10G + q10F_10H) + mu10F) * E10F + lambda10F*E10F*E10F + (q10F_10A*E10A + q10F_10B*E10B + q10F_10C*E10C + q10F_10D*E10D + q10F_10E*E10E + q10F_00F*E00F + q10F_01F*E01F + q10F_11F*E11F + q10F_10G*E10G + q10F_10H*E10H);
    
    ydot[23] = mu11F - (lambda11F + psi + (q11F_11A + q11F_11B + q11F_11C + q11F_11D + q11F_11E + q11F_00F + q11F_01F + q11F_10F + q11F_11G + q11F_11H) + mu11F) * E11F + lambda11F*E11F*E11F + (q11F_11A*E11A + q11F_11B*E11B + q11F_11C*E11C + q11F_11D*E11D + q11F_11E*E11E + q11F_00F*E00F + q11F_01F*E01F + q11F_10F*E10F + q11F_11G*E11G + q11F_11H*E11H);
    
    ydot[24] = mu00G - (lambda00G + psi + (q00G_00A + q00G_00B + q00G_00C + q00G_00D + q00G_00E + q00G_00F + q00G_01G + q00G_10G + q00G_11G + q00G_00H) + mu00G) * E00G + lambda00G*E00G*E00G + (q00G_00A*E00A + q00G_00B*E00B + q00G_00C*E00C + q00G_00D*E00D + q00G_00E*E00E + q00G_00F*E00F + q00G_01G*E01G + q00G_10G*E10G + q00G_11G*E11G + q00G_00H*E00H);
    
    ydot[25] = mu01G - (lambda01G + psi + (q01G_01A + q01G_01B + q01G_01C + q01G_01D + q01G_01E + q01G_01F + q01G_00G + q01G_10G + q01G_11G + q01G_01H) + mu01G) * E01G + lambda01G*E01G*E01G + (q01G_01A*E01A + q01G_01B*E01B + q01G_01C*E01C + q01G_01D*E01D + q01G_01E*E01E + q01G_01F*E01F + q01G_00G*E00G + q01G_10G*E10G + q01G_11G*E11G + q01G_01H*E01H);
    
    ydot[26] = mu10G - (lambda10G + psi + (q10G_10A + q10G_10B + q10G_10C + q10G_10D + q10G_10E + q10G_10F + q10G_00G + q10G_01G + q10G_11G + q10G_10H) + mu10G) * E10G + lambda10G*E10G*E10G + (q10G_10A*E10A + q10G_10B*E10B + q10G_10C*E10C + q10G_10D*E10D + q10G_10E*E10E + q10G_10F*E10F + q10G_00G*E00G + q10G_01G*E01G + q10G_11G*E11G + q10G_10H*E10H);
    
    ydot[27] = mu11G - (lambda11G + psi + (q11G_11A + q11G_11B + q11G_11C + q11G_11D + q11G_11E + q11G_11F + q11G_00G + q11G_01G + q11G_10G + q11G_11H) + mu11G) * E11G + lambda11G*E11G*E11G + (q11G_11A*E11A + q11G_11B*E11B + q11G_11C*E11C + q11G_11D*E11D + q11G_11E*E11E + q11G_11F*E11F + q11G_00G*E00G + q11G_01G*E01G + q11G_10G*E10G + q11G_11H*E11H);
    
    ydot[28] = mu00H - (lambda00H + psi + (q00H_00A + q00H_00B + q00H_00C + q00H_00D + q00H_00E + q00H_00F + q00H_00G + q00H_01H + q00H_10H + q00H_11H) + mu00H) * E00H + lambda00H*E00H*E00H + (q00H_00A*E00A + q00H_00B*E00B + q00H_00C*E00C + q00H_00D*E00D + q00H_00E*E00E + q00H_00F*E00F + q00H_00G*E00G + q00H_01H*E01H + q00H_10H*E10H + q00H_11H*E11H);
    
    ydot[29] = mu01H - (lambda01H + psi + (q01H_01A + q01H_01B + q01H_01C + q01H_01D + q01H_01E + q01H_01F + q01H_01G + q01H_00H + q01H_10H + q01H_11H) + mu01H) * E01H + lambda01H*E01H*E01H + (q01H_01A*E01A + q01H_01B*E01B + q01H_01C*E01C + q01H_01D*E01D + q01H_01E*E01E + q01H_01F*E01F + q01H_01G*E01G + q01H_00H*E00H + q01H_10H*E10H + q01H_11H*E11H);
    
    ydot[30] = mu10H - (lambda10H + psi + (q10H_10A + q10H_10B + q10H_10C + q10H_10D + q10H_10E + q10H_10F + q10H_10G + q10H_00H + q10H_01H + q10H_11H) + mu10H) * E10H + lambda10H*E10H*E10H + (q10H_10A*E10A + q10H_10B*E10B + q10H_10C*E10C + q10H_10D*E10D + q10H_10E*E10E + q10H_10F*E10F + q10H_10G*E10G + q10H_00H*E00H + q10H_01H*E01H + q10H_11H*E11H);
    
    ydot[31] = mu11H - (lambda11H + psi + (q11H_11A + q11H_11B + q11H_11C + q11H_11D + q11H_11E + q11H_11F + q11H_11G + q11H_00H + q11H_01H + q11H_10H) + mu11H) * E11H + lambda11H*E11H*E11H + (q11H_11A*E11A + q11H_11B*E11B + q11H_11C*E11C + q11H_11D*E11D + q11H_11E*E11E + q11H_11F*E11F + q11H_11G*E11G + q11H_00H*E00H + q11H_01H*E01H + q11H_10H*E10H);
    
    
    /* The D's */
    ydot[32] =  - (lambda00A + psi + (q00A_01A + q00A_10A + q00A_11A + q00A_00B + q00A_00C + q00A_00D + q00A_00E + q00A_00F + q00A_00G + q00A_00H) + mu00A) * D00A + (q00A_01A*D01A + q00A_10A*D10A + q00A_11A*D11A + q00A_00B*D00B + q00A_00C*D00C + q00A_00D*D00D + q00A_00E*D00E + q00A_00F*D00F + q00A_00G*D00G + q00A_00H*D00H);
    
    ydot[33] =  - (lambda01A + psi + (q01A_00A + q01A_10A + q01A_11A + q01A_01B + q01A_01C + q01A_01D + q01A_01E + q01A_01F + q01A_01G + q01A_01H) + mu01A) * D01A + (q01A_00A*D00A + q01A_10A*D10A + q01A_11A*D11A + q01A_01B*D01B + q01A_01C*D01C + q01A_01D*D01D + q01A_01E*D01E + q01A_01F*D01F + q01A_01G*D01G + q01A_01H*D01H);
    
    ydot[34] =  - (lambda10A + psi + (q10A_00A + q10A_01A + q10A_11A + q10A_10B + q10A_10C + q10A_10D + q10A_10E + q10A_10F + q10A_10G + q10A_10H) + mu10A) * D10A + (q10A_00A*D00A + q10A_01A*D01A + q10A_11A*D11A + q10A_10B*D10B + q10A_10C*D10C + q10A_10D*D10D + q10A_10E*D10E + q10A_10F*D10F + q10A_10G*D10G + q10A_10H*D10H);
    
    ydot[35] =  - (lambda11A + psi + (q11A_00A + q11A_01A + q11A_10A + q11A_11B + q11A_11C + q11A_11D + q11A_11E + q11A_11F + q11A_11G + q11A_11H) + mu11A) * D11A + (q11A_00A*D00A + q11A_01A*D01A + q11A_10A*D10A + q11A_11B*D11B + q11A_11C*D11C + q11A_11D*D11D + q11A_11E*D11E + q11A_11F*D11F + q11A_11G*D11G + q11A_11H*D11H);
    
    ydot[36] =  - (lambda00B + psi + (q00B_00A + q00B_01B + q00B_10B + q00B_11B + q00B_00C + q00B_00D + q00B_00E + q00B_00F + q00B_00G + q00B_00H) + mu00B) * D00B + (q00B_00A*D00A + q00B_01B*D01B + q00B_10B*D10B + q00B_11B*D11B + q00B_00C*D00C + q00B_00D*D00D + q00B_00E*D00E + q00B_00F*D00F + q00B_00G*D00G + q00B_00H*D00H);
    
    ydot[37] =  - (lambda01B + psi + (q01B_01A + q01B_00B + q01B_10B + q01B_11B + q01B_01C + q01B_01D + q01B_01E + q01B_01F + q01B_01G + q01B_01H) + mu01B) * D01B + (q01B_01A*D01A + q01B_00B*D00B + q01B_10B*D10B + q01B_11B*D11B + q01B_01C*D01C + q01B_01D*D01D + q01B_01E*D01E + q01B_01F*D01F + q01B_01G*D01G + q01B_01H*D01H);
    
    ydot[38] =  - (lambda10B + psi + (q10B_10A + q10B_00B + q10B_01B + q10B_11B + q10B_10C + q10B_10D + q10B_10E + q10B_10F + q10B_10G + q10B_10H) + mu10B) * D10B + (q10B_10A*D10A + q10B_00B*D00B + q10B_01B*D01B + q10B_11B*D11B + q10B_10C*D10C + q10B_10D*D10D + q10B_10E*D10E + q10B_10F*D10F + q10B_10G*D10G + q10B_10H*D10H);
    
    ydot[39] =  - (lambda11B + psi + (q11B_11A + q11B_00B + q11B_01B + q11B_10B + q11B_11C + q11B_11D + q11B_11E + q11B_11F + q11B_11G + q11B_11H) + mu11B) * D11B + (q11B_11A*D11A + q11B_00B*D00B + q11B_01B*D01B + q11B_10B*D10B + q11B_11C*D11C + q11B_11D*D11D + q11B_11E*D11E + q11B_11F*D11F + q11B_11G*D11G + q11B_11H*D11H);
    
    ydot[40] =  - (lambda00C + psi + (q00C_00A + q00C_00B + q00C_01C + q00C_10C + q00C_11C + q00C_00D + q00C_00E + q00C_00F + q00C_00G + q00C_00H) + mu00C) * D00C + (q00C_00A*D00A + q00C_00B*D00B + q00C_01C*D01C + q00C_10C*D10C + q00C_11C*D11C + q00C_00D*D00D + q00C_00E*D00E + q00C_00F*D00F + q00C_00G*D00G + q00C_00H*D00H);
    
    ydot[41] =  - (lambda01C + psi + (q01C_01A + q01C_01B + q01C_00C + q01C_10C + q01C_11C + q01C_01D + q01C_01E + q01C_01F + q01C_01G + q01C_01H) + mu01C) * D01C + (q01C_01A*D01A + q01C_01B*D01B + q01C_00C*D00C + q01C_10C*D10C + q01C_11C*D11C + q01C_01D*D01D + q01C_01E*D01E + q01C_01F*D01F + q01C_01G*D01G + q01C_01H*D01H);
    
    ydot[42] =  - (lambda10C + psi + (q10C_10A + q10C_10B + q10C_00C + q10C_01C + q10C_11C + q10C_10D + q10C_10E + q10C_10F + q10C_10G + q10C_10H) + mu10C) * D10C + (q10C_10A*D10A + q10C_10B*D10B + q10C_00C*D00C + q10C_01C*D01C + q10C_11C*D11C + q10C_10D*D10D + q10C_10E*D10E + q10C_10F*D10F + q10C_10G*D10G + q10C_10H*D10H);
    
    ydot[43] =  - (lambda11C + psi + (q11C_11A + q11C_11B + q11C_00C + q11C_01C + q11C_10C + q11C_11D + q11C_11E + q11C_11F + q11C_11G + q11C_11H) + mu11C) * D11C + (q11C_11A*D11A + q11C_11B*D11B + q11C_00C*D00C + q11C_01C*D01C + q11C_10C*D10C + q11C_11D*D11D + q11C_11E*D11E + q11C_11F*D11F + q11C_11G*D11G + q11C_11H*D11H);
    
    ydot[44] =  - (lambda00D + psi + (q00D_00A + q00D_00B + q00D_00C + q00D_01D + q00D_10D + q00D_11D + q00D_00E + q00D_00F + q00D_00G + q00D_00H) + mu00D) * D00D + (q00D_00A*D00A + q00D_00B*D00B + q00D_00C*D00C + q00D_01D*D01D + q00D_10D*D10D + q00D_11D*D11D + q00D_00E*D00E + q00D_00F*D00F + q00D_00G*D00G + q00D_00H*D00H);
    
    ydot[45] =  - (lambda01D + psi + (q01D_01A + q01D_01B + q01D_01C + q01D_00D + q01D_10D + q01D_11D + q01D_01E + q01D_01F + q01D_01G + q01D_01H) + mu01D) * D01D + (q01D_01A*D01A + q01D_01B*D01B + q01D_01C*D01C + q01D_00D*D00D + q01D_10D*D10D + q01D_11D*D11D + q01D_01E*D01E + q01D_01F*D01F + q01D_01G*D01G + q01D_01H*D01H);
    
    ydot[46] =  - (lambda10D + psi + (q10D_10A + q10D_10B + q10D_10C + q10D_00D + q10D_01D + q10D_11D + q10D_10E + q10D_10F + q10D_10G + q10D_10H) + mu10D) * D10D + (q10D_10A*D10A + q10D_10B*D10B + q10D_10C*D10C + q10D_00D*D00D + q10D_01D*D01D + q10D_11D*D11D + q10D_10E*D10E + q10D_10F*D10F + q10D_10G*D10G + q10D_10H*D10H);
    
    ydot[47] =  - (lambda11D + psi + (q11D_11A + q11D_11B + q11D_11C + q11D_00D + q11D_01D + q11D_10D + q11D_11E + q11D_11F + q11D_11G + q11D_11H) + mu11D) * D11D + (q11D_11A*D11A + q11D_11B*D11B + q11D_11C*D11C + q11D_00D*D00D + q11D_01D*D01D + q11D_10D*D10D + q11D_11E*D11E + q11D_11F*D11F + q11D_11G*D11G + q11D_11H*D11H);
    
    ydot[48] =  - (lambda00E + psi + (q00E_00A + q00E_00B + q00E_00C + q00E_00D + q00E_01E + q00E_10E + q00E_11E + q00E_00F + q00E_00G + q00E_00H) + mu00E) * D00E + (q00E_00A*D00A + q00E_00B*D00B + q00E_00C*D00C + q00E_00D*D00D + q00E_01E*D01E + q00E_10E*D10E + q00E_11E*D11E + q00E_00F*D00F + q00E_00G*D00G + q00E_00H*D00H);
    
    ydot[49] =  - (lambda01E + psi + (q01E_01A + q01E_01B + q01E_01C + q01E_01D + q01E_00E + q01E_10E + q01E_11E + q01E_01F + q01E_01G + q01E_01H) + mu01E) * D01E + (q01E_01A*D01A + q01E_01B*D01B + q01E_01C*D01C + q01E_01D*D01D + q01E_00E*D00E + q01E_10E*D10E + q01E_11E*D11E + q01E_01F*D01F + q01E_01G*D01G + q01E_01H*D01H);
    
    ydot[50] =  - (lambda10E + psi + (q10E_10A + q10E_10B + q10E_10C + q10E_10D + q10E_00E + q10E_01E + q10E_11E + q10E_10F + q10E_10G + q10E_10H) + mu10E) * D10E + (q10E_10A*D10A + q10E_10B*D10B + q10E_10C*D10C + q10E_10D*D10D + q10E_00E*D00E + q10E_01E*D01E + q10E_11E*D11E + q10E_10F*D10F + q10E_10G*D10G + q10E_10H*D10H);
    
    ydot[51] =  - (lambda11E + psi + (q11E_11A + q11E_11B + q11E_11C + q11E_11D + q11E_00E + q11E_01E + q11E_10E + q11E_11F + q11E_11G + q11E_11H) + mu11E) * D11E + (q11E_11A*D11A + q11E_11B*D11B + q11E_11C*D11C + q11E_11D*D11D + q11E_00E*D00E + q11E_01E*D01E + q11E_10E*D10E + q11E_11F*D11F + q11E_11G*D11G + q11E_11H*D11H);
    
    ydot[52] =  - (lambda00F + psi + (q00F_00A + q00F_00B + q00F_00C + q00F_00D + q00F_00E + q00F_01F + q00F_10F + q00F_11F + q00F_00G + q00F_00H) + mu00F) * D00F + (q00F_00A*D00A + q00F_00B*D00B + q00F_00C*D00C + q00F_00D*D00D + q00F_00E*D00E + q00F_01F*D01F + q00F_10F*D10F + q00F_11F*D11F + q00F_00G*D00G + q00F_00H*D00H);
    
    ydot[53] =  - (lambda01F + psi + (q01F_01A + q01F_01B + q01F_01C + q01F_01D + q01F_01E + q01F_00F + q01F_10F + q01F_11F + q01F_01G + q01F_01H) + mu01F) * D01F + (q01F_01A*D01A + q01F_01B*D01B + q01F_01C*D01C + q01F_01D*D01D + q01F_01E*D01E + q01F_00F*D00F + q01F_10F*D10F + q01F_11F*D11F + q01F_01G*D01G + q01F_01H*D01H);
    
    ydot[54] =  - (lambda10F + psi + (q10F_10A + q10F_10B + q10F_10C + q10F_10D + q10F_10E + q10F_00F + q10F_01F + q10F_11F + q10F_10G + q10F_10H) + mu10F) * D10F + (q10F_10A*D10A + q10F_10B*D10B + q10F_10C*D10C + q10F_10D*D10D + q10F_10E*D10E + q10F_00F*D00F + q10F_01F*D01F + q10F_11F*D11F + q10F_10G*D10G + q10F_10H*D10H);
    
    ydot[55] =  - (lambda11F + psi + (q11F_11A + q11F_11B + q11F_11C + q11F_11D + q11F_11E + q11F_00F + q11F_01F + q11F_10F + q11F_11G + q11F_11H) + mu11F) * D11F + (q11F_11A*D11A + q11F_11B*D11B + q11F_11C*D11C + q11F_11D*D11D + q11F_11E*D11E + q11F_00F*D00F + q11F_01F*D01F + q11F_10F*D10F + q11F_11G*D11G + q11F_11H*D11H);
    
    ydot[56] =  - (lambda00G + psi + (q00G_00A + q00G_00B + q00G_00C + q00G_00D + q00G_00E + q00G_00F + q00G_01G + q00G_10G + q00G_11G + q00G_00H) + mu00G) * D00G + (q00G_00A*D00A + q00G_00B*D00B + q00G_00C*D00C + q00G_00D*D00D + q00G_00E*D00E + q00G_00F*D00F + q00G_01G*D01G + q00G_10G*D10G + q00G_11G*D11G + q00G_00H*D00H);
    
    ydot[57] =  - (lambda01G + psi + (q01G_01A + q01G_01B + q01G_01C + q01G_01D + q01G_01E + q01G_01F + q01G_00G + q01G_10G + q01G_11G + q01G_01H) + mu01G) * D01G + (q01G_01A*D01A + q01G_01B*D01B + q01G_01C*D01C + q01G_01D*D01D + q01G_01E*D01E + q01G_01F*D01F + q01G_00G*D00G + q01G_10G*D10G + q01G_11G*D11G + q01G_01H*D01H);
    
    ydot[58] =  - (lambda10G + psi + (q10G_10A + q10G_10B + q10G_10C + q10G_10D + q10G_10E + q10G_10F + q10G_00G + q10G_01G + q10G_11G + q10G_10H) + mu10G) * D10G + (q10G_10A*D10A + q10G_10B*D10B + q10G_10C*D10C + q10G_10D*D10D + q10G_10E*D10E + q10G_10F*D10F + q10G_00G*D00G + q10G_01G*D01G + q10G_11G*D11G + q10G_10H*D10H);
    
    ydot[59] =  - (lambda11G + psi + (q11G_11A + q11G_11B + q11G_11C + q11G_11D + q11G_11E + q11G_11F + q11G_00G + q11G_01G + q11G_10G + q11G_11H) + mu11G) * D11G + (q11G_11A*D11A + q11G_11B*D11B + q11G_11C*D11C + q11G_11D*D11D + q11G_11E*D11E + q11G_11F*D11F + q11G_00G*D00G + q11G_01G*D01G + q11G_10G*D10G + q11G_11H*D11H);
    
    ydot[60] =  - (lambda00H + psi + (q00H_00A + q00H_00B + q00H_00C + q00H_00D + q00H_00E + q00H_00F + q00H_00G + q00H_01H + q00H_10H + q00H_11H) + mu00H) * D00H + (q00H_00A*D00A + q00H_00B*D00B + q00H_00C*D00C + q00H_00D*D00D + q00H_00E*D00E + q00H_00F*D00F + q00H_00G*D00G + q00H_01H*D01H + q00H_10H*D10H + q00H_11H*D11H);
    
    ydot[61] =  - (lambda01H + psi + (q01H_01A + q01H_01B + q01H_01C + q01H_01D + q01H_01E + q01H_01F + q01H_01G + q01H_00H + q01H_10H + q01H_11H) + mu01H) * D01H + (q01H_01A*D01A + q01H_01B*D01B + q01H_01C*D01C + q01H_01D*D01D + q01H_01E*D01E + q01H_01F*D01F + q01H_01G*D01G + q01H_00H*D00H + q01H_10H*D10H + q01H_11H*D11H);
    
    ydot[62] =  - (lambda10H + psi + (q10H_10A + q10H_10B + q10H_10C + q10H_10D + q10H_10E + q10H_10F + q10H_10G + q10H_00H + q10H_01H + q10H_11H) + mu10H) * D10H + (q10H_10A*D10A + q10H_10B*D10B + q10H_10C*D10C + q10H_10D*D10D + q10H_10E*D10E + q10H_10F*D10F + q10H_10G*D10G + q10H_00H*D00H + q10H_01H*D01H + q10H_11H*D11H);
    
    ydot[63] =  - (lambda11H + psi + (q11H_11A + q11H_11B + q11H_11C + q11H_11D + q11H_11E + q11H_11F + q11H_11G + q11H_00H + q11H_01H + q11H_10H) + mu11H) * D11H + (q11H_11A*D11A + q11H_11B*D11B + q11H_11C*D11C + q11H_11D*D11D + q11H_11E*D11E + q11H_11F*D11F + q11H_11G*D11G + q11H_00H*D00H + q11H_01H*D01H + q11H_10H*D10H);
    
}


