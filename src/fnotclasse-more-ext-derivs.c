
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <stdio.h>
#define NUMELEMENTS 380

static double params_fhinoclass[NUMELEMENTS];


void initmod_fhinoclass(void (* odeparms)(int *, double *)){
    int N = NUMELEMENTS;
    odeparms(&N, params_fhinoclass);
}


void fnotclasse_more_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
    
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
    s00A = params_fhinoclass[0],
    s11A = params_fhinoclass[1],
    s01A = params_fhinoclass[2],
    x00A = params_fhinoclass[3],
    x11A = params_fhinoclass[4],
    x01A = 0,
    d00A_11A = params_fhinoclass[5],
    d00A_01A = params_fhinoclass[6],
    d11A_00A = params_fhinoclass[7],
    d11A_01A = params_fhinoclass[8],
    d01A_00A = params_fhinoclass[9],
    d01A_11A = params_fhinoclass[10],
    
    d00A_00B = params_fhinoclass[11],
    d00A_00C = params_fhinoclass[12],
    d00A_00D = params_fhinoclass[13],
    d00A_00E = params_fhinoclass[14],
    d00A_00F = params_fhinoclass[15],
    d00A_00G = params_fhinoclass[16],
    d00A_00H = params_fhinoclass[17],
    d00A_00I = params_fhinoclass[18],
    d00A_00J = params_fhinoclass[19],
    d11A_11B = params_fhinoclass[20],
    d11A_11C = params_fhinoclass[21],
    d11A_11D = params_fhinoclass[22],
    d11A_11E = params_fhinoclass[23],
    d11A_11F = params_fhinoclass[24],
    d11A_11G = params_fhinoclass[25],
    d11A_11H = params_fhinoclass[26],
    d11A_11I = params_fhinoclass[27],
    d11A_11J = params_fhinoclass[28],
    d01A_01B = params_fhinoclass[29],
    d01A_01C = params_fhinoclass[30],
    d01A_01D = params_fhinoclass[31],
    d01A_01E = params_fhinoclass[32],
    d01A_01F = params_fhinoclass[33],
    d01A_01G = params_fhinoclass[34],
    d01A_01H = params_fhinoclass[35],
    d01A_01I = params_fhinoclass[36],
    d01A_01J = params_fhinoclass[37],
    
    s00B = params_fhinoclass[38],
    s11B = params_fhinoclass[39],
    s01B = params_fhinoclass[40],
    x00B = params_fhinoclass[41],
    x11B = params_fhinoclass[42],
    x01B = 0,
    d00B_11B = params_fhinoclass[43],
    d00B_01B = params_fhinoclass[44],
    d11B_00B = params_fhinoclass[45],
    d11B_01B = params_fhinoclass[46],
    d01B_00B = params_fhinoclass[47],
    d01B_11B = params_fhinoclass[48],
    
    d00B_00A = params_fhinoclass[49],
    d00B_00C = params_fhinoclass[50],
    d00B_00D = params_fhinoclass[51],
    d00B_00E = params_fhinoclass[52],
    d00B_00F = params_fhinoclass[53],
    d00B_00G = params_fhinoclass[54],
    d00B_00H = params_fhinoclass[55],
    d00B_00I = params_fhinoclass[56],
    d00B_00J = params_fhinoclass[57],
    d11B_11A = params_fhinoclass[58],
    d11B_11C = params_fhinoclass[59],
    d11B_11D = params_fhinoclass[60],
    d11B_11E = params_fhinoclass[61],
    d11B_11F = params_fhinoclass[62],
    d11B_11G = params_fhinoclass[63],
    d11B_11H = params_fhinoclass[64],
    d11B_11I = params_fhinoclass[65],
    d11B_11J = params_fhinoclass[66],
    d01B_01A = params_fhinoclass[67],
    d01B_01C = params_fhinoclass[68],
    d01B_01D = params_fhinoclass[69],
    d01B_01E = params_fhinoclass[70],
    d01B_01F = params_fhinoclass[71],
    d01B_01G = params_fhinoclass[72],
    d01B_01H = params_fhinoclass[73],
    d01B_01I = params_fhinoclass[74],
    d01B_01J = params_fhinoclass[75],
    
    s00C = params_fhinoclass[76],
    s11C = params_fhinoclass[77],
    s01C = params_fhinoclass[78],
    x00C = params_fhinoclass[79],
    x11C = params_fhinoclass[80],
    x01C = 0,
    d00C_11C = params_fhinoclass[81],
    d00C_01C = params_fhinoclass[82],
    d11C_00C = params_fhinoclass[83],
    d11C_01C = params_fhinoclass[84],
    d01C_00C = params_fhinoclass[85],
    d01C_11C = params_fhinoclass[86],
    
    d00C_00A = params_fhinoclass[87],
    d00C_00B = params_fhinoclass[88],
    d00C_00D = params_fhinoclass[89],
    d00C_00E = params_fhinoclass[90],
    d00C_00F = params_fhinoclass[91],
    d00C_00G = params_fhinoclass[92],
    d00C_00H = params_fhinoclass[93],
    d00C_00I = params_fhinoclass[94],
    d00C_00J = params_fhinoclass[95],
    d11C_11A = params_fhinoclass[96],
    d11C_11B = params_fhinoclass[97],
    d11C_11D = params_fhinoclass[98],
    d11C_11E = params_fhinoclass[99],
    d11C_11F = params_fhinoclass[100],
    d11C_11G = params_fhinoclass[101],
    d11C_11H = params_fhinoclass[102],
    d11C_11I = params_fhinoclass[103],
    d11C_11J = params_fhinoclass[104],
    d01C_01A = params_fhinoclass[105],
    d01C_01B = params_fhinoclass[106],
    d01C_01D = params_fhinoclass[107],
    d01C_01E = params_fhinoclass[108],
    d01C_01F = params_fhinoclass[109],
    d01C_01G = params_fhinoclass[110],
    d01C_01H = params_fhinoclass[111],
    d01C_01I = params_fhinoclass[112],
    d01C_01J = params_fhinoclass[113],
    
    s00D = params_fhinoclass[114],
    s11D = params_fhinoclass[115],
    s01D = params_fhinoclass[116],
    x00D = params_fhinoclass[117],
    x11D = params_fhinoclass[118],
    x01D = 0,
    d00D_11D = params_fhinoclass[119],
    d00D_01D = params_fhinoclass[120],
    d11D_00D = params_fhinoclass[121],
    d11D_01D = params_fhinoclass[122],
    d01D_00D = params_fhinoclass[123],
    d01D_11D = params_fhinoclass[124],
    
    d00D_00A = params_fhinoclass[125],
    d00D_00B = params_fhinoclass[126],
    d00D_00C = params_fhinoclass[127],
    d00D_00E = params_fhinoclass[128],
    d00D_00F = params_fhinoclass[129],
    d00D_00G = params_fhinoclass[130],
    d00D_00H = params_fhinoclass[131],
    d00D_00I = params_fhinoclass[132],
    d00D_00J = params_fhinoclass[133],
    d11D_11A = params_fhinoclass[134],
    d11D_11B = params_fhinoclass[135],
    d11D_11C = params_fhinoclass[136],
    d11D_11E = params_fhinoclass[137],
    d11D_11F = params_fhinoclass[138],
    d11D_11G = params_fhinoclass[139],
    d11D_11H = params_fhinoclass[140],
    d11D_11I = params_fhinoclass[141],
    d11D_11J = params_fhinoclass[142],
    d01D_01A = params_fhinoclass[143],
    d01D_01B = params_fhinoclass[144],
    d01D_01C = params_fhinoclass[145],
    d01D_01E = params_fhinoclass[146],
    d01D_01F = params_fhinoclass[147],
    d01D_01G = params_fhinoclass[148],
    d01D_01H = params_fhinoclass[149],
    d01D_01I = params_fhinoclass[150],
    d01D_01J = params_fhinoclass[151],
    
    s00E = params_fhinoclass[152],
    s11E = params_fhinoclass[153],
    s01E = params_fhinoclass[154],
    x00E = params_fhinoclass[155],
    x11E = params_fhinoclass[156],
    x01E = 0,
    d00E_11E = params_fhinoclass[157],
    d00E_01E = params_fhinoclass[158],
    d11E_00E = params_fhinoclass[159],
    d11E_01E = params_fhinoclass[160],
    d01E_00E = params_fhinoclass[161],
    d01E_11E = params_fhinoclass[162],
    
    d00E_00A = params_fhinoclass[163],
    d00E_00B = params_fhinoclass[164],
    d00E_00C = params_fhinoclass[165],
    d00E_00D = params_fhinoclass[166],
    d00E_00F = params_fhinoclass[167],
    d00E_00G = params_fhinoclass[168],
    d00E_00H = params_fhinoclass[169],
    d00E_00I = params_fhinoclass[170],
    d00E_00J = params_fhinoclass[171],
    d11E_11A = params_fhinoclass[172],
    d11E_11B = params_fhinoclass[173],
    d11E_11C = params_fhinoclass[174],
    d11E_11D = params_fhinoclass[175],
    d11E_11F = params_fhinoclass[176],
    d11E_11G = params_fhinoclass[177],
    d11E_11H = params_fhinoclass[178],
    d11E_11I = params_fhinoclass[179],
    d11E_11J = params_fhinoclass[180],
    d01E_01A = params_fhinoclass[181],
    d01E_01B = params_fhinoclass[182],
    d01E_01C = params_fhinoclass[183],
    d01E_01D = params_fhinoclass[184],
    d01E_01F = params_fhinoclass[185],
    d01E_01G = params_fhinoclass[186],
    d01E_01H = params_fhinoclass[187],
    d01E_01I = params_fhinoclass[188],
    d01E_01J = params_fhinoclass[189],
    
    s00F = params_fhinoclass[190],
    s11F = params_fhinoclass[191],
    s01F = params_fhinoclass[192],
    x00F = params_fhinoclass[193],
    x11F = params_fhinoclass[194],
    x01F = 0,
    d00F_11F = params_fhinoclass[195],
    d00F_01F = params_fhinoclass[196],
    d11F_00F = params_fhinoclass[197],
    d11F_01F = params_fhinoclass[198],
    d01F_00F = params_fhinoclass[199],
    d01F_11F = params_fhinoclass[200],
    
    d00F_00A = params_fhinoclass[201],
    d00F_00B = params_fhinoclass[202],
    d00F_00C = params_fhinoclass[203],
    d00F_00D = params_fhinoclass[204],
    d00F_00E = params_fhinoclass[205],
    d00F_00G = params_fhinoclass[206],
    d00F_00H = params_fhinoclass[207],
    d00F_00I = params_fhinoclass[208],
    d00F_00J = params_fhinoclass[209],
    d11F_11A = params_fhinoclass[210],
    d11F_11B = params_fhinoclass[211],
    d11F_11C = params_fhinoclass[212],
    d11F_11D = params_fhinoclass[213],
    d11F_11E = params_fhinoclass[214],
    d11F_11G = params_fhinoclass[215],
    d11F_11H = params_fhinoclass[216],
    d11F_11I = params_fhinoclass[217],
    d11F_11J = params_fhinoclass[218],
    d01F_01A = params_fhinoclass[219],
    d01F_01B = params_fhinoclass[220],
    d01F_01C = params_fhinoclass[221],
    d01F_01D = params_fhinoclass[222],
    d01F_01E = params_fhinoclass[223],
    d01F_01G = params_fhinoclass[224],
    d01F_01H = params_fhinoclass[225],
    d01F_01I = params_fhinoclass[226],
    d01F_01J = params_fhinoclass[227],
    
    s00G = params_fhinoclass[228],
    s11G = params_fhinoclass[229],
    s01G = params_fhinoclass[230],
    x00G = params_fhinoclass[231],
    x11G = params_fhinoclass[232],
    x01G = 0,
    d00G_11G = params_fhinoclass[233],
    d00G_01G = params_fhinoclass[234],
    d11G_00G = params_fhinoclass[235],
    d11G_01G = params_fhinoclass[236],
    d01G_00G = params_fhinoclass[237],
    d01G_11G = params_fhinoclass[238],
    
    d00G_00A = params_fhinoclass[239],
    d00G_00B = params_fhinoclass[240],
    d00G_00C = params_fhinoclass[241],
    d00G_00D = params_fhinoclass[242],
    d00G_00E = params_fhinoclass[243],
    d00G_00F = params_fhinoclass[244],
    d00G_00H = params_fhinoclass[245],
    d00G_00I = params_fhinoclass[246],
    d00G_00J = params_fhinoclass[247],
    d11G_11A = params_fhinoclass[248],
    d11G_11B = params_fhinoclass[249],
    d11G_11C = params_fhinoclass[250],
    d11G_11D = params_fhinoclass[251],
    d11G_11E = params_fhinoclass[252],
    d11G_11F = params_fhinoclass[253],
    d11G_11H = params_fhinoclass[254],
    d11G_11I = params_fhinoclass[255],
    d11G_11J = params_fhinoclass[256],
    d01G_01A = params_fhinoclass[257],
    d01G_01B = params_fhinoclass[258],
    d01G_01C = params_fhinoclass[259],
    d01G_01D = params_fhinoclass[260],
    d01G_01E = params_fhinoclass[261],
    d01G_01F = params_fhinoclass[262],
    d01G_01H = params_fhinoclass[263],
    d01G_01I = params_fhinoclass[264],
    d01G_01J = params_fhinoclass[265],
    
    s00H = params_fhinoclass[266],
    s11H = params_fhinoclass[267],
    s01H = params_fhinoclass[268],
    x00H = params_fhinoclass[269],
    x11H = params_fhinoclass[270],
    x01H = 0,
    d00H_11H = params_fhinoclass[271],
    d00H_01H = params_fhinoclass[272],
    d11H_00H = params_fhinoclass[273],
    d11H_01H = params_fhinoclass[274],
    d01H_00H = params_fhinoclass[275],
    d01H_11H = params_fhinoclass[276],
    
    d00H_00A = params_fhinoclass[277],
    d00H_00B = params_fhinoclass[278],
    d00H_00C = params_fhinoclass[279],
    d00H_00D = params_fhinoclass[280],
    d00H_00E = params_fhinoclass[281],
    d00H_00F = params_fhinoclass[282],
    d00H_00G = params_fhinoclass[283],
    d00H_00I = params_fhinoclass[284],
    d00H_00J = params_fhinoclass[285],
    d11H_11A = params_fhinoclass[286],
    d11H_11B = params_fhinoclass[287],
    d11H_11C = params_fhinoclass[288],
    d11H_11D = params_fhinoclass[289],
    d11H_11E = params_fhinoclass[290],
    d11H_11F = params_fhinoclass[291],
    d11H_11G = params_fhinoclass[292],
    d11H_11I = params_fhinoclass[293],
    d11H_11J = params_fhinoclass[294],
    d01H_01A = params_fhinoclass[295],
    d01H_01B = params_fhinoclass[296],
    d01H_01C = params_fhinoclass[297],
    d01H_01D = params_fhinoclass[298],
    d01H_01E = params_fhinoclass[299],
    d01H_01F = params_fhinoclass[300],
    d01H_01G = params_fhinoclass[301],
    d01H_01I = params_fhinoclass[302],
    d01H_01J = params_fhinoclass[303],
    
    s00I = params_fhinoclass[304],
    s11I = params_fhinoclass[305],
    s01I = params_fhinoclass[306],
    x00I = params_fhinoclass[307],
    x11I = params_fhinoclass[308],
    x01I = 0,
    d00I_11I = params_fhinoclass[309],
    d00I_01I = params_fhinoclass[310],
    d11I_00I = params_fhinoclass[311],
    d11I_01I = params_fhinoclass[312],
    d01I_00I = params_fhinoclass[313],
    d01I_11I = params_fhinoclass[314],
    
    d00I_00A = params_fhinoclass[315],
    d00I_00B = params_fhinoclass[316],
    d00I_00C = params_fhinoclass[317],
    d00I_00D = params_fhinoclass[318],
    d00I_00E = params_fhinoclass[319],
    d00I_00F = params_fhinoclass[320],
    d00I_00G = params_fhinoclass[321],
    d00I_00H = params_fhinoclass[322],
    d00I_00J = params_fhinoclass[323],
    d11I_11A = params_fhinoclass[324],
    d11I_11B = params_fhinoclass[325],
    d11I_11C = params_fhinoclass[326],
    d11I_11D = params_fhinoclass[327],
    d11I_11E = params_fhinoclass[328],
    d11I_11F = params_fhinoclass[329],
    d11I_11G = params_fhinoclass[330],
    d11I_11H = params_fhinoclass[331],
    d11I_11J = params_fhinoclass[332],
    d01I_01A = params_fhinoclass[333],
    d01I_01B = params_fhinoclass[334],
    d01I_01C = params_fhinoclass[335],
    d01I_01D = params_fhinoclass[336],
    d01I_01E = params_fhinoclass[337],
    d01I_01F = params_fhinoclass[338],
    d01I_01G = params_fhinoclass[339],
    d01I_01H = params_fhinoclass[340],
    d01I_01J = params_fhinoclass[341],
    
    s00J = params_fhinoclass[342],
    s11J = params_fhinoclass[343],
    s01J = params_fhinoclass[344],
    x00J = params_fhinoclass[345],
    x11J = params_fhinoclass[346],
    x01J = 0,
    d00J_11J = params_fhinoclass[347],
    d00J_01J = params_fhinoclass[348],
    d11J_00J = params_fhinoclass[349],
    d11J_01J = params_fhinoclass[350],
    d01J_00J = params_fhinoclass[351],
    d01J_11J = params_fhinoclass[352],
    
    d00J_00A = params_fhinoclass[353],
    d00J_00B = params_fhinoclass[354],
    d00J_00C = params_fhinoclass[355],
    d00J_00D = params_fhinoclass[356],
    d00J_00E = params_fhinoclass[357],
    d00J_00F = params_fhinoclass[358],
    d00J_00G = params_fhinoclass[359],
    d00J_00H = params_fhinoclass[360],
    d00J_00I = params_fhinoclass[361],
    d11J_11A = params_fhinoclass[362],
    d11J_11B = params_fhinoclass[363],
    d11J_11C = params_fhinoclass[364],
    d11J_11D = params_fhinoclass[365],
    d11J_11E = params_fhinoclass[366],
    d11J_11F = params_fhinoclass[367],
    d11J_11G = params_fhinoclass[368],
    d11J_11H = params_fhinoclass[369],
    d11J_11I = params_fhinoclass[370],
    d01J_01A = params_fhinoclass[371],
    d01J_01B = params_fhinoclass[372],
    d01J_01C = params_fhinoclass[373],
    d01J_01D = params_fhinoclass[374],
    d01J_01E = params_fhinoclass[375],
    d01J_01F = params_fhinoclass[376],
    d01J_01G = params_fhinoclass[377],
    d01J_01H = params_fhinoclass[378],
    d01J_01I = params_fhinoclass[379];
    


/* The E's */


ydot[0] = x00A - (s00A + (d00A_11A + d00A_01A + d00A_00B + d00A_00C + d00A_00D + d00A_00E + d00A_00F + d00A_00G + d00A_00H + d00A_00I + d00A_00J) + x00A) * E00A + s00A*E00A*E00A + (d00A_11A*E11A + d00A_01A*E01A + d00A_00B*E00B + d00A_00C*E00C + d00A_00D*E00D + d00A_00E*E00E + d00A_00F*E00F + d00A_00G*E00G + d00A_00H*E00H + d00A_00I*E00I + d00A_00J*E00J);

ydot[1] = x11A - (s11A + (d11A_00A + d11A_01A + d11A_11B + d11A_11C + d11A_11D + d11A_11E + d11A_11F + d11A_11G + d11A_11H + d11A_11I + d11A_11J) + x11A) * E11A + s11A*E11A*E11A + (d11A_00A*E00A + d11A_01A*E01A + d11A_11B*E11B + d11A_11C*E11C + d11A_11D*E11D + d11A_11E*E11E + d11A_11F*E11F + d11A_11G*E11G + d11A_11H*E11H + d11A_11I*E11I + d11A_11J*E11J);

ydot[2] = x01A - (s01A + (d01A_00A + d01A_11A + d01A_01B + d01A_01C + d01A_01D + d01A_01E + d01A_01F + d01A_01G + d01A_01H + d01A_01I + d01A_01J) + x01A) * E01A + s01A*E01A*E01A + (d01A_00A*E00A + d01A_11A*E11A + d01A_01B*E01B + d01A_01C*E01C + d01A_01D*E01D + d01A_01E*E01E + d01A_01F*E01F + d01A_01G*E01G + d01A_01H*E01H + d01A_01I*E01I + d01A_01J*E01J);

ydot[3] = x00B - (s00B + (d00B_00A + d00B_11B + d00B_01B + d00B_00C + d00B_00D + d00B_00E + d00B_00F + d00B_00G + d00B_00H + d00B_00I + d00B_00J) + x00B) * E00B + s00B*E00B*E00B + (d00B_00A*E00A + d00B_11B*E11B + d00B_01B*E01B + d00B_00C*E00C + d00B_00D*E00D + d00B_00E*E00E + d00B_00F*E00F + d00B_00G*E00G + d00B_00H*E00H + d00B_00I*E00I + d00B_00J*E00J);

ydot[4] = x11B - (s11B + (d11B_11A + d11B_00B + d11B_01B + d11B_11C + d11B_11D + d11B_11E + d11B_11F + d11B_11G + d11B_11H + d11B_11I + d11B_11J) + x11B) * E11B + s11B*E11B*E11B + (d11B_11A*E11A + d11B_00B*E00B + d11B_01B*E01B + d11B_11C*E11C + d11B_11D*E11D + d11B_11E*E11E + d11B_11F*E11F + d11B_11G*E11G + d11B_11H*E11H + d11B_11I*E11I + d11B_11J*E11J);

ydot[5] = x01B - (s01B + (d01B_01A + d01B_00B + d01B_11B + d01B_01C + d01B_01D + d01B_01E + d01B_01F + d01B_01G + d01B_01H + d01B_01I + d01B_01J) + x01B) * E01B + s01B*E01B*E01B + (d01B_01A*E01A + d01B_00B*E00B + d01B_11B*E11B + d01B_01C*E01C + d01B_01D*E01D + d01B_01E*E01E + d01B_01F*E01F + d01B_01G*E01G + d01B_01H*E01H + d01B_01I*E01I + d01B_01J*E01J);

ydot[6] = x00C - (s00C + (d00C_00A + d00C_00B + d00C_11C + d00C_01C + d00C_00D + d00C_00E + d00C_00F + d00C_00G + d00C_00H + d00C_00I + d00C_00J) + x00C) * E00C + s00C*E00C*E00C + (d00C_00A*E00A + d00C_00B*E00B + d00C_11C*E11C + d00C_01C*E01C + d00C_00D*E00D + d00C_00E*E00E + d00C_00F*E00F + d00C_00G*E00G + d00C_00H*E00H + d00C_00I*E00I + d00C_00J*E00J);

ydot[7] = x11C - (s11C + (d11C_11A + d11C_11B + d11C_00C + d11C_01C + d11C_11D + d11C_11E + d11C_11F + d11C_11G + d11C_11H + d11C_11I + d11C_11J) + x11C) * E11C + s11C*E11C*E11C + (d11C_11A*E11A + d11C_11B*E11B + d11C_00C*E00C + d11C_01C*E01C + d11C_11D*E11D + d11C_11E*E11E + d11C_11F*E11F + d11C_11G*E11G + d11C_11H*E11H + d11C_11I*E11I + d11C_11J*E11J);

ydot[8] = x01C - (s01C + (d01C_01A + d01C_01B + d01C_00C + d01C_11C + d01C_01D + d01C_01E + d01C_01F + d01C_01G + d01C_01H + d01C_01I + d01C_01J) + x01C) * E01C + s01C*E01C*E01C + (d01C_01A*E01A + d01C_01B*E01B + d01C_00C*E00C + d01C_11C*E11C + d01C_01D*E01D + d01C_01E*E01E + d01C_01F*E01F + d01C_01G*E01G + d01C_01H*E01H + d01C_01I*E01I + d01C_01J*E01J);

ydot[9] = x00D - (s00D + (d00D_00A + d00D_00B + d00D_00C + d00D_11D + d00D_01D + d00D_00E + d00D_00F + d00D_00G + d00D_00H + d00D_00I + d00D_00J) + x00D) * E00D + s00D*E00D*E00D + (d00D_00A*E00A + d00D_00B*E00B + d00D_00C*E00C + d00D_11D*E11D + d00D_01D*E01D + d00D_00E*E00E + d00D_00F*E00F + d00D_00G*E00G + d00D_00H*E00H + d00D_00I*E00I + d00D_00J*E00J);

ydot[10] = x11D - (s11D + (d11D_11A + d11D_11B + d11D_11C + d11D_00D + d11D_01D + d11D_11E + d11D_11F + d11D_11G + d11D_11H + d11D_11I + d11D_11J) + x11D) * E11D + s11D*E11D*E11D + (d11D_11A*E11A + d11D_11B*E11B + d11D_11C*E11C + d11D_00D*E00D + d11D_01D*E01D + d11D_11E*E11E + d11D_11F*E11F + d11D_11G*E11G + d11D_11H*E11H + d11D_11I*E11I + d11D_11J*E11J);

ydot[11] = x01D - (s01D + (d01D_01A + d01D_01B + d01D_01C + d01D_00D + d01D_11D + d01D_01E + d01D_01F + d01D_01G + d01D_01H + d01D_01I + d01D_01J) + x01D) * E01D + s01D*E01D*E01D + (d01D_01A*E01A + d01D_01B*E01B + d01D_01C*E01C + d01D_00D*E00D + d01D_11D*E11D + d01D_01E*E01E + d01D_01F*E01F + d01D_01G*E01G + d01D_01H*E01H + d01D_01I*E01I + d01D_01J*E01J);

ydot[12] = x00E - (s00E + (d00E_00A + d00E_00B + d00E_00C + d00E_00D + d00E_11E + d00E_01E + d00E_00F + d00E_00G + d00E_00H + d00E_00I + d00E_00J) + x00E) * E00E + s00E*E00E*E00E + (d00E_00A*E00A + d00E_00B*E00B + d00E_00C*E00C + d00E_00D*E00D + d00E_11E*E11E + d00E_01E*E01E + d00E_00F*E00F + d00E_00G*E00G + d00E_00H*E00H + d00E_00I*E00I + d00E_00J*E00J);

ydot[13] = x11E - (s11E + (d11E_11A + d11E_11B + d11E_11C + d11E_11D + d11E_00E + d11E_01E + d11E_11F + d11E_11G + d11E_11H + d11E_11I + d11E_11J) + x11E) * E11E + s11E*E11E*E11E + (d11E_11A*E11A + d11E_11B*E11B + d11E_11C*E11C + d11E_11D*E11D + d11E_00E*E00E + d11E_01E*E01E + d11E_11F*E11F + d11E_11G*E11G + d11E_11H*E11H + d11E_11I*E11I + d11E_11J*E11J);

ydot[14] = x01E - (s01E + (d01E_01A + d01E_01B + d01E_01C + d01E_01D + d01E_00E + d01E_11E + d01E_01F + d01E_01G + d01E_01H + d01E_01I + d01E_01J) + x01E) * E01E + s01E*E01E*E01E + (d01E_01A*E01A + d01E_01B*E01B + d01E_01C*E01C + d01E_01D*E01D + d01E_00E*E00E + d01E_11E*E11E + d01E_01F*E01F + d01E_01G*E01G + d01E_01H*E01H + d01E_01I*E01I + d01E_01J*E01J);

ydot[15] = x00F - (s00F + (d00F_00A + d00F_00B + d00F_00C + d00F_00D + d00F_00E + d00F_11F + d00F_01F + d00F_00G + d00F_00H + d00F_00I + d00F_00J) + x00F) * E00F + s00F*E00F*E00F + (d00F_00A*E00A + d00F_00B*E00B + d00F_00C*E00C + d00F_00D*E00D + d00F_00E*E00E + d00F_11F*E11F + d00F_01F*E01F + d00F_00G*E00G + d00F_00H*E00H + d00F_00I*E00I + d00F_00J*E00J);

ydot[16] = x11F - (s11F + (d11F_11A + d11F_11B + d11F_11C + d11F_11D + d11F_11E + d11F_00F + d11F_01F + d11F_11G + d11F_11H + d11F_11I + d11F_11J) + x11F) * E11F + s11F*E11F*E11F + (d11F_11A*E11A + d11F_11B*E11B + d11F_11C*E11C + d11F_11D*E11D + d11F_11E*E11E + d11F_00F*E00F + d11F_01F*E01F + d11F_11G*E11G + d11F_11H*E11H + d11F_11I*E11I + d11F_11J*E11J);

ydot[17] = x01F - (s01F + (d01F_01A + d01F_01B + d01F_01C + d01F_01D + d01F_01E + d01F_00F + d01F_11F + d01F_01G + d01F_01H + d01F_01I + d01F_01J) + x01F) * E01F + s01F*E01F*E01F + (d01F_01A*E01A + d01F_01B*E01B + d01F_01C*E01C + d01F_01D*E01D + d01F_01E*E01E + d01F_00F*E00F + d01F_11F*E11F + d01F_01G*E01G + d01F_01H*E01H + d01F_01I*E01I + d01F_01J*E01J);

ydot[18] = x00G - (s00G + (d00G_00A + d00G_00B + d00G_00C + d00G_00D + d00G_00E + d00G_00F + d00G_11G + d00G_01G + d00G_00H + d00G_00I + d00G_00J) + x00G) * E00G + s00G*E00G*E00G + (d00G_00A*E00A + d00G_00B*E00B + d00G_00C*E00C + d00G_00D*E00D + d00G_00E*E00E + d00G_00F*E00F + d00G_11G*E11G + d00G_01G*E01G + d00G_00H*E00H + d00G_00I*E00I + d00G_00J*E00J);

ydot[19] = x11G - (s11G + (d11G_11A + d11G_11B + d11G_11C + d11G_11D + d11G_11E + d11G_11F + d11G_00G + d11G_01G + d11G_11H + d11G_11I + d11G_11J) + x11G) * E11G + s11G*E11G*E11G + (d11G_11A*E11A + d11G_11B*E11B + d11G_11C*E11C + d11G_11D*E11D + d11G_11E*E11E + d11G_11F*E11F + d11G_00G*E00G + d11G_01G*E01G + d11G_11H*E11H + d11G_11I*E11I + d11G_11J*E11J);

ydot[20] = x01G - (s01G + (d01G_01A + d01G_01B + d01G_01C + d01G_01D + d01G_01E + d01G_01F + d01G_00G + d01G_11G + d01G_01H + d01G_01I + d01G_01J) + x01G) * E01G + s01G*E01G*E01G + (d01G_01A*E01A + d01G_01B*E01B + d01G_01C*E01C + d01G_01D*E01D + d01G_01E*E01E + d01G_01F*E01F + d01G_00G*E00G + d01G_11G*E11G + d01G_01H*E01H + d01G_01I*E01I + d01G_01J*E01J);

ydot[21] = x00H - (s00H + (d00H_00A + d00H_00B + d00H_00C + d00H_00D + d00H_00E + d00H_00F + d00H_00G + d00H_11H + d00H_01H + d00H_00I + d00H_00J) + x00H) * E00H + s00H*E00H*E00H + (d00H_00A*E00A + d00H_00B*E00B + d00H_00C*E00C + d00H_00D*E00D + d00H_00E*E00E + d00H_00F*E00F + d00H_00G*E00G + d00H_11H*E11H + d00H_01H*E01H + d00H_00I*E00I + d00H_00J*E00J);

ydot[22] = x11H - (s11H + (d11H_11A + d11H_11B + d11H_11C + d11H_11D + d11H_11E + d11H_11F + d11H_11G + d11H_00H + d11H_01H + d11H_11I + d11H_11J) + x11H) * E11H + s11H*E11H*E11H + (d11H_11A*E11A + d11H_11B*E11B + d11H_11C*E11C + d11H_11D*E11D + d11H_11E*E11E + d11H_11F*E11F + d11H_11G*E11G + d11H_00H*E00H + d11H_01H*E01H + d11H_11I*E11I + d11H_11J*E11J);

ydot[23] = x01H - (s01H + (d01H_01A + d01H_01B + d01H_01C + d01H_01D + d01H_01E + d01H_01F + d01H_01G + d01H_00H + d01H_11H + d01H_01I + d01H_01J) + x01H) * E01H + s01H*E01H*E01H + (d01H_01A*E01A + d01H_01B*E01B + d01H_01C*E01C + d01H_01D*E01D + d01H_01E*E01E + d01H_01F*E01F + d01H_01G*E01G + d01H_00H*E00H + d01H_11H*E11H + d01H_01I*E01I + d01H_01J*E01J);

ydot[24] = x00I - (s00I + (d00I_00A + d00I_00B + d00I_00C + d00I_00D + d00I_00E + d00I_00F + d00I_00G + d00I_00H + d00I_11I + d00I_01I + d00I_00J) + x00I) * E00I + s00I*E00I*E00I + (d00I_00A*E00A + d00I_00B*E00B + d00I_00C*E00C + d00I_00D*E00D + d00I_00E*E00E + d00I_00F*E00F + d00I_00G*E00G + d00I_00H*E00H + d00I_11I*E11I + d00I_01I*E01I + d00I_00J*E00J);

ydot[25] = x11I - (s11I + (d11I_11A + d11I_11B + d11I_11C + d11I_11D + d11I_11E + d11I_11F + d11I_11G + d11I_11H + d11I_00I + d11I_01I + d11I_11J) + x11I) * E11I + s11I*E11I*E11I + (d11I_11A*E11A + d11I_11B*E11B + d11I_11C*E11C + d11I_11D*E11D + d11I_11E*E11E + d11I_11F*E11F + d11I_11G*E11G + d11I_11H*E11H + d11I_00I*E00I + d11I_01I*E01I + d11I_11J*E11J);

ydot[26] = x01I - (s01I + (d01I_01A + d01I_01B + d01I_01C + d01I_01D + d01I_01E + d01I_01F + d01I_01G + d01I_01H + d01I_00I + d01I_11I + d01I_01J) + x01I) * E01I + s01I*E01I*E01I + (d01I_01A*E01A + d01I_01B*E01B + d01I_01C*E01C + d01I_01D*E01D + d01I_01E*E01E + d01I_01F*E01F + d01I_01G*E01G + d01I_01H*E01H + d01I_00I*E00I + d01I_11I*E11I + d01I_01J*E01J);

ydot[27] = x00J - (s00J + (d00J_00A + d00J_00B + d00J_00C + d00J_00D + d00J_00E + d00J_00F + d00J_00G + d00J_00H + d00J_00I + d00J_11J + d00J_01J) + x00J) * E00J + s00J*E00J*E00J + (d00J_00A*E00A + d00J_00B*E00B + d00J_00C*E00C + d00J_00D*E00D + d00J_00E*E00E + d00J_00F*E00F + d00J_00G*E00G + d00J_00H*E00H + d00J_00I*E00I + d00J_11J*E11J + d00J_01J*E01J);

ydot[28] = x11J - (s11J + (d11J_11A + d11J_11B + d11J_11C + d11J_11D + d11J_11E + d11J_11F + d11J_11G + d11J_11H + d11J_11I + d11J_00J + d11J_01J) + x11J) * E11J + s11J*E11J*E11J + (d11J_11A*E11A + d11J_11B*E11B + d11J_11C*E11C + d11J_11D*E11D + d11J_11E*E11E + d11J_11F*E11F + d11J_11G*E11G + d11J_11H*E11H + d11J_11I*E11I + d11J_00J*E00J + d11J_01J*E01J);

ydot[29] = x01J - (s01J + (d01J_01A + d01J_01B + d01J_01C + d01J_01D + d01J_01E + d01J_01F + d01J_01G + d01J_01H + d01J_01I + d01J_00J + d01J_11J) + x01J) * E01J + s01J*E01J*E01J + (d01J_01A*E01A + d01J_01B*E01B + d01J_01C*E01C + d01J_01D*E01D + d01J_01E*E01E + d01J_01F*E01F + d01J_01G*E01G + d01J_01H*E01H + d01J_01I*E01I + d01J_00J*E00J + d01J_11J*E11J);


/* The D's */

ydot[30] =  - (s00A + (d00A_11A + d00A_01A + d00A_00B + d00A_00C + d00A_00D + d00A_00E + d00A_00F + d00A_00G + d00A_00H + d00A_00I + d00A_00J) + x00A) * D00A + 2*s00A*E00A*D00A + (d00A_11A*D11A + d00A_01A*D01A + d00A_00B*D00B + d00A_00C*D00C + d00A_00D*D00D + d00A_00E*D00E + d00A_00F*D00F + d00A_00G*D00G + d00A_00H*D00H + d00A_00I*D00I + d00A_00J*D00J);

ydot[31] =  - (s11A + (d11A_00A + d11A_01A + d11A_11B + d11A_11C + d11A_11D + d11A_11E + d11A_11F + d11A_11G + d11A_11H + d11A_11I + d11A_11J) + x11A) * D11A + 2*s11A*E11A*D11A + (d11A_00A*D00A + d11A_01A*D01A + d11A_11B*D11B + d11A_11C*D11C + d11A_11D*D11D + d11A_11E*D11E + d11A_11F*D11F + d11A_11G*D11G + d11A_11H*D11H + d11A_11I*D11I + d11A_11J*D11J);

ydot[32] =  - (s01A + (d01A_00A + d01A_11A + d01A_01B + d01A_01C + d01A_01D + d01A_01E + d01A_01F + d01A_01G + d01A_01H + d01A_01I + d01A_01J) + x01A) * D01A + 2*s01A*E01A*D01A + (d01A_00A*D00A + d01A_11A*D11A + d01A_01B*D01B + d01A_01C*D01C + d01A_01D*D01D + d01A_01E*D01E + d01A_01F*D01F + d01A_01G*D01G + d01A_01H*D01H + d01A_01I*D01I + d01A_01J*D01J);

ydot[33] =  - (s00B + (d00B_00A + d00B_11B + d00B_01B + d00B_00C + d00B_00D + d00B_00E + d00B_00F + d00B_00G + d00B_00H + d00B_00I + d00B_00J) + x00B) * D00B + 2*s00B*E00B*D00B + (d00B_00A*D00A + d00B_11B*D11B + d00B_01B*D01B + d00B_00C*D00C + d00B_00D*D00D + d00B_00E*D00E + d00B_00F*D00F + d00B_00G*D00G + d00B_00H*D00H + d00B_00I*D00I + d00B_00J*D00J);

ydot[34] =  - (s11B + (d11B_11A + d11B_00B + d11B_01B + d11B_11C + d11B_11D + d11B_11E + d11B_11F + d11B_11G + d11B_11H + d11B_11I + d11B_11J) + x11B) * D11B + 2*s11B*E11B*D11B + (d11B_11A*D11A + d11B_00B*D00B + d11B_01B*D01B + d11B_11C*D11C + d11B_11D*D11D + d11B_11E*D11E + d11B_11F*D11F + d11B_11G*D11G + d11B_11H*D11H + d11B_11I*D11I + d11B_11J*D11J);

ydot[35] =  - (s01B + (d01B_01A + d01B_00B + d01B_11B + d01B_01C + d01B_01D + d01B_01E + d01B_01F + d01B_01G + d01B_01H + d01B_01I + d01B_01J) + x01B) * D01B + 2*s01B*E01B*D01B + (d01B_01A*D01A + d01B_00B*D00B + d01B_11B*D11B + d01B_01C*D01C + d01B_01D*D01D + d01B_01E*D01E + d01B_01F*D01F + d01B_01G*D01G + d01B_01H*D01H + d01B_01I*D01I + d01B_01J*D01J);

ydot[36] =  - (s00C + (d00C_00A + d00C_00B + d00C_11C + d00C_01C + d00C_00D + d00C_00E + d00C_00F + d00C_00G + d00C_00H + d00C_00I + d00C_00J) + x00C) * D00C + 2*s00C*E00C*D00C + (d00C_00A*D00A + d00C_00B*D00B + d00C_11C*D11C + d00C_01C*D01C + d00C_00D*D00D + d00C_00E*D00E + d00C_00F*D00F + d00C_00G*D00G + d00C_00H*D00H + d00C_00I*D00I + d00C_00J*D00J);

ydot[37] =  - (s11C + (d11C_11A + d11C_11B + d11C_00C + d11C_01C + d11C_11D + d11C_11E + d11C_11F + d11C_11G + d11C_11H + d11C_11I + d11C_11J) + x11C) * D11C + 2*s11C*E11C*D11C + (d11C_11A*D11A + d11C_11B*D11B + d11C_00C*D00C + d11C_01C*D01C + d11C_11D*D11D + d11C_11E*D11E + d11C_11F*D11F + d11C_11G*D11G + d11C_11H*D11H + d11C_11I*D11I + d11C_11J*D11J);

ydot[38] =  - (s01C + (d01C_01A + d01C_01B + d01C_00C + d01C_11C + d01C_01D + d01C_01E + d01C_01F + d01C_01G + d01C_01H + d01C_01I + d01C_01J) + x01C) * D01C + 2*s01C*E01C*D01C + (d01C_01A*D01A + d01C_01B*D01B + d01C_00C*D00C + d01C_11C*D11C + d01C_01D*D01D + d01C_01E*D01E + d01C_01F*D01F + d01C_01G*D01G + d01C_01H*D01H + d01C_01I*D01I + d01C_01J*D01J);

ydot[39] =  - (s00D + (d00D_00A + d00D_00B + d00D_00C + d00D_11D + d00D_01D + d00D_00E + d00D_00F + d00D_00G + d00D_00H + d00D_00I + d00D_00J) + x00D) * D00D + 2*s00D*E00D*D00D + (d00D_00A*D00A + d00D_00B*D00B + d00D_00C*D00C + d00D_11D*D11D + d00D_01D*D01D + d00D_00E*D00E + d00D_00F*D00F + d00D_00G*D00G + d00D_00H*D00H + d00D_00I*D00I + d00D_00J*D00J);

ydot[40] =  - (s11D + (d11D_11A + d11D_11B + d11D_11C + d11D_00D + d11D_01D + d11D_11E + d11D_11F + d11D_11G + d11D_11H + d11D_11I + d11D_11J) + x11D) * D11D + 2*s11D*E11D*D11D + (d11D_11A*D11A + d11D_11B*D11B + d11D_11C*D11C + d11D_00D*D00D + d11D_01D*D01D + d11D_11E*D11E + d11D_11F*D11F + d11D_11G*D11G + d11D_11H*D11H + d11D_11I*D11I + d11D_11J*D11J);

ydot[41] =  - (s01D + (d01D_01A + d01D_01B + d01D_01C + d01D_00D + d01D_11D + d01D_01E + d01D_01F + d01D_01G + d01D_01H + d01D_01I + d01D_01J) + x01D) * D01D + 2*s01D*E01D*D01D + (d01D_01A*D01A + d01D_01B*D01B + d01D_01C*D01C + d01D_00D*D00D + d01D_11D*D11D + d01D_01E*D01E + d01D_01F*D01F + d01D_01G*D01G + d01D_01H*D01H + d01D_01I*D01I + d01D_01J*D01J);

ydot[42] =  - (s00E + (d00E_00A + d00E_00B + d00E_00C + d00E_00D + d00E_11E + d00E_01E + d00E_00F + d00E_00G + d00E_00H + d00E_00I + d00E_00J) + x00E) * D00E + 2*s00E*E00E*D00E + (d00E_00A*D00A + d00E_00B*D00B + d00E_00C*D00C + d00E_00D*D00D + d00E_11E*D11E + d00E_01E*D01E + d00E_00F*D00F + d00E_00G*D00G + d00E_00H*D00H + d00E_00I*D00I + d00E_00J*D00J);

ydot[43] =  - (s11E + (d11E_11A + d11E_11B + d11E_11C + d11E_11D + d11E_00E + d11E_01E + d11E_11F + d11E_11G + d11E_11H + d11E_11I + d11E_11J) + x11E) * D11E + 2*s11E*E11E*D11E + (d11E_11A*D11A + d11E_11B*D11B + d11E_11C*D11C + d11E_11D*D11D + d11E_00E*D00E + d11E_01E*D01E + d11E_11F*D11F + d11E_11G*D11G + d11E_11H*D11H + d11E_11I*D11I + d11E_11J*D11J);

ydot[44] =  - (s01E + (d01E_01A + d01E_01B + d01E_01C + d01E_01D + d01E_00E + d01E_11E + d01E_01F + d01E_01G + d01E_01H + d01E_01I + d01E_01J) + x01E) * D01E + 2*s01E*E01E*D01E + (d01E_01A*D01A + d01E_01B*D01B + d01E_01C*D01C + d01E_01D*D01D + d01E_00E*D00E + d01E_11E*D11E + d01E_01F*D01F + d01E_01G*D01G + d01E_01H*D01H + d01E_01I*D01I + d01E_01J*D01J);

ydot[45] =  - (s00F + (d00F_00A + d00F_00B + d00F_00C + d00F_00D + d00F_00E + d00F_11F + d00F_01F + d00F_00G + d00F_00H + d00F_00I + d00F_00J) + x00F) * D00F + 2*s00F*E00F*D00F + (d00F_00A*D00A + d00F_00B*D00B + d00F_00C*D00C + d00F_00D*D00D + d00F_00E*D00E + d00F_11F*D11F + d00F_01F*D01F + d00F_00G*D00G + d00F_00H*D00H + d00F_00I*D00I + d00F_00J*D00J);

ydot[46] =  - (s11F + (d11F_11A + d11F_11B + d11F_11C + d11F_11D + d11F_11E + d11F_00F + d11F_01F + d11F_11G + d11F_11H + d11F_11I + d11F_11J) + x11F) * D11F + 2*s11F*E11F*D11F + (d11F_11A*D11A + d11F_11B*D11B + d11F_11C*D11C + d11F_11D*D11D + d11F_11E*D11E + d11F_00F*D00F + d11F_01F*D01F + d11F_11G*D11G + d11F_11H*D11H + d11F_11I*D11I + d11F_11J*D11J);

ydot[47] =  - (s01F + (d01F_01A + d01F_01B + d01F_01C + d01F_01D + d01F_01E + d01F_00F + d01F_11F + d01F_01G + d01F_01H + d01F_01I + d01F_01J) + x01F) * D01F + 2*s01F*E01F*D01F + (d01F_01A*D01A + d01F_01B*D01B + d01F_01C*D01C + d01F_01D*D01D + d01F_01E*D01E + d01F_00F*D00F + d01F_11F*D11F + d01F_01G*D01G + d01F_01H*D01H + d01F_01I*D01I + d01F_01J*D01J);

ydot[48] =  - (s00G + (d00G_00A + d00G_00B + d00G_00C + d00G_00D + d00G_00E + d00G_00F + d00G_11G + d00G_01G + d00G_00H + d00G_00I + d00G_00J) + x00G) * D00G + 2*s00G*E00G*D00G + (d00G_00A*D00A + d00G_00B*D00B + d00G_00C*D00C + d00G_00D*D00D + d00G_00E*D00E + d00G_00F*D00F + d00G_11G*D11G + d00G_01G*D01G + d00G_00H*D00H + d00G_00I*D00I + d00G_00J*D00J);

ydot[49] =  - (s11G + (d11G_11A + d11G_11B + d11G_11C + d11G_11D + d11G_11E + d11G_11F + d11G_00G + d11G_01G + d11G_11H + d11G_11I + d11G_11J) + x11G) * D11G + 2*s11G*E11G*D11G + (d11G_11A*D11A + d11G_11B*D11B + d11G_11C*D11C + d11G_11D*D11D + d11G_11E*D11E + d11G_11F*D11F + d11G_00G*D00G + d11G_01G*D01G + d11G_11H*D11H + d11G_11I*D11I + d11G_11J*D11J);

ydot[50] =  - (s01G + (d01G_01A + d01G_01B + d01G_01C + d01G_01D + d01G_01E + d01G_01F + d01G_00G + d01G_11G + d01G_01H + d01G_01I + d01G_01J) + x01G) * D01G + 2*s01G*E01G*D01G + (d01G_01A*D01A + d01G_01B*D01B + d01G_01C*D01C + d01G_01D*D01D + d01G_01E*D01E + d01G_01F*D01F + d01G_00G*D00G + d01G_11G*D11G + d01G_01H*D01H + d01G_01I*D01I + d01G_01J*D01J);

ydot[51] =  - (s00H + (d00H_00A + d00H_00B + d00H_00C + d00H_00D + d00H_00E + d00H_00F + d00H_00G + d00H_11H + d00H_01H + d00H_00I + d00H_00J) + x00H) * D00H + 2*s00H*E00H*D00H + (d00H_00A*D00A + d00H_00B*D00B + d00H_00C*D00C + d00H_00D*D00D + d00H_00E*D00E + d00H_00F*D00F + d00H_00G*D00G + d00H_11H*D11H + d00H_01H*D01H + d00H_00I*D00I + d00H_00J*D00J);

ydot[52] =  - (s11H + (d11H_11A + d11H_11B + d11H_11C + d11H_11D + d11H_11E + d11H_11F + d11H_11G + d11H_00H + d11H_01H + d11H_11I + d11H_11J) + x11H) * D11H + 2*s11H*E11H*D11H + (d11H_11A*D11A + d11H_11B*D11B + d11H_11C*D11C + d11H_11D*D11D + d11H_11E*D11E + d11H_11F*D11F + d11H_11G*D11G + d11H_00H*D00H + d11H_01H*D01H + d11H_11I*D11I + d11H_11J*D11J);

ydot[53] =  - (s01H + (d01H_01A + d01H_01B + d01H_01C + d01H_01D + d01H_01E + d01H_01F + d01H_01G + d01H_00H + d01H_11H + d01H_01I + d01H_01J) + x01H) * D01H + 2*s01H*E01H*D01H + (d01H_01A*D01A + d01H_01B*D01B + d01H_01C*D01C + d01H_01D*D01D + d01H_01E*D01E + d01H_01F*D01F + d01H_01G*D01G + d01H_00H*D00H + d01H_11H*D11H + d01H_01I*D01I + d01H_01J*D01J);

ydot[54] =  - (s00I + (d00I_00A + d00I_00B + d00I_00C + d00I_00D + d00I_00E + d00I_00F + d00I_00G + d00I_00H + d00I_11I + d00I_01I + d00I_00J) + x00I) * D00I + 2*s00I*E00I*D00I + (d00I_00A*D00A + d00I_00B*D00B + d00I_00C*D00C + d00I_00D*D00D + d00I_00E*D00E + d00I_00F*D00F + d00I_00G*D00G + d00I_00H*D00H + d00I_11I*D11I + d00I_01I*D01I + d00I_00J*D00J);

ydot[55] =  - (s11I + (d11I_11A + d11I_11B + d11I_11C + d11I_11D + d11I_11E + d11I_11F + d11I_11G + d11I_11H + d11I_00I + d11I_01I + d11I_11J) + x11I) * D11I + 2*s11I*E11I*D11I + (d11I_11A*D11A + d11I_11B*D11B + d11I_11C*D11C + d11I_11D*D11D + d11I_11E*D11E + d11I_11F*D11F + d11I_11G*D11G + d11I_11H*D11H + d11I_00I*D00I + d11I_01I*D01I + d11I_11J*D11J);

ydot[56] =  - (s01I + (d01I_01A + d01I_01B + d01I_01C + d01I_01D + d01I_01E + d01I_01F + d01I_01G + d01I_01H + d01I_00I + d01I_11I + d01I_01J) + x01I) * D01I + 2*s01I*E01I*D01I + (d01I_01A*D01A + d01I_01B*D01B + d01I_01C*D01C + d01I_01D*D01D + d01I_01E*D01E + d01I_01F*D01F + d01I_01G*D01G + d01I_01H*D01H + d01I_00I*D00I + d01I_11I*D11I + d01I_01J*D01J);

ydot[57] =  - (s00J + (d00J_00A + d00J_00B + d00J_00C + d00J_00D + d00J_00E + d00J_00F + d00J_00G + d00J_00H + d00J_00I + d00J_11J + d00J_01J) + x00J) * D00J + 2*s00J*E00J*D00J + (d00J_00A*D00A + d00J_00B*D00B + d00J_00C*D00C + d00J_00D*D00D + d00J_00E*D00E + d00J_00F*D00F + d00J_00G*D00G + d00J_00H*D00H + d00J_00I*D00I + d00J_11J*D11J + d00J_01J*D01J);

ydot[58] =  - (s11J + (d11J_11A + d11J_11B + d11J_11C + d11J_11D + d11J_11E + d11J_11F + d11J_11G + d11J_11H + d11J_11I + d11J_00J + d11J_01J) + x11J) * D11J + 2*s11J*E11J*D11J + (d11J_11A*D11A + d11J_11B*D11B + d11J_11C*D11C + d11J_11D*D11D + d11J_11E*D11E + d11J_11F*D11F + d11J_11G*D11G + d11J_11H*D11H + d11J_11I*D11I + d11J_00J*D00J + d11J_01J*D01J);

ydot[59] =  - (s01J + (d01J_01A + d01J_01B + d01J_01C + d01J_01D + d01J_01E + d01J_01F + d01J_01G + d01J_01H + d01J_01I + d01J_00J + d01J_11J) + x01J) * D01J + 2*s01J*E01J*D01J + (d01J_01A*D01A + d01J_01B*D01B + d01J_01C*D01C + d01J_01D*D01D + d01J_01E*D01E + d01J_01F*D01F + d01J_01G*D01G + d01J_01H*D01H + d01J_01I*D01I + d01J_00J*D00J + d01J_11J*D11J);

}
