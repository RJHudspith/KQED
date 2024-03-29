/**
  Kernel when x and y are non-zero
**/
#include "KQED.h"

#include "chnr_dS.h"   // chnr_dS()
#include "chnr_dT.h"   // chnr_dT()
#include "chnr_dV.h"   // chnr_dV()
#include "corr_malloc.h"

#include "getff-new.h" // set_invariants()
#include "simd.h"

// hand-unrolled direct evaluation of the kernel
static void
CONSTRUCT_FULL_KERNEL( double *kp ,
		       struct STV x )
{
  const double *vv  = (const double*)x.Vv ;
  const double *sxv = (const double*)x.Sxv ;
  const double *syv = (const double*)x.Syv ;
  const double *txv = (const double*)x.Txv ;
  const double *tyv = (const double*)x.Tyv ;
  const double sxv0 = sxv[0]/4. ;
  const double sxv1 = sxv[1]/4. ;
  const double sxv2 = sxv[2]/4. ;
  const double sxv3 = sxv[3]/4. ;
  const double sxyv0 = (sxv[0]+syv[0])/4. ;
  const double sxyv1 = (sxv[1]+syv[1])/4. ;
  const double sxyv2 = (sxv[2]+syv[2])/4. ;
  const double sxyv3 = (sxv[3]+syv[3])/4. ;
  kp[0] = (+vv[26]+vv[31]-vv[38]-vv[55]);
  kp[1] = (+vv[10]+vv[15]+vv[34]+vv[51]-2*(txv[10]+tyv[10]+txv[15]+tyv[15]));
  kp[2] = (-vv[6]-vv[18]+2*(txv[6]+tyv[6]));
  kp[3] = (-vv[7]-vv[19]+2*(txv[7]+tyv[7]));
  kp[4] = (-vv[10]-vv[15]+vv[34]+vv[51]-2*(txv[10]+txv[15])+2*(txv[10]+tyv[10]+txv[15]+tyv[15]));
  kp[5] = (+vv[26]+vv[31]+vv[38]+vv[55]-2*(txv[26]+txv[31]));
  kp[6] = (+vv[2]-vv[22]+vv[42]-vv[47]+vv[59]+vv[62]-(txv[47]+txv[62])-2*(txv[2]+tyv[2]+txv[42]+sxv2)+(txv[47]+tyv[47])-(txv[62]+tyv[62]));
  kp[7] = (+vv[3]-vv[23]+vv[43]+vv[46]-vv[58]+vv[63]-(txv[43])-2*(txv[3]+tyv[3])-(txv[43]+tyv[43]+txv[58])+(txv[58]+tyv[58])-2*(txv[63]+sxv3));
  kp[8] = (+vv[6]-vv[18]-2*(tyv[6]));
  kp[9] = (-vv[2]-vv[22]+vv[42]+vv[47]+vv[59]-vv[62]+2*(txv[2]+tyv[2]+txv[22])-(2*txv[47]+tyv[47])+(2*txv[62]+tyv[62]));   
  kp[10] = (-vv[26]+vv[31]-vv[38]-vv[55]+2*(txv[38]+sxv1));
  kp[11] = (-vv[27]-vv[30]-vv[39]+vv[54]+(2*txv[39]+tyv[39]-tyv[54]));
  kp[12] = (+vv[7]-vv[19]-2*(tyv[7]));
  kp[13] = (-vv[3]-vv[23]-vv[43]+vv[46]+vv[58]+vv[63] +2*(txv[3]+tyv[3]+txv[23])+(2*txv[43]+tyv[43])-(2*txv[58]+tyv[58]));
  kp[14] = (-vv[27]-vv[30]+vv[39]-vv[54]-(tyv[39])+(2*txv[54]+tyv[54]));
  kp[15] = (+vv[26]-vv[31]-vv[38]-vv[55]+2*(txv[55]+sxv1));
  kp[16] = (-vv[10]-vv[15]-vv[34]-vv[51]+2*(txv[10]+txv[15]));
  kp[17] = (+vv[26]+vv[31]-vv[38]-vv[55]-2*(tyv[26]+tyv[31]));
  kp[18] = (+vv[2]-vv[22]-vv[42]+vv[47]-vv[59]-vv[62]+2*(txv[22]+tyv[22]+txv[42]+sxv2)-tyv[47]+2*txv[62]+tyv[62]);
  kp[19] = (+vv[3]-vv[23]-vv[43]-vv[46]+vv[58]-vv[63]
	    +2*(txv[23]+tyv[23])+(2*txv[43]+tyv[43])-(tyv[58])+2*(txv[63]+sxv3));
  kp[20] = (-vv[26]-vv[31]-vv[38]-vv[55]+2*(txv[26]+tyv[26]+txv[31]+tyv[31]));
  kp[21] = (-vv[10]-vv[15]+vv[34]+vv[51]);
  kp[22] = (+vv[6]+vv[18]-2*(txv[18]+tyv[18]));
  kp[23] = (+vv[7]+vv[19]-2*(txv[19]+tyv[19]));
  kp[24] = (+vv[2]+vv[22]-vv[42]-vv[47]-vv[59]+vv[62]+
	    -2*(txv[2]+txv[22]+tyv[22])+(2*txv[47]+tyv[47])-(2*txv[62]+tyv[62]));
  kp[25] = (+vv[6]-vv[18]+2*tyv[18]);
  kp[26] = (+vv[10]-vv[15]+vv[34]+vv[51]-2*(txv[34]+sxv0));
  kp[27] = (+vv[11]+vv[14]+vv[35]-vv[50]-(2*txv[35]+tyv[35])+(tyv[50]));
  kp[28] = (+vv[3]+vv[23]+vv[43]-vv[46]-vv[58]-vv[63]-2*(txv[3])-(2*txv[43]+tyv[43])-2*(txv[23]+tyv[23])+(2*txv[58]+tyv[58]));
  kp[29] = (+vv[7]-vv[19]+2*(tyv[19]));
  kp[30] = (+vv[11]+vv[14]-vv[35]+vv[50]+(tyv[35])-(2*txv[50]+tyv[50]));
  kp[31] = (-vv[10]+vv[15]+vv[34]+vv[51]-2*(txv[51]+sxv0));
  kp[32] = (+vv[6]+vv[18]-2*(txv[6]));
  kp[33] = (-vv[2]+vv[22]+vv[42]+vv[47]-vv[59]+vv[62]-2*(txv[22])-2*(txv[42]+tyv[42]+sxyv2)-(tyv[47])-(2*txv[62]+tyv[62]));
  kp[34] = (+vv[26]-vv[31]-vv[38]+vv[55]+2*(tyv[38]-sxv1+sxyv1));
  kp[35] = (+vv[27]+vv[30]-vv[39]-vv[54]+(tyv[39])+(tyv[54]));
  kp[36] = (-vv[2]+vv[22]-vv[42]-vv[47]+vv[59]-vv[62]+2*(txv[2]+txv[42]+tyv[42]+sxyv2)+(tyv[47])+(2*txv[62]+tyv[62]));
  kp[37] = (-vv[6]-vv[18]+2*txv[18]);
  kp[38] = (-vv[10]+vv[15]+vv[34]-vv[51]-2*(tyv[34]-sxv0+sxyv0));
  kp[39] = (-vv[11]-vv[14]+vv[35]+vv[50]-tyv[35]-tyv[50]);
  kp[40] = (+vv[26]+vv[31]+vv[38]-vv[55]-2*(txv[38]+tyv[38]+sxyv1));
  kp[41] = (-vv[10]-vv[15]-vv[34]+vv[51]+2*(txv[34]+tyv[34]+sxyv0));
  kp[42] = (+vv[6]-vv[18]);
  kp[43] = (+vv[7]-vv[19]);
  kp[44] = (-vv[27]+vv[30]+vv[39]+vv[54]-(tyv[39])-(2*txv[54]+tyv[54]));
  kp[45] = (+vv[11]-vv[14]-vv[35]-vv[50]+tyv[35]+2*txv[50]+tyv[50]);
  kp[46] = (-vv[7]+vv[19]);
  kp[47] = (+vv[6]-vv[18]);
  kp[48] = (+vv[7]+vv[19]-2*txv[7]);
  kp[49] = (-vv[3]+vv[23]+vv[43]-vv[46]+vv[58]+vv[63]-2*(txv[23]+txv[43]+txv[63])-tyv[43]-tyv[58]-2*(tyv[63]+sxyv3));
  kp[50] = (+vv[27]+vv[30]-vv[39]-vv[54]+tyv[39]+tyv[54]);
  kp[51] = (-vv[26]+vv[31]+vv[38]-vv[55]-2*(+sxv1-tyv[55]-sxyv1));
  kp[52] = (-vv[3]+vv[23]-vv[43]+vv[46]-vv[58]-vv[63]+2*(txv[3]+txv[63]+tyv[63]+sxyv3)+(2*txv[43]+tyv[43])+(tyv[58]));
  kp[53] = (-vv[7]-vv[19]+2*(txv[19]));
  kp[54] = (-vv[11]-vv[14]+vv[35]+vv[50]-tyv[35]-tyv[50]);
  kp[55] = (+vv[10]-vv[15]-vv[34]+vv[51]-2*(tyv[51]-sxv0+sxyv0));  
  kp[56] = (+vv[27]-vv[30]+vv[39]+vv[54]-2*txv[39]-tyv[39]-tyv[54]);
  kp[57] = (-vv[11]+vv[14]-vv[35]-vv[50]+2*txv[35]+tyv[35]+tyv[50]);
  kp[58] = (+vv[7]-vv[19]);
  kp[59] = (-vv[6]+vv[18]);
  kp[60] = (+vv[26]+vv[31]-vv[38]+vv[55]-2*(txv[55]+tyv[55]+sxyv1));
  kp[61] = (-vv[10]-vv[15]+vv[34]-vv[51]+2*(txv[51]+tyv[51]+sxyv0));
  kp[62] = (+vv[6]-vv[18]);
  kp[63] = (+vv[7]-vv[19]);
  kp[64] = (-vv[25]+vv[37]+vv[47]-vv[59]);
  kp[65] = (-vv[9]-vv[33]+2*(txv[9]+tyv[9]));
  kp[66] = (+vv[5]+vv[15]+vv[17]+vv[51]-2*(txv[5]+tyv[5]+txv[15]+tyv[15]));
  kp[67] = (-vv[11]-vv[35]+2*(txv[11]+tyv[11]));
  kp[68] = (+vv[9]-vv[33]-2*tyv[9]);
  kp[69] = (-vv[25]-vv[37]+vv[47]-vv[59]+2*(txv[25]+sxv2));
  kp[70] = (-vv[1]+vv[21]+vv[31]-vv[41]+vv[55]-vv[61]+2*(txv[1]+tyv[1])-(2*txv[31]+tyv[31])+2*(txv[41])+(2*txv[61]+tyv[61]));
  kp[71] = (-vv[27]-vv[39]-vv[45]+vv[57]+(2*txv[27]+tyv[27])-tyv[57]);
  kp[72] = (-vv[5]-vv[15]+vv[17]+vv[51]+2*(tyv[5]+tyv[15]));
  kp[73] = (+vv[1]+vv[21]-vv[31]-vv[41]+vv[55]+vv[61]-2*(txv[1]+tyv[1])-2*(txv[21]+sxv1)+(tyv[31])-(2*txv[61]+tyv[61]));
  kp[74] = (+vv[25]+vv[37]+vv[47]+vv[59]-2*(txv[37]+txv[47]));
  kp[75] = (+vv[3]+vv[23]+vv[29]-vv[43]-vv[53]+vv[63]-2*(txv[3]+tyv[3])-(2*txv[23]+tyv[23])+(tyv[53])-2*(txv[63]+sxv3));
  kp[76] = (+vv[11]-vv[35]-2*tyv[11]);
  kp[77] = (+vv[27]-vv[39]-vv[45]-vv[57]-(tyv[27])+(2*txv[57]+tyv[57]));
  kp[78] = (-vv[3]-vv[23]+vv[29]-vv[43]+vv[53]+vv[63]+2*(txv[3]+tyv[3]+txv[23]+txv[43]-txv[53])+tyv[23]-tyv[53]);
  kp[79] = (-vv[25]+vv[37]-vv[47]-vv[59]+2*(txv[59]+sxv2));
  kp[80] = (+vv[9]+vv[33]-2*txv[9]);
  kp[81] = (-vv[25]+vv[37]-vv[47]+vv[59]+2*(tyv[25]-sxv2+sxyv2));
  kp[82] = (-vv[1]+vv[21]+vv[31]+vv[41]-vv[55]+vv[61]-2*(txv[21]+tyv[21]+sxyv1)-(tyv[31])-2*(txv[41])-(2*txv[61]+tyv[61]));
  kp[83] = (-vv[27]+vv[39]+vv[45]-vv[57]+(tyv[27]+tyv[57]));
  kp[84] = (+vv[25]+vv[37]+vv[47]-vv[59]-2*(txv[25]+tyv[25]+sxyv2));
  kp[85] = (+vv[9]-vv[33]);
  kp[86] = (-vv[5]-vv[17]-vv[15]+vv[51]+2*(txv[17]+tyv[17]+sxyv0));
  kp[87] = (+vv[11]-vv[35]);
  kp[88] = (-vv[1]-vv[21]-vv[31]+vv[41]+vv[55]-vv[61]+2*(txv[1])+2*(txv[21]+tyv[21]+sxyv1)+(tyv[31])+(2*txv[61]+tyv[61]));
  kp[89] = (-vv[5]+vv[15]+vv[17]-vv[51]+2*(sxv0-tyv[17]-sxyv0));
  kp[90] = (-vv[9]-vv[33]+2*txv[33]);
  kp[91] = (-vv[7]-vv[13]+vv[19]+vv[49]-(tyv[19]+tyv[49]));
  kp[92] = (+vv[27]-vv[39]+vv[45]+vv[57]-(tyv[27]+2*txv[57]+tyv[57]));
  kp[93] = (-vv[11]+vv[35]);
  kp[94] = (+vv[7]-vv[13]-vv[19]-vv[49]+tyv[19]+2*txv[49]+tyv[49]);
  kp[95] = (+vv[9]-vv[33]);
  kp[96] = (-vv[5]-vv[15]-vv[17]-vv[51]+2*(txv[5]+txv[15]));
  kp[97] = (+vv[1]-vv[21]+vv[31]-vv[41]-vv[55]-vv[61]+2*(txv[21]+sxv1)-(tyv[31])+2*(txv[41]+tyv[41])+(2*txv[61]+tyv[61]));
  kp[98] = (-vv[25]+vv[37]+vv[47]-vv[59]-2*(tyv[37]+tyv[47]));
  kp[99] = (+vv[3]-vv[23]-vv[29]-vv[43]+vv[53]-vv[63]+(2*txv[23]+tyv[23])+2*(txv[43]+tyv[43])-(tyv[53])+2*(txv[63]+sxv3));
  kp[100] = (+vv[1]-vv[21]-vv[31]+vv[41]-vv[55]+vv[61]-2*(txv[1])+(2*txv[31]+tyv[31])-2*(txv[41]+tyv[41])-(2*txv[61]+tyv[61]));
  kp[101] = (+vv[5]-vv[15]+vv[17]+vv[51]-2*(txv[17]+sxv0));
  kp[102] = (+vv[9]-vv[33]+2*tyv[33]);
  kp[103] = (+vv[7]+vv[13]+vv[19]-vv[49]-2*txv[19]-tyv[19]+tyv[49]);
  kp[104] = (-vv[25]-vv[37]-vv[47]-vv[59]+2*(txv[37]+txv[47]+tyv[37]+tyv[47]));
  kp[105] = (+vv[9]+vv[33]-2*(txv[33]+tyv[33]));
  kp[106] = (-vv[5]-vv[15]+vv[17]+vv[51]);
  kp[107] = (+vv[11]+vv[35]-2*(txv[35]+tyv[35]));
  kp[108] = (+vv[3]+vv[23]-vv[29]+vv[43]-vv[53]-vv[63]-2*(txv[3]+txv[23]+txv[43]-txv[53])-(tyv[23]+2*tyv[43]-tyv[53]));
  kp[109] = (+vv[7]+vv[13]-vv[19]+vv[49]+tyv[19]-2*txv[49]-tyv[49]);
  kp[110] = (+vv[11]-vv[35]+2*tyv[35]);
  kp[111] = (-vv[5]+vv[15]+vv[17]+vv[51]-2*(txv[51]+sxv0));
  kp[112] = (+vv[11]+vv[35]-2*(txv[11]));
  kp[113] = (-vv[27]+vv[39]+vv[45]-vv[57]+tyv[27]+tyv[57]);
  kp[114] = (-vv[3]+vv[23]-vv[29]+vv[43]+vv[53]+vv[63]-2*(txv[23]+txv[43]+txv[63])-( tyv[23]+tyv[53]+2*(tyv[63]+sxyv3) ));
  kp[115] = (+vv[25]-vv[37]+vv[47]-vv[59]+2*(tyv[59]-sxv2+sxyv2));
  kp[116] = (+vv[27]+vv[39]-vv[45]+vv[57]-(2*txv[27]+tyv[27]+tyv[57]));
  kp[117] = (+vv[11]-vv[35]);
  kp[118] = (-vv[7]+vv[13]-vv[19]-vv[49]+2*txv[19]+tyv[19]+tyv[49]);
  kp[119] = (-vv[9]+vv[33]);
  kp[120] = (-vv[3]-vv[23]+vv[29]+vv[43]-vv[53]-vv[63]+2*(txv[3]+txv[23]+txv[63])+tyv[23]+tyv[53]+2*(tyv[63]+sxyv3));
  kp[121] = (-vv[7]-vv[13]+vv[19]+vv[49]-tyv[19]-tyv[49]);
  kp[122] = (-vv[11]-vv[35]+2*(txv[35]));
  kp[123] = (+vv[5]-vv[15]-vv[17]+vv[51]-2*(tyv[51]-sxv0+sxyv0));
  kp[124] = (-vv[25]+vv[37]+vv[47]+vv[59]-2*(txv[59]+tyv[59]+sxyv2));
  kp[125] = (+vv[9]-vv[33]);
  kp[126] = (-vv[5]-vv[15]+vv[17]-vv[51]+2*(txv[51]+tyv[51]+sxyv0));
  kp[127] = (+vv[11]-vv[35]);
  kp[128] = (-vv[29]-vv[46]+vv[53]+vv[58]);
  kp[129] = (-vv[13]-vv[49]+2*(txv[13]+tyv[13]));
  kp[130] = (-vv[14]-vv[50]+2*(txv[14]+tyv[14]));
  kp[131] = (+vv[5]+vv[10]+vv[17]+vv[34]-2*(txv[5]+tyv[5]+txv[10]+tyv[10]));
  kp[132] = (+vv[13]-vv[49]-2*tyv[13]);
  kp[133] = (-vv[29]-vv[46]-vv[53]+vv[58]+2*(txv[29]+sxv3));
  kp[134] = (-vv[30]+vv[45]-vv[54]-vv[57]+2*txv[30]+tyv[30]-tyv[45]);
  kp[135] = (-vv[1]+vv[21]+vv[26]+vv[38]-vv[41]-vv[61]+2*(txv[1]-txv[26]+txv[41]+txv[61])+2*tyv[1]-tyv[26]+tyv[41] );
  kp[136] = (+vv[14]-vv[50]-2*tyv[14]);
  kp[137] = (+vv[30]-vv[45]-vv[54]-vv[57]+2*txv[45]-tyv[30]+tyv[45]);
  kp[138] = (-vv[29]-vv[46]+vv[53]-vv[58]+2*(txv[46]+sxv3));
  kp[139] = (-vv[2]-vv[22]+vv[25]+vv[37]+vv[42]-vv[62]+2*(txv[2]+txv[22]-txv[37]+txv[62])+2*tyv[2]+tyv[22]-tyv[37]);
  kp[140] = (-vv[5]-vv[10]+vv[17]+vv[34]+2*(tyv[5]+tyv[10]));
  kp[141] = (+vv[1]+vv[21]-vv[26]+vv[38]+vv[41]-vv[61]
	       -2*(txv[1]+txv[21]+sxv1+txv[41]) -2*tyv[1]+tyv[26]-tyv[41]);
  kp[142] = (+vv[2]+vv[22]+vv[25]-vv[37]+vv[42]-vv[62]
	       -2*(txv[2]+txv[22]+txv[42]+sxv2)-2*tyv[2]-tyv[22]+tyv[37]) ; 
  kp[143] = (+vv[29]+vv[46]+vv[53]+vv[58]-2*(txv[53]+txv[58]));
  kp[144] = (+vv[13]+vv[49]-2*(txv[13]));
  kp[145] = (-vv[29]+vv[46]+vv[53]-vv[58]+2*(tyv[29]-sxv3+sxyv3));
  kp[146] = (-vv[30]-vv[45]+vv[54]+vv[57]+tyv[30]+tyv[45]);
  kp[147] = (-vv[1]+vv[21]+vv[26]-vv[38]+vv[41]+vv[61]-2*(txv[21]+txv[41]+txv[61]+sxyv1)-2*tyv[21]-tyv[26]-tyv[41]);
  kp[148] = (+vv[29]-vv[46]+vv[53]+vv[58]-2*(txv[29]+tyv[29]+sxyv3));
  kp[149] = (+vv[13]-vv[49]);
  kp[150] = (+vv[14]-vv[50]);
  kp[151] = (-vv[5]-vv[10]-vv[17]+vv[34]+2*(txv[17]+tyv[17]+sxyv0));
  kp[152] = (+vv[30]+vv[45]-vv[54]+vv[57]-tyv[30]-2*txv[45]-tyv[45]);
  kp[153] = (-vv[14]+vv[50]);
  kp[154] = (+vv[13]-vv[49]);
  kp[155] = (+vv[6]-vv[9]-vv[18]-vv[33]+2*txv[33]+tyv[18]+tyv[33]);
  kp[156] = (-vv[1]-vv[21]-vv[26]+vv[38]-vv[41]+vv[61]+2*(txv[1]+txv[21]+txv[41]+sxyv1)+2*tyv[21]+tyv[26]+tyv[41]) ;
  kp[157] = (-vv[5]+vv[10]+vv[17]-vv[34]-2*(tyv[17]-sxv0+sxyv0));
  kp[158] = (-vv[6]-vv[9]+vv[18]+vv[33]-(tyv[18]+tyv[33]));
  kp[159] = (-vv[13]-vv[49]+2*(txv[49]));
  kp[160] = (+vv[14]+vv[50]-2*(txv[14]));
  kp[161] = (-vv[30]-vv[45]+vv[54]+vv[57]+tyv[30]+tyv[45]);
  kp[162] = (+vv[29]-vv[46]-vv[53]+vv[58]+2*(tyv[46]-sxv3+sxyv3));
  kp[163] = (-vv[2]+vv[22]-vv[25]+vv[37]+vv[42]+vv[62]-2*(txv[22]+txv[42]+txv[62]+sxyv2)-tyv[22]-tyv[37]-2*tyv[42] ) ;
  kp[164] = (+vv[30]+vv[45]+vv[54]-vv[57]-2*txv[30]-tyv[30]-tyv[45]);
  kp[165] = (+vv[14]-vv[50]);
  kp[166] = (-vv[13]+vv[49]);
  kp[167] = (-vv[6]+vv[9]-vv[18]-vv[33]+2*txv[18]+tyv[18]+tyv[33]);
  kp[168] = (-vv[29]+vv[46]+vv[53]+vv[58]-2*(txv[46]+tyv[46]+sxyv3));
  kp[169] = (+vv[13]-vv[49]);
  kp[170] = (+vv[14]-vv[50]);
  kp[171] = (-vv[5]-vv[10]+vv[17]-vv[34]+2*(txv[34]+tyv[34]+sxyv0));
  kp[172] = (-vv[2]-vv[22]+vv[25]-vv[37]-vv[42]+vv[62]+2*(txv[2]+txv[22]+txv[42]+sxyv2)+tyv[22]+tyv[37]+2*tyv[42] ) ;
  kp[173] = (-vv[6]-vv[9]+vv[18]+vv[33]-(tyv[18]+tyv[33]));
  kp[174] = (+vv[5]-vv[10]-vv[17]+vv[34]-2*(tyv[34]-sxv0+sxyv0));
  kp[175] = (-vv[14]-vv[50]+2*(txv[50]));
  kp[176] = (-vv[5]-vv[10]-vv[17]-vv[34]+2*(txv[5]+txv[10]));
  kp[177] = (+vv[1]-vv[21]+vv[26]-vv[38]-vv[41]-vv[61]+2*(txv[21]+txv[41]+txv[61]+sxv1)-tyv[26]+2*tyv[61]+tyv[41] ) ;
  kp[178] = (+vv[2]-vv[22]-vv[25]+vv[37]-vv[42]-vv[62]+2*(txv[22]+txv[42]+txv[62]+sxv2)+tyv[22]-tyv[37]+2*tyv[62] ) ;
  kp[179] = (-vv[29]-vv[46]+vv[53]+vv[58]-2*(tyv[53]+tyv[58]));
  kp[180] = (+vv[1]-vv[21]-vv[26]-vv[38]+vv[41]+vv[61]-2*(txv[1]-txv[26]+txv[41]+txv[61])+tyv[26]-tyv[41]-2*tyv[61]);
  kp[181] = (+vv[5]-vv[10]+vv[17]+vv[34]-2*(txv[17]+sxv0));
  kp[182] = (+vv[6]+vv[9]+vv[18]-vv[33]-(2*txv[18]+tyv[18]-tyv[33]));
  kp[183] = (+vv[13]-vv[49]+2*tyv[49]);
  kp[184] = (+vv[2]+vv[22]-vv[25]-vv[37]-vv[42]+vv[62]-2*(txv[2]+txv[22]-txv[37]+txv[62])-tyv[22]+tyv[37]-2*tyv[62] ) ;
  kp[185] = (+vv[9]+vv[6]-vv[18]+vv[33]+tyv[18]-2*txv[33]-tyv[33]);
  kp[186] = (-vv[5]+vv[10]+vv[17]+vv[34]-2*(txv[34]+sxv0));
  kp[187] = (+vv[14]-vv[50]+2*tyv[50]);
  kp[188] = (-vv[29]-vv[46]-vv[53]-vv[58]+2*(txv[53]+tyv[53]+txv[58]+tyv[58]));
  kp[189] = (+vv[13]+vv[49]-2*(txv[49]+tyv[49]));
  kp[190] = (+vv[14]+vv[50]-2*(txv[50]+tyv[50]));
  kp[191] = (-vv[5]-vv[10]+vv[17]+vv[34]);
  kp[192] = (+vv[24]-vv[36]);
  kp[193] = (+vv[8]+vv[32]+vv[47]-vv[59]-2*(txv[8]+tyv[8]+sxyv2));
  kp[194] = (-vv[4]-vv[16]-vv[31]+vv[55]+2*(txv[4]+tyv[4]+sxyv1));
  kp[195] = (+vv[27]-vv[39]);
  kp[196] = (-vv[8]+vv[32]-vv[47]+vv[59]-2*(txv[8]+sxv2-txv[8]-tyv[8]-sxyv2));
  kp[197] = (+vv[24]+vv[36]-2*txv[24]);
  kp[198] = (+vv[0]+vv[15]-vv[20]+vv[40]-vv[51]+vv[60]-2*(txv[0]+txv[40]+txv[60]+sxyv0)-2*tyv[0]-tyv[15]-tyv[60]) ;
  kp[199] = (-vv[11]+vv[35]+vv[44]-vv[56]+tyv[11]+tyv[56]);
  kp[200] = (+vv[4]-vv[16]+vv[31]-vv[55]-2*(tyv[4]-sxv1+sxyv1));
  kp[201] = (-vv[0]-vv[15]-vv[20]+vv[40]+vv[51]-vv[60]+2*(txv[0]+txv[20]+txv[60]+sxyv0)+2*tyv[0]+tyv[15]+tyv[60]);
  kp[202] = (-vv[24]-vv[36]+2*(txv[36])); 
  kp[203] = (+vv[7]-vv[19]-vv[28]+vv[52]-tyv[7]-tyv[52]);
  kp[204] = (-vv[27]+vv[39]);
  kp[205] = (+vv[11]-vv[35]+vv[44]+vv[56]-tyv[11]-2*txv[56]-tyv[56]);
  kp[206] = (-vv[7]+vv[19]-vv[28]-vv[52]+tyv[7]+2*txv[52]+tyv[52]);
  kp[207] = (+vv[24]-vv[36]);
  kp[208] = (-vv[8]-vv[32]+vv[47]-vv[59]+2*(txv[8]+sxv2));
  kp[209] = (+vv[24]-vv[36]-2*tyv[24]);
  kp[210] = (+vv[0]+vv[15]-vv[20]-vv[40]+vv[51]-vv[60]+2*(-txv[15]+txv[20]+txv[40]+txv[60])-tyv[15]+2*tyv[20]+tyv[60]);
  kp[211] = (-vv[11]-vv[35]-vv[44]+vv[56]+2*txv[11]+tyv[11]-tyv[56]);
  kp[212] = (-vv[24]-vv[36]+2*(txv[24]+tyv[24]));
  kp[213] = (-vv[8]+vv[32]+vv[47]-vv[59]);
  kp[214] = (+vv[4]+vv[16]+vv[31]+vv[55]-2*(txv[31]+tyv[31]+txv[16]+tyv[16]));
  kp[215] = (-vv[27]-vv[39]+2*(txv[27]+tyv[27]));
  kp[216] = (+vv[0]-vv[15]+vv[20]-vv[40]+vv[51]+vv[60]-2*(txv[0]+txv[20]+txv[60]+sxv0)+tyv[15]-2*tyv[20]-tyv[60]) ;
  kp[217] = (+vv[4]-vv[16]-vv[31]+vv[55]+2*(tyv[16]+tyv[31]));
  kp[218] = (+vv[8]+vv[32]+vv[47]+vv[59]-2*(txv[32]+txv[47]));
  kp[219] = (+vv[3]+vv[12]+vv[23]-vv[43]-vv[48]+vv[63]-2*(txv[3]+txv[23]+txv[63]+sxv3)-tyv[3]-2*tyv[23]+tyv[48]);
  kp[220] = (+vv[11]-vv[35]-vv[44]-vv[56]+2*txv[56]-tyv[11]+tyv[56]);
  kp[221] = (+vv[27]-vv[39]-2*tyv[27]);
  kp[222] = (-vv[3]+vv[12]-vv[23]-vv[43]+vv[48]+vv[63]+2*(txv[3]+txv[23]+txv[43]-txv[48])+tyv[3]+2*tyv[23]-tyv[48]);
  kp[223] = (-vv[8]+vv[32]-vv[47]-vv[59]+2*(txv[59]+sxv2));
  kp[224] = (+vv[4]+vv[16]-vv[31]+vv[55]-2*(txv[4]+sxv1));  
  kp[225] = (-vv[0]-vv[15]+vv[20]+vv[40]-vv[51]+vv[60]-2*(-txv[15]+txv[20]+txv[40]+txv[60])+tyv[15]-2*tyv[40]-tyv[60]);
  kp[226] = (+vv[24]-vv[36]+2*tyv[36]);
  kp[227] = (+vv[7]+vv[19]+vv[28]-vv[52]-2*txv[7]-tyv[7]+tyv[52]);
  kp[228] = (-vv[0]+vv[15]+vv[20]-vv[40]-vv[51]-vv[60]+2*(txv[0]+txv[40]+txv[60]+sxv0)-tyv[15]+2*tyv[40]+tyv[60]);
  kp[229] = (-vv[4]-vv[16]-vv[31]-vv[55]+2*(txv[16]+txv[31]));
  kp[230] = (-vv[8]+vv[32]+vv[47]-vv[59]-2*(tyv[32]+tyv[47]));
  kp[231] = (-vv[3]-vv[12]+vv[23]-vv[43]+vv[48]-vv[63]+2*(txv[3]+txv[43]+txv[63]+sxv3)+tyv[3]+2*tyv[43]-tyv[48]);
  kp[232] = (+vv[24]+vv[36]-2*(txv[36]+tyv[36]));
  kp[233] = (-vv[8]-vv[32]-vv[47]-vv[59]+2*(txv[32]+tyv[32]+txv[47]+tyv[47]));
  kp[234] = (+vv[4]-vv[16]-vv[31]+vv[55]);
  kp[235] = (+vv[27]+vv[39]-2*(txv[39]+tyv[39]));
  kp[236] = (-vv[7]+vv[19]+vv[28]+vv[52]+tyv[7]-2*txv[52]-tyv[52]);
  kp[237] = (+vv[3]-vv[12]+vv[23]+vv[43]-vv[48]-vv[63]-2*(txv[3]+txv[23]+txv[43]-txv[48])-tyv[3]-2*tyv[43]+tyv[48]);
  kp[238] = (+vv[27]-vv[39]+2*tyv[39]);
  kp[239] = (+vv[4]-vv[16]+vv[31]+vv[55]-2*(txv[55]+sxv1));
  kp[240] = (+vv[27]-vv[39]);
  kp[241] = (+vv[11]+vv[35]-vv[44]+vv[56]-2*txv[11]-tyv[11]-tyv[56]);
  kp[242] = (-vv[7]-vv[19]+vv[28]-vv[52]+2*txv[7]+tyv[7]+tyv[52]);
  kp[243] = (-vv[24]+vv[36]);
  kp[244] = (-vv[11]+vv[35]+vv[44]-vv[56]+tyv[11]+tyv[56]);
  kp[245] = (+vv[27]+vv[39]-2*txv[27]);
  kp[246] = (+vv[3]-vv[12]-vv[23]+vv[43]+vv[48]+vv[63]-2*(txv[3]+txv[43]+txv[63]+sxyv3)-tyv[3]-tyv[48]-2*tyv[63]);
  kp[247] = (+vv[8]-vv[32]+vv[47]-vv[59]+2*(tyv[59]-sxv2+sxyv2));
  kp[248] = (+vv[7]-vv[19]-vv[28]+vv[52]-tyv[7]-tyv[52]);  
  kp[249] = (-vv[3]+vv[12]-vv[23]+vv[43]-vv[48]-vv[63]+2*(txv[3]+txv[23]+txv[63]+sxyv3)+tyv[3]+tyv[48]+2*tyv[63]);
  kp[250] = (-vv[27]-vv[39]+2*txv[39]);
  kp[251] = (-vv[4]+vv[16]-vv[31]+vv[55]-2*(tyv[55]-sxv1+sxyv1));
  kp[252] = (+vv[24]-vv[36]);
  kp[253] = (-vv[8]+vv[32]+vv[47]+vv[59]-2*(txv[59]+tyv[59]+sxyv2));
  kp[254] = (+vv[4]-vv[16]-vv[31]-vv[55]+2*(txv[55]+tyv[55]+sxyv1));
  kp[255] = (+vv[27]-vv[39]);
  kp[256] = (+vv[28]-vv[52]);
  kp[257] = (+vv[12]-vv[46]+vv[48]+vv[58]-2*(txv[12]+tyv[12]+sxyv3));
  kp[258] = (+vv[30]-vv[54]);
  kp[259] = (-vv[4]-vv[16]-vv[26]+vv[38]+2*(txv[4]+tyv[4]+sxyv1));
  kp[260] = (-vv[12]+vv[46]+vv[48]-vv[58]+2*(tyv[12]-sxv3+sxyv3));
  kp[261] = (+vv[28]+vv[52]-2*txv[28]);
  kp[262] = (-vv[14]-vv[44]+vv[50]+vv[56]+tyv[14]+tyv[44]); 
  kp[263] = (+vv[0]+vv[10]-vv[20]-vv[34]+vv[40]+vv[60]-2*(txv[0]+txv[40]+txv[60]+sxyv0)-2*tyv[0]-tyv[10]-tyv[40]);
  kp[264] = (-vv[30]+vv[54]);
  kp[265] = (+vv[14]+vv[44]-vv[50]+vv[56]-tyv[14]-2*txv[44]-tyv[44]);
  kp[266] = (+vv[28]-vv[52]);
  kp[267] = (-vv[6]+vv[18]-vv[24]-vv[36]+tyv[6]+2*txv[36]+tyv[36]);
  kp[268] = (+vv[4]-vv[16]+vv[26]-vv[38]-2*(tyv[4]-sxv1+sxyv1)); 
  kp[269] = (-vv[0]-vv[10]-vv[20]+vv[34]-vv[40]+vv[60]+2*(txv[0]+txv[20]+txv[40]+sxyv0)+2*tyv[0]+tyv[10]+tyv[40]);
  kp[270] = (+vv[6]-vv[18]-vv[24]+vv[36]-tyv[6]-tyv[36]);
  kp[271] = (-vv[28]-vv[52]+2*(txv[52]));
  kp[272] = (-vv[12]-vv[46]-vv[48]+vv[58]+2*(txv[12]+sxv3));
  kp[273] = (+vv[28]-vv[52]-2*tyv[28]);
  kp[274] = (-vv[14]+vv[44]-vv[50]-vv[56]+2*txv[14]+tyv[14]-tyv[44]);
  kp[275] = (+vv[0]+vv[10]-vv[20]+vv[34]-vv[40]-vv[60]+2*(-txv[10]+txv[20]+txv[40]+txv[60])-tyv[10]+2*tyv[20]+tyv[40]);
  kp[276] = (-vv[28]-vv[52]+2*(txv[28]+tyv[28]));
  kp[277] = (-vv[12]-vv[46]+vv[48]+vv[58]);
  kp[278] = (-vv[30]-vv[54]+2*(txv[30]+tyv[30]));
  kp[279] = (+vv[26]+vv[38]+vv[4]+vv[16]-2*(txv[16]+txv[26]+tyv[16]+tyv[26]));
  kp[280] = (+vv[14]-vv[44]-vv[50]-vv[56]-tyv[14]+2*txv[44]+tyv[44]);
  kp[281] = (+vv[30]-vv[54]-2*tyv[30]);
  kp[282] = (-vv[12]-vv[46]+vv[48]-vv[58]+2*(txv[46]+sxv3));
  kp[283] = (-vv[2]+vv[8]-vv[22]+vv[32]+vv[42]-vv[62]+2*(txv[2]+txv[22]-txv[32]+txv[62])+tyv[2]+2*tyv[22]-tyv[32]);
  kp[284] = (+vv[0]-vv[10]+vv[20]+vv[34]+vv[40]-vv[60]-2*(txv[0]+txv[20]+txv[40]+sxv0)+tyv[10]-2*tyv[20]-tyv[40]);
  kp[285] = (+vv[4]-vv[16]-vv[26]+vv[38]+2*(tyv[16]+tyv[26]));  
  kp[286] = (+vv[2]+vv[8]+vv[22]-vv[32]+vv[42]-vv[62]-2*(txv[2]+txv[22]+txv[42]+sxv2)-tyv[2]-2*tyv[22]+tyv[32]);
  kp[287] = (+vv[12]+vv[46]+vv[48]+vv[58]-2*(txv[48]+txv[58]));
  kp[288] = (+vv[30]-vv[54]);
  kp[289] = (+vv[14]+vv[44]+vv[50]-vv[56]-2*txv[14]-tyv[14]-tyv[44]);
  kp[290] = (-vv[28]+vv[52]);
  kp[291] = (-vv[6]-vv[18]+vv[24]-vv[36]+2*txv[6]+tyv[6]+tyv[36]);
  kp[292] = (-vv[14]-vv[44]+vv[50]+vv[56]+tyv[14]+tyv[44]);
  kp[293] = (+vv[30]+vv[54]-2*txv[30]);
  kp[294] = (+vv[12]-vv[46]-vv[48]+vv[58]+2*(tyv[46]-sxv3+sxyv3));
  kp[295] = (+vv[2]-vv[8]-vv[22]+vv[32]+vv[42]+vv[62]-2*(txv[2]+txv[42]+txv[62]+sxyv2)-tyv[2]-tyv[32]-2*tyv[42]);
  kp[296] = (+vv[28]-vv[52]);
  kp[297] = (-vv[12]+vv[46]+vv[48]+vv[58]-2*(txv[46]+tyv[46]+sxyv3));
  kp[298] = (+vv[30]-vv[54]);
  kp[299] = (+vv[4]-vv[16]-vv[26]-vv[38]+2*(txv[38]+tyv[38]+sxyv1));
  kp[300] = (+vv[6]-vv[18]-vv[24]+vv[36]-tyv[6]-tyv[36]);
  kp[301] = (-vv[2]+vv[8]-vv[22]-vv[32]-vv[42]+vv[62]+2*(txv[2]+txv[22]+txv[42]+sxyv2)+tyv[2]+tyv[32]+2*tyv[42]);
  kp[302] = (-vv[4]+vv[16]-vv[26]+vv[38]-2*(tyv[38]-sxv1+sxyv1));
  kp[303] = (-vv[30]-vv[54]+2*(txv[54]));
  kp[304] = (+vv[4]+vv[16]-vv[26]+vv[38]-2*(txv[4]+sxv1));
  kp[305] = (-vv[0]-vv[10]+vv[20]-vv[34]+vv[40]+vv[60]-2*(-txv[10]+txv[20]+txv[40]+txv[60])+tyv[10]-tyv[40]-2*tyv[60]);
  kp[306] = (+vv[6]+vv[18]+vv[24]-vv[36]-2*txv[6]-tyv[6]+tyv[36]);
  kp[307] = (+vv[28]-vv[52]+2*tyv[52]);
  kp[308] = (-vv[0]+vv[10]+vv[20]-vv[34]-vv[40]-vv[60]+2*(txv[0]+txv[40]+txv[60]+sxv0)-tyv[10]+tyv[40]+2*tyv[60]);
  kp[309] = (-vv[4]-vv[16]-vv[26]-vv[38]+2*(txv[16]+txv[26]));
  kp[310] = (-vv[2]-vv[8]+vv[22]+vv[32]-vv[42]-vv[62]+2*(txv[2]+txv[42]+txv[62]+sxv2)+tyv[2]-tyv[32]+2*tyv[62]);
  kp[311] = (-vv[12]-vv[46]+vv[48]+vv[58]-2*(tyv[48]+tyv[58]));
  kp[312] = (-vv[6]+vv[18]+vv[24]+vv[36]+tyv[6]-2*txv[36]-tyv[36]);
  kp[313] = (+vv[2]-vv[8]+vv[22]-vv[32]-vv[42]+vv[62]-2*(txv[2]+txv[22]-txv[32]+txv[62])-tyv[2]+tyv[32]-2*tyv[62]);
  kp[314] = (+vv[4]-vv[16]+vv[26]+vv[38]-2*(txv[38]+sxv1));
  kp[315] = (+vv[30]-vv[54]+2*tyv[54]);
  kp[316] = (+vv[28]+vv[52]-2*(txv[52]+tyv[52]));
  kp[317] = (-vv[12]-vv[48]-vv[46]-vv[58]+2*(txv[48]+tyv[48]+txv[58]+tyv[58]));
  kp[318] = (+vv[30]+vv[54]-2*(txv[54]+tyv[54]));
  kp[319] = (+vv[4]-vv[16]-vv[26]+vv[38]);
  kp[320] = (+vv[44]-vv[56]);
  kp[321] = (+vv[45]-vv[57]);
  kp[322] = (+vv[12]+vv[48]-vv[29]+vv[53]-2*(txv[12]+tyv[12]+sxyv3));
  kp[323] = (-vv[8]+vv[25]-vv[32]-vv[37]+2*(txv[8]+tyv[8]+sxyv2));
  kp[324] = (-vv[45]+vv[57]);
  kp[325] = (+vv[44]-vv[56]);
  kp[326] = (+vv[13]+vv[28]-vv[49]+vv[52]-tyv[13]-2*txv[28]-tyv[28]);
  kp[327] = (-vv[9]-vv[24]+vv[33]-vv[36]+tyv[9]+2*txv[24]+tyv[24]);
  kp[328] = (-vv[12]+vv[48]+vv[29]-vv[53]+2*(tyv[12]-sxv3+sxyv3));
  kp[329] = (-vv[13]-vv[28]+vv[49]+vv[52]+tyv[13]+tyv[28]);
  kp[330] = (+vv[44]+vv[56]-2*txv[44]);
  kp[331] = (+vv[0]+vv[5]-vv[17]+vv[20]-vv[40]+vv[60]-2*(txv[0]+txv[20]+txv[60]+sxyv0)-2*tyv[0]-tyv[5]-tyv[20]);
  kp[332] = (+vv[8]-vv[25]-vv[32]+vv[37]-2*(tyv[8]-sxv2+sxyv2));
  kp[333] = (+vv[9]+vv[24]-vv[33]-vv[36]-tyv[9]-tyv[24]);
  kp[334] = (-vv[0]-vv[5]+vv[17]-vv[20]-vv[40]+vv[60]+2*(txv[0]+txv[20]+txv[40]+sxyv0)+2*tyv[0]+tyv[5]+tyv[20]);
  kp[335] = (-vv[44]-vv[56]+2*txv[56]);
  kp[336] = (+vv[45]-vv[57]);
  kp[337] = (-vv[44]+vv[56]);
  kp[338] = (+vv[13]+vv[28]+vv[49]-vv[52]-2*txv[13]-tyv[13]-tyv[28]);
  kp[339] = (-vv[9]-vv[24]-vv[33]+vv[36]+2*txv[9]+tyv[9]+tyv[24]);
  kp[340] = (+vv[44]-vv[56]);
  kp[341] = (+vv[45]-vv[57]);
  kp[342] = (-vv[12]+vv[29]+vv[48]+vv[53]-2*(txv[29]+tyv[29]+sxyv3));
  kp[343] = (+vv[8]-vv[25]-vv[32]-vv[37]+2*(txv[25]+tyv[25]+sxyv2));
  kp[344] = (-vv[13]-vv[28]+vv[49]+vv[52]+tyv[13]+tyv[28]);
  kp[345] = (+vv[12]-vv[29]-vv[48]+vv[53]+2*(tyv[29]-sxv3+sxyv3));
  kp[346] = (+vv[45]+vv[57]-2*(txv[45]));
  kp[347] = (+vv[1]-vv[4]+vv[16]+vv[21]-vv[41]+vv[61]-2*(txv[1]+txv[21]+txv[61]+sxyv1)-tyv[1]-tyv[16]-2*tyv[21]);
  kp[348] = (+vv[9]+vv[24]-vv[33]-vv[36]-tyv[9]-tyv[24]);
  kp[349] = (-vv[8]+vv[25]+vv[32]-vv[37]-2*(tyv[25]-sxv2+sxyv2));
  kp[350] = (-vv[1]+vv[4]-vv[16]-vv[21]-vv[41]+vv[61]+2*(txv[1]+txv[21]+txv[41]+sxyv1)+tyv[1]+tyv[16]+2*tyv[21]);
  kp[351] = (-vv[45]-vv[57]+2*txv[57]);
  kp[352] = (-vv[12]-vv[29]-vv[48]+vv[53]+2*(txv[12]+sxv3));
  kp[353] = (-vv[13]+vv[28]-vv[49]-vv[52]+2*txv[13]+tyv[13]-tyv[28]);
  kp[354] = (+vv[44]-vv[56]-2*tyv[44]);
  kp[355] = (+vv[0]+vv[5]+vv[17]-vv[20]-vv[40]-vv[60]+2*(-txv[5]+txv[20]+txv[40]+txv[60])-tyv[5]+tyv[20]+2*tyv[40]);
  kp[356] = (+vv[13]-vv[28]-vv[49]-vv[52]+2*txv[28]-tyv[13]+tyv[28]);
  kp[357] = (-vv[12]-vv[29]+vv[48]-vv[53]+2*(txv[29]+sxv3));
  kp[358] = (+vv[45]-vv[57]-2*tyv[45]);
  kp[359] = (-vv[1]+vv[4]+vv[16]+vv[21]-vv[41]-vv[61]+2*(txv[1]-txv[16]+txv[41]+txv[61])+tyv[1]-tyv[16]+2*tyv[41]);
  kp[360] = (-vv[44]-vv[56]+2*(txv[44]+tyv[44]));
  kp[361] = (-vv[45]-vv[57]+2*(txv[45]+tyv[45]));
  kp[362] = (-vv[12]-vv[29]+vv[48]+vv[53]);
  kp[363] = (+vv[8]+vv[25]+vv[32]+vv[37]-2*(txv[32]+tyv[32]+txv[37]+tyv[37]));
  kp[364] = (+vv[0]-vv[5]+vv[17]+vv[20]+vv[40]-vv[60]-2*(txv[0]+txv[20]+txv[40]+sxv0)+tyv[5]-tyv[20]-2*tyv[40]);
  kp[365] = (+vv[1]+vv[4]-vv[16]+vv[21]+vv[41]-vv[61]-2*(txv[1]+txv[21]+txv[41]+sxv1)-tyv[1]+tyv[16]-2*tyv[41]);
  kp[366] = (+vv[8]+vv[25]-vv[32]-vv[37]+2*(tyv[32]+tyv[37]));
  kp[367] = (+vv[12]+vv[29]+vv[48]+vv[53]-2*(txv[48]+txv[53]));
  kp[368] = (+vv[8]+vv[25]+vv[32]-vv[37]-2*(txv[8]+sxv2));
  kp[369] = (+vv[9]-vv[24]+vv[33]+vv[36]-2*txv[9]-tyv[9]+tyv[24]);
  kp[370] = (-vv[0]-vv[5]-vv[17]+vv[20]+vv[40]+vv[60]-2*(-txv[5]+txv[20]+txv[40]+txv[60])+tyv[5]-tyv[20]-2*tyv[60]);
  kp[371] = (+vv[44]-vv[56]+2*(tyv[56]));
  kp[372] = (-vv[9]+vv[24]+vv[33]+vv[36]+tyv[9]-2*txv[24]-tyv[24]);
  kp[373] = (+vv[8]+vv[25]-vv[32]+vv[37]-2*(txv[25]+sxv2));
  kp[374] = (+vv[1]-vv[4]-vv[16]-vv[21]+vv[41]+vv[61]-2*(txv[1]-txv[16]+txv[41]+txv[61])-tyv[1]+tyv[16]-2*tyv[61]);
  kp[375] = (+vv[45]-vv[57]+2*tyv[57]);
  kp[376] = (-vv[0]+vv[5]-vv[17]-vv[20]+vv[40]-vv[60]
	       +2*(txv[0]+txv[20]+txv[60]+sxv0)-tyv[5]+tyv[20]+2*tyv[60]) ;
  kp[377] = (-vv[1]-vv[4]+vv[16]-vv[21]+vv[41]-vv[61]
	       +2*(txv[1]+txv[21]+txv[61]+sxv1)+tyv[1]-tyv[16]+2*tyv[61]);
  kp[378] = (-vv[8]-vv[25]-vv[32]-vv[37]+2*(txv[32]+txv[37]));
  kp[379] = (-vv[12]-vv[29]+vv[48]+vv[53]-2*(tyv[48]+tyv[53]));
  kp[380] = (+vv[44]+vv[56]-2*(txv[56]+tyv[56]));
  kp[381] = (+vv[45]+vv[57]-2*(txv[57]+tyv[57]));
  kp[382] = (-vv[12]-vv[29]-vv[48]-vv[53]+2*(txv[48]+tyv[48]+txv[53]+tyv[53]));
  kp[383] = (+vv[8]+vv[25]-vv[32]-vv[37]);
  return ;
}

// returns 1 if something bad happened
static int
init_STV( const double xv[4] ,
	  const double yv[4] ,
	  struct QED_kernel_temps t ,
	  struct STV *k )
{
  const struct invariants Inv = set_invariants( xv , yv , t.Grid ) ;

  memset( k -> Sxv , 0 , 4*sizeof( double ) ) ;
  memset( k -> Syv , 0 , 4*sizeof( double ) ) ;
  memset( k -> Txv , 0 , 64*sizeof( double ) ) ;
  memset( k -> Tyv , 0 , 64*sizeof( double ) ) ;
  memset( k -> Vv  , 0 , 64*sizeof( double ) ) ;

  // precompuations for the interpolations
  const size_t ix1 = Inv.INVx.idx ;
  size_t ix2 = ix1+1 ;
  if( ix2 >= (size_t)t.Grid.nstpx ) {
    ix2 = ix1 ;
  }

#if (defined HAVE_IMMINTRIN_H) && (defined __AVX__)

  const size_t thread = (size_t)omp_get_thread_num() ;
  const size_t toff = t.Grid.Nffa*thread ;
  
  const size_t nx1 = (size_t)t.Grid.nfx[ix1];
  const size_t nx2 = (size_t)t.Grid.nfx[ix2];
  const size_t nmin = (size_t)((nx1<nx2)?nx1:nx2);

  // y edge case
  const size_t iy1 = Inv.INVy.idx ;
  size_t iy2 = iy1 + 1 ;
  if( iy2 >= (size_t)t.Grid.nstpy ) {
    iy2 = iy1 ;
  }
 
  size_t i ;
  for( i = 0 ; i < (size_t)t.Grid.Nffa ; i++ ) {
    
    t.Grid.PC[i+toff].nmax = (nx1<nx2)?nx2:nx1 ;
    
    const float *Fm1 = t.Grid.Ffm[i][ix1][iy1] ;
    const float *Fm2 = t.Grid.Ffm[i][ix1][iy2] ;
    const float *Fm3 = t.Grid.Ffm[i][ix2][iy1] ;
    const float *Fm4 = t.Grid.Ffm[i][ix2][iy2] ;
    const float *Fp1 = t.Grid.Ffp[i][ix1][iy1] ;
    const float *Fp2 = t.Grid.Ffp[i][ix1][iy2] ;
    const float *Fp3 = t.Grid.Ffp[i][ix2][iy1] ;
    const float *Fp4 = t.Grid.Ffp[i][ix2][iy2] ;

    __m256d *pFm = (__m256d*)t.Grid.PC[i+toff].Fm ;
    __m256d *pFp = (__m256d*)t.Grid.PC[i+toff].Fp ;
    
    size_t j ;
    for( j = 0 ; j < nmin ; j++ ) {
      pFm[j] = _mm256_setr_pd( *Fm1 , *Fm2 , *Fm3 , *Fm4 ) ;
      pFp[j] = _mm256_setr_pd( *Fp1 , *Fp2 , *Fp3 , *Fp4 ) ;
      Fm1++ ; Fm2++ ; Fm3++ ; Fm4++ ;
      Fp1++ ; Fp2++ ; Fp3++ ; Fp4++ ;	   
    }
    if( nx1 < nx2) {
      for( ; j < nx2 ; j++ ) {
	pFm[j] = _mm256_setr_pd( 0 , 0 , *Fm3 , *Fm4 ) ;
	pFp[j] = _mm256_setr_pd( 0 , 0 , *Fp3 , *Fp4 ) ;
	Fm3++ ; Fm4++ ;
	Fp3++ ; Fp4++ ;    
      }
    } else {
      for( ; j < nx1 ; j++ ) {
	pFm[j] = _mm256_setr_pd( *Fm1 , *Fm2 , 0 , 0 ) ;
	pFp[j] = _mm256_setr_pd( *Fp1 , *Fp2 , 0 , 0 ) ;
	Fm1++ ; Fm2++ ;
	Fp1++ ; Fp2++ ;	
      }
    }
  }

  if( chnr_dS( xv, yv, Inv, t.Grid, t.Grid.PC+toff, k->Sxv, k->Syv ) ||
      chnr_dT( xv, yv, Inv, t.Grid, t.Grid.PC+toff, k->Txv, k->Tyv ) ||
      chnr_dV( xv, yv, Inv, t.Grid, t.Grid.PC+toff, k->Vv ) ) {
    return 1 ;
  }
  
#else

  if( chnr_dS( xv, yv, Inv, t.Grid, t.Grid.PC, k->Sxv, k->Syv ) ||
      chnr_dT( xv, yv, Inv, t.Grid, t.Grid.PC, k->Txv, k->Tyv ) ||
      chnr_dV( xv, yv, Inv, t.Grid, t.Grid.PC, k->Vv ) ) {
    return 1 ;
  }
#endif

  return 0 ;
}

// computes the kernel [rho - sigma];[mu][nu][lambda]
int
kernelQED( const double xv[4] ,
	   const double yv[4] ,
	   const struct QED_kernel_temps t ,
	   double kerv[6][4][4][4] )
{
  // FFs of QED kernel and their derivatives wrt x,cb,y
  struct STV x_y ;
  if( init_STV( xv, yv, t, &x_y ) ) return 1 ;

  // point out the kernel again
  double *kp = (double*)kerv ;
  CONSTRUCT_FULL_KERNEL( kp , x_y ) ;

  return 0 ;
}
