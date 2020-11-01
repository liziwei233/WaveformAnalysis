#include "global.h"

bool gBinaryFlag=1;
int gNsample = 1024;
float ghorizontal_interval = 0.2; // ns
//float gvertical_gain_ch = 1;
//float gvertical_gain_tr = 1;
float gvertical_gain_ch[32] = {0.2586,0.2590,0.2608,0.2588,0.2579,0.2612,0.2610,0.2618,0.2584,0.2582,0.2587,0.2598,0.2590,0.2607,0.2624,0.2601};
float gvertical_gain_tr[2] = {0.5567, 0.5532};
float vertical_offset = 0;
float horizontal_offset = 0;
bool filecreate=0;