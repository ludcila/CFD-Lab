#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#include <math.h>
#define Q 			19
#define D 			3
#define W_12_36		12.0/36.0
#define W_2_36		2.0/36.0
#define W_1_36		1.0/36.0
#define SQRT3		1.73205080756887729353

	static const int LATTICEVELOCITIES[19][3] = {
	
		{ 0, -1, -1},
		{-1,  0, -1},
		{ 0,  0, -1},
		{ 1,  0, -1},
		{ 0,  1, -1},
		
		{-1, -1,  0},
		{ 0, -1,  0},
		{ 1, -1,  0},
		{-1,  0,  0},
		{ 0,  0,  0},
		
		{ 1,  0,  0},
		{-1,  1,  0},
		{ 0,  1,  0},
		{ 1,  1,  0},
		{ 0, -1,  1},
		
		{-1,  0,  1},
		{ 0,  0,  1},
		{ 1,  0,  1},
		{ 0,  1,  1}
		
	};
	
	static const double LATTICEWEIGHTS[19] = {
	
		W_1_36,
		W_1_36,
		W_2_36,
		W_1_36,
		W_1_36,
		
		W_1_36,
		W_2_36,
		W_1_36,
		W_2_36,
		W_12_36,
		
		W_2_36,
		W_1_36,
		W_2_36,
		W_1_36,
		W_1_36,
		
		W_1_36,
		W_2_36,
		W_1_36,
		W_1_36		
		
	};
	
	static const double C_S = 1.0 / SQRT3;

#endif

