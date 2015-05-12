#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#include <math.h>
#define Q 		19
#define D 		3
#define SQRT2 	1.41421356237309504880
#define SQRT3	1.73205080756887729353

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
	
		SQRT2,
		SQRT2,
		1,
		SQRT2,
		SQRT2,
		
		SQRT2,
		1,
		SQRT2,
		1,
		0,
		
		1,
		SQRT2,
		1,
		SQRT2,
		SQRT2,
		
		SQRT2,
		1,
		SQRT2,
		SQRT2		
		
	};
	
	static const double C_S = 1.0 / SQRT3;

#endif

