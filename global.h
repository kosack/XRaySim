/*****************************************************************
 * g  l  o  b  a  l  .  h                                        *
 * Constants and global vars                                     *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/
#ifndef __GLOBAL_H__
#define __GLOBAL_H__


/////////////////////////////////////////////////////////////////////////
// CONSTANTS:


#define		VERSION			0.99				// Version 


extern bool DEBUG;
extern bool GRAPHICS;

extern double PI;							// PI constant defined in xraysim.cpp

typedef unsigned char Byte;
typedef int Voxel;
#define		EMPTY_VOXEL		-1

#ifndef USE_OPENGL
#define USE_OPENGL
#endif

#define		debug			if(DEBUG) cout << "DEBUG: "	
#define		MATRIX_MODES	2				// number of matrix modes
#define		MAXGRAINS		255				// max number of grains 
#define		TOLER			1e-5			// intersect tolerance
#define		INFINITY		1.0e8			// really big number
#define		MAX_MILLER_INDEX 5				// miller indices to compute G's for

#define		RADTODEG		180.0/PI		// angle conversion
#define		DEGTORAD		PI/180.0		// angle conversion


enum matrix_modes {							// list of matrix modes
	MODEL_MATRIX_MODE=0,
	VIEW_MATRIX_MODE=1
};

enum axes {
	X = 0,
	Y = 1,
	Z = 2,
	W = 3
};

enum GrainType { SC, FCC, BCC };




#endif //__GLOBAL_H__