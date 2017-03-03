///////////////////////////////////////////////////////////////////
// s c a n s . h
// scanning procedures - boundary mapping, etc.                  
//                                                               
// started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          
///////////////////////////////////////////////////////////////////

#include <iostream.h>
/////////////////////////////////////////////////////////
//  s c a n s . h
//
//  include file for scan procedures
/////////////////////////////////////////////////////////

#define X	0
#define Y	1
#define Z	2
#define OMEGA	0
#define CHI		0
#define PHI		0
#define PI		(4*atan(1))


extern	bool toggle[];
extern	double value[];
extern	int voxel_disp;
extern	int boundaryfileext;
extern	char *boundaryfile;
extern	char *detectordatafile;
extern	char *logfile;
extern  ExpWindow *expwin;
extern  DetWindow *detwin;


void autorotate(void);
void find_peak(void);
int  get_bragg_intensity(void);
void step(double,double,bool);
void rotstep(double, double&, double&);
void boundary_mark(ofstream);
bool boundary_search(ofstream);
void map_grain(void);
void map_crystal(void);
