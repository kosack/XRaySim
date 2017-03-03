/*****************************************************************
 * i n t e r f a c e . h                                         *
 * graphical user interface using Fltk                           *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/
#ifndef __INTERFACE_H__
#define __INTERFACE_H__

#include <iostream.h>
#include <stdio.h>
#include <FL/Fl.H>
#include <FL/gl.h>
#include <GL/glu.h>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Hor_Slider.H>
#include <FL/Fl_Roller.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Menu_Item.H> 
#include <FL/Fl_Choice.H> 
#include <FL/Fl_Menu_Bar.H> 
#include <FL/fl_file_chooser.H> 
#include <FL/Fl_Ask.H> 
#include <FL/fl_draw.H> 
#include <FL/math.h>
#include "global.h"
#include "libsvl/SVL.h"
#include "xraysim.h"


/////////////////////////////////////////////////////////////////////
// CONSTANTS AND ENUMS:

const int NUMTOGGLES = 8;
enum toggles { 
	T_UPDATE_DETECTOR, T_UPDATE_VIEW, T_BEAM_ON,
	T_SAMPLE_REF_FRAME, T_OMEGA, T_CHI, T_PHI, T_LOG
};	


const int NUMVALUES = 29;
enum values {
	V_OMEGA, V_CHI, V_PHI, 
	V_OMEGASTEP, V_CHISTEP, V_PHISTEP,
	V_OMEGASTOP, V_CHISTOP, V_PHISTOP,
	V_TX, V_TY, V_TZ,
	V_LOOKAT_X, V_LOOKAT_Y, V_LOOKAT_Z,
	V_LOOKFROM_X, V_LOOKFROM_Y, V_LOOKFROM_Z,
	V_UP_X, V_UP_Y, V_UP_Z, 
	V_DETDIST,
	V_R, V_THETA, V_PSI,

	V_BOUNDARY_TOLER,
	V_BOUNDARY_TURN,
	V_BOUNDARY_STEPSIZE

};


enum window_layout {
	W_WIDTH = 800,
	W_HEIGHT = 600,
	W_XGRID = 8,
	W_YGRID = 8
};



enum commands {
	C_EXIT, C_AUTOROTATE, C_FIND_PEAK, C_FIREONE, C_MAP_GRAIN,
	C_MAP_CRYSTAL
};


enum menuoptions { 
	M_BEAMON, M_BEAMOFF, M_AUTOOMEGA, M_AUTOCHI, M_AUTOPHI, M_EXIT,
	M_SAVEDET, M_RESETDET, M_NULL, M_STATUS, M_RF_SAMP, M_RF_DET,
	M_DET_UPDATE_AUTO, M_DET_UPDATE_MAN, M_DIS_UPDATE_AUTO, M_DIS_UPDATE_MAN
};


enum VoxelDisplayTypes {
	VOXEL_NONE,
	VOXEL_POINTS,
	VOXEL_BOUNDARY_POINTS
};

#define POS true;
#define NEG false;


class XRayScope;
extern XRayScope scope;

///////////////////////////////////////////////////////////////////////
// Prototypes:

void interface_init(void);
int interface_start(int argc, char **argv);		// start the user interface
void callbackConstrainSel(Fl_Widget *o, void*);


///////////////////////////////////////////////////////////////////////
// CLASSES

class ExpWindow : public Fl_Gl_Window {
	// Experiment view window

  public:
	ExpWindow(int x, int y, int w, int h, const char *l, XRayScope *s);
	void draw(void);
	void init(void);

  private:
	XRayScope *scope;

};


class DetWindow : public Fl_Gl_Window {
	// Detector View window

  public:
	DetWindow(int x, int y, int w, int h, const char *l, XRayScope *s);
	void draw();
	void draw_overlay(void);
	void init();
	void eraseBox();
	int handle(int);

	int sx, sy, sw, sh, ix, iy;

  private:
	XRayScope *scope;

};


#endif //__INTERFACE_H__