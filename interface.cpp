/*****************************************************************
 * i n t e r f a c e . c p p                                     *
 * graphical user interface using Fltk                           *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/

#include <iostream.h>
#include "interface.h"
#include "scans.h"


/////////////////////////////////////////////////////////////////////
// DATA FILES AND NAMES
char *boundaryfile;
char *detectordatafile;
char *logfile;

/////////////////////////////////////////////////////////////////////
// GLOBALS:
bool toggle[NUMTOGGLES];
double value[NUMVALUES];
int voxel_disp = VOXEL_NONE;
int boundaryfileext = 0;

ExpWindow *expwin;
DetWindow *detwin;


ExpWindow::ExpWindow(int x,int y,int w,int h,const char *l, XRayScope *s) : Fl_Gl_Window(x,y,w,h,l) {

	scope = s;

}




void
ExpWindow::draw() {
	// Display main viewing pane

	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);

	// SET THE VIEWING POSITION:
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();	// modify the default projection matrix (pop when done displaying)
	gluLookAt(value[V_LOOKFROM_X], value[V_LOOKFROM_Y], value[V_LOOKFROM_Z],
			  value[V_LOOKAT_X], value[V_LOOKAT_Y], value[V_LOOKAT_Z],
			  value[V_UP_X], value[V_UP_Y], value[V_UP_Z]);

	// IF WE'RE NOT IN THE SAMPLE REFERENCE FRAME, THEN TRANSFORM
	// PROJECTION INTO DETECTOR'S FRAME (inverse Euler rotation)
	if (!toggle[T_SAMPLE_REF_FRAME]){
		glTranslated(value[V_TX], value[V_TY], value[V_TZ]); 
		glRotated(-value[V_PHI]  *RADTODEG, 0,0,1);
		glRotated(-value[V_CHI]*RADTODEG, 1,0,0);
		glRotated(-value[V_OMEGA]  *RADTODEG, 0,0,1);
	}
	glMatrixMode(GL_MODELVIEW);

	scope->rotateSample(value[V_OMEGA], value[V_CHI], value[V_PHI]);
	scope->translateSample(value[V_TX], value[V_TY], value[V_TZ]);
	scope->display();
	if (toggle[T_BEAM_ON]){
		scope->activate();		
	}

	// Draw some axes:
	glBegin(GL_LINES);
		glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(10,0,0);
		glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,10,0);
		glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,10);
	glEnd();



	glMatrixMode(GL_PROJECTION);
	glPopMatrix();	// restore the default projection matrix
	glMatrixMode(GL_MODELVIEW);

}


void
ExpWindow::init(void){

	make_current();
	mode(FL_DOUBLE|FL_RGB);

	glClearColor(0,0,0.4,1);

	// SET UP OpenGL RENDERING OPTIONS:	
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glShadeModel(GL_FLAT);
	glEnable(GL_DEPTH_TEST); // Z-Buffer


	// SETUP OpenGL LIGHTING:
	GLfloat white8[] = {0.8, 0.8, 0.8, 1.};
	GLfloat white2[] = {0.2, 0.2, 0.2, 1.};
	GLfloat black[] = {0., 0., 0., 1.};
	GLfloat mat_shininess[] = {5064.};              /* Phong exponent */
	GLfloat light0_position[] = {1., 1., 5., 0.}; /* directional light (w=0) */
	GLfloat white[] = {1., 1., 1., 1.};
	GLfloat light1_position[] = {-3., 1., -1., 0.};
	GLfloat red[] = {1., .3, .3, 1.};
	GLfloat light2_position[] = {1., 8., -2., 0.};
	GLfloat blue[] = {.3, .4, 1., 1.};

	glMatrixMode(GL_MODELVIEW);	  
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, white2);   /* no ambient */
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, white8);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	  
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
	glLightfv(GL_LIGHT0, GL_SPECULAR, white);
	glEnable(GL_LIGHT0);
		  
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, white);
	glLightfv(GL_LIGHT1, GL_SPECULAR, white);
	glEnable(GL_LIGHT1);
	  
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, white);
	glLightfv(GL_LIGHT2, GL_SPECULAR, white);
	glEnable(GL_LIGHT2);
	  
	//glEnable(GL_NORMALIZE);       /* normalize normal vectors */
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);      /* two-sided lighting*/
	glEnable(GL_LIGHTING);

	// SET UP DEFAULT PROJECTION MATRIX:
 	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w(),h());
	gluPerspective(50, 1.25, 0.01, 50.0);
	glMatrixMode(GL_MODELVIEW);

}

////////////////////////////////////////////////////////////////////////////////
// DETECTOR WINDOW METHODS

DetWindow::DetWindow(int x,int y,int w,int h,const char *l, XRayScope *s) : Fl_Gl_Window(x,y,w,h,l) {
	scope = s;
	sx=sy=ix=iy=0;
	sw=1; sh=1;
}

void
DetWindow::draw(void){

	make_current();
	glClearColor(0.4,0.4,0.4,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	scope->detector.image->display();
	
}

void
DetWindow::draw_overlay(void) {
	
	make_overlay_current();
	ortho();
	gl_color(FL_RED);
	glBegin(GL_LINE_LOOP);
		glVertex2d(sx,sy);
		glVertex2d(sx,sy+sh);
		glVertex2d(sx+sw,sy+sh);
		glVertex2d(sx+sw,sy);
	glEnd();

}


void
DetWindow::eraseBox(void){
	window()->make_current();
	fl_overlay_clear();
}

int 
DetWindow::handle(int event) {
	static int dragged;
	static int button;
	int x2,y2;
	switch (event) {
		case FL_PUSH:
			ix = Fl::event_x(); if (ix<0) ix=0; if (ix>=w()) ix=w()-1;
			iy = h()-Fl::event_y(); if (iy<0) iy=0; if (iy>=h()) iy=h()-1;
			if(ix >=scope->detector.getXRes()) ix = scope->detector.getXRes()-1;
			if(iy >=scope->detector.getYRes()) iy = scope->detector.getYRes()-1;
			dragged = 0;
			button = Fl::event_button();
			redraw_overlay();
			return 1;
		case FL_DRAG:
			dragged = 1;
			x2 = Fl::event_x(); if (x2<0) x2=0; if (x2>=w()) x2=w()-1;
			y2 = h()-Fl::event_y(); if (y2<0) y2=0; if (y2>=h()) y2=h()-1;
			if(x2 >=scope->detector.getXRes()) x2 = scope->detector.getXRes()-1;
			if(y2 >=scope->detector.getYRes()) y2 = scope->detector.getYRes()-1;
			if (button != 1) {ix = x2; iy = y2; return 1;}
			if (ix < x2) {sx = ix; sw = x2-ix;} else {sx = x2; sw = ix-x2;}
			if (iy < y2) {sy = iy; sh = y2-iy;} else {sy = y2; sh = iy-y2;}
			redraw_overlay();
			return 1;
		case FL_RELEASE:
			callbackConstrainSel(NULL, NULL);
			break;

	}
	return 0;

}


void
DetWindow::init(void){

	make_current();
	glLoadIdentity();
	ortho(); 

}





////////////////////////////////////////////////////////////////////////////////
// CALLBACKS:


void
callbackAbout(Fl_Widget *o, void *){

	fl_message("XRaySim - by Karl Kosack <kosack@andrew.cmu.edu> 1998");

}

void
callbackToggle(Fl_Widget *o, long arg){

	toggle[(int)arg] = (bool)(((Fl_Check_Button*)o)->value());

	if((int)arg==T_BEAM_ON) expwin->redraw();

}

void
callbackWindow(Fl_Widget *o, long arg){
	cout << "Bye!" << endl;
	exit(1);
}


void
callbackCommand(Fl_Widget *o, long arg){

	switch ((int)arg){

		case C_EXIT: exit(0); break;

		case C_AUTOROTATE: autorotate(); break;

		case C_FIND_PEAK: find_peak(); break;

		case C_FIREONE: 
			scope.rotateSample(value[V_OMEGA], value[V_CHI], value[V_PHI]);
			scope.translateSample(value[V_TX], value[V_TY], value[V_TZ]);
			scope.activate();
			expwin->redraw();
			detwin->redraw();
			break;

		case C_MAP_GRAIN:
			map_grain();
			break;

		case C_MAP_CRYSTAL:
			map_crystal();
			break;

	}


}


void callbackViewPos(Fl_Widget *o, long arg){


	switch(arg){
		case V_R:
			value[(int)arg] = (double)(((Fl_Value_Input*)o)->value());
		case V_THETA:
			value[(int)arg] = (double)(((Fl_Value_Input*)o)->value());
		case V_PSI:
			value[(int)arg] = (double)(((Fl_Value_Input*)o)->value());
	}

	value[V_LOOKFROM_X] = value[V_R] * sin(value[V_THETA]) * sin(value[V_PSI]);	
	value[V_LOOKFROM_Y] = value[V_R] * sin(value[V_THETA]) * cos(value[V_PSI]);	
	value[V_LOOKFROM_Z] = value[V_R] * cos(value[V_THETA]);
	expwin->redraw();

}

void
callbackValue(Fl_Widget *o, long arg){

	value[(int)arg] = (double)(((Fl_Value_Input*)o)->value());
	//if (toggle[T_UPDATE_VIEW]) 
	expwin->redraw();

}


void
callbackExit(Fl_Widget *o, void*){

	exit(0);

}

void
callbackClear(Fl_Widget *o, void*){
	scope.detector.image->clear();
	detwin->redraw();
}


void
callbackLoad(Fl_Widget *o, void*){
	char *filename;

	filename = fl_file_chooser( "Choose sample datafile", "*.xrs", "sample.xrs" );
	if(filename == NULL) return;
	else scope.loadSample(filename);

}

void
callbackSave(Fl_Widget *o, void*){
	scope.writeLogHeader("Long detector exposure - ignore orientation info");
	scope.writeLogEntry();
	fl_message( "Detector output was saved to disk.");
}

void
callbackShowDet(Fl_Widget *o, void*){
	//detwin->init();
	detwin->show();
	detwin->redraw();
}


void
callbackOutput(Fl_Widget *o, long arg){

	switch((int)((Fl_Choice*)o)->value()){
		case 0:
			toggle[T_LOG] = false;
			break;

		case 1:
			toggle[T_LOG] = true;
			break;

	}


}

void
callbackVoxelDisp(Fl_Widget *o, long arg){

	switch((int)((Fl_Choice*)o)->value()){
		case VOXEL_NONE:
			scope.setVoxelDisp(VOXEL_NONE);
			break;
		case VOXEL_POINTS:
			scope.setVoxelDisp(VOXEL_POINTS);
			break;
		case VOXEL_BOUNDARY_POINTS:
			scope.setVoxelDisp(VOXEL_BOUNDARY_POINTS);
			break;
	}

	expwin->redraw();

}

void
callbackBeamDiv(Fl_Widget *o, long arg){
	scope.setBeamDiv((double)(((Fl_Value_Input*)o)->value()));
}

void
callbackDetDist(Fl_Widget *o, long arg){
	scope.setDetectorDist((double)(((Fl_Value_Input*)o)->value()));
	expwin->redraw();
}


void
callbackGenerate(Fl_Widget *o, void*){
	scope.generateSample();
	expwin->redraw();
}

void
callbackSampVolMin(Fl_Widget *o, long arg){
	double min=scope.getProbeVolumeMin(), max=scope.getProbeVolumeMax();
	scope.setProbeVolumeSize((double)(((Fl_Value_Input*)o)->value()), max);
}


void
callbackSampVolMax(Fl_Widget *o, long arg){
	double min=scope.getProbeVolumeMin(), max=scope.getProbeVolumeMax();
	scope.setProbeVolumeSize(min, (double)(((Fl_Value_Input*)o)->value()));
}


void 
callbackConstrainSel(Fl_Widget *o, void*){
	int fullval, val=0;
	int tx0 = detwin->sx;
	int ty0 = detwin->sy;
	int tx1 = detwin->sw + detwin->sx;
	int ty1 = detwin->sh + detwin->sy;
	int vx0 = 1, vy0=1, vx1=-1, vy1=-1;

	fullval = scope.watchDetectorSpot( tx0, ty0, tx1, ty1);
	val = fullval;

	while(1){


		tx0+=vx0; if (tx0>=tx1) {tx0=tx1-vx0; vx0=0;}
		val=scope.watchDetectorSpot(tx0,ty0,tx1,ty1);
		if (val != fullval) { tx0 -= 1; vx0=0; }

		ty0+=vy0; if (ty0>=ty1) {ty0=ty1-vy0; vy0=0;}
		val=scope.watchDetectorSpot(tx0,ty0,tx1,ty1);
		if (val != fullval) { ty0 -= 1; vy0=0; }

		tx1+=vx1; if (tx1<=tx0) {tx1=tx0-vx1; vx1=0;}
		val=scope.watchDetectorSpot(tx0,ty0,tx1,ty1);
		if (val != fullval) { tx1 -= -1; vx1=0; }

		ty1+=vy1; if (ty1<=ty0) {ty1=ty0-vy1; vy1=0;}
		val=scope.watchDetectorSpot(tx0,ty0,tx1,ty1);
		if (val != fullval) { ty1 -= -1; vy1=0; }

		if(abs(tx1-tx0) <= 1 && abs(ty1-ty0) <=1) break;
		if(vx0==0 && vx1==0 && vy0==0 && vy1==0) break;

		detwin->sx = tx0;
		detwin->sy = ty0;
		detwin->sw = tx1-tx0;
		detwin->sh = ty1-ty0;
		detwin->draw_overlay();

	}

	if (tx0 >= 2) tx0 -=2;
	if (ty0 >= 2) ty0 -=2;
	if (tx1 <= scope.detector.getXRes()-2) tx1 +=2;
	if (ty1 <= scope.detector.getYRes()-2) ty1 +=2;


	detwin->sx = tx0;
	detwin->sy = ty0;
	detwin->sw = tx1-tx0;
	detwin->sh = ty1-ty0;
	detwin->redraw_overlay();


}






////////////////////////////////////////////////////////////////////////////////
// INTERFACE FUNCTIONS

void
interface_init(void){
	// SET DEFAULTS:

	int i;

	for(i=0; i<NUMVALUES; i++){
		value[i] = 0;
	}

	for(i=0; i<NUMTOGGLES; i++){
		toggle[i] = false;
	}

	value[V_UP_Z] = 1.0;
	value[V_R] = 15;
	value[V_THETA] = 1.4;
	value[V_PSI] = -1.3;
	value[V_LOOKFROM_X] = -14.243;
	value[V_LOOKFROM_Y] = 3.954;
	value[V_LOOKFROM_Z] = 2.549;

	value[V_OMEGASTOP] = 2*_PI;
	value[V_CHISTOP]   = 2*_PI;
	value[V_PHISTOP]   = 2*_PI;

	value[V_BOUNDARY_STEPSIZE] = 0.0025;
	value[V_BOUNDARY_TOLER] = 0;
	value[V_BOUNDARY_TURN] = 0.3;

}

int
interface_start(int argc, char **argv){

	int y=0,x=0;

	Fl_Menu_Item menutable[] = {
		{"&Sample",0,0,0,FL_SUBMENU},
			{"&Load Sample", 0, callbackLoad, 0},
			{"&Randomize Volume", 0, callbackGenerate, 0, FL_MENU_DIVIDER},
			{"E&xit",	FL_CTRL+'c', callbackExit, 0},
			{0},
		{"&Detector",0,0,0,FL_SUBMENU},
			{"&Reset/Clear", 0, callbackClear, 0},
			{"&Save output", 0 , callbackSave, 0,FL_MENU_DIVIDER},
			{"&Clear Selection", FL_CTRL+'a', 0, 0},
			{"Con&strain Selection", FL_CTRL+'s', callbackConstrainSel, 0},
			{0},
		{"&Window",0,0,0,FL_SUBMENU},
			{"&Show Detector",FL_CTRL+'d', callbackShowDet, 0, FL_MENU_DIVIDER},
			{"&About",0,callbackAbout,0},
			{0},
		{0},
	{0}
	};

	cout << "User interface: " << flush;

	interface_init();
	
	//===========================================================================
	// CONTROL WINDOW
	Fl_Window ctrlwin(0,30,W_WIDTH,W_HEIGHT,"XRaySim - 3D x-ray microscope simulator");
	ctrlwin.callback(callbackWindow);
	Fl_Menu_Bar menubar(0,0,W_WIDTH,W_YGRID*4, "Menubar"); menubar.menu(menutable);


	//EXPERIMENT CONTROLS
	y = W_YGRID*8;
	Fl_Group expctrls(W_XGRID,y,W_XGRID*47,W_YGRID*65,"Experiment Controls");
		expctrls.box(FL_ENGRAVED_BOX); 
		y += W_YGRID*4;
		
		//----- OMEGA
		Fl_Check_Button omegaon(W_XGRID*7, y, W_XGRID*5, W_YGRID*3, "omega");
		omegaon.align(FL_ALIGN_LEFT);
		omegaon.callback(callbackToggle, FL_WHEN_CHANGED);
		omegaon.value(toggle[T_OMEGA]);
		omegaon.argument(T_OMEGA);

		Fl_Value_Input omega(W_XGRID*11, y, W_XGRID*10, W_YGRID*3, "Start");
		omega.align(FL_ALIGN_TOP);
		omega.range(-4*_PI, 4*_PI);
		omega.step(0.001);
		omega.callback(callbackValue, FL_WHEN_CHANGED);
		omega.argument(V_OMEGA);

		Fl_Value_Input omegastep(W_XGRID*22, y, W_XGRID*10, W_YGRID*3, "Step");
		omegastep.align(FL_ALIGN_TOP);
		omegastep.range(0,1); omegastep.step(0.0001);
		omegastep.callback(callbackValue, FL_WHEN_CHANGED);
		omegastep.argument(V_OMEGASTEP);

		Fl_Value_Input omegastop(W_XGRID*33, y, W_XGRID*10, W_YGRID*3, "Stop");
		omegastop.align(FL_ALIGN_TOP);
		omegastop.step(0.0001);
		omegastop.range(-4*_PI, 4*_PI);
		omegastop.callback(callbackValue, FL_WHEN_CHANGED);
		omegastop.value(value[V_OMEGASTOP]);
		omegastop.argument(V_OMEGASTOP);


		y += W_YGRID*4;


		//----- CHI
		Fl_Check_Button chion(W_XGRID*7, y, W_XGRID*5, W_YGRID*3, "chi");
		chion.align(FL_ALIGN_LEFT);
		chion.callback(callbackToggle, FL_WHEN_CHANGED);
		chion.argument(T_CHI);

		Fl_Value_Input chi(W_XGRID*11, y, W_XGRID*10, W_YGRID*3, 0);
		chi.align(FL_ALIGN_TOP);
		chi.step(0.001);
		chi.range(-4*_PI, 4*_PI);
		chi.callback(callbackValue, FL_WHEN_CHANGED);
		chi.argument(V_CHI);

		Fl_Value_Input chistep(W_XGRID*22, y, W_XGRID*10, W_YGRID*3, 0);
		chistep.align(FL_ALIGN_TOP);
		chistep.range(0,1); chistep.step(0.0001);
		chistep.callback(callbackValue, FL_WHEN_CHANGED);
		chistep.argument(V_CHISTEP);

		Fl_Value_Input chistop(W_XGRID*33, y, W_XGRID*10, W_YGRID*3, 0);
		chistop.align(FL_ALIGN_TOP);
		chistop.range(-4*_PI, 4*_PI);
		chistop.step(0.0001);
		chistop.value(value[V_CHISTOP]);
		chistop.callback(callbackValue, FL_WHEN_CHANGED);
		chistop.argument(V_CHISTOP);
		
		y += W_YGRID*4;

		//----- PHI
		Fl_Check_Button phion(W_XGRID*7, y, W_XGRID*5, W_YGRID*3, "phi");
		phion.align(FL_ALIGN_LEFT);
		phion.callback(callbackToggle, FL_WHEN_CHANGED);
		phion.argument(T_PHI);

		Fl_Value_Input phi(W_XGRID*11, y, W_XGRID*10, W_YGRID*3, 0);
		phi.align(FL_ALIGN_TOP);
		phi.step(0.001);
		phi.range(-4*_PI, 4*_PI);
		phi.callback(callbackValue, FL_WHEN_CHANGED);
		phi.argument(V_PHI);

		Fl_Value_Input phistep(W_XGRID*22, y, W_XGRID*10, W_YGRID*3, 0);
		phistep.align(FL_ALIGN_TOP);
		phistep.range(0,1); phistep.step(0.0001);
		phistep.callback(callbackValue, FL_WHEN_CHANGED);
		phistep.argument(V_PHISTEP);

		Fl_Value_Input phistop(W_XGRID*33, y, W_XGRID*10, W_YGRID*3, 0);
		phistop.align(FL_ALIGN_TOP);
		phistop.step(0.0001);
		phistop.range(-4*_PI, 4*_PI);
		phistop.value(value[V_PHISTOP]);
		phistop.callback(callbackValue, FL_WHEN_CHANGED);
		phistop.argument(V_PHISTOP);
		y += W_YGRID*6;

		// TX
		Fl_Value_Input tx(W_XGRID*7, y, W_XGRID*10, W_YGRID*3, "tx");
		tx.align(FL_ALIGN_LEFT);
		tx.step(0.0001);
		tx.range(-5, 5);
		tx.callback(callbackValue, FL_WHEN_CHANGED);
		tx.argument(V_TX);

		// TY
		Fl_Value_Input ty(W_XGRID*21, y, W_XGRID*10, W_YGRID*3, "ty");
		ty.align(FL_ALIGN_LEFT);
		ty.step(0.0001);
		ty.range(-5, 5);
		ty.callback(callbackValue, FL_WHEN_CHANGED);
		ty.argument(V_TY);

		// TZ
		Fl_Value_Input tz(W_XGRID*35, y, W_XGRID*10, W_YGRID*3, "tz");
		tz.align(FL_ALIGN_LEFT);
		tz.step(0.0001);
		tz.range(-5, 5);
		tz.callback(callbackValue, FL_WHEN_CHANGED);
		tz.argument(V_TZ);
		y += W_YGRID*4;




		//beam
		Fl_Light_Button beam_toggle(W_XGRID*2,y,W_XGRID*10,W_YGRID*3,"Beam");
		beam_toggle.callback(callbackToggle);
		beam_toggle.argument(T_BEAM_ON);

		// Autorotate
		Fl_Button autorot(W_XGRID *15, y, W_XGRID*15, W_YGRID*3, "AutoRotate");
		autorot.callback(callbackCommand);
		autorot.color2(FL_CYAN);
		autorot.argument(C_AUTOROTATE);

		// fire one
		Fl_Button fireone(W_XGRID *30, y, W_XGRID*15, W_YGRID*3, "Fire One");
		fireone.callback(callbackCommand);
		fireone.color2(FL_RED);
		fireone.argument(C_FIREONE);
		y += W_YGRID*3;

		// map grain
		Fl_Button boundary(W_XGRID *15, y, W_XGRID*15, W_YGRID*3, "Map Grain");
		boundary.callback(callbackCommand);
		boundary.color2(FL_GREEN);
		boundary.argument(C_MAP_GRAIN);

		// find peak
		Fl_Button findpeak(W_XGRID *30, y, W_XGRID*15, W_YGRID*3, "Find Peak");
		findpeak.callback(callbackCommand);
		findpeak.color2(FL_YELLOW);
		findpeak.argument(C_FIND_PEAK);
		y += W_YGRID*4;

		// map crystal
		Fl_Button mapcrystal(W_XGRID *15, y, W_XGRID*30, W_YGRID*3, "Map Crystal");
		mapcrystal.callback(callbackCommand);
		mapcrystal.color2(FL_YELLOW);
		mapcrystal.argument(C_MAP_CRYSTAL);
		y += W_YGRID*4;

		Fl_Value_Input beamdiv(W_XGRID*22, y, W_XGRID*10, W_YGRID*3, "Beam divergence");
		beamdiv.value(scope.getBeamDiv());
		beamdiv.range(0,1); beamdiv.step(0.0001);
		beamdiv.callback(callbackBeamDiv, FL_WHEN_CHANGED);
		y += W_YGRID*4;

		Fl_Value_Input detdist(W_XGRID*22, y, W_XGRID*10, W_YGRID*3, "Detector dist");
		detdist.range(0.5,50); detdist.step(0.5);
		detdist.value(scope.getDetectorDist());
		detdist.callback(callbackDetDist, FL_WHEN_CHANGED);
		y += W_YGRID*4;

		// sample vol min
		Fl_Value_Input svolmin(W_XGRID*22, y, W_XGRID*10, W_YGRID*3, "Probe Vol min/max");
		svolmin.range(0.5,50); detdist.step(0.5);
		svolmin.value(scope.getProbeVolumeMin());
		svolmin.range(-10,0);
		svolmin.step(0.001);
		svolmin.argument(0);
		svolmin.callback(callbackSampVolMin, FL_WHEN_CHANGED);

		// sample vol max
		Fl_Value_Input svolmax(W_XGRID*32, y, W_XGRID*10, W_YGRID*3, 0);
		svolmax.range(0.5,50); detdist.step(0.5);
		svolmax.value(scope.getProbeVolumeMax());
		svolmax.range(0,10);
		svolmax.step(0.001);
		svolmax.argument(1);
		svolmax.callback(callbackSampVolMax, FL_WHEN_CHANGED);
		y += W_YGRID*4;

		// output to
		Fl_Menu_Item writetochoices[] = {{"Display",0,0,0},{"Log file",0,0,0},{0}};
		Fl_Choice writeto(W_XGRID*22,y,W_XGRID*10,W_YGRID*3,"Output to:");
		writeto.menu(writetochoices);
		writeto.callback(callbackOutput);
		y+=W_YGRID*10;

		x = W_XGRID *10;
		//boundary toler
		Fl_Value_Input boundtoler(x, y, W_XGRID*10, W_YGRID*3, "Tolerance");
		boundtoler.range(0,10); 
		boundtoler.align(FL_ALIGN_TOP);
		boundtoler.step(1);
		boundtoler.value(value[V_BOUNDARY_TOLER]);
		boundtoler.callback(callbackValue, V_BOUNDARY_TOLER);

		//boundary step
		Fl_Value_Input boundstep(x+W_XGRID*11, y, W_XGRID*10, W_YGRID*3, "Step size");
		boundstep.range(0.001,1); 
		boundstep.step(0.0001);
		boundstep.align(FL_ALIGN_TOP);
		boundstep.value(value[V_BOUNDARY_STEPSIZE]);
		boundstep.callback(callbackValue, V_BOUNDARY_STEPSIZE);

		//boundary turn
		Fl_Value_Input boundturn(x+W_XGRID*22, y, W_XGRID*10, W_YGRID*3, "Turn");
		boundturn.range(0.001,1); 
		boundturn.step(0.0001);
		boundturn.align(FL_ALIGN_TOP);
		boundturn.value(value[V_BOUNDARY_TURN]);
		boundturn.callback(callbackValue, V_BOUNDARY_TURN);



	expctrls.end();

	// DISPLAY CONTROLS
	x = W_XGRID * 50;
	y = W_YGRID * 48;
	Fl_Group dispctrls(x,y,W_XGRID*47,W_YGRID*20,"Display Controls");
		dispctrls.box(FL_ENGRAVED_BOX);
		y += W_YGRID*3;

		Fl_Check_Button update_disp(x+W_XGRID*2,y,W_XGRID*25,W_YGRID*3,"Display update");
		update_disp.callback(callbackToggle, FL_WHEN_CHANGED); 
		update_disp.color2(FL_BLUE);
		update_disp.argument(T_UPDATE_VIEW); 
		y += W_YGRID*3;

		Fl_Check_Button update_det(x+W_XGRID*2,y,W_XGRID*25,W_YGRID*3,"Detector update");
		update_det.callback(callbackToggle, FL_WHEN_CHANGED); 
		update_det.color2(FL_BLUE);
		update_det.argument(T_UPDATE_DETECTOR); 
		y += W_YGRID*3;

		Fl_Check_Button sample_ref(x+W_XGRID*2,y,W_XGRID*25,W_YGRID*3,"Sample Reference Frame");
		sample_ref.callback(callbackToggle, FL_WHEN_CHANGED); 
		sample_ref.color2(FL_BLUE);
		sample_ref.argument(T_SAMPLE_REF_FRAME); 
		y += W_YGRID*3;

		//R ROLLER
		y -= W_YGRID*9;
		Fl_Roller xroller(x+W_XGRID*35,y,W_XGRID*10,W_YGRID*3,"r");
		xroller.align(FL_ALIGN_LEFT);
		xroller.type(FL_HORIZONTAL);
		xroller.range(0.05,25); xroller.step(0.01);
		xroller.value(value[V_R]);
		xroller.callback(callbackViewPos);
		xroller.argument(V_R);
		y += W_YGRID*3;

		

		//THETA ROLLER
		Fl_Roller yroller(x+W_XGRID*35,y,W_XGRID*10,W_YGRID*3,"theta");
		yroller.align(FL_ALIGN_LEFT);
		yroller.type(FL_HORIZONTAL);
		yroller.range(0.005,_PI-0.005); yroller.step(0.005);
		yroller.value(value[V_THETA]);
		yroller.callback(callbackViewPos);
		yroller.argument(V_THETA);
		y += W_YGRID*3;

		//PSI ROLLER
		Fl_Roller zroller(x+W_XGRID*35,y,W_XGRID*10,W_YGRID*3,"psi");
		zroller.align(FL_ALIGN_LEFT);
		zroller.type(FL_HORIZONTAL);
		zroller.range(-2*_PI,2*_PI); zroller.step(0.005);
		zroller.value(value[V_PSI]);
		zroller.callback(callbackViewPos);
		zroller.argument(V_PSI);
		y += W_YGRID*4;

		// output to
		Fl_Menu_Item voxeldispchoices[]= {{"None",0,0,0},
										 {"Points",0,0,0},
										 {0}};
		Fl_Choice voxeldisp(x+W_XGRID*22,y,W_XGRID*17,W_YGRID*3,"Voxel display:");
		voxeldisp.menu(voxeldispchoices);
		voxeldisp.callback(callbackVoxelDisp);
		y+=W_YGRID*6;


	dispctrls.end();


	//===========================================================================
	// EXPERIMENT WINDOW
		
	Fl_Group expborder(397,61,381,306, 0);
	expborder.box(FL_DOWN_BOX);
		expwin = new ExpWindow(400,64,375,300, "Experiment Display", &scope);
		expwin->box(FL_DOWN_BOX);
		expwin->mode(FL_DOUBLE|FL_RGB);
		expwin->end();
	expborder.end();

	ctrlwin.end();




	//===========================================================================
	// DETECTOR WINDOW
	detwin = new DetWindow(450,400,scope.detector.getXRes(),scope.detector.getYRes(), "Detector Display", &scope);
	detwin->size_range(50,50,1000,1000,1,1,1);
	detwin->resizable();
	detwin->end();


	//===========================================================================
	// SHOW WINDOWS
	ctrlwin.show(argc, argv);
	expwin->show(argc, argv);
	detwin->show(argc, argv);

	// INIT WINDOWS:
	expwin->init();
	detwin->init();
	expwin->redraw();
	detwin->redraw();

	cout << " active." << endl;

	return Fl::run();
}


