///////////////////////////////////////////////////////////////////
// s c a n s . c p p                                             
// scanning procedures - boundary mapping, etc.                  
//                                                               
// started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          
///////////////////////////////////////////////////////////////////

#include "interface.h"
#include "scans.h"

//////////////////////////////////////////////////////////////////
// GLOBALS:

double centery, centerz;	// for keeping track of centering info

void translate(void){
	scope.translateSample(value[V_TX], value[V_TY], value[V_TZ]);
}

void rotate(void){
	scope.rotateSample(value[V_OMEGA], value[V_CHI], value[V_PHI]);
}


//////////////////////////////////////////////////////////////////
//
//  AUTOROTATE: Perform sample rotation                          
//
//////////////////////////////////////////////////////////////////

void 
autorotate(void){
	
	double tomega=value[V_OMEGA],tchi=value[V_CHI],tphi=value[V_PHI];

	// CHECK CONDITIONS - IS ROTATION POSSIBLE?
	if((!toggle[T_OMEGA]) && (!toggle[T_CHI]) && (!toggle[T_PHI])){
		fl_alert("You need to activate at least one rotation.");
		return;
	}

	if(toggle[T_OMEGA] && (value[V_OMEGASTEP]==0)) {
		fl_alert("Please specify a non-zero OMEGA step.");
		return;
	}
	else if(toggle[T_CHI] && (value[V_CHISTEP]==0)) {
		fl_alert("Please specify a non-zero CHI step.");
		return;
	}
	else if(toggle[T_PHI] && (value[V_PHISTEP]==0)) {
		fl_alert("Please specify a non-zero PHI step.");
		return;
	}
	
	if (toggle[T_LOG] && toggle[T_UPDATE_DETECTOR]) {
		fl_alert("You can't write a logfile AND update the detector simultaneously");
		return;
	}


	//DO THE ROTATION
	if(toggle[T_LOG]) scope.writeLogHeader("AutoRotate");
	while (1){
	
		if(toggle[T_OMEGA]){  if(value[V_OMEGA] >= value[V_OMEGASTOP]) break;}
		if(toggle[T_CHI]){    if(value[V_CHI]   >= value[V_CHISTOP]) break;}
		if(toggle[T_PHI]){    if(value[V_PHI]   >= value[V_PHISTOP]) break;}			

		//reset detector between frames if the logfile is being written...
		if (toggle[T_LOG]) scope.detector.image->clear();

		// rotate the sample					
		rotate();

		if(toggle[T_UPDATE_VIEW]){ expwin->draw(); expwin->swap_buffers();}
		if(toggle[T_UPDATE_DETECTOR]){ detwin->draw(); detwin->swap_buffers();}
		if(! (toggle[T_BEAM_ON] && toggle[T_UPDATE_VIEW]) )
			scope.activate();  

		// write the log if needed:
		if (toggle[T_LOG]) scope.writeLogEntry();

		if(toggle[T_OMEGA]) value[V_OMEGA]  += value[V_OMEGASTEP];
		if(toggle[T_CHI])	value[V_CHI]    += value[V_CHISTEP];
		if(toggle[T_PHI])   value[V_PHI]    += value[V_PHISTEP];
	}
	
	value[V_OMEGA] = tomega;
	value[V_CHI] = tchi;
	value[V_PHI] = tphi;
	expwin->redraw();
	detwin->redraw();

}


//////////////////////////////////////////////////////////////////
//
//	FIND_PEAK  - rotate the sample until a sufficient 
//		bragg peak is located 
//
//////////////////////////////////////////////////////////////////
void find_peak(void) {

	double intensity;

	//=====================================================
	// CHECK CONDITIONS - IS ROTATION POSSIBLE?

	if((!toggle[T_OMEGA]) && (!toggle[T_CHI]) && (!toggle[T_PHI])){
		fl_alert("You need to activate at least one rotation.");
		return;
	}

	if(toggle[T_OMEGA] && (value[V_OMEGASTEP]==0)) {
		fl_alert("Please specify a non-zero OMEGA step.");
		return;
	}
	else if(toggle[T_CHI] && (value[V_CHISTEP]==0)) {
		fl_alert("Please specify a non-zero CHI step.");
		return;
	}
	else if(toggle[T_PHI] && (value[V_PHISTEP]==0)) {
		fl_alert("Please specify a non-zero PHI step.");
		return;
	}
	
	if (toggle[T_LOG] && toggle[T_UPDATE_DETECTOR]) {
		fl_alert("You can't write a logfile AND update the detector simultaneously");
		return;
	}

	//=====================================================
	//DO THE ROTATION


	intensity = scope.watchDetectorSpot(0,0,scope.detector.getXRes()-1,
												scope.detector.getYRes()-1);

	if(intensity ==0 ){

		if(toggle[T_LOG]) scope.writeLogHeader("Peak Search");
		while (1){
	
			if(toggle[T_OMEGA]){  if(value[V_OMEGA] >= value[V_OMEGASTOP]) break;}
			if(toggle[T_CHI]){    if(value[V_CHI]   >= value[V_CHISTOP]) break;}
			if(toggle[T_PHI]){    if(value[V_PHI]   >= value[V_PHISTOP]) break;}			
	
			// move one step...
			if(toggle[T_OMEGA]) value[V_OMEGA]  += value[V_OMEGASTEP];
			if(toggle[T_CHI])	value[V_CHI]    += value[V_CHISTEP];
			if(toggle[T_PHI])   value[V_PHI]    += value[V_PHISTEP];


			scope.detector.image->clear();			//reset detector
			rotate();

			if(toggle[T_UPDATE_VIEW]){ expwin->draw(); expwin->swap_buffers();}
			if(toggle[T_UPDATE_DETECTOR]){ detwin->draw(); detwin->swap_buffers();}
			if(! (toggle[T_BEAM_ON] && toggle[T_UPDATE_VIEW]) )
				scope.activate();  


			// GET THE INTEGRATED INTENSITY : is it nonzero?
		
			intensity = scope.watchDetectorSpot(0,0,scope.detector.getXRes()-1,
												scope.detector.getYRes()-1);
		
			if (intensity != 0){
				//fl_alert("Peak found: please select it.");
				break;	
			}

			if(toggle[T_OMEGA]) value[V_OMEGA]  += value[V_OMEGASTEP];
			if(toggle[T_CHI])	value[V_CHI]    += value[V_CHISTEP];
			if(toggle[T_PHI])   value[V_PHI]    += value[V_PHISTEP];

			// BAIL OUT IF NO PEAKS ARE APPEARING IN THE DETECTOR 
			// (SAMPLE MAY BE OUT OF THE BEAM)
			if ( (value[V_OMEGA]+value[V_CHI]+value[V_PHI] > 6*PI)) return;

		}
	}

	//TODO: fix this
	// Just a hack to automatically grab the peak:
	// (won't work if there are peaks from multiple grains 
	//  on the detector)
	detwin->sx = 0; detwin->sy=0;
	detwin->sw = scope.detector.getXRes()-1;
	detwin->sh = scope.detector.getYRes()-1;
	callbackConstrainSel(NULL,NULL);

	
	expwin->redraw();
	detwin->redraw();

}


//////////////////////////////////////////////////////////////////
//
//  Helper functions for the boundary search scans:
//
//////////////////////////////////////////////////////////////////

int get_bragg_intensity(void){
// helper function for boundary_search() - 
// returns the intensity of the peak within the
// selection box
	scope.detector.image->clear();
	scope.activate();
	return scope.watchDetectorSpot(detwin->sx,
									detwin->sy,
									detwin->sw+detwin->sx,
									detwin->sh+detwin->sy);
}


void step(double hstep, double vstep, bool move){
	value[V_TY] += hstep;
	value[V_TZ] += vstep;
	if(move) translate();
}

void rotstep(double theta, double &y, double &z){
	double yp = y*cos(theta) - z*sin(theta);
	double zp = y*sin(theta) + z*cos(theta);
	y=yp;
	z=zp;
}

void boundary_mark(ofstream outbound) { 
	double vty, vtz, vtx;
	MatrixStack M;
	Vec4 p;

	// TRANSFORM THE COORDINATES INTO SAMPLE REF FRAME AND OUTPUT...
	M.eulerRot(scope.getOmega(), scope.getChi(), scope.getPhi(), true);
	p = M.getMatrix() * Vec4(value[V_TX], value[V_TY], value[V_TZ],1);
	outbound  << p[X] << "\t" << p[Y] << "\t" << p[Z] << endl;
}


//////////////////////////////////////////////////////////////////
//
//  BOUNDARY_SEARCH: Search grain boundary from present 
//		starting conditions (maps one ring)
//
//////////////////////////////////////////////////////////////////

bool
boundary_search(ofstream outbound){

	double starty = value[V_TY];
	double startz = value[V_TZ];
	double inity, initz, d;

	double hstep = value[V_BOUNDARY_STEPSIZE];
	double vstep = 0;
	double turn = value[V_BOUNDARY_TURN];
	double spread;
	int curint, lastint;
	int halfint;
	int i, hits=0, misses=0;
	bool origsign, cursign;
	int newdir;
	bool fullloop=false;

	double the=0, lastthe=0, the0=0;


	expwin->make_current();

	//===================================================================
	// FIRST, CHECK IF THERE IS A BRAGG PEAK INSIDE THE SELECTION:

	if( (curint=get_bragg_intensity()) == 0) return false;
	halfint = curint/2;
	if (halfint == 0) return false;

	cout << "Boundary Search in progress..."<<endl;
	cout << "        file ext = " << boundaryfileext-1 << endl;
	cout << "  half-intensity = " << halfint << endl;

	//===================================================================
	// LOCATE THE BOUNDARY EDGE

	while(curint >= halfint){
		step(hstep,0, true);
		curint=get_bragg_intensity();
	}
	boundary_mark(outbound); 

	centery=0; centerz=0;

	lastint=curint;
	inity = value[V_TY];
	initz = value[V_TZ];

	//===================================================================
	// STEP AROUND THE BOUNDARY:
	//	start moving tangent to the boundary
	//	if value decreases, move back, rotate ccw and try again
	//	if value increases, move back, rotate cw and try again
	//	keep going until you are back at the initial boundary position again

	hstep = 0; vstep = value[V_BOUNDARY_STEPSIZE];
	the0 = atan2( initz-startz, inity-starty );
	lastthe=the0;


	for(i=0; i<2000; i++){	

		lastint = curint;
		lastthe = the;
		step(hstep,vstep, true);
		curint = get_bragg_intensity();


		expwin->draw();
		expwin->swap_buffers();


		// Where are we in polar angle?
		the = atan2( value[V_TZ]-startz, value[V_TY]-starty );

		if ( ((lastthe-the0)<0) && ((the-the0) > 0) ) {
			cout << "  traced entire boundary." << endl;
			fullloop=true;
			break;
		}
		
		if ( abs(curint-halfint) <= value[V_BOUNDARY_TOLER] ) {  
			// near enough to half int boundary
			boundary_mark(outbound); hits++;
			centery += value[V_TY];	centerz += value[V_TZ];
			continue;
		}
		else if (curint<halfint && lastint>halfint ){
			// fell somewhere in between
			// should turn, but less than normal
			boundary_mark(outbound); hits++;
			centery += value[V_TY];	centerz += value[V_TZ];
			rotstep(turn/2, hstep, vstep);	
		}
		else if (curint>halfint && lastint<halfint ){
			//fell somwehere in between 
			boundary_mark(outbound); hits++;
			centery += value[V_TY];	centerz += value[V_TZ];
			rotstep(-turn/2, hstep, vstep);	

		}
		else if (curint < halfint){
			// move back and rotate ccw
			misses++;
			step(-hstep, -vstep, false);
			rotstep(turn, hstep, vstep);	

		}
		else if (curint > halfint){
			//rotate cw
			misses++;
			step(-hstep, -vstep, false);
			rotstep(-turn, hstep, vstep);	
		}

	} 

	outbound << endl;

	centery = centery/hits;
	centerz = centerz/hits;

	//===================================================================
	// RESTORE TRANSLATIONS
	value[V_TY] = starty;
	value[V_TZ] = startz;
	translate();
	cout << "  HITS: "<<hits<<"  MISSES: "<< misses 
		 << "  EFFICIENCY: " << (double)hits/(hits+misses) << endl;
	if (!fullloop) {
		cout<<"  NOTE: boundary search was unable to close the boundary, so there" << endl
			<<"     may be points missing, or the boundary may have been traced" <<endl
			<<"     multiple times." << endl;
	}		
	cout << "Done." << endl;

	return true;

}


//////////////////////////////////////////////////////////////////
//
//  MAP_GRAIN: map the entire boundary of the grain which is
//		currently diffracting (produces full 3-d map)
//
//////////////////////////////////////////////////////////////////
void 
map_grain(void){

	double startx = value[V_TX];
	double starty = value[V_TY];
	double startz = value[V_TZ];
	ofstream outbound;
	char filename[255];
	double cy, cz;
	int curint, halfint;


	sprintf( filename,"boundary.%3.3d",boundaryfileext++ );
	outbound.open( filename, ios::out );

	//===================================================================
	// FIRST, CHECK IF THERE IS A BRAGG PEAK INSIDE THE SELECTION:

	if( (curint=get_bragg_intensity()) == 0) {
		fl_alert("There doesn't appear to be a Bragg peak inside the selection. Either select one now or use FIND PEAK to locate one");
		outbound.close();
		return;
	}
	halfint = curint/2;
	if (halfint == 0) {
		fl_alert("The currently selected peak is not bright enough, please choose a more intense peak");
		outbound.close();
		return;
	}


	//=================================================================
	// SET THE STARTING POSITION AND MAP IN POSITIVE DIRECTION: 
	//   (if there is a prev mapped boundary, use its center position

	while (boundary_search(outbound)){

		cout << "Recentered to: (y,z)=" << value[V_TY] <<","<<value[V_TZ]<<endl;
		value[V_TY] = centery;
		value[V_TZ] = centerz;
		value[V_TX] += 0.002;
		translate();
		outbound << endl;
	}

	//=================================================================
	// NOW, GO BACK AND DO THE NEGATIVE DIRECTION...	
		
	value[V_TX] = startx;
	value[V_TY] = starty;
	value[V_TZ] = startz;
	translate();

	while (boundary_search(outbound)){
		value[V_TY] = centery;
		value[V_TZ] = centerz;
		value[V_TX] -= 0.002;
		translate();
		outbound << endl;
	}

	//=================================================================
	// RESET STARTING VALUES...

	value[V_TX] = startx;
	value[V_TY] = starty;
	value[V_TZ] = startz;
	translate();

	outbound.close();

	cout << "GRAIN MAPPING COMPLETE." << endl;

}





//////////////////////////////////////////////////////////////////
//
//	MAP_CRYSTAL - map central grain and its neighbors
//	
//////////////////////////////////////////////////////////////////

void
map_crystal(void) {
//TODO: fix find_peak so that it returns if no peak found
//			(i.e. if sample leaves the beam...)

	int dir[6][3]={{0,0,1},				// DIRECTION TO MAP IN
				   {0,0,-1},
				   {0,1,0},
				   {0,-1,0},
				   {1,0,0},
				   {-1,0,0}};
	int i,j;						
	double center[3] = {0,0,0};			// STARTING CENTER POSITION
	double rot[3]    = {0,0,0};			// STARTING ROTATATION
	MatrixStack M;						
	Vec4 direct;						// DIRECTION FOR EACH STEP

	//=================================================================
	// FIRST, ORIENT THE CRYSTAL SO THAT IT IS DIFFRACTING:
	//  (and remember the starting positions)

	find_peak();

	center[X] = value[V_TX];
	center[Y] = value[V_TY];
	center[Z] = value[V_TZ];
	rot[OMEGA]= value[V_OMEGA];
	rot[CHI]  = value[V_CHI];
	rot[PHI]  = value[V_PHI];

	//=================================================================
	// MAP THE CENTRAL GRAIN:

	cout << "Mapping central grain..." << endl;
	map_grain();


	//=================================================================
	// FOR EACH DIRECTION FROM THE CENTRAL GRAIN 
	//	(up, down, left, right, in, out)
	//   MOVE TO THE NEW GRAIN AND MAP IT...

	for (i=0; i<6; i++) {

		cout << "Mapping border grain "<<i<<" of 6..." << endl;

		// FIRST, MOVE TO THE CENTER GRAIN AND DIFFRACT...

		value[V_TX] = center[X];
		value[V_TY] = center[Y];
		value[V_TZ] = center[Z];
		translate();
		value[V_OMEGA] = rot[OMEGA];
		value[V_CHI]   = rot[CHI];
		value[V_PHI]   = rot[PHI];		
		rotate();
		find_peak();


		cout << "!!!!!: intensity="<<get_bragg_intensity() << endl;

		// COMPUTE THE MOVEMENT STEPS ACCOUNTING FOR THE 
		// ROTATION OF THE SAMPLE:

		M.eulerRot( value[V_OMEGA],value[V_CHI],value[V_PHI],true );
		direct = M.getMatrix() * 
				(Vec4(dir[i][X],dir[i][Y],dir[i][Z],1)
					*value[V_BOUNDARY_STEPSIZE]);

		// NOW MOVE IN CURRENT DIRECTION UNTIL DIFFRACTION 
		// IS GONE...
		
		cout << "MOVING TO NEW GRAIN: in direction" << direct << endl;

		while ( get_bragg_intensity() > 0 ) {

			value[V_TX] += direct[X];
			value[V_TY] += direct[Y];
			value[V_TZ] += direct[Z];
			translate();
			cout <<"MOVING: "<<value[V_TX]<<","<<value[V_TY]<<"," <<value[V_TZ]<<","<<endl; 

		}

		// TAKE A FEW MORE STEPS TO ENSURE BEAM IS FULLY IN THE GRAIN:

		for (j=0; j<2; j++){

			value[V_TX] += direct[X];
			value[V_TY] += direct[Y];
			value[V_TZ] += direct[Z];
			translate();

		}

		// MAP THE NEW GRAIN:

		scope.display();
		cout << "\tLooking for peak..." << endl;
		find_peak();
		cout << "\tMapping..." << endl;
		map_grain();


	}

}

