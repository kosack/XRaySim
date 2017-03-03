/*****************************************************************
 * x  r  a  y  s  i  m  .  c  p  p                               *
 *                                                               *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/

///////////////////////////////////////////////////////////////////
// INFO:
//
// This program simulates a 3-D X-ray microscope.  In such a setup
// an x-ray beam is passed through a crystal sample which is made
// up of various crystal "grains" of differing orientation. As the
// beam travels through the sample, it is diffracted by the grains,
// producing a pattern on a 2-d detector which is behind the sample.
// 
// In this simulation, a crystal grain is represented by a group of 
// cubic volumetric pixels (voxels), through which the path of a ray
// is traced. See manual for more info.
//
// REQUIREMENTS:
// 
// To compile this program, you will need to have the OpenGL and Fltk
// libraries installed if you want 3D visualization/GUI. This program also 
// uses the Simple Vector Library (SVL) by Andrew Wilmott, which is 
// included with the xraysim source. The source for the Fltk library is
// also included.   
//
// TO DISABLE OPENGL IF NOT NEEDED: comment out the line:
//             #define USE_OPENGL
//		in xraysim.h and all opengl code will not be compiled in.
//		Also, remove the opengl libraries from the project.
//
// TO DISABLE THE GUI:
//	 To completely remove it, don't link with the FlTk libraries, 
//	 remove "interface.cpp" and "interface.h" from the project, and
//	 do not call interface_start() in main.
//
//	 Otherwise, a simpler method is to just remove the line which 
//	 calls interface_start() in main and the program will not use the GUI.
//	 (though it will still be compiled in)

/////////////////////////////////////////////////////////////////////////
// INCLUDES:

#include "xraysim.h"

bool DEBUG = 0;
bool GRAPHICS = 1;
double PI;


/////////////////////////////////////////////////////////////////////////
// GLOBALS:

XRayScope scope;
//GLUquadricObj *quadric;

/////////////////////////////////////////////////////////////////////////
// Image METHODS:

Image::Image(int w, int h){
	width = w; height=h;
	if(!(data = new int[width*height])){				// 3 bits per pixel

		cerr << "Image::Image: ERROR: couldn't allocate image buffer!"<< endl;
		exit(1);

	}
	this->clear();
}


void
Image::display(void){

	#ifdef USE_OPENGL // only do this if OPENGL interface is needed
	glPixelZoom(2,2);
	glPixelTransferi(GL_RED_SCALE, 100000000); 
	glPixelTransferi(GL_GREEN_SCALE, 100000000); 
	glPixelTransferi(GL_BLUE_SCALE, 100000000); 
	glDrawPixels( width, height, GL_LUMINANCE, GL_INT, data );
	#endif

}

void
Image::clear(void){
	int i;
	for( i=0; i< (width*height); i++){			// Clear the image
		data[i] = 0;
	}
	debug << "Image was cleared." << endl;
}

void 
Image::setPixel(int x, int y){
	// average a color into a pixel - (color range should be 0..1)
	double a = 0.3;  

	if (x >= width || x<0 || y >= height || y <0){
		cout << "Image::setPixel: Index ("<<x<<","<<y<<") out of range!"<<endl;
		return;
	}

	data[x+y*width]++;

}

int
Image::getPixel(int x, int y){
	return data[x+y*width];
}


bool
Image::write(char *file){
	// save image to a file
	// output intensity and position

	return false;
}

/////////////////////////////////////////////////////////////////////////
// Detector METHODS:

bool
Detector::intersect(const Ray &ray, Isect &isect){
// Calculate intersection point of ray and detector
// don't really care about returning the isect, since the
// detector autorecords intersection

	double d, tden, t;
	Vec4 point;
	int i;

	// Does the ray intersect the plane which contains the detector?
	tden = (nor[X]*ray.dir[X] + 
			nor[Y]*ray.dir[Y] + 
			nor[Z]*ray.dir[Z]);

	//debug << "Detector: tden=" << tden << endl;
	if (-TOLER<tden && tden< TOLER) return false;  

	// Where does the ray intersect?
	d =	  -(dot(nor,points[0]));
	t =   -(dot(nor,ray.orig) + d) / tden; 

	//debug << "Detector: t=" << t << " d="<<d<< endl;
	if (t<=TOLER) return false;
	
	point = ray.orig + t * ray.dir;

	bool test=1;
	for (i=0; i<3; i++){  // iterate over veritices
		test &= (dot((cross(points[i+1]-points[i], nor,Vec4(0,0,0,1))),
					(point - points[i])) < 0);
	}

	if (!test) return false;

	//debug << "Ray hit detector at: " << point << endl;

	// Now figure out the 2-d coordinate relative to the detector
	Vec4 r, u,v;
	int x,y;
	r = point - points[3]; r[3] = 0;
	v = points[0] - points[3]; v[3]=0;
	u = points[2] - points[3]; u[3]=0;

	x = (int) ((dot(r,u)/width/width)  * (float)image->getWidth());
	y = (int) ((dot(r,v)/height/height)  * (float)image->getHeight());

	debug << "Detector::intersect: x,y = " << x<<","<<y<<endl;

	if (x<0 || y<0) return false;


	image->setPixel(x,y);
	isect.x = x;
	isect.y = y;

	return true;


}

/////////////////////////////////////////////////////////////////////////
// Grain METHODS:

bool
Grain::even(int i){
	if (i < 0) i *= -1;
	if (i==0) return true;
	if (floor(i/2) == floor((i+1)/2)) return true;
	else return false;
}

Grain::Grain(Mat4 m, double ta, GrainType gt, XRayScope *parent){
	orientation = m;	// Set rotation matrix
	a = ta;				// set lattice constant	
	type = gt;			// set type
	scope = parent;		// set parent pointer

	Vec4 B[3];			// basis vectors
	Vec4 *tG;
	double tpa=2*_PI/a;
	int h,k,l;
	bool skip =0;

	//set a pseudocolor:
	pcolor = Vec3(((double)rand()/RAND_MAX),((double)rand()/RAND_MAX),((double)rand()/RAND_MAX));


	// CALCULATE A LIST OF G VALUES:
	//	G = h*B1 + k*B2 + l*B3      where Bi are the basis vectors

	//	Bi = (2*PI/a) * ihat ... multiplied by orientation matrix
	B[X] = orientation * Vec4(tpa,0,0,0);
	B[Y] = orientation * Vec4(0,tpa,0,0);
	B[Z] = orientation * Vec4(0,0,tpa,0);
	
	// iterate over miller indices
	// (apply selection rules to the G's)

	for (h=-MAX_MILLER_INDEX; h<=MAX_MILLER_INDEX; h++){
		for (k=-MAX_MILLER_INDEX; k<=MAX_MILLER_INDEX; k++){
			for (l=-MAX_MILLER_INDEX; l<=MAX_MILLER_INDEX; l++){

				skip=false;
				if (h==0&&k==0&&l==0) skip=true;
				switch (type) {
					case SC: 
						// Selection rule: all
						break;
					case FCC:
						// Selection rule:  hkl all even or all odd
						if (! ((even(h) && even(k) && even(l)) ||
							  (!(even(h)) && (!even(k)) && (!even(l))))){
							skip=true;		
						}
						break;
					case BCC:
						// Selection rule: h+k+l even
						if (! even(h+k+l) ) skip = true; 
						break;				

				}

				if(!skip){

					// Add G to the GList linked list...
					tG = new Vec4();
					*tG = h*B[X] + k*B[Y] +l*B[Z];
					glist.addItem(tG);

					tG = new Vec4();
					*tG = -(h*B[X] + k*B[Y] +l*B[Z]);
					glist.addItem(tG);
//					cout << "Found a G at: ["<<h<<k<<l<<"]" << endl;

				}
			}
		}
	}

	debug << "Number of G's in list = " << glist.getLength() << endl;
}


bool
Grain::getDiffraction(const Ray &r, List<Vec4> &klist){
// Figure out if diffraction occurs and in which direction
// should return a list of diffracted rays...

	Vec4 *g;
	Vec4 *kf;
	Vec4 ki;
	bool diff=false;
	double test;
	double diverge;

	if(glist.getLength() > 0){

		diverge = (2*_PI/r.lambda);
		ki = r.dir * diverge;

		diverge = diverge * diverge * scope->getBeamDiv();
		
		glist.beginIteration();
		while (g = glist.getItem()){
		
			// DIFFRACT OR NOT?

			test = abs(dot(*g, *g)+2*dot(ki,*g));

			if (test <= diverge) {  // condition for scattering
				debug << "Grain::getDiffraction: Ray was diffracted!  " <<  ki+*g << endl;

				kf = new Vec4();
				*kf = norm(ki + *g); 
				klist.addItem(kf);

				diff = true;
			}

		}

	}
	else { 
		cout << "Grain::getDiffraction(): No G's found"<<endl;
		return false;
	}

	if (diff)	return true;
	else return false;

}



/////////////////////////////////////////////////////////////////////////
// Volume METHODS:

int toint(double a){
	if( (a-(double)floor(a)) <= 0.5 ) return floor(a);
	else return ceil(a);
}


Volume::Volume(double txlen,double tylen,double tzlen,double tvoxsize){
	voxel = NULL;
	voxel_disp = 0;
	resize(txlen, tylen, tzlen, tvoxsize);
}


Volume::~Volume(void){
	int i;

	delete [] voxel;

	for(i=0; i<numgrains; i++){
		delete grainlist[i];
	}

}


int
Volume::resize(double txlen, double tylen, double tzlen, double tvoxsize){
	int i;

	xlen=txlen; ylen=tylen; zlen=tzlen; voxsize=tvoxsize;
	xvox = toint(xlen/voxsize);
	yvox = toint(ylen/voxsize);
	zvox = toint(zlen/voxsize);	
	totalvox = xvox*yvox*zvox;

	posx = -(xlen/2);
	posy = -(ylen/2);
	posz = -(zlen/2);


	if(voxel != NULL) delete [] voxel;

	if(!(voxel = new Voxel[totalvox])){
		cerr << "FAILED TO ALLOCATE ENOUGH MEMORY FOR VOLUME BUFFER." << endl;
		cerr << "---------------------------------------------------" << endl;
		cerr << "     TOTAL VOXELS: " <<totalvox<<endl;
		cerr << "  MEMORY REQUIRED: "<<totalvox*sizeof(Voxel)/1024<<" KB"<<endl;
		exit(1);
	}

	if(voxel == NULL) {
		cerr << "HEY! What happened?" << endl;
		exit(1);
	}

	for (i=0; i<totalvox; i++)
		voxel[i] = 0;
	
	cout << "CREATED NEW VOLUME: " << totalvox << " voxels"
		 << "; "<< totalvox*sizeof(Voxel)/1024 << " KB"<<endl;

	return totalvox;

}



Vec4
Volume::centerOfVoxel(int i, int j, int k) {
	return Vec4(posx+i*voxsize+voxsize/2,
				posy+j*voxsize+voxsize/2,
				posz+k*voxsize+voxsize/2, 1);
}


Voxel 
Volume::voxelAt(int i, int j, int k){

	assert( (i+j*xvox+k*xvox*yvox) < totalvox);
	return voxel[i + j*xvox + k*xvox*yvox];
}


bool
Volume::voxelIndexOf( Vec4 p, int &i, int &j, int &k){
	
	Vec3 pp = Vec3(p[X] - posx, p[Y]-posy, p[Z]-posz);
	
	i = (int)(pp[X] / voxsize);
	j = (int)(pp[Y] / voxsize);
	k = (int)(pp[Z] / voxsize);		

	if(i>xvox || j>yvox || k>zvox || i<0 || j<0 || k<0){
		i=j=k=0;
		cerr << "Volume::voxelIndexOf(): VOXEL INDEX OUT OF RANGE!" << endl;
		return false;
	}

	return true;
}

void
Volume::setVoxel(int i, int j, int k, Voxel vox){
	voxel[i + j*xvox + k*xvox*yvox] = vox;
}


void
Volume::display(void){

#ifdef USE_OPENGL // SHOW DIFFRACTED RAY
	
	double minx, miny, minz, maxx, maxy, maxz;
	int i,j,k;

	minx = posx;
	miny = posy;
	minz = posz;

	maxx = posx + xlen;
	maxy = posy + ylen;
	maxz = posz + zlen;

	glColor4f(0.0,0.6,0.9,1);

	glBegin(GL_LINE_LOOP); // xy, high z
		glNormal3d( 0,0,1);
		glVertex3d( maxx,miny,maxz );
		glVertex3d( maxx,maxy,maxz );
		glVertex3d( minx,maxy,maxz );
		glVertex3d( minx,miny,maxz );
	glEnd();
	glBegin(GL_LINE_LOOP); // xy, low z
		glNormal3d( 0,0,-1);
		glVertex3d( minx,miny  ,minz   );
		glVertex3d( minx,maxy  ,minz   );
		glVertex3d( maxx,maxy  ,minz   );
		glVertex3d( maxx,miny  ,minz   );
	glEnd();
	glBegin(GL_LINE_LOOP); // xz, high y
		glNormal3d( 0,1,0);
		glVertex3d( maxx  ,maxy  ,minz   );
		glVertex3d( maxx  ,maxy  ,maxz   );
		glVertex3d( minx  ,maxy  ,maxz   );
		glVertex3d( minx  ,maxy  ,minz   );
	glEnd();
	glBegin(GL_LINE_LOOP); // xz, low y
		glNormal3d( 0,1,0);
		glVertex3d( minx  ,miny  ,minz   );
		glVertex3d( maxx  ,miny  ,minz   );
		glVertex3d( maxx  ,miny  ,maxz   );
		glVertex3d( minx  ,miny  ,maxz   );
	glEnd();
	glBegin(GL_LINE_LOOP); // yz, high x
		glNormal3d( -1,0,0);
		glVertex3d( maxx  ,miny  ,minz   );
		glVertex3d( maxx  ,maxy  ,minz   );
		glVertex3d( maxx  ,maxy  ,maxz   );
		glVertex3d( maxx  ,miny  ,maxz   );
	glEnd();
	glBegin(GL_LINE_LOOP); // yz, low x
		glNormal3d( -1,0,0);
		glVertex3d( minx  ,maxy  ,minz   );
		glVertex3d( minx  ,maxy  ,maxz   );
		glVertex3d( minx  ,miny  ,maxz   );
		glVertex3d( minx  ,miny  ,minz   );
	glEnd();

	// SHOW THE INTERNAL VOXELS:

	if(voxel_disp != 0) {

		glBegin(GL_POINTS);

		if(voxel_disp == 1){ // show all points
			for (k=1; k<zvox; k++){
				for (j=1; j<yvox; j++){
					for (i=1; i<xvox; i++){
						glColor3dv(grainlist[voxelAt(i,j,k)]->pcolor.Ref());
						glVertex3d(posx+i*voxsize,posy+j*voxsize, posz+k*voxsize);
					}
				}
			}
		}
		else { // just show boundary

		}


		glEnd();
	}	

#endif	// end opengl

}


#define TEST(AXIS, LIMIT, OPER, TFUN1, TFUN2){							\
		if ((cur_t = (LIMIT-ray.orig[AXIS])/ray.dir[AXIS]) OPER t){		\
			p=ray.point(cur_t - TOLER);									\
			if(TFUN1(p) && TFUN2(p))									\
			t=cur_t;													\
		}																\
	}																	\
;


#define within_x(P)	( (P[X] >= minx ) && (P[X] < maxx) )

#define within_y(P)	( (P[Y] >= miny ) && (P[Y] < maxy) )	

#define within_z(P) ( (P[Z] >= minz ) && (P[Z] < maxz) )


bool
Volume::intersect(int i, int j, int k, const Ray &ray, List<Isect> *isectlist){

	double t = INFINITY;
	double cur_t;
	Vec4 p;
	Isect *isect;
	double minx, miny, minz, maxx, maxy, maxz;

	if( voxelAt(i,j,k) == EMPTY_VOXEL ) return false; 

	minx = posx + i*voxsize;		// lower left corner
	miny = posy + j*voxsize;
	minz = posz + k*voxsize;

	maxx = posx + (i+1)*voxsize;	// upper right corner
	maxy = posy + (j+1)*voxsize;
	maxz = posz + (k+1)*voxsize;

	
	// find the entry point:
	
	// check Z faces:
	if (ray.dir[Z]){
		TEST(Z, maxz, <, within_x, within_y); 
		TEST(Z, minz, <, within_x, within_y); 
	}

	// check Y faces:
	if (ray.dir[Y]){
		TEST(Y, maxy, <, within_x, within_z); 
		TEST(Y, miny, <, within_x, within_z); 
	}

	// check X faces:
	if (ray.dir[X]){
		TEST(X, maxx, <, within_y, within_z); 
		TEST(X, minx, <, within_y, within_z); 
	}

	if (t>=INFINITY) return false;  // didn't hit

	else { // RAY HIT VOXEL...

		isect = new Isect();
		isect->t = t;
		isect->i = i; isect->j=j; isect->k=k;
		isect->grainid = voxelAt(i,j,k);
		isect->enter = true;
		isectlist->addItem(isect);

		assert( isect->grainid >= 0);

		// SHOW THE BEAMPATH
		glColor4f(1,1,0,1);
		glBegin(GL_POINTS);
		glVertex4dv(ray.point(t).Ref());
		glEnd();

		return true;
	}


}


bool
Volume::boundingVoxels(const Ray &ray,int &i0,int &j0,int &k0, 
                       int &i1,int &j1,int &k1, double bmin ,double bmax, double detdist){
// Computes the entry and exit point in the volume and returns
// which voxel is at those points.  (i0, j0, k0) is the entry point
// and (i1, k1, j1) is the exit point.
// bmin and bmax are the probe volume bounds (only look at voxels within them)
// TODO: merge this code with the intersect method

	double t = INFINITY, t0=0, t1=0, bt0=0, bt1=0;
	double cur_t;
	Vec4 p;
	Isect *isect;
	double minx, miny, minz, maxx, maxy, maxz, deltat;

	minx = posx;		// lower left corner
	miny = posy;
	minz = posz;
	maxx = posx + xlen;	// upper right corner
	maxy = posy + ylen;
	maxz = posz + zlen;

	//===========================================================
	// FIND ENTRY POINT:

	// check Z faces:
	if (ray.dir[Z]){
		TEST(Z, maxz, <, within_x, within_y); 
		TEST(Z, minz, <, within_x, within_y); 
	}

	// check Y faces:
	if (ray.dir[Y]){
		TEST(Y, maxy, <, within_x, within_z); 
		TEST(Y, miny, <, within_x, within_z); 
	}

	// check X faces:
	if (ray.dir[X]){
		TEST(X, maxx, <, within_y, within_z); 
		TEST(X, minx, <, within_y, within_z); 
	}
	if (t==INFINITY) return false;  // didn't hit
	t0 = t;
	
	//===========================================================
	// FIND EXIT POINT:

	// check Z faces:
	if (ray.dir[Z]){
		TEST(Z, maxz, >, within_x, within_y); 
		TEST(Z, minz, >, within_x, within_y); 
	}

	// check Y faces:
	if (ray.dir[Y]){
		TEST(Y, maxy, >, within_x, within_z); 
		TEST(Y, miny, >, within_x, within_z); 
	}

	// check X faces:
	if (ray.dir[X]){
		TEST(X, maxx, >, within_y, within_z); 
		TEST(X, minx, >, within_y, within_z); 
	}
	t1 = t;

	
	//===========================================================	
	// APPLY PROBE VOLUME BOUNDS TO THE INTERSECTION

	bt0 = bt1 = detdist; 
	bt0 += bmin;
	bt1 += bmax;

	if (bt1 < t0 || bt0 > t1) return false; // outside volume
	if (bt0 < t0 || bt0 > t1) bt0 = t0;
	if (bt1 > t1 || bt1 < t0) bt1 = t1;


	//===========================================================	
	// CALCULATE VOXEL AT ENTRY/EXIT
	
	deltat = voxsize/16;
	if(!voxelIndexOf(ray.point(bt0+deltat), i0,j0,k0)){
		return false;
	}
	if(!voxelIndexOf(ray.point(bt1+deltat), i1,j1,k1)){
		return false;
	}

	return true;
}



void
Volume::load(char *file){
// Loads volume data from a file
// assumes that grainlist has already been built by XRayScope
// (each voxel stores an index to the grain array)

	
	// TODO: write me.  
	// see Volume::save(..)

}

void
Volume::save(char *file){

	ofstream outfile;

	//TODO: fix me I'm broken
	outfile.open(file);	
	outfile << "XRAYSIM-BINARY-VOLUME-DATA";
	outfile.write((char *)totalvox, sizeof(int));
	outfile.write((char *)voxel,sizeof(Voxel)*totalvox);
	outfile.close();

}

void
Volume::generate(VolumeType type){
	// Generate a sample polycrystal

	int i,j,k;
	double minx, miny, minz, maxx, maxy, maxz;
	Voxel n, nearest;
	Voxel g;
	Vec3 *rgrain;
	double neardist,tmp;

	//if (voxel != NULL) delete [] voxel;

	switch(type) {

		case VOLUME_CUBE:
	
			g=0;
			for(k=0; k<zvox; k++){				
				for(j=0; j<yvox; j++){
					for(i=0; i<xvox; i++){
						if (g>=(numgrains-1)) g=0;
						cout << "g=" << g<<endl;
						setVoxel(i,j,k, g);
						g++;
					}
				}
			}				

			
			break;

		case VOLUME_RANDOM:
		// for random algorithm:

		// For each grain, place a random point in a 3-d volume
		// then, build a set of voxels that fill the bounding box, and
		// set the grain for each to be the grain of the random point
		// closest to it.  This generates an approximately physical random
		// polycrystal.
			cout << "Building random polycrystal..." << endl;
			rgrain = new Vec3[numgrains];
			minx = posx;
			miny = posy;
			minz = posz;
			maxx = posx + xlen;
			maxy = posy + ylen;
			maxz = posz + zlen;			
			//srand(10);
			rgrain[0] = Vec3(0,0,0);
			for(n=1; n<numgrains; n++){
				rgrain[n] =Vec3(  (((double)rand()/RAND_MAX)*(maxx-minx) + minx),
								  (((double)rand()/RAND_MAX)*(maxy-miny) + miny),
								  (((double)rand()/RAND_MAX)*(maxz-minz) + minz));
			}
					

			for(k=0; k<zvox; k++){				
				for(j=0; j<yvox; j++){
					for(i=0; i<xvox; i++){

						// find nearest random point
						nearest = 0;
						neardist = 1000000;
						for(n=0; n<numgrains; n++){
							if ((tmp=len(rgrain[n]-centerOfVoxel(i,j,k))) < neardist){
								neardist = tmp;
								nearest = n;
							}
						}
						
						// set it
						setVoxel(i,j,k, nearest);

					}
				}
				if (k % 10 == 0) cout << "\r   " << (double)k/zvox*100 <<"%"<< flush;
			}
			delete [] rgrain;
			cout << "\r   " << (double)k/zvox*100 <<"%"<< endl;
			cout << "done." << endl;
			break;
	
	}

}

/////////////////////////////////////////////////////////////////////////
// XRayScope METHODS:

XRayScope::XRayScope(void){
	tx=ty=tz=0;
	omega=chi=phi=0;
	debug << "XRayScope was created." << endl;


	sample = new Volume(.05, .05, .05, 0.05);
	sample->numgrains=0;
	boundmin = -0.005;
	boundmax = 0.005;

	Mat4 m = vl_1;		// initialize matrix
	beamdiv = 0.01;		// set default beam divergence

	strcpy(logfilename, "scope.log");
	logfile.open(logfilename, ios::out);
	if (logfile.fail()) {
		cout << "Couldn't open logfile." << endl;
		exit(1);
	}
	logfile << "# XRAYSIM-LOGFILE "<< VERSION << endl;

}

XRayScope::~XRayScope(void){
	
	logfile.close();
	debug << "XRayScope was destroyed." << endl;

}


void
XRayScope::setLogFilename(const char *filename) {
// set a new logfilename;

	// close the current log and change the filename:
	logfile.close();
	strcpy( logfilename, filename );

	// open the new file
	logfile.open(logfilename);
	if (logfile.fail()){
		cerr << "ERROR: couldn't open log file: " << logfilename << endl;
	}

	logfile << "# XRAYSIM-LOGFILE "<< VERSION << endl;

}


void
XRayScope::loadSample(char *file){
	// load in sample datafile

	ifstream infile;
	char s[255];
	double version;
	int i;
	GrainType type;
	double omega=0,chi=0,phi=0, a=0;
	double xl=0, yl=0, zl=0, vs=0;

	char string[255];
	char delimiters[] = " ,=\t";
	char *token, *key;
	double value;
	bool startflag=true;	
	
	debug << "XRayScope::loadSample() - not totally implemented yet."<<endl;

	infile.open(file);
	if(infile.fail()){
		cout << "XRayScope::loadSample: couldn't open '" << file <<"' for input."<< endl;
		exit(1);
	}

	// LOAD HEADER INFO:
	cout << "Loading data from: " << file << endl;
	infile >> s;
	if (strcmp(s,"XRAYSIM-SAMPLE")){
		cout << "'"<<file<<"' does not appear to be an XRaySim datafile!" << endl;
		cout << s << endl;
		exit(1);
	}
	infile >> version;
	if(version > VERSION){
		cerr << "WARNING: the sample data you are loading appears to be from" << endl
			 << "   a newer version of XRaySim than the one that you are" << endl
			 << "   currently runnning: ("<<VERSION<<").  Please use at least"<<endl
			 <<	"   version " << version<< " to load this datafile." << endl << endl; 
		exit(1);
	}

	// DELETE ANY GRAINS THAT CURRENTLY EXIST...
	
	for (i=0; i<sample->numgrains; i++)
		if (sample->grainlist[i] != NULL) delete sample->grainlist[i];
	sample->numgrains = 0;	

	// NOW START PROCESSING TOKENS:
	//   for each line...
	while(! infile.eof() ){
		infile.getline(string, 255, '\n');
		startflag = true;

		if(startflag) {token = strtok( string, delimiters ); startflag=false;}
		else token = strtok( NULL, delimiters );

		// COMMENT:
		if (!strcmp(string,"")) continue;
		if (token[0] == '#') continue;

		
		// KEY/VALUE PAIR:
		else if (!strcmp(token, ">")) {

			key = strtok( NULL, delimiters );
			value = atof( strtok( NULL, delimiters ) );		

			// process key:
			if(!strcmp(key,"lambda")) emitter.ray.lambda = value;
			else if(!strcmp(key,"width")) xl=value;
			else if(!strcmp(key,"height")) yl=value;
			else if(!strcmp(key,"depth")) zl=value;
			else if(!strcmp(key,"voxsize")) vs=value;
			else if(!strcmp(key,"numgrains")) sample->numgrains=value;

			cout << "VALUE: " << key << " = " << value << endl;
			continue;

		}

		//GRAIN INFO
		else {
			omega = atof(token);
			chi	  = atof(strtok( NULL, delimiters ));
			phi	  = atof(strtok( NULL, delimiters ));
			a	  = atof(strtok( NULL, delimiters ));
			token = strtok(NULL, delimiters);

			cout << "grain " << sample->numgrains << " is "<< token << " lattice, a="<<a
				 << " rotated at: " << omega<<","<<chi<<","<<phi<<endl;

			if (strcmp(token, "FCC")==0) type = FCC; 
			else if (strcmp(token, "BCC")==0) type = BCC; 
			else if (strcmp(token, "SC")==0) type = SC; 

			matrix.push();
			matrix.loadIdentity();
			matrix.eulerRot(omega*DEGTORAD,chi*DEGTORAD,phi*DEGTORAD,false);	// rotate sample
			sample->grainlist[sample->numgrains] = new Grain(matrix.getMatrix(), a, type, this); 
			matrix.pop();

			sample->numgrains++;

		}

	}


	if(xl==0 || yl==0 || zl==0 || vs==0){
		cerr << "PROBLEM: Sample has a dimension which is 0 in length! " << endl
			 << "    Sample NOT loaded" << endl;
		return;
	}
	if (a==0) {
		cerr << "PROBLEM: No lattice constant was specified in the datafile." << endl
			 << "    Sample NOT loaded." << endl;
	}
	if (emitter.ray.lambda==0) {
		cerr << "PROBLEM: No wavelength was specified in the datafile." << endl
			 << "    Sample NOT loaded." << endl;
	}

	sample->resize(xl, yl, zl, vs);


}

void
XRayScope::rotateSample(double tomega, double tchi, double tphi){
	// set sample rotation

	omega = tomega;
	chi = tchi;
	phi   = tphi;

	redraw();

}

void
XRayScope::translateSample(double ttx, double tty, double ttz){
	// set sample translation

	tx =ttx; ty=tty; tz=ttz;

	redraw();

}

void
XRayScope::setDetectorDist(double dist){ detdist=dist; redraw(); }

void
XRayScope::setProbeVolumeSize(double min, double max) {
// set size of area where beam has effect
	boundmin = min;
	boundmax = max;
}

void
XRayScope::setDetectorSize(double w, double h){
	detector.setSize(w,h);
	redraw();
}

void
XRayScope::setDetectorRes(int x, int y){
	detector.setRes(x,y);
	redraw();
}


void
XRayScope::redraw(void){
	// Recompute the positions of everything:

	int i;
	Mat4 M;

	debug << "Recomputing Scene..." << endl;

	double hw=detector.getWidth()/2;
	double hh=detector.getHeight()/2;

	// first reset everything to defaults:
	detector.points[0] = Vec4(detdist,-hw, hh ,1); // UL
	detector.points[1] = Vec4(detdist, hw, hh ,1); // UR
	detector.points[2] = Vec4(detdist, hw,-hh ,1); // LR
	detector.points[3] = Vec4(detdist,-hw,-hh ,1); // LL
	detector.nor       = Vec4(-1,0,0,1);		 // normal

	emitter.ray.orig = Vec4(-detdist,0,0,1); // origin
	emitter.ray.dir  = Vec4(1,0,0,0); // xhat direction

	// reset the matrix:
	matrix.loadIdentity();

	// set up the "sample" rotation and translation matrix
	// (actually rotate and translate the world, leaving sample fixed)
	matrix.push();
	matrix.eulerRot(omega, chi, phi, true);  // perform inverse Euler rotation

	matrix.push();
	matrix.translate(-tx,-ty,-tz);
	

	// transform emitter/detector:
	M=matrix.getMatrix();	
	for(i=0; i<4; i++){
		detector.points[i] = M * detector.points[i];
	}
	emitter.ray.orig = M * emitter.ray.orig;

	matrix.pop(); // don't want to translate the direction
	M=matrix.getMatrix();	
	emitter.ray.dir = M * emitter.ray.dir;  
	detector.nor    = M * detector.nor;
	
	matrix.pop();	// back to identity;

}


void
XRayScope::writeLogHeader(char *desc){
	logfile << "# ------------------------------------------------------------------"<<endl;
	logfile << "# DESCRIPTION: " << desc << endl;
	logfile << "# KEYS/VALUES:"<< endl;
	logfile	<< "> xres=" << detector.image->getWidth() << endl;
	logfile	<< "> yres=" << detector.image->getHeight()<< endl;
	logfile	<< "> width=" << detector.getWidth() << endl;
	logfile	<< "> height=" << detector.getHeight() << endl;
	logfile	<< "> dist=" << this->detdist << endl;
	logfile	<< "> lambda=" << emitter.ray.lambda << endl;
	logfile	<< "> a=" << sample->grainlist[0]->a << endl;
	logfile << "# DATA:" << endl;
	logfile	<< "#omega\tchi\tphi\ttx\tty\ttz\tydet\tzdet\tintens" << endl;
}


void
XRayScope::writeLogEntry(void){
	// Write the state of the detector to the logfile

	int i,j,pix;
	for(j=0; j<detector.image->getHeight(); j++){
		for(i=0; i<detector.image->getWidth(); i++){
			pix = detector.image->getPixel(i,j);
			if (pix > 0){	
				logfile << omega<<"\t"<<chi<<"\t"<<phi<<"\t"						
						<< tx <<"\t"<<ty <<"\t"<<tz<<"\t"
						<< i <<"\t"<<j<<"\t"
						<< pix
						<< endl;
			}
		}
	}
}

void
XRayScope::display(void){
	// Display the parts of the detector in OpenGL

#ifdef USE_OPENGL // ONLY COMPILE THIS IF USING OPENGL

	int i;

	glLoadIdentity();

	// draw the detector:
	glColor4f(0.2,0.2,0.2,0.8);
	glBegin(GL_QUADS);
		for(i=0; i<4; i++){
			glVertex4dv(detector.points[i].Ref());
		}
	glEnd();

	//draw the emitter:
	glPushMatrix();
	glTranslatef(emitter.ray.orig[X],emitter.ray.orig[Y],emitter.ray.orig[Z]);
	glColor3f(1,1,0.5);
	//gluSphere(quadric, 0.25, 10, 10);
	glPopMatrix();
	glBegin(GL_LINES);
		glVertex4dv(emitter.ray.orig.Ref());
		glVertex4dv((emitter.ray.orig+emitter.ray.dir).Ref());
	glEnd();


	// draw the sample:
	sample->display();

#endif // end OPENGL

}

void
XRayScope::activate(void){
	// Turn on the beam - raytrace it


	// Start by tracing the beam starting at the emitter:
	debug << "Raytracing the beam..." << endl;
	trace(emitter.ray);
	
}


void
XRayScope::trace(Ray &ray){
	// intersect with everything:
	
	Isect isect;
	Isect *is;
	List<Isect> *isectlist;
	List<Vec4>  *klist;
	Ray diffray;
	Vec4 *kf;
	int i, j, k;
	int x1, y1, z1, x2, y2, z2, h, dx, dy, dz, l, m, n, xinc, yinc, zinc;
	int	err1, err2, dx2, dy2, dz2;
	int count = 0;

	// first check detector 
		//debug << "Intersecting with detector..." << endl;
		//ray.pcolor = Vec3(1,0,0);
		//detector.intersect(ray, isect);
		//ray.pcolor = Vec3(1,1,1);
	
	// now intersect with the sample
	debug << "Intersecting with sample..." << endl;
	
	//===============================================================
	//ITERATE THROUGH VOXELS AND BUILD INTERSECTION LIST...
	// this is done with a 3-D Bresenham algorithm to minimize
	// intersection testing

	// first, compute the entry and exit voxels:	
	if( !(sample->boundingVoxels(ray, x1, y1, z1, x2, y2, z2, boundmin-tx,boundmax-tx, detdist+tx))) 
		return;
	
	//now, using the Bresenham algorithm, step between them and check
	//the corresponding voxels...

	isectlist = new List<Isect>();

	i = x1; j = y1; k = z1; 
	dx=x2-x1;  dy=y2-y1;  dz=z2-z1; 

	xinc = (dx<0) ? -1 : 1; l = abs(dx);
	yinc = (dy<0) ? -1 : 1; m = abs(dy);
	zinc = (dz<0) ? -1 : 1; n = abs(dz);
	dx2  = (l<<1);
	dy2  = (m<<1);
	dz2  = (n<<1);

	if ((l>=m) && (l>=n)) {
		err1 = dy2 - l;
		err2 = dz2 - l;
		for(h=0; h<l; h++){
			if(i<sample->getXvox() && j<sample->getYvox() && k<sample->getZvox())
				sample->intersect( i,j,k,ray,isectlist );	
			if (err1 > 0) {
				j += yinc;
				err1 -= dx2;
			}
			if (err2 > 0){
				k += zinc;
				err2 -= dx2;
			}
			err1 += dy2;
			err2 += dz2;
			i += xinc;
		}
	}
	else if ((m>=l) && (m>=n)){
		err1 = dx2 - m;
		err2 = dz2 - m;
		for (h=0; h<m; h++){
			if(i<sample->getXvox() && j<sample->getYvox() && k<sample->getZvox())
				sample->intersect( i,j,k,ray,isectlist );	
			if (err1 > 0){
				i += xinc;
				err1 -= dy2;
			}
			if (err2 > 0) {
				k += zinc;
				err2 -= dy2;
			}
			err1 += dx2;
			err2 += dz2;
			j += yinc;
		}
	}
	else {
		err1 = dy2 - n;
		err2 = dx2 - n;
		for(h=0; h<n; h++){
			sample->intersect( i,j,k,ray,isectlist );	
			if (err1>0){
				j += yinc;
				err1 -= dz2;
			}
			if(err2 > 0){
				i += xinc;
				err2 -= dz2;
			}
			err1 += dy2;
			err2 += dx2;
			k += zinc;
		}
	}
	if(i<sample->getXvox() && j<sample->getYvox() && k<sample->getZvox())
		sample->intersect( i,j,k,ray,isectlist );	


	//===============================================================
	// NOW, LOOK AT EACH INTERSECTION AND CALCULATE DIFFRACTED RAYS
	// for each isect in the isectlist...
	if (isectlist->isEmpty())  return;

	isectlist->beginIteration();
	while(is = isectlist->getItem()){

		// for each sample intersection, calculate the 
		// diffracted beam(s) and send them off to the detector:
		// Diffracted beam origin is from center of voxel:

		if (is->grainid >= sample->numgrains){
			cout << "ERROR: Grainid " << is->grainid << " is out of bounds!"<< endl;
			continue;
		}
	
		// IF THERE IS DIFFRACTION, GET A LIST OF DIFFRACTED WAVEVECTORS:
		klist = new List<Vec4>();

		if(sample->grainlist[is->grainid]->getDiffraction(ray,*klist)){  
			
			diffray.orig = sample->centerOfVoxel(is->i, is->j, is->k);

			// FOR EACH DIFFRACTED RAY...
			klist->beginIteration();
			while (kf = klist->getItem()){
				diffray.pcolor = Vec3(1,1,1);
				diffray.dir = *kf;


				#ifdef USE_OPENGL // SHOW DIFFRACTED RAY
				glColor4f(1,0,1,1);
				glBegin(GL_LINES);
					glVertex4dv(diffray.orig.Ref());
					glVertex4dv((diffray.orig+diffray.dir*4).Ref());	
				glEnd();
				#endif			  // END USE_OPENGL..


				// FIRE IT OFF TO THE DETECTOR
				detector.intersect(diffray, isect);

			}
		}
		delete klist;
	}
	delete isectlist;

}

void
XRayScope::generateSample(void){
// Assign grains to each voxel in the sample
// based on a particular algorithm
	
	sample->generate(VOLUME_RANDOM);

}

int
XRayScope::watchDetectorSpot(int x0, int y0, int x1, int y1){
// Return the integrated value of a "spot" on the detector
// by integrating values over the rectangle defined by (x0,y0)..(x1,y1) 

	int i,j;	
	int value  = 0;
	int width  = detector.image->getWidth();
	int height = detector.image->getHeight();

	if(x1<x0 || y1<y0 || y0<0 || y1<0 || x0<0 || x1<0 || x0>=width || x1>=width || y0>=height || y1>=height){ 
		cerr << "watchDetector(): INCORRECT BOUNDS!"<< endl;
		return 0;
	}


	//Calculate integrated intensity...
	for (i=x0; i<x1; i++){
		for (j=y0; j<y1; j++){
			value += detector.image->getPixel(i,j);
		}
	}

	return value;

}


/////////////////////////////////////////////////////////////////////////
// MAIN:

int
main(int argc, char *argv[]){

	int i;
	float t;

	PI = acos(-1.0);				// Define PI constant from global.h	

	cout << "PI="<<PI;
	
	cout << endl;
	cout << "3D X-Ray Microscope Simulator   <kosack@andrew.cmu.edu>" << endl;
	cout << "-------------------------------------------------------" << endl;

	// PROCESS COMMAND LINE ARGUMENTS:
	if(argc==1) cout << "type: "<<argv[0]<<" -h for list of options\n"<<endl;

	for(i=1; i<argc; i++){
		if (strcmp(argv[i],"-d")==0){
			cout << "Debugging info: on" << endl;
			DEBUG = 1;
		}
		else if (strcmp(argv[i],"-t")==0){
			GRAPHICS = 0;
			cout << "Graphic Display/User interface: disabled" << endl;
		}
		else if(strcmp(argv[i],"-h")==0) {
			cout	<< "SWITCHES: " << endl
					<< "   -d    turn on debugging info (verbose)" << endl
					<< "   -t	 text mode (no graphics)" << endl
					<< "   -h    display this help message" << endl;
						
			exit(0);
		}
	}


	scope.loadSample("sample.xrs");	


	//SET UP THE MACHINE:
	scope.setDetectorDist(5.5);
	scope.setDetectorRes(500,500);
	scope.setDetectorSize(5.5,5.5);
	scope.setBeamDiv(0.001);
	scope.setProbeVolumeSize(-0.005, 0.005);
	scope.generateSample();

	scope.rotateSample(0,0,0);			
	scope.translateSample(0,0,0);

	// INITIALIZE GRAPHICS DISPLAY AND START THE USER INTERFACE:
	if (GRAPHICS){
		interface_start(argc, argv);
	}
	else {

		// here's an example scan without the GUI interface:
		// (this will run if you specify the -t option)

		cout << "Spinning in OMEGA" << endl;
		for(t=0; t<2*_PI; t+=0.0016){
			scope.rotateSample(t,0,0);
			scope.activate();
		}
		scope.writeLogHeader("final autorotate output - ignore orientation");
		scope.writeLogEntry();
	}

	return 0;
}




