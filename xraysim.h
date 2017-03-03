/*****************************************************************
 * x  r  a  y  s  i  m  .  h                                     *
 *                                                               *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/

#ifndef __XRAYSIM_H__
#define __XRAYSIM_H__

/////////////////////////////////////////////////////////////////////////
// INCLUDES:

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "libsvl/SVL.h"
#include "global.h"
#include "interface.h"
#include "MatrixStack.h"
#include "list.h"

class Grain;
class Ray;
class Isect;
class List<Isect>;


/////////////////////////////////////////////////////////////////////////
// 3DXDM PARTS:

class XRayScope;

class Grain {
	// cubic crystal grain

  public:
	Grain(Mat4 rm, double ta, GrainType gt, XRayScope *parent);
	bool getDiffraction( const Ray &r, List<Vec4> &klist );
	double		a;									// lattice constant
	Vec3		pcolor;								// pseudocolor for display

  private:
	Mat4		orientation;						// Rotation matrix;
	GrainType	type;								// fcc, bcc, sc, etc...
	List<Vec4>	glist;								// list of g's
	bool		even(int i);						// support func
	XRayScope   *scope;								// scope that this grain is in

};

class Ray {
  public:
	Ray(){pcolor=Vec3(1,1,1); orig=Vec4(0,0,0,1); dir=Vec4(1,0,0,0); lambda=1.54e-8;}
	Vec4	point(double t) const {return orig+t*dir;}
	Vec4	orig;					// starting point
	Vec4	dir;					// direction
	double	intensity;				// intensity of ray
	double	lambda;					// wavelength
	Vec3	pcolor;					// pseudocolor
};

class Isect {
	// Intersection of ray and object
  public:
	double t;						// ray parameter of intersection
	bool enter;						// entering or leaving?
	Voxel  grainid;					// intersected grain
	int x,y;
	int i,j,k;						// which voxel (if voxel)
};

class Emitter {
	// An X-ray source
  public:
	Ray ray;						// ray to fire

  private:
};


class Image {
	// Detector image buffer

  public:
	Image(int w, int h);
	~Image(){ delete [] data; }
	int getWidth(void){ return width;}		  
	int getHeight(void){ return height;}
	void clear(void);							// clear the image to black
	void setPixel(int x, int y);				// Add data to pixel at x,y 
	int getPixel(int x, int y);					// get the color of pixel at x,y
	bool write(char *file);						// save image to file
	void display(void);							// OpenGL display

  private:
	int width, height;
	int *data;

	
};

class Intersectable {
  public:
	virtual bool intersect(const Ray &ray, Isect &isect){return false;};
};



class Detector : public Intersectable {
// An X-ray detector

	friend class XRayScope;
	
  public:
	Detector(){xres=25;yres=25; width=1; height=1;image = new Image(xres,yres);}
	~Detector(){delete image;}
	void setRes(int x, int y){ delete image;xres=x;yres=y;
				image=new Image(abs(xres),abs(yres)); };
	void setSize(double w, double h){ width=w; height=h;}
	double getWidth(void){return width;}
	double getHeight(void){return height;}
	int getXRes(){return xres;}
	int getYRes(){return yres;}
	bool intersect(const Ray &ray, Isect &isect);
	Image *image;


  protected:
  	Vec4 points[4];					// Upperleft, upright, lowright, lowleft
	Vec4 nor;						// normal to plane

  private:
	double width, height;
	int xres, yres;

};


enum VolumeType {
	VOLUME_CUBE,
	VOLUME_RANDOM
};



class Volume {

  public:

	Volume(double txlen, double tylen, double tzlen, double tvoxsize);
	~Volume(void);
	
	bool	intersect(int i, int j, int k, const Ray &ray, List<Isect> *isectlist);
	void	load(char *file);
	void	save(char *file);					
	void	display(void);						// OpenGL display
	void	generate(VolumeType type);
	Voxel	voxelAt(int, int, int);				// return grainid of voxel at i,j,k
	bool	voxelIndexOf(Vec4,int&,int&,int&);	// get indices of voxel of point p
	void	setVoxel(int, int, int, Voxel);		// set grainid of voxel
	int		getXvox(){return xvox;}			
	int		getYvox(){return yvox;}
	int		getZvox(){return zvox;}
	Vec4	centerOfVoxel(int,int,int);
	bool	boundingVoxels(const Ray&,int&,int&,int&,int&,int&,int&,double,double,double);
	int		resize(double txlen, double tylen, double tzlen, double tvoxsize);

	Grain		*grainlist[5000];				// array of grains in sample
	int			numgrains;						// number of grains loaded
	int			voxel_disp;

  private:



	Voxel		*voxel;					// voxel array
	double		xlen, ylen, zlen;		// dimensions of volume box
	double		voxsize;				// size of voxels to fill the box
	double		posx, posy, posz;		// position of  voxel (0,0,0)
	int			totalvox;				// total number of voxels in box
	int			xvox, yvox, zvox;		// index of max voxel in each dimension
	int			readDataEntry(ifstream infile, char *key, double &value);


};

/////////////////////////////////////////////////////////////////////////
// THE X-RAY MICROSCOPE ITSELF:

class XRayScope {

  public:
	
	XRayScope();
	~XRayScope();

	void loadSample(char *file);				// load a sample datafile
	void generateSample();						// Make a test sample
	void rotateSample(double,double,double);	// Euler rotation of sample
	void translateSample(double,double,double);	// Sample position in beam
	void activate();							// raytrace the beam
	void display();								// OpenGL support
	void setDetectorDist(double);				// set detector/emitter distance
	void setDetectorSize(double w, double h);	// set width/height of detector
	void setDetectorRes(int x, int y);			// set resolution of detector
	void setBeamDiv(double bd){beamdiv=bd;}		// get beam divergence
	void setProbeVolumeSize(double,double);		// set probe volume min,max
	double getBeamDiv(){return beamdiv;}	    // get beam divergence
	double getDetectorDist(){return detdist;}	// get detector distance
	double getProbeVolumeMin(){return boundmin;}
	double getProbeVolumeMax(){return boundmax;}
	double getOmega(){return omega;}
	double getChi(){return chi;}
	double getPhi(){return phi;}
	void setLogFilename(const char *filename);	// set current log file 
	void writeLogHeader(char*);					// write log header
	void writeLogEntry();						// write log entry
	int watchDetectorSpot(int,int,int,int );    // return avg intensity of spot

	void setVoxelDisp(int d){sample->voxel_disp = d;}


	Detector	detector;						// 2-d x-ray detector


  private:

	void redraw(void);							// redraw the apparatus
	void trace(Ray &ray);						// intersect with everything
	double		omega,chi,phi,tx,ty,tz;			// sample position (rot and trans)
	double		detdist;						// detector/beam distance from sample
	MatrixStack matrix;							// model transformation matrix
	Emitter		emitter;						// x-ray source
	Volume		*sample;						// sample volume (voxel octree)
	double		beamdiv;						// beam divergence
	ofstream	logfile;						// output logfile
	char		logfilename[255];				// log filename
	double		boundmin;						// minimum for sampling boundary
	double		boundmax;						// maximum for sampling boundary


};


#endif //__XRAYSIM_H__