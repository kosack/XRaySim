/*****************************************************************
 * M  a  t  r  i  x  S  t  a  c  k  .  c  p  p                   *
 * Matrix Operations                                             *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/
#ifndef __MATRIXSTACK_H__
#define __MATRIXSTACK_H__

#include <iostream.h>
#include "libsvl/SVL.h"
#include "global.h"


struct MatrixNode {
	// Part of MatrixStack
	Mat4 mat;
	MatrixNode *next;
};

class MatrixStack {
	// Viewing and model matrices

  public:
	MatrixStack(void);
	~MatrixStack(void);
	void loadIdentity(void);				// Load the identity matrix
	void eulerRot(double,double,double,bool);// rotate with euler angles
	void translate(double,double,double);	// mult by translation matrix
	void push(void);						// Push copy of matrix onto stack
	void pop(void);							// Pop off top matrix on stack
	Mat4 getMatrix(void);					// get the current matrix
	void mult(Mat4 m);						// multiply current matrix by m
	void print(void);						// print the current matrix
	void setMatrixMode(int);				// set the current matrix mode
	int  getMatrixMode(void);				// get current matrix mode
	int	 getLength(void);					// get the length of the stack

  private:
	int curmode;							// Current matrix mode
	MatrixNode *mode[MATRIX_MODES];			// pointers to each stack mode...
	MatrixNode *top;						// stack of MatrixNodes

};

#endif //__MATRIXSTACK_H__