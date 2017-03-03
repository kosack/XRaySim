/*****************************************************************
 * M  a  t  r  i  x  S  t  a  c  k  .  c  p  p                   *
 * Matrix Operations                                             *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/

#include "MatrixStack.h"

MatrixStack::MatrixStack(void){
		
	int i;

	for (i=0; i<MATRIX_MODES; i++){
		mode[i] = new MatrixNode;
		mode[i]->next = NULL;
		mode[i]->mat  = vl_1;  // ident matrix
	}

	top = mode[MODEL_MATRIX_MODE];
	curmode = MODEL_MATRIX_MODE;

}

MatrixStack::~MatrixStack(void){
	
	while(top->next != NULL){
		pop();
	}	

	delete top;

}


void
MatrixStack::loadIdentity(void){
	
	top->mat = vl_1;  // defined by SVL - 

}

void
MatrixStack::push(void){

	MatrixNode *temp = new MatrixNode;
	int i,j;

	for(i=0; i<4; i++){
		for(j=0; j<4; j++){
			temp->mat[i][j] = top->mat[i][j];
		}
	}

	temp->next = top;
	top = temp;

}


void 
MatrixStack::pop(void){

	if (top->next != NULL){
		MatrixNode *temp = top;
		top = top->next;
		delete temp;
	}
	else {
		cout << "================" << endl;
		cout << "\nMatrixStack: YOU POPPED ONCE TO MANY! Check your code." << endl;
		print();
		cout << "================" << endl;
	}
	
}

Mat4
MatrixStack::getMatrix(void){
	return top->mat;
}

void
MatrixStack::mult(Mat4 m){
	
	top->mat = top->mat * m;  // TODO: ORDER OF OPS CORRECT?

}

void
MatrixStack::print(void){
	cout << "CURRENT MATRIX: (" ;
	switch (curmode) {
		case MODEL_MATRIX_MODE: cout << "MODEL_MATRIX_MODE) "; break;
		case VIEW_MATRIX_MODE : cout << "VIEW_MATRIX_MODE) " ; break;
	}

	cout << "(stack size = " << getLength() << ") " << endl;
	cout << top->mat << endl;
}

void
MatrixStack::setMatrixMode(int m){
	if (m<0 || m>MATRIX_MODES){
		cerr << "INVALID MATRIX MODE: " << m << "  (mode not changed)" <<endl;
	}
	else {
		curmode = m;
		top = mode[m]; // set the current mode
	}
}

int
MatrixStack::getMatrixMode(void){
	return curmode;
}

int
MatrixStack::getLength(void){
	MatrixNode *temp;
	int i=0;

	temp = top;
	while (temp != NULL){
		temp = temp->next;
		i++;
	}

	return i;
}

void
MatrixStack::eulerRot(double phi, double theta, double psi, bool inverse = false){
	
	double cosphi = cos(phi), sinphi = sin(phi);
	double cospsi = cos(psi), sinpsi = sin(psi);
	double costheta = cos(theta), sintheta=sin(theta);

	Mat4 eulermat = vl_1;

	eulermat[0][0] = cospsi*cosphi - costheta*sinphi*sinpsi;
	eulermat[1][0] =-sinpsi*cosphi - costheta*sinphi*cospsi;
	eulermat[2][0] = sintheta*sinphi;  
	eulermat[0][1] = cospsi*sinphi + costheta*cosphi*sinpsi;
	eulermat[1][1] =-sinpsi*sinphi + costheta*cosphi*cospsi;
	eulermat[2][1] =-sintheta*cosphi;
	eulermat[0][2] = sinpsi*sintheta;
	eulermat[1][2] = cospsi*sintheta;
	eulermat[2][2] = costheta;


	if (inverse) eulermat = inv(eulermat); 

	mult(eulermat);

}

void
MatrixStack::translate(double tx, double ty, double tz){

	Mat4 trans = vl_1;
	trans[0][3] = tx;
	trans[1][3] = ty;
	trans[2][3] = tz;
	mult(trans);

}