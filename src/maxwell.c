#include "model.h"
#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#define _PI 3.1415926535897932384626

void MaxwellNumFlux2dTM(double wL[],double wR[],double* vnorm,double* flux){
	double n1,n2,r;
	n1=vnorm[0];
	n2=vnorm[1];
	r = sqrt(n1*n1 +n2*n2);
	
	/*
	* Elements of maxtrix Ai.ni(-)
	*/
	
	double 	am11,am12,am13,
			am21,am22,am23,
			am31,am32,am33;
	
	am11 = -n2*n2/r;
	am12 =  n1*n2/r;
	am13 = -n2;
	am21 =  n1*n2;
	am22 = -n1*n1;
	am23 =  n1;
	am31 = -n2;
	am32 =  n1;
	am33 = -r;
	
	/*
	* Elements of matrix Ai.ni(+)
	*/
	double 	ap11,ap12,ap13,
			ap21,ap22,ap23,
			ap31,ap32,ap33;
	
	ap11 = -am11;
	ap12 = -am12;
	ap13 =  am13;
	ap21 = -am21;
	ap22 = -am22;
	ap23 =  am23;
	ap31 =  am31;
	ap32 =  am32;
	ap33 = -am33;
	
	/*
	*	UPWIND FLUX
	*/
	flux[0] = 0.5*(ap11*wL[0] + ap12*wL[1] + ap13*wL[2])
			+ 0.5*(am11*wR[0] + am12*wR[1] + am13*wL[2]);
	flux[1] = 0.5*(ap21*wL[0] + ap22*wL[1] + ap23*wL[2])
			+ 0.5*(am21*wR[0] + am22*wR[1] + am23*wL[2]);
	flux[2] = 0.5*(ap31*wL[0] + ap32*wL[1] + ap33*wL[2]) 
			+ 0.5*(am31*wR[0] + am32*wR[1] + am33*wL[2]);
	
	/*
	*	CENTERED FLUX F(wL,wR,n) = 0.5*A1n1(wL + wR) + 0.5*A2n2(wL + wR)
	*
	*	flux[0] = 0.5*(-n2*(wL[3] + wR[3]));
	*	flux[1] = 0.5*(n1*(wL[3] + wR[3]));
	*	flux[2] = 0.5*( -n2*(wL[1] +wR[1]) + n1*(wL[2] +wR[2]) );
	*/
	 //assert(fabs(vnorm[3])<1e-8);
};

void MaxwellMetalBoundary2DTM(double x[3],double t,double wL[],double* vnorm,
						 double* flux){
  flux[0]=0;
  flux[1]=0;
  flux[2]=0;
};

//Here we dont use silver-muller boundary condition...

void MaxwellBoundary2dTM(double x[3],double t,double wL[],double* vnorm,
						 double* flux){
  double wR[2];
  MaxwellImposedData2dTM(x,t,wR);
  MaxwellNumFlux2dTM(wL,wR,vnorm,flux);
};


void MaxwellInitData2dTM(double x[3],double w[]){

  double t=0;
  MaxwellImposedData2dTM(x,t,w);

};

void MaxwellImposedData2dTM(double x[3], double t, double* W){
	
	double pi=_PI;
    double k=2*pi;
    double theta=15;
    double onde=cos(k*(cos(theta)*x[0]+sin(theta)*x[1]-t) + pi/2.0 );
    W[0] = -sin(theta)*onde;
    W[1] =  cos(theta)*onde;
    W[2] =  onde;
	/* Constant solution
	double cst = 4;
	W[0] = cst;
    W[1] = cst;
    W[2] = cst;
	*/
	
	/*Cavity periodic function in unit square. Mode m,n
	double c=1;
	int n = 1;
	int m = 1;
	int omega = 1;
	double u = m*pi*x[0];
	double v = n*pi*x[1]
	double w = omega*t
	W[0] = -c*c*(n*pi/omega)*cos(u)*sin(v)*sin(w);
    W[1] =  c*c*(m*pi/omega)*sin(u)*cos(v)*sin(w);
    W[2] =  cos(u)*cos(v)*cos(w);
	*/
	
};
