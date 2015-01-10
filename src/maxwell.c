#include "model.h"
#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#define _PI 3.1415926535897932384626
#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define _M 3;

void MaxwellNumFlux2dTM(double wL[],double wR[],double* vnorm,double* flux){
	double n1,n2,n3,r,r2,k,l,m,f1,f2,f3;
	double n1c,n2c;
	double a,b,theta,c,nx,ny,norme;
    int i,sgnn1,sgnn2;
	// double* wL=WaL;
	// double* wR=WaR;	
	
	// double flux2[3];
	
	
 // printf("%f\n",vnorm[1]/vnorm[0]);

	
	
	// if(vnorm[1]<1e-3){
	// theta=0;
	// }
	
  // norme = sqrt(vnorm[0]*vnorm[0]+ vnorm[1]*vnorm[1]);
	
	r=sqrt(vnorm[0]*vnorm[0] + vnorm[1]*vnorm[1]);
	// printf("r1 = %f\n",vnorm[2]);
	  
  // if(nx<1e-3){
	// theta= _PI*0.5 ;
  // }
  
	theta=0;
	n1=0;
	n2=0;
	
	
	if( abs(vnorm[0])>1e-8 ){
		if(vnorm[0]<0) sgnn1=-1;
		if(vnorm[0]>0) sgnn1=1;
		if(vnorm[1]<0) sgnn2=-1;
		if(vnorm[1]>0) sgnn2=1;
		
		theta = atan(vnorm[1]/vnorm[0]);
		n1=sgnn1*r*cos(theta);
		// n2=sgnn2*r*sin(theta);
		 if(sgnn1==1)
		 n2=1*r*sin(theta);
		 if(sgnn1==-1) 
		 n2=-1*r*sin(theta);
		 
		  assert(fabs(n1-vnorm[0])<1e-8);
		  assert(fabs(n2-vnorm[1])<1e-8);
		// printf("%f, %f, %f \n",n2-vnorm[1],n2,vnorm[1]);
		 // if(sgnn2==1 && sgnn1==1) printf("%f, %f, %f \n",n2+vnorm[1],n2,vnorm[1]);
		 // if(sgnn2==-1 && sgnn1==1) printf("%f, %f, %f \n",n2+vnorm[1],n2,vnorm[1]);
	}
	
		if( abs(vnorm[1])<1e-8 && vnorm[1] > 1e-8){
		if(vnorm[0]<0) sgnn1=-1;
		if(vnorm[0]>0) sgnn1=1;
		if(vnorm[1]<0) sgnn2=-1;
		if(vnorm[1]>0) sgnn2=1;
		
		theta = _PI*0.5 - atan(vnorm[0]/vnorm[1]);
		// n1=sgnn1*r*cos(theta);
		n2=sgnn2*r*sin(theta);
		if(sgnn2==1)
		 n1=1*r*cos(theta);
		 if(sgnn2==-1) 
		 n1=-1*r*cos(theta);
		 // printf("%f, %f, %f \n",n1-vnorm[0],n1,vnorm[0]);
		  assert(fabs(n1-vnorm[0])<1e-8);
		  assert(fabs(n2-vnorm[1])<1e-8);
		// printf("%f, %f, %f \n",n2-vnorm[1],n2,vnorm[1]);
		 // if(sgnn2==1 && sgnn1==1) printf("%f, %f, %f \n",n2+vnorm[1],n2,vnorm[1]);
		 // if(sgnn2==-1 && sgnn1==1) printf("%f, %f, %f \n",n2+vnorm[1],n2,vnorm[1]);
	}
	
	
	
	// if(abs(vnorm[1])>1e-8){
		// theta = _PI*0.5 - atan(vnorm[0]/vnorm[1]);
    // }
	
	// if(abs(vnorm[0])<1e-8 && abs(vnorm[1])<1e-8){
	// theta=0;
	// n1=0;
	// n2=0;
	// r=0;
	// }
	
	// if(abs(vnorm[1])>abs(vnorm[0])){
	  // theta =  _PI*0.5 - theta;
	// }
  // }

	 // n1=r*cos(theta);
	 // printf("%f\n",n2-nx);
	
	 
	 
	 // n2=r*sin(theta);
 
	 a=cos(theta);
	 b=sin(theta);
	 
	 assert(fabs(b*n2)<100);
	 assert(fabs(a*n2)<100);
	 assert(fabs(a*n1)<100);
	 assert(r<100);
			
	 // printf("%f\n",theta);
	
	
	/*
	*	RUSANOV  F(wL,wR,n) = 0.5*A1n1(wL + wR) + 0.5*A2n2(wL + wR) - (r/2)*(WR - WL)
	*
	*/	
	
	// flux[0] = -0.5*vnorm[1]*( wL[2] + wR[2] ) -0.5*r*( wR[0]-wL[0] );
	// flux[1] =  0.5*vnorm[0]*( wL[2] + wR[2] ) -0.5*r*( wR[1]-wL[1] ) ;
	// flux[2] = -0.5*vnorm[1]*( wL[0] + wR[0] )  + 0.5*vnorm[0]*( wL[1] +wR[1] )  - 0.5*r*( wR[2]-wL[2] );
	
	/*
	*	CENTERED FLUX F(wL,wR,n) = 0.5*A1n1(wL + wR) + 0.5*A2n2(wL + wR)
	*
	*/	
	
	// flux[0] = -(1./2)*vnorm[1]*( wL[2] + wR[2] );
	// flux[1] =  (1./2)*vnorm[0]*( wL[2] + wR[2] ) ;
	// flux[2] = -(1/2.)*vnorm[1]*( wL[0] + wR[0] )  + (1/2.)*vnorm[0]*( wL[1] +wR[1] ) ;
	
	
	
	/*
	*	UPWIND F(wL,wR,n) = Aini-.wR + Aini+.wL
	*
	*/	

		
	flux[0]   =- b*n2*wR[0]  + a*n2*wR[1] -  n2*wR[2] +  b*n2*wL[0]  -  a*n2*wL[1] -  n2*wL[2];
	flux[1]   =  a*n2*wR[0]   - a*n1*wR[1] + n1*wR[2]  -  a*n2*wL[0]  + a*n1*wL[1] +  n1*wL[2];
	flux[2]   =-   n2*wR[0]  +   n1*wR[1]  -      r*wR[2]   -    n2*wL[0]  +   n1*wL[1] +  r*wL[2]; 
	
   flux[0] =(1./2)*flux[0];
   flux[1] =(1./2)*flux[1];
   flux[2] =(1./2)*flux[2];
   assert(r<100);
   
	 assert(fabs( flux[0] )<100);
	 assert(fabs( flux[1] )<100);
	 assert(fabs( flux[2] )<100);
};


//Here we dont use silver-muller boundary condition...

void MaxwellBoundary2dTM(double x[3],double t,double wL[],double* vnorm,
						 double* flux){
  
  double wR[3];
 

  // double eps = 0.1;
  // double distxy =  sqrt( (x[0]-0)*(x[0]-0) + (x[1]-0)*(x[1]-0) );
  // MaxwellImposedData2dTM(x,t,wR);
  // MaxwellNumFlux2dTM(wL,wR,vnorm,flux);
  // if(distxy > 0.3 &&  distxy  <= 0.3+eps){
  // wR[0]=0;
  // wR[1]=0;
  // wR[2]=0;
  // }
   // if( distxy < 0.3 ){
  // wR[0]=0;
  // wR[1]=0;
  // wR[2]=0;
  // }
  // wR[0]=wL[0];
  // wR[1]=wL[1];
  // wR[2]=wL[2];
  // MaxwellMetalBoundary2DTM(x,t,wL,vnorm,flux);
  // else{
  MaxwellImposedData2dTM(x,t,wR);
  // }
  MaxwellNumFlux2dTM(wL,wR,vnorm,flux);
};


void MaxwellInitData2dTM(double x[3],double w[]){

  double t=0;
  MaxwellImposedData2dTM(x,t,w);

};

void MaxwellImposedData2dTM(double x[3], double t, double* w){
	
	
	
	double pi=_PI;
    // double k=2*pi;
    // double theta=pi/2;
    // double onde=cos(k*(cos(theta)*x[0]+sin(theta)*x[1]-t) + pi/2.0 );
    // W[0] = -sin(theta)*onde;
    // W[1] =  cos(theta)*onde;
    // W[2] =  onde;
	// double distxy =  sqrt( (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) +(x[2]-0)*(x[2]-0) );
  // double distz =  sqrt((x[2]-0)*(x[2]-0));	

	w[0] = 0;
	w[1] = cos(2*pi*(x[0]-t));
	w[2] = cos(2*pi*(x[0]-t));
	
	// w[0]=exp(-t/2);
	// w[1]=exp(-t/2);
	// w[2]=exp(-t/2);
	
	
	
	

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
