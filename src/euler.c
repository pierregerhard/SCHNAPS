#include "euler.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

const double Euler_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double Euler_v[] = {1,0,0};

const double Euler_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};

double WL [] = {0, 0, 0, 0} ;
double WR [] = {0, 0, 0, 0} ;

void EulerNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    Euler_v[0] * vnorm[0] +
    Euler_v[1] * vnorm[1] +
    Euler_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   //flux[0] = vnp * wL[0] + vnm * wR[0];
   /*
   flux[0] = vnorm[0] ; 
   flux[1] = vnorm[0];
   flux[2] = vnorm[0] ;
   flux[3] = vnorm[0] ;
   */
   flux[0] = vnp * wL[0] + vnm * wR[0];
   flux[1] = vnp * wL[1] + vnm * wR[1];
   flux[2] = vnp * wL[2] + vnm * wR[2];
   flux[3] = vnp * wL[3] + vnm * wR[3];

};

void EulerNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){
  WL[0] = wL[0] ;
  WL[1] = wL[1] ;
  WL[2] = wL[2] ;
  WL[3] = wL[3] ;
  WR[0] = wR[0] ;
  WR[1] = wR[1] ;
  WR[2] = wR[2] ;
  WR[3] = wR[3] ;
  double lambda_max = sqrt(wL[1]*wL[1]/wL[0]/wL[0]) + sqrt (GAMMA*wL[0]) ;
   
   flux[0] = (wL[1]*vnorm[0] + wL[2]*vnorm[1] 
   		   + wR[1]*vnorm[0] + wR[2] *vnorm[1])/2 
		   - lambda_max*(wR[0]-wL[0])  ; 

   flux[1] = ((wL[1]*wL[1]/wL[0] + wL[0]*wL[0])*vnorm[0] + wL[1]*wL[2]/wL[0]*vnorm[1] 
   		   +  (wR[1]*wR[1]/wR[0] + wR[0]*wR[0])*vnorm[0] + wR[1]*wR[2]/wR[0]*vnorm[1])/2 
		   - lambda_max*(wR[1]-wL[1])   ;

   flux[2] = (wL[1]*wL[2]/wL[0]*vnorm[0] + (wL[2]*wL[2]/wL[0] + wL[0]*wL[0])*vnorm[1] 
   		   + wR[1]*wR[2]/wR[0]*vnorm[0] + (wR[2]*wR[2]/wR[0] + wR[0]*wR[0])*vnorm[1])/2 
		   - lambda_max*(wR[2]-wL[2])   ;

   flux[3] = ((wL[1]*wL[3]/wL[0] + wL[0]*wL[0]*wL[1]/wL[0])*vnorm[0] + (wL[2]*wL[3]/wL[0] + wL[0]*wL[0]*wL[1]/wL[0])*vnorm[1]
   		   + (wR[1]*wR[3]/wR[0] + wR[0]*wR[0]*wR[1]/wR[0])*vnorm[0] + (wR[2]*wR[3]/wR[0] + wR[0]*wR[0]*wR[1]/wR[0])*vnorm[1])/2
		   - lambda_max*(wR[3]-wL[3])   ;
};

void EulerBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[4];
  EulerImposedData(x,t,wR);
  EulerNumFlux(wL,wR,vnorm,flux);
};

void EulerBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[4];
  EulerImposedData2d(x,t,wR);
  EulerNumFlux2d(wL,wR,vnorm,flux);
};

void EulerInitData(double x[3],double w[]){

  double t=0;
  EulerImposedData(x,t,w);

};

void EulerInitData2d(double x[3],double w[]){

  double t=0;
  EulerImposedData2d(x,t,w);

};

void riemann (double * wL, double * wR, double xi, double * w) {
	doublereal rL,uL,vL,pL,rR,uR,vR,pR,r,u,v,p,y,um ;
	rL=wL[0] ;
	uL=wL[1]/rL ;
	vL=wL[2]/rL ;
	pL=(wL[3]-0.5*(wL[0]*wL[1]+wL[0]*wL[2]))*(GAMMA-1) ;

	rR=wR[0] ;
	uR=wR[1]/rR ;
	vR=wR[2]/rR ;
	pR=(wR[3]-0.5*(wR[0]*wR[1]+wR[0]*wR[2]))*(GAMMA-1) ;

	//riemann77_(&rL,&uL,&pL,&rR,&uR,&pR,&xi,&r,&u,&p,&um) ;

	if (um > xi) 
		v=vL ;
	else 
		v=vR ;
	w[0]=r ;
	w[1]=r*u ; 
	w[2]=r*v ;
	w[3]=p/(GAMMA-1)+0.5*r*(u*u+v*v) ;
}

void EulerImposedData(double x[3],double t,double w[]){
  double xi = x[0]/t ;
  riemann (WL, WR, xi, w) ;
};

void EulerImposedData2d(double x[3],double t,double w[]){
  double xi = x[0]/t ;
  riemann (WL, WR, xi, w) ;
};

void TestEulerBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestEulerImposedData(x,t,wR);
  EulerNumFlux(wL,wR,vnorm,flux);
};

void TestEulerBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[4];
  TestEulerImposedData2d(x,t,wR);
  EulerNumFlux2d(wL,wR,vnorm,flux);
};

void TestEulerInitData(double x[3],double w[]){

  double t=0;
  TestEulerImposedData(x,t,w);

};


void TestEulerInitData2d(double x[3],double w[]){

  double t=0;
  TestEulerImposedData2d(x,t,w);

};

void TestEulerImposedData(double x[3],double t,double w[]){

  double vx =
    Euler_v[0] * x[0] +
    Euler_v[1] * x[1] +
    Euler_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestEulerImposedData2d(double x[3],double t,double w[]){

  double vx =
    Euler_v2d[0] * x[0] +
    Euler_v2d[1] * x[1] +
    Euler_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};
