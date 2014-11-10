#include "shallow.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

const double transport_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double transport_v[] = {1,0,0};

const double transport_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};


void ShallowNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    transport_v[0] * vnorm[0] +
    transport_v[1] * vnorm[1] +
    transport_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];

};

void ShallowNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    transport_v2d[0] * vnorm[0] +
    transport_v2d[1] * vnorm[1] +
    transport_v2d[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];
   /* if (fabs(vnorm[2])>1e-6){ */
   /*   printf("vnds %lf %lf %lf \n",vnorm[0],vnorm[1],vnorm[2]); */
   /* } */
   // verify that 2d computations are actually
   // activated
   assert(fabs(vnorm[2])<1e-8);


};

void ShallowBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  ShallowImposedData(x,t,wR);
  ShallowNumFlux(wL,wR,vnorm,flux);
};

void ShallowBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  ShallowImposedData2d(x,t,wR);
  ShallowNumFlux2d(wL,wR,vnorm,flux);
};

void ShallowInitData(double x[3],double w[]){

  double t=0;
  ShallowImposedData(x,t,w);

};

void ShallowInitData2d(double x[3],double w[]){

  double t=0;
  ShallowImposedData2d(x,t,w);

};


void ShallowImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void ShallowImposedData2d(double x[3],double t,double w[]){

  double vx =
    transport_v2d[0] * x[0] +
    transport_v2d[1] * x[1] +
    transport_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void TestShallowBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestShallowImposedData(x,t,wR);
  ShallowNumFlux(wL,wR,vnorm,flux);
};

void TestShallowBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestShallowImposedData2d(x,t,wR);
  ShallowNumFlux2d(wL,wR,vnorm,flux);
};

void TestShallowInitData(double x[3],double w[]){

  double t=0;
  TestShallowImposedData(x,t,w);

};


void TestShallowInitData2d(double x[3],double w[]){

  double t=0;
  TestShallowImposedData2d(x,t,w);

};

void TestShallowImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestShallowImposedData2d(double x[3],double t,double w[]){

  double vx =
    transport_v2d[0] * x[0] +
    transport_v2d[1] * x[1] +
    transport_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};


