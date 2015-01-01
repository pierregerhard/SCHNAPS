#include "stvenantfp.h"
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


void StvenantFpNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    transport_v[0] * vnorm[0] +
    transport_v[1] * vnorm[1] +
    transport_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];

};

void StvenantFpNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){
  
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

void StvenantFpBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  StvenantFpImposedData(x,t,wR);
  StvenantFpNumFlux(wL,wR,vnorm,flux);
};

void StvenantFpBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  StvenantFpImposedData2d(x,t,wR);
  StvenantFpNumFlux2d(wL,wR,vnorm,flux);
};

void StvenantFpInitData(double x[3],double w[]){

  double t=0;
  StvenantFpImposedData(x,t,w);
};

void StvenantFpInitData2d(double x[3],double w[]){

  double t=0;
  StvenantFpImposedData2d(x,t,w);

};


void StvenantFpImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void StvenantFpImposedData2d(double x[3],double t,double w[]){

  double vx =
    transport_v2d[0] * x[0] +
    transport_v2d[1] * x[1] +
    transport_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void TestStvenantFpBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestStvenantFpImposedData(x,t,wR);
  StvenantFpNumFlux(wL,wR,vnorm,flux);
};

void TestStvenantFpBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestStvenantFpImposedData2d(x,t,wR);
  StvenantFpNumFlux2d(wL,wR,vnorm,flux);
};

void TestStvenantFpInitData(double x[3],double w[]){

  double t=0;
  TestStvenantFpImposedData(x,t,w);

};


void TestStvenantFpInitData2d(double x[3],double w[]){

  double t=0;
  TestStvenantFpImposedData2d(x,t,w);

};

void TestStvenantFpImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestStvenantFpImposedData2d(double x[3],double t,double w[]){

  double vx =
    transport_v2d[0] * x[0] +
    transport_v2d[1] * x[1] +
    transport_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};
