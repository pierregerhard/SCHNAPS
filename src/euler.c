#include "euler.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

const double euler_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double euler_v[] = {1,0,0};

const double euler_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};


void EulerNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    euler_v[0] * vnorm[0] +
    euler_v[1] * vnorm[1] +
    euler_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];

};

void EulerNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    euler_v2d[0] * vnorm[0] +
    euler_v2d[1] * vnorm[1] +
    euler_v2d[2] * vnorm[2];

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

void EulerBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  eulerImposedData(x,t,wR);
  eulerNumFlux(wL,wR,vnorm,flux);
};

void EulerBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  eulerImposedData2d(x,t,wR);
  eulerNumFlux2d(wL,wR,vnorm,flux);
};

void EulerInitData(double x[3],double w[]){

  double t=0;
  eulerImposedData(x,t,w);

};

void EulerInitData2d(double x[3],double w[]){

  double t=0;
  eulerImposedData2d(x,t,w);

};


void EulerImposedData(double x[3],double t,double w[]){

  double vx =
    euler_v[0] * x[0] +
    euler_v[1] * x[1] +
    euler_v[2] * x[2];

  double xx = vx - t;

  //w[0]=cos(xx);
  w[0] = 5. ;
};

void EulerImposedData2d(double x[3],double t,double w[]){

  double vx =
    euler_v2d[0] * x[0] +
    euler_v2d[1] * x[1] +
    euler_v2d[2] * x[2];

  double xx = vx - t;

  //w[0]=cos(xx);
  w[0] = 5. ;
};

void TestEulerBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TesteulerImposedData(x,t,wR);
  eulerNumFlux(wL,wR,vnorm,flux);
};

void TestEulerBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TesteulerImposedData2d(x,t,wR);
  eulerNumFlux2d(wL,wR,vnorm,flux);
};

void TestEulerInitData(double x[3],double w[]){

  double t=0;
  TesteulerImposedData(x,t,w);

};


void TestEulerInitData2d(double x[3],double w[]){

  double t=0;
  TesteulerImposedData2d(x,t,w);

};

void TestEulerImposedData(double x[3],double t,double w[]){

  double vx =
    euler_v[0] * x[0] +
    euler_v[1] * x[1] +
    euler_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestEulerImposedData2d(double x[3],double t,double w[]){

  double vx =
    euler_v2d[0] * x[0] +
    euler_v2d[1] * x[1] +
    euler_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};
