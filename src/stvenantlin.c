#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

const double stvenantlin_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double stvenantlin_v[] = {1,0,0};

const double stvenantlin_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};


void StVenantLinNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    stvenantlin_v[0] * vnorm[0] +
    stvenantlin_v[1] * vnorm[1] +
    stvenantlin_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];

};

void StVenantLinNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    stvenantlin_v2d[0] * vnorm[0] +
    stvenantlin_v2d[1] * vnorm[1] +
    stvenantlin_v2d[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

  // flux[0] = vnp * wL[0] + vnm * wR[0];
   flux[0]=0;
   flux[1]=0;
   flux[2]=0;
   /* if (fabs(vnorm[2])>1e-6){ */
   /*   printf("vnds %lf %lf %lf \n",vnorm[0],vnorm[1],vnorm[2]); */
   /* } */
   // verify that 2d computations are actually
   // activated
   assert(fabs(vnorm[2])<1e-8);


};

void StVenantLinBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  StVenantLinImposedData(x,t,wR);
  StVenantLinNumFlux(wL,wR,vnorm,flux);
};

void StVenantLinBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  StVenantLinImposedData2d(x,t,wR);
  StVenantLinNumFlux2d(wL,wR,vnorm,flux);
};

void StVenantLinInitData(double x[3],double w[]){

  double t=0;
  StVenantLinImposedData(x,t,w);

};

void StVenantLinInitData2d(double x[3],double w[]){

  double t=0;
  StVenantLinImposedData2d(x,t,w);

};


void StVenantLinImposedData(double x[3],double t,double w[]){

  double vx =
    stvenantlin_v[0] * x[0] +
    stvenantlin_v[1] * x[1] +
    stvenantlin_v[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
 };

void StVenantLinImposedData2d(double x[3],double t,double w[]){

  double vx =
    stvenantlin_v2d[0] * x[0] +
    stvenantlin_v2d[1] * x[1] +
    stvenantlin_v2d[2] * x[2];

  double xx = vx - t;

 // w[0]=cos(xx);
 w[0]=1.;
 w[1]=1.;
 w[2]=1.;
};

void TestStVenantLinBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestStVenantLinImposedData(x,t,wR);
  StVenantLinNumFlux(wL,wR,vnorm,flux);
};

void TestStVenantLinBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestStVenantLinImposedData2d(x,t,wR);
  StVenantLinNumFlux2d(wL,wR,vnorm,flux);
};

void TestStVenantLinInitData(double x[3],double w[]){

  double t=0;
  TestStVenantLinImposedData(x,t,w);

};


void TestStVenantLinInitData2d(double x[3],double w[]){

  double t=0;
  TestStVenantLinImposedData2d(x,t,w);

};

void TestStVenantLinImposedData(double x[3],double t,double w[]){

  double vx =
    stvenantlin_v[0] * x[0] +
    stvenantlin_v[1] * x[1] +
    stvenantlin_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestStVenantLinImposedData2d(double x[3],double t,double w[]){

  double vx =
    stvenantlin_v2d[0] * x[0] +
    stvenantlin_v2d[1] * x[1] +
    stvenantlin_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};
