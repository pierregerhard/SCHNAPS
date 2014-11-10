#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

const double waves_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double transport_v[] = {1,0,0};

const double waves_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};


void WavesNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    waves_v[0] * vnorm[0] +
    waves_v[1] * vnorm[1] +
    waves_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];

};

void WavesNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    waves_v2d[0] * vnorm[0] +
    waves_v2d[1] * vnorm[1] +
    waves_v2d[2] * vnorm[2];

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

void WavesBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  WavesImposedData(x,t,wR);
  WavesNumFlux(wL,wR,vnorm,flux);
};

void WavesBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  WavesImposedData2d(x,t,wR);
  WavesNumFlux2d(wL,wR,vnorm,flux);
};

void WavesInitData(double x[3],double w[]){

  double t=0;
  WavesImposedData(x,t,w);

};

void WavesInitData2d(double x[3],double w[]){

  double t=0;
  WavesImposedData2d(x,t,w);

};


void WavesImposedData(double x[3],double t,double w[]){

  double vx =
    waves_v[0] * x[0] +
    waves_v[1] * x[1] +
    waves_v[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void WavesImposedData2d(double x[3],double t,double w[]){

  double vx =
    waves_v2d[0] * x[0] +
    waves_v2d[1] * x[1] +
    waves_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void TestWavesBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestWavesImposedData(x,t,wR);
  WavesNumFlux(wL,wR,vnorm,flux);
};

void TestWavesBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestWavesImposedData2d(x,t,wR);
  WavesNumFlux2d(wL,wR,vnorm,flux);
};

void TestWavesInitData(double x[3],double w[]){

  double t=0;
  TestWavesImposedData(x,t,w);

};


void TestWavesInitData2d(double x[3],double w[]){

  double t=0;
  TestWavesImposedData2d(x,t,w);

};

void TestWavesImposedData(double x[3],double t,double w[]){

  double vx =
    waves_v[0] * x[0] +
    waves_v[1] * x[1] +
    waves_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestWavesImposedData2d(double x[3],double t,double w[]){

  double vx =
    waves_v2d[0] * x[0] +
    waves_v2d[1] * x[1] +
    waves_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};
