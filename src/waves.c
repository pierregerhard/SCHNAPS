#include "waves.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)
#define PI  (3.14159265359)

const double waves_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double transport_v[] = {1,0,0};

const double waves_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};

// Initialisation de la variable globale pour stocker n1 et n2
double waves_norm[] = {
    0,
    0,
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
  
    double n1,n2,s;
 
    n1=vnorm[0];
    n2=vnorm[1];

    // remplit waves_norm par n1 et n2
    waves_norm[0]=vnorm[0];
    waves_norm[1]=vnorm[1];
    s=sqrt(n1*n1+n2*n2);
     
    // Matrix A.n


    double A[3][3];
    A[0][0]=0;
    A[0][1]=-n1;
    A[0][2]=-n2;
    A[1][0]=-n1;
    A[1][1]=0;
    A[1][2]=0;
    A[2][0]=-n2;
    A[2][1]=0;
    A[2][2]=0;

    //Rusanov numerical flux

   
    flux[0]=(A[0][0]*(wL[0]+wR[0])+A[0][1]*(wL[1]+wR[1])+A[0][2]*(wL[2]+wR[2]))/2 - s/2*(wR[0]-wL[0]);
    flux[1]=(A[1][0]*(wL[0]+wR[0])+A[1][1]*(wL[1]+wR[1])+A[1][2]*(wL[2]+wR[2]))/2 - s/2*(wR[1]-wL[1]);
    flux[2]=(A[2][0]*(wL[0]+wR[0])+A[2][1]*(wL[1]+wR[1])+A[2][2]*(wL[2]+wR[2]))/2 - s/2*(wR[2]-wL[2]);
   

    // Matrix A.n(-)

    double AM[3][3];
    AM[0][0]=s/2;
    AM[0][1]=n1/2;
    AM[0][2]=n2/2;
    AM[1][0]=n1/2;
    AM[1][1]=n1*n1/(2*s);
    AM[1][2]=n1*n2/(2*s);
    AM[2][0]=n2/2;
    AM[2][1]=n1*n2/(2*s);
    AM[2][2]=n2*n2/(2*s);

    // Matrix A.n(+)

    double AP[3][3];
    AP[0][0]=-s/2;
    AP[0][1]=n1/2;
    AP[0][2]=n2/2;
    AP[1][0]=n1/2;
    AP[1][1]=-n1*n1/(2*s);
    AP[1][2]=-n1*n2/(2*s);
    AP[2][0]=n2/2;
    AP[2][1]=-n1*n2/(2*s);
    AP[2][2]=-n2*n2/(2*s);


// Godunov numerical flux

/*
    flux[0]=AM[0][0]*wR[0]+AM[0][1]*wR[1]+AM[0][2]*wR[2]+AP[0][0]*wL[0]+AP[0][1]*wL[1]+AP[0][2]*wL[2];
    flux[1]=AM[1][0]*wR[0]+AM[1][1]*wR[1]+AM[1][2]*wR[2]+AP[1][0]*wL[0]+AP[1][1]*wL[1]+AP[1][2]*wL[2];
    flux[2]=AM[2][0]*wR[0]+AM[2][1]*wR[1]+AM[2][2]*wR[2]+AP[2][0]*wL[0]+AP[2][1]*wL[1]+AP[2][2]*wL[2];

    assert(fabs(vnorm[2])<1e-8);
*/

};

void WavesBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[3];
  WavesImposedData(x,t,wR);
  WavesNumFlux(wL,wR,vnorm,flux);
};

void WavesBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[3];
  WavesImposedData2d(x,t,wR);
  if(x[0]*x[0]+x[1]*x[1]==1)
  {
      wR[0]=wL[0];
      wR[1]=wL[1]-2*(wL[1]*vnorm[0]+wL[2]*vnorm[1])*vnorm[0];
      wR[2]=wL[2]-2*(wL[1]*vnorm[0]+wL[2]*vnorm[1])*vnorm[1];
  }
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

// Onde vers x
/*
 w[0]=(x[0]-t);
 w[1]=-(x[0]-t);
 w[2]=0;//cos(4*x[1]-4*t);
*/

// Onde en biais

 w[0]=cos(4*(x[0]-x[1])/sqrt(2)+4*t);
 w[1]=cos(4*(x[0]-x[1])/sqrt(2)+4*t)/sqrt(2);
 w[2]=-cos(4*(x[0]-x[1])/sqrt(2)+4*t)/sqrt(2);


/*
 w[0]=cos(x[0]+x[1]-t)/sqrt(2);
 w[1]=-waves_norm[0]*cos(x[0]+x[1]-t)/sqrt(2);
 w[2]=-waves_norm[1]*cos(x[0]+x[1]-t)/sqrt(2);
*/

/*
 w[0]=(3*(x[0]-x[1])/sqrt(2)+3*t)+(3*(x[0]-x[1])/sqrt(2)-3*t);
 w[1]=(3*(x[0]-x[1])/sqrt(2)+3*t)/sqrt(2)+(3*(x[0]-x[1])/sqrt(2)-3*t)/sqrt(2);
 w[2]=-(3*(x[0]-x[1])/sqrt(2)+3*t)/sqrt(2)-(3*(x[0]-x[1])/sqrt(2)-3*t)/sqrt(2);
*/

/*
 w[0]=1;
 w[1]=1;
 w[2]=1;

*/

};

void TestWavesBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[3];
  TestWavesImposedData(x,t,wR);
  WavesNumFlux(wL,wR,vnorm,flux);
};

void TestWavesBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[3];
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

  w[0]=0;
  w[1]=(-1*waves_norm[1]/waves_norm[0])*cos(waves_norm[0]*x[0]+waves_norm[1]*x[1]-t);
  w[2]=cos(waves_norm[0]*x[0]+waves_norm[1]*x[1]-t);




};
