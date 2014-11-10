#ifndef _EULER_H
#define _EULER_H
#include <math.h>
#include <stdio.h>

void EulerNumFlux(double wL[],double wR[],double vn[3],double* flux);

void EulerNumFlux2d(double wL[],double wR[],double vn[3],double* flux);

void EulerBoundaryFlux(double* x,double t,double* wL,double* vn,
			   double* flux);

void EulerBoundaryFlux2d(double* x,double t,double* wL,double* vn,
			     double* flux);

void EulerInitData(double* x,double* w);

void EulerInitData2d(double* x,double* w);

void EulerImposedData(double* x,double t,double* w);

void EulerImposedData2d(double* x,double t,double* w);

void TestEulerBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux);

void TestEulerBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux);

void TestEulerInitData(double x[3],double w[]) ;

void TestEulerInitData2d(double x[3],double w[]);

void TestEulerImposedData(double x[3],double t,double w[]);

void TestEulerImposedData2d(double x[3],double t,double w[]);

#endif
