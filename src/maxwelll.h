#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <model.h>

void MaxwellNumFlux2dTM(double wL[],double wR[],double* vnorm,double* flux);
void MaxwellMetalBoundary2DTM(double x[3],double t,double wL[],double* vnorm,double* flux);
void MaxwellBoundary2dTM(double x[3],double t,double wL[],double* vnorm, double* flux);
void MaxwellInitData2dTM(double x[3],double w[]);
void MaxwellImposedData2dTM(double x[3], double t, double* W);