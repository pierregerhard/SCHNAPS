#ifndef _EULER_H
#define _EULER_H
#include "f2c.h"
#include "gasdyn77.h"

#define GAMMA 2

//! \brief particular flux for the Euler model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void EulerNumFlux(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular flux for the 2d Euler model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void EulerNumFlux2d(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular boundary flux for the Euler model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void EulerBoundaryFlux(double* x,double t,double* wL,double* vn,
			   double* flux);
//! \brief particular boundary flux for the 2d Euler model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void EulerBoundaryFlux2d(double* x,double t,double* wL,double* vn,
			     double* flux);
//! \brief particular init data for the Euler model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void EulerInitData(double* x,double* w);
//! \brief particular init data for the 2d Euler model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void EulerInitData2d(double* x,double* w);
//! \brief particular imposed data for the Euler model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void EulerImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d Euler model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void EulerImposedData2d(double* x,double t,double* w);

//! \brief particular flux for testing the Euler model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestEulerBoundaryFlux(double* x,double t,double* wL,double* vn,
			       double* flux);
//! \brief particular flux for testing the 2d Euler model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestEulerBoundaryFlux2d(double* x,double t,double* wL,double* vn,
				 double* flux);
//! \brief particular init data for the Euler model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestEulerInitData(double* x,double* w);
//! \brief particular init data for the 2d Euler model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestEulerInitData2d(double* x,double* w);
//! \brief particular imposed data for the Euler model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestEulerImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d Euler model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestEulerImposedData2d(double* x,double t,double* w);


void riemann (double * wL, double * wR, double xi, double * w) ;
#endif
