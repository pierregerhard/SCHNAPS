#ifndef _STVENANTLIN_H
#define _STVENANTLIN_H

//essai
#define const_g 9.81
#define H0 0.01 


//! \brief particular flux for the StVenantLin model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void StVenantLinNumFlux(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular flux for the 2d StVenantLin model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void StVenantLinNumFlux2d(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular boundary flux for the StVenantLin model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void StVenantLinBoundaryFlux(double* x,double t,double* wL,double* vn,
			   double* flux);
//! \brief particular boundary flux for the 2d StVenantLin model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void StVenantLinBoundaryFlux2d(double* x,double t,double* wL,double* vn,
			     double* flux);
//! \brief particular init data for the StVenantLin model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void StVenantLinInitData(double* x,double* w);
//! \brief particular init data for the 2d StVenantLin model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void StVenantLinInitData2d(double* x,double* w);
//! \brief particular imposed data for the StVenantLin model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void StVenantLinImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d StVenantLin model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void StVenantLinImposedData2d(double* x,double t,double* w);

//! \brief particular flux for testing the StVenantLin model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestStVenantLinBoundaryFlux(double* x,double t,double* wL,double* vn,
			       double* flux);
//! \brief particular flux for testing the 2d StVenantLin model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestStVenantLinBoundaryFlux2d(double* x,double t,double* wL,double* vn,
				 double* flux);
//! \brief particular init data for the StVenantLin model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestStVenantLinInitData(double* x,double* w);
//! \brief particular init data for the 2d StVenantLin model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestStVenantLinInitData2d(double* x,double* w);
//! \brief particular imposed data for the StVenantLin model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestStVenantLinImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d StVenantLin model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestStVenantLinImposedData2d(double* x,double t,double* w);


#endif
