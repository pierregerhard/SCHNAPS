#ifndef _MODEL_H
#define _MODEL_H

//! \brief a unified framework for all physical models
typedef struct Model{
  //! number of conservative variables
  int m;
  //! \brief a pointer to the numflux function
  //! \param[in] wL,wR : left and right states
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*NumFlux)(double wL[],double wR[],double vn[3],double flux[]);

  //! \brief a pointer to the boundary flux function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] wL : left state
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*BoundaryFlux)(double x[3],double t,double wL[],double vn[3],double flux[]);

  //! \brief a pointer to the init data function
  // !\param[in] x : space position
  //! \param[out] w : init state at point x
  void (*InitData)(double x[3],double w[]);

  //! \brief a pointer to the imposed data function
  //!\param[in] x,t : space and time position
  //! \param[out] w : imposed state at point x and time t
  void (*ImposedData)(double x[3],double t,double w[]);

} Model;



//! \brief particular flux for the 2d stvenant model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vecto
//! \param[out] flux : the flux
void SvenantFpNumFlux2d(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular boundary flux for the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void StvenantFpBoundaryFlux(double* x,double t,double* wL,double* vn,double *flux);
//! \brief particular boundary flux for the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void StvenantFpBoundaryFlux2d(double* x,double t,double* wL,double* vn,
			    double* flux);
//! \brief particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransportInitData(double* x,double* w);
//! \brief particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void StvenantFpInitData(double* x,double* w);
//! \brief particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void SvenantFpInitData2d(double* x,double* w);
//! \brief particular imposed data for the 2d transport model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void SvenantFpImposedData2d(double* x,double t,double* w);

//! \brief particular flux for testing the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestSvenantFpBoundaryFlux2d(double* x,double t,double* wL,double* vn,
				 double* flux);
//! \brief particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestSvenantFpInitData(double* x,double* w);
//! \brief particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestSvenantFpInitData2d(double* x,double* w);
//! \brief particular imposed data for the transport model
//! \param[in] x,t : space and time position
//! \param[out] wt : imposed state at point x and time t
void TestStvenantFImposeData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d transport model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestSvenantFpImposedData2d(double* x,double t,double* w);


#endif
