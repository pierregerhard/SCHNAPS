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


//! \brief particular flux for the waves model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void WavesNumFlux(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular flux for the 2d waves model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void WavesNumFlux2d(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular boundary flux for the waves model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void WavesBoundaryFlux(double* x,double t,double* wL,double* vn,
			   double* flux);
//! \brief particular boundary flux for the 2d waves model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void WavesBoundaryFlux2d(double* x,double t,double* wL,double* vn,
			     double* flux);
//! \brief particular init data for the waves model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void WavesInitData(double* x,double* w);
//! \brief particular init data for the 2d waves model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void WavesInitData2d(double* x,double* w);
//! \brief particular imposed data for the waves model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void WavesImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d waves model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void WavesImposedData2d(double* x,double t,double* w);

//! \brief particular flux for testing the waves model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestWavesBoundaryFlux(double* x,double t,double* wL,double* vn,
			       double* flux);
//! \brief particular flux for testing the 2d waves model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestWavesBoundaryFlux2d(double* x,double t,double* wL,double* vn,
				 double* flux);
//! \brief particular init data for the waves model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestWavesInitData(double* x,double* w);
//! \brief particular init data for the 2d waves model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestWavesInitData2d(double* x,double* w);
//! \brief particular imposed data for the waves model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestWavesImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the 2d waves model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestWavesImposedData2d(double* x,double t,double* w);


#endif
