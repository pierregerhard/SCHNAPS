#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "interpolation.h"
#include "model.h"

//! \brief Data structure for managing a  discrete vector field
//! solution of a DG approximation 
typedef struct Field{
  //! underlying mesh 
  MacroMesh macromesh;
  //! physical and numerical model
  Model model;
  //! interpolation used for each component of the field
  Interpolation interp;
  //! a copy of the interpolation parameters
  int interp_param[8];
  //! current time
  double tnow;
  //! cfl parameter min_i (vol_i / surf_i)
  double hmin;
  //! time step
  //! dt has to be smaller than hmin / vmax
  double dt;

  //! activate or not 2D computations
  bool is2d;

  //! fields at time steps n 
  double* wn;
  //! fields at time steps n+1   
  double* wnp1;
  //! time derivative of the field
  double* dtwn;

  //! \brief memory arrangement of field components
  //! \param[in] param interpolation parameters
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  int (*varindex)(int* param, int elem, int ipg, int iv);

} Field;


//! \brief a simple struct for packing a field
//! and a cell range.  To be passed to a thread
//! as a void* pointer.
typedef struct MacroCell{
  int first_cell; //!< first cell index
  int last_cell_p1;  //!< last cell index + 1
  Field* field; //! pointer to a the field
} MacroCell;

//! \brief memory arrangement of field components.
//! Generic implementation.
//! \param[in] param interpolation parameters
//! \param[in] elem macro element index
//! \param[in] ipg glop index
//! \param[in] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn. 
int GenericVarindex(int* param, int elem, int ipg, int iv);

//! field initialization. Computation of the initial
//! at each glop.
//! \param[inout] f a field
void InitField(Field* f);

//! \brief compute the Discontinuous Galerkin volume terms
//! slow version (no optimization using the tensors products)
//! \param[inout] f a field
void DGVolumeSlow(Field* f);

//! \brief compute the Discontinuous Galerkin inter-macrocells boundary terms.
//! Slow version  (no optimization using the tensors products)
//! \param[inout] f a field
void DGMacroCellInterfaceSlow(Field* f);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. One subcell implementation.
//! \param[inout] f a field
void dtFieldSlow(Field* f);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. Works with several subcells.
//! Fast version: multithreaded and with tensor products optimizations 
//! \param[inout] f a field
void dtField(Field* f);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void* DGMacroCellInterface(void* mcell);
//! \brief compute the Discontinuous Galerkin volume terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void* DGVolume(void* mcell);
//! \brief compute the Discontinuous Galerkin inter-subcells terms
//! \param[inout] mcell a MacroCell
void* DGSubCellInterface(void* mcell);
//! \brief  apply the DG mass term
//! \param[inout] mcell a MacroCell
void* DGMass(void* mcell);

//! \brief time integration by a second order Runge-Kutta algorithm 
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2(Field* f,double tmax);

void RK2StVenantLin(Field* f,double tmax);

//!  \brief time integration by a second order Runge-Kutta algorithm.
//! slow version
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2Copy(Field* f,double tmax);

//! \brief save the results in the gmsh format
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! \param[in] f a field
//! with the analytical solution
//! \param[in] filename the path to the gmsh visualization file.
void PlotField(int typplot,int compare,Field* f,char* filename);

//! \brief  display the field on screen
//! \param[in] f the field.
void DisplayField(Field* f);

//!  \brief compute the normalized L2 distance with the imposed data
//! \param[in] f the field.
//! \returns the error. 
double L2error(Field* f);


#endif
