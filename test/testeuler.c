//#include "macromesh.h"
//#include "geometry.h"
//#include "interpolation.h"
#include "test.h"
#include "euler.h"
#include "field.h"
#include "f2c.h"
#include "gasdyn77.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(void) {
  
  // unit tests
    
  int resu=TestEuler();
	 

  if (resu) printf("Model test OK !\n");
  else printf("Model test failed !\n");

  return !resu;
} 

int TestEuler(void){
  int test = true ;
  Field f;
  f.model.m=4; // only one conservative variable
  f.model.NumFlux=EulerNumFlux2d;
  f.model.BoundaryFlux=EulerBoundaryFlux2d;
  f.model.InitData=EulerInitData2d;
  f.model.ImposedData=EulerImposedData2d;
  f.varindex=GenericVarindex;


  f.interp.interp_param[0]=4;  // _M
  f.interp.interp_param[1]=1;  // x direction degree
  f.interp.interp_param[2]=1;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=2;  // x direction refinement
  f.interp.interp_param[5]=2;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement


  // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"disque.msh");
  // try to detect a 2d mesh
  bool is2d=Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  InitField(&f);
  f.is2d=true;


  // prudence...
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);


  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  RK2(&f,0.2);
 
  // save the results and the error
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);
  return test;


};
