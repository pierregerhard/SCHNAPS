//#include "macromesh.h"
//#include "geometry.h"
//#include "interpolation.h"
#include "test.h"
#include "stvenantlin.h"
#include "field.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(void) {
  
  // unit tests
    
  int resu=TestStvenantFp();
	 

  if (resu) printf("St venant test OK !\n");
  else printf("St venant test failed !\n");

  return !resu;
} 

int TestStvenantFp(void){


  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=StvenantFpNumFlux;
  f.model.BoundaryFlux=StvenantFpBoundaryFlux;
  f.model.InitData=StvenantFpInitData;
  f.model.ImposedData=StvenantFpImposedData;
  f.varindex=GenericVarindex;


  f.interp.interp_param[0]=3;  // _M
  f.interp.interp_param[1]=3;  // x direction degree
  f.interp.interp_param[2]=3;  // y direction degree
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
  RK2(&f,1.);
 
  // save the results and the error
  PlotField(0,(1==0),&f,"dgvisustvenantlin.msh");
  PlotField(0,(1==1),&f,"dgerrorstvenantlin.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);
  return 0;

};


