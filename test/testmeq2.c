#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"



int main(void) {
  
  // unit tests
    
  int resu=TestmEq2();
	 
  if (resu) printf("m greater than 1 test OK !\n");
  else printf("m greater than 1 test failed !\n");

  return !resu;
} 


int TestmEq2(void) {

  bool test=true;

  Field f;

  int vec=2;
  
  if (vec==2) {
    f.model.m=2; // num of conservative variables
    f.model.NumFlux=VecTransNumFlux2d;
    f.model.BoundaryFlux=VecTransBoundaryFlux2d;
    f.model.InitData=VecTransInitData2d;
    f.model.ImposedData=VecTransImposedData2d;
    f.varindex=GenericVarindex;
    
    
    f.interp.interp_param[0]=2;  // _M
    f.interp.interp_param[1]=2;  // x direction degree
    f.interp.interp_param[2]=2;  // y direction degree
    f.interp.interp_param[3]=0;  // z direction degree
    f.interp.interp_param[4]=4;  // x direction refinement
    f.interp.interp_param[5]=4;  // y direction refinement
    f.interp.interp_param[6]=1;  // z direction refinement
  }
  else {
    f.model.m=1; // num of conservative variables
    f.model.NumFlux=TransportNumFlux2d;
    f.model.BoundaryFlux=TestTransportBoundaryFlux2d;
    f.model.InitData=TestTransportInitData2d;
    f.model.ImposedData=TestTransportImposedData2d;
    f.varindex=GenericVarindex;
    
    
    f.interp.interp_param[0]=1;  // _M
    f.interp.interp_param[1]=2;  // x direction degree
    f.interp.interp_param[2]=2;  // y direction degree
    f.interp.interp_param[3]=0;  // z direction degree
    f.interp.interp_param[4]=4;  // x direction refinement
    f.interp.interp_param[5]=4;  // y direction refinement
    f.interp.interp_param[6]=1;  // z direction refinement
  }
  // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
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

  // time derivative
  //dtField(&f);
  //DisplayField(&f);
  //assert(1==2);
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  RK2(&f,0.1);
 
  // save the results and the error
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);
  test= test && (dd<1e-7);
  return test;

};


