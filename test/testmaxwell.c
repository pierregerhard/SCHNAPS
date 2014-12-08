#include "maxwell.h"
#include "test.h"
#include "field.h"
#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  int resu=TestMaxwell();
  if (resu) printf("Model test OK !\n");
  else printf("Model test failed !\n");
  return resu;
}  

int TestMaxwell(void){

  int test = (1==1);

  Field f;
  f.model.m=3;
  f.model.NumFlux=MaxwellNumFlux2dTM;
  f.model.BoundaryFlux=MaxwellBoundary2dTM;
  f.model.InitData=MaxwellInitData2dTM;
  f.model.ImposedData=MaxwellImposedData2dTM;
  f.varindex=GenericVarindex;

  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=1;  // x direction degree
  f.interp.interp_param[2]=1;  // y direction degree
  f.interp.interp_param[3]=1;  // z direction degree
  f.interp.interp_param[4]=1;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  bool is2d=Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);
  BuildConnectivity(&(f.macromesh));

  InitField(&f);
  f.is2d=true;
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);
  RK2(&f,1.);
  
  PlotField(0,(1==0),&f,"maxwelldg_field.msh");
  PlotField(0,(1==1),&f,"maxwelldg_error.msh");
  double dd=L2error(&f);
  printf("erreur L2=%f \n",dd);
  return test;


};
