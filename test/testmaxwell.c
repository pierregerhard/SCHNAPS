#include "maxwell.h"
#include "test.h"
#include "field.h"
#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

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

  f.interp.interp_param[0]=3;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=3;  // x direction refinement
  f.interp.interp_param[5]=3;  // y direction refinement
  f.interp.interp_param[6]=3;  // z direction refinement

  ReadMacroMesh(&(f.macromesh),"test/disque2d.msh");
  bool is2d=Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);
  BuildConnectivity(&(f.macromesh));

  InitField(&f);
  f.is2d=true;
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);
  printf("cfl param =%f\n",f.hmin);
  
   RK2(&f,0.1);
  
  PlotField(0,(1==0),&f,"maxwelldg_field.msh");
  // PlotField(0,(1==1),&f,"maxwelldg_error.msh");
  // double dd=L2error(&f);
  // printf("erreur L2=%f \n",dd);
  return 1;


};
