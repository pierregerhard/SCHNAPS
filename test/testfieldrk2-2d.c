#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
int main(void) {
  
  // unit tests
    
  int resu=TestFieldRK2_2D();
	 
  if (resu) printf("Field RK2 2D test OK !\n");
  else printf("Field RK2 2D test failed !\n");

  return !resu;
} 



int TestFieldRK2_2D(void) {

  bool test=true;

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux2d;
  f.model.BoundaryFlux=TransportBoundaryFlux2d;
  f.model.InitData=TransportInitData2d;
  f.model.ImposedData=TransportImposedData2d;
  f.varindex=GenericVarindex;


  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=1;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement


  ReadMacroMesh(&(f.macromesh),"test/testdisque2d.msh");
  bool is2d=Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  InitField(&f);
  // require a 2d computation
  f.is2d=true;


  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);


  RK2(&f,0.2);
 
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);

  test = test && (dd<0.01);

  return test;

};
