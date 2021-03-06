#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
int main(void) {
  
  // unit tests
    
  int resu=TestFieldSubCellDGVol();
	 
  if (resu) printf("Field DG Subcell Vol test OK !\n");
  else printf("Field DG Subcell Vol test failed !\n");

  return !resu;
} 




int TestFieldSubCellDGVol(void){

  int test = (1==1);

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TestTransportBoundaryFlux;
  f.model.InitData=TestTransportInitData;
  f.model.ImposedData=TestTransportImposedData;
  f.varindex=GenericVarindex;

  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=2;  // z direction degree
  f.interp.interp_param[4]=2;  // x direction refinement
  f.interp.interp_param[5]=2;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement


  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testdisque.msh");
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));
  //AffineMapMacroMesh(&(f.macromesh));
  PrintMacroMesh(&(f.macromesh));

  InitField(&f);
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);


  MacroCell mcell[f.macromesh.nbelems];

  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    mcell[ie].field=&f;
    mcell[ie].first_cell=ie;
    mcell[ie].last_cell_p1=ie+1;
  }

  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    DGMacroCellInterface((void*) (mcell+ie));
  }
  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    DGSubCellInterface((void*) (mcell+ie));
  }
  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    DGVolume((void*) (mcell+ie));
  }
  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    DGMass((void*) (mcell+ie));
  }


  /* DGMacroCellInterface(&f); */
  /* DGSubCellInterface(&f); */

  /* DGVolume(&f); */

  /* DGMass(&f); */
  
  DisplayField(&f);  

  int yes_compare = 1;
  int no_compare = 0;

  PlotField(0,no_compare,&f,"visu.msh");
  PlotField(0,yes_compare,&f,"error.msh");

  // test the time derivative with the exact solution
  for(int i=0;i<f.model.m * f.macromesh.nbelems * 
	NPG(f.interp.interp_param+1);i++){
    test = test && fabs(4*f.wn[i]-pow(f.dtwn[i],2))<1e-2;
    assert(test);
  }
  
  return test;



};
