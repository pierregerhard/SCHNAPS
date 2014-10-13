#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  
  // unit tests
    
  int resu=TestField();
	 
  if (resu) printf("Field test OK !\n");
  else printf("Field test failed !\n");

  return !resu;
} 




int TestField(void){

  int test = (1==1);

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TestTransportBoundaryFlux;
  f.model.InitData=TestTransportInitData;
  f.model.ImposedData=TestTransportImposedData;
  f.varindex=GenericVarindex;

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));

  InitField(&f);

  PlotField(0,(1==1),&f,"testvisufield.msh");
  
  return test;


};