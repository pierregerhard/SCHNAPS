#include "macromesh.h"

#define _GNU_SOURCE  // for avoiding a compiler warning
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "geometry.h"
#include "interpolation.h"
#include <math.h>

void ReadMacroMesh(MacroMesh* m,char* filename){

  m->is2d=false;

  FILE* f=NULL;

  char* line=NULL;
  size_t linesize=0;
  size_t ret;

  printf("Read mesh file %s\n",filename);

  f=fopen(filename,"r");
  assert(f != NULL);

  

  do {
    ret=getline(&line,&linesize,f);
  }
  while(strcmp(line,"$Nodes\n") != 0);   

  // read the nodes data
  ret=getline(&line,&linesize,f);
  m->nbnodes=atoi(line);
  printf("nbnodes=%d\n",m->nbnodes);

  m->node=malloc(3 * m->nbnodes * sizeof(double));
  assert(m->node);

  for(int i=0;i<m->nbnodes;i++){    
    ret=getdelim(&line,&linesize,(int) ' ',f); // node number
    ret=getdelim(&line,&linesize,(int) ' ',f); // x
    m->node[3*i+0]=atof(line);
    ret=getdelim(&line,&linesize,(int) ' ',f); // y
    m->node[3*i+1]=atof(line);
    ret=getline(&line,&linesize,f); // z (end of the line)
    m->node[3*i+2]=atof(line);
    /* printf("Node %d x=%f y=%f z=%f\n",i, */
    /* 	   m->node[3*i+0], */
    /* 	   m->node[3*i+1], */
    /* 	   m->node[3*i+2]); */
  }    
  ret=getline(&line,&linesize,f);
  //printf("%s",line);
  // check that we have reached the end of nodes
  assert(strcmp(line,"$EndNodes\n") == 0);
  
  // Now read all the elements of the mesh
  do {
    ret=getline(&line,&linesize,f);
  }
  while(strcmp(line,"$Elements\n") != 0);   

  // size of the gmsh elems list 
  ret=getline(&line,&linesize,f);
  int nball=atoi(line);
  // allocate to a too big size
  m->elem2node=malloc(20 * sizeof(int) * nball);
  assert(m->elem2node);

  // now count only the H20 elems (code=17)
  m->nbelems=0;
  int countnode=0;
  for(int i=0;i<nball;i++){
    ret=getdelim(&line,&linesize,(int) ' ',f); // elem number
    ret=getdelim(&line,&linesize,(int) ' ',f); // elem type
    int elemtype=atoi(line);
    if (elemtype != 17) {
      ret=getline(&line,&linesize,f);
    }
    else {
      m->nbelems++;
      ret=getdelim(&line,&linesize,(int) ' ',f); //useless code
      ret=getdelim(&line,&linesize,(int) ' ',f); //useless code
      ret=getdelim(&line,&linesize,(int) ' ',f); //useless code
      for(int j=0;j<19;j++){
	ret=getdelim(&line,&linesize,(int) ' ',f);
	//printf("%d ",atoi(line));
	m->elem2node[countnode]=atoi(line)-1;
	countnode++;
      }
      ret=getline(&line,&linesize,f);
      //printf("%d\n",atoi(line));
      m->elem2node[countnode]=atoi(line)-1;
      countnode++;
    }
  }	
  ret=getline(&line,&linesize,f);
  //printf("%s",line);
  // check that we have reached the end of nodes
  assert(strcmp(line,"$EndElements\n") == 0);
  printf("nbelems=%d\n",m->nbelems);
  m->elem2node=realloc(m->elem2node,20 * sizeof(int) * m->nbelems);
  assert(m->elem2node);

  m->elem2elem=NULL;

}


void AffineMap(double* x){

  //double A[3][3]={{1,2,1},{0,-1,4},{7,8,-5}};
  double A[3][3]={{0,-1,0},{-2,0,0},{0,0,-1}};
  //double A[3][3]={1,0,0,0,2,0,0,0,1};
  double x0[3]={0,0,1};
  //double x0[3]={0,0,0};

  double newx[3];

  for(int i=0;i<3;i++){
    newx[i]=x0[i];
    for(int j=0;j<3;j++){
      newx[i]+=A[i][j]*x[j];
    }
  }
  x[0]=newx[0];
  x[1]=newx[1];
  x[2]=newx[2];

}


void AffineMapMacroMesh(MacroMesh* m){

    for(int ino=0;ino<m->nbnodes;ino++){
      AffineMap(&(m->node[ino*3]));
    }
}



// display macromesh data on standard output
void PrintMacroMesh(MacroMesh* m){
  printf("Print macromesh...\n");
  int start=1;
  printf("nbnodes=%d\n",m->nbnodes);
  for(int i=0;i<m->nbnodes;i++){
    printf("node %d x=%f y=%f z=%f\n",i+start,
	   m->node[3*i+0],
	   m->node[3*i+1],
	   m->node[3*i+2]);
  }
  printf("nbelems=%d\n",m->nbelems);
  for(int i=0;i<m->nbelems;i++){
    printf("elem %d -> ",i+start);
    for(int j=0;j<20;j++){
      printf("%d ",m->elem2node[20*i+j]+start);
    }
    printf("\n");
  }
  if (m->elem2elem !=0) {
    for(int i=0;i<m->nbelems;i++){
      printf("elem %d voisins: ",i+start);
      for(int j=0;j<6;j++){
	printf("%d ",m->elem2elem[6*i+j]+start);
      }
      printf("\n");
    }
  }
}

// build other connectivity arrays
void BuildConnectivity(MacroMesh* m){

  printf("Build connectivity...\n");

  assert(m->elem2elem==NULL);
  
  // build a list of faces
  // each face is made of four corners of 
  // the hexaedron mesh
  Face4Sort* face;
  Face4Sort* f;

  face=malloc(6*sizeof(Face4Sort)*m->nbelems);

  assert(face);

  int face2locnode[6][4]={
    {0,1,5,4},
    {1,2,6,5},
    {2,3,7,6},
    {0,4,7,3},
    {5,6,7,4},
    {0,3,2,1},
  };


  for(int ie=0;ie<m->nbelems;ie++){
    for(int ifa=0;ifa<6;ifa++){
      f=face+ifa+6*ie;
      for(int ino=0;ino<4;ino++){
	f->node[ino]=m->elem2node[face2locnode[ifa][ino]+20*ie];
      }
      f->left=ie;
      f->locfaceleft=ifa;
      f->right=-1;
      f->locfaceright=-1;
      OrderFace4Sort(f);
      /* printf("elem=%d ifa=%d left=%d nodes %d %d %d %d\n",ie,ifa, */
      /* 	     f->left,f->node[0], */
      /* 	     f->node[1],f->node[2],f->node[3]); */
    }
  }

  // now sort the list of faces
  qsort(face,6*m->nbelems,sizeof(Face4Sort),CompareFace4Sort);
  // check
  /* for(int ie=0;ie<m->nbelems;ie++){ */
  /*   for(int ifa=0;ifa<6;ifa++){ */
  /*     f=face+ifa+6*ie; */
  /*     printf("left=%d right=%d, nodes %d %d %d %d\n", */
  /* 	     f->left,f->right,f->node[0], */
  /* 	     f->node[1],f->node[2],f->node[3]); */
  /*   } */
  /* } */

  // allocate element connectivity array
  m->elem2elem=malloc(6 * m->nbelems * sizeof(int));
  for(int i=0;i<6 * m->nbelems;i++){
    m->elem2elem[i]=-1;
  }

  // now, two successive equal faces
  // correspond to two neighbours in the 
  // element list
  for(int ifa=0;ifa<6*m->nbelems-1;ifa++){
    Face4Sort* f1=face+ifa;
    Face4Sort* f2=face+ifa+1;
    if (CompareFace4Sort(f1,f2)==0){
      int ie1=f1->left;
      int if1=f1->locfaceleft;
      int ie2=f2->left;
      int if2=f2->locfaceleft;
	m->elem2elem[if1+6*ie1]=ie2;
	m->elem2elem[if2+6*ie2]=ie1;
    }
  }

  free(face);

  // check
  /* for(int ie=0;ie<m->nbelems;ie++){ */
  /*   for(int ifa=0;ifa<6;ifa++){ */
  /*     printf("elem=%d face=%d, voisin=%d\n", */
  /* 	     ie,ifa,m->elem2elem[ifa+6*ie]); */
  /*   } */
  /* } */

  
}

// compare two integers
int CompareInt(const void* a,const void* b){
  return(*(int*)a - *(int*)b);
}

// sort the nodes list of the face
void OrderFace4Sort(Face4Sort* f){
  qsort ( f->node, 4, sizeof(int), CompareInt);
}

// compare two ordered four-corner faces
// lexicographical order
int CompareFace4Sort(const void* a,const void* b){
  Face4Sort* f1= (Face4Sort*)a; 
  Face4Sort* f2= (Face4Sort*)b;

  int r= f1->node[0]-f2->node[0];
  if (r==0) r=f1->node[1]-f2->node[1];
  if (r==0) r=f1->node[2]-f2->node[2];
  if (r==0) r=f1->node[3]-f2->node[3];
  return r;
  
};

void CheckMacroMesh(MacroMesh* m,int* param){

  Geom g;

  double face_centers[6][3]={
    {0.5,0.0,0.5},
    {1.0,0.5,0.5},
    {0.5,1.0,0.5},
    {0.0,0.5,0.5},
    {0.5,0.5,1.0},
    {0.5,0.5,0.0},
  };

  /* double refnormal[6][3]={{0,-1,0},{1,0,0}, */
  /* 			  {0,1,0},{-1,0,0}, */
  /* 			  {0,0,1},{0,0,-1}}; */


  for(int ie=0;ie<m->nbelems;ie++){
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=m->elem2node[20*ie+inoloc];
      g.physnode[inoloc][0]=m->node[3*ino+0];
      g.physnode[inoloc][1]=m->node[3*ino+1];
      g.physnode[inoloc][2]=m->node[3*ino+2];
    }
    
    // test that the ref_ipg function
    // is compatible with ref_pg_vol
    //int param[7]={_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
    for(int ipg=0;ipg<NPG(param);ipg++){
      double xref1[3],xref_in[3];
      double wpg;
      ref_pg_vol(param,ipg,xref1,&wpg,xref_in);
      memcpy(g.xref,xref1,sizeof(g.xref));

      g.ifa=0;
      GeomRef2Phy(&g);
      GeomPhy2Ref(&g);

      // if (param[4]==1 && param[5]==1 && param[6]==1){ 
      //printf("ipg %d ipg2 %d xref %f %f %f\n",ipg,
      //	     ref_ipg(param,xref_in),xref_in[0],xref_in[1],xref_in[2]);
      assert(ipg==ref_ipg(param,xref_in));
	//}
    }

    // middle of the element
    g.xref[0]=0.5;
    g.xref[1]=0.5;
    g.xref[2]=0.5;

    GeomRef2Phy(&g);
    double xphym[3];
    memcpy(xphym,g.xphy,sizeof(xphym));
 
    for(int ifa=0;ifa<6;ifa++){
      // middle of the face
      memcpy(g.xref,face_centers[ifa],sizeof(g.xref));
      g.ifa=ifa;
      GeomRef2Phy(&g);
      // check volume  orientation
      assert(g.det >0);
      
      double vec[3]={g.xphy[0]-xphym[0],
		     g.xphy[1]-xphym[1],g.xphy[2]-xphym[2]};
      
      // check face orientation
      assert(g.vnds[0]*vec[0]+g.vnds[1]*vec[1]+g.vnds[2]*vec[2] > 0);

      // check compatibility between face and volume numbering
        for(int ipgf=0;ipgf<NPGF(param,ifa);ipgf++){
          double xpgref[3],wpg;
          // get the coordinates of the Gauss point
          ref_pg_face(param,ifa,ipgf,xpgref,&wpg,NULL);
          // recover the volume gauss point from
          // the face index
          int ipgv=param[6];
          double xpgref2[3],wpg2;
          ref_pg_vol(param,ipgv,xpgref2,&wpg2,NULL);
	  // in 2D do not check upper and lower face
	  if (m->is2d){
	    if (ifa !=4 && ifa!=5) {
	      assert(Dist(xpgref,xpgref2)<1e-11);
	    }
	  }
	  // in 3D check all faces
	  else {
	    assert(Dist(xpgref,xpgref2)<1e-11);
	  }
	  	      
        }


    }


  }


  // check that the faces are defined by the same mapping
  // with opposite normals
   for (int ie=0;ie<m->nbelems;ie++){
     //int param[8]={1,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=m->elem2node[20*ie+inoloc];
      physnode[inoloc][0]=m->node[3*ino+0];
      physnode[inoloc][1]=m->node[3*ino+1];
      physnode[inoloc][2]=m->node[3*ino+2];
    }

    // loop on the 6 faces
    for(int ifa=0;ifa<6;ifa++){
      // get the right elem or the boundary id
      int ieR=m->elem2elem[6*ie+ifa];
      double physnodeR[20][3];
      if (ieR >= 0) {
      	for(int inoloc=0;inoloc<20;inoloc++){
      	  int ino=m->elem2node[20*ieR+inoloc];
      	  physnodeR[inoloc][0]=m->node[3*ino+0];
      	  physnodeR[inoloc][1]=m->node[3*ino+1];
      	  physnodeR[inoloc][2]=m->node[3*ino+2];
      	}
      }
      
      // loop on the glops (numerical integration)
      // of the face ifa
      for(int ipgf=0;ipgf<NPGF(param,ifa);ipgf++){
  	double xpgref[3],xpgref_in[3],wpg;
  	//double xpgref2[3],wpg2;
  	// get the coordinates of the Gauss point
  	ref_pg_face(param,ifa,ipgf,xpgref,&wpg,xpgref_in);
	//printf("xref_in=%f %f %f\n",xpgref_in[0],xpgref_in[1],xpgref_in[2]);

  	// recover the volume gauss point from
  	// the face index
  	//int ipg=param[6];

	
  	// get the left value of w at the gauss point
  	// the basis functions is also the gauss point index
  	// normal vector at gauss point ipg
  	//double dpsiref[3];dpsi[3];
  	double dtau[3][3],codtau[3][3],xpg[3],xpg_in[3];
  	double vnds[3];

	// compute the "slightly inside" position
  	Ref2Phy(physnode,
  		xpgref_in,
  		NULL,ifa, // dpsiref,ifa
  		xpg_in,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds
	//printf("ie=%d ifa=%d xrefL=%f %f %f\n",ie,
	//     ifa,xpgref_in[0],xpgref_in[1],xpgref_in[2]);

	// compute the exact ref position
 	Ref2Phy(physnode,
  		xpgref,
  		NULL,ifa, // dpsiref,ifa
  		xpg,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds

  	if (ieR >=0) {  // the right element exists
  	  // find the corresponding point in the right elem
  	  double xref[3];
	  //printf("xpg_in=%f %f %f\n",xpg_in[0],xpg_in[1],xpg_in[2]);
	  Phy2Ref(physnodeR,xpg_in,xref);
          int ifaR=0;
          while (m->elem2elem[6*ieR+ifaR] != ie) ifaR++;
          assert(ifaR<6);
	  //printf("ieR=%d ifaR=%d xref=%f %f %f\n",ieR,
	  //	 ifaR,xref[0],xref[1],xref[2]);
  	  int ipgR=ref_ipg(param,xref);
	  //printf("ok\n");
	  double xpgR[3],xrefR[3],wpgR;
	  ref_pg_vol(param, ipgR, xrefR, &wpgR,NULL);
          double dtauR[3][3],codtauR[3][3];double vndsR[3];
	  Ref2Phy(physnodeR,
		  xrefR,
		  NULL,ifaR, // dphiref,ifa
		  xpgR,dtauR,  
		  codtauR,NULL,vndsR); // codtau,dphi,vnds
	  // printf("x1=%f %f %f x2=%f %f %f\n",xpg[0],xpg[1],xpg[2],
	  //xpgR[0],xpgR[1],xpgR[2]);
          assert(Dist(xpg,xpgR)<1e-11);
          assert(fabs(vnds[0]+vndsR[0])<1e-11);
          assert(fabs(vnds[1]+vndsR[1])<1e-11);
          assert(fabs(vnds[1]+vndsR[1])<1e-11);
        }
      }
    }
   }
  


};

// detect if the mesh is 2D
// and then permut the nodes so that
// the z direction coincides in the reference
// or physical frame
bool Detect2DMacroMesh(MacroMesh* m){

  m->is2d= true;


  // do not permut the node if the connectivity
  // is already built
  if (m->elem2elem != NULL)
    printf("Cannot permut nodes before building connectivity\n");
  assert(m->elem2elem == 0);

  for(int ie=0;ie<m->nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=m->elem2node[20*ie+inoloc];
      physnode[inoloc][0]=m->node[3*ino+0];
      physnode[inoloc][1]=m->node[3*ino+1];
      physnode[inoloc][2]=m->node[3*ino+2];
    }

    // we decide that the mesh is 2D if the 
    // middles of the elements have a constant z 
    // coordinate equal to 0.5
    double zmil=0;
    for(int inoloc=0;inoloc<20;inoloc++){
      zmil+=physnode[inoloc][2];
    }
    zmil/=20;
    // the mesh is not 2d
    if (fabs(zmil-0.5)>1e-6) {
      m->is2d=false;
      return m->is2d;
    }
  }

  printf("Detection of a 2D mesh\n");
  for(int ie=0;ie<m->nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=m->elem2node[20*ie+inoloc];
      physnode[inoloc][0]=m->node[3*ino+0];
      physnode[inoloc][1]=m->node[3*ino+1];
      physnode[inoloc][2]=m->node[3*ino+2];
    }
    // if the mesh is 2d permut the nodes
    // in order that the z^ and z axis are the 
    // same

    double face_centers[6][3]={
      {0.5,0.0,0.5},
      {1.0,0.5,0.5},
      {0.5,1.0,0.5},
      {0.0,0.5,0.5},
      {0.5,0.5,1.0},
      {0.5,0.5,0.0},
    };

    // rotation of the cube around the origin
    // at most two rotations are needed to put the cube
    // in a correct position
    for(int irot=0;irot<2;irot++){
      // compute the normal to face 4
      double vnds[3],dtau[3][3],codtau[3][3];
      Ref2Phy(physnode,
	      face_centers[4],
	      NULL,4, // dphiref,ifa
	      NULL,dtau,
	      codtau,NULL,vnds); // codtau,dphi,vnds

      double d=sqrt(vnds[0]*vnds[0]+vnds[1]*vnds[1]+vnds[2]*vnds[2]);
      // if the normal is not up or down
      // we have to permut the nodes
      if (fabs(vnds[2]/d)<0.9){
	printf("irot=%d rotating the element %d\n",irot,ie);
	int oldnum[20];
	int newnum[20]={1,5,6,2,4,8,7,3,11,9,10,17,18,13,19,12,16,14,20,15};
	for(int inoloc=0;inoloc<20;inoloc++) {
	  newnum[inoloc]--;
	  oldnum[inoloc]=m->elem2node[20*ie+inoloc];
	}
	// rotate the node numbering
	for(int inoloc=0;inoloc<20;inoloc++) {
	  m->elem2node[20*ie+inoloc]=oldnum[newnum[inoloc]];
	}
	// get the rotated nodes coordinates
	for(int inoloc=0;inoloc<20;inoloc++){
	  int ino=m->elem2node[20*ie+inoloc];
	  physnode[inoloc][0]=m->node[3*ino+0];
	  physnode[inoloc][1]=m->node[3*ino+1];
	  physnode[inoloc][2]=m->node[3*ino+2];
	}

      }
    }


  }


  return m->is2d;

};




  

