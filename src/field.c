#include "field.h"
#include "geometry.h"
#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "global.h"
#include <math.h>

// param[0] = M
// param[1] = deg x
// param[2] = deg y
// param[3] = deg z
// param[4] = raf x
// param[5] = raf y
// param[6] = raf z
int GenericVarindex(int* param, int elem, int ipg, int iv){

  int npg= (param[1]+1)*(param[2]+1)*(param[3]+1);

  return iv + param[0] * ( ipg + npg * elem);

}

void InitField(Field* f){

  int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ};
  double w[f->model.m];
  double xpg[3];
  double xref[3],omega;
  double dtau[9];
  double physnode[20*3];

  int nmem=f->model.m * f->macromesh.nbelems * 
    NPG(param+1);
  printf("allocate %d doubles\n",nmem);
  f->wn=malloc(nmem * sizeof(double));
  assert(f->wn);	       
  f->wnp1=malloc(nmem * sizeof(double));
  assert(f->wnp1);	       
  f->dtwn=malloc(nmem * sizeof(double));
  assert(f->dtwn);	       

  f->tnow=0;

  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc*3+0]=f->macromesh.node[3*ino+0];
      physnode[inoloc*3+1]=f->macromesh.node[3*ino+1];
      physnode[inoloc*3+2]=f->macromesh.node[3*ino+2];
    }
    for(int ipg=0;ipg<NPG(param+1);ipg++){
      ref_pg_vol(param+1, ipg, xref, &omega);
      Ref2Phy(physnode,
	      xref,
	      0,0, // dphiref,ifa
              xpg,dtau,  
	      NULL,NULL,NULL); // codtau,dphi,vnds
      // check the reverse transform at all the GLOPS
      double xref2[3];
      Phy2Ref(physnode,xpg,xref2);
      assert(sqrt((xref2[0]-xref[0])*(xref2[0]-xref[0])+
      		(xref2[0]-xref[0])*(xref2[0]-xref[0])+
      		  (xref2[0]-xref[0])*(xref2[0]-xref[0])) < 1e-8);
      
      f->model.InitData(xpg,w);
      //printf("xpg %f %f %f w=%f\n",xpg[0],xpg[1],xpg[2],w[0]);
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(param,ie,ipg,iv);
	f->wn[imem]=w[iv];
	//printf("imem=%d\n",imem);
      }
    }
  }

};

// display the field on screen
void DisplayField(Field* f){
  int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ};
  //int param[8]={f->model.m,2,2,2,1,1,1};

  printf("Display field...\n");
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    printf("elem %d\n",ie);
    for(int ipg=0;ipg<NPG(param+1);ipg++){
    printf("Gauss point %d\n",ipg);
    printf("w= ");
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(param,ie,ipg,iv);
	printf("%f ",f->wn[imem]);
      }
      printf("\n");
    printf("dtw= ");
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(param,ie,ipg,iv);
	printf("%f ",f->dtwn[imem]);
      }
      printf("\n");
    }
  }



};



// save the results in the gmsh format
void PlotField(Field* f,char* filename){

  const int hexa64ref[3*64]={
    0,0,3,
    3,0,3,
    3,3,3,
    0,3,3,
    0,0,0,3,0,0,3,3,0,0,3,0,
    1,0,3,2,0,3,0,1,3,0,2,3,0,0,2,0,0,1,3,1,3,3,2,3,
    3,0,2,3,0,1,2,3,3,1,3,3,3,3,2,3,3,1,0,3,2,0,3,1,
    1,0,0,2,0,0,0,1,0,0,2,0,3,1,0,3,2,0,2,3,0,1,3,0,
    1,1,3,1,2,3,2,2,3,2,1,3,1,0,2,2,0,2,2,0,1,1,0,1,
    0,1,2,0,1,1,0,2,1,0,2,2,3,1,2,3,2,2,3,2,1,3,1,1,
    2,3,2,1,3,2,1,3,1,2,3,1,1,1,0,2,1,0,2,2,0,1,2,0,
    1,1,2,2,1,2,2,2,2,1,2,2,1,1,1,2,1,1,2,2,1,1,2,1};

  int* elem2nodes = f->macromesh.elem2node;
  double* node = f->macromesh.node;


  FILE * gmshfile;
  gmshfile = fopen( filename, "w" );

  // data plots
  int mw = f->model.m;
  int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ};
  int degre = param[1];
  int npgf = NPGF(param+1,0);
  int npgv = NPG(param+1);
  int nnodes = 20;

  double Xn[3*nnodes];
  double Xr[3];
  double Xphy[3];
  double Vnds[3];  // normal vector times surface element dS
  double DX[3*3];  // jacobian matrix
  double coDX[3*3];  // comatrix of the jacobian matrix


  int vindex_param[] = {npgv, mw};

  // header
  fprintf(gmshfile,"$MeshFormat\n2.2 0 %d\n",(int) sizeof(double));
  //int one=1;
  //fwrite((char*) &one,sizeof(int),1,gmshfile);
  fprintf(gmshfile,"$EndMeshFormat\n$Nodes\n%d\n",f->macromesh.nbelems*64);
  // nodes
  for(int i=0;i<f->macromesh.nbelems;i++){
    // get the nodes of element L
    for(int ino=0;ino<nnodes;ino++){
      int numnoe=elem2nodes[nnodes*i+ino];
      for(int ii=0;ii<3;ii++){
        Xn[3*ino+ii]=node[3*numnoe+ii];
      }
    }
    for(int ino=0;ino<64;ino++){
      int nnoe=64*i+ino+1;
      Xr[0]=(double) (hexa64ref[3*ino+0]) / 3;
      Xr[1]=(double) (hexa64ref[3*ino+1]) / 3;
      Xr[2]=(double) (hexa64ref[3*ino+2]) / 3;

      Ref2Phy(Xn,
	      Xr,
	      NULL,
	      0,
	      Xphy,
	      DX,
	      coDX,
	      NULL,
	      Vnds);

      double Xplot[3];
      Xplot[0]=Xphy[0];
      Xplot[1]=Xphy[1];
      Xplot[2]=Xphy[2];

      // fwrite((char*) &nnoe,sizeof(int),1,gmshfile);
      // fwrite((char*) &(Xplot[0]),sizeof(double),1,gmshfile);
      // fwrite((char*) &(Xplot[1]),sizeof(double),1,gmshfile);
      // fwrite((char*) &(Xplot[2]),sizeof(double),1,gmshfile);
      fprintf(gmshfile,"%d %f %f %f\n",nnoe,Xplot[0],Xplot[1],Xplot[2]);
      
    }
  }
  // elements

  fprintf(gmshfile,"$EndNodes\n");
  fprintf(gmshfile,"$Elements\n");
  fprintf(gmshfile,"%d\n",f->macromesh.nbelems);



  int elm_type=92;
  int num_elm_follow=f->macromesh.nbelems;
  int num_tags=0;

  // fwrite((char*) &elm_type,sizeof(int),1,gmshfile);
  // fwrite((char*) &num_elm_follow,sizeof(int),1,gmshfile);
  // fwrite((char*) &num_tags,sizeof(int),1,gmshfile);

  for(int i=0;i<f->macromesh.nbelems;i++){
    
    int numelem=i+1;
    //fwrite((char*) &numelem,sizeof(int),1,gmshfile);
    fprintf(gmshfile,"%d ",numelem);
    fprintf(gmshfile,"%d ",elm_type);
    fprintf(gmshfile,"%d ",num_tags);
    
    
    for(int ii=0;ii<64;ii++){
      int numnoe=64*i + ii +1;
      //fwrite((char*) &numnoe,sizeof(int),1,gmshfile);
      fprintf(gmshfile,"%d ",numnoe);
    }
    fprintf(gmshfile,"\n");
  }
  
  fprintf(gmshfile,"$EndElements\n");
  
  
  // now display data
  for(int typplot=0;typplot<mw;typplot++){
    
    fprintf(gmshfile,"$NodeData\n");
    fprintf(gmshfile,"1\n");
    fprintf(gmshfile,"\"Field %d\"\n",typplot);

    double t = 0;

    // f << 1 <<std::endl;
    // f << t << std::endl;
    // f << 3 << std::endl;
    // f << 0 << std::endl;
    // f << 1 << std::endl;
    fprintf(gmshfile,"1\n%f\n3\n0\n1\n",t);

    //f << 64 * zf->zone_mesh()->nb_elems() <<std::endl;
    fprintf(gmshfile,"%d\n",64*f->macromesh.nbelems);

    
    for(int i=0;i<f->macromesh.nbelems;i++){
      for(int ino=0;ino<20;ino++){
	int numnoe=elem2nodes[nnodes*i+ino];
	for(int ii=0;ii<3;ii++){
	  Xn[3*ino+ii]=node[3*numnoe+ii];
	}
      }
      
      // data at the eight nodes
      for(int ii=0;ii<64;ii++){
	int nodenumber=64*i + ii +1;
	
	Xr[0]=(double) (hexa64ref[3*ii+0]) / 3;
	Xr[1]=(double) (hexa64ref[3*ii+1]) / 3;
	Xr[2]=(double) (hexa64ref[3*ii+2]) / 3;
	
	Ref2Phy(Xn,
		Xr,
		NULL,
		0,
		Xphy,
		DX,
		NULL,
		NULL,
		NULL);
	

	double value=0;
	for(int ib=0;ib<npgv;ib++){
	  double psi;
	  psi_ref(param+1, ib, Xr, &psi, NULL);
	  
	  int vi = f->varindex(param, i, ib, typplot);
	  //printf("i=%d psi=%f w=%f\n",vi,psi,f->wn[vi]);
	  value += psi * f->dtwn[vi];

	}

	// Exact solution
	// clac::real* w = new clac::real[zf->model()->nb_vars()];
	// zf->model()->boundary_data(Xphy, 2.5, w);
	// value = w[typplot];
	// delete[] w;

	//fwrite(const void *ptr, size_t size_of_elements, size_t number_of_elements, FILE *a_file);
	//fwrite((char*) &nodenumber, sizeof(int),1,gmshfile);
	//fwrite((char*) &value, sizeof(double),1,gmshfile);
	//fprintf(gmshfile,"%d %f\n",nodenumber,value);
	fprintf(gmshfile,"%d %f\n",nodenumber,value);
      }

    }

    fprintf(gmshfile,"\n$EndNodeData\n");

  }

    
  fclose(gmshfile);
  

}

// apply the Discontinuous Galerkin approximation for computing
// the time derivative of the field
void dtField(Field* f){

  // interpolation params
  // warning: this is ugly, but the last
  // parameter is used for computing the volume
  // GLOP index from the face GLOP index...
  // ugly too: the first parameter is not used by all
  // utilities. we have sometimes to jump over : pass param+1
  // instead of param...
  int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ};

  // init to zero the time derivative
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int ipg=0;ipg<NPG(param+1);ipg++){
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(param,ie,ipg,iv);
	f->dtwn[imem]=0;
      }
    }
  }

  // 
  // assembly of the surface terms
  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20*3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc*3+0]=f->macromesh.node[3*ino+0];
      physnode[inoloc*3+1]=f->macromesh.node[3*ino+1];
      physnode[inoloc*3+2]=f->macromesh.node[3*ino+2];
    }

    // loop on the 6 faces
    for(int ifa=0;ifa<6;ifa++){
      // get the right elem or the boundary id
      int ieR=f->macromesh.elem2elem[6*ie+ifa];
      double physnodeR[20*3];
      // if we are not at a boundary
      // get the nodes of the right elem
      if (ieR >= 0) {
  	for(int inoloc=0;inoloc<20;inoloc++){
  	  int ino=f->macromesh.elem2node[20*ieR+inoloc];
  	  physnodeR[inoloc*3+0]=f->macromesh.node[3*ino+0];
  	  physnodeR[inoloc*3+1]=f->macromesh.node[3*ino+1];
  	  physnodeR[inoloc*3+2]=f->macromesh.node[3*ino+2];
  	}
      }
      
      // loop on the glops (numerical integration)
      // of the face ifa
      for(int ipgf=0;ipgf<NPGF(param+1,ifa);ipgf++){
  	double xpgref[3],wpg;
  	double xpgref2[3],wpg2;
  	// get the coordinates of the Gauss point
  	ref_pg_face(param+1,ifa,ipgf,xpgref,&wpg);

  	// recover the volume gauss point from
  	// the face index
  	int ipg=param[7];
  	/* ref_pg_vol(param+1,ipg,xpgref2,&wpg2); */
	/* assert(fabs(xpgref2[0]-xpgref[0]) */
	/*        +fabs(xpgref2[1]-xpgref[1]) */
	/*        +fabs(xpgref2[2]-xpgref[2])<1e-10); */
  	// get the left value of w at the gauss point
  	double wL[3],wR[3];
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(param,ie,ipg,iv);
  	  wL[iv]=f->wn[imem];
  	}
  	// the basis functions is also the gauss point index
  	int ib=ipg;
  	// normal vector at gauss point ipg
  	//double dpsiref[3];dpsi[3];
  	double dtau[3*3],codtau[3*3],xpg[3];
  	double vnds[3];
  	Ref2Phy(physnode,
  		xpgref,
  		NULL,ifa, // dpsiref,ifa
  		xpg,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds
	//printf("ifa=%d ipg=%d vnds=%f %f %f \n",ifa,ipg,
	// vnds[0],vnds[1],vnds[2]);
  	double flux[f->model.m];
  	if (ieR >=0) {  // the right element exists
  	  // find the corresponding point in the right elem
  	  double xref[3];
  	  Phy2Ref(physnodeR,xpg,xref);
  	  int ipgR=ref_ipg(param+1,xref);
  	  for(int iv=0;iv<f->model.m;iv++){
  	    int imem=f->varindex(param,ieR,ipgR,iv);
  	    wR[iv]=f->wn[imem];
  	  }
  	  // int_dL F(wL,wR,grad phi_ib )
  	  f->model.NumFlux(wL,wR,vnds,flux);
	  assert(1==2);
  	}
  	else { //the right element does not exist
  	  f->model.BoundaryFlux(xpg,f->tnow,wL,vnds,flux);
  	}
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(param,ie,ib,iv);
  	  f->dtwn[imem]+=flux[iv]*wpg;
  	}
	
      }

    }
  }

  // assembly of the volume terms
  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20*3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc*3+0]=f->macromesh.node[3*ino+0];
      physnode[inoloc*3+1]=f->macromesh.node[3*ino+1];
      physnode[inoloc*3+2]=f->macromesh.node[3*ino+2];
    }

    // mass matrix
    double masspg[NPG(param+1)];
    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(param+1);ipg++){
      double xpgref[3],wpg;
      // get the coordinates of the Gauss point
      ref_pg_vol(param+1,ipg,xpgref,&wpg);

      // get the value of w at the gauss point
      double w[3];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // loop on the basis functions
      for(int ib=0;ib<NPG(param+1);ib++){
	// gradient of psi_ib at gauss point ipg
	double dpsiref[3],dpsi[3];
	double dtau[3*3],codtau[3*3];//,xpg[3];
	grad_psi_pg(param+1,ib,ipg,dpsiref);
	Ref2Phy(physnode, // phys. nodes
		xpgref,  // xref
		dpsiref,NULL, // dpsiref,ifa
		NULL,dtau,  // xphy,dtau
		codtau,dpsi,NULL); // codtau,dpsi,vnds
	// remember the diagonal mass term
	if (ib == ipg){
	  double det=dtau[0]*codtau[0]+dtau[1]*codtau[1]+dtau[2]*codtau[2];
	  masspg[ipg]=wpg*det;
	}
	// int_L F(w,w,grad phi_ib )
	double flux[f->model.m];
	f->model.NumFlux(w,w,dpsi,flux);
	for(int iv=0;iv<f->model.m;iv++){
	  int imem=f->varindex(param,ie,ib,iv);
	  f->dtwn[imem]-=flux[iv]*wpg;
	}
      }
    }
    for(int ipg=0;ipg<NPG(param+1);ipg++){
      // apply the inverse of the diagonal mass matrix
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(param,ie,ipg,iv);
	f->dtwn[imem]/=masspg[ipg];
	//	printf("masspg=%f\n",masspg);	
      }
    }

  }

};




