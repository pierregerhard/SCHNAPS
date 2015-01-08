#ifndef _GASDYN77_H_
#define _GASDYN77_H_
#include "f2c.h"
doublereal phia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl) ;

doublereal dphia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl) ;

doublereal ha_(doublereal *pinfa, doublereal *ga, doublereal *ta, doublereal *
	pa, doublereal *pl) ;

doublereal dha_(doublereal *pinfa, doublereal *ga, doublereal *ta, doublereal 
	*pa, doublereal *pl) ;

doublereal powa_(doublereal *a, doublereal *b) ;

doublereal psia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl) ;

doublereal dpsia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl) ;

doublereal pp_(doublereal *pinfa, doublereal *ga, doublereal *ta, doublereal *
	pa, doublereal *psi) ;

doublereal xhia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl) ;

doublereal dxhia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl) ;

int riemann77_(doublereal *rg, doublereal *ug, doublereal *
	pg, doublereal *rd, doublereal *ud, doublereal *pd, doublereal *xi, 
	doublereal *r__, doublereal *u, doublereal *p, doublereal *um) ;

#endif

