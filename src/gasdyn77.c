/* gasdyn77.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include "f2c.h"
/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/* $$$c234567 */
/* $$$      program euler */
/* $$$      implicit double precision (a-h,o-z) */
/* $$$ */
/* $$$      nmax=1000 */
/* $$$      a=-1 */
/* $$$      b=1 */
/* $$$      h=(b-a)/nmax */
/* $$$ */
/* $$$      rg=2 */
/* $$$      ug=0. */
/* $$$      pg=2 */
/* $$$      yg=1 */
/* $$$ */
/* $$$      rd=1 */
/* $$$      ud=0 */
/* $$$      pd=1 */
/* $$$      yd=0 */
/* $$$ */
/* $$$      time=0.5d0 */
/* $$$ */
/* $$$      open(1,file='rho') */
/* $$$      open(2,file='u') */
/* $$$      open(3,file='p') */
/* $$$      open(4,file='y') */
/* $$$      do i=0,nmax */
/* $$$         xi=(a+i*h)/time */
/* $$$         call riemann77(rg,ug,pg,yg,rd,ud,pd,yd,xi,r,u,p,y) */
/* $$$         write(1,*) a+i*h,r */
/* $$$         write(2,*) a+i*h,u */
/* $$$         write(3,*) a+i*h,p */
/* $$$         write(4,*) a+i*h,y */
/* $$$      enddo */
/* $$$ */
/* $$$      close(1) */
/* $$$      close(2) */
/* $$$      close(3) */
/* $$$      close(4) */
/* $$$ */
/* $$$      return */
/* $$$      end */
/* ce sous-programme renvoie les coefficients de la loi */
/* d'état en fonction de la fraction yv. */
/* $$$      subroutine eos(gam,pinf,yv) */
/* $$$ */
/* $$$      implicit double precision(a-h,o-z) */
/* $$$ */
/* $$$      gam1=1.4d0 */
/* $$$      gam2=1.4d0 */
/* $$$      pinf1=0 */
/* $$$      pinf2=0 */
/* $$$ */
/* $$$      t0=yv/(gam1-1.d0)+(1.d0-yv)/(gam2-1.d0) */
/* $$$      gam=1.d0+1.d0/t0 */
/* $$$      t0=yv*gam1*pinf1/(gam1-1.d0)+(1.d0-yv)*gam2*pinf2/(gam2-1.d0) */
/* $$$      pinf=(gam-1.d0)*t0/gam */
/* $$$ */
/* $$$      return */
/* $$$      end */
/* la suite concerne la résolution du problème de Riemann */
/* ... */
/* courbe de Hugoniot pour les chocs */
/* 234567 */
doublereal phia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal t0;
    extern doublereal ha_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    t0 = sqrt((d__1 = (*pl - *pa) * (*ta - ha_(pinfa, ga, ta, pa, pl)), abs(
	    d__1)));
    if (*pl <= *pa) {
	t0 = -t0;
    }
    ret_val = t0;
    return ret_val;
} /* phia_ */

/* dérivée par rapport à pl de la fonction précédente */
/* 234567 */
doublereal dphia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal t0, pi, pi0;
    extern doublereal dha_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    pi = *pl + *pinfa;
    pi0 = *pa + *pinfa;
    t0 = sqrt(*ta * 2. / ((*ga - 1.) * pi0 + (*ga + 1.) * pi)) * .5 - sqrt(((*
	    ga - 1.) * pi0 + (*ga + 1.) * pi) / 2. / *ta) * .5 * dha_(pinfa, 
	    ga, ta, pa, pl);
    if (*pl <= *pa) {
	t0 = -t0;
    }
    ret_val = t0;
    return ret_val;
} /* dphia_ */

/* calcul de tau=1/rho du côté (l) d'une 1 ou 3 onde en fonction de pl */
doublereal ha_(doublereal *pinfa, doublereal *ga, doublereal *ta, doublereal *
	pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal t, pi, pi0;

    pi = *pl + *pinfa;
    pi0 = *pa + *pinfa;
/* cas du choc */
    if (*pl > *pa) {
	t = *ta * ((*ga + 1.) * pi0 + (*ga - 1.) * pi) / ((*ga + 1.) * pi + (*
		ga - 1.) * pi0);
/* cas de la détente */
    } else {
	d__1 = pi0 / pi;
	d__2 = 1. / *ga;
	t = pow_dd(&d__1, &d__2) * *ta;
    }
    ret_val = t;
    return ret_val;
} /* ha_ */

/* dérivée par rapport à pl de la fonction précédente */
doublereal dha_(doublereal *pinfa, doublereal *ga, doublereal *ta, doublereal 
	*pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal t, pi, pi0;

    pi = *pl + *pinfa;
    pi0 = *pa + *pinfa;
/* choc */
    if (*pl > *pa) {
/* Computing 2nd power */
	d__1 = pi * (*ga + 1.) + pi0 * (*ga - 1.);
	t = *ga * -4. * *ta * pi0 / (d__1 * d__1);
/*     détente (ne sert qu'au déboguage) */
    } else {
	d__1 = 1. / *ga;
	d__2 = -(*ga + 1.) / *ga;
	t = -(*ta) * pow_dd(&pi0, &d__1) / *ga * pow_dd(&pi, &d__2);
    }
    ret_val = t;
    return ret_val;
} /* dha_ */

/* fonction puissance */
doublereal powa_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    ret_val = pow_dd(a, b);
    return ret_val;
} /* powa_ */

/* courbes isentropiques pour les détentes */
doublereal psia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c0, t1;
    extern doublereal powa_(doublereal *, doublereal *);

    c0 = sqrt(*ga * (*pa + *pinfa) * *ta);
/*      write(*,*) 'powa',pl,pa,pinfa,ga */
    d__1 = (*pl + *pinfa) / (*pa + *pinfa);
    d__2 = (*ga - 1.) / 2. / *ga;
    t1 = c0 * 2. / (*ga - 1.) * (powa_(&d__1, &d__2) - 1.);
    ret_val = t1;
    return ret_val;
} /* psia_ */

/* dérivée par rapport à pl de la fonction précédente */
doublereal dpsia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c0, t0;
    extern doublereal powa_(doublereal *, doublereal *);

    c0 = sqrt(*ga * (*pa + *pinfa) * *ta);
    d__1 = *pa + *pinfa;
    d__2 = (1. - *ga) / 2. / *ga;
    d__3 = *pl + *pinfa;
    d__4 = -(*ga + 1.) / 2. / *ga;
    t0 = c0 / *ga * powa_(&d__1, &d__2) * powa_(&d__3, &d__4);
    ret_val = t0;
    return ret_val;
} /* dpsia_ */

/* calcul de la pression en fonction de la vitesse dans une détente */
doublereal pp_(doublereal *pinfa, doublereal *ga, doublereal *ta, doublereal *
	pa, doublereal *psi)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, c0, t0;
    extern doublereal powa_(doublereal *, doublereal *);

    c0 = sqrt(*ga * (*pa + *pinfa) * *ta);
    c__ = (*ga - 1.) / (*ga + 1.) * (*psi + 2. / (*ga - 1.) * c0);
    d__1 = *pa + *pinfa;
    d__2 = 1. / *ga;
    t0 = c__ * c__ / *ga / *ta / powa_(&d__1, &d__2);
    d__1 = *ga / (*ga - 1.);
    t0 = powa_(&t0, &d__1) - *pinfa;
    ret_val = t0;
    return ret_val;
} /* pp_ */

/* courbes mixtes chocs-détentes */
doublereal xhia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal t0;
    extern doublereal phia_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), psia_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    if (*pl > *pa) {
	t0 = phia_(pinfa, ga, ta, pa, pl);
    } else {
	t0 = psia_(pinfa, ga, ta, pa, pl);
    }
    ret_val = t0;
    return ret_val;
} /* xhia_ */

/* dérivée par rapport à pl de la fonction précédente */
doublereal dxhia_(doublereal *pinfa, doublereal *ga, doublereal *ta, 
	doublereal *pa, doublereal *pl)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal t0;
    extern doublereal dphia_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dpsia_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    if (*pl - *pa > 0.) {
/*         write(*,*) 'choc 1',pl,pa,pl-pa */
	t0 = dphia_(pinfa, ga, ta, pa, pl);
/*         write(*,*) 'choc 2' */
    } else {
/*         write(*,*) 'detente1' */
	t0 = dpsia_(pinfa, ga, ta, pa, pl);
/*         write(*,*) 'detente2' */
    }
    ret_val = t0;
    return ret_val;
} /* dxhia_ */

/* fonction de résolution du problème de Riemann */
/* Subroutine */ int riemann77_(doublereal *rg, doublereal *ug, doublereal *
	pg, doublereal *rd, doublereal *ud, doublereal *pd, doublereal *xi, 
	doublereal *r__, doublereal *u, doublereal *p, doublereal *um)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer iterriem;
    static doublereal f, r1, r2;
    extern doublereal ha_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal fd, gd, ff, fg, gg, df, dp;
    static integer it;
    static doublereal pn, pm;
    extern doublereal pp_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal um1, um2, eps, err, psi, vit[5], gamd, gamg;
    extern doublereal xhia_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal pinf;
    extern doublereal psia_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal psid, crit, psig, pmin;
    extern doublereal dxhia_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal pinfd;
    extern doublereal dpsia_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal pinfg;
    static integer newton;
    static doublereal errriem;

    /* Fortran I/O blocks */
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };


    psig = 1.;
    psid = 1.;
/*      common/iter/iterriem/err/errriem */
/* initialisation */
/*      call eos(gamg,pinfg,psig) */
/*      call eos(gamd,pinfd,psid) */
    gamg = 1.4;
    gamd = 1.4;
    pinfg = 0.;
    pinfd = 0.;
    fg = gamg - 1.;
    fd = gamd - 1.;
/* initialisation */
    pn = 0.;
    err = 1.;
    eps = 1e-14;
    it = 0;
    pmin = -pinfg;
    if (pmin < -pinfd) {
	pmin = -pinfd;
    }
    gg = fg + 1.;
    gd = fd + 1.;
/* critère d'apparition du vide */
/*      write(*,*) 'crit1',pd,eps*pd,pg,eps*pd */
    d__1 = fd + 1.;
    d__2 = 1. / *rd;
    d__3 = pmin + eps * *pd;
    d__4 = fg + 1.;
    d__5 = 1. / *rg;
    d__6 = pmin + eps * *pg;
    crit = *ug - *ud - xhia_(&pinfd, &d__1, &d__2, pd, &d__3) - xhia_(&pinfg, 
	    &d__4, &d__5, pg, &d__6);
/*      write(*,*) 'crit2' */
    if (crit <= 0.) {
	s_wsle(&io___36);
	do_lio(&c__9, &c__1, "apparition du vide", (ftnlen)18);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
/* valeur initiale newton */
/*      pn=pg */
/*      if (pn.lt.pd) pn=pd */
    pn = (*pg + *pd) * .5;
/* nombre max d'itérations */
    newton = 100;
    while(err > eps && it < newton) {
	++it;
/*        write(*,*) it,err,'riem=',rg,ug,pg,rd,ud,pd */
/*     terme de gauche (voir Rouy) */
	ff = *ug - *ud;
	df = 0.;
	d__1 = 1. / *rg;
	ff -= xhia_(&pinfg, &gg, &d__1, pg, &pn);
	d__1 = 1. / *rg;
	df -= dxhia_(&pinfg, &gg, &d__1, pg, &pn);
/*     terme de droite (voir Rouy) */
/*         write(*,*) 'coucou1' */
	d__1 = 1. / *rd;
	ff -= xhia_(&pinfd, &gd, &d__1, pd, &pn);
/*         write(*,*) 'coucou2',rd,pd,pn */
	d__1 = 1. / *rd;
	df -= dxhia_(&pinfd, &gd, &d__1, pd, &pn);
/*         write(*,*) 'coucou3' */
	dp = ff / df;
	pn -= dp;
/*         if (dabs(dp).lt.dabs(ff)) then */
/*            err = dabs(ff) */
/*         else */
	err = (d__1 = dp / pn, abs(d__1));
/*         endif */
    }
    if (it > iterriem) {
	iterriem = it;
    }
    if (err > errriem) {
	errriem = err;
    }
    if (it == newton) {
	s_wsle(&io___43);
	do_lio(&c__9, &c__1, "non convergence ", (ftnlen)16);
	e_wsle();
	s_wsle(&io___44);
	do_lio(&c__9, &c__1, "ff,dp,pn", (ftnlen)8);
	do_lio(&c__5, &c__1, (char *)&ff, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&pn, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___45);
	do_lio(&c__9, &c__1, "rg,ug,pg,rd,ud,pd,psid,psig", (ftnlen)27);
	do_lio(&c__5, &c__1, (char *)&(*rg), (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&(*ug), (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&(*pg), (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&(*rd), (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&(*ud), (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&(*pd), (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&psid, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&psig, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    pm = pn;
    d__1 = 1. / *rd;
    r2 = 1. / ha_(&pinfd, &gd, &d__1, pd, &pn);
    d__1 = 1. / *rg;
    r1 = 1. / ha_(&pinfg, &gg, &d__1, pg, &pn);
/* vitesse de la discontinuité de contact */
    d__1 = 1. / *rg;
    um1 = *ug - xhia_(&pinfg, &gg, &d__1, pg, &pn);
    d__1 = 1. / *rd;
    um2 = *ud + xhia_(&pinfd, &gd, &d__1, pd, &pn);
    *um = (um1 + um2) * .5;
/*     vitesses caractéristiques */
/*     1-détente */
    if (pm <= *pg) {
	d__1 = 1. / *rg;
	vit[0] = *ug - 1. / dpsia_(&pinfg, &gg, &d__1, pg, pg) / *rg;
	d__1 = 1. / *rg;
	vit[1] = um1 - 1. / dpsia_(&pinfg, &gg, &d__1, pg, &pm) / r1;
/*     1-choc */
    } else {
	vit[0] = *ug - sqrt(((gg + 1.) * (pm + pinfg) + (gg - 1.) * (*pg + 
		pinfg)) * .5 / *rg);
	vit[1] = vit[0];
    }
/*     contact */
    vit[2] = *um;
/*     3-détente */
    if (pm <= *pd) {
	d__1 = 1. / *rd;
	vit[3] = um2 + 1. / dpsia_(&pinfd, &gd, &d__1, pd, &pm) / r2;
	d__1 = 1. / *rd;
	vit[4] = *ud + 1. / dpsia_(&pinfd, &gd, &d__1, pd, pd) / *rd;
/*     3-choc */
    } else {
	vit[3] = *ud + sqrt(((gd + 1.) * (pm + pinfd) + (gd - 1.) * (*pd + 
		pinfd)) * .5 / *rd);
	vit[4] = vit[3];
    }
/* cette fonction calcule les variables primitives */
/* de la solution du problème de Riemann en xi=x/t */
/* état gauche */
    if (*xi < vit[0]) {
	*r__ = *rg;
	*p = *pg;
	*u = *ug;
	f = fg;
	pinf = pinfg;
	psi = psig;
    }
/*     calcul dans la 1-détente */
    if (*xi >= vit[0] && *xi < vit[1]) {
	d__1 = fg + 1.;
	d__2 = 1. / *rg;
	d__3 = *ug - *xi;
	*p = pp_(&pinfg, &d__1, &d__2, pg, &d__3);
	d__1 = fg + 1.;
	d__2 = 1. / *rg;
	*r__ = 1. / ha_(&pinfg, &d__1, &d__2, pg, p);
	f = fg;
	pinf = pinfg;
	psi = psig;
	d__1 = fg + 1.;
	d__2 = 1. / *rg;
	*u = *ug - psia_(&pinfg, &d__1, &d__2, pg, p);
    }
/*     état (I) (voir Rouy) */
    if (*xi >= vit[1] && *xi < vit[2]) {
	*r__ = r1;
	*p = pm;
	*u = um1;
	f = fg;
	pinf = pinfg;
	psi = psig;
    }
/*     état (II) (voir Rouy) */
    if (*xi >= vit[2] && *xi <= vit[3]) {
	*r__ = r2;
	*p = pm;
	*u = um2;
	f = fd;
	pinf = pinfd;
	psi = psid;
    }
/* calcul dans la 3-détente */
    if (*xi > vit[3] && *xi < vit[4]) {
	d__1 = fd + 1.;
	d__2 = 1. / *rd;
	d__3 = *xi - *ud;
	*p = pp_(&pinfd, &d__1, &d__2, pd, &d__3);
	d__1 = fd + 1.;
	d__2 = 1. / *rd;
	*r__ = 1. / ha_(&pinfd, &d__1, &d__2, pd, p);
	f = fd;
	pinf = pinfd;
	psi = psid;
	d__1 = fd + 1.;
	d__2 = 1. / *rd;
	*u = *ud + psia_(&pinfd, &d__1, &d__2, pd, p);
    }
/*     état droit */
    if (*xi >= vit[4]) {
	*r__ = *rd;
	*p = *pd;
	*u = *ud;
	f = fd;
	pinf = pinfd;
	psi = psid;
    }
    return 0;
} /* riemann77_ */

