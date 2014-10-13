/* heating.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;

/* Subroutine */ int heating_(integer *iuheat, integer *nlyr, integer *nza, 
	real *umu0_1__, doublereal *wnmin, doublereal *wnmax, real *p, real *
	t, real *alt, doublereal *aid_lr__, doublereal *dirsoflx, doublereal *
	dnsoflx, doublereal *upsoflx, doublereal *dnthflx, doublereal *
	upthflx, doublereal *soheat, doublereal *thheat)
{
    /* Initialized data */

    static real pi = 3.141593f;

    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    alist al__1;

    /* Builtin functions */
    integer f_rew(alist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    double acos(doublereal);

    /* Local variables */
    static real thheat_d__, soheat_d__;
    static integer k, nz;
    static real pbar, tbar, altkm, solang;
    static doublereal thflxn[70], soflxn[70];

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, "(/,/,1a)", 0 };
    static cilist io___8 = { 0, 0, 0, "(/,/,1a,2(1pe14.6,1a),/,1a,1pe13.5,1a"
	    ",/,/,3a,/,3a)", 0 };
    static cilist io___10 = { 0, 0, 0, "(3(1pe12.4),24x,5(1pe13.5),2(1pe13.4"
	    "))", 0 };
    static cilist io___15 = { 0, 0, 0, "(5(1pe12.4),5(1pe13.5),2(1pe13.4))", 
	    0 };



/* ccccccccccccccccccccccccccc  h e a t i n g  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes flux divergences adn heating rates     cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      nlyr - number of of layers in the model atmosphere            cc */
/* c       nza - number of solar zenith angles for which heating rates  cc */
/* c             are needed.                                            cc */
/* c     flxdn - spectrally-integrated downward flux (W/m**2).          cc */
/* c     flxup - spectrally-integrated upward flux (W/m**2).            cc */
/* c    aid_lr - aidiabatic lapse rate (g/cp) in mks units              cc */
/* c         p - pressure at each level (pascals)                       cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      soheat - solar radiative heating rate (K/sec)                 cc */
/* c      thheat - thermal radiative heating rate (K/sec)               cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  h e a t i n g  ccccccccccccccccccccccccccc */




/* ccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc */
/* c								      cc */
/* c    p u r p o s e :                                                 cc */
/* c								      cc */
/* c    this include file specifies the dimensions of variables used    cc */
/* c    in the program smt_do.                                          cc */
/* c								      cc */
/* c    q u a n t i t i e s :                                           cc */
/* c								      cc */
/* c      nsp - Maximum number of spectral intervals 		      cc */
/* c       kp - Maximum number of vertical levels.		      cc */
/* c    nmode - Maximum number of aerosol particle modes.               cc */
/* c     ngas - Maximum number of absorbing gases.		      cc */
/* c     ngrp - Maximum number of SMT grouping bins.		      cc */
/* c    mxumu - Maximum number of user output zenith angles in discrete cc */
/* c            ordinate model.					      cc */
/* c    mxphi - Maximum number of output azimuth angles in the discrete cc */
/* c            ordinate model.					      cc */
/* c    mxmom - Maximum number of legendre polynomial moments in the    cc */
/* c            discrete ordinate method.				      cc */
/* c    mxulv - Maximum number of output levels in the discrete 	      cc */
/* c            ordinate model.					      cc */
/* c    mxrad - Maximum number of output radiance levels.		      cc */
/* c     nsol - Maximum number of solar zenith angles.		      cc */
/* c      nex - Maximum number of spectrally-dependent optical          cc */
/* c            (gases, aerosol modes, surface albedos, and solar flux. cc */
/* c     mxpd - maximum number of radiance partial derivatives          cc */
/* c								      cc */
/* ccccccccccccccccccccccccccccc  p a r a m  ccccccccccccccccccccccccccccc */




/* ****   atmospheric structure variables */


/* ****   double precision variables for heating rate integration */



    /* Parameter adjustments */
    --thheat;
    soheat -= 71;
    --upthflx;
    --dnthflx;
    upsoflx -= 71;
    dnsoflx -= 71;
    dirsoflx -= 71;
    --aid_lr__;
    --alt;
    --t;
    --p;
    --umu0_1__;

    /* Function Body */

/* ****     define solar heating rates at each solar zenith angle */

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {

/* ****      define the net flux at each level: */

	i__2 = *nlyr + 1;
	for (k = 1; k <= i__2; ++k) {
	    soflxn[k - 1] = dirsoflx[k + nz * 70] + dnsoflx[k + nz * 70] - 
		    upsoflx[k + nz * 70];
/* L1001: */
	}

/* ****       define the heating rate in each layer */

	i__2 = *nlyr;
	for (k = 1; k <= i__2; ++k) {

	    soheat[k + nz * 70] = aid_lr__[k] * (soflxn[k - 1] - soflxn[k]) / 
		    (p[k + 1] - p[k]);
/* L1021: */
	}
/* L1201: */
    }

/* ****    define the net flux at each level: */

    i__1 = *nlyr + 1;
    for (k = 1; k <= i__1; ++k) {
	thflxn[k - 1] = dnthflx[k] - upthflx[k];
/* L2001: */
    }

/* ***** define thermal cooling rates (independent of zenith angle) */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	thheat[k] = aid_lr__[k] * (thflxn[k] - thflxn[k - 1]) / (p[k + 1] - p[
		k]);
/* L2021: */
    }
    al__1.aerr = 0;
    al__1.aunit = *iuheat;
    f_rew(&al__1);
    io___6.ciunit = *iuheat;
    s_wsfe(&io___6);
    do_fio(&c__1, " R a d i a t i v e     H e a t i n g    R a t e s ", (
	    ftnlen)50);
    e_wsfe();

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	solang = acos(umu0_1__[nz]) * 180.f / pi;
	io___8.ciunit = *iuheat;
	s_wsfe(&io___8);
	do_fio(&c__1, "Spectral Range =", (ftnlen)16);
	do_fio(&c__1, (char *)&(*wnmin), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " -", (ftnlen)2);
	do_fio(&c__1, (char *)&(*wnmax), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " wavenumbers", (ftnlen)12);
	do_fio(&c__1, "Solar Zenith Angle =", (ftnlen)20);
	do_fio(&c__1, (char *)&solang, (ftnlen)sizeof(real));
	do_fio(&c__1, " degrees", (ftnlen)8);
	do_fio(&c__1, "   pressure  temperature   altitude     solar Q ", (
		ftnlen)48);
	do_fio(&c__1, "   thermal Q  dir sol flx  dn sol flx   up sol flx", (
		ftnlen)50);
	do_fio(&c__1, "   dn th flx   up th flx     p(flux)     t(flux)", (
		ftnlen)48);
	do_fio(&c__1, "     (bar)       (K)         (km)      (K/day) ", (
		ftnlen)47);
	do_fio(&c__1, "    (K/day)   (watts/m*m)  (watts/m*m)  (watts/m*m)", (
		ftnlen)51);
	do_fio(&c__1, "   (watts/m*m)  (watts/m*m)     (bar)       (K)", (
		ftnlen)47);
	e_wsfe();

	pbar = p[1] * 1e-5f;
	io___10.ciunit = *iuheat;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&pbar, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&t[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&alt[1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dirsoflx[nz * 70 + 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&dnsoflx[nz * 70 + 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&upsoflx[nz * 70 + 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&dnthflx[1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&upthflx[1], (ftnlen)sizeof(doublereal));
	r__1 = p[1] * 1e-5f;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&t[1], (ftnlen)sizeof(real));
	e_wsfe();
	i__2 = *nlyr;
	for (k = 1; k <= i__2; ++k) {
	    pbar = (p[k] + p[k + 1]) * 5e-6f;
	    tbar = (t[k] + t[k + 1]) * .5f;
	    altkm = (alt[k] + alt[k + 1]) * .5f;
	    soheat_d__ = (real) soheat[k + nz * 70] * 86400.f;
	    thheat_d__ = (real) thheat[k] * 86400.f;
	    io___15.ciunit = *iuheat;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&pbar, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&tbar, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&altkm, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&soheat_d__, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&thheat_d__, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dirsoflx[k + 1 + nz * 70], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&dnsoflx[k + 1 + nz * 70], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&upsoflx[k + 1 + nz * 70], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&dnthflx[k + 1], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&upthflx[k + 1], (ftnlen)sizeof(doublereal))
		    ;
	    r__1 = p[k + 1] * 1e-5f;
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&t[k + 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L6201: */
	}
/* L6221: */
    }

    return 0;
} /* heating_ */

