/* add_tflx.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int add_tflx__(integer *nlyr, doublereal *tl, doublereal *rl,
	 doublereal *fup, doublereal *fdn, doublereal *urt, doublereal *drb, 
	doublereal *uft, doublereal *dfb, doublereal *ft, doublereal *fb, 
	doublereal *flxu, doublereal *flxd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k, l;
    static doublereal dd[70], rb[70], du[70], rt[70];
    static integer nlev;


/* cccccccccccccccccccccccccc  a d d _ t f l x   ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this program uses the layer adding method for finding thermal   cc */
/* c    fluxes in a vertically inhomogeneous scattering, absorbing      cc */
/* c    atmosphere.                                                     cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr - number of computational model layers                  cc */
/* c         rl - diffuse flux reflectance of each layer                cc */
/* c         tl - diffuse flux transmittance of each homogeneous layer  cc */
/* c        fup - f+, upward flux at the top of each homogeneous layer  cc */
/* c        fdn - f-, downward flux at the base of each layer           cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      flxu - upward flux at nlyr+1 layer boundries.                 cc */
/* c             (flxu(l) refers to the upward flux at the top          cc */
/* c              of layer l)                                           cc */
/* c      flxd - downward flux at nlyr+1 layer boundries.               cc */
/* c             (flxd(l) refers to the downward flux at the bottom     cc */
/* c              of layer l-1)                                         cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  a d d _ t f l x   ccccccccccccccccccccccccc */




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



/* ****   layer flux transmittances, absorptances, and fluxes */


/* ****   variables for flux layer adding method */


/* ****    quantitities used on radiance adding */


/* ****   output */


    /* Parameter adjustments */
    --flxd;
    --flxu;
    --fb;
    --ft;
    --dfb;
    --uft;
    --drb;
    --urt;
    --fdn;
    --fup;
    --rl;
    --tl;

    /* Function Body */
    nlev = *nlyr + 1;

/* ****  find layer-integrated reflectivities.  the quantities are: */
/*      rt(k): reflectance of an inhomogeneous layer extendind from */
/*             level 1 to level k */
/*      rb(k): reflectance of an inhomogeneous layer extendind from */
/*             the surface (nlev) to level k */

    rt[0] = 0.f;
    rb[nlev - 1] = rl[nlev];
    urt[nlev] = 0.;
    drb[1] = 0.;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	dd[k - 1] = 1. / (1. - rt[k - 1] * rl[k]);
	urt[k] = tl[k] * dd[k - 1];
/* Computing 2nd power */
	d__1 = tl[k];
	rt[k] = rl[k] + d__1 * d__1 * rt[k - 1] * dd[k - 1];
	du[nlev - k - 1] = 1. / (1. - rl[nlev - k] * rb[nlev - k]);
	drb[nlev - k + 1] = tl[nlev - k] * du[nlev - k - 1];
/* Computing 2nd power */
	d__1 = tl[nlev - k];
	rb[nlev - k - 1] = rl[nlev - k] + d__1 * d__1 * rb[nlev - k] * du[
		nlev - k - 1];
/* L1201: */
    }

/* ****   use adding method to find upward and downward fluxes */
/*       for combined layers.  start at top, adding homogeneous */
/*       layers to the base of the existing inhomogeneous layer */

    ft[1] = 0.;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	ft[k + 1] = fdn[k + 1] + tl[k] * (ft[k] + rt[k - 1] * fup[k]) * dd[k 
		- 1];
/* L1401: */
    }

/* ****   use adding method to find upward and downward fluxes */
/*       for combined layers.  start at bottom. */

    fb[nlev] = fup[nlev];
    i__1 = *nlyr;
    for (l = 1; l <= i__1; ++l) {
	k = nlev - l;
	fb[k] = fup[k] + tl[k] * (fb[k + 1] + rb[k] * fdn[k + 1]) * du[k - 1];
/* L1601: */
    }

/* ****  find the total upward and downward fluxes at interfaces */
/*      between inhomogeneous layers. */

    flxd[1] = fdn[1];
    i__1 = nlev;
    for (k = 2; k <= i__1; ++k) {
	flxd[k] = (ft[k] + rt[k - 1] * fb[k]) / (1. - rt[k - 1] * rb[k - 1]);
/* L2001: */
    }

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	flxu[k] = fb[k] + rb[k - 1] * flxd[k];
/* L2021: */
    }

/* ****   find the variables dft and ufb for radiance adding */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	uft[k] = (fup[k] + rl[k] * ft[k]) * dd[k - 1];
	dfb[k + 1] = (fdn[k + 1] + rl[k] * fb[k + 1]) * du[k - 1];
/* L3001: */
    }

    return 0;
} /* add_tflx__ */

