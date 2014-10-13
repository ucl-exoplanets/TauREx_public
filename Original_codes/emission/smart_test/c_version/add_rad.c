/* add_rad.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int add_rad__(integer *nlyr, integer *nzdn, integer *nzup, 
	integer *numu, integer *nphi, doublereal *trnrad, doublereal *refrad, 
	doublereal *brdf, doublereal *rad_src0__, doublereal *ft, doublereal *
	fb, doublereal *urt, doublereal *drb, doublereal *uft, doublereal *
	dfb, doublereal *rad0)
{
    /* Initialized data */

    static doublereal pi = 3.141592654;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k;
    static doublereal rb[17920]	/* was [16][16][70] */, rt[17920]	/* 
	    was [16][16][70] */, fdn[17920]	/* was [16][16][70] */;
    static integer naz;
    static doublereal fup[17920]	/* was [16][16][70] */;
    static integer nze, nlev;
    static doublereal rad_ib__[17920]	/* was [16][16][70] */, rad_it__[
	    17920]	/* was [16][16][70] */;


/* cccccccccccccccccccccccccc  a d d _ r a d   ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this program uses the layer adding method for finding fluxes    cc */
/* c    in a vertically inhomogeneous scattering, absorbing atmosphere. cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr - number of computational model layers                  cc */
/* c         rl - diffuse flux reflectance of each layer                cc */
/* c         tl - diffuse flux transmittance of each homogeneous layer  cc */
/* c        fup - f-, upward flux at the top of each homogeneous layer  cc */
/* c        fdn - f+, downward flux at the base of each layer           cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      radu - upward radiance at nlyr+1 layer boundries.             cc */
/* c             (flxu(l) refers to the upward flux at the top          cc */
/* c              of layer l)                                           cc */
/* c      radd - downward radiance at nlyr+1 layer boundries.           cc */
/* c             (flxd(l) refers to the downward flux at the bottom     cc */
/* c              of layer l-1)                                         cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  a d d _ r a d   ccccccccccccccccccccccccccc */




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


/* ****   surface bi-directional reflection function */


/* ****    radiance source terms */


/* ****    flux quantitities used on radiance adding */



/* ****   variables for flux layer adding method */



/* ****   output upward and downward radiances */



    /* Parameter adjustments */
    rad0 -= 273;
    --dfb;
    --uft;
    --drb;
    --urt;
    --fb;
    --ft;
    rad_src0__ -= 273;
    brdf -= 17;
    refrad -= 17;
    trnrad -= 17;

    /* Function Body */

    nlev = *nlyr + 1;

/* ****  find layer-integrated reflectivities.  the quantities are: */
/*      rt(k): reflectance of an inhomogeneous layer extendind from */
/*             level 1 to level k */
/*      rb(k): reflectance of an inhomogeneous layer extendind from */
/*             the surface (nlev) to level k */

    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {

/* ****       reflectances of layers extending from the top of the */
/*           atmosphere to each level for downward streams */

	i__2 = *nzdn;
	for (nze = 1; nze <= i__2; ++nze) {
	    rt[nze + (naz + 16 << 4) - 273] = refrad[nze + 16];
	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {
		rt[nze + (naz + (k + 1 << 4) << 4) - 273] = refrad[nze + (k <<
			 4)] + trnrad[nze + (k << 4)] * rt[nze + (naz + (k << 
			4) << 4) - 273] * urt[k] / pi;
/* L1201: */
	    }

/* L1221: */
	}

/* ****       reflectances of layers extending from the surface to */
/*           each level of the atmosphere */

	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {
	    rb[nze + (naz + (nlev << 4) << 4) - 273] = brdf[nze + (naz << 4)];
	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {
		rb[nze + (naz + (nlev - k << 4) << 4) - 273] = refrad[nze + (
			nlev - k << 4)] + trnrad[nze + (nlev - k << 4)] * rb[
			nze + (naz + (nlev - k + 1 << 4) << 4) - 273] * drb[
			nlev - k + 1] / pi;
/* L1241: */
	    }
/* L1261: */
	}
/* L1281: */
    }

/* ****   rescale source functions with current values */
/*       note: the upward radiances at the surface (k=nlev) are */
/*             not scaled since tl(nlev) = 0.0 */

    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {

/* ****       scale downward radiance source terms */

	i__2 = *nzdn;
	for (nze = 1; nze <= i__2; ++nze) {
	    fdn[nze + (naz + 16 << 4) - 273] = rad_src0__[nze + (naz + 16 << 
		    4)];
	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {
		fdn[nze + (naz + (k + 1 << 4) << 4) - 273] = rad_src0__[nze + 
			(naz + (k + 1 << 4) << 4)] - trnrad[nze + (k << 4)] * 
			rt[nze + (naz + (k << 4) << 4) - 273] * uft[k] / pi;
/*                  fdn(nze,naz,k+1) = rad_src0(nze,naz,k+1) */
/* L1401: */
	    }
/* L1421: */
	}

/* ****       scale upward radiance source terms */

	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {
	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {
		fup[nze + (naz + (k << 4) << 4) - 273] = rad_src0__[nze + (
			naz + (k << 4) << 4)] - trnrad[nze + (k << 4)] * rb[
			nze + (naz + (k + 1 << 4) << 4) - 273] * dfb[k + 1] / 
			pi;
/*                  fup(nze,naz,k) = rad_src0(nze,naz,k) */
/* L1441: */
	    }
	    fup[nze + (naz + (nlev << 4) << 4) - 273] = rad_src0__[nze + (naz 
		    + (nlev << 4) << 4)];
/* L1461: */
	}
/* L1481: */
    }

    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {

/* ****       use adding method to find downward radinaces */
/*           for combined layers. */

	i__2 = *nzdn;
	for (nze = 1; nze <= i__2; ++nze) {

/*              start at top, adding homogeneous layers to the base */
/*               of the existing inhomogeneous layer */

	    rad_it__[nze + (naz + 16 << 4) - 273] = fdn[nze + (naz + 16 << 4) 
		    - 273];
	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {
		rad_it__[nze + (naz + (k + 1 << 4) << 4) - 273] = fdn[nze + (
			naz + (k + 1 << 4) << 4) - 273] + trnrad[nze + (k << 
			4)] * (rad_it__[nze + (naz + (k << 4) << 4) - 273] + 
			rt[nze + (naz + (k << 4) << 4) - 273] * uft[k] / pi);
/* L2001: */
	    }
/* L2021: */
	}

	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {

/* ****            use adding method to find upward and downward radiances */
/*                for combined layers starting at bottom of atmosphere. */

	    rad_ib__[nze + (naz + (nlev << 4) << 4) - 273] = fup[nze + (naz + 
		    (nlev << 4) << 4) - 273];
	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {
		rad_ib__[nze + (naz + (nlev - k << 4) << 4) - 273] = fup[nze 
			+ (naz + (nlev - k << 4) << 4) - 273] + trnrad[nze + (
			nlev - k << 4)] * (rad_ib__[nze + (naz + (nlev - k + 
			1 << 4) << 4) - 273] + rb[nze + (naz + (nlev - k + 1 
			<< 4) << 4) - 273] * dfb[nlev - k + 1] / pi);
/* L2041: */
	    }
/* L2061: */
	}
/* L2081: */
    }

/* ****  find the total upward and downward fluxes at interfaces */
/*      between inhomogeneous layers. */

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {

/* ****           use adding method to find downward radinaces */
/*               at layer interfaces */

	    i__3 = *nzdn;
	    for (nze = 1; nze <= i__3; ++nze) {
		rad0[nze + (naz + (k << 4) << 4)] = rad_it__[nze + (naz + (k 
			<< 4) << 4) - 273] + rt[nze + (naz + (k << 4) << 4) - 
			273] * fb[k] / pi;
/* L3001: */
	    }

/* ****           use adding method to find downward radinaces */
/*               at layer interfaces */

	    i__3 = *numu;
	    for (nze = *nzup; nze <= i__3; ++nze) {
		rad0[nze + (naz + (k << 4) << 4)] = rad_ib__[nze + (naz + (k 
			<< 4) << 4) - 273] + rb[nze + (naz + (k << 4) << 4) - 
			273] * ft[k] / pi;
/* L3021: */
	    }
/* L3041: */
	}
/* L3061: */
    }

    return 0;
} /* add_rad__ */

