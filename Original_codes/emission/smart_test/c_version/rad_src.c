/* rad_src.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rad_src__(integer *ng0, integer *l0, integer *nlyr, 
	integer *nphi, integer *nzdn, integer *numu, doublereal *trnrad_b__, 
	doublereal *rad_rt__, doublereal *rad_rb__, doublereal *flx_ft__, 
	doublereal *flx_fb__, doublereal *flx_dfb__, doublereal *flx_uft__, 
	real *rfldn, real *flup, real *uu, doublereal *radsrc)
{
    /* Initialized data */

    static doublereal pi = 3.141592654;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, l, ibn, naz, nze, nlev;
    static doublereal rad_ib__[70], rad_it__[70];


/* ccccccccccccccccccccccccccc  r a d _ s r c   cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine used a simplified flux adding method to find    cc */
/* c    the effective radiance source term in each layer, given the     cc */
/* c    radiance and flux field from DISORT.                            cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        ng0 - index  of current spectral bin                        cc */
/* c       nlyr - number of computational model layers                  cc */
/* c   trnrad_b - scattering transmittance of an isolate layer for bin  cc */
/* c     rad_rt - radiance reflectance of layer over this level         cc */
/* c     rad_rb - radiance reflectance of layer below this level        cc */
/* c         uu - angle dependent radiance at interface of 2 layers     cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       radsrc: radiance source for an isolated layer                cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  r a d _ s r c   cccccccccccccccccccccccccc */




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



/* ****   DISORT output */


/* ****   binned layer transmittances and absorptances */


/* ****   variables for radiance layer adding method */


/* ****  quantities for radiance calculation: downward and upward fluxes */
/*      at layer interfaces */


/* ****   internal binned source function variables */



/* ****   define pi */

    /* Parameter adjustments */
    radsrc -= 18193;
    uu -= 1137;
    --flup;
    --rfldn;
    flx_uft__ -= 71;
    flx_dfb__ -= 71;
    flx_fb__ -= 71;
    flx_ft__ -= 71;
    rad_rb__ -= 18193;
    rad_rt__ -= 18193;
    trnrad_b__ -= 6737;

    /* Function Body */

    ibn = *ng0;
    l = *l0;
    nlev = *nlyr + 1;

/* ****    use the modified flux adding method to find radiance sources */

    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {

	i__2 = *nzdn;
	for (nze = 1; nze <= i__2; ++nze) {

	    i__3 = nlev;
	    for (k = 1; k <= i__3; ++k) {
		rad_it__[k - 1] = uu[nze + (k + naz * 70 << 4)] - rad_rt__[
			nze + (naz + (k + l * 70 << 4) << 4)] * flx_fb__[k + 
			l * 70] / pi;
		rad_ib__[k - 1] = uu[*numu - nze + 1 + (k + naz * 70 << 4)] - 
			rad_rb__[*numu - nze + 1 + (naz + (k + l * 70 << 4) <<
			 4)] * flx_ft__[k + l * 70] / pi;
/* L1001: */
	    }

/* ****            find downward radiance sources, starting at top */

	    radsrc[nze + (naz + (l * 70 + 1 << 4) << 4)] = rad_it__[0];

	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {

/* ****               find downward radiances */

		radsrc[nze + (naz + (k + 1 + l * 70 << 4) << 4)] = rad_it__[k]
			 - trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] * (
			rad_it__[k - 1] + rad_rt__[nze + (naz + (k + l * 70 <<
			 4) << 4)] * flx_uft__[k + l * 70] / pi);
/* L1041: */
	    }

	    i__3 = *nlyr;
	    for (k = 1; k <= i__3; ++k) {

/* ****               find the upward radiances */

		radsrc[*numu - nze + 1 + (naz + (k + l * 70 << 4) << 4)] = 
			rad_ib__[k - 1] - trnrad_b__[*numu - nze + 1 + (k + (
			l + ibn * 5) * 70 << 4)] * (rad_ib__[k] + rad_rb__[*
			numu - nze + 1 + (naz + (k + 1 + l * 70 << 4) << 4)] *
			 flx_dfb__[k + 1 + l * 70] / pi);
/* L1061: */
	    }

	    radsrc[*numu - nze + 1 + (naz + (nlev + l * 70 << 4) << 4)] = 
		    rad_ib__[nlev - 1];

/* L1221: */
	}
/* L1241: */
    }

    return 0;
} /* rad_src__ */

