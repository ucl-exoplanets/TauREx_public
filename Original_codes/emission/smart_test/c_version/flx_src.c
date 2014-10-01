/* flx_src.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int flx_src__(integer *ng0, integer *l0, integer *nlyr, 
	doublereal *trnflx_b__, doublereal *refflx_b__, doublereal *flx_rt__, 
	doublereal *flx_rb__, doublereal *flx_du__, doublereal *flx_dd__, 
	real *rfldn, real *flup, doublereal *flx_fb__, doublereal *flx_ft__, 
	doublereal *flx_dft__, doublereal *flx_uft__, doublereal *flx_dfb__, 
	doublereal *flx_ufb__, doublereal *upflxsrc, doublereal *dnflxsrc)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k, l, ibn, nlev;


/* cccccccccccccccccccccccccc  f l x _ s r c   ccccccccccccccccccccccccccc */
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
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       upflxsrc: f- (Crisp 1986 eq 7) upward flux at base           cc */
/* c                     of an isolated layer.                          cc */
/* c       dnflxsrc: f+ (Crisp 1986 eq 8) downward flux at top          cc */
/* c                    of isolated  layer.                             cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  f l x _ s r c   ccccccccccccccccccccccccccc */




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


/* ****   binned layer flux transmittances and absorptances */


/* ****   variables for radiance layer adding method */


/* ****  quantities for radiance calculation: downward and upward fluxes */
/*      at layer interfaces */


/* ****   output variables */


    /* Parameter adjustments */
    dnflxsrc -= 71;
    upflxsrc -= 71;
    flx_ufb__ -= 71;
    flx_dfb__ -= 71;
    flx_uft__ -= 71;
    flx_dft__ -= 71;
    flx_ft__ -= 71;
    flx_fb__ -= 71;
    --flup;
    --rfldn;
    flx_dd__ -= 71;
    flx_du__ -= 71;
    flx_rb__ -= 71;
    flx_rt__ -= 71;
    refflx_b__ -= 421;
    trnflx_b__ -= 421;

    /* Function Body */
    nlev = *nlyr + 1;
    l = *l0;
    ibn = *ng0;

/* ****     use the flux adding method to solve for flx_ft, the */
/*         downward diffuse flux at the base of each */
/*         inhomogeneous layer extending from layers 1 to k */
/*         and flx_fb, the upward flux at the top of each */
/*         inhomogeneous layer extending from the surface to k */

/* ****   find the downward flux at the base of layer k */
/*                  (level k+1) (dflux) eq 22 Crisp 1986) */

    flx_ft__[l * 70 + 1] = 0.;
    i__1 = nlev;
    for (k = 2; k <= i__1; ++k) {
	flx_ft__[k + l * 70] = rfldn[k] - flx_rt__[k + l * 70] * flup[k];
/* L1001: */
    }

/* ****   find the upward flux at the top of an inhomogeneous */
/*       layer extending from level k to the surface (uflux) */
/*       (simplifed version of Eq 23 of Crisp 1986) */

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	flx_fb__[k + l * 70] = flup[k] - flx_rb__[k + l * 70] * rfldn[k];
/* L1201: */
    }

/* ****   find the fluxes at the top of each homogeneous layer */
/*       (upflxsrc) and the bottom of each layer */

    dnflxsrc[l * 70 + 1] = 0.f;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

/*            solve for the upward flux at the top of each */
/*            homogenous layer,upflxsrc (ufl in delted). */
/*            To derive this equation, I started with Eqs. 18 and */
/*            21 of Crisp 1986, solved for f+ and f-, respectively, */
/*            and then substituted the value of f+ into the expression */
/*            for f- */

/*            upflxsrc = f+ */

/* Computing 2nd power */
	d__1 = trnflx_b__[k + (l + ibn * 5) * 70];
	upflxsrc[k + l * 70] = (flx_fb__[k + l * 70] - trnflx_b__[k + (l + 
		ibn * 5) * 70] * (flx_fb__[k + 1 + l * 70] + flx_rb__[k + 1 + 
		l * 70] * (flx_ft__[k + 1 + l * 70] - trnflx_b__[k + (l + ibn 
		* 5) * 70] * flx_ft__[k + l * 70] * flx_dd__[k + l * 70])) * 
		flx_du__[k + l * 70]) / (1. - d__1 * d__1 * flx_rb__[k + 1 + 
		l * 70] * flx_rt__[k + l * 70] * flx_dd__[k + l * 70] * 
		flx_du__[k + l * 70]);

/*             solve for the downward flux at the base of each */
/*                   homogenous layer,f- = dnflxsrc */

	dnflxsrc[k + 1 + l * 70] = flx_ft__[k + 1 + l * 70] - trnflx_b__[k + (
		l + ibn * 5) * 70] * (flx_ft__[k + l * 70] + flx_rt__[k + l * 
		70] * upflxsrc[k + l * 70]) * flx_dd__[k + l * 70];

/* L2001: */
    }
    upflxsrc[nlev + l * 70] = flx_fb__[nlev + l * 70];

/* ****   find the variables flx_dft and flx_ufb for radiance adding */

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	flx_uft__[k + l * 70] = (upflxsrc[k + l * 70] + refflx_b__[k + (l + 
		ibn * 5) * 70] * flx_ft__[k + l * 70]) * flx_dd__[k + l * 70];
	flx_dft__[k + l * 70] = (flx_ft__[k + l * 70] + flx_rt__[k + l * 70] *
		 upflxsrc[k + l * 70]) * flx_dd__[k + l * 70];
/* L3001: */
    }

    flx_dfb__[l * 70 + 1] = 0.f;

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	flx_dfb__[k + 1 + l * 70] = (dnflxsrc[k + 1 + l * 70] + refflx_b__[k 
		+ (l + ibn * 5) * 70] * flx_fb__[k + 1 + l * 70]) * flx_du__[
		k + l * 70];
	flx_ufb__[k + l * 70] = flx_fb__[k + l * 70] + flx_rb__[k + l * 70] * 
		flx_dfb__[k + l * 70];
/* L3021: */
    }
    flx_ufb__[nlev + l * 70] = flx_fb__[nlev + l * 70] + flx_rb__[nlev + l * 
	    70] * flx_dfb__[nlev + l * 70];

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*       write(*,'(1a,16(1pe12.4))') 'rfldn   ', */
/*     - (rfldn(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'rfldn   ',(rfldn(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flup    ',(flup(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flx_ft  ',(flx_ft(k,l),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flx_fb  ',(flx_fb(k,l),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flx_rt  ',(flx_rt(k,l),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flx_rb  ',(flx_rb(k,l),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flx_rb  ',(flx_dd(k,l),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'flx_rb  ',(flx_du(k,l),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'refflx_b', */
/*     -                             (refflx_b(k,l,ibn),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'trnflx_b', */
/*     -                             (trnflx_b(k,l,ibn),k=1,nlev) */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    return 0;
} /* flx_src__ */

