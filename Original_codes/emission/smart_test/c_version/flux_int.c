/* flux_int.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int flux_int__(logical *lsolar, logical *lplanck, integer *
	nlyr, integer *nz0, integer *ninter, doublereal *delnu, real *
	dir_s_flx__, real *dn_s_flx__, real *up_s_flx__, real *dn_t_flx__, 
	real *up_t_flx__, doublereal *dirsoflx, doublereal *dnsoflx, 
	doublereal *upsoflx, doublereal *dnthflx, doublereal *upthflx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, nz;


/* ccccccccccccccccccccccccc   f l u x _ i n t   ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine integrates level-dependent solar and thermal    cc */
/* c    fluxes over wavenumber.  These spectrally-integrated values     cc */
/* c    are used to find heating rates.                                 cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c       nlyr - number of computational model layers                  cc */
/* c        nzo - index of current solar zenith angle                   cc */
/* c      delnu - current spectral resolution (cm**-1)                  cc */
/* c     ninter - counter that records first time routine is called     cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    dir_s_flx - direct downward solar flux at wavenumber wn         cc */
/* c     dn_s_flx - total downward solar flux at wavenumber wn          cc */
/* c     up_s_flx - upward solar flux (irradiance) at wavenumber wn     cc */
/* c     dn_t_flx - downward thermal flux (irradiance) at wavenumber wn cc */
/* c     up_t_flx - upward thermal flux (irradiance) at wavenumber wn   cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc   f l u x _ i n t   ccccccccccccccccccccccccc */




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




/* ****   spectrally-dependent output flux and radiance */


/* ****   spectrally dependent fluxes */


/* ****   double precision integration variables for heating rates */


/* ****   specify the solar zenith angle index */

    /* Parameter adjustments */
    --upthflx;
    --dnthflx;
    upsoflx -= 71;
    dnsoflx -= 71;
    dirsoflx -= 71;
    up_t_flx__ -= 71;
    dn_t_flx__ -= 71;
    up_s_flx__ -= 771;
    dn_s_flx__ -= 771;
    dir_s_flx__ -= 771;

    /* Function Body */
    nz = *nz0;

    if (*ninter > 1) {

/* ****    add contributions to solar heating rates: */

	if (*lsolar) {

/* ****         o u t p u t    l e v e l   l o o p */

	    i__1 = *nlyr + 1;
	    for (k = 1; k <= i__1; ++k) {

/* ****          define spectrally-integrated solar fluxes for */
/*              heating rate calculation. */

		dirsoflx[k + nz * 70] += *delnu * .5f * (dir_s_flx__[k + (nz 
			+ 20) * 70] + dir_s_flx__[k + (nz + 10) * 70]);
		dir_s_flx__[k + (nz + 10) * 70] = dir_s_flx__[k + (nz + 20) * 
			70];
		dnsoflx[k + nz * 70] += *delnu * .5f * (dn_s_flx__[k + (nz + 
			20) * 70] + dn_s_flx__[k + (nz + 10) * 70]);
		dn_s_flx__[k + (nz + 10) * 70] = dn_s_flx__[k + (nz + 20) * 
			70];
		upsoflx[k + nz * 70] += *delnu * .5f * (up_s_flx__[k + (nz + 
			20) * 70] + up_s_flx__[k + (nz + 10) * 70]);
		up_s_flx__[k + (nz + 10) * 70] = up_s_flx__[k + (nz + 20) * 
			70];

/* L7001: */
	    }

	}

/* ****    add contributions to thermal cooling rates: */

	if (*lplanck && nz == 1) {

/* ****      define spectrally-integrated thermal fluxes for */
/*          heating rate calculation. */

	    i__1 = *nlyr + 1;
	    for (k = 1; k <= i__1; ++k) {
		dnthflx[k] += *delnu * .5f * (dn_t_flx__[k + 140] + 
			dn_t_flx__[k + 70]);
		dn_t_flx__[k + 70] = dn_t_flx__[k + 140];
		upthflx[k] += *delnu * .5f * (up_t_flx__[k + 140] + 
			up_t_flx__[k + 70]);
		up_t_flx__[k + 70] = up_t_flx__[k + 140];
/* L4541: */
	    }
	}
    } else {

/* *****     initialize flux integration arrays */

	i__1 = *nlyr + 1;
	for (k = 1; k <= i__1; ++k) {
	    up_s_flx__[k + (nz + 10) * 70] = up_s_flx__[k + (nz + 20) * 70];
	    dn_s_flx__[k + (nz + 10) * 70] = dn_s_flx__[k + (nz + 20) * 70];
	    dir_s_flx__[k + (nz + 10) * 70] = dir_s_flx__[k + (nz + 20) * 70];
/* L4601: */
	}
	i__1 = *nlyr + 1;
	for (k = 1; k <= i__1; ++k) {
	    dn_t_flx__[k + 70] = dn_t_flx__[k + 140];
	    up_t_flx__[k + 70] = up_t_flx__[k + 140];
/* L4621: */
	}
    }

    return 0;
} /* flux_int__ */

