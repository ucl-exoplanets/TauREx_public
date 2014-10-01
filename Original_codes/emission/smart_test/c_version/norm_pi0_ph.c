/* norm_pi0_ph.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int norm_pi0_ph__(integer *nmommx, integer *nlay, integer *
	nt_pd__, real *small_tau__, real *dtauex, real *dtausc, real *
	co_pi0__, real *g, real *phmom)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, l, mom;


/* cccccccccccccccccccccccc  n o r m _ pi0_ph  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine normalizes the single scattering albedo and     cc */
/* c    scattering phase function in each atmospheric layer.            cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlay - number of layers in the model atmosphere.             cc */
/* c     nmommx - maximum number of legendre polynomial moments         cc */
/* c      nt_pd - number of temperature profiles needed                 cc */
/* c              1 - radiances only, 2 - partial derivaties            cc */
/* c  small_tau - very small optical depth (small_tau*small_tau)        cc */
/* c     dtauex - differential extintion optical depth in each layer    cc */
/* c     dtausc - differential scattering optical depth in each layer   cc */
/* c     co_pi0 - single scattering co-albedo in each layer             cc */
/* c      phmom - scattering phase function moments in each layer.      cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      phmom - scattering phase function moments in each layer.      cc */
/* c     co_pi0 - single scattering co-albedo in each layer             cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccc  n o r m _ pi0_ph  ccccccccccccccccccccccccccc */




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




/* *****   monochromatic optical properties */


/* ****   normalize the effective asymmetry parameter */
/*       and moments of phase function. */

/*       NOTE:  The values of these parameters must be normalized by */
/*              dtausc rather than dtauex because the rayleigh scattering */
/*              contribution is explicitly included in these values. */

    /* Parameter adjustments */
    phmom -= 201;
    --g;
    co_pi0__ -= 71;
    --dtausc;
    dtauex -= 71;

    /* Function Body */
    i__1 = *nlay;
    for (k = 1; k <= i__1; ++k) {
	if (dtausc[k] > *small_tau__) {
	    g[k] /= dtausc[k];
	    i__2 = *nmommx;
	    for (mom = 0; mom <= i__2; ++mom) {
		phmom[mom + k * 201] /= dtausc[k];
/* L4401: */
	    }
	} else {
	    g[k] = 0.f;
	    phmom[k * 201] = 1.f;
	    i__2 = *nmommx;
	    for (mom = 1; mom <= i__2; ++mom) {
		phmom[mom + k * 201] = 0.f;
/* L4421: */
	    }
	}
/* L4441: */
    }

/* ****   define the single scattering CO-albedo (1. - pi0) and */
/*       normalize the asymmetry parameter.  Note: the asymmetry */
/*       parameter is normalized by */

    i__1 = *nt_pd__;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *nlay;
	for (k = 1; k <= i__2; ++k) {
	    if (dtauex[k + l * 70] > *small_tau__) {
		co_pi0__[k + l * 70] = 1.f - dtausc[k] / dtauex[k + l * 70];
	    } else {
		co_pi0__[k + l * 70] = 1.f;
	    }
/* L4601: */
	}
/* L4621: */
    }

    return 0;
} /* norm_pi0_ph__ */

