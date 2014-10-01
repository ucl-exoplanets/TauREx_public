/* init_bin.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int init_bin__(integer *nlyr, integer *nref, integer *nza, 
	integer *ngroup, integer *iwngrp, integer *ismterr, integer *itau1cnt,
	 integer *itaucnt, integer *ipi0cnt, integer *igcnt, integer *isurcnt,
	 integer *levtau1, integer *nmomgrp, doublereal *dnugrp, doublereal *
	wngrp, real *surfgrp, real *albgrp, real *taugrp, real *pi0grp, real *
	ggrp, real *pmomgrp, real *dtau_b__, real *copi0_b__, real *g_b__, 
	real *pmom_b__, real *sur_b__, real *alb_b__, real *dx_b_i__, real *
	dalb_b_i__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, l, m, n, nr, nz;


/* cccccccccccccccccccccc    i n i t _ b i n     ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c     p u r p o s e :                                                cc */
/* c                                                                    cc */
/* c    this subroutine initializes values set in the binning routine.  cc */
/* c                                                                    cc */
/* c     i n p u t                                                      cc */
/* c                                                                    cc */
/* c       nlyr - number of layers in nominal atmosphere                cc */
/* c       nref - number of surface optical properties specified at     cc */
/* c              each wavelength.                                      cc */
/* c     ngroup - number of spectral bins (0 to ngrp)                   cc */
/* c     iwngrp - index of spectral bin                                 cc */
/* c    ismterr - error flag: ngroup > ngrp                             cc */
/* c   itau1cnt - number of bins rejected by poor fit at tau=1          cc */
/* c    itaucnt - number of bins rejected by poor fit away from tau=1   cc */
/* c    ipi0cnt - number of bins rejected by poor pi0 fit               cc */
/* c      igcnt - number of bins rejected by poor g fit                 cc */
/* c    isurcnt - number of bins rejected by poor albedo fit            cc */
/* c    levtau1 - index of tau=1 level                                  cc */
/* c    nmomgrp - maximum number of phase function moments in bin       cc */
/* c     taugrp - mean, min, and max optical depth in each bin          cc */
/* c     piogrp - mean, min, and max single scattering albedo in  bin   cc */
/* c    pmomgrp - mean, min, and max phase function in each bin         cc */
/* c       ggrp - mean, min, and max asymmetry parmeter in each bin     cc */
/* c     dtau_b - perturbed differential optical depth in each bin      cc */
/* c      pi0_b - perturbed single scattering albedo in each bin        cc */
/* c     pmom_b - perturbed scattering phase function in each bin       cc */
/* c        g_b - perturbed scattering asymmetery parametr in each bin  cc */
/* c      sur_b - perturbed surface albedo in each bin                  cc */
/* c                                                                    cc */
/* c     o u t p u t                                                    cc */
/* c                                                                    cc */
/* c    zeroed variables.                                               cc */
/* c                                                                    cc */
/* cccccccccccccccccccccc    i n i t _ b i n     ccccccccccccccccccccccccc */




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






/* ****   initialize spectral binning parameters */

    /* Parameter adjustments */
    dalb_b_i__ -= 11;
    dx_b_i__ -= 421;
    alb_b__ -= 61;
    sur_b__ -= 25;
    pmom_b__ -= 84621;
    g_b__ -= 421;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    pmomgrp -= 7218111;
    ggrp -= 35911;
    pi0grp -= 35911;
    taugrp -= 35911;
    albgrp -= 5633;
    surfgrp -= 2561;
    wngrp -= 513;
    --dnugrp;
    --nmomgrp;
    levtau1 -= 513;

    /* Function Body */
    *ismterr = 0;
    *itau1cnt = 0;
    *itaucnt = 0;
    *ipi0cnt = 0;
    *igcnt = 0;
    *isurcnt = 0;
    *iwngrp = 0;
    *ngroup = 0;

    for (n = 1; n <= 512; ++n) {
	levtau1[n + 512] = *nlyr;
	levtau1[n + 1024] = *nlyr;
	nmomgrp[n] = 0;
	dnugrp[n] = 0.;

	for (l = 1; l <= 3; ++l) {
	    wngrp[n + (l << 9)] = 0.;
	    i__1 = *nref;
	    for (nr = 1; nr <= i__1; ++nr) {
		surfgrp[n + (nr + (l << 2) << 9)] = 0.f;
/* L2001: */
	    }
	    i__1 = *nza;
	    for (nz = 1; nz <= i__1; ++nz) {
		albgrp[n + (nz + l * 10 << 9)] = 0.f;
/* L2201: */
	    }
	    for (k = 1; k <= 70; ++k) {
		taugrp[k + (n + (l << 9)) * 70] = 0.f;
		pi0grp[k + (n + (l << 9)) * 70] = 0.f;
		ggrp[k + (n + (l << 9)) * 70] = 0.f;
		pmomgrp[(k + (n + (l << 9)) * 70) * 201] = 1.f;
		for (m = 1; m <= 200; ++m) {
		    pmomgrp[m + (k + (n + (l << 9)) * 70) * 201] = 0.f;
/* L2401: */
		}
/* L2601: */
	    }
/* L2801: */
	}

	for (l = 1; l <= 5; ++l) {
	    for (k = 1; k <= 70; ++k) {
		dtau_b__[k + (l + n * 5) * 70] = 0.f;
		copi0_b__[k + (l + n * 5) * 70] = 0.f;
		g_b__[k + (l + n * 5) * 70] = 0.f;
		pmom_b__[(k + (l + n * 5) * 70) * 201] = 1.f;
		dx_b_i__[k + (l + n * 5) * 70] = 0.f;
		for (m = 1; m <= 200; ++m) {
		    pmom_b__[m + (k + (l + n * 5) * 70) * 201] = 0.f;
/* L3001: */
		}
/* L3021: */
	    }
/* L3061: */
	}

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    dalb_b_i__[nz + n * 10] = 0.f;
/* L3101: */
	}

	for (l = 1; l <= 5; ++l) {
	    i__1 = *nref;
	    for (nr = 1; nr <= i__1; ++nr) {
		sur_b__[nr + (l + n * 5 << 2)] = 0.f;
/* L3201: */
	    }
	    i__1 = *nza;
	    for (nz = 1; nz <= i__1; ++nz) {
		alb_b__[nz + (l + n * 5) * 10] = 0.f;
/* L3221: */
	    }
/* L3241: */
	}
/* L3281: */
    }

    return 0;
} /* init_bin__ */

