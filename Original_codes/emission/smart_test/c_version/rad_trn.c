/* rad_trn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rad_trn__(logical *usrang, logical *lplanck, integer *
	l_1__, integer *l_2__, integer *ng0, integer *nstr, integer *numu, 
	real *umu, real *umu_f__, real *gwt_f__, real *dtau, real *copi0, 
	real *g0, doublereal *trnrad, doublereal *absrad, doublereal *
	trn_cmu__, doublereal *abs_cmu__)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer k, nze;
    static doublereal dtauabs, dtauexe;


/* cccccccccccccccccccccccccc  r a d _ t r n   ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the effective layer radiance and flux  cc */
/* c    transmittances and absorptances for each spectral bin.  These   cc */
/* c    quantities are used in the spectral mapping and jacobian        cc */
/* c    calculations.                                                   cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c     usrang: logical variable.  If false, use computational angles, cc */
/* c             if true, find radiances for user specified angles.     cc */
/* c        l_1: index of first layer for which values are needed       cc */
/* c        l_2: index of last layer for which values are needed        cc */
/* c        ng0: spectral bin number: 0=> single scat, >0=> mult scat   cc */
/* c       nstr: number of zenith angle streams.                        cc */
/* c       numu: number of output zenith angles.                        cc */
/* c       dtau: layer optical depth                                    cc */
/* c      copi0: layer single scattering albedo                         cc */
/* c         g0: layer asymmetry parameter                              cc */
/* c        umu: zenith angle of each stream.                           cc */
/* c      umu_f: gaussian point for each stream.                        cc */
/* c      gwt_f: gaussian weight for each stream.                       cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    trn_cmu: layer transmittance for each radiance stream.          cc */
/* c    abs_cmu: layer absorptance for each rediance stream             cc */
/* c     trnrad: layer transmittance for each user radiance stream.     cc */
/* c     absrad: layer absorptance for each user radiance stream.       cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  r a d _ t r n   ccccccccccccccccccccccccccc */




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





/* ***    set up direction cosines for flux integration */


/* *****   dummy state vector optical properties. */


/* ****    layer transmission values for simplified adding method */


/* ****    find the layer transittances and reflectances for each layer */

    /* Parameter adjustments */
    abs_cmu__ -= 17;
    trn_cmu__ -= 17;
    absrad -= 17;
    trnrad -= 17;
    --g0;
    --copi0;
    --dtau;
    --gwt_f__;
    --umu_f__;
    --umu;

    /* Function Body */
    i__1 = *l_2__;
    for (k = *l_1__; k <= i__1; ++k) {

/* ****       define the delta-m scaled extinction optical depth */

	dtauexe = dtau[k];
/*          dtauexe = dtau(k)*(1.0 - (1.0 - copi0(k))*g0(k)*g0(k)) */

/* ****       define the absorption optical depth */

	dtauabs = dtau[k] * copi0[k];

/*           define the direct transmission and absorption */
/*           along each computational stream */

	i__2 = *nstr / 2;
	for (nze = 1; nze <= i__2; ++nze) {
	    trn_cmu__[nze + *nstr / 2 + (k << 4)] = exp(-dtauexe / umu_f__[
		    nze + *nstr / 2]);
	    trn_cmu__[*nstr / 2 - nze + 1 + (k << 4)] = trn_cmu__[nze + *nstr 
		    / 2 + (k << 4)];
/* L2001: */
	}
	if (*ng0 == 0 || *lplanck) {
	    i__2 = *nstr / 2;
	    for (nze = 1; nze <= i__2; ++nze) {
		abs_cmu__[nze + *nstr / 2 + (k << 4)] = copi0[k] * (1. - exp(
			-dtauabs / umu_f__[nze + *nstr / 2]));
		abs_cmu__[*nstr / 2 - nze + 1 + (k << 4)] = abs_cmu__[nze + *
			nstr / 2 + (k << 4)];
/* L2021: */
	    }
	}

	if (! (*usrang)) {

/*             define the direct transmission along each stream */

	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		trnrad[nze + (k << 4)] = trn_cmu__[nze + (k << 4)];
		absrad[nze + (k << 4)] = abs_cmu__[nze + (k << 4)];
/* L2201: */
	    }

	} else {

/* ****         usrang = .true.  Recompute radiances on the */
/*             computational grid and integrate these over */
/*             angle to get fluxes */

/*             define the direct transmission along each stream */

	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		trnrad[nze + (k << 4)] = exp(-dtauexe / (r__1 = umu[nze], 
			dabs(r__1)));
		absrad[nze + (k << 4)] = copi0[k] * (1. - exp(-dtauabs / (
			r__1 = umu[nze], dabs(r__1))));
/* L2401: */
	    }

	}

/* L2621: */
    }

    return 0;
} /* rad_trn__ */

