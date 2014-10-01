/* skip_wn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int skip_wn__(integer *iwnmode, integer *ij0, integer *nlay, 
	integer *iwnflg, integer *ntau_pd__, integer *numu, integer *irad, 
	integer *nza, doublereal *dwnmx, real *small_tau__, real *small_pi0__,
	 real *small_g__, real *small_alb__, real *tauerr, real *pi0err, real 
	*phferr, real *surferr, doublereal *wn, real *dtauex, real *dtausc, 
	real *co_pi0__, real *g, real *phmom, real *surf_opt__, real *alb, 
	doublereal *wni, real *dtauexi, real *dtausci, real *co_pi0i__, real *
	gi, real *phmomi, real *surf_opti__, real *albi, real *tau_ext__, 
	real *tau_sca__, real *g_sca__, real *tau_exti__, real *tau_scai__, 
	real *g_scai__, real *p_ray_0__, real *p_ray__, real *tau_ray__, real 
	*p_ray_0i__, real *p_rayi__, real *tau_rayi__, real *p_gas_0__, real *
	p_gas__, real *tau_gas__, real *p_gas_0i__, real *p_gasi__, real *
	tau_gasi__, real *p_aer_0__, real *p_aer__, real *tau_aer__, real *
	p_aer_0i__, real *p_aeri__, real *tau_aeri__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, l, m, n, nz, mom, nze;
    static doublereal diff, dwni;


/* cccccccccccccccccccccccc    s k i p _ w n    cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine determines whether the optical properties at    cc */
/* c    a spectral grid point will make a contribution to the spectrum, cc */
/* c    when compared to its neighboring points.  If not, it is skipped.cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                      nza                              cc */
/* c    iwnmode - operating mode:                                       cc */
/* c              (0) setup                                             cc */
/* c              (1) check to see if point x(2) is needed              cc */
/* c              (2) all 3 points are good                             cc */
/* c              (3) overwrite point x(2) with point x(3)              cc */
/* c        ij0 - wn index counter (1,2,3)                              cc */
/* c       nlay - number of levels (nlev or nlev+1 if pressure is a     cc */
/* c              variable part of the state vector                     cc */
/* c     iwnflg - flag stating if this point is needed:                 cc */
/* c              (0) no, (1) yes                                       cc */
/* c      dwnmx -  maximum spectral interval between output points      cc */
/* c  small_tau - small optical depth                                   cc */
/* c  small_pi0 - small single scattering albedo                        cc */
/* c    small_g - small asymmetry parameter                             cc */
/* c     tauerr - optical depth relative binning error (0. to ~0.8)     cc */
/* c     pi0err - co-single scattering albedo absolute binning error    cc */
/* c     phferr - asymmetry factor absolute binning error               cc */
/* c    surferr - surface optical property binning error                cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    input values, reset to skip values at specified wavenumber      cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccc    s k i p _ w n    cccccccccccccccccccccccccc */




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



/* ****   atmospheric optical properties */


/* ****   atmospheric scattering optical depth and surface properties */



/* ****    pressure of tau = 1 and column optical depths */


/* ****    optical depths of each variable component of state vector */







/* ****    check to see if this is a setup step (iwnmode = 0), */
/*        a wavelength skip test is to be performed (iwnmode = 1), */
/*        or a wavelength skip step is performed (iwnmode = 2) */

    /* Parameter adjustments */
    --tau_aeri__;
    p_aeri__ -= 17;
    --p_aer_0i__;
    --p_aer__;
    --tau_gasi__;
    p_gasi__ -= 17;
    --p_gas_0i__;
    --p_gas__;
    --tau_rayi__;
    p_rayi__ -= 17;
    --p_ray_0i__;
    --p_ray__;
    g_scai__ -= 771;
    tau_scai__ -= 771;
    tau_exti__ -= 2381;
    g_sca__ -= 71;
    tau_sca__ -= 71;
    tau_ext__ -= 281;
    albi -= 11;
    surf_opti__ -= 5;
    phmomi -= 14271;
    gi -= 71;
    co_pi0i__ -= 211;
    dtausci -= 71;
    dtauexi -= 211;
    --wni;
    --alb;
    --surf_opt__;
    phmom -= 201;
    --g;
    co_pi0__ -= 71;
    --dtausc;
    dtauex -= 71;

    /* Function Body */
    if (*iwnmode == 1) {

/* ****       determine if the optical properties are changing */
/*           rapidly enough to justify an RT calculation for this point */

	wni[*ij0] = *wn;

	if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	    p_ray_0i__[*ij0] = *p_ray_0__;
	    p_gas_0i__[*ij0] = *p_gas_0__;
	    p_aer_0i__[*ij0] = *p_aer_0__;
	    tau_rayi__[*ij0] = *tau_ray__;
	    tau_gasi__[*ij0] = *tau_gas__;
	    tau_aeri__[*ij0] = *tau_aer__;
	    i__1 = *numu;
	    for (nze = 1; nze <= i__1; ++nze) {
		p_rayi__[nze + (*ij0 << 4)] = p_ray__[nze];
		p_gasi__[nze + (*ij0 << 4)] = p_gas__[nze];
		p_aeri__[nze + (*ij0 << 4)] = p_aer__[nze];
/* L2001: */
	    }

	}

	i__1 = *nlay;
	for (k = 1; k <= i__1; ++k) {
	    dtauexi[k + ((*ij0 << 1) + 1) * 70] = dtauex[k + 70];
	    dtauexi[k + ((*ij0 << 1) + 2) * 70] = dtauex[k + 140];
	    dtausci[k + *ij0 * 70] = dtausc[k];
	    co_pi0i__[k + ((*ij0 << 1) + 1) * 70] = co_pi0__[k + 70];
	    co_pi0i__[k + ((*ij0 << 1) + 2) * 70] = co_pi0__[k + 140];
	    gi[k + *ij0 * 70] = g[k];
	    for (mom = 0; mom <= 200; ++mom) {
		phmomi[mom + (k + *ij0 * 70) * 201] = phmom[mom + k * 201];
/* L2021: */
	    }
	    i__2 = *ntau_pd__;
	    for (n = 1; n <= i__2; ++n) {
		tau_exti__[k + ((n + *ij0 * 10) * 3 + 1) * 70] = tau_ext__[k 
			+ (n * 3 + 1) * 70];
		tau_exti__[k + ((n + *ij0 * 10) * 3 + 2) * 70] = tau_ext__[k 
			+ (n * 3 + 2) * 70];
		tau_exti__[k + ((n + *ij0 * 10) * 3 + 3) * 70] = tau_ext__[k 
			+ (n * 3 + 3) * 70];
		tau_scai__[k + (n + *ij0 * 10) * 70] = tau_sca__[k + n * 70];
		g_scai__[k + (n + *ij0 * 10) * 70] = g_sca__[k + n * 70];
/* L2041: */
	    }
/* L2061: */
	}
	for (l = 1; l <= 4; ++l) {
	    surf_opti__[l + (*ij0 << 2)] = surf_opt__[l];
/* L2081: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    albi[nz + *ij0 * 10] = alb[nz];
/* L2091: */
	}

	if (*ij0 < 3) {
	    return 0;
	}

	*iwnflg = 0;

	if (wni[3] - wni[1] > *dwnmx) {
	    *iwnflg = 1;
	}
	dwni = (wni[2] - wni[1]) / (wni[3] - wni[1]);
	i__1 = *nlay;
	for (k = 1; k <= i__1; ++k) {
	    diff = dtauexi[k + 350] - (dtauexi[k + 210] + (dtauexi[k + 490] - 
		    dtauexi[k + 210]) * dwni);
	    if (abs(diff) > *tauerr * .01f * (dtauexi[k + 350] + *small_tau__)
		    ) {
		*iwnflg = 1;
	    }

	    diff = dtauexi[k + 420] - (dtauexi[k + 280] + (dtauexi[k + 560] - 
		    dtauexi[k + 280]) * dwni);
	    if (abs(diff) > *tauerr * .01f * (dtauexi[k + 420] + *small_tau__)
		    ) {
		*iwnflg = 1;
	    }

	    diff = dtausci[k + 140] - (dtausci[k + 70] + (dtausci[k + 210] - 
		    dtausci[k + 70]) * dwni);
	    if (abs(diff) > *tauerr * .01f * (dtausci[k + 140] + *small_tau__)
		    ) {
		*iwnflg = 1;
	    }

	    diff = co_pi0i__[k + 350] - (co_pi0i__[k + 210] + (co_pi0i__[k + 
		    490] - co_pi0i__[k + 210]) * dwni);
	    if (abs(diff) > *pi0err * .01f * (co_pi0i__[k + 350] + *
		    small_pi0__)) {
		*iwnflg = 1;
	    }

	    diff = co_pi0i__[k + 420] - (co_pi0i__[k + 280] + (co_pi0i__[k + 
		    560] - co_pi0i__[k + 280]) * dwni);
	    if (abs(diff) > *pi0err * .01f * (co_pi0i__[k + 420] + *
		    small_pi0__)) {
		*iwnflg = 1;
	    }

	    diff = gi[k + 140] - (gi[k + 70] + (gi[k + 210] - gi[k + 70]) * 
		    dwni);
	    if (abs(diff) > *phferr * .01f * (gi[k + 140] + *small_g__)) {
		*iwnflg = 1;
	    }

	    for (mom = 0; mom <= 200; ++mom) {
		diff = phmomi[mom + (k + 140) * 201] - (phmomi[mom + (k + 70) 
			* 201] + (phmomi[mom + (k + 210) * 201] - phmomi[mom 
			+ (k + 70) * 201]) * dwni);
		if (abs(diff) > *phferr * .01f * (phmomi[mom + (k + 140) * 
			201] + *small_g__)) {
		    *iwnflg = 1;
		}
/* L3001: */
	    }
/* L3021: */
	}
	for (l = 1; l <= 4; ++l) {
	    surf_opti__[l + (*ij0 << 2)] = surf_opt__[l];
	    diff = surf_opti__[l + 8] - (surf_opti__[l + 4] + (surf_opti__[l 
		    + 12] - surf_opti__[l + 4]) * dwni);
	    if (abs(diff) > *surferr * .01f * (surf_opti__[l + 8] + *
		    small_alb__)) {
		*iwnflg = 1;
	    }
/* L3041: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    albi[nz + *ij0 * 10] = alb[nz];
	    diff = albi[nz + 20] - (albi[nz + 10] + (albi[nz + 30] - albi[nz 
		    + 10]) * dwni);
	    if (abs(diff) > *surferr * .01f * (albi[l + 20] + *small_alb__)) {
		*iwnflg = 1;
	    }
/* L3061: */
	}

    } else {

	if (*iwnmode == 2) {

/* *****       All 3 points are good. put x(1) into working arrays */
/*            and process that point, and then pack x(2) into x(1) */
/*            and x(3) into x(2) */

	    *wn = wni[1];
	    wni[1] = wni[2];
	    wni[2] = wni[3];

	    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
		*p_ray_0__ = p_ray_0i__[1];
		p_ray_0i__[1] = p_ray_0i__[2];
		p_ray_0i__[2] = p_ray_0i__[3];
		*p_gas_0__ = p_gas_0i__[1];
		p_gas_0i__[1] = p_gas_0i__[2];
		p_gas_0i__[2] = p_gas_0i__[3];
		*p_aer_0__ = p_aer_0i__[1];
		p_aer_0i__[1] = p_aer_0i__[2];
		p_aer_0i__[2] = p_aer_0i__[3];
		*tau_ray__ = tau_rayi__[1];
		tau_rayi__[1] = tau_rayi__[2];
		tau_rayi__[2] = tau_rayi__[3];
		*tau_gas__ = tau_gasi__[1];
		tau_gasi__[1] = tau_gasi__[2];
		tau_gasi__[2] = tau_gasi__[3];
		*tau_aer__ = tau_aeri__[1];
		tau_aeri__[1] = tau_aeri__[2];
		tau_aeri__[2] = tau_aeri__[3];
		i__1 = *numu;
		for (nze = 1; nze <= i__1; ++nze) {
		    p_ray__[nze] = p_rayi__[nze + 16];
		    p_rayi__[nze + 16] = p_rayi__[nze + 32];
		    p_rayi__[nze + 32] = p_rayi__[nze + 48];
		    p_gas__[nze] = p_gasi__[nze + 16];
		    p_gasi__[nze + 16] = p_gasi__[nze + 32];
		    p_gasi__[nze + 32] = p_gasi__[nze + 48];
		    p_aer__[nze] = p_aeri__[nze + 16];
		    p_aeri__[nze + 16] = p_aeri__[nze + 32];
		    p_aeri__[nze + 32] = p_aeri__[nze + 48];
/* L4001: */
		}

	    }

	    i__1 = *nlay;
	    for (k = 1; k <= i__1; ++k) {
		dtauex[k + 70] = dtauexi[k + 210];
		dtauexi[k + 210] = dtauexi[k + 350];
		dtauexi[k + 350] = dtauexi[k + 490];
		dtauex[k + 140] = dtauexi[k + 280];
		dtauexi[k + 280] = dtauexi[k + 420];
		dtauexi[k + 420] = dtauexi[k + 560];
		dtausc[k] = dtausci[k + 70];
		dtausci[k + 70] = dtausci[k + 140];
		dtausci[k + 140] = dtausci[k + 210];
		co_pi0__[k + 70] = co_pi0i__[k + 210];
		co_pi0i__[k + 210] = co_pi0i__[k + 350];
		co_pi0i__[k + 350] = co_pi0i__[k + 490];
		co_pi0__[k + 140] = co_pi0i__[k + 280];
		co_pi0i__[k + 280] = co_pi0i__[k + 420];
		co_pi0i__[k + 420] = co_pi0i__[k + 560];
		g[k] = gi[k + 70];
		gi[k + 70] = gi[k + 140];
		gi[k + 140] = gi[k + 210];
		for (mom = 0; mom <= 200; ++mom) {
		    phmom[mom + k * 201] = phmomi[mom + (k + 70) * 201];
		    phmomi[mom + (k + 70) * 201] = phmomi[mom + (k + 140) * 
			    201];
		    phmomi[mom + (k + 140) * 201] = phmomi[mom + (k + 210) * 
			    201];
/* L4021: */
		}
		i__2 = *ntau_pd__;
		for (n = 1; n <= i__2; ++n) {
		    for (m = 1; m <= 3; ++m) {
			tau_ext__[k + (m + n * 3) * 70] = tau_exti__[k + (m + 
				(n + 10) * 3) * 70];
			tau_exti__[k + (m + (n + 10) * 3) * 70] = tau_exti__[
				k + (m + (n + 20) * 3) * 70];
			tau_exti__[k + (m + (n + 20) * 3) * 70] = tau_exti__[
				k + (m + (n + 30) * 3) * 70];
/* L4031: */
		    }
		    tau_sca__[k + n * 70] = tau_scai__[k + (n + 10) * 70];
		    tau_scai__[k + (n + 10) * 70] = tau_scai__[k + (n + 20) * 
			    70];
		    tau_scai__[k + (n + 20) * 70] = tau_scai__[k + (n + 30) * 
			    70];
		    g_sca__[k + n * 70] = g_scai__[k + (n + 10) * 70];
		    g_scai__[k + (n + 10) * 70] = g_scai__[k + (n + 20) * 70];
		    g_scai__[k + (n + 20) * 70] = g_scai__[k + (n + 30) * 70];
/* L4041: */
		}
/* L4061: */
	    }
	    for (l = 1; l <= 4; ++l) {
		surf_opt__[l] = surf_opti__[l + 4];
		surf_opti__[l + 4] = surf_opti__[l + 8];
		surf_opti__[l + 8] = surf_opti__[l + 12];
/* L4081: */
	    }
	    i__1 = *nza;
	    for (nz = 1; nz <= i__1; ++nz) {
		alb[nz] = albi[nz + 10];
		albi[nz + 10] = albi[nz + 20];
		albi[nz + 20] = albi[nz + 30];
/* L4091: */
	    }

	} else {

/* ****        iwnmode = 3: skip point x(2).  Load x(3) into x(2) */

	    wni[2] = wni[3];

	    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
		p_ray_0i__[2] = p_ray_0i__[3];
		p_gas_0i__[2] = p_gas_0i__[3];
		p_aer_0i__[2] = p_aer_0i__[3];
		tau_rayi__[2] = tau_rayi__[3];
		tau_gasi__[2] = tau_gasi__[3];
		tau_aeri__[2] = tau_aeri__[3];
		i__1 = *numu;
		for (nze = 1; nze <= i__1; ++nze) {
		    p_rayi__[nze + 32] = p_rayi__[nze + 48];
		    p_gasi__[nze + 32] = p_gasi__[nze + 48];
		    p_aeri__[nze + 32] = p_aeri__[nze + 48];
/* L5001: */
		}
	    }

	    i__1 = *nlay;
	    for (k = 1; k <= i__1; ++k) {
		dtauexi[k + 350] = dtauexi[k + 490];
		dtauexi[k + 420] = dtauexi[k + 560];
		dtausci[k + 140] = dtausci[k + 210];
		co_pi0i__[k + 350] = co_pi0i__[k + 490];
		co_pi0i__[k + 420] = co_pi0i__[k + 560];
		gi[k + 140] = gi[k + 210];
		for (mom = 0; mom <= 200; ++mom) {
		    phmomi[mom + (k + 140) * 201] = phmomi[mom + (k + 210) * 
			    201];
/* L5021: */
		}
		i__2 = *ntau_pd__;
		for (n = 1; n <= i__2; ++n) {
		    tau_exti__[k + ((n + 20) * 3 + 1) * 70] = tau_exti__[k + (
			    (n + 30) * 3 + 1) * 70];
		    tau_exti__[k + ((n + 20) * 3 + 2) * 70] = tau_exti__[k + (
			    (n + 30) * 3 + 2) * 70];
		    tau_exti__[k + ((n + 20) * 3 + 3) * 70] = tau_exti__[k + (
			    (n + 30) * 3 + 3) * 70];
		    tau_scai__[k + (n + 20) * 70] = tau_scai__[k + (n + 30) * 
			    70];
		    g_scai__[k + (n + 20) * 70] = g_scai__[k + (n + 30) * 70];
/* L5041: */
		}
/* L5061: */
	    }
	    for (l = 1; l <= 4; ++l) {
		surf_opti__[l + 8] = surf_opti__[l + 12];
/* L5081: */
	    }
	    i__1 = *nza;
	    for (nz = 1; nz <= i__1; ++nz) {
		albi[nz + 20] = albi[nz + 30];
/* L5091: */
	    }

	}

    }

    return 0;
} /* skip_wn__ */

