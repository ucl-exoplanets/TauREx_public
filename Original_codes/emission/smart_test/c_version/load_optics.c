/* load_optics.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int load_optics__(logical *lamber, integer *ng0, integer *
	nlyr, integer *nza, integer *nstate, integer *istate, integer *
	n_rad__, integer *nref, integer *iref, integer *nmomgrp, real *umu0, 
	real *tauerr, real *pi0err, real *phferr, real *surferr, real *taumn, 
	real *surfgrp, real *albgrp, real *taugrp, real *pi0grp, real *ggrp, 
	real *pmomgrp, real *alb_b__, real *sur_b__, real *phiw, real *ws, 
	real *dtau_b__, real *copi0_b__, real *g_b__, real *pmom_b__, real *
	dx_b_i__, real *dalb_b_i__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static real deltau_b__, ss_range__;
    static integer k, n;
    static real tau_range__, deltaumin;
    static integer nr, nz, ibn, mom, ist;
    extern doublereal dref_(real *, real *, integer *);
    static real gmin, dref0, dsurf, delg_b__, dref_i__[10], ssamin, taumin, 
	    g_range__, delpi0_b__, surf_pr__[4], surfmin;


/* cccccccccccccccccccccc  l o a d _ o p t i c s   ccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this routine loads optical depths, single scattering albedos,   cc */
/* c    surface reflectance and surface pressures into arrays from use  cc */
/* c    in the routine sm_eq_trn.                                       cc */
/* c    NOTE: this version implements a perturbation that is            cc */
/* c          related to the actual range in each bin.                  cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        ng0 - index  of current spectral bin                        cc */
/* c       nlyr - number of computational model layers                  cc */
/* c     nstate - number of variable elements in the state vecror       cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c              1 - surface pressure                                  cc */
/* c              2 - surface/atmospheric temperature                   cc */
/* c              3 - gas absorption coeffient                          cc */
/* c              4 - cloud/aerosol optical depth                       cc */
/* c              5 - surface albedo                                    cc */
/* c      n_rad - radiance calculations flag for each bin               cc */
/*       n_rad = 1 : radiances for nominal state structure */
/*             = 2 : radiances for perturbed optical depths */
/*             = 3 : radiances for perturbed single scattering albedos */
/*             = 4 : radiances for perturbed phase functions */
/*             = 5 : radiances for perturbed surface reflectances */
/* c       nref - number of surface optical properties specified at     cc */
/* c              each wavelength.                                      cc */
/* c     tauerr - optical depth relative binning error (0. to ~0.8)     cc */
/* c     pi0err - co-single scattering albedo absolute binning error    cc */
/* c    surferr - surface optical property binning error                cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      sur_b - surface albedo for this bin                           cc */
/* c     dtau_b - differential optical depth in each layer for this bin cc */
/* c    copi0_b - single scattering co-albedo in each layer in this bin cc */
/* c        g_b - asymmetry parameter in each layer for this bin        cc */
/* c                                                                    cc */
/* cccccccccccccccccccccc  l o a d _ o p t i c s   ccccccccccccccccccccccc */




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





/* *****   state vector variables. */


/* ****   internal variables */



/* ****    set group counter */

    /* Parameter adjustments */
    dalb_b_i__ -= 11;
    dx_b_i__ -= 421;
    pmom_b__ -= 84621;
    g_b__ -= 421;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    sur_b__ -= 25;
    alb_b__ -= 61;
    pmomgrp -= 7218111;
    ggrp -= 35911;
    pi0grp -= 35911;
    taugrp -= 35911;
    albgrp -= 5633;
    surfgrp -= 2561;
    --umu0;
    --nmomgrp;
    n_rad__ -= 6;
    --istate;

    /* Function Body */
    ibn = *ng0;

/* ****    define a minimum change for delta-tau */

/*       deltaumin = 1.e-5 */
    deltaumin = .001f;

/* ****   set optical depths for the nomical case and for each variable */
/*       component of the optical depth structure and set */
/*       the radiance calculation flag, n_rad, where */
/*       n_rad = 1 : radiances for nominal state structure */
/*             = 2 : radiances for perturbed optical depths */
/*             = 3 : radiances for perturbed single scattering albedos */
/*             = 4 : radiances for perturbed phase functions */
/*             = 5 : radiances for perturbed surface reflectances */

/* ****       define the column-integrated optcal depth, single scatting */
/*           albedo, and particle phase function ranges in each bin */

    tau_range__ = 0.f;
    ss_range__ = 0.f;
    g_range__ = 0.f;

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

/* ****       define the column-integrated optcal depth, single scatting */
/*           albedo, and particle phase function ranges */

	tau_range__ = tau_range__ + taugrp[k + (ibn + 1536) * 70] - taugrp[k 
		+ (ibn + 1024) * 70];
	ss_range__ = ss_range__ + pi0grp[k + (ibn + 1536) * 70] - pi0grp[k + (
		ibn + 1024) * 70];
	g_range__ = g_range__ + ggrp[k + (ibn + 1536) * 70] - ggrp[k + (ibn + 
		1024) * 70];
/* L1001: */
    }

    if (istate[1] == 1) {

/* ****      add the difference between the nominal surface layer */
/*          optical depth and the optical depth for the layer with */
/*          a perturbed surface pressure. */

	tau_range__ += (r__1 = taugrp[*nlyr + 1 + (ibn + 512) * 70] - taugrp[*
		nlyr + (ibn + 512) * 70], dabs(r__1));
	ss_range__ += (r__1 = pi0grp[*nlyr + 1 + (ibn + 512) * 70] - pi0grp[*
		nlyr + (ibn + 512) * 70], dabs(r__1));
	g_range__ += (r__1 = ggrp[*nlyr + 1 + (ibn + 512) * 70] - ggrp[*nlyr 
		+ (ibn + 512) * 70], dabs(r__1));

    }

/* ****     n o m i n a l    c a s e */

    n_rad__[ibn * 5 + 1] = 1;

/* ****    define the optical properties for the nominal state vector */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	dtau_b__[k + (ibn * 5 + 1) * 70] = taugrp[k + (ibn + 512) * 70];
	copi0_b__[k + (ibn * 5 + 1) * 70] = pi0grp[k + (ibn + 512) * 70];
	g_b__[k + (ibn * 5 + 1) * 70] = ggrp[k + (ibn + 512) * 70];
	pmom_b__[(k + (ibn * 5 + 1) * 70) * 201] = 1.f;
	i__2 = nmomgrp[ibn];
	for (mom = 1; mom <= i__2; ++mom) {
	    pmom_b__[mom + (k + (ibn * 5 + 1) * 70) * 201] = pmomgrp[mom + (k 
		    + (ibn + 512) * 70) * 201];
/* L1011: */
	}
/* L1021: */
    }

/* ****   define the nominal surface albedo */

    i__1 = *nref;
    for (nr = 1; nr <= i__1; ++nr) {
	sur_b__[nr + (ibn * 5 + 1 << 2)] = surfgrp[ibn + (nr + 4 << 9)];
/* L1041: */
    }
    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	alb_b__[nz + (ibn * 5 + 1) * 10] = albgrp[ibn + (nz + 10 << 9)];
/* L1061: */
    }

/* ****     o p t i c a l   p r o p e r t y    p e r t u r b a t i o n s */

/* ***      check to see if optical depth or single scattering albedo */
/*         are variable parts of the state structure */

    ist = 0;
    i__1 = *nstate;
    for (n = 1; n <= i__1; ++n) {
	if ((i__2 = istate[n], abs(i__2)) >= 1 && (i__3 = istate[n], abs(i__3)
		) <= 4) {
	    ist = 1;
	}
/* L1301: */
    }

/* ****   determine if optical depth varies across the bin */

    taumin = *tauerr * .001f * *tauerr + *taumn;
/*      if(ibn .eq. 1) write(*,'(/,1a,2(1pe12.4))') */
/*     - 'load_optics: taumin,tau_range: ', */
/*     - taumin,tau_range */
    if (tau_range__ > taumin || ist == 1) {

	n_rad__[ibn * 5 + 2] = 1;

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {

/* ****         note: add at least deltaumin to this quantity to */
/*                   ensure that radiance difference is no-zero, and */
/*                   that the jacobian denominator is not too small */

	    dtau_b__[k + (ibn * 5 + 2) * 70] = deltaumin + (*tauerr * *tauerr 
		    + 1.f) * taugrp[k + (ibn + 512) * 70] + (taugrp[k + (ibn 
		    + 1536) * 70] - taugrp[k + (ibn + 1024) * 70]) * .25f;
	    copi0_b__[k + (ibn * 5 + 2) * 70] = pi0grp[k + (ibn + 512) * 70];
	    g_b__[k + (ibn * 5 + 2) * 70] = ggrp[k + (ibn + 512) * 70];
	    pmom_b__[(k + (ibn * 5 + 2) * 70) * 201] = 1.f;
	    i__2 = nmomgrp[ibn];
	    for (mom = 1; mom <= i__2; ++mom) {
		pmom_b__[mom + (k + (ibn * 5 + 2) * 70) * 201] = pmomgrp[mom 
			+ (k + (ibn + 512) * 70) * 201];
/* L1321: */
	    }

/* L1351: */
	}
/*      if(ibn .eq. 1) then */
/*     - write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dtau_b(1)   ', */
/*     - (dtau_b(k,1,ibn),k=1,nlyr) */
/*     - write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dtau_b(2)   ', */
/*     - (dtau_b(k,2,ibn),k=1,nlyr) */
/*       endif */

	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    sur_b__[nr + (ibn * 5 + 2 << 2)] = surfgrp[ibn + (nr + 4 << 9)];
/* L1341: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb_b__[nz + (ibn * 5 + 2) * 10] = albgrp[ibn + (nz + 10 << 9)];
/* L1361: */
	}

    } else {

	n_rad__[ibn * 5 + 2] = 0;

    }

/* ****     s i n g l e   s c a t t e r i n g   a l b e d o */

/* ****   check single scattering co-albedo range in this bin */

    ssamin = *pi0err * .001f * *pi0err + 1e-5f;

    if (ss_range__ > ssamin || ist == 1) {

	n_rad__[ibn * 5 + 3] = 1;

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    dtau_b__[k + (ibn * 5 + 3) * 70] = taugrp[k + (ibn + 512) * 70];

	    copi0_b__[k + (ibn * 5 + 3) * 70] = (*pi0err * *pi0err + 1.f) * 
		    pi0grp[k + (ibn + 512) * 70] + .01f + (pi0grp[k + (ibn + 
		    1536) * 70] - pi0grp[k + (ibn + 1024) * 70]) * .25f;

	    if (copi0_b__[k + (ibn * 5 + 3) * 70] > 1.f) {
		copi0_b__[k + (ibn * 5 + 3) * 70] = pi0grp[k + (ibn + 512) * 
			70] / (*pi0err * *pi0err + 1.f) - .01f - (pi0grp[k + (
			ibn + 1536) * 70] - pi0grp[k + (ibn + 1024) * 70]) * 
			.25f;
	    }

	    g_b__[k + (ibn * 5 + 3) * 70] = ggrp[k + (ibn + 512) * 70];
	    pmom_b__[(k + (ibn * 5 + 3) * 70) * 201] = 1.f;
	    i__2 = nmomgrp[ibn];
	    for (mom = 1; mom <= i__2; ++mom) {
		pmom_b__[mom + (k + (ibn * 5 + 3) * 70) * 201] = pmomgrp[mom 
			+ (k + (ibn + 512) * 70) * 201];
/* L1401: */
	    }
/* L1411: */
	}

	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    sur_b__[nr + (ibn * 5 + 3 << 2)] = surfgrp[ibn + (nr + 4 << 9)];
/* L1421: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb_b__[nz + (ibn * 5 + 3) * 10] = albgrp[ibn + (nz + 10 << 9)];
/* L1441: */
	}

    } else {

	n_rad__[ibn * 5 + 3] = 0;

    }

/* ****     s c a t t e r i n g   p h a s e    f u n c t i o n */

/* *****      NOTE: phase functions cannot be arbitrarily multiplied by */
/*           a constant, or they will lose their normaization. */
/*           replace the values by the min and max binned values */

/* ****   check phase function range in this bin */

    gmin = *phferr * .001f * *phferr + 1e-5f;
    if (g_range__ > gmin || ist == 1) {

	n_rad__[ibn * 5 + 4] = 1;

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    dtau_b__[k + (ibn * 5 + 4) * 70] = taugrp[k + (ibn + 512) * 70];
	    copi0_b__[k + (ibn * 5 + 4) * 70] = pi0grp[k + (ibn + 512) * 70];
	    g_b__[k + (ibn * 5 + 4) * 70] = ggrp[k + (ibn + 1536) * 70];
	    pmom_b__[(k + (ibn * 5 + 4) * 70) * 201] = 1.f;
	    i__2 = nmomgrp[ibn];
	    for (mom = 1; mom <= i__2; ++mom) {
		pmom_b__[mom + (k + (ibn * 5 + 4) * 70) * 201] = pmomgrp[mom 
			+ (k + (ibn + 1536) * 70) * 201];
/* L1521: */
	    }
/* L1541: */
	}

	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    sur_b__[nr + (ibn * 5 + 4 << 2)] = surfgrp[ibn + (nr + 4 << 9)];
/* L1561: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb_b__[nz + (ibn * 5 + 4) * 10] = albgrp[ibn + (nz + 10 << 9)];
/* L1581: */
	}

    } else {

	n_rad__[ibn * 5 + 4] = 0;

    }

/* ****    s u r f a c e    o p t i c a l   p r o p e r t i e s */

/* ****   in this version, the surface reflectance is always */
/*       a variable part of the state structure if any variable changes */
/*       because it affects the upward radiances at the surface */

    ist = 0;
    i__1 = *nstate;
    for (n = 1; n <= i__1; ++n) {
	if ((i__2 = istate[n], abs(i__2)) == 5) {
	    ist = 1;
	}
/* L1601: */
    }

/* *****  set the optical depths and single scattering albedos to nominal */
/*       values */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	dtau_b__[k + (ibn * 5 + 5) * 70] = taugrp[k + (ibn + 512) * 70];
	copi0_b__[k + (ibn * 5 + 5) * 70] = pi0grp[k + (ibn + 512) * 70];
	g_b__[k + (ibn * 5 + 5) * 70] = ggrp[k + (ibn + 512) * 70];
	pmom_b__[(k + (ibn * 5 + 5) * 70) * 201] = 0.f;
	i__2 = nmomgrp[ibn];
	for (mom = 1; mom <= i__2; ++mom) {
	    pmom_b__[mom + (k + (ibn * 5 + 5) * 70) * 201] = pmomgrp[mom + (k 
		    + (ibn + 512) * 70) * 201];
/* L1611: */
	}
/* L1621: */
    }

/* ****   determine the surface reflectance variations in this bin */

    n_rad__[ibn * 5 + 5] = 0;
    surfmin = *surferr * .001f * *surferr + 1e-5f;
    i__1 = *nref;
    for (nr = 1; nr <= i__1; ++nr) {

	dsurf = surfgrp[ibn + (nr + 12 << 9)] - surfgrp[ibn + (nr + 8 << 9)];

	if (dsurf > surfmin || ist == 1) {

/* ****        radiances must be calculated for 2 surface reflectances */

	    n_rad__[ibn * 5 + 5] = 1;

	    sur_b__[nr + (ibn * 5 + 5 << 2)] = dsurf * .25f + .01f + (*
		    surferr * *surferr + 1.f) * surfgrp[ibn + (nr + 4 << 9)];
	    if (sur_b__[nr + (ibn * 5 + 5 << 2)] > 1.f) {
		sur_b__[nr + (ibn * 5 + 5 << 2)] = surfgrp[ibn + (nr + 4 << 9)
			] * .99f - dsurf * .25f;
	    }

	} else {

	    sur_b__[nr + (ibn * 5 + 5 << 2)] = surfgrp[ibn + (nr + 4 << 9)];

	}

/* L1661: */
    }

    if (*lamber) {

/* ****     set the perturbed albedo equal to the first perturbed sur_b */
/*         value */

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb_b__[nz + (ibn * 5 + 5) * 10] = sur_b__[(ibn * 5 + 5 << 2) + 1]
		    ;
/* L1671: */
	}

    } else {

/* ****   use the disort routine dref to calculate the value of alb */
/*       that corresponds to this difference in the parameters */

	i__1 = *nref;
	for (nr = 1; nr <= i__1; ++nr) {
	    surf_pr__[nr - 1] = sur_b__[nr + (n * 5 + 5 << 2)];
/* L1701: */
	}
	if (*iref == 4) {

/* ****      set wind speed and direction for Cox/Munk model */

	    surf_pr__[2] = *ws;
	    surf_pr__[3] = *phiw;
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb_b__[nz + (ibn * 5 + 5) * 10] = dref_(&umu0[nz], surf_pr__, 
		    iref);
/* L1721: */
	}

    }

/* ****    compute denomenators for flux and radiance partial derivatives */
/*        optical depth, single scattering co-albedo, and surface albedo */

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	dref0 = alb_b__[nz + (ibn * 5 + 5) * 10] - alb_b__[nz + (ibn * 5 + 1) 
		* 10];
	if (dabs(dref0) > surfmin) {
	    dref_i__[nz - 1] = 1.f / dref0;
	} else {
	    dref_i__[nz - 1] = 0.f;
	}
/* L1741: */
    }

/* ****   define the amplitudes of the perturbations */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

/* ****       nominal calse */

	dx_b_i__[k + (ibn * 5 + 1) * 70] = 0.f;

/* ****       optical depth: n_rad(2,ibn) = 1 */

	deltau_b__ = dtau_b__[k + (ibn * 5 + 2) * 70] - dtau_b__[k + (ibn * 5 
		+ 1) * 70];
	if (deltau_b__ > *tauerr * .001f * *tauerr * dtau_b__[k + (ibn * 5 + 
		1) * 70] + *taumn) {
	    dx_b_i__[k + (ibn * 5 + 2) * 70] = 1.f / deltau_b__;
	} else {
	    dx_b_i__[k + (ibn * 5 + 2) * 70] = 0.f;
	}

/* ****       single scattering albedo: n_rad(3,ibn) = 1 */

	delpi0_b__ = copi0_b__[k + (ibn * 5 + 3) * 70] - copi0_b__[k + (ibn * 
		5 + 1) * 70];
	if (dabs(delpi0_b__) > ssamin) {
	    dx_b_i__[k + (ibn * 5 + 3) * 70] = 1.f / delpi0_b__;
	} else {
	    dx_b_i__[k + (ibn * 5 + 3) * 70] = 0.f;
	}

/* *****        phase function: n_rad(4,ibn) = 1 */

	delg_b__ = g_b__[k + (ibn * 5 + 4) * 70] - g_b__[k + (ibn * 5 + 1) * 
		70];
	if (dabs(delg_b__) > gmin) {
	    dx_b_i__[k + (ibn * 5 + 4) * 70] = 1.f / delg_b__;
	} else {
	    dx_b_i__[k + (ibn * 5 + 4) * 70] = 0.f;
	}

/* L2021: */
    }

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	dalb_b_i__[nz + ibn * 10] = dref_i__[nz - 1];
/* L2101: */
    }

    return 0;
} /* load_optics__ */

