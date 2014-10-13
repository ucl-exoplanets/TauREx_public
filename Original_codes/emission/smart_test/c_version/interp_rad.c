/* interp_rad.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int interp_rad__(logical *lsolar, logical *lplanck, integer *
	ng0, integer *l_1__, integer *l_2__, integer *nlev, integer *nz0, 
	integer *nphi, integer *numu, integer *nzup, integer *nzdn, 
	doublereal *wn_io__, real *umu0nz, real *umu, real *alb_b__, real *
	dtau_b__, real *copi0_b__, real *g_b__, doublereal *trnflx_b__, 
	doublereal *dtrnflxdx, doublereal *refflx_b__, doublereal *drefflxdx, 
	doublereal *absflx_b__, doublereal *dabsflxdx, doublereal *refrad_b__,
	 doublereal *drefraddx, doublereal *absrad_b__, doublereal *dabsraddx,
	 doublereal *brdf_b__, doublereal *dbrdfdx, real *tsurf, real *tatm, 
	real *alb0, real *dtau, real *copi0, real *g0, real *dalb, real *
	deltau, real *delpi0, real *delg, doublereal *dnsflxsrc_b__, 
	doublereal *ddnsflxdx, doublereal *upsflxsrc_b__, doublereal *
	dupsflxdx, doublereal *dntflxsrc_b__, doublereal *ddntflxdx, 
	doublereal *uptflxsrc_b__, doublereal *duptflxdx, doublereal *
	sradsrc_b__, doublereal *dsraddx, doublereal *tradsrc_b__, doublereal 
	*dtraddx, doublereal *bb, doublereal *bb_sur__, doublereal *
	bb_flx_up__, doublereal *bb_flx_dn__, doublereal *bb_rad__, 
	doublereal *dn_s_src__, doublereal *up_s_src__, doublereal *
	rad_s_src__, doublereal *dn_t_src__, doublereal *up_t_src__, 
	doublereal *rad_t_src__, doublereal *trndir, doublereal *trnflx, 
	doublereal *refflx, doublereal *absflx, doublereal *trnrad, 
	doublereal *refrad, doublereal *absrad, real *upsflx, real *dnsflx, 
	real *dirsflx, real *sol_rad__, real *uptflx, real *dntflx, real *
	th_rad__)
{
    /* Initialized data */

    static doublereal pi = 3.141592654;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int add_sflx__(integer *, doublereal *, 
	    doublereal *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    add_tflx__(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer k;
    static doublereal fb[70], ft[70], rl[70], tl[70];
    extern /* Subroutine */ int interp_src__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, real *, real *, real *, real *, real *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal dfb[70], fdn[70];
    static integer ibn;
    static doublereal drb[70], dir[70];
    static integer naz;
    static doublereal fup[70];
    static integer nze;
    static doublereal uft[70], urt[70], rad0[17920]	/* was [16][16][70] */
	    ;
    static integer lev1;
    static doublereal brdf[256]	/* was [16][16] */, flxd[70];
    static real wnio;
    static doublereal flxu[70];
    static integer nlyr;
    static doublereal dtau0;
    static integer l2tau;
    static real tsurf0[1];
    extern /* Subroutine */ int planck_(integer *, integer *, integer *, 
	    integer *, real *, real *, doublereal *);
    static real solflx;
    extern /* Subroutine */ int add_rad__(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal rad_src0__[17920]	/* was [16][16][70] */;


/* cccccccccccccccccccccc   i n t e r p _ r a d   cccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine linearly interpolates the binned radiance and   cc */
/* c    flux source functions for each layer to the optical depth,      cc */
/* c    single scatterng albedo, asymmetry parameter, surface pressure, cc */
/* c    and surface albedo at the wavelength of interest.               cc */
/* c                                                                    cc */
/* c    notes:                                                          cc */
/* c           This version uses a full surface bdrf, but only the      cc */
/* c           first element in the surface optical property vector     cc */
/* c           is used in the interpolation.  The assumption is that    cc */
/* c           variation in the other elements are correlated with      cc */
/* c           first element.                                           cc */
/* c           This version uses optical depth interpolation.           cc */
/* c                                                                    cc */
/* c    NOTE: this version has been modified to produce layer dependent cc */
/* c    flux source functions and their partial derivatives, as well    cc */
/* c    as the corresponding flux values.                               cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c         l1 - starting vertical index for source calculation        cc */
/* c         l2 - ending vertical index for source calculation          cc */
/* c        ng0 - spectral mapping bin index for this spectral point    cc */
/* c       nlev - number of model vertical levels                       cc */
/* c        nz0 - index of this solar zenith angle (1 to nsol)          cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c      wn_io - current wavenumber (cm**-1)                           cc */
/* c       nzup - number of upward steams for radiances                 cc */
/* c       nzdn - number of downward streams for radiances              cc */
/* c       numu - total number of streams for radiances                 cc */
/* c       nphi - total number of azimuths for radiances                cc */
/* c dnsflxsrc_b - downward flux source at each level for this bin      cc */
/* c upsflxsrc_b - downward flux source at each level for this bin      cc */
/* c   sradsrc_b - angle-dependent radiance source at each level for    cc */
/* c              this bin                                              cc */
/* c  ddnsflxdx - downward flux jacobian at each level for this bin     cc */
/* c  dupsflxdx - upward flux jacobian at each level for this bin       cc */
/* c    dsraddx - radiance jacobian at each level for this bin          cc */
/* c dntflxsrc_b - downward flux source at each level for this bin      cc */
/* c uptflxsrc_b - downward flux source at each level for this bin      cc */
/* c   tradsrc_b - angle-dependent radiance source at each level for    cc */
/* c              this bin                                              cc */
/* c  ddntflxdx - downward flux jacobian at each level for this bin     cc */
/* c  duptflxdx - upward flux jacobian at each level for this bin       cc */
/* c    dtraddx - radiance jacobian at each level for this bin          cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     uptflx - wn-dependent upward thermal flux                      cc */
/* c     dntflx - wn-dependent downward thermal flux                    cc */
/* c     th_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c     upsflx - wn-dependent upward solar flux                        cc */
/* c     dnsflx - wn-dependent downward direct+diffuse solar flux       cc */
/* c    dirsflx - wn-dependent downward direct solar flux               cc */
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c                                                                    cc */
/* cccccccccccccccccccccc   i n t e r p _ r a d   cccccccccccccccccccccccc */




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





/* ****   binned optical properties */


/* *****   state vector variables. */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   binned layer flux transmittances and absorptances */


/* ****   layer flux transmittances, absorptances, and fluxes */


/* ****   surface bi-directional reflection function */



/* ****   output upward and downward radiances */


/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */





/* ****   vertical structure variables at this wavelength */


/* ****    layer transmission values for simplified adding method */



/* ****   interpolated solar and thermal source functions */



/* ****    quantitities used on radiance adding */



/* ****   planck function thermal source variables */


/* ****   spectrally-dependent output flux and radiance */


/* ****   define pi */

    /* Parameter adjustments */
    th_rad__ -= 273;
    --dntflx;
    --uptflx;
    sol_rad__ -= 273;
    --dirsflx;
    --dnsflx;
    --upsflx;
    absrad -= 17;
    refrad -= 17;
    trnrad -= 17;
    --absflx;
    --refflx;
    --trnflx;
    --trndir;
    rad_t_src__ -= 273;
    --up_t_src__;
    --dn_t_src__;
    rad_s_src__ -= 273;
    --up_s_src__;
    --dn_s_src__;
    bb_rad__ -= 273;
    --bb_flx_dn__;
    --bb_flx_up__;
    --bb;
    dtraddx -= 89873;
    tradsrc_b__ -= 18193;
    dsraddx -= 89873;
    sradsrc_b__ -= 18193;
    duptflxdx -= 351;
    uptflxsrc_b__ -= 71;
    ddntflxdx -= 351;
    dntflxsrc_b__ -= 71;
    dupsflxdx -= 351;
    upsflxsrc_b__ -= 71;
    ddnsflxdx -= 351;
    dnsflxsrc_b__ -= 71;
    --delg;
    --delpi0;
    --deltau;
    --g0;
    --copi0;
    --dtau;
    --tatm;
    dbrdfdx -= 273;
    brdf_b__ -= 1553;
    dabsraddx -= 4497;
    absrad_b__ -= 6737;
    drefraddx -= 4497;
    refrad_b__ -= 6737;
    dabsflxdx -= 281;
    absflx_b__ -= 421;
    drefflxdx -= 281;
    refflx_b__ -= 421;
    dtrnflxdx -= 281;
    trnflx_b__ -= 421;
    g_b__ -= 421;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    alb_b__ -= 61;
    --umu;

    /* Function Body */

/* ****   specify the bin index */

    ibn = *ng0;
    nlyr = *nlev - 1;

/* ****  set the level range counters */

    if (*l_1__ == 1) {
	lev1 = 2;
    } else {
	lev1 = *l_1__;
    }
    if (*l_2__ < *nlev) {
	l2tau = *l_2__;
    } else {
	l2tau = *nlev - 1;
    }

/* ****    find the surface optical properties differences. */
/*        for surface optical properties, use only the first element */
/*        of the surface optical property vector for interpolation. */

/*      dalb = 0.0 */
    *dalb = *alb0 - alb_b__[*nz0 + (ibn * 5 + 1) * 10];

/* ****     find the optical depth and single scattering albedo */
/*         and asymmetry parameter differences. */

    i__1 = l2tau;
    for (k = *l_1__; k <= i__1; ++k) {
	deltau[k] = dtau[k] - dtau_b__[k + (ibn * 5 + 1) * 70];
	delpi0[k] = copi0[k] - copi0_b__[k + (ibn * 5 + 1) * 70];
	delg[k] = g0[k] - g_b__[k + (ibn * 5 + 1) * 70];
/*          deltau(k) = 0.0 */
/*          delpi0(k) = 0.0 */
/*          delg(k) = 0.0 */
/* L1001: */
    }

/* ****     Layer transmittances are used in mapping fluxes and */
/*         radiances back to the high resolution spectrum. */
/*         find the transmission through each layer for each stream. */

    i__1 = l2tau;
    for (k = *l_1__; k <= i__1; ++k) {

/* ****            find the flux transmittances and absorptances */

	trnflx[k] = trnflx_b__[k + (ibn * 5 + 1) * 70] + dtrnflxdx[k + (ibn * 
		3 + 1) * 70] * deltau[k] + dtrnflxdx[k + (ibn * 3 + 2) * 70] *
		 delpi0[k] + dtrnflxdx[k + (ibn * 3 + 3) * 70] * delg[k];
	refflx[k] = refflx_b__[k + (ibn * 5 + 1) * 70] + drefflxdx[k + (ibn * 
		3 + 1) * 70] * deltau[k] + drefflxdx[k + (ibn * 3 + 2) * 70] *
		 delpi0[k] + drefflxdx[k + (ibn * 3 + 3) * 70] * delg[k];

/* ****            find the radiance transmittances and absorptances */

	i__2 = *numu;
	for (nze = 1; nze <= i__2; ++nze) {
	    dtau0 = dtau[k] / (r__1 = umu[nze], dabs(r__1));
	    trnrad[nze + (k << 4)] = exp(-dtau0);
	    refrad[nze + (k << 4)] = refrad_b__[nze + (k + (ibn * 5 + 1) * 70 
		    << 4)] + drefraddx[nze + (k + (ibn * 3 + 1) * 70 << 4)] * 
		    deltau[k] + drefraddx[nze + (k + (ibn * 3 + 2) * 70 << 4)]
		     * delpi0[k] + drefraddx[nze + (k + (ibn * 3 + 3) * 70 << 
		    4)] * delg[k];
/* L1201: */
	}
/* L1221: */
    }

    refflx[*nlev] = *alb0;
    trnflx[*nlev] = 0.f;
    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {
	i__2 = *numu;
	for (nze = 1; nze <= i__2; ++nze) {
	    brdf[nze + (naz << 4) - 17] = brdf_b__[nze + (naz + (ibn * 5 + 1 
		    << 4) << 4)] + dbrdfdx[nze + (naz + (ibn << 4) << 4)] * *
		    dalb;
/* L1241: */
	}
/* L1261: */
    }

    if (*lplanck) {
	i__1 = l2tau;
	for (k = *l_1__; k <= i__1; ++k) {

/* ****            find the flux transmittances and absorptances */

	    absflx[k] = absflx_b__[k + (ibn * 5 + 1) * 70] + dabsflxdx[k + (
		    ibn * 3 + 1) * 70] * deltau[k] + dabsflxdx[k + (ibn * 3 + 
		    2) * 70] * delpi0[k] + dabsflxdx[k + (ibn * 3 + 3) * 70] *
		     delg[k];
	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {

/* ****            find the radiance transmittances and absorptances */

		absrad[nze + (k << 4)] = absrad_b__[nze + (k + (ibn * 5 + 1) *
			 70 << 4)] + dabsraddx[nze + (k + (ibn * 3 + 1) * 70 
			<< 4)] * deltau[k] + dabsraddx[nze + (k + (ibn * 3 + 
			2) * 70 << 4)] * delpi0[k] + dabsraddx[nze + (k + (
			ibn * 3 + 3) * 70 << 4)] * delg[k];
/* L1401: */
	    }
/* L1421: */
	}

    }

    if (*lsolar) {

/* ****    find the direct fluxes, dir */

	solflx = 1.f;
	dir[0] = *umu0nz * solflx;
	i__1 = nlyr;
	for (k = 1; k <= i__1; ++k) {
	    dir[k] = dir[k - 1] * trndir[k];
/* L1601: */
	}

    }

/* ****    i n t e r p o l a t e   s o u r c e   f u n c t i o n s */

    interp_src__(lsolar, lplanck, ng0, nlev, l_1__, l_2__, &lev1, &l2tau, 
	    nzup, nzdn, numu, nphi, &dnsflxsrc_b__[71], &upsflxsrc_b__[71], &
	    sradsrc_b__[18193], &ddnsflxdx[351], &dupsflxdx[351], &dsraddx[
	    89873], &dntflxsrc_b__[71], &uptflxsrc_b__[71], &tradsrc_b__[
	    18193], &ddntflxdx[351], &duptflxdx[351], &dtraddx[89873], dalb, &
	    deltau[1], &delpi0[1], &delg[1], alb0, wn_io__, &up_s_src__[1], &
	    dn_s_src__[1], &rad_s_src__[273], &up_t_src__[1], &dn_t_src__[1], 
	    &rad_t_src__[273]);

    if (*lsolar) {

/* ********************************************************************** */

/* *****         f i n d   s o l a r    r a d i a n c e s */

/* ********************************************************************** */

/* ****      initialize solar radiances and fluxes at all levels */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    upsflx[k] = 0.f;
	    dnsflx[k] = 0.f;
	    dirsflx[k] = 0.f;

/* ****        initialize output radiances */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = 1; nze <= i__3; ++nze) {
		    sol_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L2001: */
		}
/* L2021: */
	    }
/* L2041: */
	}

/* ****     use the flux adding method to find the */
/*         upward and downward fluxes */

	fdn[0] = dn_s_src__[1];
	i__1 = nlyr;
	for (k = 1; k <= i__1; ++k) {
	    rl[k - 1] = refflx[k];
	    tl[k - 1] = trnflx[k];
	    fdn[k] = dn_s_src__[k + 1];
	    fup[k - 1] = up_s_src__[k];
/* L2201: */
	}

	rl[*nlev - 1] = *alb0;
	tl[*nlev - 1] = 0.;

/* *****    define the surface flux source function - for fluxes */
/*         the reflected solar gives and exact solution (everything */
/*         else is accouted for in the other terms). */

	fup[*nlev - 1] = *alb0 * dir[*nlev - 1];

	add_sflx__(&nlyr, tl, rl, umu0nz, dir, fup, fdn, urt, drb, uft, dfb, 
		ft, fb, flxu, flxd);

/* ****         set direct and diffuse fluxes at all levels */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    upsflx[k] = (real) flxu[k - 1];
	    dnsflx[k] = (real) flxd[k - 1];
	    dirsflx[k] = (real) dir[k - 1];
/* L2441: */
	}

/* ****     find upward and downward radiances at all levels */

/* ****     load the radiance source terms */

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {
	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		i__3 = *nlev;
		for (k = 1; k <= i__3; ++k) {
		    rad_src0__[nze + (naz + (k << 4) << 4) - 273] = 
			    rad_s_src__[nze + (naz + (k << 4) << 4)];
/* L2601: */
		}
/* L2621: */
	    }
/* L2641: */
	}

	add_rad__(&nlyr, nzdn, nzup, numu, nphi, &trnrad[17], &refrad[17], 
		brdf, rad_src0__, ft, fb, urt, drb, uft, dfb, rad0);

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {
	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		i__3 = *nlev;
		for (k = 1; k <= i__3; ++k) {
		    sol_rad__[nze + (naz + (k << 4) << 4)] = (real) rad0[nze 
			    + (naz + (k << 4) << 4) - 273];
/* L2801: */
		}
/* L2821: */
	    }
/* L2841: */
	}

/* ****    exit solar block */

    }

/* *********************************************************************** */

/* ****           f i n d    t h e r m a l    r a d i a n c e s */

/* *********************************************************************** */

    if (*lplanck && *nz0 == 1) {

/* ****     find the surface planck source at this wavenumber */

	wnio = (real) (*wn_io__);
	if (*l_2__ == *nlev) {
	    tsurf0[0] = *tsurf;
	    planck_(&c__1, &c__1, &c__1, &c__1, &wnio, tsurf0, &bb[1]);
	    *bb_sur__ = bb[1];
	    bb_flx_up__[*nlev] = pi * *bb_sur__ * (1. - *alb0);
	}

/* ****    find the atmospheric planck functions at this wavenumber */

	planck_(&c__1, l_2__, &c__1, &c__1, &wnio, &tatm[1], &bb[1]);

/* ****         define the black body source functions for flux and */
/*             radiance, assuming a linear-in-tau formulation. */

	if (*l_1__ == 1) {
	    bb_flx_dn__[1] = pi * bb[1] * copi0[1] * absflx[1];
	}
	i__1 = l2tau;
	for (k = *l_1__; k <= i__1; ++k) {
	    bb_flx_dn__[k + 1] = pi * .5 * (bb[k] + bb[k + 1]) * copi0[k] * 
		    absflx[k];
	    bb_flx_up__[k] = pi * .5 * (bb[k] + bb[k + 1]) * copi0[k] * 
		    absflx[k];
/* L3001: */
	}

/* ****    thermal sources for downward radiance streams at top of */
/*        atmosphere and upward radiance streams at surface */

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {
	    if (*l_1__ == 1) {

/* ****          the following creates a non-zero scaling factor for */
/*              the top layer */

		i__2 = *nzdn;
		for (nze = 1; nze <= i__2; ++nze) {
		    bb_rad__[nze + (naz + 16 << 4)] = bb[1] * copi0[1] * 
			    absrad[nze + 16];
/* L3201: */
		}

	    }

	    if (*l_2__ == *nlev) {

/* ****           set upward thermal source at surface */

		i__2 = *numu;
		for (nze = *nzup; nze <= i__2; ++nze) {
		    bb_rad__[nze + (naz + (*nlev << 4) << 4)] = *bb_sur__ * (
			    1. - brdf[nze + (naz << 4) - 17]);
/* L3221: */
		}

	    }
/* L3241: */
	}

/* ****             thermal source for downward streams */

	i__1 = l2tau;
	for (k = *l_1__; k <= i__1; ++k) {
	    naz = 1;
	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {
		bb_rad__[nze + (naz + (k + 1 << 4) << 4)] = (bb[k] + bb[k + 1]
			) * .5 * copi0[k] * absrad[nze + (k << 4)];
/* L3401: */
	    }

/* ****         set values at all other azimuth angles */
/*             (no azimuth dependence above surface) */

	    i__2 = *nphi;
	    for (naz = 2; naz <= i__2; ++naz) {
		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
		    bb_rad__[nze + (naz + (k + 1 << 4) << 4)] = bb_rad__[nze 
			    + ((k + 1 << 4) + 1 << 4)];
/* L3421: */
		}
/* L3441: */
	    }
/* L3461: */
	}

/* ****       thermal source for upward streams */

	i__1 = l2tau;
	for (k = *l_1__; k <= i__1; ++k) {
	    naz = 1;
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		bb_rad__[nze + (naz + (k << 4) << 4)] = (bb[k] + bb[k + 1]) * 
			.5 * copi0[k] * absrad[nze + (k << 4)];
/* L3601: */
	    }
	    i__2 = *nphi;
	    for (naz = 2; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
		    bb_rad__[nze + (naz + (k << 4) << 4)] = bb_rad__[nze + ((
			    k << 4) + 1 << 4)];
/* L3621: */
		}
/* L3641: */
	    }
/* L3661: */
	}

/* ****     initialize thermal fluxes and radiances at all levels */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    uptflx[k] = 0.f;
	    dntflx[k] = 0.f;

/* ****     initialize thermal radiances at all levels */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = 1; nze <= i__3; ++nze) {
		    th_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L4201: */
		}
/* L4221: */
	    }
/* L4241: */
	}

/* ****  use the flux adding method to find the upward and downward fluxes */

	i__1 = nlyr;
	for (k = 1; k <= i__1; ++k) {
	    rl[k - 1] = refflx[k];
	    tl[k - 1] = trnflx[k];
	    fdn[k - 1] = bb_flx_dn__[k] * dn_t_src__[k] + bb_flx_dn__[k];
	    fup[k - 1] = bb_flx_up__[k] * up_t_src__[k] + bb_flx_up__[k];
/* L4401: */
	}
	k = *nlev;
	rl[k - 1] = *alb0;
	tl[k - 1] = 0.;
	fdn[k - 1] = bb_flx_dn__[k] * dn_t_src__[k] + bb_flx_dn__[k];
	fup[k - 1] = bb_flx_up__[k] * up_t_src__[k] + bb_flx_up__[k];

	add_tflx__(&nlyr, tl, rl, fup, fdn, urt, drb, uft, dfb, ft, fb, flxu, 
		flxd);

/* ****         set downward direct and diffuse fluxes at all levels */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    uptflx[k] = (real) flxu[k - 1];
	    dntflx[k] = (real) flxd[k - 1];
/* L4421: */
	}
/*      if(nz0 .eq. 1 .and.  (l_1 .eq. 1 .or. l_1 .eq. 10)) then */
/*       write(*,'(1a,16(1pe12.4))') 'dn_t_src:      ', */
/*     - (dn_t_src(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'up_t_src:      ', */
/*     - (up_t_src(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'dntflx:      ', */
/*     - (dntflx(k),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'uptflx:      ', */
/*     - (uptflx(k),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'bb_flx_dn    ', */
/*     - (bb_flx_dn(k),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'bb_flx_up    ', */
/*     - (bb_flx_up(k),k=1,nlev) */
/*      do 7 nze=1,numu */
/* 7      write(*,'(1a,i3,16(1pe12.4))') 'bb_rad ',nze, */
/*     - (bb_rad(nze,1,k),k=1,nlev) */
/*      endif */

/* ****     rescale the radiance source terms */

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {

/* ****        scale downward thermal sources at top of atmosphere */

	    k = 1;
	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {
/*                rad_src0(nze,naz,k) = bb_rad(nze,naz,k)* */
/*     -                                  rad_t_src(nze,naz,k) + */
/*     -                                  bb_rad(nze,naz,k) */
		rad_src0__[nze + (naz + (k << 4) << 4) - 273] = rad_t_src__[
			nze + (naz + (k << 4) << 4)] + bb_rad__[nze + (naz + (
			k << 4) << 4)];
/* L4601: */
	    }

/* ****        scale upward thermal sources at surface */

	    k = *nlev;
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
/*                rad_src0(nze,naz,k) = bb_rad(nze,naz,k)* */
/*     -                                rad_t_src(nze,naz,k) + */
/*     -                                bb_rad(nze,naz,k) */
		rad_src0__[nze + (naz + (k << 4) << 4) - 273] = rad_t_src__[
			nze + (naz + (k << 4) << 4)] + bb_rad__[nze + (naz + (
			k << 4) << 4)];
/* L4621: */
	    }
/* L4641: */
	}

/* *****    scale upward and downward thermal sources at other levels */

	i__1 = nlyr;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
/*                        rad_src0(nze,naz,k+1) = bb_rad(nze,naz,k)* */
/*     -                                         rad_t_src(nze,naz,k+1) + */
/*     -                                          bb_rad(nze,naz,k) */
		    rad_src0__[nze + (naz + (k + 1 << 4) << 4) - 273] = 
			    rad_t_src__[nze + (naz + (k + 1 << 4) << 4)] + 
			    bb_rad__[nze + (naz + (k << 4) << 4)];
/* L4801: */
		}
		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
/*                        rad_src0(nze,naz,k) = bb_rad(nze,naz,k)* */
/*     -                                        rad_t_src(nze,naz,k) + */
/*     -                                        bb_rad(nze,naz,k) */
		    rad_src0__[nze + (naz + (k << 4) << 4) - 273] = 
			    rad_t_src__[nze + (naz + (k << 4) << 4)] + 
			    bb_rad__[nze + (naz + (k << 4) << 4)];
/* L4821: */
		}
/* L4841: */
	    }
/* L4861: */
	}

	add_rad__(&nlyr, nzdn, nzup, numu, nphi, &trnrad[17], &refrad[17], 
		brdf, rad_src0__, ft, fb, urt, drb, uft, dfb, rad0);

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {
	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		i__3 = *nlev;
		for (k = 1; k <= i__3; ++k) {
		    th_rad__[nze + (naz + (k << 4) << 4)] = (real) rad0[nze + 
			    (naz + (k << 4) << 4) - 273];
/* L4901: */
		}
/* L4921: */
	    }
/* L4941: */
	}

    }

    return 0;
} /* interp_rad__ */

