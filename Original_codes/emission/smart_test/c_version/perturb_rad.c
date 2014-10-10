/* perturb_rad.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int perturb_rad__(logical *lsolar, logical *lplanck, logical 
	*lamber, logical *usrang, integer *ni0, integer *nz0, integer *nlyr, 
	integer *iref, integer *nlout, integer *levout, integer *k_out__, 
	integer *ibin, integer *nphi, integer *nstr, integer *numu, integer *
	nzup, integer *nzdn, integer *npd, integer *npd1, integer *ipd, 
	integer *igs, integer *modepd, doublereal *wn_io__, doublereal *dvi, 
	real *umu_f__, real *gwt_f__, real *umu, real *phi, real *umu0nz, 
	real *phi0nz, real *pd_pert__, real *dx, real *alb, real *ts, real *t,
	 real *p, real *dp_dp__, real *rmix, real *dtauaer, real *tau_ext__, 
	real *tau_sca__, real *g_sca__, real *dtauex, real *co_pi0__, real *g,
	 real *alb_b__, real *dtau_b__, real *copi0_b__, real *g_b__, 
	doublereal *trnflx_b__, doublereal *dtrnflxdx, doublereal *refflx_b__,
	 doublereal *drefflxdx, doublereal *absflx_b__, doublereal *dabsflxdx,
	 doublereal *refrad_b__, doublereal *drefraddx, doublereal *
	absrad_b__, doublereal *dabsraddx, doublereal *brdf_b__, doublereal *
	dbrdfdx, doublereal *dnsflxsrc_b__, doublereal *upsflxsrc_b__, 
	doublereal *sradsrc_b__, doublereal *dntflxsrc_b__, doublereal *
	uptflxsrc_b__, doublereal *tradsrc_b__, doublereal *ddntflxdx, 
	doublereal *duptflxdx, doublereal *dtraddx, doublereal *ddnsflxdx, 
	doublereal *dupsflxdx, doublereal *dsraddx, real *dalb0, real *
	deltau0, real *delpi00, real *delg0, doublereal *bb0, doublereal *
	bb_sur0__, doublereal *bb_flx_up0__, doublereal *bb_flx_dn0__, 
	doublereal *bb_rad0__, doublereal *dn_s_src0__, doublereal *
	up_s_src0__, doublereal *rad_s_src0__, doublereal *dn_t_src0__, 
	doublereal *up_t_src0__, doublereal *rad_t_src0__, doublereal *
	trndir0, doublereal *trnflx0, doublereal *refflx0, doublereal *
	absflx0, doublereal *trnrad0, doublereal *refrad0, doublereal *
	absrad0, doublereal *trndir1_i__, doublereal *trnflx1_i__, doublereal 
	*refflx1_i__, doublereal *absflx1_i__, doublereal *dtrndir1_i__, 
	doublereal *dtrnflx1_i__, doublereal *drefflx1_i__, doublereal *
	dabsflx1_i__, doublereal *dn_s_src1__, doublereal *up_s_src1__, 
	doublereal *dup_s_src1__, doublereal *ddn_s_src1__, doublereal *
	dn_t_src1__, doublereal *up_t_src1__, doublereal *dup_t_src1__, 
	doublereal *ddn_t_src1__, doublereal *sol_rad1__, doublereal *
	dsol_rad1__, doublereal *th_rad1__, doublereal *dth_rad1__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal dn_s_src__[70], dn_t_src__[70], up_s_src__[70], 
	    up_t_src__[70];
    static integer k, n;
    static doublereal bb_flx_dn__[70], rad_s_src__[17920]	/* was [16][
	    16][70] */, bb_flx_up__[70], rad_t_src__[17920]	/* was [16][
	    16][70] */;
    static real g0[70];
    static integer l1, l2;
    static doublereal bb[70];
    static integer ni, nz, l_1__, l_2__, ng0;
    extern /* Subroutine */ int interp_rad__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, real *, real *, real *, real *
	    , real *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static integer lpd;
    static real dns[3];
    static integer naz, nze;
    static real ups[3], alb0;
    extern /* Subroutine */ int int_rad_lev__(logical *, logical *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *), 
	    grey_eq_trn__(logical *, logical *, logical *, logical *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static real dalb, delg[70], dtau[70], dnth[3], dirs[3];
    static integer nlev;
    static real tatm[70], upth[3], copi0[70], fbeam, rad_s__[768]	/* 
	    was [16][16][3] */, btemp, temis, ttemp, tsurf, delpi0[70];
    static doublereal bb_rad__[17920]	/* was [16][16][70] */, absrad[1120]	
	    /* was [16][70] */;
    static real rad_th__[768]	/* was [16][16][3] */;
    static doublereal refrad[1120]	/* was [16][70] */;
    static real th_rad__[17920]	/* was [16][16][70] */;
    static doublereal bb_sur__;
    static real deltau[70];
    static doublereal absflx[70], refflx[70], trnrad[1120]	/* was [16][
	    70] */;
    static real dnsflx[70], dntflx[70];
    static doublereal trndir[70], trnflx[70];
    static real upsflx[70], uptflx[70], dtau_l1__, dtau_l2__;
    static doublereal abs_cmu__[1120]	/* was [16][70] */;
    static real sol_rad__[17920]	/* was [16][16][70] */, dtausca;
    extern /* Subroutine */ int rad_trn__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static real dtaumin;
    static doublereal trn_cmu__[1120]	/* was [16][70] */;
    static real dirsflx[70], tautest, dtausca0;


/* ccccccccccccccccccccccc  p e r t u r b _ r a d   cccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine finds evaluates the partial derivatives of the  cc */
/* c    radiances with respect to each variable member of the state     cc */
/* c    vector.                                                         cc */
/* c                                                                    cc */
/* c    NOTE: this version has been modified to produce layer dependent cc */
/* c    flux source functions and their partial derivatives, rather     cc */
/* c    than the corresponding flux values.                             cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        nz0 - index of this solar zenith angle (1 to nsol)          cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       numu - number of zenith angles in input file                 cc */
/* c     uptflx - wn-dependent upward thermal flux                      cc */
/* c     dntflx - wn-dependent downward thermal flux                    cc */
/* c     th_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c     upsflx - wn-dependent upward solar flux                        cc */
/* c     dnsflx - wn-dependent downward diffuse + direct solar flux     cc */
/* c    dirsflx - wn-dependent downward direct solar flux               cc */
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c     th_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccc  p e r t u r b _ r a d   cccccccccccccccccccccc */




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




/* ****   counters used in map_back and write_mono */


/* ****   spectral binning parameters */


/* ****   gaussian angles for flux integration */


/* ****   spectral binning parameters */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   binned layer flux transmittances and absorptances */


/* ****    layer radiance transmission values for simplified adding method */


/* ****    layer flux transmission values for simplified adding method */



/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */





/* ****   vertical structure variables at this wavelength */




/* *****   variables passed to int_rad_lev */


/* ****   spectrally-dependent output flux and radiance */


/* ****    optical property differences for this wn */


/* ****    planck functions used in interp_rad */


/* ****   interpolated solar and thermal source functions */


/* ****    baseline optical property differences for this wn */


/* ****    baseline planck functions used in interp_rad */


/* ****    baseline layer transmission values for simplified adding method */


/* ****   baseline flux and radiance source functions */


/* ****     perturrbed monochormatic radiance and wn interpolation values */


/* ****   flux source functions at wavenumber wn_io for perturbed state */


/* ****   flux source functions wavenumber derivatives */


/* ****    flux transmission and absorption needed for simplified */
/*        adding method at wavenumber wn_io */



/* ****    specify the solar zenith angle and gas index */

    /* Parameter adjustments */
    dth_rad1__ -= 54801;
    th_rad1__ -= 592401;
    dsol_rad1__ -= 592401;
    sol_rad1__ -= 1667601;
    ddn_t_src1__ -= 71;
    dup_t_src1__ -= 71;
    up_t_src1__ -= 771;
    dn_t_src1__ -= 771;
    ddn_s_src1__ -= 771;
    dup_s_src1__ -= 771;
    up_s_src1__ -= 2171;
    dn_s_src1__ -= 2171;
    dabsflx1_i__ -= 71;
    drefflx1_i__ -= 71;
    dtrnflx1_i__ -= 71;
    dtrndir1_i__ -= 771;
    absflx1_i__ -= 771;
    refflx1_i__ -= 771;
    trnflx1_i__ -= 771;
    trndir1_i__ -= 2171;
    absrad0 -= 17;
    refrad0 -= 17;
    trnrad0 -= 17;
    --absflx0;
    --refflx0;
    --trnflx0;
    --trndir0;
    rad_t_src0__ -= 273;
    --up_t_src0__;
    --dn_t_src0__;
    rad_s_src0__ -= 273;
    --up_s_src0__;
    --dn_s_src0__;
    bb_rad0__ -= 273;
    --bb_flx_dn0__;
    --bb_flx_up0__;
    --bb0;
    --delg0;
    --delpi00;
    --deltau0;
    dsraddx -= 89873;
    dupsflxdx -= 351;
    ddnsflxdx -= 351;
    dtraddx -= 89873;
    duptflxdx -= 351;
    ddntflxdx -= 351;
    tradsrc_b__ -= 18193;
    uptflxsrc_b__ -= 71;
    dntflxsrc_b__ -= 71;
    sradsrc_b__ -= 18193;
    upsflxsrc_b__ -= 71;
    dnsflxsrc_b__ -= 71;
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
    --g;
    co_pi0__ -= 71;
    dtauex -= 71;
    g_sca__ -= 71;
    tau_sca__ -= 71;
    tau_ext__ -= 281;
    dtauaer -= 11;
    rmix -= 71;
    --dp_dp__;
    --p;
    --t;
    --alb;
    dx -= 71;
    --pd_pert__;
    --phi;
    --umu;
    --gwt_f__;
    --umu_f__;
    --modepd;
    --igs;
    --ipd;
    --k_out__;

    /* Function Body */
    nz = *nz0;
    ni = *ni0;
    nlev = *nlyr + 1;

/* ****       define a small value of delta tau.  If the optical depth */
/*           perturbation is smaller than this, but greater than taumn, */
/*           set the value to this value */

    dtaumin = 5e-4f;

/* ****       evaluate the radiances and partial derivatives */
/*           at this wavenumber for each optical property */
/*           that can vary */

    i__1 = *npd;
    for (n = 1; n <= i__1; ++n) {

/* ****        determine the number of levels at which the state */
/*            vector variable chages (T=nlev, tau=nlyr, alb=1) */

	if ((i__2 = ipd[n], abs(i__2)) == 1) {

/* ****         surface pressure */

	    l1 = nlev;
	    l2 = nlev;

	} else {

	    if ((i__2 = ipd[n], abs(i__2)) == 2) {

/* ****            atmospheric/surface temperature */

		l1 = 1;
		l2 = nlev;
	    } else {

		if ((i__2 = ipd[n], abs(i__2)) == 3) {

/* ****             gas optical depth */

		    l1 = 1;
		    l2 = nlev;

		} else {

		    if ((i__2 = ipd[n], abs(i__2)) == 4) {

/* ****               aerosol optical depth */

			l1 = 1;
			l2 = nlev;

		    } else {

/* ****               surface albedo */

			l1 = nlev;
			l2 = nlev;

		    }
		}
	    }
	}

/* ****      enter the loop over perturbation level */

	i__2 = l2;
	for (lpd = l1; lpd <= i__2; ++lpd) {

/* ****           reload baseline arrays for optical properties, */
/*               planck functions, and source functions */

	    alb0 = alb[nz];
	    dalb = *dalb0;
	    i__3 = nlev - 1;
	    for (k = 1; k <= i__3; ++k) {
		dtau[k - 1] = dtauex[k + 70];
		copi0[k - 1] = co_pi0__[k + 70];
		g0[k - 1] = g[k];
		deltau[k - 1] = deltau0[k];
		delpi0[k - 1] = delpi00[k];
		delg[k - 1] = delg0[k];
/* L1001: */
	    }

	    i__3 = nlev - 1;
	    for (k = 1; k <= i__3; ++k) {
		trndir[k - 1] = trndir0[k];
		trnflx[k - 1] = trnflx0[k];
		refflx[k - 1] = refflx0[k];
		absflx[k - 1] = absflx0[k];
		i__4 = *numu;
		for (nze = 1; nze <= i__4; ++nze) {
		    trnrad[nze + (k << 4) - 17] = trnrad0[nze + (k << 4)];
		    refrad[nze + (k << 4) - 17] = refrad0[nze + (k << 4)];
		    absrad[nze + (k << 4) - 17] = absrad0[nze + (k << 4)];
/* L1201: */
		}
/* L1261: */
	    }

	    if (*lsolar) {
		i__3 = nlev;
		for (k = 1; k <= i__3; ++k) {
		    dn_s_src__[k - 1] = dn_s_src0__[k];
		    up_s_src__[k - 1] = up_s_src0__[k];
		    i__4 = *nphi;
		    for (naz = 1; naz <= i__4; ++naz) {
			i__5 = *numu;
			for (nze = 1; nze <= i__5; ++nze) {
			    rad_s_src__[nze + (naz + (k << 4) << 4) - 273] = 
				    rad_s_src0__[nze + (naz + (k << 4) << 4)];
/* L1441: */
			}
/* L1461: */
		    }
/* L1481: */
		}
	    }

	    if (*lplanck) {
		tsurf = *ts;
		bb_sur__ = *bb_sur0__;
		i__3 = nlev;
		for (k = 1; k <= i__3; ++k) {
		    tatm[k - 1] = t[k];
		    bb[k - 1] = bb0[k];
		    bb_flx_up__[k - 1] = bb_flx_up0__[k];
		    bb_flx_dn__[k - 1] = bb_flx_dn0__[k];
		    dn_t_src__[k - 1] = dn_t_src0__[k];
		    up_t_src__[k - 1] = up_t_src0__[k];
		    i__4 = *nphi;
		    for (naz = 1; naz <= i__4; ++naz) {
			i__5 = *numu;
			for (nze = 1; nze <= i__5; ++nze) {
			    bb_rad__[nze + (naz + (k << 4) << 4) - 273] = 
				    bb_rad0__[nze + (naz + (k << 4) << 4)];
			    rad_t_src__[nze + (naz + (k << 4) << 4) - 273] = 
				    rad_t_src0__[nze + (naz + (k << 4) << 4)];
/* L1641: */
			}
/* L1661: */
		    }
/* L1681: */
		}
	    }

/* ****           define the perturbation increments, dx */

	    if ((i__3 = ipd[n], abs(i__3)) == 1) {

/*                define the surface pressure change */

		dx[lpd + n * 70] = pd_pert__[n] * p[nlev];

/* *****           calculate radiances for perturbed surface pressure: */
/*                substitute dtauex(nlev,1) for dtauex(nlyr,1) */

/* ****             replace the differential optical depth and */
/*                 single scattering co-albedo in the lowest layer */

		dtau[*nlyr - 1] = dtauex[nlev + 70];
		copi0[*nlyr - 1] = co_pi0__[nlev + 70];

/* ****             set the levels for recalculating source functions */

		l_1__ = nlev - 1;
		l_2__ = nlev;

	    } else {

		if ((i__3 = ipd[n], abs(i__3)) == 2) {

/* ****               define atmospheric or surface temperature change */

		    if (lpd <= nlev) {
			dx[lpd + n * 70] = pd_pert__[n] * t[lpd];
		    } else {
			dx[lpd + n * 70] = pd_pert__[n] * *ts;
		    }

/* *****             calculate radiances for perturbed temperatures */

		    if (lpd <= nlev) {
			tatm[lpd - 1] = (pd_pert__[n] + 1.f) * t[lpd];
		    } else {
			tsurf = (pd_pert__[n] + 1.f) * *ts;
		    }

/* ****             set the levels for recalculating source functions */

		    if (lpd == 1) {
			l_1__ = 1;
		    } else {
			l_1__ = lpd - 1;
		    }
		    l_2__ = lpd;

/* ****               modify optical depths and single scattering abledos */
/*                   for temperature dependence */

		    dtau[l_1__ - 1] = dtauex[l_1__ + 140];
		    copi0[l_1__ - 1] = co_pi0__[l_1__ + 140];
		    if (l_2__ < nlev) {
			dtau[l_2__ - 1] = dtauex[l_2__ + 140];
			copi0[l_2__ - 1] = co_pi0__[l_2__ + 140];
		    }

		} else {

		    if ((i__3 = ipd[n], abs(i__3)) == 3) {

/* ****               gas mixing ratio perturbation */

/* ****               set the levels for recalculating source functions */

			if (lpd == 1) {
			    l_1__ = 1;
			} else {
			    l_1__ = lpd - 1;
			}

			l_2__ = lpd;

/* ****                 define the gas optical depth change */

			dx[lpd + n * 70] = pd_pert__[n] * rmix[lpd + igs[n - *
				npd1] * 70];

/* *****               find optical depths for perturbed gas mixing ratios */

/* ****                     note: there should be no perturbation in */
/*                               g(k) since the asymmetry parameter */
/*                               for a gas, g_sca, should be 0.0 */

/* ****                    a change in mixing ratio at level lpd affects */
/*                        the optical depth in layer l_1 and l_2 */

/* ****                    modify gas rmix at bottom of layer l_2. */
/*                        subtract off the nominal optical depth of */
/*                        the gas, tau_ext(l_1,1,n-npd1), and add back */
/*                        the perturbed value tau_ext(l_1,3,n-npd1) */

			dtau_l1__ = -tau_ext__[l_1__ + ((n - *npd1) * 3 + 1) *
				 70] + tau_ext__[l_1__ + ((n - *npd1) * 3 + 3)
				 * 70];

/* ****                    modify gas mixing ratio at top of layer l_2 */
/*                        subtract off the nominal optical depth of */
/*                        the gas, tau_ext(l_3,1,n-npd1), and add back */
/*                        the perturbed value tau_ext(l_3,3,n-npd1) */

			if (l_2__ < nlev) {
			    dtau_l2__ = -tau_ext__[l_2__ + ((n - *npd1) * 3 + 
				    1) * 70] + tau_ext__[l_2__ + ((n - *npd1) 
				    * 3 + 2) * 70];
			} else {
			    dtau_l2__ = 0.f;
			}

/* ****                  define a delta-tau test variable that will be */
/*                      used to define dx in the case where dtau_l1 and */
/*                      dtau_l2 are very small */

			tautest = (dtau_l1__ + dtau_l2__) * .5f;

			if (tautest >= dtaumin) {

/* ****                   adjust optical depths in layers l_1 and l_2 */
/*                       by an amount consistent with */
/*                       pd_pert(n)*rmix(lpd) */

			    dtau[l_1__ - 1] = dtauex[l_1__ + 70] + dtau_l1__;
			    if (l_2__ < nlev) {
				dtau[l_2__ - 1] = dtauex[l_2__ + 70] + 
					dtau_l2__;
			    }

			} else {

/* ****                     don't adjust the optical depth */

			    dx[lpd + n * 70] = 0.f;
			    dtau[l_1__ - 1] = dtauex[l_1__ + 70];
			    if (l2 < nlev) {
				dtau[l_2__ - 1] = dtauex[l_2__ + 70];
			    }

			}

			if (dtau[l_1__ - 1] != 0.f) {

/* ****                  find the total scattering optical depth of the */
/*                      unperturbed layer, l_1 */

			    dtausca0 = (1.f - co_pi0__[l_1__ + 70]) * dtauex[
				    l_1__ + 70];

/* ****                  then next statement assumes that perturbing the */
/*                      mixing ratio does not change the scattering */
/*                      optical depth, while it does change the */
/*                      single scattering albedo. */

			    copi0[l_1__ - 1] = 1.f - dtausca0 / dtau[l_1__ - 
				    1];

			} else {

			    copi0[l_1__ - 1] = 0.f;

			}

			if (dtau[l_2__ - 1] != 0.f && l_2__ < nlev) {

/* ****                  find the total scattering optical depth of the */
/*                      unperturbed layer l_2 */

			    dtausca0 = (1.f - co_pi0__[l_2__ + 70]) * dtauex[
				    l_2__ + 70];

			    copi0[l_2__ - 1] = 1.f - dtausca0 / dtau[l_2__ - 
				    1];

			} else {

			    copi0[l_2__ - 1] = 0.f;

			}
/*      if(wn_io .gt. 5307.45 .and. wn_io .lt. 5307.5 .and. n .eq. 1 */
/*     - .and. lpd .eq. 12 .and. nz .eq. 1) */
/*     - write(*,'(1a,i5,1pe14.6,14(1pe12.4))') */
/*     - 'perturb_rad: ',ibin,wn_io, */
/*     - tau_ext(l_1,1,n-npd1),tau_ext(l_1,3,n-npd1), */
/*     - dtauex(l_1,1),dtau(l_1), */
/*     - tau_ext(l_2,1,n-npd1),tau_ext(l_2,2,n-npd1), */
/*     - dtauex(l_2,1),dtau(l_2), */
/*     - copi0(l_1),copi0(l_2),tautest,dx(lpd,n) */

		    } else {

			if ((i__3 = ipd[n], abs(i__3)) == 4) {

/* ****                   define the cloud/aerosol optical depth change */

			    dx[lpd + n * 70] = pd_pert__[n] * dtauaer[modepd[
				    n - *npd1] + lpd * 10];

/* ****                   set levels for recalculating source functions */

			    l_1__ = lpd;
			    l_2__ = lpd;

/* *****                  calculate radiances for perturbed aerosol */
/*                       optical depth */

			    if (lpd < nlev) {

				dtau[lpd - 1] = dtauex[lpd + 70] + pd_pert__[
					n] * tau_ext__[lpd + ((n - *npd1) * 3 
					+ 1) * 70];

				if (dtau[lpd - 1] != 0.f) {

/* ****                       define the nominal scattering optical depth */

				    dtausca0 = (1.f - co_pi0__[lpd + 70]) * 
					    dtauex[lpd + 70];

/* ****                       perturb the scattering optical depth */

				    dtausca = dtausca0 + pd_pert__[n] * 
					    tau_sca__[lpd + (n - *npd1) * 70];

/* ****                      recompute perturbed single scattering albedo */

				    copi0[lpd - 1] = 1.f - dtausca / dtau[lpd 
					    - 1];

/* ****                        calculate  perturbed asymmetry parameter */

				    if (dtausca != 0.f) {
					g0[lpd - 1] = (g[lpd] * dtausca0 + 
						pd_pert__[n] * tau_sca__[lpd 
						+ (n - *npd1) * 70] * g_sca__[
						lpd + (n - *npd1) * 70]) / 
						dtausca;
				    } else {
					g0[lpd - 1] = 0.f;
				    }

				} else {
				    copi0[lpd - 1] = 0.f;
				    g0[lpd - 1] = 0.f;
				}

			    }

			} else {

/* ****                  set levels for recalculating source functions */

			    l_1__ = nlev;
			    l_2__ = nlev;

/* ****                   define the surface albedo change */

			    dx[lpd + n * 70] = pd_pert__[n] * alb[nz];

/* ***                   calculate radiances for perturbed albedos */

			    alb0 = (pd_pert__[n] + 1.f) * alb[nz];

			}

		    }

		}

	    }

/* ****           Layer transmittances are used in mapping fluxes and */
/*               radiances back to the high resolution spectrum. */
/*               find the transmission through each layer for each */
/*               stream. */

	    ng0 = *ibin;

	    if (*lsolar) {

/* ****             find the direct beam transmission for sun */

		i__3 = *nlyr;
		for (k = 1; k <= i__3; ++k) {
		    trndir[k - 1] = exp(-dtau[k - 1] / *umu0nz);
/* L3011: */
		}

	    }

/*      if(wn_io .gt. 6320.5 .and. wn_io .lt. 6320.6 .and. n .eq. 1 */
/*     - .and. lpd .eq. 10 .and. nz .eq. 1) */
/*     - write(*,'(/,1a,i5,1pe14.6,14(1pe12.4))') */
/*     - 'perturb_rad: ',ibin,wn_io, */
/*     - dtau(lpd),copi0(lpd),g0(lpd) */

/* ****           calculate radiances for this case */

	    if (*ibin > 0) {

		interp_rad__(lsolar, lplanck, &ng0, &l_1__, &l_2__, &nlev, 
			nz0, nphi, numu, nzup, nzdn, wn_io__, umu0nz, &umu[1],
			 &alb_b__[61], &dtau_b__[421], &copi0_b__[421], &
			g_b__[421], &trnflx_b__[421], &dtrnflxdx[281], &
			refflx_b__[421], &drefflxdx[281], &absflx_b__[421], &
			dabsflxdx[281], &refrad_b__[6737], &drefraddx[4497], &
			absrad_b__[6737], &dabsraddx[4497], &brdf_b__[1553], &
			dbrdfdx[273], &tsurf, tatm, &alb0, dtau, copi0, g0, &
			dalb, deltau, delpi0, delg, &dnsflxsrc_b__[71], &
			ddnsflxdx[351], &upsflxsrc_b__[71], &dupsflxdx[351], &
			dntflxsrc_b__[71], &ddntflxdx[351], &uptflxsrc_b__[71]
			, &duptflxdx[351], &sradsrc_b__[18193], &dsraddx[
			89873], &tradsrc_b__[18193], &dtraddx[89873], bb, &
			bb_sur__, bb_flx_up__, bb_flx_dn__, bb_rad__, 
			dn_s_src__, up_s_src__, rad_s_src__, dn_t_src__, 
			up_t_src__, rad_t_src__, trndir, trnflx, refflx, 
			absflx, trnrad, refrad, absrad, upsflx, dnsflx, 
			dirsflx, sol_rad__, uptflx, dntflx, th_rad__);

/* ****            initialize flux transmision values for simplified */
/*                adding sourse function integration */

		trnflx1_i__[lpd + (n + ni * 10) * 70] = trnflx[lpd - 1];
		refflx1_i__[lpd + (n + ni * 10) * 70] = refflx[lpd - 1];
		absflx1_i__[lpd + (n + ni * 10) * 70] = absflx[lpd - 1];
		if (*lsolar) {
		    if (lpd == 1) {
			trndir1_i__[lpd + (n + (ni + (nz << 1)) * 10) * 70] = 
				trndir[lpd - 1];
		    } else {
			trndir1_i__[lpd + (n + (ni + (nz << 1)) * 10) * 70] = 
				trndir[lpd - 1] * trndir1_i__[lpd - 1 + (n + (
				ni + (nz << 1)) * 10) * 70];
		    }
		}

	    } else {

/* ****             use the pure absorption model to calculate radiances */

		fbeam = 1.f;
		temis = 1.f;
		btemp = tsurf;
		ttemp = 0.f;

		rad_trn__(usrang, lplanck, &l_1__, &l_2__, &ng0, nstr, numu, &
			umu[1], &umu_f__[1], &gwt_f__[1], dtau, copi0, g0, 
			trnrad, absrad, trn_cmu__, abs_cmu__);

		grey_eq_trn__(lsolar, lplanck, usrang, lamber, &nlev, nz0, 
			nstr, numu, nphi, nzup, nzdn, iref, wn_io__, dtau, 
			copi0, tatm, &umu[1], &umu_f__[1], &gwt_f__[1], &phi[
			1], umu0nz, phi0nz, &fbeam, &alb0, &btemp, &ttemp, &
			temis, trndir, trnflx, absflx, trn_cmu__, abs_cmu__, 
			trnrad, absrad, upsflx, dnsflx, dirsflx, sol_rad__, 
			uptflx, dntflx, th_rad__, dn_s_src__, up_s_src__, 
			dn_t_src__, up_t_src__);

/* ****            initialize flux transmision values for simplified */
/*                adding sourse function integration */

		trnflx1_i__[lpd + (n + ni * 10) * 70] = trnflx[lpd - 1];
		refflx1_i__[lpd + (n + ni * 10) * 70] = refflx[lpd - 1];
		absflx1_i__[lpd + (n + ni * 10) * 70] = absflx[lpd - 1];

		if (*lsolar) {
		    if (lpd == 1) {
			trndir1_i__[lpd + (n + (ni + (nz << 1)) * 10) * 70] = 
				trndir[lpd - 1];
		    } else {
			trndir1_i__[lpd + (n + (ni + (nz << 1)) * 10) * 70] = 
				trndir[lpd - 1] * trndir1_i__[lpd - 1 + (n + (
				ni + (nz << 1)) * 10) * 70];
		    }
		}

	    }

/* ****           find the spectral gradients in the flux transmission and */
/*               absorption values used for the simplified adding method */

	    if (*lsolar) {
		dtrndir1_i__[lpd + (n + nz * 10) * 70] = (trndir1_i__[lpd + (
			n + ((nz << 1) + 2) * 10) * 70] - trndir1_i__[lpd + (
			n + ((nz << 1) + 1) * 10) * 70]) * *dvi;
	    }
	    dtrnflx1_i__[lpd + n * 70] = (trnflx1_i__[lpd + (n + 20) * 70] - 
		    trnflx1_i__[lpd + (n + 10) * 70]) * *dvi;
	    drefflx1_i__[lpd + n * 70] = (refflx1_i__[lpd + (n + 20) * 70] - 
		    refflx1_i__[lpd + (n + 10) * 70]) * *dvi;
	    dabsflx1_i__[lpd + n * 70] = (absflx1_i__[lpd + (n + 20) * 70] - 
		    absflx1_i__[lpd + (n + 10) * 70]) * *dvi;

/* ****             interpolate flux and radiance values to their */
/*                 output levels */

	    int_rad_lev__(lsolar, lplanck, nphi, numu, &nz, nlyr, nlout, 
		    levout, &k_out__[1], &dp_dp__[1], upsflx, dnsflx, dirsflx,
		     uptflx, dntflx, sol_rad__, th_rad__, ups, dns, dirs, 
		    upth, dnth, rad_s__, rad_th__);

/* ****             find the spectral gradients in the flux and radiance */


	    if (*lsolar) {

/* ****             define layer dependent solar source terms at wn_io */

		dn_s_src1__[lpd + (n + (ni + (nz << 1)) * 10) * 70] = 
			dn_s_src__[lpd - 1];
		up_s_src1__[lpd + (n + (ni + (nz << 1)) * 10) * 70] = 
			up_s_src__[lpd - 1];

/* ****           define solar source term wavelength derivatives */

		dup_s_src1__[lpd + (n + nz * 10) * 70] = (up_s_src1__[lpd + (
			n + ((nz << 1) + 2) * 10) * 70] - up_s_src1__[lpd + (
			n + ((nz << 1) + 1) * 10) * 70]) * *dvi;
		ddn_s_src1__[lpd + (n + nz * 10) * 70] = (dn_s_src1__[lpd + (
			n + ((nz << 1) + 2) * 10) * 70] - dn_s_src1__[lpd + (
			n + ((nz << 1) + 1) * 10) * 70]) * *dvi;

/* ****              define solar radiance wavelength derivatives */
/*                  at the output levels, k_out */

		i__3 = *nlout;
		for (k = 1; k <= i__3; ++k) {
		    i__4 = *nphi;
		    for (naz = 1; naz <= i__4; ++naz) {
			i__5 = *numu;
			for (nze = 1; nze <= i__5; ++nze) {

/* ****                         define solar radiances and wavelength */
/*                             derivatives at output levels, k_out */

			    sol_rad1__[nze + (naz + (k + (lpd + (n + (ni + (
				    nz << 1)) * 10) * 70) * 3 << 4) << 4)] = 
				    rad_s__[nze + (naz + (k << 4) << 4) - 273]
				    ;
			    dsol_rad1__[nze + (naz + (k + (lpd + (n + nz * 10)
				     * 70) * 3 << 4) << 4)] = (sol_rad1__[nze 
				    + (naz + (k + (lpd + (n + ((nz << 1) + 2) 
				    * 10) * 70) * 3 << 4) << 4)] - sol_rad1__[
				    nze + (naz + (k + (lpd + (n + ((nz << 1) 
				    + 1) * 10) * 70) * 3 << 4) << 4)]) * *dvi;
/* L3221: */
			}
/* L3241: */
		    }
/* L3261: */
		}
/*      if(wn_io .gt. 6320.5 .and. wn_io .lt. 6320.6 .and. n .eq. 1 */
/*     - .and. lpd .eq. 10 .and. nz .eq. 1) then */
/*      write(*,'(1a,16(1pe12.4))') 'sol_rad1(1)', */
/*     - (sol_rad1(numu,1,1,k,n,1,nz),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'sol_rad1(2)', */
/*     - (sol_rad1(numu,1,1,k,n,2,nz),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'dsol_rad1: ', */
/*     - (dsol_rad1(numu,1,1,k,n,nz),k=1,nlev) */
/*      endif */

	    }

	    if (*lplanck && nz == 1) {

/* ****             define layer dependent solar source terms at wn_io */

		dn_t_src1__[lpd + (n + ni * 10) * 70] = dn_t_src__[lpd - 1];
		up_t_src1__[lpd + (n + ni * 10) * 70] = up_t_src__[lpd - 1];

/* ****           define solar source term wavelength derivatives */

		dup_t_src1__[lpd + n * 70] = (up_t_src1__[lpd + (n + 20) * 70]
			 - up_t_src1__[lpd + (n + 10) * 70]) * *dvi;
		ddn_t_src1__[lpd + n * 70] = (dn_t_src1__[lpd + (n + 20) * 70]
			 - dn_t_src1__[lpd + (n + 10) * 70]) * *dvi;

/* ****              define thermal radiance wavelength derivatives */
/*                  at the output levels, k_out */

		i__3 = *nlout;
		for (k = 1; k <= i__3; ++k) {
		    i__4 = *nphi;
		    for (naz = 1; naz <= i__4; ++naz) {
			i__5 = *numu;
			for (nze = 1; nze <= i__5; ++nze) {

/* ****                         define solar radiances and wavelength */
/*                             derivatives at output levels, k_out */

			    th_rad1__[nze + (naz + (k + (lpd + (n + ni * 10) *
				     70) * 3 << 4) << 4)] = rad_th__[nze + (
				    naz + (k << 4) << 4) - 273];
			    dth_rad1__[nze + (naz + (k + (lpd + n * 70) * 3 <<
				     4) << 4)] = (th_rad1__[nze + (naz + (k + 
				    (lpd + (n + 20) * 70) * 3 << 4) << 4)] - 
				    th_rad1__[nze + (naz + (k + (lpd + (n + 
				    10) * 70) * 3 << 4) << 4)]) * *dvi;
/* L3421: */
			}
/* L3441: */
		    }
/* L3461: */
		}

	    }

/* L4001: */
	}

/* ****      exit loop over output levels */

/* L5001: */
    }

/* ****   exit loop over constituent */

    return 0;
} /* perturb_rad__ */

