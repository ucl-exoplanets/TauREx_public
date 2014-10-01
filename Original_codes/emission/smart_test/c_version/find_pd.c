/* find_pd.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int find_pd__(logical *lsolar, logical *lplanck, integer *
	nz0, integer *nlyr, integer *nlout, integer *nphi, integer *numu, 
	integer *npd, integer *ipd, doublereal *wn, doublereal *d_wn__, real *
	dx, real *solflx, real *rad, real *trn_dir__, real *trn_flx__, real *
	ref_flx__, real *abs_flx__, real *dns_src__, real *ups_src__, real *
	dnt_src__, real *upt_src__, doublereal *trndir1_i__, doublereal *
	trnflx1_i__, doublereal *refflx1_i__, doublereal *absflx1_i__, 
	doublereal *dtrndir1_i__, doublereal *dtrnflx1_i__, doublereal *
	drefflx1_i__, doublereal *dabsflx1_i__, doublereal *dn_s_src1__, 
	doublereal *up_s_src1__, doublereal *ddn_s_src1__, doublereal *
	dup_s_src1__, doublereal *dn_t_src1__, doublereal *up_t_src1__, 
	doublereal *ddn_t_src1__, doublereal *dup_t_src1__, doublereal *
	sol_rad1__, doublereal *dsol_rad1__, doublereal *th_rad1__, 
	doublereal *dth_rad1__, real *pd_trndir__, real *pd_trnflx__, real *
	pd_refflx__, real *pd_absflx__, real *pd_dns_src__, real *
	pd_ups_src__, real *pd_dnt_src__, real *pd_upt_src__, real *pd_rad__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static doublereal ups_src1__;
    static real upt_src1__;
    static integer k, n, nz, l_1__, l_2__, lpd;
    static real dxi[700]	/* was [70][10] */;
    static integer naz;
    static real dwn;
    static integer nze;
    static real rad1[256]	/* was [16][16] */;
    static integer nlev;
    static real rad_s1__;
    static doublereal absflx, refflx, trndir, trnflx;
    static real rad_th1__;
    static doublereal dns_src1__;
    static real dnt_src1__;


/* cccccccccccccccccccccccccc  f i n d _ p d  cccccccccccccccccccccccccccc */
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
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c   up_s_src - layer upward flux source function                     cc */
/* c   dn_s_src - layer downward flux source function                   cc */
/* c     th_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  f i n d _ p d  cccccccccccccccccccccccccccc */




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




/* ****    flux transmission and absorption functions at wavenuber, wn */



/* ****   flux source functions at wavenumber wn_io for perturbed state */


/* ****   flux source functions wavenumber derivatives */


/* ****     monochormatic radiance for basic state */


/* ****    flux transmission and absorption needed for simplified */
/*        adding method at wavenumber wnio */


/* ****     perturrbed monochormatic radiance and wn interpolation values */


/* ****    flux transmission and absorption partial */
/*        derivatives for simplified adding method at wavenumber wn */


/* ****   output flux and radiance partial derivatives */






/* ****    specify the solar zenith angle and gas index */

    /* Parameter adjustments */
    pd_rad__ -= 54801;
    pd_upt_src__ -= 71;
    pd_dnt_src__ -= 71;
    pd_ups_src__ -= 71;
    pd_dns_src__ -= 71;
    pd_absflx__ -= 71;
    pd_refflx__ -= 71;
    pd_trnflx__ -= 71;
    pd_trndir__ -= 71;
    dth_rad1__ -= 54801;
    th_rad1__ -= 592401;
    dsol_rad1__ -= 592401;
    sol_rad1__ -= 1667601;
    dup_t_src1__ -= 71;
    ddn_t_src1__ -= 71;
    up_t_src1__ -= 771;
    dn_t_src1__ -= 771;
    dup_s_src1__ -= 771;
    ddn_s_src1__ -= 771;
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
    --upt_src__;
    --dnt_src__;
    --ups_src__;
    --dns_src__;
    --abs_flx__;
    --ref_flx__;
    --trn_flx__;
    --trn_dir__;
    rad -= 273;
    dx -= 71;
    --ipd;

    /* Function Body */
    nz = *nz0;
    nlev = *nlyr + 1;
    dwn = (real) (*d_wn__);

/* ****       evaluate the radiances and partial derivatives */
/*           at this wavenumber for each optical property */
/*           that can vary */

    i__1 = *npd;
    for (n = 1; n <= i__1; ++n) {

/* ****        determine the number of levels at which the state */
/*            vector variable chages (T=nlev, tau=nlyr, alb=1) */

	if ((i__2 = ipd[n], abs(i__2)) == 1) {

/* ****         surface pressure */

	    l_1__ = nlev;
	    l_2__ = nlev;

	} else {

	    if ((i__2 = ipd[n], abs(i__2)) == 2) {

/* ****            atmospheric/surface temperature */

		l_1__ = 1;
		l_2__ = nlev;
	    } else {

		if ((i__2 = ipd[n], abs(i__2)) == 3 || (i__3 = ipd[n], abs(
			i__3)) == 4) {

/* ****             gas or aerosol optical depth */

		    l_1__ = 1;
		    l_2__ = nlev;

		} else {

/* ****             this is surface albedo */

		    l_1__ = nlev;
		    l_2__ = nlev;
		}
	    }
	}

/* ****      enter the loop over perturbation level */

	i__2 = l_2__;
	for (lpd = l_1__; lpd <= i__2; ++lpd) {

/* *****         define the inverse of the state vector change */

	    if (dx[lpd + n * 70] != 0.f) {
		dxi[lpd + n * 70 - 71] = 1.f / dx[lpd + n * 70];
	    } else {
		dxi[lpd + n * 70 - 71] = 0.f;
	    }

/* ****           interpolate flux transmision and absorption functions */
/*               and jacobians to this wavenumber */

	    if (lpd <= nlev) {
		if (ipd[n] < 0) {
		    trnflx = trnflx1_i__[lpd + (n + 10) * 70] + dtrnflx1_i__[
			    lpd + n * 70] * *d_wn__;
		    refflx = refflx1_i__[lpd + (n + 10) * 70] + drefflx1_i__[
			    lpd + n * 70] * *d_wn__;
		    absflx = absflx1_i__[lpd + (n + 10) * 70] + dabsflx1_i__[
			    lpd + n * 70] * *d_wn__;

		    pd_trnflx__[lpd + n * 70] = (real) ((trnflx - trn_flx__[
			    lpd]) * dxi[lpd + n * 70 - 71]);
		    pd_refflx__[lpd + n * 70] = (real) ((refflx - ref_flx__[
			    lpd]) * dxi[lpd + n * 70 - 71]);
		    pd_absflx__[lpd + n * 70] = (real) ((absflx - abs_flx__[
			    lpd]) * dxi[lpd + n * 70 - 71]);
		} else {
		    pd_trnflx__[lpd + n * 70] = 0.f;
		    pd_refflx__[lpd + n * 70] = 0.f;
		    pd_absflx__[lpd + n * 70] = 0.f;
		}
	    } else {
		pd_trnflx__[lpd + n * 70] = 0.f;
		pd_refflx__[lpd + n * 70] = 0.f;
		pd_absflx__[lpd + n * 70] = 0.f;
	    }

	    if (l_2__ <= nlev) {
		if (*lsolar) {

		    trndir = trndir1_i__[lpd + (n + ((nz << 1) + 1) * 10) * 
			    70] + dtrndir1_i__[lpd + (n + nz * 10) * 70] * *
			    d_wn__;

		    pd_trndir__[lpd + n * 70] = (real) ((trndir - trn_dir__[
			    lpd]) * dxi[lpd + n * 70 - 71]);

		    if (ipd[n] < 0) {

/* *****                  interpolate perturbed layer dependent solar */
/*                       source terms to wn */

			dns_src1__ = dn_s_src1__[lpd + (n + ((nz << 1) + 1) * 
				10) * 70] + ddn_s_src1__[lpd + (n + nz * 10) *
				 70] * *d_wn__;
			ups_src1__ = up_s_src1__[lpd + (n + ((nz << 1) + 1) * 
				10) * 70] + dup_s_src1__[lpd + (n + nz * 10) *
				 70] * *d_wn__;

/* *****              find jacobians for solar source functions */

			pd_dns_src__[lpd + n * 70] = (real) (dns_src1__ - 
				dns_src__[lpd]) * dxi[lpd + n * 70 - 71];
			pd_ups_src__[lpd + n * 70] = (real) (ups_src1__ - 
				ups_src__[lpd]) * dxi[lpd + n * 70 - 71];

		    } else {

			pd_dns_src__[lpd + n * 70] = 0.f;
			pd_ups_src__[lpd + n * 70] = 0.f;

		    }

		}

		if (*lplanck) {

		    if (ipd[n] < 0) {

/* *****                  interpolate perturbed layer dependent solar */
/*                       source terms to wn */

			dnt_src1__ = (real) (dn_t_src1__[lpd + (n + 10) * 70] 
				+ ddn_t_src1__[lpd + n * 70] * *d_wn__);
			upt_src1__ = (real) (up_t_src1__[lpd + (n + 10) * 70] 
				+ dup_t_src1__[lpd + n * 70] * *d_wn__);

/* *****              find jacobians for thermal source functions */

			pd_dnt_src__[lpd + n * 70] = (dnt_src1__ - dnt_src__[
				lpd]) * dxi[lpd + n * 70 - 71];
			pd_upt_src__[lpd + n * 70] = (upt_src1__ - upt_src__[
				lpd]) * dxi[lpd + n * 70 - 71];

		    } else {

			pd_dnt_src__[lpd + n * 70] = 0.f;
			pd_upt_src__[lpd + n * 70] = 0.f;

		    }

		}

	    } else {

		if (*lsolar) {

		    pd_dns_src__[lpd + n * 70] = 0.f;
		    pd_ups_src__[lpd + n * 70] = 0.f;
		    pd_trndir__[lpd + n * 70] = 0.f;

		}
		if (*lplanck) {

		    pd_dns_src__[lpd + n * 70] = 0.f;
		    pd_ups_src__[lpd + n * 70] = 0.f;

		}

	    }

/* ****          interpolate perturbed radiances to wavenumber, wn, */
/*              and find partial derivatives */

	    if (ipd[n] > 0) {
		i__3 = *nlout;
		for (k = 1; k <= i__3; ++k) {
		    i__4 = *nphi;
		    for (naz = 1; naz <= i__4; ++naz) {
			i__5 = *numu;
			for (nze = 1; nze <= i__5; ++nze) {

/* *****                         interpolate solar radiance to wn */

			    rad_s1__ = 0.f;
			    if (*lsolar) {
				rad_s1__ = (real) (*solflx * (sol_rad1__[nze 
					+ (naz + (k + (lpd + (n + ((nz << 1) 
					+ 1) * 10) * 70) * 3 << 4) << 4)] + 
					dsol_rad1__[nze + (naz + (k + (lpd + (
					n + nz * 10) * 70) * 3 << 4) << 4)] * 
					dwn));
			    }

/* *****                       interpolate thermal radiance to wn */

			    rad_th1__ = 0.f;
			    if (*lplanck) {
				rad_th1__ = (real) (th_rad1__[nze + (naz + (k 
					+ (lpd + (n + 10) * 70) * 3 << 4) << 
					4)] + dth_rad1__[nze + (naz + (k + (
					lpd + n * 70) * 3 << 4) << 4)] * dwn);
			    }

/* ****                         define the total (solar+thermal) radiance */

			    rad1[nze + (naz << 4) - 17] = rad_s1__ + 
				    rad_th1__;

/* ****                         find radiance jacobian */

			    pd_rad__[nze + (naz + (k + (lpd + n * 70) * 3 << 
				    4) << 4)] = (rad1[nze + (naz << 4) - 17] 
				    - rad[nze + (naz + (k << 4) << 4)]) * dxi[
				    lpd + n * 70 - 71];
/* L4021: */
			}
/* L4041: */
		    }
/* L4061: */
		}

	    }

/* L4981: */
	}

/* ****      exit loop over output levels */

/* L5001: */
    }

/* ****   exit loop over constituent */

    return 0;
} /* find_pd__ */

