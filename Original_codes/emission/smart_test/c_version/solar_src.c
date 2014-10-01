/* solar_src.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int solar_src__(integer *ng0, integer *l0, integer *nz0, 
	integer *nlyr, integer *nzdn, integer *nzup, integer *numu, integer *
	nphi, integer *n_rad__, real *alb0, real *umu0, real *rfldir, real *
	rfldn, real *dx_b_i__, real *dalb_b_i__, doublereal *trnrad_b__, 
	doublereal *rad_rt__, doublereal *rad_rb__, doublereal *flx_ft__, 
	doublereal *flx_fb__, doublereal *trnflx_b__, doublereal *flx_rt__, 
	doublereal *flx_rb__, doublereal *flx_du__, doublereal *flx_dd__, 
	doublereal *flx_dft__, doublereal *flx_uft__, doublereal *flx_dfb__, 
	doublereal *flx_ufb__, doublereal *dnflxsrc, doublereal *upflxsrc, 
	doublereal *radsrc, doublereal *dnsflxsrc, doublereal *upsflxsrc, 
	doublereal *sradsrc, doublereal *dnsflxsrc_b__, doublereal *
	upsflxsrc_b__, doublereal *sradsrc_b__, doublereal *ddnsflxdx, 
	doublereal *dupsflxdx, doublereal *dsraddx)
{
    /* Initialized data */

    static doublereal pi = 3.141592654;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, l, nz, ibn, naz, nze;
    static real dflx;
    static integer nlev;
    static real uflx;
    static doublereal test1[70], test2[70];


/* ccccccccccccccccccccccccc  s o l a r _ s r c  ccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine scales binned solar radiance and flux sources   cc */
/* c    and then loads them and and their jacobians into the arrays     cc */
/* c    used by interp_rad                                              cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      rfldir - downward direct solar flux at each level for bin     cc */
/* c dnsflxsrc - downward flux source at each level for this bin        cc */
/* c upsflxsrc - downward flux source at each level for this bin        cc */
/* c   sradsrc - angle-dependent radiance source at each level for      cc */
/* c               this bin                                             cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c dnsflxsrc_b - downward solar flux source at each level for bin     cc */
/* c upsflxsrc_b - downward solar flux source at each level for bin     cc */
/* c   sradsrc_b - angle-dependent solar radiance source at each level  cc */
/* c              for this bin                                          cc */
/* c  ddnsflxdx - downward solar flux jacobian at each level for  bin   cc */
/* c  dupsflxdx - upward solar flux jacobian at each level for bin      cc */
/* c    dsraddx - solar radiance jacobian at each level for bin         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  s o l a r _ s r c  ccccccccccccccccccccccccc */




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



/* ****   internal binned source function variables */


/* ****   binned solar source function variables and jacobians */


/* ****   source_fcn output variables */


/* ****   scaled solar source function variables and jacobians */



/* ****  quantities for radiance calculation: downward and upward fluxes */
/*      at layer interfaces */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   variables for radiance layer adding method */




/* ****   define pi */

    /* Parameter adjustments */
    dsraddx -= 89873;
    dupsflxdx -= 351;
    ddnsflxdx -= 351;
    sradsrc_b__ -= 18193;
    upsflxsrc_b__ -= 71;
    dnsflxsrc_b__ -= 71;
    sradsrc -= 18193;
    upsflxsrc -= 71;
    dnsflxsrc -= 71;
    radsrc -= 18193;
    upflxsrc -= 71;
    dnflxsrc -= 71;
    flx_ufb__ -= 71;
    flx_dfb__ -= 71;
    flx_uft__ -= 71;
    flx_dft__ -= 71;
    flx_dd__ -= 71;
    flx_du__ -= 71;
    flx_rb__ -= 71;
    flx_rt__ -= 71;
    trnflx_b__ -= 421;
    flx_fb__ -= 71;
    flx_ft__ -= 71;
    rad_rb__ -= 18193;
    rad_rt__ -= 18193;
    trnrad_b__ -= 6737;
    dalb_b_i__ -= 11;
    dx_b_i__ -= 421;
    --rfldn;
    --rfldir;
    n_rad__ -= 6;

    /* Function Body */

    l = *l0;
    ibn = *ng0;
    nz = *nz0;
    nlev = *nlyr + 1;

    if (l == 1 || n_rad__[l + ibn * 5] != 0) {

/* ****     set values of fluxes at top and bottom of atmosphere. */

	dnsflxsrc[l * 70 + 1] = 0.f;
	upsflxsrc[nlev + l * 70] = *alb0 * rfldir[nlev];

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {

/* ****         define the flux scaling factor */

	    dflx = rfldir[k] / *umu0 + flx_dft__[k + l * 70];

	    if (dflx > 0.f) {

		test1[k] = dnflxsrc[k + 1 + l * 70] / (rfldir[k] / *umu0 + (
			flx_ft__[k + l * 70] + flx_rt__[k + l * 70] * 
			upflxsrc[k + l * 70]) * flx_dd__[k + l * 70]);
		uflx = rfldir[k + 1] + flx_ufb__[k + 1 + l * 70];
		test2[k] = dnflxsrc[k + 1 + l * 70] / (rfldir[k] / *umu0 + (
			flx_ft__[k + l * 70] + flx_rt__[k + l * 70] * 
			upflxsrc[k + l * 70] / uflx) * flx_dd__[k + l * 70]);
		dnsflxsrc[k + 1 + l * 70] = dnflxsrc[k + 1 + l * 70] / dflx;

	    } else {

		dnsflxsrc[k + 1 + l * 70] = 0.f;

	    }

	    uflx = rfldir[k + 1] + flx_ufb__[k + 1 + l * 70];

	    if (uflx != 0.f) {

		upsflxsrc[k + l * 70] = upflxsrc[k + l * 70] / uflx;

	    } else {

		upsflxsrc[k + l * 70] = 0.f;

	    }

/* L1201: */
	}

/* ****      scale radiances */
/*             note: the upward radiances at the surface (k=nlev) are */
/*                   not scaled since trnrad_b(nlev) = 0.0 */

/* ****     set values of radiances at top and bottom of atmosphere */
/*         where they are not scaled */

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {
	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {
		sradsrc[nze + (naz + (l * 70 + 1 << 4) << 4)] = radsrc[nze + (
			naz + (l * 70 + 1 << 4) << 4)];
/* L2001: */
	    }
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		sradsrc[nze + (naz + (nlev + l * 70 << 4) << 4)] = radsrc[nze 
			+ (naz + (nlev + l * 70 << 4) << 4)];
/* L2021: */
	    }
/* L2041: */
	}

/* ****     scale values at other layers */

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
		    sradsrc[nze + (naz + (k + 1 + l * 70 << 4) << 4)] = 
			    radsrc[nze + (naz + (k + 1 + l * 70 << 4) << 4)] 
			    + trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] 
			    * rad_rt__[nze + (naz + (k + l * 70 << 4) << 4)] *
			     flx_uft__[k + l * 70] / pi;
/* L2201: */
		}

		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
		    sradsrc[nze + (naz + (k + l * 70 << 4) << 4)] = radsrc[
			    nze + (naz + (k + l * 70 << 4) << 4)] + 
			    trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] * 
			    rad_rb__[nze + (naz + (k + 1 + l * 70 << 4) << 4)]
			     * flx_dfb__[k + 1 + l * 70] / pi;
/* L2221: */
		}
/* L2241: */
	    }
/* L2261: */
	}

    }

/* ****  store flux and radiance source variables in bin arrays */

    if (*l0 == 1) {

/* ****           load nominal flux and radiance source terms */

	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {

	    dnsflxsrc_b__[k + ibn * 70] = dnsflxsrc[k + 70];
	    upsflxsrc_b__[k + ibn * 70] = upsflxsrc[k + 70];

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = 1; nze <= i__3; ++nze) {
		    sradsrc_b__[nze + (naz + (k + ibn * 70 << 4) << 4)] = 
			    sradsrc[nze + (naz + (k + 70 << 4) << 4)];
/* L2601: */
		}
/* L2621: */
	    }
/* L2641: */
	}

    } else {

	if (n_rad__[l + ibn * 5] != 0) {

/* ****         flux derivatives at top of atmosphere and surface */

	    ddnsflxdx[(l - 1 + (ibn << 2)) * 70 + 1] = 0.;

	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {

/* ****         downward radiance derivatives at top of atmosphere */

		i__2 = *nzdn;
		for (nze = 1; nze <= i__2; ++nze) {
		    dsraddx[nze + (naz + ((l - 1 + (ibn << 2)) * 70 + 1 << 4) 
			    << 4)] = 0.f;
/* L2721: */
		}

/* ****           upward radiance streams derivatives at the surface */

		if (l < 5) {
		    i__2 = *numu;
		    for (nze = *nzup; nze <= i__2; ++nze) {
			dsraddx[nze + (naz + (nlev + (l - 1 + (ibn << 2)) * 
				70 << 4) << 4)] = 0.f;
/* L2741: */
		    }

		} else {

		    i__2 = *numu;
		    for (nze = *nzup; nze <= i__2; ++nze) {
			dsraddx[nze + (naz + (nlev + (l - 1 + (ibn << 2)) * 
				70 << 4) << 4)] = dalb_b_i__[nz + ibn * 10] * 
				(sradsrc[nze + (naz + (nlev + l * 70 << 4) << 
				4)] - sradsrc[nze + (naz + (nlev + 70 << 4) <<
				 4)]);
/* L2761: */
		    }
		}
/* L2781: */
	    }

/* ***        find flux and radiance derivatives at other levels */

	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {

/* ****          flux derivatives each level */

		ddnsflxdx[k + 1 + (l - 1 + (ibn << 2)) * 70] = dx_b_i__[k + (
			l + ibn * 5) * 70] * (dnsflxsrc[k + 1 + l * 70] - 
			dnsflxsrc[k + 71]);
		dupsflxdx[k + (l - 1 + (ibn << 2)) * 70] = dx_b_i__[k + (l + 
			ibn * 5) * 70] * (upsflxsrc[k + l * 70] - upsflxsrc[k 
			+ 70]);

		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {

/* ****             downward radiance derivatives at each level */

		    i__3 = *nzdn;
		    for (nze = 1; nze <= i__3; ++nze) {
			dsraddx[nze + (naz + (k + 1 + (l - 1 + (ibn << 2)) * 
				70 << 4) << 4)] = dx_b_i__[k + (l + ibn * 5) *
				 70] * (sradsrc[nze + (naz + (k + 1 + l * 70 
				<< 4) << 4)] - sradsrc[nze + (naz + (k + 71 <<
				 4) << 4)]);
/* L2801: */
		    }

/* ****                 upward radiance streams derivatives at  surface */

		    i__3 = *numu;
		    for (nze = *nzup; nze <= i__3; ++nze) {
			dsraddx[nze + (naz + (k + (l - 1 + (ibn << 2)) * 70 <<
				 4) << 4)] = dx_b_i__[k + (l + ibn * 5) * 70] 
				* (sradsrc[nze + (naz + (k + l * 70 << 4) << 
				4)] - sradsrc[nze + (naz + (k + 70 << 4) << 4)
				]);
/* L2821: */
		    }
/* L2841: */
		}
/* L2861: */
	    }

	} else {

/* ****       n_rad = 0, not a variable part of the state vector */

	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {
		ddnsflxdx[k + (l - 1 + (ibn << 2)) * 70] = 0.;
		dupsflxdx[k + (l - 1 + (ibn << 2)) * 70] = 0.;
		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__3 = *numu;
		    for (nze = 1; nze <= i__3; ++nze) {
			dsraddx[nze + (naz + (k + (l - 1 + (ibn << 2)) * 70 <<
				 4) << 4)] = 0.;
/* L2901: */
		    }
/* L2921: */
		}
/* L2941: */
	    }

	}

    }

    return 0;
} /* solar_src__ */

