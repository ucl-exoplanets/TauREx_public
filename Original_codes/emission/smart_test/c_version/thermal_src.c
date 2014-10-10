/* thermal_src.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int thermal_src__(integer *ng0, integer *l0, integer *nz0, 
	integer *nlyr, integer *nzdn, integer *nzup, integer *numu, integer *
	nphi, integer *n_rad__, real *wng0, real *alb0, doublereal *brdf_b__, 
	real *ts, real *t, real *dx_b_i__, real *dalb_b_i__, doublereal *
	trnrad_b__, doublereal *rad_rt__, doublereal *rad_rb__, doublereal *
	flx_uft__, doublereal *flx_dfb__, doublereal *dnflxsrc, doublereal *
	upflxsrc, doublereal *radsrc, doublereal *dntflxsrc, doublereal *
	uptflxsrc, doublereal *tradsrc, real *copi0_b__, doublereal *
	absflx_b__, doublereal *absrad_b__, doublereal *dntflxsrc_b__, 
	doublereal *uptflxsrc_b__, doublereal *tradsrc_b__, doublereal *
	ddntflxdx, doublereal *duptflxdx, doublereal *dtraddx, doublereal *
	bb_rad__)
{
    /* Initialized data */

    static doublereal pi = 3.141592654;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, l;
    static doublereal bb_flx_dn__[70], bb_flx_up__[70], bb[70];
    static integer nz, ibn, naz, nze, nlev;
    static real tsurf0[1];
    extern /* Subroutine */ int planck_(integer *, integer *, integer *, 
	    integer *, real *, real *, doublereal *);
    static doublereal bb_sur__;
    static integer iunits;


/* cccccccccccccccccccccccc  t h e r m a l _ s r c  cccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine scales binned thermal radiance and flux sources cc */
/* c    and then loads them and and their jacobians into the arrays     cc */
/* c    used by interp_rad                                              cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c         ng0 - spectral bin index                                   cc */
/* c   dntflxsrc - downward flux source at each level for this bin      cc */
/* c   uptflxsrc - downward flux source at each level for this bin      cc */
/* c     tradsrc - angle-dependent radiance source at each level for    cc */
/* c               this bin                                             cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c dntflxsrc_b - downward thermal flux source at each level for bin   cc */
/* c uptflxsrc_b - downward thermal flux source at each level for bin   cc */
/* c   tradsrc_b - angle-dependent thermal radiance source at each      cc */
/* c              level for this bin                                    cc */
/* c  ddntflxdx - downward thermal flux jacobian at each level for bin  cc */
/* c  duptflxdx - upward thermal flux jacobian at each level for bin    cc */
/* c    dtraddx - thermal radiance jacobian at each level for bin       cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccc  t h e r m a l _ s r c  cccccccccccccccccccccc */




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




/* ****   source_fcn output variables */


/* ****   scaled thermal source function variables */


/* ****   internal black body sources */


/* ****   binned thermal source function variables and jacobians */



/* ****  quantities for radiance calculation: downward and upward fluxes */
/*      at layer interfaces */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   variables for radiance layer adding method */





    /* Parameter adjustments */
    bb_rad__ -= 273;
    dtraddx -= 89873;
    duptflxdx -= 351;
    ddntflxdx -= 351;
    tradsrc_b__ -= 18193;
    uptflxsrc_b__ -= 71;
    dntflxsrc_b__ -= 71;
    absrad_b__ -= 6737;
    absflx_b__ -= 421;
    copi0_b__ -= 421;
    tradsrc -= 18193;
    uptflxsrc -= 71;
    dntflxsrc -= 71;
    radsrc -= 18193;
    upflxsrc -= 71;
    dnflxsrc -= 71;
    flx_dfb__ -= 71;
    flx_uft__ -= 71;
    rad_rb__ -= 18193;
    rad_rt__ -= 18193;
    trnrad_b__ -= 6737;
    dalb_b_i__ -= 11;
    dx_b_i__ -= 421;
    --t;
    brdf_b__ -= 1553;
    n_rad__ -= 6;

    /* Function Body */

    l = *l0;
    ibn = *ng0;
    nz = *nz0;
    nlev = *nlyr + 1;
    iunits = 1;

    if (l == 1 || n_rad__[l + ibn * 5] != 0) {

/* ****   evaluate planck function at the surface */

	tsurf0[0] = *ts;

	planck_(&c__1, &c__1, &c__1, &iunits, wng0, tsurf0, bb);

	bb_sur__ = bb[0];

/* ****    find the atmospheric planck functions at this wavenumber */

	planck_(&c__1, &nlev, &c__1, &iunits, wng0, &t[1], bb);

/* ****         define the black body source functions for flux and */
/*             radiance, assuming a linear-in-tau formulation. */

	bb_flx_dn__[0] = pi * bb[0] * copi0_b__[(l + ibn * 5) * 70 + 1] * 
		absflx_b__[(l + ibn * 5) * 70 + 1];
	bb_flx_up__[nlev - 1] = pi * bb_sur__ * (1. - *alb0);

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    bb_flx_dn__[k] = pi * .5 * (bb[k - 1] + bb[k]) * copi0_b__[k + (l 
		    + ibn * 5) * 70] * absflx_b__[k + (l + ibn * 5) * 70];
	    bb_flx_up__[k - 1] = pi * .5 * (bb[k - 1] + bb[k]) * copi0_b__[k 
		    + (l + ibn * 5) * 70] * absflx_b__[k + (l + ibn * 5) * 70]
		    ;
/* L1001: */
	}

/* ****    thermal sources for downward radiance streams at top of */
/*        atmosphere and upward radiance streams at surface */

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {
	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {

/* ****            define a non-zero thermal source for scaling top level */

		bb_rad__[nze + (naz + 16 << 4)] = bb[0] * copi0_b__[(l + ibn *
			 5) * 70 + 1] * absrad_b__[nze + ((l + ibn * 5) * 70 
			+ 1 << 4)];
/* L1201: */
	    }

/* ****         define the surface blackbody emission */

	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		bb_rad__[nze + (naz + (nlev << 4) << 4)] = bb_sur__ * (1. - 
			brdf_b__[nze + (naz + (l + ibn * 5 << 4) << 4)]);
/* L1221: */
	    }
/* L1241: */
	}

/* ****             thermal source for downward streams */

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    naz = 1;
	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {
		bb_rad__[nze + (naz + (k + 1 << 4) << 4)] = (bb[k - 1] + bb[k]
			) * .5 * copi0_b__[k + (l + ibn * 5) * 70] * 
			absrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)];
/* L1401: */
	    }

/* ****         set values at all other azimuth angles */
/*             (no azimuth dependence above surface) */

	    i__2 = *nphi;
	    for (naz = 2; naz <= i__2; ++naz) {
		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
		    bb_rad__[nze + (naz + (k + 1 << 4) << 4)] = bb_rad__[nze 
			    + ((k + 1 << 4) + 1 << 4)];
/* L1421: */
		}
/* L1441: */
	    }
/* L1461: */
	}

/* ****       thermal source for upward streams */

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    naz = 1;
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		bb_rad__[nze + (naz + (k << 4) << 4)] = (bb[k - 1] + bb[k]) * 
			.5 * copi0_b__[k + (l + ibn * 5) * 70] * absrad_b__[
			nze + (k + (l + ibn * 5) * 70 << 4)];
/* L1601: */
	    }
	    i__2 = *nphi;
	    for (naz = 2; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
		    bb_rad__[nze + (naz + (k << 4) << 4)] = bb_rad__[nze + ((
			    k << 4) + 1 << 4)];
/* L1621: */
		}
/* L1641: */
	    }
/* L1661: */
	}

/* ****     normalize the fluxes and radiances by the black body */
/*         emission in the layer. */

	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {

/* ****          scale downward fluxes */

	    if (bb_flx_dn__[k - 1] > 0.) {
		dntflxsrc[k + l * 70] = (dnflxsrc[k + l * 70] - bb_flx_dn__[k 
			- 1]) / bb_flx_dn__[k - 1];
	    } else {
		dntflxsrc[k + l * 70] = 0.;
	    }
/*              dntflxsrc(k,l) = dnflxsrc(k,l) */

/* ****         scale upward fluxes */

	    if (bb_flx_up__[k - 1] > 0.) {
		uptflxsrc[k + l * 70] = (upflxsrc[k + l * 70] - bb_flx_up__[k 
			- 1]) / bb_flx_up__[k - 1];
	    } else {
		uptflxsrc[k + l * 70] = 0.;
	    }
/*            uptflxsrc(k,l) = upflxsrc(k,l) */
/* L2001: */
	}

/* ****     scale radiance source terms at top of atmosphere and surface */

	i__1 = *nphi;
	for (naz = 1; naz <= i__1; ++naz) {

/* ****          scale downward radiance sources at top of atmosphere */

	    k = 1;
	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {
		tradsrc[nze + (naz + (k + l * 70 << 4) << 4)] = radsrc[nze + (
			naz + (k + l * 70 << 4) << 4)] - bb_rad__[nze + (naz 
			+ (k << 4) << 4)];
/* L2021: */
	    }

/* ****              scale upward radiance sources at the surface */

	    k = nlev;
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		tradsrc[nze + (naz + (k + l * 70 << 4) << 4)] = radsrc[nze + (
			naz + (k + l * 70 << 4) << 4)] - bb_rad__[nze + (naz 
			+ (k << 4) << 4)];
/* L2041: */
	    }
/* L2061: */
	}

/* ****     scale values at other levels */

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {

/* ****              scale downward radiance sources */

		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
		    tradsrc[nze + (naz + (k + 1 + l * 70 << 4) << 4)] = 
			    radsrc[nze + (naz + (k + 1 + l * 70 << 4) << 4)] 
			    + trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] 
			    * rad_rt__[nze + (naz + (k + l * 70 << 4) << 4)] *
			     flx_uft__[k + l * 70] / pi - bb_rad__[nze + (naz 
			    + (k << 4) << 4)];
/* L2101: */
		}

/* ****              scale upward radiance sources */

		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
		    tradsrc[nze + (naz + (k + l * 70 << 4) << 4)] = radsrc[
			    nze + (naz + (k + l * 70 << 4) << 4)] + 
			    trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] * 
			    rad_rb__[nze + (naz + (k + 1 + l * 70 << 4) << 4)]
			     * flx_dfb__[k + 1 + l * 70] / pi - bb_rad__[nze 
			    + (naz + (k << 4) << 4)];
/* L2121: */
		}
/* L2141: */
	    }

/* L2161: */
	}

    }

    if (*l0 == 1) {

/* ****           load nominal flux and radiance source terms */

	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {

	    dntflxsrc_b__[k + ibn * 70] = dntflxsrc[k + 70];
	    uptflxsrc_b__[k + ibn * 70] = uptflxsrc[k + 70];

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = 1; nze <= i__3; ++nze) {
		    tradsrc_b__[nze + (naz + (k + ibn * 70 << 4) << 4)] = 
			    tradsrc[nze + (naz + (k + 70 << 4) << 4)];
/* L2201: */
		}
/* L2221: */
	    }
/* L2241: */
	}

    } else {

	if (n_rad__[l + ibn * 5] != 0) {

/* ****         flux derivatives at top of atmosphere and surface */

	    ddntflxdx[(l - 1 + (ibn << 2)) * 70 + 1] = 0.;
	    duptflxdx[nlev + (l - 1 + (ibn << 2)) * 70] = dx_b_i__[nlev + (l 
		    + ibn * 5) * 70] * (uptflxsrc[nlev + l * 70] - uptflxsrc[
		    nlev + 70]);

	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {

/* ****         downward radiance derivatives at top of atmosphere */

		i__2 = *nzdn;
		for (nze = 1; nze <= i__2; ++nze) {
		    dtraddx[nze + (naz + ((l - 1 + (ibn << 2)) * 70 + 1 << 4) 
			    << 4)] = 0.f;
/* L2401: */
		}

/* ****           upward radiance streams derivatives at the surface */

		if (l < 5) {
		    i__2 = *numu;
		    for (nze = *nzup; nze <= i__2; ++nze) {
			dtraddx[nze + (naz + (nlev + (l - 1 + (ibn << 2)) * 
				70 << 4) << 4)] = 0.f;
/* L2421: */
		    }

		} else {

		    i__2 = *numu;
		    for (nze = *nzup; nze <= i__2; ++nze) {
			dtraddx[nze + (naz + (nlev + (l - 1 + (ibn << 2)) * 
				70 << 4) << 4)] = dalb_b_i__[nz + ibn * 10] * 
				(tradsrc[nze + (naz + (nlev + l * 70 << 4) << 
				4)] - tradsrc[nze + (naz + (nlev + 70 << 4) <<
				 4)]);
/* L2441: */
		    }

		}

/* L2461: */
	    }

/* ***        find flux and radiance derivatives at other levels */

	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {

/* ****          flux derivatives each level */

		ddntflxdx[k + 1 + (l - 1 + (ibn << 2)) * 70] = dx_b_i__[k + (
			l + ibn * 5) * 70] * (dntflxsrc[k + 1 + l * 70] - 
			dntflxsrc[k + 71]);
		duptflxdx[k + (l - 1 + (ibn << 2)) * 70] = dx_b_i__[k + (l + 
			ibn * 5) * 70] * (uptflxsrc[k + l * 70] - uptflxsrc[k 
			+ 70]);

		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {

/* ****             downward radiance derivatives at each level */

		    i__3 = *nzdn;
		    for (nze = 1; nze <= i__3; ++nze) {
			dtraddx[nze + (naz + (k + 1 + (l - 1 + (ibn << 2)) * 
				70 << 4) << 4)] = dx_b_i__[k + (l + ibn * 5) *
				 70] * (tradsrc[nze + (naz + (k + 1 + l * 70 
				<< 4) << 4)] - tradsrc[nze + (naz + (k + 71 <<
				 4) << 4)]);
/* L2601: */
		    }

/* ****                 upward radiance streams derivatives at  surface */

		    i__3 = *numu;
		    for (nze = *nzup; nze <= i__3; ++nze) {
			dtraddx[nze + (naz + (k + (l - 1 + (ibn << 2)) * 70 <<
				 4) << 4)] = dx_b_i__[k + (l + ibn * 5) * 70] 
				* (tradsrc[nze + (naz + (k + l * 70 << 4) << 
				4)] - tradsrc[nze + (naz + (k + 70 << 4) << 4)
				]);
/* L2621: */
		    }
/* L2641: */
		}
/* L2661: */
	    }

	} else {

/* ****       n_rad = 0, not a variable part of the state vector */

	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {
		ddntflxdx[k + (l - 1 + (ibn << 2)) * 70] = 0.;
		duptflxdx[k + (l - 1 + (ibn << 2)) * 70] = 0.;
		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__3 = *numu;
		    for (nze = 1; nze <= i__3; ++nze) {
			dtraddx[nze + (naz + (k + (l - 1 + (ibn << 2)) * 70 <<
				 4) << 4)] = 0.;
/* L2801: */
		    }
/* L2821: */
		}
/* L2841: */
	    }

	}

    }

    return 0;
} /* thermal_src__ */

