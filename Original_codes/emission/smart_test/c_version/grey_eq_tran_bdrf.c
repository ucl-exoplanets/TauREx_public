/* grey_eq_tran_bdrf.f -- translated by f2c (version 20100827).
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

static integer c__8 = 8;
static integer c__1 = 1;

/* Subroutine */ int grey_eq_trn__(integer *nlyr, integer *nstr, integer *
	numu, integer *nphi, integer *nzup, integer *nzdn, integer *iref, 
	real *wnio, real *dtau, real *tatm, real *umu, real *phi, real *
	umu0nz, real *phi0nz, real *fbeam, real *surf_opt__, real *btemp, 
	real *ttemp, real *temis, logical *lsolar, logical *lplanck, logical *
	usrang, logical *lamber, real *upsflx, real *dnsflx, real *dirsflx, 
	real *sol_rad__, real *uptflx, real *dntflx, real *th_rad__)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    double asin(doublereal), exp(doublereal);

    /* Local variables */
    static real b[70];
    static integer k;
    static doublereal trans_flx__[14000]	/* was [70][200] */;
    static real fa, bs, bt, pi;
    static doublereal tr[14000]	/* was [70][200] */;
    static real bs0[1], bt0[1];
    static doublereal tr1[14000]	/* was [70][200] */, rad[14000]	/* 
	    was [70][200] */;
    static real bdr;
    static integer naz, nze;
    static doublereal uth[14000]	/* was [200][70] */, u0_s__[1120]	
	    /* was [16][70] */;
    static real dbdt[70];
    extern doublereal dref_(real *, real *, integer *);
    static real dphi, c1nu3;
    extern doublereal bdref_(real *, real *, real *, real *, integer *);
    static real g_phi__[16];
    static doublereal dtauk, trans[14000]	/* was [70][200] */;
    static real twopi;
    static doublereal trsol[70];
    static real btemp0[1], ttemp0[1];
    extern /* Subroutine */ int planck_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *);
    static doublereal tr_flx__[14000]	/* was [70][200] */, uu_flx__[17920]	
	    /* was [16][70][16] */;
    static real umu_ze__[200], gwt_ze__[200];
    extern /* Subroutine */ int qgausn_(integer *, real *, real *);
    static doublereal tr1_flx__[14000]	/* was [70][200] */, rad_flx__[14000]	
	    /* was [70][200] */;
    static real bdr_flx__, gwt_phi__[16];
    static doublereal uth_flx__[14000]	/* was [200][70] */;


/* cccccccccccccccccccccc  g r e y _ e q _ t r n  cccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine solves the equation of transfer for radiances   cc */
/* c    and fluxes in non-scattering atmospheres.  It assumes that the  cc */
/* c    planck function varies linearly with optical depth throughout   cc */
/* c    each level.                                                     cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       dtau - differential extintion optical depth in each layer    cc */
/* c       tatm - temperature of each atmospheric level                 cc */
/* c         wn - wavelnumber                                           cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c     usrang - output radiances at user angles? (logical: T/F)       cc */
/* c       numu - number of output zenith angles used in D/O code       cc */
/* c        umu - emission zenith angle cosines                         cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c        phi - emission azimuths read from input header              cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angle cosines read from header          cc */
/* c      fbeam - intensity of collimated flux at top of atmosphere     cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c                 note: this version assumes this                    cc */
/* c   surf_opt - surface albedo and other surface optical properties   cc */
/* c      btemp - surface temperature                                   cc */
/* c      ttemp - temperature of uppermost model layer (space)          cc */
/* c      temis - emissivity of uppermost model layer (space)           cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     dirsflx - downward direct solar flux at each model level       cc */
/* c      dnsflx - downward diffuse solar flux at each model level      cc */
/* c      upsflx - upward (diffuse) solar flux at each model level      cc */
/* c     sol_rad - solar radiance at each model level                   cc */
/* c      dntflx - downward diffuse thermal flux at each model level    cc */
/* c      uptflx - upward (diffuse) thermal flux at each model level    cc */
/* c      th_rad - thermal radiance at each model level                 cc */
/* c                                                                    cc */
/* cccccccccccccccccccccc  g r e y _ e q _ t r n  cccccccccccccccccccccccc */




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







/*       define local variables */





/* ****   spectrally-dependent output flux and radiance */

/*     .. */
    /* Parameter adjustments */
    th_rad__ -= 273;
    --dntflx;
    --uptflx;
    sol_rad__ -= 273;
    --dirsflx;
    --dnsflx;
    --upsflx;
    --surf_opt__;
    --phi;
    --umu;
    --tatm;
    --dtau;

    /* Function Body */

    if (pass1) {
	pass1 = FALSE_;
	pi = asin(1.f) * 2.f;
	twopi = pi * 2.f;

/* ****      find the gaussian angles needed to perform the integration */
/*          over azimuth */

	qgausn_(&c__8, g_phi__, gwt_phi__);

	for (k = 1; k <= 8; ++k) {
	    g_phi__[k + 7] = -g_phi__[k - 1];
	    gwt_phi__[k + 7] = gwt_phi__[k - 1];
/* L1001: */
	}

/* ****   find the gaussian angles needed to perform the flux integration */

	i__1 = *nstr / 2;
	qgausn_(&i__1, umu_ze__, gwt_ze__);

/* ****      restack the values into order used by disort */

	i__1 = *nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {
	    gwt_ze__[nze + *nstr / 2 - 1] = twopi * gwt_ze__[nze - 1] * 
		    umu_ze__[nze - 1];
	    umu_ze__[nze + *nstr / 2 - 1] = umu_ze__[nze - 1];
/* L1101: */
	}
	i__1 = *nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {
	    gwt_ze__[nze - 1] = gwt_ze__[*nstr - nze];
	    umu_ze__[nze - 1] = -umu_ze__[*nstr - nze];
/* L1121: */
	}

    }

/* ****   find the layer transmission for each downward zenith angle */

    i__1 = *nzdn;
    for (nze = 1; nze <= i__1; ++nze) {
	i__2 = *nlyr;
	for (k = 1; k <= i__2; ++k) {
	    dtauk = dtau[k] / umu[nze];
	    tr[k + nze * 70 - 71] = exp(dtauk);
	    tr1[k + nze * 70 - 71] = 1. - tr[k + nze * 70 - 71];
/* L1402: */
	}
    }

/* ****   find the layer transmission for each upward zenith angle */

    i__2 = *numu;
    for (nze = *nzup; nze <= i__2; ++nze) {
	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    dtauk = dtau[k] / umu[nze];
	    tr[k + nze * 70 - 71] = exp(-dtauk);
	    tr1[k + nze * 70 - 71] = 1. - tr[k + nze * 70 - 71];
/* L1422: */
	}
    }

    if (*usrang) {

/* ****    calculate the transmittances for use in finding solar fluxes */

/* ****    find the layer transmission for each downward zenith angle */

	i__1 = *nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {
	    i__2 = *nlyr;
	    for (k = 1; k <= i__2; ++k) {
		dtauk = dtau[k] / umu_ze__[nze - 1];
		tr_flx__[k + nze * 70 - 71] = exp(dtauk);
		tr1_flx__[k + nze * 70 - 71] = 1. - tr_flx__[k + nze * 70 - 
			71];
/* L1502: */
	    }
	}

/* ****     find the layer transmission for each upward zenith angle */

	i__2 = *nstr;
	for (nze = *nstr / 2 + 1; nze <= i__2; ++nze) {
	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {
		dtauk = dtau[k] / umu_ze__[nze - 1];
		tr_flx__[k + nze * 70 - 71] = exp(-dtauk);
		tr1_flx__[k + nze * 70 - 71] = 1. - tr_flx__[k + nze * 70 - 
			71];
/* L1522: */
	    }
	}

    }

/* ****  initialize radiances and fluxes */

    i__1 = *nlyr + 1;
    for (k = 1; k <= i__1; ++k) {
	dnsflx[k] = 0.f;
	dirsflx[k] = 0.f;
	upsflx[k] = 0.f;
	dntflx[k] = 0.f;
	uptflx[k] = 0.f;

	for (nze = 1; nze <= 16; ++nze) {
	    u0_s__[nze + (k << 4) - 17] = 0.f;
	    for (naz = 1; naz <= 16; ++naz) {
		sol_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
		th_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L1602: */
	    }
	}

	if (*usrang) {
	    i__2 = *nstr;
	    for (nze = 1; nze <= i__2; ++nze) {
		for (naz = 1; naz <= 16; ++naz) {
		    uu_flx__[nze + (k + naz * 70 << 4) - 1137] = 0.f;
/* L1622: */
		}
	    }

	}

/* L1641: */
    }

/* ****             s o l a r    r a d i a n c e s */

    if (*lsolar) {

/* ****      find the transmission along the solar zenith angle */

	trsol[0] = 1.f;
	i__1 = *nlyr + 1;
	for (k = 2; k <= i__1; ++k) {
	    dtauk = dtau[k - 1] / *umu0nz;
	    trsol[k - 1] = trsol[k - 2] * exp(-dtauk);
/* L2001: */
	}

/* ****     find the transmission between the top of the atmosphere and */
/*         each model level */

	i__1 = *nzdn;
	for (nze = 1; nze <= i__1; ++nze) {
	    trans[nze * 70 - 70] = 1.;
	    i__2 = *nlyr + 1;
	    for (k = 2; k <= i__2; ++k) {
		trans[k + nze * 70 - 71] = trans[k - 1 + nze * 70 - 71] * tr[
			k - 1 + nze * 70 - 71];
/* L2022: */
	    }
	}

/* ****      find the transmission along each stream from the top of */
/*          the atmosphere, to the surface, and back up. */

	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {
	    trans[*nlyr + 1 + nze * 70 - 71] = trsol[*nlyr];
	    for (k = *nlyr; k >= 1; --k) {
		trans[k + nze * 70 - 71] = trans[k + 1 + nze * 70 - 71] * tr[
			k + nze * 70 - 71];
/* L2042: */
	    }
	}

	if (*usrang) {

/* ****       find the transmission between the top of the atmosphere and */
/*           each model level */

	    i__2 = *nzdn;
	    for (nze = 1; nze <= i__2; ++nze) {
		trans_flx__[nze * 70 - 70] = 1.;
		i__1 = *nlyr + 1;
		for (k = 2; k <= i__1; ++k) {
		    trans_flx__[k + nze * 70 - 71] = trans_flx__[k - 1 + nze *
			     70 - 71] * tr_flx__[k - 1 + nze * 70 - 71];
/* L2122: */
		}
	    }

/* ****      find the transmission along each stream from the top of */
/*          the atmosphere, to the surface, and back up. */

	    i__1 = *nstr;
	    for (nze = *nzup; nze <= i__1; ++nze) {
		trans_flx__[*nlyr + 1 + nze * 70 - 71] = trsol[*nlyr];
		for (k = *nlyr; k >= 1; --k) {
		    trans_flx__[k + nze * 70 - 71] = trans_flx__[k + 1 + nze *
			     70 - 71] * tr_flx__[k + nze * 70 - 71];
/* L2142: */
		}
	    }
	} else {

/* ****        use the transmission values along the usual gaussian angles */

	    i__1 = *nstr;
	    for (nze = *nzup; nze <= i__1; ++nze) {
		i__2 = *nlyr + 1;
		for (k = 1; k <= i__2; ++k) {
		    trans_flx__[k + nze * 70 - 71] = trans[k + nze * 70 - 71];
/* L2202: */
		}
	    }
	}

/* ****     find the direct solar irradiances at each level */
/*         Note: the downward diffuse radiance should be zero */
/*         at each zenith angle and azimuth for the non-scattering case */

	i__2 = *nlyr + 1;
	for (k = 1; k <= i__2; ++k) {
	    dirsflx[k] = *umu0nz * *fbeam * trsol[k - 1];
/* L2501: */
	}

/* ****     find the upward scattered solar radiances at each level. */

	if (*lamber) {

/* ****        u p w a r d   r a d i a n c e s - l a m b e r t */

	    fa = *umu0nz * *fbeam * surf_opt__[1] / pi;
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		i__1 = *nlyr + 1;
		for (k = 1; k <= i__1; ++k) {
		    sol_rad__[nze + ((k << 4) + 1 << 4)] = fa * trans[k + nze 
			    * 70 - 71];
		    i__3 = *nphi;
		    for (naz = 2; naz <= i__3; ++naz) {
			sol_rad__[nze + (naz + (k << 4) << 4)] = sol_rad__[
				nze + ((k << 4) + 1 << 4)];
/* L2601: */
		    }
/* L2621: */
		}
/* L2641: */
	    }

/* ****        u p w a r d   s o l a r    f l u x - l a m b e r t */

/* ****      compute the azimuthally-averaged radiances */
/*          along the models computational streams */

/* ****      This version assumes the surface is lambertian */

	    i__2 = *nstr;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		i__1 = *nlyr + 1;
		for (k = 1; k <= i__1; ++k) {
		    u0_s__[nze + (k << 4) - 17] = *umu0nz * *fbeam * 
			    surf_opt__[1] * trans_flx__[k + nze * 70 - 71] / 
			    pi;
/* L2661: */
		}
/* L2681: */
	    }

	} else {

/* ****       u p w a r d    r a d i a n c e s - b r d f */

	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		i__1 = *nphi;
		for (naz = 1; naz <= i__1; ++naz) {
		    dphi = pi * (r__1 = phi[naz] - *phi0nz, dabs(r__1)) / 
			    180.f;
		    bdr = *umu0nz * *fbeam * bdref_(&umu[nze], umu0nz, &dphi, 
			    &surf_opt__[1], iref) / pi;
		    i__3 = *nlyr + 1;
		    for (k = 1; k <= i__3; ++k) {
			sol_rad__[nze + (naz + (k << 4) << 4)] = bdr * trans[
				k + nze * 70 - 71];
/* L2701: */
		    }
/* L2721: */
		}
/* L2741: */
	    }

/* ****        u p w a r d   s o l a r    f l u x  - b r d f */

/* ****      compute the azimuthally-averaged radiances */
/*          along the models computational streams */

/* ****      This version assumes the surface is lambertian */

	    bdr_flx__ = *umu0nz * *fbeam * dref_(umu0nz, &surf_opt__[1], iref)
		     / pi;

	    i__2 = *nstr;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		i__1 = *nlyr + 1;
		for (k = 1; k <= i__1; ++k) {
		    u0_s__[nze + (k << 4) - 17] = bdr * trans_flx__[k + nze * 
			    70 - 71];
/* L2761: */
		}
/* L2781: */
	    }

	}

/* ****      use gaussian quadrature to compute the upward fluxes */

	i__2 = *nstr;
	for (nze = *nstr / 2 + 1; nze <= i__2; ++nze) {
	    i__1 = *nlyr + 1;
	    for (k = 1; k <= i__1; ++k) {
		upsflx[k] += gwt_ze__[nze - 1] * u0_s__[nze + (k << 4) - 17];
/* L2902: */
	    }
	}

    }

/* ****    t h e r m a l   r a d i a n c e s   a n d    f l u x e s */

    if (*lplanck) {

/* ****   evaluate planck function at surface and at each model level */

/* ****           s u r f a c e    t h e r m a l    e m i s s i o n */

	btemp0[0] = *btemp;

	planck_(&c__1, &c__1, &c__1, &c__1, wnio, btemp0, bs0, &c1nu3);

	bs = bs0[0];

/* ****           t o p    o f   a t m o s p h e r e    e m i s s i o n */

	if (*ttemp > 0.f && *temis > 0.f) {
	    ttemp0[0] = *ttemp;

	    planck_(&c__1, &c__1, &c__1, &c__1, wnio, ttemp0, bt0, &c1nu3);

	    bt = bt0[0];
	} else {
	    bt = 0.f;
	}

/* ****        a t m o s p h e r i c     t h e r m a l    e m i s s i o n */

	planck_(&c__1, &c__1, &c__1, &c__1, wnio, &tatm[1], b, &c1nu3);

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {

	    planck_(&c__1, &c__1, &c__1, &c__1, wnio, &tatm[k], &b[k], &c1nu3)
		    ;

	    if (dtau[k] != 0.f) {
		dbdt[k - 1] = (b[k] - b[k - 1]) / dtau[k];
	    } else {
		dbdt[k - 1] = 0.f;
	    }

/* L3021: */
	}

/* ****     evaluate angle-dependent radiance contribution by each layer */


/* ****     downward radiance at base of this layer */

	i__1 = *nzdn;
	for (nze = 1; nze <= i__1; ++nze) {
	    i__2 = *nlyr;
	    for (k = 1; k <= i__2; ++k) {
		rad[k + nze * 70 - 71] = b[k - 1] * tr1[k + nze * 70 - 71] + 
			b[k] - b[k - 1] + umu[nze] * dbdt[k - 1] * tr1[k + 
			nze * 70 - 71];
/* L3102: */
	    }
	}

/* ****      upward radiance at top of this layer */

	i__2 = *nzup;
	for (nze = *nzdn + 1; nze <= i__2; ++nze) {
	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {
		rad[k + nze * 70 - 71] = b[k - 1] * tr1[k + nze * 70 - 71] - (
			b[k] - b[k - 1]) * tr[k + nze * 70 - 71] + umu[nze] * 
			dbdt[k - 1] * tr1[k + nze * 70 - 71];
/* L3122: */
	    }
	}

	if (*usrang) {

/* ****     find downward and upward radiances on the gaussian */
/*         computational grid for use in the flx calculations */

/* ****     downward radiance at base of this layer */

	    i__1 = *nstr / 2;
	    for (nze = 1; nze <= i__1; ++nze) {
		i__2 = *nlyr;
		for (k = 1; k <= i__2; ++k) {
		    rad_flx__[k + nze * 70 - 71] = b[k - 1] * tr1_flx__[k + 
			    nze * 70 - 71] + b[k] - b[k - 1] + umu_ze__[nze - 
			    1] * dbdt[k - 1] * tr1_flx__[k + nze * 70 - 71];
/* L3142: */
		}
	    }

/* ****        upward radiance at top of this layer */

	    i__2 = *nstr;
	    for (nze = *nstr / 2 + 1; nze <= i__2; ++nze) {
		i__1 = *nlyr;
		for (k = 1; k <= i__1; ++k) {
		    rad_flx__[k + nze * 70 - 71] = b[k - 1] * tr1_flx__[k + 
			    nze * 70 - 71] - (b[k] - b[k - 1]) * tr_flx__[k + 
			    nze * 70 - 71] + umu_ze__[nze - 1] * dbdt[k - 1] *
			     tr1_flx__[k + nze * 70 - 71];
/* L3162: */
		}
	    }
	} else {
	    i__1 = *nstr;
	    for (nze = 1; nze <= i__1; ++nze) {
		i__2 = *nlyr;
		for (k = 1; k <= i__2; ++k) {
		    rad_flx__[k + nze * 70 - 71] = rad[k + nze * 70 - 71];
/* L3182: */
		}
	    }
	}

/* ****     d o w n w a r d     t h e r m a l     r a d i a n c e */

/* ****     find the downward diffuse radiances at each level and angle */
/*         by integrating downward from the top of the atmosphere. */

	i__2 = *nzdn;
	for (nze = 1; nze <= i__2; ++nze) {

/* ****         set radiance at top of atmosphere to top BC */

	    uth[nze - 1] = *temis * bt;

/* ****             specify downward radiances at all other levels */

	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {
		uth[nze + (k + 1) * 200 - 201] = rad[k + nze * 70 - 71] + tr[
			k + nze * 70 - 71] * uth[nze + k * 200 - 201];
/* L3201: */
	    }

/* L3221: */
	}

/* ****      find the downward irradiance (flux) at each level */

	if (*usrang) {

/* ****       find downward radiances on the gaussian computational grid */

	    i__2 = *nstr / 2;
	    for (nze = 1; nze <= i__2; ++nze) {

/* ****           set radiance at top of atmosphere to top BC */

		uth_flx__[nze - 1] = *temis * bt;

/* ****               specify downward radiances at all other levels */

		i__1 = *nlyr;
		for (k = 1; k <= i__1; ++k) {
		    uth_flx__[nze + (k + 1) * 200 - 201] = rad_flx__[k + nze *
			     70 - 71] + tr_flx__[k + nze * 70 - 71] * 
			    uth_flx__[nze + k * 200 - 201];
/* L3301: */
		}

/* L3321: */
	    }
	} else {

	    i__2 = *nstr / 2;
	    for (nze = 1; nze <= i__2; ++nze) {
		i__1 = *nlyr;
		for (k = 1; k <= i__1; ++k) {
		    uth_flx__[nze + (k + 1) * 200 - 201] = uth[nze + (k + 1) *
			     200 - 201];
/* L3342: */
		}
	    }

	}

/* ****      use gaussian quadrature to integrate over the upper */
/*          hemisphere (note, this explicitly assumes that the */
/*          downward thermal radiances are azimuthally uniform */

	i__1 = *nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {
	    i__2 = *nlyr + 1;
	    for (k = 1; k <= i__2; ++k) {
		dntflx[k] += gwt_ze__[nze - 1] * uth_flx__[nze + k * 200 - 
			201];
/* L3422: */
	    }
	}

/* ****     u p w a r d     t h e r m a l    r a d i a n c e */

/* ****     find the upward diffuse radiances at each level and angle */
/*         by integrating upward from the surface. */

	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {

/* ****         specify the surface boundary condition */
/*             note: This version still assumes a lamber surface */
/*             at thermal wavelengths */

/* ****           assume that the the surface is lambertian */

	    uth[nze + (*nlyr + 1) * 200 - 201] = (1.f - surf_opt__[1]) * bs + 
		    surf_opt__[1] * dntflx[*nlyr + 1] / pi;

/* ****         find radiance radiances at all other levels */

	    for (k = *nlyr; k >= 1; --k) {
		uth[nze + k * 200 - 201] = rad[k + nze * 70 - 71] + tr[k + 
			nze * 70 - 71] * uth[nze + (k + 1) * 200 - 201];
/* L3601: */
	    }

/* L3621: */
	}

/* ****      find the upward irradiance (flux) at each level */

/* ****      use gaussian quadrature to integrate over the lower hemisphere */

	i__2 = *nstr;
	for (nze = *nzup; nze <= i__2; ++nze) {
	    i__1 = *nlyr + 1;
	    for (k = 1; k <= i__1; ++k) {
		uptflx[k] += gwt_ze__[nze - 1] * uth_flx__[nze + k * 200 - 
			201];
/* L3822: */
	    }
	}

    }

    return 0;
} /* grey_eq_trn__ */

