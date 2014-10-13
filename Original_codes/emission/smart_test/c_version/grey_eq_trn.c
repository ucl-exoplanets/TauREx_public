/* grey_eq_trn.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int grey_eq_trn__(logical *lsolar, logical *lplanck, logical 
	*usrang, logical *lamber, integer *nlev, integer *nz0, integer *nstr, 
	integer *numu, integer *nphi, integer *nzup, integer *nzdn, integer *
	iref, doublereal *wn_io__, real *dtau, real *copi0, real *tatm, real *
	umu, real *umu_f__, real *gwt_f__, real *phi, real *umu0nz, real *
	phi0nz, real *fbeam, real *alb0, real *btemp, real *ttemp, real *
	temis, doublereal *trndir, doublereal *trnflx, doublereal *absflx, 
	doublereal *trn_cmu__, doublereal *abs_cmu__, doublereal *trnrad, 
	doublereal *absrad, real *upsflx, real *dnsflx, real *dirsflx, real *
	sol_rad__, real *uptflx, real *dntflx, real *th_rad__, doublereal *
	dn_s_src__, doublereal *up_s_src__, doublereal *dn_t_src__, 
	doublereal *up_t_src__)
{
    /* Initialized data */

    static integer ipass1 = 0;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double acos(doublereal);

    /* Local variables */
    static real b[70];
    static integer k;
    static doublereal bb_flx_dn__[70], bb_flx_up__[70], trans_flx__[1120]	
	    /* was [70][16] */, bb[70];
    static real fa, bs;
    static integer kk;
    static real bt, pi;
    static integer naz, nze;
    static real dbdt[70], wnio;
    static integer nlyr;
    static doublereal thrad[1120]	/* was [16][70] */, trans[1120]	/* 
	    was [70][16] */;
    static real twopi;
    static doublereal trsol[70];
    static real btemp0[1], ttemp0[1];
    extern /* Subroutine */ int planck_(integer *, integer *, integer *, 
	    integer *, real *, real *, doublereal *);
    static real rad_flx__[1120]	/* was [70][16] */;
    static doublereal dth_flx__[1120]	/* was [16][70] */, uth_flx__[1120]	
	    /* was [16][70] */;


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
/* c       alb0 - surface albedo and other surface optical properties   cc */
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



/* ****   input integer variables */



/* ****   input real variables */



/*       define local real variables */



/* ****    layer transmission values for simplified adding method */





/* ****   spectrally-dependent output flux and radiance */


/* ****   interpolated solar and thermal source functions */



/* ****   define the speed of light (mks) */

    /* Parameter adjustments */
    --up_t_src__;
    --dn_t_src__;
    --up_s_src__;
    --dn_s_src__;
    th_rad__ -= 273;
    --dntflx;
    --uptflx;
    sol_rad__ -= 273;
    --dirsflx;
    --dnsflx;
    --upsflx;
    absrad -= 17;
    trnrad -= 17;
    abs_cmu__ -= 17;
    trn_cmu__ -= 17;
    --absflx;
    --trnflx;
    --trndir;
    --phi;
    --gwt_f__;
    --umu_f__;
    --umu;
    --tatm;
    --copi0;
    --dtau;

    /* Function Body */

    nlyr = *nlev - 1;
    if (ipass1 == 0) {

	ipass1 = 1;
	pi = acos(-1.f);
	twopi = pi * 2.f;

    }

/* ****         find the flux transmittance and absorptance by */
/*             integrating transmittance over zenith angle using */
/*             gaussian quadrature */

    i__1 = nlyr;
    for (k = 1; k <= i__1; ++k) {
	trnflx[k] = 0.f;
	absflx[k] = 0.f;
	i__2 = *nstr / 2;
	for (nze = 1; nze <= i__2; ++nze) {
	    trnflx[k] += gwt_f__[nze] * 2.f * umu_f__[nze] * trn_cmu__[nze + (
		    k << 4)];
	    absflx[k] += gwt_f__[nze] * 2.f * umu_f__[nze] * abs_cmu__[nze + (
		    k << 4)];
/* L1201: */
	}
/* L1221: */
    }

/* *********************************************************************** */

/* ****             s o l a r    r a d i a n c e s */

/* *********************************************************************** */

    if (*lsolar) {

/* ****  initialize radiances and fluxes */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    dn_s_src__[k] = 0.f;
	    up_s_src__[k] = 0.f;
	    dnsflx[k] = 0.f;
	    dirsflx[k] = 0.f;
	    upsflx[k] = 0.f;

	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		rad_flx__[k + nze * 70 - 71] = 0.f;
		i__3 = *nphi;
		for (naz = 1; naz <= i__3; ++naz) {
		    sol_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L1801: */
		}
/* L1821: */
	    }
/* L1841: */
	}

/* ****      find the transmission along the solar zenith angle */

	trsol[0] = 1.;
	i__1 = *nlev;
	for (k = 2; k <= i__1; ++k) {
	    trsol[k - 1] = trsol[k - 2] * trndir[k];
/* L2001: */
	}

/* ****      find the transmission along each stream from the top of */
/*          the atmosphere, to the surface, and back up. */

	i__1 = *numu;
	for (nze = *nzup; nze <= i__1; ++nze) {
	    trans[*nlev + nze * 70 - 71] = trsol[*nlev - 1];
	    for (k = nlyr; k >= 1; --k) {
		trans[k + nze * 70 - 71] = trans[k + 1 + nze * 70 - 71] * 
			trnrad[nze + (k << 4)];
/* L2021: */
	    }
/* L2041: */
	}

/* ****      find the transmission along each upward stream from */
/*          surface to the top of the atmosphere */

	i__1 = *nstr;
	for (nze = *nstr / 2 + 1; nze <= i__1; ++nze) {
	    trans_flx__[*nlev + nze * 70 - 71] = trsol[*nlev - 1];
	    for (k = nlyr; k >= 1; --k) {
		trans_flx__[k + nze * 70 - 71] = trans_flx__[k + 1 + nze * 70 
			- 71] * trn_cmu__[nze + (k << 4)];
/* L2061: */
	    }
/* L2081: */
	}

/* ****     find the direct downward solar irradiances at each level */
/*         Note: the downward diffuse radiance should be zero */
/*         at each zenith angle and azimuth for the non-scattering case */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    dirsflx[k] = *umu0nz * *fbeam * (real) trsol[k - 1];
/* L2101: */
	}

/* ****     find the upward scattered solar radiances at each level. */

/* ****        u p w a r d   r a d i a n c e s - l a m b e r t */

	fa = *umu0nz * *fbeam * *alb0 / pi;
	i__1 = *numu;
	for (nze = *nzup; nze <= i__1; ++nze) {
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		sol_rad__[nze + ((k << 4) + 1 << 4)] = fa * (real) trans[k + 
			nze * 70 - 71];
		i__3 = *nphi;
		for (naz = 2; naz <= i__3; ++naz) {
		    sol_rad__[nze + (naz + (k << 4) << 4)] = sol_rad__[nze + (
			    (k << 4) + 1 << 4)];
/* L2201: */
		}
/* L2221: */
	    }
/* L2241: */
	}

/* ****        u p w a r d   s o l a r    f l u x - l a m b e r t */

/* ****      compute the azimuthally-averaged radiances */
/*          along the models computational streams */

	i__1 = *nstr;
	for (nze = *nstr / 2 + 1; nze <= i__1; ++nze) {
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		rad_flx__[k + nze * 70 - 71] = fa * (real) trans_flx__[k + 
			nze * 70 - 71];
/* L2261: */
	    }
/* L2281: */
	}

/* ****      use gaussian quadrature to compute the upward fluxes */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nstr;
	    for (nze = *nstr / 2 + 1; nze <= i__2; ++nze) {
		upsflx[k] += twopi * umu_f__[nze] * gwt_f__[nze] * rad_flx__[
			k + nze * 70 - 71];
/* L2301: */
	    }
/* L2321: */
	}

/* ****      initialize the effective 'source' for each layer */

/*         note: the downward diffuse solar fluxes are zero by definition */

/* ****     initialize the upward fluxes at the surface */

	up_s_src__[*nlev] = upsflx[*nlev];

	i__1 = *nlev;
	for (k = 2; k <= i__1; ++k) {

/* ****         find upward flux source at top of each layer */

	    kk = *nlev - k + 1;
	    up_s_src__[kk] = upsflx[kk] - upsflx[kk + 1] * trnflx[kk];

/* L2401: */
	}

    }
/* *********************************************************************** */

/* ****    t h e r m a l   r a d i a n c e s   a n d    f l u x e s */

/* *********************************************************************** */
    if (*lplanck && *nz0 == 1) {

/* ****  initialize radiances and fluxes */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    dn_t_src__[k] = 0.f;
	    up_t_src__[k] = 0.f;
	    dntflx[k] = 0.f;
	    uptflx[k] = 0.f;

	    i__2 = *numu;
	    for (nze = 1; nze <= i__2; ++nze) {
		i__3 = *nphi;
		for (naz = 1; naz <= i__3; ++naz) {
		    th_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L3001: */
		}
/* L3021: */
	    }
/* L3041: */
	}

/* ****   evaluate planck function at surface and at each model level */

/* ****           s u r f a c e    t h e r m a l    e m i s s i o n */

	wnio = (real) (*wn_io__);
	btemp0[0] = *btemp;

	planck_(&c__1, &c__1, &c__1, &c__1, &wnio, btemp0, bb);

	bs = (real) bb[0];

/* ****           t o p    o f   a t m o s p h e r e    e m i s s i o n */

	if (*ttemp > 0.f && *temis > 0.f) {
	    ttemp0[0] = *ttemp;

	    planck_(&c__1, &c__1, &c__1, &c__1, &wnio, ttemp0, bb);

	    bt = (real) bb[0];
	} else {
	    bt = 0.f;
	}

/* ****        a t m o s p h e r i c     t h e r m a l    e m i s s i o n */

	planck_(&c__1, nlev, &c__1, &c__1, &wnio, &tatm[1], bb);

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    b[k - 1] = (real) bb[k - 1];
/* L3101: */
	}

	i__1 = nlyr;
	for (k = 1; k <= i__1; ++k) {

	    if (copi0[k] * dtau[k] != 0.f) {
		dbdt[k - 1] = (b[k] - b[k - 1]) / (copi0[k] * dtau[k]);
	    } else {
		dbdt[k - 1] = 0.f;
	    }

/* L3121: */
	}

/* ****     d o w n w a r d     t h e r m a l     r a d i a n c e */

/* ****     find the downward diffuse radiances at each level and angle */
/*         by integrating downward from the top of the atmosphere. */

	i__1 = *nzdn;
	for (nze = 1; nze <= i__1; ++nze) {

/* ****         set radiance at top of atmosphere to top BC */

	    thrad[nze - 1] = *temis * bt;

/* ****             specify downward radiances at all other levels */

	    i__2 = nlyr;
	    for (k = 1; k <= i__2; ++k) {
		thrad[nze + (k + 1 << 4) - 17] = (real) (b[k] - b[k - 1] + (b[
			k - 1] + umu[nze] * dbdt[k - 1]) * absrad[nze + (k << 
			4)] + trnrad[nze + (k << 4)] * thrad[nze + (k << 4) - 
			17]);
/* Danie Liang */
		if (dtau[k] == 0.f) {
		    thrad[nze + (k + 1 << 4) - 17] = thrad[nze + (k << 4) - 
			    17];
		}
/* L3201: */
	    }

/* L3221: */
	}

/* ****       find downward radiances on the gaussian computational grid */

	i__1 = *nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {

/* ****         set radiance at top of atmosphere to top BC */

	    dth_flx__[nze - 1] = *temis * bt;

/* ****             specify downward radiances at all other levels */

	    i__2 = nlyr;
	    for (k = 1; k <= i__2; ++k) {
		dth_flx__[nze + (k + 1 << 4) - 17] = b[k] - b[k - 1] + (b[k - 
			1] + umu_f__[nze] * dbdt[k - 1]) * abs_cmu__[nze + (k 
			<< 4)] + trn_cmu__[nze + (k << 4)] * dth_flx__[nze + (
			k << 4) - 17];
/* Danie Liang */
		if (dtau[k] == 0.f) {
		    dth_flx__[nze + (k + 1 << 4) - 17] = dth_flx__[nze + (k <<
			     4) - 17];
		}
		if (dth_flx__[nze + (k + 1 << 4) - 17] < 0.) {
		    dth_flx__[nze + (k + 1 << 4) - 17] = 0.f;
		}
/* L3301: */
	    }

/* L3321: */
	}

/* ****      use gaussian quadrature to integrate over the upper */
/*          hemisphere (note, this explicitly assumes that the */
/*          downward thermal radiances are azimuthally uniform */

	i__1 = *nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
/* Danie Liang, change size of -twopi to +twopi */
		dntflx[k] += twopi * umu_f__[nze] * gwt_f__[nze] * (real) 
			dth_flx__[nze + (k << 4) - 17];
/* L3341: */
	    }
/* L3351: */
	}

/* ****     u p w a r d     t h e r m a l    r a d i a n c e */

/* ****     find the upward diffuse radiances at each level and angle */
/*         by integrating upward from the surface. */

/* ****        use a lambertian surface */

	i__1 = *numu;
	for (nze = *nzup; nze <= i__1; ++nze) {

/* ****           specify the surface boundary condition */
/*               note: This version still assumes a lamber surface */
/*               at thermal wavelengths */

/* ****             assume that the the surface is lambertian */

	    thrad[nze + (*nlev << 4) - 17] = (1.f - *alb0) * bs + *alb0 * 
		    dntflx[*nlev] / pi;

/* ****           find radiance radiances at all other levels */

	    for (k = nlyr; k >= 1; --k) {
/* Danie Liang */
		thrad[nze + (k << 4) - 17] = b[k - 1] * absrad[nze + (k << 4)]
			 - (b[k] - b[k - 1]) * trnrad[nze + (k << 4)] + umu[
			nze] * dbdt[k - 1] * absrad[nze + (k << 4)] + trnrad[
			nze + (k << 4)] * thrad[nze + (k + 1 << 4) - 17];
/*                thrad(nze,k) = b(k+1) - b(k) + (b(k) - */
/*     -                           umu(nze)*dbdt(k))*absrad(nze,k) + */
/*     -                           trnrad(nze,k)*thrad(nze,k+1) */
/* L3401: */
	    }

/* L3421: */
	}

/* *****    load output arrays */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = 1; nze <= i__3; ++nze) {
		    th_rad__[nze + (naz + (k << 4) << 4)] = (real) thrad[nze 
			    + (k << 4) - 17];
/* L3601: */
		}
/* L3621: */
	    }
/* L3641: */
	}

/* ****     u p w a r d     t h e r m a l    f l u x */

/* ****       find downward radiances on gaussian computational grid */

	i__1 = *nstr;
	for (nze = *nstr / 2 + 1; nze <= i__1; ++nze) {

/* ****         set radiance at top of atmosphere to top BC */

	    uth_flx__[nze + (*nlev << 4) - 17] = (1.f - *alb0) * bs + *alb0 * 
		    dntflx[*nlev] / pi;

/* ****             specify downward radiances at all other levels */

	    for (k = nlyr; k >= 1; --k) {
/*                uth_flx(nze,k) = b(k+1) - b(k) + (b(k) - */
/*     -                           umu_f(nze)*dbdt(k))*abs_cmu(nze,k) + */
/*     -                           trn_cmu(nze,k)*uth_flx(nze,k+1) */
/* Danie Liang */
		uth_flx__[nze + (k << 4) - 17] = b[k - 1] * absrad[nze + (k <<
			 4)] - (b[k] - b[k - 1]) * trn_cmu__[nze + (k << 4)] 
			+ umu_f__[nze] * dbdt[k - 1] * abs_cmu__[nze + (k << 
			4)] + trn_cmu__[nze + (k << 4)] * uth_flx__[nze + (k 
			+ 1 << 4) - 17];
/* L3701: */
	    }

/* L3721: */
	}

/* ****      use gaussian quadrature to integrate over the upper */
/*          hemisphere (note, this explicitly assumes that the */
/*          downward thermal radiances are azimuthally uniform */

	i__1 = *nstr;
	for (nze = *nstr / 2 + 1; nze <= i__1; ++nze) {
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		uptflx[k] += twopi * umu_f__[nze] * gwt_f__[nze] * (real) 
			uth_flx__[nze + (k << 4) - 17];
/* L3741: */
	    }
/* L3751: */
	}

/* ****    t h e r m a l   s o u r c e    f u n c t i o n s */

/* ****           define the black body source functions for flux and */
/*               radiance, assuming a linear-in-tau formulation. */

	bb_flx_dn__[0] = 0.f;
	bb_flx_up__[*nlev - 1] = pi * bs * (1. - *alb0);
	i__1 = nlyr;
	for (k = 1; k <= i__1; ++k) {
	    bb_flx_up__[k - 1] = pi * .5f * (b[k - 1] + b[k]) * absflx[k];
	    bb_flx_dn__[k] = pi * .5f * (b[k - 1] + b[k]) * absflx[k];
/* L3801: */
	}

/* ****     compute the effective 'source function' for flux in each */
/*         layer for use in climate modeling */

	dn_t_src__[1] = 0.f;
	up_t_src__[*nlev] = uptflx[*nlev] - bb_flx_up__[*nlev - 1];

	i__1 = *nlev;
	for (k = 2; k <= i__1; ++k) {

/* ****         find downward source at bottom of each layer */

	    dn_t_src__[k] = dntflx[k] - dntflx[k - 1] * trnflx[k - 1] - 
		    bb_flx_dn__[k - 1];

/* ****         find upward source at top of each layer */

	    kk = *nlev - k + 1;
	    up_t_src__[kk] = uptflx[kk] - bb_flx_up__[kk - 1] - uptflx[kk + 1]
		     * trnflx[kk];

/* L3821: */
	}

    }

    return 0;
} /* grey_eq_trn__ */

