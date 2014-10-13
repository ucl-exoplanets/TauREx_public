/* ray_tau.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ray_tau__(integer *ne, integer *nlyr, integer *ncomp, 
	integer *icomp, real *volmix, real *p, real *grav, real *wgtatm, 
	integer *nzup, integer *nstr, integer *iflext, integer *nsiext, 
	integer *istate, real *pd_frac__, doublereal *wnext, doublereal *wn, 
	real *umu, real *dtauex, real *dtausc, real *phmom, real *tauray, 
	real *p_ray__, real *p_ray_0__)
{
    /* Initialized data */

    static real a0 = 6.02297e26f;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static real dtau_tst__;
    static integer k, l;
    static real wl;
    static integer nze, nlay, ktau1, nt_pd__;
    static real dtray;
    extern doublereal raylei_(real *, integer *, integer *, real *);
    static real sigray, ptauray;


/* cccccccccccccccccccccccccc  r a y _ t a u  cccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine finds the Rayleigh scattering cross sections    cc */
/* c    at each model level, and computes the effective monochromatic   cc */
/* c    Rayleigh scattering optical depth for each atmospheric layer.   cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr - number of layers in the model atmosphere         .    cc */
/* c      ncomp - number of rayleigh-scattering constituentes           cc */
/* c      icomp - index of each rayleigh scattering constituent         cc */
/* c         a0 - avogodro's number.                                    cc */
/* c     volmix - volume mixing ratio of each rayleigh scatterer        cc */
/* c       nzup - index of first upward radiance stream.                cc */
/* c          p - pressure at each model level (Pascals).               cc */
/* c       grav - gravitational acceleration at each level (m/s^2).     cc */
/* c       nstr - number of computational zenith angles (streams).      cc */
/* c        umu - cosine of each upward radiance stream.                cc */
/* c         wn - current wavenumber (cm**-1).                          cc */
/* c     iflext - index for extinction sources                          cc */
/* c     nsiext - index counter for extinction sources                  cc */
/* c      wnext - wavenumber of the extincition source                  cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     tauray - column-integrated rayleigh optical depth.             cc */
/* c     dtauex - monochromatic extinction optical depth in each layer. cc */
/* c     dtausc - monochromatic scattering optical depth in each layer. cc */
/* c      phmom - scattering phase function moments in each layer.      cc */
/* c      p_ray - pressure level of rayliegh tau=1 for each upward      cc */
/* c             stream (bars).                                         cc */
/* c    p_ray_0 - pressure of Rayleigh tau=1 for a vertical stream.     cc */
/* c     iflext - index for extinction sources                          cc */
/* c     nsiext - index counter for extinction sources                  cc */
/* c      wnext - wavenumber of the extincition source                  cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  r a y _ t a u  cccccccccccccccccccccccccccc */




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



/* *****   wavenumber grid variables */




/* ****   partial derivative fractional change */


/* ****  rayleigh scattering optical properties */


/* ****    define avagadro's number (mks units: molecules/kmole) */

    /* Parameter adjustments */
    --volmix;
    --icomp;
    --p;
    --grav;
    --iflext;
    nsiext -= 3;
    --istate;
    --pd_frac__;
    wnext -= 3;
    --umu;
    dtauex -= 71;
    --dtausc;
    phmom -= 201;
    --tauray;
    --p_ray__;

    /* Function Body */

/* *****   determine whether temperature is a variable part of */
/*        the state vector.  if it is, extinction optical depth must be */
/*        computed for the background and perterbed temperature profile */

    if (istate[2] == 0) {
	nt_pd__ = 1;
    } else {
	nt_pd__ = 2;
    }

/* ****   define the local wavelength (microns) */

    wl = (real) (1e4 / *wn);

/* ****       find the rayleigh scattering cross section at wl */

    sigray = raylei_(&wl, ncomp, &icomp[1], &volmix[1]);
    tauray[1] = 0.f;

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

/* ****         find rayleigh scattering optical depth */

	dtray = a0 * sigray * (p[k + 1] - p[k]) / (*wgtatm * .5f * (grav[k] + 
		grav[k + 1]));
	tauray[k + 1] = tauray[k] + dtray;
	dtausc[k] += dtray;
	i__2 = nt_pd__;
	for (l = 1; l <= i__2; ++l) {
	    dtauex[k + l * 70] += dtray;
/* L4001: */
	}

/* ****        define rayleigh scattering phase function moments */

	phmom[k * 201] += dtray;
	phmom[k * 201 + 1] += dtray * .1f;
/* L4021: */
    }

/* ****    if surface pressure is a variable part of the state vector, */
/*        find the optical depth for a layer below the base of the */
/*        atmosphere that is 1% thicker than the nominal layer */

    if (istate[1] != 0) {
	nlay = *nlyr + 1;
	k = nlay;
	dtray = a0 * sigray * ((pd_frac__[1] + 1.f) * p[k] - p[k - 1]) / (*
		wgtatm * .5f * (grav[k] + grav[k - 1]));
	dtausc[k] += dtray;
	dtauex[k + 70] += dtray;
/*      write(*,'(1a,i5,3(1pe14.6))') 'ray_tau: nlay,pd_frac,dtauex:', */
/*     -      nlay,pd_frac(1),dtauex(nlyr,1),dtauex(nlay,1) */

/* ****        define rayleigh scattering phase function moments */

	phmom[k * 201] += dtray;
	phmom[k * 201 + 1] += dtray * .1f;
    }

/* ****       find pressure of rayleigh scattering optical depth unity: */

    ktau1 = *nlyr + 1;
    ptauray = 0.f;
    i__1 = *nstr;
    for (nze = *nzup; nze <= i__1; ++nze) {
	i__2 = *nlyr;
	for (k = 1; k <= i__2; ++k) {
	    if (ptauray == 0.f && tauray[k + 1] / umu[nze] >= 1.f) {
		ktau1 = k + 1;
		ptauray = p[k];
	    }
/* L4101: */
	}
	dtau_tst__ = (tauray[ktau1] - tauray[ktau1 - 1]) / umu[nze];
	if (dtau_tst__ != 0.f) {
	    p_ray__[nze] = (p[ktau1 - 1] + (1.f - tauray[ktau1 - 1] / umu[nze]
		    ) * (p[ktau1] - p[ktau1 - 1]) / dtau_tst__) * 1e-5f;
	} else {
	    p_ray__[nze] = p[ktau1] * 1e-5f;
	}
/* L4121: */
    }

/* ****       find pressure of rayleigh scattering optical depth unity */
/*           for a vertical path: */

    ktau1 = *nlyr + 1;
    ptauray = 0.f;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	if (ptauray == 0.f && tauray[k + 1] >= 1.f) {
	    ktau1 = k + 1;
	    ptauray = p[k];
	}
/* L4141: */
    }
    dtau_tst__ = tauray[ktau1] - tauray[ktau1 - 1];
    if (dtau_tst__ != 0.f) {
	*p_ray_0__ = (p[ktau1 - 1] + (1.f - tauray[ktau1 - 1]) * (p[ktau1] - 
		p[ktau1 - 1]) / dtau_tst__) * 1e-5f;
    } else {
	*p_ray_0__ = p[ktau1] * 1e-5f;
    }

/* ****       update extinction counter and rayleigh flag */

    ++(*ne);
    iflext[*ne] = 0;

/* *****       find the next required rayleigh scattering wavenumber */

    nsiext[(*ne << 1) + 2] = 2;
    if (tauray[ktau1] < .01f) {
	wnext[(*ne << 1) + 2] = *wn * 1.05;
    } else {
	if (tauray[ktau1] < .1f) {
	    wnext[(*ne << 1) + 2] = *wn * 1.01;
	} else {
	    wnext[(*ne << 1) + 2] = *wn * 1.001;
	}
    }

    return 0;
} /* ray_tau__ */

