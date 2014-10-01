/* gasptau.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int gasptau_(integer *ne, integer *iuabc, integer *ng0, 
	integer *ngt0, integer *igpab, integer *nlev, integer *nsiext, 
	integer *nlay, integer *istate, integer *nt_pd__, integer *ntaug, 
	doublereal *wn, doublereal *wnmin, doublereal *wnmax, doublereal *
	wnext, real *p, real *t, real *alt, real *rmix, real *d_pres__, real *
	d_temp__, real *d_rmix__, real *dtaug, real *dtaugdv, doublereal *
	wneof, integer *io_end__, integer *io_err__)
{
    /* Initialized data */

    static real rgas = 8314.f;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Local variables */
    extern /* Subroutine */ int xyinterp_(real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer k, l, m, n, ng, ni, mm, ns;
    static real tt, tt0[1], pab[70], dvi;
    static integer ngt;
    static real rho[210]	/* was [70][3] */, pab0[700000]	/* was [70][
	    10000] */, pab1[70];
    static integer npab[15]	/* was [3][5] */, nlyr, npab0;
    static real tpab0[70];
    static integer ngpab, ntpab;
    static doublereal wnpab[150000]	/* was [10000][3][5] */;
    static integer nsipr[15]	/* was [3][5] */, nxmax, nymax, iuabc0;
    static doublereal wnpab0[10000];
    static real convrt;
    extern /* Subroutine */ int readpab_(integer *, doublereal *, doublereal *
	    , integer *, integer *, real *, doublereal *, real *, integer *, 
	    integer *, doublereal *);
    static real dtaupab[63000000]	/* was [70][10000][3][5][2][3] */;


/* ccccccccccccccccccccccccc  g a s p t a u  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads gas p-induced absorption coefficients &   cc */
/* c    interpolates them to the appropriate pressure, and temperature  cc */
/* c    grids.                                                          cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c         ne - index of extinction source                            cc */
/* c        ng0 - gas index (1 - ngases)                                cc */
/* c      nt_pd - number of temperature values where absorption         cc */
/* c              coefficients are needed (2 are needed if temperature  cc */
/* c      wnmin - minimum wavenumber (cm**-1)                           cc */
/* c      wnmax - maximum wavenumber (cm**-1)                           cc */
/* c     d_pres - perturbed surface pressure used when pressure is a    cc */
/* c              variable part of the state vector                     cc */
/* c     d_temp - fractional temperature perturbation used when         cc */
/* c              temperature is a variable componen of the state       cc */
/* c              vector.                                               cc */
/* c     d_rmix - perturbed gas mixing ratio used when rmix is a        cc */
/* c              variable part of the state vector                     cc */
/* c      ntaug - number of optical depths needed for this gas:         cc */
/* c              =1 if this gas is not a variable part of state vector cc */
/* c              =3 if this gas is a variable part of state vector     cc */
/* c              (nominal, rmix perturbed at top of layer, rmix        cc */
/* c               perturbed at bottom of layer).                       cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     io_end - end-of file flag                                      cc */
/* c     io_err - flag specifying error in file                         cc */
/* c      wneof - wavenumber at end of file                             cc */
/* c      wnpab - wavenumber of each pressure-induced optical depth     cc */
/* c        pab - pressure-induced gas optical depth at each input wn   cc */
/* c       npab - number of wavenumbers were pab's are specified        cc */
/* c      dtaug - gas optical depth in each layer and input wavelength  cc */
/* c    dtaugdv - wavenumber gradient of each gas optical depth         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  g a s p t a u  ccccccccccccccccccccccccccccc */




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






/* ****   atmospheric structure variables */


/* ****   gas absorption coefficients */


/* ****   partial derivative fractional change */


    /* Parameter adjustments */
    --io_err__;
    --io_end__;
    --wneof;
    dtaugdv -= 23801;
    dtaug -= 47671;
    d_rmix__ -= 3;
    --d_temp__;
    --d_pres__;
    rmix -= 71;
    --alt;
    --t;
    --p;
    wnext -= 3;
    --istate;
    nsiext -= 3;
    igpab -= 4;
    iuabc -= 4;

    /* Function Body */

    nlyr = *nlev - 1;
    ng = *ng0;
    ngt = *ngt0;
    ngpab = igpab[ngt + ng * 3];

    if (nsiext[(*ne << 1) + 2] == 0) {

/* ****     no pressure induced absorption coefficents have been */
/*         read in yet.  define conversion from cgs to cm-amagat */
/*         (1 amagat = density at stp = 0.044616 kmoles/m**3) */

	convrt = 22.413483951945491f;

/* ****     find the density (amagats) and absorber amount (cm-amagat**2) */

/*            write(*,*) '08/18, output rmix to see if it is mmr' */

	i__1 = *nt_pd__;
	for (l = 1; l <= i__1; ++l) {
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		tt = d_temp__[l] * t[k];
		rho[k + l * 70 - 71] = convrt * rmix[k + ng * 70] * p[k] / (
			rgas * tt);
/* cccccccccccc */
/*            write(*,*) rmix(k,ng) */
/* cccccccccccc */
/* L1001: */
	    }

/* L1021: */
	}

	if (istate[1] != 0) {

/* ****        pressure is a variable part of the state vector: */
/*            create a perturbed layer at  base of the atmosphere with */
/*            perturbed pressure, but the same temperature. */
	    rho[*nlev] = convrt * rmix[*nlev + ng * 70] * d_pres__[2] * p[*
		    nlev] / (rgas * t[*nlev]);
	}

/* ****     read wavelength-dependent pressure-induced */
/*         absorption coefficient for each gas. */

	iuabc0 = iuabc[ngt + ng * 3];

	readpab_(&iuabc0, wnmin, wnmax, &npab0, &ntpab, tpab0, wnpab0, pab0, &
		io_end__[*ne], &io_err__[*ne], &wneof[1]);

	npab[ngt + ngpab * 3 - 4] = npab0;
	if (npab[ngt + ngpab * 3 - 4] == 0) {

/* ****      there is no continuum absorption in this spectral range */

	    io_end__[*ne] = 1;
	    nsiext[(*ne << 1) + 1] = 1;
	    nsiext[(*ne << 1) + 2] = 2;
	    wneof[1] = -9999.;
	    wneof[2] = -9999.;
	    wnext[(*ne << 1) + 1] = *wnmin;
	    wnext[(*ne << 1) + 2] = *wnmax;

/* ****      inialize gas optical depth - use nlyr+1 layers just in case */
/*          istate(1) = 1. */

	    i__1 = *ntaug;
	    for (m = 1; m <= i__1; ++m) {
		i__2 = *nlay;
		for (k = 1; k <= i__2; ++k) {
		    dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 1) * 70] = 
			    0.f;
		    dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 2) * 70] = 
			    0.f;
		    dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
		    dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 1) * 70] = 
			    0.f;
		    dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 2) * 70] = 
			    0.f;
		    dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L1041: */
		}
/* L1061: */
	    }

	    return 0;

	}

/* ****     interpolate absorption coefficients to the model atmosphere */
/*         temperature grid */

	if (ntpab > 1) {

/* ****      pressured induced absorption coefficients are defined at */
/*          ntpab temperatures */

	    nxmax = 70;
	    nymax = 1;

	    i__1 = npab[ngt + ngpab * 3 - 4];
	    for (n = 1; n <= i__1; ++n) {
		i__2 = *nt_pd__;
		for (l = 1; l <= i__2; ++l) {
		    i__3 = *nlev;
		    for (k = 2; k <= i__3; ++k) {
			i__4 = ntpab;
			for (mm = 1; mm <= i__4; ++mm) {
			    pab1[mm - 1] = pab0[mm + n * 70 - 71];
/* L1201: */
			}
			tt0[0] = d_temp__[l] * t[k];

			xyinterp_(tpab0, pab1, tt0, &pab[k - 1], &nxmax, &
				nymax, &ntpab, &nymax, &nymax);

/* L1221: */
		    }

		    i__3 = nlyr;
		    for (k = 1; k <= i__3; ++k) {
/* Computing 2nd power */
			r__1 = rho[k + 1 + l * 70 - 71];
/* Computing 2nd power */
			r__2 = rho[k + l * 70 - 71];
			dtaupab[k + (n + (ngt + (ngpab + (l + 2) * 5) * 3) * 
				10000) * 70 - 34300071] = (r__1 * r__1 * pab[
				k] + r__2 * r__2 * pab[k - 1]) * 5e4f * (alt[
				k] - alt[k + 1]);

			if (*ntaug > 1) {

/* ****                    find the (rho**2*pab)dz for mixing ratio */
/*                        perturbation at top of layer k (level k) */

/* Computing 2nd power */
			    r__1 = rho[k + 1 + l * 70 - 71];
/* Computing 2nd power */
			    r__2 = d_rmix__[(ng << 1) + 2] * rho[k + l * 70 - 
				    71];
			    dtaupab[k + (n + (ngt + (ngpab + (l + 4) * 5) * 3)
				     * 10000) * 70 - 34300071] = (r__1 * r__1 
				    * pab[k] + r__2 * r__2 * pab[k - 1]) * 
				    5e4f * (alt[k] - alt[k + 1]);

/* ****                     find the (rho**2)dz for mixing ratio */
/*                         perturbation at bottom of layer k (level k+1) */

/* Computing 2nd power */
			    r__1 = d_rmix__[(ng << 1) + 2] * rho[k + 1 + l * 
				    70 - 71];
/* Computing 2nd power */
			    r__2 = rho[k + l * 70 - 71];
			    dtaupab[k + (n + (ngt + (ngpab + (l + 6) * 5) * 3)
				     * 10000) * 70 - 34300071] = (r__1 * r__1 
				    * pab[k] + r__2 * r__2 * pab[k - 1]) * 
				    5e4f * (alt[k] - alt[k + 1]);
			}

/* L1241: */
		    }

/* L1261: */
		}

		if (istate[1] != 0) {

/* ****             pressure is a variable part of the state vector. */
/*                 find optical depth for perturbed surface pressure. */
/*                 (use pab values for the nomial surface temperature) */

/* Computing 2nd power */
		    r__1 = rho[*nlev + 1 + l * 70 - 71];
/* Computing 2nd power */
		    r__2 = rho[*nlev - 1 + l * 70 - 71];
		    dtaupab[*nlay + (n + (ngt + (ngpab + 15) * 3) * 10000) * 
			    70 - 34300071] = (r__1 * r__1 * pab[*nlev - 1] + 
			    r__2 * r__2 * pab[*nlev - 2]) * 5e4f * (alt[*nlev 
			    - 1] - alt[*nlev + 1]);
		}

/* L1281: */
	    }

	} else {

	    i__1 = npab[ngt + ngpab * 3 - 4];
	    for (n = 1; n <= i__1; ++n) {
		i__2 = *ntaug;
		for (m = 1; m <= i__2; ++m) {
		    i__3 = *nt_pd__;
		    for (l = 1; l <= i__3; ++l) {
			i__4 = nlyr;
			for (k = 1; k <= i__4; ++k) {
/* Computing 2nd power */
			    r__1 = rho[k + 1 + l * 70 - 71];
/* Computing 2nd power */
			    r__2 = rho[k + l * 70 - 71];
			    dtaupab[k + (n + (ngt + (ngpab + (l + (m << 1)) * 
				    5) * 3) * 10000) * 70 - 34300071] = pab0[
				    n * 70 - 70] * 5e4f * (r__1 * r__1 + r__2 
				    * r__2) * (alt[k] - alt[k + 1]);
/* L1601: */
			}

/* L1621: */
		    }

/* L1631: */
		}

		if (istate[1] != 0) {

/* ****             pressure is a variable part of the state vector. */
/*                 find optical depth for perturbed surface pressure. */
/*                 (use pab values for the nomial surface temperature) */

/* Computing 2nd power */
		    r__1 = rho[*nlev + 1 + l * 70 - 71];
/* Computing 2nd power */
		    r__2 = rho[*nlev - 1 + l * 70 - 71];
		    dtaupab[*nlay + (n + (ngt + (ngpab + 15) * 3) * 10000) * 
			    70 - 34300071] = pab0[n * 70 - 70] * 5e4f * (r__1 
			    * r__1 + r__2 * r__2) * (alt[*nlev - 1] - alt[*
			    nlev + 1]);
		}

/* L1641: */
	    }

	}

	i__1 = npab[ngt + ngpab * 3 - 4];
	for (n = 1; n <= i__1; ++n) {
	    wnpab[n + (ngt + ngpab * 3) * 10000 - 40001] = wnpab0[n - 1];
/* L1801: */
	}

	wneof[1] = wnpab[(ngt + ngpab * 3) * 10000 - 40000];
	if (npab[ngt + ngpab * 3 - 4] > 0) {
	    wneof[2] = wnpab[npab[ngt + ngpab * 3 - 4] + (ngt + ngpab * 3) * 
		    10000 - 40001];
	} else {
	    wneof[2] = -99999.;
	}

/* ****     initialize spectral counters */

	nsiext[(*ne << 1) + 1] = 1;
	nsiext[(*ne << 1) + 2] = 2;
	nsipr[ngt + ngpab * 3 - 4] = 0;

/* ****   find the spectral index for the first wn */

	if (wnpab[(ngt + ngpab * 3) * 10000 - 40000] > *wnmin) {
	    wnext[(*ne << 1) + 1] = 0.;
	    wnext[(*ne << 1) + 2] = *wnmin;

/* ****       initialize optical depths. Define nlyr+1 just in case */
/*           pressure is a variable part of the state vector. */

	    i__1 = *ntaug;
	    for (m = 1; m <= i__1; ++m) {
		i__2 = *nlay;
		for (k = 1; k <= i__2; ++k) {
		    dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 1) * 70] = 
			    0.f;
		    dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 2) * 70] = 
			    0.f;
		    dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
		    dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 1) * 70] = 
			    0.f;
		    dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 2) * 70] = 
			    0.f;
		    dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L2001: */
		}
/* L2011: */
	    }

	} else {

/* ****           find the first coefficient within wnmin - wnmax */

	    i__1 = npab[ngt + ngpab * 3 - 4];
	    for (ns = 1; ns <= i__1; ++ns) {
		if (wnpab[ns + (ngt + ngpab * 3) * 10000 - 40001] > *wnmin) {
		    goto L2041;
		}
		nsipr[ngt + ngpab * 3 - 4] = ns;
/* L2021: */
	    }

L2041:
	    wnext[(*ne << 1) + 2] = wnpab[nsipr[ngt + ngpab * 3 - 4] + (ngt + 
		    ngpab * 3) * 10000 - 40001];

	    i__1 = *ntaug;
	    for (m = 1; m <= i__1; ++m) {
		i__2 = *nt_pd__;
		for (l = 1; l <= i__2; ++l) {
		    i__3 = *nlay;
		    for (k = 1; k <= i__3; ++k) {
			dtaug[k + ((*ne + (l + (m << 1)) * 113 << 1) + 2) * 
				70] = dtaupab[k + (nsipr[ngt + ngpab * 3 - 4] 
				+ (ngt + (ngpab + (l + (m << 1)) * 5) * 3) * 
				10000) * 70 - 34300071];
/* L2061: */
		    }

/* L2081: */
		}

/* L2091: */
	    }

	}

    }

/* ****    increment pressure-induced absorption counter. */

L3001:
    ++nsipr[ngt + ngpab * 3 - 4];

/* ****    swap spectral interval counters */

    ni = nsiext[(*ne << 1) + 1];
    nsiext[(*ne << 1) + 1] = nsiext[(*ne << 1) + 2];
    nsiext[(*ne << 1) + 2] = ni;

/* ****  find optical depth at next input wavenumber. */

    if (nsipr[ngt + ngpab * 3 - 4] <= npab[ngt + ngpab * 3 - 4]) {
	wnext[ni + (*ne << 1)] = wnpab[nsipr[ngt + ngpab * 3 - 4] + (ngt + 
		ngpab * 3) * 10000 - 40001];
	dvi = (real) (1. / (wnext[(*ne << 1) + 2] - wnext[(*ne << 1) + 1]));
	i__1 = *ntaug;
	for (m = 1; m <= i__1; ++m) {
	    i__2 = *nt_pd__;
	    for (l = 1; l <= i__2; ++l) {
		i__3 = *nlay;
		for (k = 1; k <= i__3; ++k) {
		    dtaug[k + (ni + (*ne + (l + (m << 1)) * 113 << 1)) * 70] =
			     dtaupab[k + (nsipr[ngt + ngpab * 3 - 4] + (ngt + 
			    (ngpab + (l + (m << 1)) * 5) * 3) * 10000) * 70 - 
			    34300071];
		    dtaugdv[k + (*ne + (l + (m << 1)) * 113) * 70] = (dtaug[k 
			    + ((*ne + (l + (m << 1)) * 113 << 1) + 2) * 70] - 
			    dtaug[k + ((*ne + (l + (m << 1)) * 113 << 1) + 1) 
			    * 70]) * dvi;
/* L3021: */
		}
/* L3041: */
	    }
/* L3051: */
	}
	wneof[2] = wnext[nsiext[(*ne << 1) + 2] + (*ne << 1)];
	if (wnext[ni + (*ne << 1)] < *wn) {
	    goto L3001;
	}

    } else {

/* ****     you've run out of coefficients.  set further values to zero. */

	io_end__[*ne] = 1;
	io_err__[*ne] = 0;
	wneof[2] = wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)];

	wnext[ni + (*ne << 1)] = *wnmax;

/* ****    set all optical depths to zero. */

	i__1 = *ntaug;
	for (m = 1; m <= i__1; ++m) {
	    i__2 = *nlay;
	    for (k = 1; k <= i__2; ++k) {
		dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 1) * 70] = 0.f;
		dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 2) * 70] = 0.f;
		dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
		dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 1) * 70] = 0.f;
		dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 2) * 70] = 0.f;
		dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L3061: */
	    }
/* L3071: */
	}

    }

    return 0;
} /* gasptau_ */

