/* gasetau.f -- translated by f2c (version 20100827).
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

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__4 = 4;

/* Subroutine */ int gasetau_(integer *ne, integer *iuabc, integer *ng0, 
	integer *ngt0, integer *nlev, integer *nlay, integer *nsiext, integer 
	*nt_pd__, integer *istate, integer *ntaug, doublereal *wn, doublereal 
	*wnmin, doublereal *wnmax, doublereal *wnext, real *p, real *t, real *
	alt, real *rmix, real *ratm, real *d_pres__, real *d_temp__, real *
	d_rmix__, real *dtaug, real *dtaugdv, doublereal *wneof, integer *
	io_end__, integer *io_err__)
{
    /* Initialized data */

    static real a0 = 6.02e26f;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    extern /* Subroutine */ int readxsec_(integer *, integer *, real *, 
	    doublereal *, doublereal *, doublereal *, real *, integer *, 
	    integer *);
    static integer k, l, m, n, ng, ni, ns;
    static real dvi;
    static integer ngt;
    static real xsec[2100000]	/* was [14000][3][50] */;
    static integer nsix[150]	/* was [3][50] */, nxmx;
    static real xsec0[14000];
    static integer nxsec[150]	/* was [3][50] */;
    static real atoms[140]	/* was [70][2] */, rhodz[21000]	/* was [70][
	    50][2][3] */;
    static integer iquit;
    static real wgtgs;
    static doublereal wlxsec[14000], wnxsec[2100000]	/* was [14000][3][50] 
	    */;
    static real convert;

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 96, 0, 0, 0 };



/* ccccccccccccccccccccccccc  g a s e t a u  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads wavelength-dependent gas absorption       cc */
/* c    cross-sections for electronic transitions, computes the gas     cc */
/* c    optical depth, and interpolates these values to the             cc */
/* c    appropriate pressure, and temperature grids.                    cc */
/* c                                                                    cc */
/* c    note: this version does not include temperature dependence.     cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c         ne - index of extinction source                            cc */
/* c        ng0 - gas index (1 - ngases)                                cc */
/* c      nt_pd - number of temperature values where absorption         cc */
/* c              coefficients are needed (2 are needed if temperature  cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c       nlay - number of layers were optical depths are needed:      cc */
/* c              nlev-1 if pressure is not a varible part of the state cc */
/* c              vector, nlev if it is (istate(1) = 1)                 cc */
/* c      nt_pd - number of temperature values where absorption         cc */
/* c              coefficients are needed (2 are needed if temperature  cc */
/*               is a variable component of the state vector.          cc */
/* c      wnmin - minimum wavenumber (cm**-1)                           cc */
/* c      wnmax - maximum wavenumber (cm**-1)                           cc */
/* c     d_pres - perturbed surface pressure used when pressure is a    cc */
/* c              variable part of the state vector                     cc */
/* c     d_temp - fractional temperature perturbation used when         cc */
/* c              temperature is a variable componen of the state       cc */
/* c              vector.                                               cc */
/* c     d_rmix - perturbed gas mixing ratio used when rmix is a        cc */
/* c              variable part of the state vector                     cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     io_end - end-of file flag                                      cc */
/* c     io_err - flag specifying error in file                         cc */
/* c      wneof - wavenumber at end of file                             cc */
/* c     wnxsec - wavenumber of each absorption cross section           cc */
/* c       xsec - cross-section at each input wn                        cc */
/* c      nxsec - number of wavenumbers were xsec's are specified       cc */
/* c      dtaug - gas optical depth in each layer and input wavelength  cc */
/* c    dtaugdv - wavenumber gradient of each gas optical depth         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  g a s e t a u  ccccccccccccccccccccccccccccc */




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


/* ****    uv gas absorption cross-sections */


/* ****   partial derivative fractional change */


/* ****    specify avagodro's number (kg/kmole). To convert from */
/*        density to number density - multiply by a0/wgt */

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
    iuabc -= 4;

    /* Function Body */

    ng = *ng0;
    ngt = *ngt0;
    nxmx = 14000;
    iquit = 0;

    if (nsiext[(*ne << 1) + 2] == 0) {

/* ****     read wavelength-dependent gas absorption cross-sections */

	readxsec_(&iuabc[ngt + ng * 3], &nxsec[ngt + ng * 3 - 4], &wgtgs, 
		wnmin, wnmax, wlxsec, xsec0, &nxmx, &iquit);

/* ****     check ordering of absorption coefficients.  Reverse order */
/*         if necessary */

	if (wlxsec[0] < wlxsec[nxsec[ngt + ng * 3 - 4] - 1]) {
	    i__1 = nxsec[ngt + ng * 3 - 4];
	    for (n = 1; n <= i__1; ++n) {
		if (wlxsec[nxsec[ngt + ng * 3 - 4] - n] > 0.) {
		    wnxsec[n + (ngt + ng * 3) * 14000 - 56001] = 1e4 / wlxsec[
			    nxsec[ngt + ng * 3 - 4] - n];
		} else {
		    wnxsec[n + (ngt + ng * 3) * 14000 - 56001] = 1e8;
		}
		xsec[n + (ngt + ng * 3) * 14000 - 56001] = xsec0[nxsec[ngt + 
			ng * 3 - 4] - n];
/* L1001: */
	    }

	} else {

	    i__1 = nxsec[ngt + ng * 3 - 4];
	    for (n = 1; n <= i__1; ++n) {
		if (wlxsec[n - 1] > 0.) {
		    wnxsec[n + (ngt + ng * 3) * 14000 - 56001] = 1e4 / wlxsec[
			    n - 1];
		} else {
		    wnxsec[n + (ngt + ng * 3) * 14000 - 56001] = 1e8;
		}
		xsec[n + (ngt + ng * 3) * 14000 - 56001] = xsec0[n - 1];
/* L1021: */
	    }
	}

/* ****     set start-of-buffer wavenumber */

	if (wneof[1] <= 0. || wnxsec[n + (ngt + ng * 3) * 14000 - 56001] <= 
		wneof[1]) {
	    wneof[1] = wnxsec[n + (ngt + ng * 3) * 14000 - 56001];
	}

/* ****      define factor to convert  rho(kg/m**3)*dz(km) to */
/*          N(atoms/cm**3)*dz(cm) = (a0/wgt)*(1e-6m**3/cm**3)*1e5cm/km */
/*          or 0.1*a0/wgt - note wgt divides out of rho and convert. */

	convert = a0 * .1f / (*ratm * wgtgs);
	i__1 = *nt_pd__;
	for (l = 1; l <= i__1; ++l) {
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		atoms[k - 1] = convert * rmix[k + ng * 70] * p[k] / (d_temp__[
			l] * t[k]);
		if (*ntaug > 1) {
		    atoms[k + 69] = convert * d_rmix__[(ng << 1) + 2] * rmix[
			    k + ng * 70] * p[k] / (d_temp__[l] * t[k]);
		}
/* L2001: */
	    }

	    if (istate[1] != 0) {

/* ****         if pressure is a variable part of the state vector */
/*             create a perturbed layer at the base of the atmosphere */
/*             with 1% more pressure than the surface layer, but the */
/*             same temperature. (needed only for l=1) */

		k = *nlev;
		atoms[k] = convert * d_pres__[2] * rmix[k + ng * 70] * p[k] / 
			(d_temp__[l] * t[k]);
	    }

/* ****         find rho*dz in each layer */

	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		rhodz[k + (ng + (l + 2) * 50) * 70 - 10571] = (atoms[k] + 
			atoms[k - 1]) * .5f * (alt[k] - alt[k + 1]);
		if (*ntaug > 1) {

/* ****                  find the sigma*rho*dz for mixing ratio */
/*                      perturbation at top of layer k (level k) */

		    rhodz[k + (ng + (l + 4) * 50) * 70 - 10571] = (atoms[k] + 
			    atoms[k + 69]) * .5f * (alt[k] - alt[k + 1]);

/* ****                  find the sigma*rho*dz for mixing ratio */
/*                      perturbation at bottom of layer k (level k+1) */

		    rhodz[k + (ng + (l + 6) * 50) * 70 - 10571] = (atoms[k + 
			    70] + atoms[k - 1]) * .5f * (alt[k] - alt[k + 1]);
		}
/* L2201: */
	    }

	    if (istate[1] != 0) {

/* ****           find sigma*rho*dz for lower layer with peturbed pressure */

		rhodz[*nlev + 1 + (ng + (l + 2) * 50) * 70 - 10571] = (atoms[*
			nlev] + atoms[*nlev - 2]) * .5f * (alt[*nlev - 1] - 
			alt[*nlev + 1]);
		if (*ntaug > 1) {
		    rhodz[*nlev + 1 + (ng + (l + 4) * 50) * 70 - 10571] = 
			    rhodz[*nlev + 1 + (ng + (l + 2) * 50) * 70 - 
			    10571];
		    rhodz[*nlev + 1 + (ng + (l + 6) * 50) * 70 - 10571] = 
			    rhodz[*nlev + 1 + (ng + (l + 2) * 50) * 70 - 
			    10571];
		}

	    }

/* L2221: */
	}

/* ****     initialize spectral counters */

	nsiext[(*ne << 1) + 1] = 1;
	nsiext[(*ne << 1) + 2] = 2;
	nsix[ngt + ng * 3 - 4] = 0;

/* ****      find the spectral index for the first wn */

	if (wnxsec[(ngt + ng * 3) * 14000 - 56000] > *wnmin) {
	    wnext[(*ne << 1) + 1] = 0.;
	    wnext[(*ne << 1) + 2] = *wnmin;

/* ****      inialize gas optical depth - use nlev layers just in case */
/*          istate(1) = +/-1. */

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
/* L2401: */
		}
/* L2421: */
	    }

	} else {

/* ****           find the first coefficient within wnmin - wnmax */

	    i__1 = nxsec[ngt + ng * 3 - 4];
	    for (ns = 1; ns <= i__1; ++ns) {
		if (wnxsec[ns + (ngt + ng * 3) * 14000 - 56001] > *wnmin) {
		    goto L2801;
		}
		nsix[ngt + ng * 3 - 4] = ns;
/* L2601: */
	    }

L2801:
	    wnext[(*ne << 1) + 2] = wnxsec[nsix[ngt + ng * 3 - 4] + (ngt + ng 
		    * 3) * 14000 - 56001];
	    i__1 = *nt_pd__;
	    for (l = 1; l <= i__1; ++l) {
		i__2 = *ntaug;
		for (m = 1; m <= i__2; ++m) {
		    i__3 = *nlay;
		    for (k = 1; k <= i__3; ++k) {
			s_wsle(&io___21);
			do_lio(&c__5, &c__1, (char *)&(*wn), (ftnlen)sizeof(
				doublereal));
			do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(
				integer));
			do_lio(&c__3, &c__1, (char *)&ng, (ftnlen)sizeof(
				integer));
			do_lio(&c__3, &c__1, (char *)&l, (ftnlen)sizeof(
				integer));
			do_lio(&c__3, &c__1, (char *)&m, (ftnlen)sizeof(
				integer));
			do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(
				integer));
			do_lio(&c__4, &c__1, (char *)&rhodz[k + (ng + (l + (m 
				<< 1)) * 50) * 70 - 10571], (ftnlen)sizeof(
				real));
			e_wsle();
			dtaug[k + ((*ne + (l + (m << 1)) * 113 << 1) + 2) * 
				70] = xsec[nsix[ngt + ng * 3 - 4] + (ngt + ng 
				* 3) * 14000 - 56001] * rhodz[k + (ng + (l + (
				m << 1)) * 50) * 70 - 10571];
/* L2821: */
		    }
/* L2841: */
		}
/* L2861: */
	    }
	}

    }

/* ****    increment continuum absorption counter. */

L3001:
    ++nsix[ngt + ng * 3 - 4];

/* ****    swap spectral interval counters */

    ni = nsiext[(*ne << 1) + 1];
    nsiext[(*ne << 1) + 1] = nsiext[(*ne << 1) + 2];
    nsiext[(*ne << 1) + 2] = ni;

/* ****  find optical depth at next input wavenumber */

    if (nsix[ngt + ng * 3 - 4] <= nxsec[ngt + ng * 3 - 4]) {
	wnext[ni + (*ne << 1)] = wnxsec[nsix[ngt + ng * 3 - 4] + (ngt + ng * 
		3) * 14000 - 56001];
	dvi = (real) (1. / (wnext[(*ne << 1) + 2] - wnext[(*ne << 1) + 1]));
	i__1 = *nt_pd__;
	for (l = 1; l <= i__1; ++l) {

/* ****         find an additional value if pressure is a variable */
/*             part of the state vector */

	    i__2 = *ntaug;
	    for (m = 1; m <= i__2; ++m) {
		i__3 = *nlay;
		for (k = 1; k <= i__3; ++k) {
		    dtaug[k + (ni + (*ne + (l + (m << 1)) * 113 << 1)) * 70] =
			     xsec[nsix[ngt + ng * 3 - 4] + (ngt + ng * 3) * 
			    14000 - 56001] * rhodz[k + (ng + (l + (m << 1)) * 
			    50) * 70 - 10571];
		    dtaugdv[k + (*ne + (l + (m << 1)) * 113) * 70] = (dtaug[k 
			    + ((*ne + (l + (m << 1)) * 113 << 1) + 2) * 70] - 
			    dtaug[k + ((*ne + (l + (m << 1)) * 113 << 1) + 1) 
			    * 70]) * dvi;
/* L3021: */
		}
/* L3031: */
	    }
/* L3041: */
	}
	wneof[2] = wnext[nsiext[(*ne << 1) + 2] + (*ne << 1)];
	if (wnext[ni + (*ne << 1)] < *wn) {
	    goto L3001;
	}

    } else {

/* ****     you've run out of coefficients.  punt! */

	io_end__[*ne] = 1;
	io_err__[*ne] = 0;
	wneof[2] = wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)];

	wnext[ni + (*ne << 1)] = *wnmax;

	i__1 = *nlay;
	for (k = 1; k <= i__1; ++k) {
	    dtaug[k + (ni + (*ne + ((m << 1) + 1) * 113 << 1)) * 70] = 0.f;
	    dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
	    dtaug[k + (ni + (*ne + ((m << 1) + 2) * 113 << 1)) * 70] = 0.f;
	    dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L3061: */
	}
    }

    return 0;
} /* gasetau_ */

