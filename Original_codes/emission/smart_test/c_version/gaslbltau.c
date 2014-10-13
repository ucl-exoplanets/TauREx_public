/* gaslbltau.f -- translated by f2c (version 20100827).
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
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__5 = 5;

/* Subroutine */ int gaslbltau_(integer *ne, integer *iuabc, integer *ng0, 
	integer *ngt0, integer *nlev, integer *npg, integer *ntg, integer *
	nbgas, integer *nsiext, integer *nlay, integer *istate, integer *
	nt_pd__, integer *ntaug, doublereal *wn, doublereal *wnmin, 
	doublereal *wnmax, doublereal *wnext, real *p, real *t, real *z__, 
	real *grav, real *rmix, real *d_pres__, real *d_temp__, real *
	d_rmix__, real *dtaug, real *dtaugdv, doublereal *wneof, integer *
	io_end__, integer *io_err__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void);
    double log(doublereal);
    integer i_nint(real *);
    double exp(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    extern /* Subroutine */ int xyinterp_(real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer i__, k, l, m, n1, n2;
    static real p0[70], p5, z1[70];
    static integer ng;
    static real at[70];
    static integer ni;
    static real bt[70], rg[3500]	/* was [70][50] */;
    static doublereal dv;
    static real tt, dt0, dt02, dvi;
    static integer itl, ngt;
    static real tol, stm, abc0[420]	/* was [70][3][2] */, stp, wns, abc1[
	    210]	/* was [70][3] */, abc2[420]	/* was [70][2][3] */, 
	    dnu0, delt, ppas[10500]	/* was [70][3][50] */, tatm[210]	
	    /* was [70][3] */;
    static integer nlyr;
    static real tatm1[31500]	/* was [70][3][3][50] */, sigma;
    static integer int_p__[150]	/* was [3][50] */, int_t__[150]	/* was [3][50]
	     */, itlev[10500]	/* was [70][3][50] */;
    static real pstar[1];
    static integer nxmax, nymax;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, 0, 0 };
    static cilist io___9 = { 0, 0, 0, 0, 0 };
    static cilist io___12 = { 0, 0, 0, 0, 0 };
    static cilist io___26 = { 1, 0, 1, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };



/* ccccccccccccccccccccccc  g a s l b l t a u  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads gas absorption coefficients, inter-       cc */
/* c    polates them to the appropriate pressure and temperature        cc */
/* c    grids, and computes the gas optical depth, and it derivative    cc */
/* c    with respect to wavenumber                                      cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    ntemp - Maximum number of temperature profiles for gas          cc */
/* c            absorption coefficients.				      cc */
/* c         ne - index of extinction source                            cc */
/* c         ng - gas index (1 - ngases)                                cc */
/* c        ngt - number of different absorption coefficeint types      cc */
/* c              for this gas.                                         cc */
/* c       nlev - number of levels in the model atmosphere (<kp).       cc */
/* c       nlay - number of layers were optical depths are needed:      cc */
/* c              nlev-1 if pressure is not a varible part of the state cc */
/* c              vector, nlev if it is (istate(1) = 1)                 cc */
/* c      nt_pd - number of temperature values where absorption         cc */
/* c              coefficients are needed (2 are needed if temperature  cc */
/*               is a variable component of the state vector.          cc */
/* c      ntaug - number of optical depth profile of each gas:          cc */
/* c                - if istate not equal to 3, ntaug = 1               cc */
/* c                - if istate equal to 3, ntaug = 3                   cc */
/* c      wnmin - minimum wavenumber (cm**-1)                           cc */
/* c      wnmax - maximum wavenumber (cm**-1)                           cc */
/* c         p0 - pressures of input absorption coefficients            cc */
/* c         t0 - temperatures of input absorption coefficients         cc */
/* c      rmix0 - mixing ratios of input absorption coefficients        cc */
/* c         wn - wavenumber of each gas absorption coefficient         cc */
/* c        abs - input absorption coefficient                          cc */
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
/* c      dtaug - gas optical depth in each layer and input wavelength  cc */
/* c    dtaugdv - wavenumber gradient of each gas optical depth         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccc  g a s l b l t a u  ccccccccccccccccccccccccccc */





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







/* ****    atmospheric structure and gas mixing ratios */


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
    --grav;
    --z__;
    --t;
    --p;
    wnext -= 3;
    --istate;
    nsiext -= 3;
    nbgas -= 4;
    ntg -= 4;
    npg -= 4;
    iuabc -= 4;

    /* Function Body */
    tol = 1e-7f;
    p5 = .5f;
    nlyr = *nlev - 1;
    ng = *ng0;
    ngt = *ngt0;

    if (nsiext[(*ne << 1) + 2] == 0) {

/* ****    initialize gas spectral interval counters */

	nsiext[(*ne << 1) + 1] = 1;
	nsiext[(*ne << 1) + 2] = 2;
	wnext[(*ne << 1) + 1] = *wnmin;
	wnext[(*ne << 1) + 2] = *wnmin;

/* ****   open read the atmospheric structure information */

	io___6.ciunit = iuabc[ngt + ng * 3];
	s_rsue(&io___6);
	i__1 = npg[ngt + ng * 3];
	for (k = 1; k <= i__1; ++k) {
	    do_uio(&c__1, (char *)&p0[k - 1], (ftnlen)sizeof(real));
	}
	e_rsue();
	io___9.ciunit = iuabc[ngt + ng * 3];
	s_rsue(&io___9);
	i__1 = npg[ngt + ng * 3];
	for (k = 1; k <= i__1; ++k) {
	    i__2 = nbgas[ngt + ng * 3];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_uio(&c__1, (char *)&rg[k + i__ * 70 - 71], (ftnlen)sizeof(
			real));
	    }
	}
	e_rsue();
	io___12.ciunit = iuabc[ngt + ng * 3];
	s_rsue(&io___12);
	i__2 = ntg[ngt + ng * 3];
	for (l = 1; l <= i__2; ++l) {
	    i__1 = npg[ngt + ng * 3];
	    for (k = 1; k <= i__1; ++k) {
		do_uio(&c__1, (char *)&tatm[k + l * 70 - 71], (ftnlen)sizeof(
			real));
	    }
	}
	e_rsue();

/* *****    determine if absorption coeficients must be interpolated */
/*         to a new pressure grid: */

	if (npg[ngt + ng * 3] != *nlev) {

/* ****      the input atmosphere has a different number of levels. */
/*          a pressure interpolation is required. */

	    int_p__[ngt + ng * 3 - 4] = 1;
	} else {

/* ****     the input atmosphere has the same number of levels as the */
/*         absorption coefficient file.  Check each level */

	    int_p__[ngt + ng * 3 - 4] = 0;
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {
		if ((r__1 = p0[k - 1] - p[k], dabs(r__1)) > p[k] * .01f) {
		    int_p__[ngt + ng * 3 - 4] = 1;
		}
/* L1201: */
	    }
	}

	if (int_p__[ngt + ng * 3 - 4] != 0) {

/* ****      interpolate temperatures to the output p grid. */

	    i__1 = npg[ngt + ng * 3];
	    for (k = 1; k <= i__1; ++k) {
		ppas[k + (ngt + ng * 3) * 70 - 281] = p0[k - 1];
		if (ppas[k + (ngt + ng * 3) * 70 - 281] != 0.f) {
		    z1[k - 1] = log(ppas[k + (ngt + ng * 3) * 70 - 281]);
		} else {
		    z1[k - 1] = -80.f;
		}
/* L2421: */
	    }

/* ****       interpolate absorption coefficients to output p grid. */

	    nxmax = 70;
	    nymax = 3;

	    xyinterp_(z1, tatm, &z__[1], &tatm1[((ngt + ng * 3) * 3 + 1) * 70 
		    - 910], &nxmax, &nymax, &npg[ngt + ng * 3], nlev, &ntg[
		    ngt + ng * 3]);
	} else {
	    i__1 = ntg[ngt + ng * 3];
	    for (l = 1; l <= i__1; ++l) {
		i__2 = *nlev;
		for (k = 1; k <= i__2; ++k) {
		    tatm1[k + (l + (ngt + ng * 3) * 3) * 70 - 911] = tatm[k + 
			    l * 70 - 71];
/* L2441: */
		}
/* L2461: */
	    }
	}

/* ****    determine if a temperature interpolation is needed. */

	if (ntg[ngt + ng * 3] > 1) {
	    int_t__[ngt + ng * 3 - 4] = 0;
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {

/* ****            find the index of the temperature that is closest */

		dt0 = tatm1[k + ((ngt + ng * 3) * 3 + 2) * 70 - 911] - tatm1[
			k + ((ngt + ng * 3) * 3 + 1) * 70 - 911];
		r__1 = (t[k] - tatm1[k + ((ngt + ng * 3) * 3 + 1) * 70 - 911])
			 / dt0;
		itlev[k + (ngt + ng * 3) * 70 - 281] = i_nint(&r__1) + 1;

/* ****           check to see if itlev is with range of tatm */

		if (itlev[k + (ngt + ng * 3) * 70 - 281] < 1) {
		    itlev[k + (ngt + ng * 3) * 70 - 281] = 1;
		} else {
		    if (itlev[k + (ngt + ng * 3) * 70 - 281] > ntg[ngt + ng * 
			    3]) {
			itlev[k + (ngt + ng * 3) * 70 - 281] = ntg[ngt + ng * 
				3];
		    }
		}

/* ****           determine if a temperature interpolation is needed */

		if ((r__1 = tatm1[k + (itlev[k + (ngt + ng * 3) * 70 - 281] + 
			(ngt + ng * 3) * 3) * 70 - 911] - t[k], dabs(r__1)) > 
			t[k] * 5e-5f || istate[2] != 0) {
		    int_t__[ngt + ng * 3 - 4] = 1;
		}
/* L2501: */
	    }
	} else {
	    int_t__[ngt + ng * 3 - 4] = 0;
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {
		itlev[k + (ngt + ng * 3) * 70 - 281] = 1;
/* L2521: */
	    }
	}

/* ****   initialize absorption coefficients for this gas */

	ni = nsiext[(*ne << 1) + 2];
	wnext[(*ne << 1) + 2] = *wnmin;
	i__1 = *ntaug;
	for (m = 1; m <= i__1; ++m) {
	    i__2 = *nlay;
	    for (k = 1; k <= i__2; ++k) {
		dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 1) * 70] = 0.f;
		dtaug[k + ((*ne + ((m << 1) + 1) * 113 << 1) + 2) * 70] = 0.f;
		dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 1) * 70] = 0.f;
		dtaug[k + ((*ne + ((m << 1) + 2) * 113 << 1) + 2) * 70] = 0.f;
/* L2601: */
	    }
/* L2621: */
	}

    }

/* ****    swap spectral interval counters */

L2801:
    ni = nsiext[(*ne << 1) + 1];
    nsiext[(*ne << 1) + 1] = nsiext[(*ne << 1) + 2];
    nsiext[(*ne << 1) + 2] = ni;

/* ****     read the next absorption coefficient for this gas */

    io___26.ciunit = iuabc[ngt + ng * 3];
    i__1 = s_rsue(&io___26);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&wns, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&dnu0, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L100001;
    }
    i__2 = ntg[ngt + ng * 3];
    for (l = 1; l <= i__2; ++l) {
	i__3 = npg[ngt + ng * 3];
	for (k = 1; k <= i__3; ++k) {
	    i__1 = do_uio(&c__1, (char *)&abc0[k + (l + ni * 3) * 70 - 281], (
		    ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100001;
	    }
	}
    }
    i__1 = e_rsue();
L100001:
    if (i__1 < 0) {
	goto L4001;
    }
    if (i__1 > 0) {
	goto L4201;
    }

/* ****   determine if this segment is in desired spectral window */

    wnext[ni + (*ne << 1)] = wns;
    if (wnext[ni + (*ne << 1)] <= *wnmin) {
	wneof[1] = wnext[ni + (*ne << 1)];
    }
    if (wnext[ni + (*ne << 1)] < *wn) {
	goto L2801;
    }

    wneof[2] = wnext[ni + (*ne << 1)];

/* ****    interpolate absorption coefficients to appropriate output grid. */

    if (wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)] <= *wnmin) {
	n1 = 1;
    } else {
	n1 = 2;
    }
    n2 = 2;

    nxmax = 70;
    nymax = 3;

    i__1 = n2;
    for (ni = n1; ni <= i__1; ++ni) {
	if (int_p__[ngt + ng * 3 - 4] == 1) {

/* ****    interpolate absorption coefficients to p grid */

	    xyinterp_(&ppas[(ngt + ng * 3) * 70 - 280], &abc0[(nsiext[ni + (*
		    ne << 1)] * 3 + 1) * 70 - 280], &p[1], abc1, &nxmax, &
		    nymax, &npg[ngt + ng * 3], nlev, &ntg[ngt + ng * 3]);

	} else {

/* ****        no p interpolation needed.  Just load abc1 array */

	    i__3 = ntg[ngt + ng * 3];
	    for (l = 1; l <= i__3; ++l) {
		i__2 = *nlay + 1;
		for (k = 1; k <= i__2; ++k) {
		    abc1[k + l * 70 - 71] = abc0[k + (l + nsiext[ni + (*ne << 
			    1)] * 3) * 70 - 281];
/* L3101: */
		}
/* L3121: */
	    }

	}

	if (istate[1] != 0) {

/* ****         interpolate abc to layer (1.0+p_pert(1))*p(nlev) */

	    pstar[0] = d_pres__[2] * p[*nlev];
	    xyinterp_(&ppas[(ngt + ng * 3) * 70 - 280], &abc0[(nsiext[ni + (*
		    ne << 1)] * 3 + 1) * 70 - 280], pstar, &abc1[*nlev], &
		    nxmax, &nymax, &npg[ngt + ng * 3], &c__1, &ntg[ngt + ng * 
		    3]);

/* ****         define tatm1 at nlev+1 */

	    i__3 = ntg[ngt + ng * 3];
	    for (l = 1; l <= i__3; ++l) {
		tatm1[*nlev + 1 + (l + (ngt + ng * 3) * 3) * 70 - 911] = 
			tatm1[*nlev + (l + (ngt + ng * 3) * 3) * 70 - 911];
/* L3141: */
	    }

	}

/* ****       interpolate absorption coefficients to output T grid. */
/*           (assume that log(abc) is quadratic in t.) */
/*           Weight results by mass-mixing-ratio/gravity (cgs units) */

	if (int_t__[ngt + ng * 3 - 4] >= 1) {
	    i__3 = *nlay + 1;
	    for (k = 1; k <= i__3; ++k) {
		dt0 = tatm1[k + ((ngt + ng * 3) * 3 + 2) * 70 - 911] - tatm1[
			k + ((ngt + ng * 3) * 3 + 1) * 70 - 911];
		dt02 = dt0 * dt0;
		i__2 = ntg[ngt + ng * 3] - 1;
		for (l = 2; l <= i__2; ++l) {
		    stp = log(abc1[k + (l + 1) * 70 - 71] / abc1[k + l * 70 - 
			    71]);
		    stm = log(abc1[k + (l - 1) * 70 - 71] / abc1[k + l * 70 - 
			    71]);
		    at[l - 1] = p5 * (stp - stm) / dt0;
		    bt[l - 1] = p5 * (stp + stm) / dt02;
/* L3201: */
		}

		i__2 = *nt_pd__;
		for (l = 1; l <= i__2; ++l) {
		    tt = d_temp__[l] * t[k];
		    r__1 = (tt - tatm1[k + ((ngt + ng * 3) * 3 + 1) * 70 - 
			    911]) / dt0;
		    itl = i_nint(&r__1) + 1;
		    if (itl < 2) {
			itl = 2;
		    }
		    if (itl > ntg[ngt + ng * 3] - 1) {
			itl = ntg[ngt + ng * 3] - 1;
		    }
		    delt = tt - tatm1[k + (itl + (ngt + ng * 3) * 3) * 70 - 
			    911];
		    sigma = abc1[k + itl * 70 - 71] * .1f * exp(at[itl - 1] * 
			    delt + bt[itl - 1] * delt * delt);
		    abc2[k + (l + 2) * 70 - 211] = sigma * rmix[k + ng * 70] /
			     grav[k];
		    if (*ntaug > 1) {
			abc2[k + (l + 4) * 70 - 211] = sigma * d_rmix__[(ng <<
				 1) + 2] * rmix[k + ng * 70] / grav[k];
		    }
/* L3221: */
		}

/* L3241: */
	    }

	} else {

/* ****       values are only specified at one temperature */

	    i__3 = *nlay + 1;
	    for (k = 1; k <= i__3; ++k) {
		abc2[k - 1] = abc1[k + itlev[k + (ngt + ng * 3) * 70 - 281] * 
			70 - 71] * .1f * rmix[k + ng * 70] / grav[k];

		if (*ntaug > 1) {

/* ****               find the absorption for perturbed profile */

		    abc2[k + 139] = abc1[k + itlev[k + (ngt + ng * 3) * 70 - 
			    281] * 70 - 71] * .1f * d_rmix__[(ng << 1) + 2] * 
			    rmix[k + ng * 70] / grav[k];
		}
/* L3401: */
	    }

	}

/* ****         find the gas optical depth in each layer. */

	i__3 = *nt_pd__;
	for (l = 1; l <= i__3; ++l) {
	    i__2 = nlyr;
	    for (k = 1; k <= i__2; ++k) {
		dtaug[k + (nsiext[ni + (*ne << 1)] + (*ne + (l + 2) * 113 << 
			1)) * 70] = p5 * (abc2[k + 1 + (l + 2) * 70 - 211] + 
			abc2[k + (l + 2) * 70 - 211]) * (p[k + 1] - p[k]);

		if (*ntaug > 1) {

/* ****                  find the optical depth for mixing ratio */
/*                      perturbation at top of layer k (level k) */

		    dtaug[k + (nsiext[ni + (*ne << 1)] + (*ne + (l + 4) * 113 
			    << 1)) * 70] = p5 * (abc2[k + 1 + (l + 2) * 70 - 
			    211] + abc2[k + (l + 4) * 70 - 211]) * (p[k + 1] 
			    - p[k]);

/* ****                  find the optical depth for mixing ratio */
/*                      perturbation at bottom of layer k (level k+1) */

		    dtaug[k + (nsiext[ni + (*ne << 1)] + (*ne + (l + 6) * 113 
			    << 1)) * 70] = p5 * (abc2[k + 1 + (l + 4) * 70 - 
			    211] + abc2[k + (l + 2) * 70 - 211]) * (p[k + 1] 
			    - p[k]);
		}
/* L3421: */
	    }

/* L3441: */
	}

	if (istate[1] != 0) {

/* ****         surface pressure is a variable component of the */
/*             state vector - find perturbation optical depth */
/*             for the lowest pressure level. */

	    k = *nlev;
	    dtaug[k + (nsiext[ni + (*ne << 1)] + (*ne + 339 << 1)) * 70] = p5 
		    * (abc2[k] + abc2[k - 2]) * (d_pres__[2] * p[k] - p[k - 1]
		    );
	}

/* L3481: */
    }

/* ****     compute the derivative of the optical depth w.r.t. wavenumber */

    dv = wnext[(*ne << 1) + 2] - wnext[(*ne << 1) + 1];
    if (abs(dv) > tol * wnext[(*ne << 1) + 1]) {
	dvi = (real) (1. / dv);
	i__1 = *ntaug;
	for (m = 1; m <= i__1; ++m) {
	    i__3 = *nt_pd__;
	    for (l = 1; l <= i__3; ++l) {
		i__2 = *nlay;
		for (k = 1; k <= i__2; ++k) {
		    dtaugdv[k + (*ne + (l + (m << 1)) * 113) * 70] = (dtaug[k 
			    + ((*ne + (l + (m << 1)) * 113 << 1) + 2) * 70] - 
			    dtaug[k + ((*ne + (l + (m << 1)) * 113 << 1) + 1) 
			    * 70]) * dvi;
/* L3801: */
		}
/* L3821: */
	    }
/* L3841: */
	}

    } else {

	i__1 = *ntaug;
	for (m = 1; m <= i__1; ++m) {
	    i__3 = *nlay;
	    for (k = 1; k <= i__3; ++k) {
		dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
		dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L3861: */
	    }
/* L3881: */
	}
    }

    return 0;

/* ***   end of file encountered */

L4001:
    io_end__[*ne] = 1;

    wneof[2] = wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)];

/* ****     set all gas optical depths and derivatives to zero */

    wnext[nsiext[(*ne << 1) + 2] + (*ne << 1)] = *wnmax * 1.01;
    i__1 = *ntaug;
    for (m = 1; m <= i__1; ++m) {
	i__3 = *nlay;
	for (k = 1; k <= i__3; ++k) {
	    dtaug[k + (nsiext[(*ne << 1) + 1] + (*ne + ((m << 1) + 1) * 113 <<
		     1)) * 70] = 0.f;
	    dtaug[k + (nsiext[(*ne << 1) + 2] + (*ne + ((m << 1) + 1) * 113 <<
		     1)) * 70] = 0.f;
	    dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
	    dtaug[k + (nsiext[(*ne << 1) + 1] + (*ne + ((m << 1) + 2) * 113 <<
		     1)) * 70] = 0.f;
	    dtaug[k + (nsiext[(*ne << 1) + 2] + (*ne + ((m << 1) + 2) * 113 <<
		     1)) * 70] = 0.f;
	    dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L4021: */
	}
/* L4041: */
    }

    return 0;

/* ****   file read error encountered */

L4201:
    io_err__[*ne] = 1;
    wneof[2] = wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)];
    s_wsle(&io___46);
    do_lio(&c__9, &c__1, "gaslbltau: error reading unit iuabc=", (ftnlen)36);
    do_lio(&c__3, &c__1, (char *)&iuabc[ngt + ng * 3], (ftnlen)sizeof(integer)
	    );
    do_lio(&c__9, &c__1, " at wavenumber", (ftnlen)14);
    do_lio(&c__5, &c__1, (char *)&wneof[2], (ftnlen)sizeof(doublereal));
    e_wsle();

    wnext[nsiext[(*ne << 1) + 2] + (*ne << 1)] = *wnmax * 1.01;

/* ****     set all gas optical depths and derivatives to zero */

    i__1 = *ntaug;
    for (m = 1; m <= i__1; ++m) {
	i__3 = *nlay;
	for (k = 1; k <= i__3; ++k) {
	    dtaug[k + (nsiext[(*ne << 1) + 1] + (*ne + ((m << 1) + 1) * 113 <<
		     1)) * 70] = 0.f;
	    dtaug[k + (nsiext[(*ne << 1) + 2] + (*ne + ((m << 1) + 1) * 113 <<
		     1)) * 70] = 0.f;
	    dtaugdv[k + (*ne + ((m << 1) + 1) * 113) * 70] = 0.f;
	    dtaug[k + (nsiext[(*ne << 1) + 1] + (*ne + ((m << 1) + 2) * 113 <<
		     1)) * 70] = 0.f;
	    dtaug[k + (nsiext[(*ne << 1) + 2] + (*ne + ((m << 1) + 2) * 113 <<
		     1)) * 70] = 0.f;
	    dtaugdv[k + (*ne + ((m << 1) + 2) * 113) * 70] = 0.f;
/* L4221: */
	}
/* L4241: */
    }

    return 0;
} /* gaslbltau_ */

