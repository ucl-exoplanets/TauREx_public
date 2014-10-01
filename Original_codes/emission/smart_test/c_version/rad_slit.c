/* rad_slit.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rad_slit__(integer *iuout, integer *nz0, integer *ifrm, 
	integer *islit, integer *iscwt, integer *iord, integer *nspt, 
	doublereal *wnmin, doublereal *wnmax, doublereal *dwn, doublereal *
	width, real *points, doublereal *wnz, real *z__, integer *iquit)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void),
	     f_clos(cllist *);

    /* Local variables */
    static doublereal dwn_dist__, a[15284682]	/* was [5162][21][141] */;
    static integer j, n;
    static doublereal w[2961]	/* was [21][141] */;
    static integer i0, i1;
    static real z0[727842]	/* was [5162][141] */;
    static integer nc, in, nn;
    static real wl;
    static integer nw, nz;
    static doublereal wt;
    static integer in0[141], in1[141], in2[141];
    static doublereal wn0[141], wn1[141];
    static real zm1[727842]	/* was [5162][141] */, fcn[2];
    static integer in01;
    static doublereal dnu;
    static integer nwn;
    static doublereal dnu0[141], wni0[141], wnm1[141];
    static integer nfcn;
    static doublereal dnui, hwhm, dist;
    static integer ndwn;
    static real dwnz;
    static integer nout[141], i0new[141], i1new[141], i2new[141];
    static real dfdnu[2];
    static doublereal delnu[2], delwn[141];
    static real wni_o__[141], wn_io__[2961]	/* was [21][141] */, spect[
	    5162];
    static doublereal dzdwn[5162], dnuwt;
    static integer nwnwt;
    static doublereal wnmin0[141], wnmax0[141];
    static real spectc[5162];

    /* Fortran I/O blocks */
    static cilist io___30 = { 0, 6, 0, "(1x,1a,/,1x,1a,/,1x,1a/,1x,1a)", 0 };
    static cilist io___31 = { 0, 6, 0, "(1a,2i10,3(1pe13.5))", 0 };
    static cilist io___46 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___47 = { 0, 0, 0, 0, 0 };
    static cilist io___48 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___49 = { 0, 0, 0, 0, 0 };
    static cilist io___50 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___51 = { 0, 0, 0, 0, 0 };
    static cilist io___52 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___53 = { 0, 0, 0, 0, 0 };
    static cilist io___55 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___56 = { 0, 0, 0, 0, 0 };
    static cilist io___57 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___58 = { 0, 0, 0, 0, 0 };
    static cilist io___65 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___66 = { 0, 0, 0, 0, 0 };
    static cilist io___67 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___68 = { 0, 0, 0, 0, 0 };
    static cilist io___69 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___70 = { 0, 0, 0, 0, 0 };
    static cilist io___71 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___72 = { 0, 0, 0, 0, 0 };



/* ccccccccccccccccccccccccc  r a d _  s l i t  cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine uses a weighted-average scheme to map randomly  cc */
/* c    spaced data onto a regular, but arbitrarily-spaced output grid. cc */
/* c                                                                    cc */
/* c    the output of this routine is actually the weighted sum of all  cc */
/* c    contributions at each point on the output grid, and the sum of  cc */
/* c    all weights that contribute at that point.  the weighted sum at cc */
/* c    each point must be divided by the sum of all weights at that    cc */
/* c    point by the calling program after all input points are         cc */
/* c    included.                                                       cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c     ifrm - output file format: 1) ascii, 2) binary.                cc */
/* c    iord - ordering of input grid: 1) wavelength, 2) wavenumber.    cc */
/* c    iscwt - spectral resolution flag: 0) constant,                  cc */
/* c                                      1) linear with wn             cc */
/* c    width - half-width at half max for weighting function.          cc */
/* c            these quantities are also used to normalize the         cc */
/* c            spatial/temporal domain of the grid.                    cc */
/* c      wn1 - wavenumber of each input data point.                    cc */
/* c            (only the first 2 are used for the output map)          cc */
/* c        z - vector containing values of input quantity at wn.       cc */
/* c     nspt - number of z values at each wavenumber.                  cc */
/* c    wn_io - output wavenumber grid.                                 cc */
/* c    wnmin - minimum wavenumber in input and output grid             cc */
/* c    wnmax - maximum wavenumber in input and output grid             cc */
/* c      dwn - spacing of output spectral grid (cm**-1)                cc */
/* c                                                                    cc */
/* c     islit - index specifying type of weighting function.           cc */
/* c                 1) boxcar (constant) weighting within r0, ry, r0   cc */
/* c                 2) triangular (approximate slit spectrometer)      cc */
/* c     nfcn - number of values in fcn table.                          cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c        wn - wavenumber at each output sample point.                cc */
/* c     spect - value of each ouput spectral quantity at point wn.     cc */
/* c                                                                    cc */
/* c    note: this version of slit is designed to accomodate multiple   cc */
/* c          streams and solar zenith angles.                          cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  r a d _  s l i t  cccccccccccccccccccccccccc */




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

















    /* Parameter adjustments */
    --z__;
    --points;

    /* Function Body */
    nz = *nz0;

    if (*iquit < 0) {
	return 0;
    }

    if (points[nz] <= 1.f) {

/* ****    initialize index pointers for wrap-around output buffer. */

	in0[nz - 1] = 1;
	in2[nz - 1] = 0;
	wni0[nz - 1] = *wnmin;
	nout[nz - 1] = 0;

/* ****     find index of last contributing element.  Assume the ouput */
/*         grid is equally-spaced with spacing, dwn. */

	nwn = (integer) ((*wnmax - *wnmin) / *dwn) + 1;
	in1[nz - 1] = nwn;
	if (in1[nz - 1] > 21) {
	    in1[nz - 1] = 21;
	}
	delwn[nz - 1] = 0.f;
	wn0[nz - 1] = -9999.f;

/* ****      initialize output arrays for range extending from i0 to i1: */

	i__1 = in1[nz - 1];
	for (n = in0[nz - 1]; n <= i__1; ++n) {
	    wn_io__[n + nz * 21 - 22] = (real) (*wnmin + delwn[nz - 1]);
	    delwn[nz - 1] += *dwn;
	    w[n + nz * 21 - 22] = 0.f;
	    i__2 = *nspt;
	    for (j = 1; j <= i__2; ++j) {
		a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L1001: */
	    }
/* L1021: */
	}

/* ****     initialize internal variables */

	wn0[nz - 1] = *wnz;
	wnm1[nz - 1] = *wnz;
	i__1 = *nspt;
	for (nc = 1; nc <= i__1; ++nc) {
	    z0[nc + nz * 5162 - 5163] = z__[nc];
	    zm1[nc + nz * 5162 - 5163] = z__[nc];
/* L1041: */
	}

/* ****          you can't do any more, so return */

/*          write(*,'(/,1a,i5)') */
/*     -      'Returning from rad_slit after initialing unit: ',iuout */

	return 0;

    } else {

/* ****      This value's contribution to the output spectrum depends on */
/*          the spectral interval width that it represents.  Because this */
/*          width is not uniform, it must be estimated from the */
/*          distance between adjacent intervals. */

	if (points[nz] == 2.f) {

/* ****         now you know the width of interval 1.  express all other */
/*             interval widths relative to this width. */

	    dnu0[nz - 1] = (d__1 = *wnz - wn0[nz - 1], abs(d__1));
	    dnuwt = 1.f;
	} else {

	    if (dnu0[nz - 1] != 0.f) {
		dnuwt = (*wnz - wnm1[nz - 1]) * .5f / dnu0[nz - 1];
	    } else {
		dnuwt = 0.f;
	    }

	}

    }

/* ****     find the distance between point wn0 and the */
/*         previous point, wnm1. */

    dwn_dist__ = wn0[nz - 1] - wnm1[nz - 1];

/* ****    determine if the slit function must be re-defined */

    if (points[nz] <= 2.f || *iscwt == 1 && dwn_dist__ > wn0[nz - 1] * .005f) 
	    {

/* ****     define the half-width of the weighting function */

	if (*iscwt == 1) {
	    if (wn0[nz - 1] > 0.f) {
		hwhm = *width * .5f * wn0[nz - 1];
	    } else {
		hwhm = *width * .5f * *wnmin;
	    }
	} else {
	    hwhm = *width;
	}

/* ****           d e f i n e    s l i t    f u n c t i o n */

/* ****     define weighting functions and their derivative. */

	if (*islit == 1) {

/* ****          b o x c a r    f u n c t i o n */

	    nfcn = 2;
	    dnu = hwhm;
	    dnui = 1. / dnu;
	    delnu[0] = 0.f;
	    fcn[0] = 1.f;
	    dfdnu[0] = 0.f;
	    delnu[1] = hwhm;
	    fcn[1] = 1.f;
	}

	if (*islit == 2) {

/* ****          t r i a n g u l a r   f u n c t i o n */

	    nfcn = 2;
	    dnu = hwhm * 2.f;
	    dnui = 1. / dnu;
	    delnu[0] = 0.f;
	    fcn[0] = 1.f;
	    dfdnu[0] = (real) (-1.f / dnu);
	    delnu[1] = hwhm * 2.f;
	    fcn[1] = 0.f;
	}

/* ****    determine if output grid is wide enough to */
/*        include weighting function. */

	nwnwt = (integer) (delnu[nfcn - 1] * 2.f / *dwn);

	if (nwnwt > 21) {
	    s_wsfe(&io___30);
	    do_fio(&c__1, "Number of output grid points needed to resolve ", (
		    ftnlen)47);
	    do_fio(&c__1, "weighting function exceeds dimension bound (nv0).",
		     (ftnlen)49);
	    do_fio(&c__1, "Use a smaller number of points to resolve ", (
		    ftnlen)42);
	    do_fio(&c__1, "weighting function.", (ftnlen)19);
	    e_wsfe();
	    s_wsfe(&io___31);
	    do_fio(&c__1, "nwnwt, nfcn, wnz, delnu, dwn", (ftnlen)28);
	    do_fio(&c__1, (char *)&nwnwt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nfcn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*wnz), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&delnu[nfcn - 1], (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&(*dwn), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    s_stop("", (ftnlen)0);
	}
    }

/* *****     s p e c t r a l    c o n v o l u t i o n   l o o p */

/* ****    find the number of output intervals needed between */
/*        wn0 and wnm1 */

    ndwn = (integer) (abs(dwn_dist__) / *dwn) + 1;

/* ****    if more than one spectral interval is needed, find */
/*        the spectral gradient */

    if (ndwn > 1) {
	i__1 = *nspt;
	for (nc = 1; nc <= i__1; ++nc) {
	    dzdwn[nc - 1] = (z0[nc + nz * 5162 - 5163] - zm1[nc + nz * 5162 - 
		    5163]) / dwn_dist__;
/* L2001: */
	}

    } else {
	i__1 = *nspt;
	for (nc = 1; nc <= i__1; ++nc) {
	    dzdwn[nc - 1] = 0.f;
/* L2021: */
	}
    }

/* ****    enter loop over output spectral intervals */

    i__1 = ndwn;
    for (nw = 1; nw <= i__1; ++nw) {

	dwnz = (real) ((real) (ndwn - nw) * *dwn);
	wn1[nz - 1] = wn0[nz - 1] - dwnz;
	i__2 = *nspt;
	for (nc = 1; nc <= i__2; ++nc) {
	    spect[nc - 1] = (real) (z0[nc + nz * 5162 - 5163] - dzdwn[nc - 1] 
		    * dwnz);
/* L2101: */
	}

	if (wn1[nz - 1] >= *wnmax || *iquit == 1) {
	    goto L3002;
	}

	if (points[nz] <= 1.f) {

/* ****        if the input spectral grid does not extend beyond end of */
/*            output grid, set wn1 to wnmin - delnu.  This is a fudge */
/*              to fix a computational problem. */

	    if (wn1[nz - 1] > *wnmin - delnu[nfcn - 1]) {
		wn1[nz - 1] = *wnmin - delnu[nfcn - 1];
	    }
	}

/* ****       find min and max wavenumber of the output grid where this */
/*           input value contributes. */

	wnmin0[nz - 1] = wn1[nz - 1] - delnu[nfcn - 1];
	wnmax0[nz - 1] = wn1[nz - 1] + delnu[nfcn - 1];

/* ****      find index of the output grid that corresponds to wnmin0. */
/*          (note: This approach assumes an equally-spaced output grid.) */

	i0 = (integer) ((wnmin0[nz - 1] - wni0[nz - 1]) / *dwn);

/* ****        determine if this interval is beyond shortwave limit of */
/*            output wrap-around buffer. */

	if (i0 < 0) {
	    i0 = 0;
	} else {

/* ****      move the initial point in the wrap-around output buffer by */
/*          i0 points.  write any completed output intervals to disk. */

	    if (nout[nz - 1] == 0) {
		in01 = in0[nz - 1];
		nout[nz - 1] = 1;
	    } else {
		in01 = in0[nz - 1] + 1;
	    }

	    i0new[nz - 1] = in0[nz - 1] + i0;
	    if (i0new[nz - 1] <= 21) {
		if (i0 > 1) {
		    i__2 = i0new[nz - 1];
		    for (n = in01; n <= i__2; ++n) {
			if (w[n + nz * 21 - 22] != 0.f) {
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				spectc[j - 1] = (real) (a[j + (n + nz * 21) * 
					5162 - 113565] / w[n + nz * 21 - 22]);
				a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2201: */
			    }
			    w[n + nz * 21 - 22] = 0.f;
			} else {
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				spectc[j - 1] = 0.f;
				a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2221: */
			    }
			}
			wl = 1e4f / wn_io__[n + nz * 21 - 22];
			wni_o__[nz - 1] = wn_io__[n + nz * 21 - 22];
			if (*iord == 1) {
			    if (*ifrm == 1) {
				io___46.ciunit = *iuout;
				s_wsfe(&io___46);
				do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 
					22], (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(
					real));
				i__3 = *nspt;
				for (j = 1; j <= i__3; ++j) {
				    do_fio(&c__1, (char *)&spectc[j - 1], (
					    ftnlen)sizeof(real));
				}
				e_wsfe();
			    } else {
				io___47.ciunit = *iuout;
				s_wsue(&io___47);
				do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 
					22], (ftnlen)sizeof(real));
				do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(
					real));
				i__3 = *nspt;
				for (j = 1; j <= i__3; ++j) {
				    do_uio(&c__1, (char *)&spectc[j - 1], (
					    ftnlen)sizeof(real));
				}
				e_wsue();
			    }
			} else {
			    if (*ifrm == 1) {
				io___48.ciunit = *iuout;
				s_wsfe(&io___48);
				do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 
					22], (ftnlen)sizeof(real));
				i__3 = *nspt;
				for (j = 1; j <= i__3; ++j) {
				    do_fio(&c__1, (char *)&spectc[j - 1], (
					    ftnlen)sizeof(real));
				}
				e_wsfe();
			    } else {
				io___49.ciunit = *iuout;
				s_wsue(&io___49);
				do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(
					real));
				do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 
					22], (ftnlen)sizeof(real));
				i__3 = *nspt;
				for (j = 1; j <= i__3; ++j) {
				    do_uio(&c__1, (char *)&spectc[j - 1], (
					    ftnlen)sizeof(real));
				}
				e_wsue();
			    }
			}
/* L2241: */
		    }

/* ****            reset the in0 pointer to in0 + i0 */

		    in0[nz - 1] = i0new[nz - 1];
		    wni0[nz - 1] = wn_io__[in0[nz - 1] + nz * 21 - 22];
		}

	    } else {

/* ****          write all values between in0 and nv0 to disk.  wrap the */
/*              starting interval around the end of the output buffer */
/*              and write values between 1 and i0new - nv0 to disk. */

		for (n = in01; n <= 21; ++n) {
		    if (w[n + nz * 21 - 22] != 0.f) {
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    spectc[j - 1] = (real) (a[j + (n + nz * 21) * 
				    5162 - 113565] / w[n + nz * 21 - 22]);
			    a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2301: */
			}
			w[n + nz * 21 - 22] = 0.f;
		    } else {
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    spectc[j - 1] = 0.f;
			    a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2321: */
			}
		    }
		    wl = 1e4 / wn_io__[n + nz * 21 - 22];
		    wni_o__[nz - 1] = wn_io__[n + nz * 21 - 22];
		    if (*iord == 1) {
			if (*ifrm == 1) {
			    io___50.ciunit = *iuout;
			    s_wsfe(&io___50);
			    do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    i__2 = *nspt;
			    for (j = 1; j <= i__2; ++j) {
				do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsfe();
			} else {
			    io___51.ciunit = *iuout;
			    s_wsue(&io___51);
			    do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    i__2 = *nspt;
			    for (j = 1; j <= i__2; ++j) {
				do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsue();
			}
		    } else {
			if (*ifrm == 1) {
			    io___52.ciunit = *iuout;
			    s_wsfe(&io___52);
			    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    i__2 = *nspt;
			    for (j = 1; j <= i__2; ++j) {
				do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsfe();
			} else {
			    io___53.ciunit = *iuout;
			    s_wsue(&io___53);
			    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    i__2 = *nspt;
			    for (j = 1; j <= i__2; ++j) {
				do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsue();
			}
		    }

/* L2341: */
		}

/* ****            write values between beginning of the buffer and in0 */

		in0[nz - 1] = i0new[nz - 1] - 21;
		i__2 = in0[nz - 1];
		for (nn = 1; nn <= i__2; ++nn) {
		    if (nn <= 21) {
			n = nn;
			if (w[n + nz * 21 - 22] != 0.f) {
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				spectc[j - 1] = (real) (a[j + (n + nz * 21) * 
					5162 - 113565] / w[n + nz * 21 - 22]);
				a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2401: */
			    }
			    w[n + nz * 21 - 22] = 0.f;
			} else {
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				spectc[j - 1] = 0.f;
				a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2421: */
			    }
			}

		    } else {

			n = 1;
			wn_io__[n + nz * 21 - 22] = (real) (*wnmin + delwn[nz 
				- 1]);
			delwn[nz - 1] += *dwn;
			i__3 = *nspt;
			for (j = 1; j <= i__3; ++j) {
			    spectc[j - 1] = 0.f;
			    a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2461: */
			}
		    }

		    wl = 1e4f / wn_io__[n + nz * 21 - 22];
		    wni_o__[nz - 1] = wn_io__[n + nz * 21 - 22];
		    if (*iord == 1) {
			if (*ifrm == 1) {
			    io___55.ciunit = *iuout;
			    s_wsfe(&io___55);
			    do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsfe();
			} else {
			    io___56.ciunit = *iuout;
			    s_wsue(&io___56);
			    do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsue();
			}
		    } else {
			if (*ifrm == 1) {
			    io___57.ciunit = *iuout;
			    s_wsfe(&io___57);
			    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsfe();
			} else {
			    io___58.ciunit = *iuout;
			    s_wsue(&io___58);
			    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			    do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], 
				    (ftnlen)sizeof(real));
			    i__3 = *nspt;
			    for (j = 1; j <= i__3; ++j) {
				do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
					sizeof(real));
			    }
			    e_wsue();
			}
		    }

/* L2481: */
		}

/* ****           reset end-of-buffer counters */

		if (in0[nz - 1] <= 21) {
		    in1[nz - 1] = in2[nz - 1];
		} else {
		    in0[nz - 1] = 1;
		    in1[nz - 1] = 0;
		}
		in2[nz - 1] = 0;

		wni0[nz - 1] = wn_io__[in0[nz - 1] + nz * 21 - 22];

	    }
	}

/* ****       find the index of the maximum output wavenumber bin affected */
/*           by this input quantity. */

	i1 = (integer) ((wnmax0[nz - 1] - wni0[nz - 1]) / *dwn);

/* ****       determine if interval is beyond end of wrap-around buffer. */

	i1new[nz - 1] = in0[nz - 1] + i1;
	if (i1new[nz - 1] <= 21) {

/* ****      create new output wavenumbers and initialize output arrays. */

	    i__2 = i1new[nz - 1];
	    for (n = in1[nz - 1] + 1; n <= i__2; ++n) {
		wn_io__[n + nz * 21 - 22] = (real) (*wnmin + delwn[nz - 1]);
		delwn[nz - 1] += *dwn;
		w[n + nz * 21 - 22] = 0.f;
		i__3 = *nspt;
		for (j = 1; j <= i__3; ++j) {
		    a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2501: */
		}
/* L2521: */
	    }

	    if (i1new[nz - 1] > in1[nz - 1]) {
		in1[nz - 1] = i1new[nz - 1];
	    }
	    in2[nz - 1] = 0;

	} else {

/* ****         initialize wavenumbers and spectral quantities between */
/*             current in1 and end of wrap-around buffer. */

	    for (n = in1[nz - 1] + 1; n <= 21; ++n) {
		wn_io__[n + nz * 21 - 22] = (real) (*wnmin + delwn[nz - 1]);
		delwn[nz - 1] += *dwn;
		w[n + nz * 21 - 22] = 0.f;
		i__2 = *nspt;
		for (j = 1; j <= i__2; ++j) {
		    a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2601: */
		}
/* L2621: */
	    }

/* ****         define new termination index (in2) and initialize values */
/*             between the beginning of the buffer and this point. */

	    in1[nz - 1] = 21;
	    i2new[nz - 1] = i1new[nz - 1] - 21;
	    if (i2new[nz - 1] > 21) {
		i2new[nz - 1] = 21;
	    }
	    i__2 = i2new[nz - 1];
	    for (n = in2[nz - 1] + 1; n <= i__2; ++n) {
		wn_io__[n + nz * 21 - 22] = (real) (*wnmin + delwn[nz - 1]);
		delwn[nz - 1] += *dwn;
		w[n + nz * 21 - 22] = 0.f;
		i__3 = *nspt;
		for (j = 1; j <= i__3; ++j) {
		    a[j + (n + nz * 21) * 5162 - 113565] = 0.f;
/* L2641: */
		}
/* L2661: */
	    }
	    in2[nz - 1] = i2new[nz - 1];
	    if (in2[nz - 1] >= in0[nz - 1]) {
		in2[nz - 1] = in0[nz - 1];
	    }

	}

/* ****      add distance-weighted average contribution to each output bin */

	if (nout[nz - 1] == 0) {
	    in01 = in0[nz - 1];
	    nout[nz - 1] = 1;
	} else {
	    in01 = in0[nz - 1] + 1;
	}

	i__2 = in1[nz - 1];
	for (n = in01; n <= i__2; ++n) {
	    dist = (d__1 = wn1[nz - 1] - wn_io__[n + nz * 21 - 22], abs(d__1))
		    ;
	    in = (integer) (dnui * dist) + 1;
	    if (in < 2) {
		wt = (fcn[in - 1] + dfdnu[in - 1] * (dist - delnu[in - 1])) * 
			dnuwt;
		if (wt != 0.f) {
		    w[n + nz * 21 - 22] += wt;
		    i__3 = *nspt;
		    for (j = 1; j <= i__3; ++j) {
			a[j + (n + nz * 21) * 5162 - 113565] += wt * spect[j 
				- 1];
/* L2701: */
		    }
		}
	    }
/* L2721: */
	}

	i__2 = in2[nz - 1];
	for (n = 1; n <= i__2; ++n) {
	    dist = (d__1 = wn1[nz - 1] - wn_io__[n + nz * 21 - 22], abs(d__1))
		    ;
	    in = (integer) (dnui * dist) + 1;
	    if (in < 2) {
		wt = (fcn[in - 1] + dfdnu[in - 1] * (dist - delnu[in - 1])) * 
			dnuwt;
		if (wt != 0.f) {
		    w[n + nz * 21 - 22] += wt;
		    i__3 = *nspt;
		    for (j = 1; j <= i__3; ++j) {
			a[j + (n + nz * 21) * 5162 - 113565] += wt * spect[j 
				- 1];
/* L2741: */
		    }
		}
	    }
/* L2761: */
	}

	wnm1[nz - 1] = wn0[nz - 1];

/* L2801: */
    }

/* ****   reset old wavenumber and spectral variables */

    wnm1[nz - 1] = wn0[nz - 1];
    wn0[nz - 1] = *wnz;
    i__1 = *nspt;
    for (nc = 1; nc <= i__1; ++nc) {
	zm1[nc + nz * 5162 - 5163] = z0[nc + nz * 5162 - 5163];
	z0[nc + nz * 5162 - 5163] = z__[nc];
/* L2901: */
    }

/* ****      g e t    t h e    n e x t    i n p u t    r e c o r d : */

    if (wn1[nz - 1] <= *wnmax && *iquit == 0) {
	return 0;
    }

/* ****      flush remaining buffers and quit */

L3002:

/*      write(*,'(1x,1a,i5,1a,i5,1a,1pe14.6)') */
/*     - 'rad_slit: flushing buffers and closing unit: ', */
/*     - iuout,' for variable: ',nz0,' at wavenumber: ',wnz */

    i__1 = in1[nz - 1];
    for (n = in0[nz - 1] + 1; n <= i__1; ++n) {
	if (wn_io__[n + nz * 21 - 22] <= *wnmax && wn_io__[n + nz * 21 - 22] 
		> wni_o__[nz - 1]) {
	    if (w[n + nz * 21 - 22] != 0.f) {
		i__2 = *nspt;
		for (j = 1; j <= i__2; ++j) {
		    spectc[j - 1] = (real) (a[j + (n + nz * 21) * 5162 - 
			    113565] / w[n + nz * 21 - 22]);
/* L3201: */
		}
		wni_o__[nz - 1] = wn_io__[n + nz * 21 - 22];
		wl = 1e4 / wn_io__[n + nz * 21 - 22];
		if (*iord == 1) {
		    if (*ifrm == 1) {
			io___65.ciunit = *iuout;
			s_wsfe(&io___65);
			do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsfe();
		    } else {
			io___66.ciunit = *iuout;
			s_wsue(&io___66);
			do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsue();
		    }
		} else {
		    if (*ifrm == 1) {
			io___67.ciunit = *iuout;
			s_wsfe(&io___67);
			do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsfe();
		    } else {
			io___68.ciunit = *iuout;
			s_wsue(&io___68);
			do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsue();
		    }
		}
	    }
	}

/* L3221: */
    }

    i__1 = in2[nz - 1];
    for (n = 1; n <= i__1; ++n) {
	if (wn_io__[n + nz * 21 - 22] <= *wnmax && wn_io__[n + nz * 21 - 22] 
		> wni_o__[nz - 1]) {
	    if (w[n + nz * 21 - 22] != 0.f) {
		i__2 = *nspt;
		for (j = 1; j <= i__2; ++j) {
		    spectc[j - 1] = (real) (a[j + (n + nz * 21) * 5162 - 
			    113565] / w[n + nz * 21 - 22]);
/* L3241: */
		}
		wl = 1e4 / wn_io__[n + nz * 21 - 22];
		wni_o__[nz - 1] = wn_io__[n + nz * 21 - 22];
		if (*iord == 1) {
		    if (*ifrm == 1) {
			io___69.ciunit = *iuout;
			s_wsfe(&io___69);
			do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsfe();
		    } else {
			io___70.ciunit = *iuout;
			s_wsue(&io___70);
			do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsue();
		    }
		} else {
		    if (*ifrm == 1) {
			io___71.ciunit = *iuout;
			s_wsfe(&io___71);
			do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_fio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsfe();
		    } else {
			io___72.ciunit = *iuout;
			s_wsue(&io___72);
			do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
			do_uio(&c__1, (char *)&wn_io__[n + nz * 21 - 22], (
				ftnlen)sizeof(real));
			i__2 = *nspt;
			for (j = 1; j <= i__2; ++j) {
			    do_uio(&c__1, (char *)&spectc[j - 1], (ftnlen)
				    sizeof(real));
			}
			e_wsue();
		    }
		}
	    }
	}

/* L3261: */
    }

    *iquit = -1;
    cl__1.cerr = 0;
    cl__1.cunit = *iuout;
    cl__1.csta = 0;
    f_clos(&cl__1);

    return 0;
} /* rad_slit__ */

