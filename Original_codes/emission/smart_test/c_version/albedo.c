/* albedo.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int albedo_(integer *ne, integer *iusur, integer *iflext, 
	integer *nsiext, integer *iref, integer *nref, integer *nza, logical *
	lamber, real *umu0, doublereal *wn, doublereal *wnmin, doublereal *
	wnmax, doublereal *wnext, real *surf_0__, real *dsurfdv, real *ws, 
	real *phiw, real *surf_opt__, real *alb, doublereal *wn_eof__, 
	integer *io_end__, integer *io_err__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void)
	    ;

    /* Local variables */
    static integer n, ni;
    static doublereal dv;
    static integer nz;
    static doublereal dvi;
    extern doublereal dref_(real *, real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 1, 0, 1, 0, 0 };
    static cilist io___7 = { 0, 6, 0, "(/,1a,1pe14.6)", 0 };
    static cilist io___8 = { 0, 6, 0, "(/,1a,1pe12.4)", 0 };



/* ccccccccccccccccccccccccc   a l b e d o   cccccccccccccccccccccccccccc */
/* c                                                                   cc */
/* c    p u r p o s e :                                                cc */
/* c                                                                   cc */
/* c    this subroutine interpolates surface properties to the         cc */
/* c    common wavenumber grid used by line-by-line programs.          cc */
/* c                                                                   cc */
/* c     note: this version calculates the lamber albedo when a        cc */
/* c           non-lambertian albedo is selected.                      cc */
/* c                                                                   cc */
/* c    i n p u t :                                                    cc */
/* c                                                                   cc */
/* c       ne: optical property index counter (prior to albedo)        cc */
/* c   iflext: spectral flag: 1 => find next albedo point              cc */
/* c   nsiext: spectral index counter                                  cc */
/* c    iusur: input unit number for surface albedo file.              cc */
/* c     nref: number of wavelength dependent surface properties       cc */
/* c       wn: current wavenumber (cm-1)                               cc */
/* c    wnext: wavenumbers where albedos are specified (cm-1)          cc */
/* c    wnmax: maximum wavenumber (cm-1)                               cc */
/* c                                                                   cc */
/* c    o u t p u t :                                                  cc */
/* c                                                                   cc */
/* c   surf_0: surface properties at spectral point, wnext(ni,ne)      cc */
/* c  surf_opt: interpolated surface albedo at spectral point wn       cc */
/* c                                                                   cc */
/* ccccccccccccccccccccccccc   a l b e d o   cccccccccccccccccccccccccccc */





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



/* *****   albedo i/o flags and counters */


/* ****   output albedos and spectral gradients */




/* ****    variables for the output wavenumber grid */


/* ****   increment index counter and check new-value flag */

    /* Parameter adjustments */
    --io_err__;
    --io_end__;
    wn_eof__ -= 3;
    --alb;
    --surf_opt__;
    --dsurfdv;
    surf_0__ -= 3;
    wnext -= 3;
    --umu0;
    nsiext -= 3;
    --iflext;

    /* Function Body */
    ++(*ne);

/* ****   initialize values if this is the first interval */

    if (nsiext[(*ne << 1) + 2] == 0) {

/* ****      initialize spectral interval counter for this aerosols mode */

	nsiext[(*ne << 1) + 1] = 1;
	nsiext[(*ne << 1) + 2] = 2;

/* ****      initialize the aerosol optical properties */

	wnext[(*ne << 1) + 2] = 0.;
	i__1 = *nref;
	for (n = 1; n <= i__1; ++n) {
	    surf_0__[(n << 1) + 2] = 0.f;
/* L1001: */
	}

    }

    if (iflext[*ne] == 1) {

/* ****    swap spectral interval counters */

L2001:
	ni = nsiext[(*ne << 1) + 1];
	nsiext[(*ne << 1) + 1] = nsiext[(*ne << 1) + 2];
	nsiext[(*ne << 1) + 2] = ni;

/* ****     read the next surface optical property value for this gas */

	io___3.ciunit = *iusur;
	i__1 = s_rsue(&io___3);
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = do_uio(&c__1, (char *)&wnext[ni + (*ne << 1)], (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L100001;
	}
	i__2 = *nref;
	for (n = 1; n <= i__2; ++n) {
	    i__1 = do_uio(&c__1, (char *)&surf_0__[ni + (n << 1)], (ftnlen)
		    sizeof(real));
	    if (i__1 != 0) {
		goto L100001;
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

/* ****       determine if this segment is in desired spectral window */

	if (wnext[ni + (*ne << 1)] <= *wnmin) {
	    wn_eof__[(*ne << 1) + 1] = wnext[ni + (*ne << 1)];
	}

	if (wnext[ni + (*ne << 1)] <= *wn) {
	    goto L2001;
	}

	wn_eof__[(*ne << 1) + 2] = wnext[ni + (*ne << 1)];

/* ****    find the derivative of each quantity wrt wavenumber */

	dvi = 1. / (wnext[(*ne << 1) + 2] - wnext[(*ne << 1) + 1]);
	i__1 = *nref;
	for (n = 1; n <= i__1; ++n) {
	    dsurfdv[n] = (real) ((surf_0__[(n << 1) + 2] - surf_0__[(n << 1) 
		    + 1]) * dvi);
/* L2021: */
	}

	iflext[*ne] = 0;

    }

/* ****       interpolate  surface properties to this wn */

    dv = *wn - wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)];
    i__1 = *nref;
    for (n = 1; n <= i__1; ++n) {
	surf_opt__[n] = (real) (surf_0__[nsiext[(*ne << 1) + 1] + (n << 1)] + 
		dsurfdv[n] * dv);
/* L2041: */
    }

    if (*lamber) {

/* ****     a lambert albedo is used */

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb[nz] = surf_opt__[1];
/* L2061: */
	}

    } else {

	if (*iref == 4) {

/* ****      set wind speed and direction for Cox/Munk model */

	    surf_opt__[3] = *ws;
	    surf_opt__[4] = *phiw;
	}

/* ****    use the disort dref routine to find the flux albedo */
/*        for each solar zenith angle */

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    alb[nz] = dref_(&umu0[nz], &surf_opt__[1], iref);
/* L2081: */
	}

    }

    return 0;

/* ****   end of surface optical property file: */

L4001:
    wnext[ni + (*ne << 1)] = *wnmax;
    i__1 = *nref;
    for (n = 1; n <= i__1; ++n) {
	surf_0__[ni + (n << 1)] = 0.f;
	dsurfdv[n] = 0.f;
/* L4021: */
    }
    iflext[*ne] = 0;
    s_wsfe(&io___7);
    do_fio(&c__1, "End of surface albedo file at wavenumber ", (ftnlen)41);
    do_fio(&c__1, (char *)&(*wn), (ftnlen)sizeof(doublereal));
    e_wsfe();
    io_end__[*ne] = 1;

    return 0;

/* ****   error in surface optical property file */

L4201:
    s_wsfe(&io___8);
    do_fio(&c__1, "Error in surface albedo file after wavenumber", (ftnlen)45)
	    ;
    do_fio(&c__1, (char *)&wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)], (
	    ftnlen)sizeof(doublereal));
    e_wsfe();
    io_err__[*ne] = 1;

    return 0;
} /* albedo_ */

