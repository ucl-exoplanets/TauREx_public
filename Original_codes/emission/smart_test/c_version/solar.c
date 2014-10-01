/* solar.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int solar_(integer *iusol, integer *ne, integer *nsivar, 
	integer *iend, doublereal *wnmax, doublereal *wn, doublereal *wnvar, 
	real *sol0, real *dsoldv)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void)
	    ;

    /* Local variables */
    static integer ni;

    /* Fortran I/O blocks */
    static cilist io___2 = { 1, 0, 1, 0, 0 };
    static cilist io___3 = { 0, 6, 0, "(/,1a,1pe12.4)", 0 };
    static cilist io___4 = { 0, 6, 0, "(/,2a,1pe12.4)", 0 };



/* cccccccccccccccccccccccccccc  s o l a r   ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine interpolates solar fluxes  to the standard      cc */
/* c    standard wavenumber grid used by line-by-line programs.         cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      iusol - unit number of file with solar fluxes                 cc */
/* c      wnmax - maximum wavenumber                                    cc */
/* c     nsivar - spectral interval initialization flag                 cc */
/* c      wnvar - wavenumber at which solar fluxes are specified        cc */
/* c         wn - current wavenumber                                    cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      wnvar - wavenumber at which solar fluxes are specified        cc */
/* c       sol0 - solar flux at wavenumber wnvar                        cc */
/* c     dsoldv - derivative of solar flux with respect to wavenumber   cc */
/* c       iend - end of file flag                                      cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccccc  s o l a r   ccccccccccccccccccccccccccccc */




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
    --sol0;
    --wnvar;
    --iend;
    --nsivar;

    /* Function Body */
    if (nsivar[2] == 0) {

/* ****    initialize spectral interval counter for this aerosols mode */

	nsivar[1] = 1;
	nsivar[2] = 2;

/* ****     initialize the aerosol optical properties */

	wnvar[2] = 0.;
	sol0[2] = 0.f;

    }

/* ****    swap spectral interval counters */

L2001:
    ni = nsivar[1];
    nsivar[1] = nsivar[2];
    nsivar[2] = ni;

/* ****     read the next solar flux value */

    io___2.ciunit = *iusol;
    i__1 = s_rsue(&io___2);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&wnvar[ni], (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&sol0[ni], (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = e_rsue();
L100001:
    if (i__1 < 0) {
	goto L4001;
    }
    if (i__1 > 0) {
	goto L4201;
    }

/* ****     determine if this segment is in desired spectral window */

    if (wnvar[ni] <= *wn) {
	goto L2001;
    }

/* ****    find the derivative of each quantity wrt wavenumber */

    *dsoldv = (real) ((sol0[2] - sol0[1]) / (wnvar[2] - wnvar[1]));

    return 0;

L4001:
    s_wsfe(&io___3);
    do_fio(&c__1, "End of solar flux file at wavenumber", (ftnlen)36);
    do_fio(&c__1, (char *)&wnvar[nsivar[1]], (ftnlen)sizeof(doublereal));
    e_wsfe();

    wnvar[ni] = *wnmax * 1.0001;
    sol0[ni] = 0.f;
    iend[*ne] = 1;

    return 0;

L4201:
    s_wsfe(&io___4);
    do_fio(&c__1, "Error in solar flux file ", (ftnlen)25);
    do_fio(&c__1, "after wavenumber", (ftnlen)16);
    do_fio(&c__1, (char *)&wnvar[nsivar[1]], (ftnlen)sizeof(doublereal));
    e_wsfe();

    wnvar[ni] = *wnmax * 1.01;
    sol0[ni] = 0.f;
    iend[*ne] = 1;

    return 0;
} /* solar_ */

