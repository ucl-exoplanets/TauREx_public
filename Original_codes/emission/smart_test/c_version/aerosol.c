/* aerosol.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int aerosol_(integer *ne, integer *mode, integer *iu, 
	doublereal *wn, doublereal *wnmin, doublereal *wnmax, real *aerqext, 
	real *aerqsca, real *aerg0, real *aerpmom, real *dpmomdv, real *
	dqextdv, real *dqscadv, real *dg0dv, integer *nmomaer, integer *
	nmom_mx__, integer *nsiext, doublereal *wnext, doublereal *wn_eof__, 
	integer *io_end__, integer *io_err__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer ni;
    static real dvi;
    static integer mom;

    /* Fortran I/O blocks */
    static cilist io___3 = { 1, 0, 1, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };



/* ccccccccccccccccccccccccc  a e r o s o l  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine interpolates aerosol optical properties to the  cc */
/* c    standard wavenumber grid used by line-by-line programs.         cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c         ne : index of radiatively active component                 cc */
/* c       mode : aerosol particle mode index                           cc */
/* c       nmom : number of legendre polynomial coefficients for the    cc */
/* c              scattering phase function expansion.                  cc */
/* c          iu: unit for aerosol optical property i/o                 cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    aerqext(ni,mode): aerosol extincition efficiency                cc */
/* c    aerqsca(ni,mode): aerosol scattering efficiency                 cc */
/* c      aerg0(ni,mode): scattering asymmetry parameter                cc */
/* c    dqextdv(mode): rate of change of qext with wavenumber           cc */
/* c    dqscadv(mode): rate of change of qsca with wavenumber           cc */
/* c      dg0dv(mode): rate of change of g with wavenumber              cc */
/* c    aerpmom(ni,mom,mode) : legendre polynomial momemt of aerosol    cc */
/* c                         scattering phase function                  cc */
/* c    dpmomdv(0,mode) : rate of change of pmom with wavenumber        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  a e r o s o l  ccccccccccccccccccccccccccccc */



/* ****   aerosol moment variables */


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



/* ****    variables for the output wavenumber grid */




/* ****   input aerosol wavelength-dependent variables */


    /* Parameter adjustments */
    --io_err__;
    --io_end__;
    wn_eof__ -= 3;
    wnext -= 3;
    nsiext -= 3;
    --nmom_mx__;
    nmomaer -= 3;
    --dg0dv;
    --dqscadv;
    --dqextdv;
    dpmomdv -= 201;
    aerpmom -= 403;
    aerg0 -= 3;
    aerqsca -= 3;
    aerqext -= 3;

    /* Function Body */
    if (nsiext[(*ne << 1) + 2] == 0) {
/*        write(*,*) 'intializing aerosol properties for mode: ',mode */

/* ****    initialize spectral interval counter for this aerosols mode */

	nsiext[(*ne << 1) + 1] = 1;
	nsiext[(*ne << 1) + 2] = 2;

/* ****     initialize the aerosol optical properties */

	wnext[(*ne << 1) + 2] = 0.;
	aerqext[(*mode << 1) + 2] = 0.f;
	aerqsca[(*mode << 1) + 2] = 0.f;
	aerg0[(*mode << 1) + 2] = 0.f;
	aerpmom[*mode * 402 + 2] = 1.f;
	nmomaer[(*mode << 1) + 2] = 1;
	for (mom = 1; mom <= 200; ++mom) {
	    aerpmom[(mom + *mode * 201 << 1) + 2] = 0.f;
/* L1201: */
	}

    }

/* ****    swap spectral interval counters */

L2001:
    ni = nsiext[(*ne << 1) + 1];
    nsiext[(*ne << 1) + 1] = nsiext[(*ne << 1) + 2];
    nsiext[(*ne << 1) + 2] = ni;

/* ****      read aerosol extinction efficiencies and phase functions */

    io___3.ciunit = *iu;
    i__1 = s_rsue(&io___3);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&nmomaer[ni + (*mode << 1)], (ftnlen)sizeof(
	    integer));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&wnext[ni + (*ne << 1)], (ftnlen)sizeof(
	    doublereal));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&aerqext[ni + (*mode << 1)], (ftnlen)sizeof(
	    real));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&aerqsca[ni + (*mode << 1)], (ftnlen)sizeof(
	    real));
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_uio(&c__1, (char *)&aerg0[ni + (*mode << 1)], (ftnlen)sizeof(
	    real));
    if (i__1 != 0) {
	goto L100001;
    }
    i__2 = nmomaer[ni + (*mode << 1)];
    for (mom = 0; mom <= i__2; ++mom) {
	i__1 = do_uio(&c__1, (char *)&aerpmom[ni + (mom + *mode * 201 << 1)], 
		(ftnlen)sizeof(real));
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

/* ****     determine if this segment is in desired spectral window */

    if (wnext[ni + (*ne << 1)] <= *wnmin) {
	wn_eof__[(*ne << 1) + 1] = wnext[ni + (*ne << 1)];
    }
    if (wnext[ni + (*ne << 1)] <= *wn) {
	goto L2001;
    }
    wn_eof__[(*ne << 1) + 2] = wnext[ni + (*ne << 1)];

/* ****    find the derivative of each quantity wrt wavenumber */

    dvi = (real) (1. / (wnext[(*ne << 1) + 2] - wnext[(*ne << 1) + 1]));
    dqextdv[*mode] = (aerqext[(*mode << 1) + 2] - aerqext[(*mode << 1) + 1]) *
	     dvi;
    dqscadv[*mode] = (aerqsca[(*mode << 1) + 2] - aerqsca[(*mode << 1) + 1]) *
	     dvi;
    dg0dv[*mode] = (aerg0[(*mode << 1) + 2] - aerg0[(*mode << 1) + 1]) * dvi;

/* ****   find the mean values of the phase function moments */
/*        - the number of phase function moments is the maximum */
/*          of the two at ni = 1,2 */

    nmom_mx__[*mode] = nmomaer[(*mode << 1) + 1];
    if (nmomaer[(*mode << 1) + 2] > nmom_mx__[*mode]) {
	nmom_mx__[*mode] = nmomaer[(*mode << 1) + 2];
    }

    i__1 = nmom_mx__[*mode];
    for (mom = 1; mom <= i__1; ++mom) {
	dpmomdv[mom + *mode * 201] = (aerpmom[(mom + *mode * 201 << 1) + 2] - 
		aerpmom[(mom + *mode * 201 << 1) + 1]) * dvi;
/* L3101: */
    }

    return 0;

L4001:
    s_wsle(&io___5);
    do_lio(&c__9, &c__1, "End of aerosol file for mode", (ftnlen)28);
    do_lio(&c__3, &c__1, (char *)&(*mode), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, "after wavenumber", (ftnlen)16);
    do_lio(&c__5, &c__1, (char *)&wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)], 
	    (ftnlen)sizeof(doublereal));
    e_wsle();

    io_end__[*ne] = 1;
    wnext[ni + (*ne << 1)] = *wnmax * 1.01;
    aerqext[ni + (*mode << 1)] = 0.f;
    aerqsca[ni + (*mode << 1)] = 0.f;
    aerg0[ni + (*mode << 1)] = 0.f;
    dqextdv[*mode] = 0.f;
    dqscadv[*mode] = 0.f;
    dg0dv[*mode] = 0.f;
    aerpmom[ni + *mode * 402] = 1.f;
    nmomaer[ni + (*mode << 1)] = 1;
    dpmomdv[*mode * 201] = 0.f;
    for (mom = 1; mom <= 200; ++mom) {
	aerpmom[ni + (mom + *mode * 201 << 1)] = 0.f;
	dpmomdv[mom + *mode * 201] = 0.f;
/* L4101: */
    }

    return 0;

L4201:
    s_wsle(&io___6);
    do_lio(&c__9, &c__1, "Error in aerosol file for mode", (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*mode), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, "after wavenumber", (ftnlen)16);
    do_lio(&c__5, &c__1, (char *)&wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)], 
	    (ftnlen)sizeof(doublereal));
    e_wsle();

    io_err__[*ne] = 1;
    wnext[ni + (*ne << 1)] = *wnmax * 1.01;
    aerqext[ni + (*mode << 1)] = 0.f;
    aerqsca[ni + (*mode << 1)] = 0.f;
    aerg0[ni + (*mode << 1)] = 0.f;
    dqextdv[*mode] = 0.f;
    dqscadv[*mode] = 0.f;
    dg0dv[*mode] = 0.f;
    aerpmom[ni + *mode * 402] = 1.f;
    nmomaer[ni + (*mode << 1)] = 1;
    dpmomdv[*mode * 201] = 0.f;
    for (mom = 1; mom <= 200; ++mom) {
	aerpmom[ni + (mom + *mode * 201 << 1)] = 0.f;
	dpmomdv[mom + *mode * 201] = 0.f;
/* L4221: */
    }

    return 0;
} /* aerosol_ */

