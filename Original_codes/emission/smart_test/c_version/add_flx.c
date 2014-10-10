/* add_flx.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int add_flx__(integer *nlyr, doublereal *tl, doublereal *rl, 
	doublereal *fup, doublereal *fdn, doublereal *urt, doublereal *drb, 
	doublereal *uft, doublereal *dfb, doublereal *ft, doublereal *fb, 
	doublereal *flxu, doublereal *flxd)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer k, l;
    static doublereal dd[70], rb[70], du[70], rt[70];
    static integer nlev;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___9 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___10 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___11 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___12 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___13 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___14 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___15 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___16 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___17 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___18 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___19 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___20 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___21 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };
    static cilist io___22 = { 0, 6, 0, "(1a,16(1pe12.4))", 0 };



/* cccccccccccccccccccccccccc  a d d _ f l x   ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this program uses the layer adding method for finding fluxes    cc */
/* c    in a vertically inhomogeneous scattering, absorbing atmosphere. cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       nlyr - number of computational model layers                  cc */
/* c         rl - diffuse flux reflectance of each layer                cc */
/* c         tl - diffuse flux transmittance of each homogeneous layer  cc */
/* c        fup - f+, upward flux at the top of each homogeneous layer  cc */
/* c        fdn - f-, downward flux at the base of each layer           cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      flxu - upward flux at nlyr+1 layer boundries.                 cc */
/* c             (flxu(l) refers to the upward flux at the top          cc */
/* c              of layer l)                                           cc */
/* c      flxd - downward flux at nlyr+1 layer boundries.               cc */
/* c             (flxd(l) refers to the downward flux at the bottom     cc */
/* c              of layer l-1)                                         cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  a d d _ f l x   ccccccccccccccccccccccccccc */




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



/* ****   layer flux transmittances, absorptances, and fluxes */


/* ****   variables for flux layer adding method */


/* ****    quantitities used on radiance adding */


/* ****   output */


    /* Parameter adjustments */
    --flxd;
    --flxu;
    --fb;
    --ft;
    --dfb;
    --uft;
    --drb;
    --urt;
    --fdn;
    --fup;
    --rl;
    --tl;

    /* Function Body */
    nlev = *nlyr + 1;

/* ****  find layer-integrated reflectivities.  the quantities are: */
/*      rt(k): reflectance of an inhomogeneous layer extendind from */
/*             level 1 to level k */
/*      rb(k): reflectance of an inhomogeneous layer extendind from */
/*             the surface (nlev) to level k */

    rt[0] = 0.f;
    rb[nlev - 1] = rl[nlev];
    urt[nlev] = 0.;
    drb[1] = 0.;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	dd[k - 1] = 1. / (1. - rt[k - 1] * rl[k]);
	urt[k] = tl[k] * dd[k - 1];
/* Computing 2nd power */
	d__1 = tl[k];
	rt[k] = rl[k] + d__1 * d__1 * rt[k - 1] * dd[k - 1];
	du[nlev - k - 1] = 1. / (1. - rl[nlev - k] * rb[nlev - k]);
	drb[nlev - k + 1] = tl[nlev - k] * du[nlev - k - 1];
/* Computing 2nd power */
	d__1 = tl[nlev - k];
	rb[nlev - k - 1] = rl[nlev - k] + d__1 * d__1 * rb[nlev - k] * du[
		nlev - k - 1];
/* L1201: */
    }

/* ****   use adding method to find upward and downward fluxes */
/*       for combined layers.  start at top, adding homogeneous */
/*       layers to the base of the existing inhomogeneous layer */

    ft[1] = 0.;
    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	ft[k + 1] = fdn[k + 1] + tl[k] * (ft[k] + rt[k - 1] * fup[k]) * dd[k 
		- 1];
/* L1401: */
    }

/* ****   use adding method to find upward and downward fluxes */
/*       for combined layers.  start at bottom. */

    fb[nlev] = fup[nlev];
    i__1 = *nlyr;
    for (l = 1; l <= i__1; ++l) {
	k = nlev - l;
	fb[k] = fup[k] + tl[k] * (fb[k + 1] + rb[k] * fdn[k + 1]) * du[k - 1];
/* L1601: */
    }

/* ****  find the total upward and downward fluxes at interfaces */
/*      between inhomogeneous layers. */

    flxd[1] = fdn[1];
    i__1 = nlev;
    for (k = 2; k <= i__1; ++k) {
	flxd[k] = (ft[k] + rt[k - 1] * fb[k]) / (1. - rt[k - 1] * rb[k - 1]);
/* L2001: */
    }

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	flxu[k] = fb[k] + rb[k - 1] * flxd[k];
/* L2021: */
    }

/* ****   find the variables dft and ufb for radiance adding */

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	uft[k] = (fup[k] + rl[k] * ft[k]) * dd[k - 1];
	dfb[k + 1] = (fdn[k + 1] + rl[k] * fb[k + 1]) * du[k - 1];
/* L3001: */
    }
    s_wsfe(&io___8);
    do_fio(&c__1, "in add_flx: ", (ftnlen)12);
    e_wsfe();
    s_wsfe(&io___9);
    do_fio(&c__1, "rl:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&rl[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___10);
    do_fio(&c__1, "tl:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&tl[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___11);
    do_fio(&c__1, "dd:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&dd[k - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___12);
    do_fio(&c__1, "du:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&du[k - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___13);
    do_fio(&c__1, "rt:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&rt[k - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___14);
    do_fio(&c__1, "rb:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&rb[k - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___15);
    do_fio(&c__1, "ft:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&ft[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___16);
    do_fio(&c__1, "fb:    ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&fb[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___17);
    do_fio(&c__1, "fdn:   ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&fdn[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___18);
    do_fio(&c__1, "fup:   ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&fup[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___19);
    do_fio(&c__1, "uft:   ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&uft[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___20);
    do_fio(&c__1, "dfb:   ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&dfb[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___21);
    do_fio(&c__1, "flxd:  ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&flxd[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___22);
    do_fio(&c__1, "flxu:  ", (ftnlen)7);
    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	do_fio(&c__1, (char *)&flxu[k], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

    return 0;
} /* add_flx__ */

