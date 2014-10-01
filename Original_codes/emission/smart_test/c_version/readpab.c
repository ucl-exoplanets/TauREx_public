/* readpab.f -- translated by f2c (version 20100827).
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
static integer c__255 = 255;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__10000 = 10000;

/* Subroutine */ int readpab_(integer *iu, doublereal *wnmin, doublereal *
	wnmax, integer *npab0, integer *ntpab, real *tpab0, doublereal *
	wnpab0, real *pab0, integer *ioend, integer *ioerr, doublereal *wneof)
{
    /* System generated locals */
    integer i__1, i__2;
    cllist cl__1;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_rsle(cilist *), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, n, ic;
    static doublereal wn0;
    static integer nlb, len, ntb;
    static char str[255], str1[1*255];
    static integer itst;
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, "(1a)", 0 };
    static icilist io___9 = { 1, str, 0, 0, 255, 1 };
    static cilist io___11 = { 0, 0, 0, "(/,/)", 0 };
    static cilist io___12 = { 0, 6, 0, "(/,1a,i5,1a)", 0 };
    static cilist io___14 = { 1, 0, 1, 0, 0 };
    static cilist io___16 = { 0, 6, 0, "(/,1a,1pe14.6)", 0 };
    static cilist io___17 = { 0, 6, 0, "(/,1a,1pe14.6)", 0 };
    static cilist io___18 = { 0, 6, 0, "(/,1a,1pe14.6,/2a,i5)", 0 };



/* cccccccccccccccccccccccccc  r e a d p a b  cccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads wavelength-dependent pressure-induced     cc */
/* c    absorption coefficients of gases at one or more temperatures    cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        iu - unit number for absorption coefficeints.               cc */
/* c     wnpab - wavenumbers where absorption coefficients are given    cc */
/* c       pab - absorption coefficient at each wavenumber              cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      npab0 - number of wavenumbers where absorption coefficients   cc */
/* c             are specified                                          cc */
/* c     wnpab - wavenumbers where absorption coefficients are given    cc */
/* c       pab - absorption coefficient at each wavenumber              cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  r e a d p a b  cccccccccccccccccccccccccccc */




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






/* ****    read the number of temperatures */

    /* Parameter adjustments */
    --wneof;
    pab0 -= 71;
    --wnpab0;
    --tpab0;

    /* Function Body */
    ic = 0;
L1021:
    if (ic < 100) {
	++ic;

/* ****     skip over any ascii header at the top of the file, and */
/*         find the first record with numbers in it */

	io___2.ciunit = *iu;
	s_rsfe(&io___2);
	do_fio(&c__1, str, (ftnlen)255);
	e_rsfe();

/* *****   find the first non-blank character */

	charsp_(str, str1, &len, &c__255, &nlb, &ntb, (ftnlen)255, (ftnlen)1);

	itst = *(unsigned char *)&str1[0];
	if (itst >= 48 && itst <= 57 && len > 0) {

/*         this record includes numerical values */

	    i__1 = s_rsli(&io___9);
	    if (i__1 != 0) {
		goto L1021;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*ntpab), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L1021;
	    }
	    i__2 = *ntpab;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_lio(&c__4, &c__1, (char *)&tpab0[i__], (ftnlen)
			sizeof(real));
		if (i__1 != 0) {
		    goto L1021;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1021;
	    }
	    io___11.ciunit = *iu;
	    s_rsfe(&io___11);
	    e_rsfe();
	} else {
	    goto L1021;
	}
    } else {
	s_wsfe(&io___12);
	do_fio(&c__1, " Error reading unit", (ftnlen)19);
	do_fio(&c__1, (char *)&(*iu), (ftnlen)sizeof(integer));
	do_fio(&c__1, " in readpab", (ftnlen)11);
	e_wsfe();
	s_stop("", (ftnlen)0);
    }

    n = 0;
L2001:
    io___14.ciunit = *iu;
    i__1 = s_rsle(&io___14);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_lio(&c__5, &c__1, (char *)&wn0, (ftnlen)sizeof(doublereal));
    if (i__1 != 0) {
	goto L100001;
    }
    i__2 = *ntpab;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_lio(&c__4, &c__1, (char *)&pab0[i__ + (n + 1) * 70], (
		ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L100001;
	}
    }
    i__1 = e_rsle();
L100001:
    if (i__1 < 0) {
	goto L2021;
    }
    if (i__1 > 0) {
	goto L2041;
    }

    wnpab0[n + 1] = wn0;
    if (wnpab0[n + 1] >= *wnmin) {
	++n;
    }
    if (n > 10000) {
	goto L2061;
    }

    if (wn0 <= *wnmax) {
	goto L2001;
    }

L2021:
    cl__1.cerr = 0;
    cl__1.cunit = *iu;
    cl__1.csta = 0;
    f_clos(&cl__1);
    *npab0 = n;

    *ioend = 1;
    *ioerr = 0;
    wneof[1] = wn0;
    s_wsfe(&io___16);
    do_fio(&c__1, "EOF in readpab after wn = ", (ftnlen)26);
    do_fio(&c__1, (char *)&wnpab0[n], (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;

L2041:
    *npab0 = n;
    *ioend = 0;
    *ioerr = 1;
    wneof[1] = wn0;
    s_wsfe(&io___17);
    do_fio(&c__1, "Error in readpab after wn = ", (ftnlen)28);
    do_fio(&c__1, (char *)&wnpab0[n], (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;

L2061:
    *npab0 = n;
    *ioend = 0;
    *ioerr = 1;
    wneof[1] = wn0;
    s_wsfe(&io___18);
    do_fio(&c__1, "Error in readpab after wn = ", (ftnlen)28);
    do_fio(&c__1, (char *)&wnpab0[n], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, "Number of elements in array exceeds dimension bound:", (
	    ftnlen)52);
    do_fio(&c__1, " nsp = ", (ftnlen)7);
    do_fio(&c__1, (char *)&c__10000, (ftnlen)sizeof(integer));
    e_wsfe();

    return 0;
} /* readpab_ */

