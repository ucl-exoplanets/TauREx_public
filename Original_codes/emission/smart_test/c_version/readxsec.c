/* readxsec.f -- translated by f2c (version 20100827).
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

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__5 = 5;

/* Subroutine */ int readxsec_(integer *iu, integer *nxsec, real *wgtgs, 
	doublereal *wnmin, doublereal *wnmax, doublereal *wl, real *xsec, 
	integer *nxmx, integer *iquit)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsle(cilist *), e_rsle(void), do_lio(integer *, integer *, char 
	    *, ftnlen), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), 
	    e_wsfe(void);

    /* Local variables */
    static integer n;
    static doublereal wl0;
    static real xsec0;
    static integer nskip;
    static doublereal wlmin, wlmax;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, 0, 0 };
    static cilist io___6 = { 0, 0, 0, 0, 0 };
    static cilist io___7 = { 0, 0, 0, 0, 0 };
    static cilist io___8 = { 1, 0, 1, 0, 0 };
    static cilist io___11 = { 0, 6, 0, "(/,1a,1pe14.6,/2a,i5)", 0 };
    static cilist io___12 = { 0, 6, 0, "(/,1a,1pe14.6)", 0 };



/* ccccccccccccccccccccccccc  r e a d x s e c  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads wavelength-dependent gas absorption       cc */
/* c    cross-sections.  This version of the code                       cc */
/* c    assumes no pressure or temperature dependence for these         cc */
/* c    coefficients.                                                   cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        iu - unit number for absorption coefficeints.               cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     nxsec - number of wavenumbers where absorption coefficients    cc */
/* c             are specified                                          cc */
/* c        wl - wavenumbers where absorption coefficients are given    cc */
/* c      xsec - absorption coefficient at each wavenumber              cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  r e a d x s e c  ccccccccccccccccccccccccccc */




/* ****    skip over the file header */

    /* Parameter adjustments */
    --xsec;
    --wl;

    /* Function Body */
    wlmin = 1e4 / *wnmax;
    wlmax = 1e4 / *wnmin;

    nskip = 2;
    i__1 = nskip;
    for (n = 1; n <= i__1; ++n) {
	io___5.ciunit = *iu;
	s_rsle(&io___5);
	e_rsle();
/* L1001: */
    }

/* ****    read the gas molecular weight (kg/kmole) */

    io___6.ciunit = *iu;
    s_rsle(&io___6);
    do_lio(&c__4, &c__1, (char *)&(*wgtgs), (ftnlen)sizeof(real));
    e_rsle();

    nskip = 4;
    i__1 = nskip;
    for (n = 1; n <= i__1; ++n) {
	io___7.ciunit = *iu;
	s_rsle(&io___7);
	e_rsle();
/* L1021: */
    }

/* ****   read the gas absorption cross-section (cm**2) */

    n = 0;
L2001:
    if (n < *nxmx && *iquit == 0) {
	io___8.ciunit = *iu;
	i__1 = s_rsle(&io___8);
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&wl0, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&xsec0, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = e_rsle();
L100001:
	if (i__1 < 0) {
	    goto L2022;
	}
	if (i__1 > 0) {
	    goto L2042;
	}

	if (wl0 >= wlmin && wl0 <= wlmax) {
	    ++n;
	    wl[n] = wl0;
	    xsec[n] = xsec0;
	} else {
	    if (n < 2) {
		n = 1;
		wl[n] = wl0;
		xsec[n] = xsec0;
	    } else {
		++n;
		wl[n] = wl0;
		xsec[n] = xsec0;
		goto L2022;
	    }
	}

	goto L2001;
    } else {
	if (n > *nxmx) {
	    s_wsfe(&io___11);
	    do_fio(&c__1, "Error in readxsec after wl = ", (ftnlen)29);
	    do_fio(&c__1, (char *)&wl[n], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, "Number of elements in array exceeds dimension bou"
		    "nd:", (ftnlen)52);
	    do_fio(&c__1, " nxmx = ", (ftnlen)8);
	    do_fio(&c__1, (char *)&(*nxmx), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

L2022:
    *nxsec = n;

    *iquit = 1;
    return 0;

L2042:
    *iquit = -1;
    s_wsfe(&io___12);
    do_fio(&c__1, "Error in readxsec after wl = ", (ftnlen)29);
    do_fio(&c__1, (char *)&wl[n], (ftnlen)sizeof(doublereal));
    e_wsfe();

    return 0;
} /* readxsec_ */

