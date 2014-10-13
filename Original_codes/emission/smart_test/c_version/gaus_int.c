/* gaus_int.f -- translated by f2c (version 20100827).
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

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_rsle(cilist *), e_rsle(void), s_wsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double cos(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer n;
    static real pi, gmu[16], gwt[16];
    static integer nstr;
    static real gsum1;
    extern /* Subroutine */ int qgausn_(integer *, real *, real *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static cilist io___3 = { 0, 5, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, "(/,/,1a,/,1a)", 0 };
    static cilist io___10 = { 0, 6, 0, "(i5,2e12.4)", 0 };
    static cilist io___11 = { 0, 6, 0, "(/,/,1a,1pe12.4)", 0 };



/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                     cc */
/* c    p u r p o s e :                                                  cc */
/* c                                                                     cc */
/* c    this program uses gaussian quadrature schemes to evaulate        cc */
/* c    integrals of sample functions.                                   cc */
/* c                                                                     cc */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */



    pi = 3.14159f;

    s_wsle(&io___2);
    do_lio(&c__9, &c__1, "enter number of gaussian points: ", (ftnlen)33);
    e_wsle();
    s_rsle(&io___3);
    do_lio(&c__3, &c__1, (char *)&nstr, (ftnlen)sizeof(integer));
    e_rsle();

/* ****   call standard gaussian quadrature routine: */

    qgausn_(&nstr, gmu, gwt);

    s_wsfe(&io___7);
    do_fio(&c__1, " qgausn output ", (ftnlen)15);
    do_fio(&c__1, "  n     gmu      gwt", (ftnlen)20);
    e_wsfe();
    gsum1 = 0.f;
    i__1 = nstr;
    for (n = 1; n <= i__1; ++n) {
/*          gmu(nstr-n+1) = -gmu(n) */
/*          gwt(nstr-n+1) = gwt(n) */
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&gmu[n - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&gwt[n - 1], (ftnlen)sizeof(real));
	e_wsfe();
	gsum1 += gwt[n - 1] * cos(pi * gmu[n - 1]);
/* L1001: */
    }

    s_wsfe(&io___11);
    do_fio(&c__1, "gsum1 =", (ftnlen)7);
    do_fio(&c__1, (char *)&gsum1, (ftnlen)sizeof(real));
    e_wsfe();

    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int gaus_int__ () { MAIN__ (); return 0; }
