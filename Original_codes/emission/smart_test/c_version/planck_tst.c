/* planck_tst.f -- translated by f2c (version 20100827).
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
static integer c__5 = 5;
static integer c__4 = 4;

/* Subroutine */ int planck_(integer *itype, integer *nlev, integer *iuin, 
	integer *iuout, real *w, real *t, doublereal *b)
{
    /* Initialized data */

    static doublereal h__ = 6.6262e-34;
    static doublereal bk = 1.3806e-23;
    static doublereal c__ = 2.998e8;
    static doublereal twopi = 6.283185;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double exp(doublereal);

    /* Local variables */
    static integer j;
    static doublereal c2, fv, ang, arg, c1fv, c2fv;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };



/* ccccccccccccccccccccccccccc  p l a n c k  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine evaluates the planck function at 'nlev'         cc */
/* c    temperatures at a specific wavelenght or wavenumber.            cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c   itype - calculation type: 1) radinace, 2) flux                   cc */
/* c    nlev - number of levels.                                        cc */
/* c    iuin - input frequency/wavelength unit:                         cc */
/* c           1: wavenumber (cm**-1)                                   cc */
/* c           2: wavelength (microns)                                  cc */
/* c           3: wavelength (nanometers)                               cc */
/* c           4: wavelength (Angstroms)                                cc */
/* c   iuout - output units desired for planck function:                cc */
/* c           1: Watts/m**2/cm**-1                                     cc */
/* c           2: Watts/m**2/micron                                     cc */
/* c           3: Watts/m**2/nanometer                                  cc */
/* c           4: Watts/m**2/Angstrom                                   cc */
/* c           5: Watts/m**2/Hz                                         cc */
/* c       w - wavenumber (cm**-1) or wavelength                        cc */
/* c       t - temperature at each level                                cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c       b - planck function (units specified by iuout)               cc */
/* c    c1fv - numerator of planck function in desired units            cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  p l a n c k cccccccccccccccccccccccccccccc */





/* ****   define Planks constant, h (J s), Boltzmann's constant, k (J/K) */
/*       and the speed of light, c (m/sec). */

    /* Parameter adjustments */
    --b;
    --t;

    /* Function Body */


/* ****  evaluate the 2nd radiation constant and convert */
/*      input wavenumber or wavelength to Hz */

    s_wsle(&io___5);
    do_lio(&c__9, &c__1, "in planck: iuin,h, bk", (ftnlen)21);
    do_lio(&c__3, &c__1, (char *)&(*iuin), (ftnlen)sizeof(integer));
    do_lio(&c__5, &c__1, (char *)&h__, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&bk, (ftnlen)sizeof(doublereal));
    e_wsle();
    c2 = h__ / bk;
    fv = c__ * 100. * *w;
    c2fv = c2 * fv;

/* ****   specify the output units of the Planck function. */

    if (*itype == 1) {
	ang = 2.;
    } else {
	ang = twopi;
    }
    c1fv = ang * 100. * h__ * fv * fv * fv / c__;
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "c2fv, c1fv", (ftnlen)10);
    do_lio(&c__5, &c__1, (char *)&c2fv, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&c1fv, (ftnlen)sizeof(doublereal));
    e_wsle();

/* compute the Planck function based on the temperature profile. */

    i__1 = *nlev;
    for (j = 1; j <= i__1; ++j) {
	arg = c2fv / t[j];
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, "j,c2fv,t(j),arg=", (ftnlen)16);
	do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_lio(&c__5, &c__1, (char *)&c2fv, (ftnlen)sizeof(doublereal));
	do_lio(&c__4, &c__1, (char *)&t[j], (ftnlen)sizeof(real));
	do_lio(&c__5, &c__1, (char *)&arg, (ftnlen)sizeof(doublereal));
	e_wsle();
	if (arg > 0.f && arg < 200.) {
	    b[j] = c1fv / (exp(arg) - 1.);
	} else {

/* ****        use the Rayleigh-Jeans limit */

	    b[j] = c1fv * exp(-arg);

	}
/* L4201: */
    }

    return 0;
} /* planck_ */

