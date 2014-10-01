/* planck.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int planck_(integer *itype, integer *nlev, integer *iuin, 
	integer *iunits, real *w, real *t, doublereal *b)
{
    /* Initialized data */

    static doublereal h__ = 6.6262e-34;
    static doublereal bk = 1.3806e-23;
    static doublereal c__ = 2.998e8;
    static doublereal twopi = 6.283185;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer k;
    static doublereal c2, fv, ang, arg, wlm, c1fv, c2fv;


/* ccccccccccccccccccccccccccc  p l a n c k  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine evaluates the planck function at 'nlev'         cc */
/* c    temperatures at a specific wavelenght or wavenumber.            cc */
/* c    note: if w = 0.0, this version returns a unit value for b       cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c   itype - calculation type: 1) radiance, 2) flux                   cc */
/* c    nlev - number of levels.                                        cc */
/* c    iuin - input frequency/wavelength unit:                         cc */
/* c           1: wavenumber (cm**-1)                                   cc */
/* c           2: wavelength (microns)                                  cc */
/* c           3: wavelength (nanometers)                               cc */
/* c           4: wavelength (Angstroms)                                cc */
/* c   iunits - output units desired for planck function:               cc */
/* c           0: unit flux: b(k) = 1.0                                 cc */
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
/* c       b - planck function (units specified by iunits)               cc */
/* c    c1fv - numerator of planck function in desired units            cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  p l a n c k  ccccccccccccccccccccccccccccc */





/* ****   define Planks constant, h (J s), Boltzmann's constant, k (J/K) */
/*       and the speed of light, c (m/sec). */

    /* Parameter adjustments */
    --b;
    --t;

    /* Function Body */


/* ****  evaluate the 2nd radiation constant and convert */
/*      input wavenumber or wavelength to Hz */

    if (*iunits == 0) {

/* ****   return a unit flux */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    b[k] = 1.f;
/* L1001: */
	}

	return 0;

    }

    c2 = h__ / bk;
    if (*iuin == 1) {
	fv = c__ * 100. * *w;
	wlm = .01 / *w;
    } else {
	if (*iuin == 2) {
	    fv = c__ * 1e6 / *w;
	    wlm = *w * .01;
	} else {
	    if (*iuin == 3) {
		fv = c__ * 1e9 / *w;
		wlm = *w * 1e-9;
	    } else {
		fv = c__ * 1e10f / *w;
		wlm = *w * 1e-10f;
	    }
	}
    }
    c2fv = c2 * fv;

/* ****   specify the output units of the Planck function. */

    if (*itype == 1) {
	ang = 2.;
    } else {
	ang = twopi;
    }
    if (*iunits == 1) {
	c1fv = ang * 100. * h__ * fv * fv * fv / c__;
    } else {
	if (*iunits == 2) {
/* Computing 5th power */
	    d__1 = wlm, d__2 = d__1, d__1 *= d__1;
	    c1fv = ang * 1e-6 * h__ * c__ * c__ / (d__2 * (d__1 * d__1));
	} else {
	    if (*iunits == 3) {
/* Computing 5th power */
		d__1 = wlm, d__2 = d__1, d__1 *= d__1;
		c1fv = ang * 1e-9 * h__ * c__ * c__ / (d__2 * (d__1 * d__1));
	    } else {
		if (*iunits == 4) {
/* Computing 5th power */
		    d__1 = wlm, d__2 = d__1, d__1 *= d__1;
		    c1fv = ang * 1e-10 * h__ * c__ * c__ / (d__2 * (d__1 * 
			    d__1));
		} else {
		    c1fv = ang * h__ * fv * fv * fv / (c__ * c__);
		}
	    }
	}
    }

/* compute the Planck function based on the temperature profile. */

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	arg = c2fv / t[k];
	if (arg < 200.) {
	    b[k] = c1fv / (exp(arg) - 1.);
	} else {

/* ****        use the Rayleigh-Jeans limit */

	    b[k] = c1fv * exp(-arg);

	}
/* L4201: */
    }

    return 0;
} /* planck_ */

