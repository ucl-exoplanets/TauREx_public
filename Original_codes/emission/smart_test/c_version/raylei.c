/* raylei.f -- translated by f2c (version 20100827).
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

doublereal raylei_(real *wleff, integer *ncomp, integer *icomp, real *volmix)
{
    /* Initialized data */

    static doublereal delta[6] = { .0279,.078,.021,.058,0.,0. };
    static doublereal a[6] = { 2.871e-4,4.39e-4,2.906e-4,2.663e-4,1.358e-4,
	    3.48e-5 };
    static doublereal b[6] = { .00567,.0064,.0077,.00507,.00752,.0023 };
    static doublereal cnst = 1.031e-24;

    /* System generated locals */
    integer i__1;
    real ret_val;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal r__, r2, sum, wl2i, aniso;


/* ccccccccccccccccccccccccccccc  r a y l e i  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the rayleigh scattering cross section  cc */
/* c    per molecule (meters**-2)  for atmospheres composed of (1) air, cc */
/* c    (2) co2, (3) n2, (4) o2, (5) h2, (6) he, or any combination of  cc */
/* c    these gases.                                                    cc */
/* c    to find the rayleigh scattering optical depth at any level of   cc */
/* c    the atmosphere, these cross sections must be multiplied by      cc */
/* c    the pathlength-integrated number density, n(z)*dz.  in a        cc */
/* c    hydrostatic atmosphere, n(z)*dz = a0*dp/(u0*grav), where u0     cc */
/* c    is the mean molecular mass (kg/mole), and a0 is avagadro's      cc */
/* c    number(kg/kmole)                                                cc */
/* c                                                                    cc */
/* c    r e f e r e n c e s :                                           cc */
/* c                                                                    cc */
/* c    e.j. mccartney, optics of the atmosphere, wiley, p. 187-215,    cc */
/* c          1976.                                                     cc */
/* c    a.t. young, revised depolarization corrections for atmospheric  cc */
/* c          extinction, appl. opt. 19, 3427-3428, 1980.               cc */
/* c    c.w. allen, astrophysical quantities, athlone press, p. 87,     cc */
/* c         1964.                                                      cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    ncomp = number of major atmosphere constituents                 cc */
/* c    icomp = atmosphere constituent index                            cc */
/* c            (1) air  (2) co2  (3) n2  (4) o2 (5) h2 (6) he          cc */
/* c   volmix = volume mixing ratio of each constituent                 cc */
/* c    wleff = effective wavelength (microns)                          cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    raylei - raleigh scattering cross section per molecule (m**2)   cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccc  r a y l e i  ccccccccccccccccccccccccccc */




/*       depolarization factors for air, co2, n2, o2 (young, 1980) */

    /* Parameter adjustments */
    --volmix;
    --icomp;

    /* Function Body */

/* ***    wavelength dependence coefficients for the refractive index */
/*       (allen, 1964) (note: wavelengths must be in microns) */


/* ****   Define the constant multiplier, cnst.  this */
/* 	constant is equal to 24*pi**3/(1.e-24*L**2), */
/* 	where L = loschmidt's number (mks units) , */
/* 	(L = 2.687e25 molecules/m3) and the factor 1.e-24 is */
/* 	needed to convert wl**4 from microns to meters. */


    wl2i = 1. / (*wleff * *wleff);
    sum = 0.f;

    i__1 = *ncomp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__ = a[icomp[i__] - 1] * (b[icomp[i__] - 1] * wl2i + 1.) + 1.;
	aniso = (delta[icomp[i__] - 1] * 3. + 6.) / (6. - delta[icomp[i__] - 
		1] * 7.);
	r2 = r__ * r__ + 2.;
/* Computing 2nd power */
	d__1 = (r__ * r__ - 1.f) / r2;
	sum += volmix[i__] * aniso * (d__1 * d__1);
/* L1201: */
    }

    ret_val = (real) (cnst * wl2i * wl2i * sum);
/* Danie Liang */
    ret_val = 0.f;

    return ret_val;
} /* raylei_ */

