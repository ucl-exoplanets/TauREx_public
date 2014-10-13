/* find_units.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int find_units__(integer *iuin, integer *iunits, doublereal *
	wn, real *units)
{
    /* Initialized data */

    static doublereal hc = 1.9864421e-25;


/* ccccccccccccccccccccc  f i n d _ u n i t s  ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine sets the scaling factor for converting radiance cc */
/* c    units between several choices, specified by the index, iunits.  cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c       iuin - index of input radiance units:                        cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) ergs/s/cm**2/sr/cm-1                               cc */
/* c              5) photons/s/m**2/sr/micron                           cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) ergs/s/cm**2/sr/cm-1                               cc */
/* c              5) photons/s/m**2/sr/micron                           cc */
/* c         wn - wavenumber (cm-1)                                     cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      units - conversion factor needed to convert units from        cc */
/* c              iuin to iunits.                                       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccc  f i n d _ u n i t s  ccccccccccccccccccccccccccc */






/* ****   define the product of the speed of light, c (2.99792458e+8 m/s) */
/*       and Plancks Constant, h (6.626068e-34 m**2 kg/s) */


/* ****   scale output units to W/m**2/sr/cm**-1 */

    if (*iuin != *iunits && *wn != 0.) {

	if (*iuin == 1) {

/* ****       input units are W/m**2/cm**-1 */

	    if (*iunits == 2) {

/* ****         convert from W/m**2/cm**-1 to W/m**2/micron */

		*units = (real) (*wn * 1e-4 * *wn);

	    } else {

		if (*iunits == 3) {

/* ****           convert from W/m**2/cm**-1 to W/m**2/nm */

		    *units = (real) (*wn * 1e-7 * *wn);

		} else {

		    if (*iunits == 4) {

/* ****           convert from W/m**2/cm**-1 to erg/s/cm**2/cm**-1 */

			*units = 1e3f;

		    } else {

/* ****           convert from W/m**2/cm**-1 to photons/s/m**2/micron */

			*units = (real) (*wn * 1e-6 / hc);

		    }

		}
	    }

	} else {

	    if (*iuin == 2) {

/* ****          input units are W/m**2/micron */

		if (*iunits == 1) {

/* ****           convert from W/m**2/micron to W/m**2/cm**-1 */

		    *units = (real) (1e4 / (*wn * *wn));

		} else {

		    if (*iunits == 3) {

/* ****             convert from W/m**2/micron to W/m**2/nm */

			*units = .001f;

		    } else {

			if (*iunits == 4) {

/* ****               convert from W/m**2/micron to erg/s/cm**2/cm**-1 */

			    *units = (real) (1e7 / (*wn * *wn));

			} else {

/* ****               convert from W/m**2/micron to photons/s/m**2/micron */

			    *units = (real) (.01 / (*wn * hc));

			}

		    }

		}

	    } else {

		if (*iuin == 3) {

/* ****           the input units are W/m**2/nm */

		    if (*iunits == 1) {

/* ****             convert from W/m**2/nm to W/m**2/cm**-1 */

			*units = (real) (1e7 / (*wn * *wn));

		    } else {

			if (*iunits == 2) {

/* ****               convert from W/m**2/nm to W/m**2/micron */

			    *units = 1e3f;

			} else {

			    if (*iunits == 4) {

/* ****                 convert from W/m**2/nm to erg/s/cm**2/cm**-1 */

				*units = (real) (1e10 / (*wn * *wn));

			    } else {

/* ****                 convert from W/m**2/nm to photons/s/m**2/micron */

				*units = (real) (10. / (*wn * hc));

			    }

			}

		    }

		} else {

		    if (*iuin == 4) {

/* ****            the input units are erg/s/cm**2/cm**-1 */

			if (*iunits == 1) {

/* ****               convert from erg/s/cm**2/cm**-1 to W/m**2/cm**-1 */

			    *units = .001f;

			} else {

			    if (*iunits == 2) {

/* ****                convert from erg/s/cm**2/cm**-1 to W/m**2/micron */

				*units = (real) (*wn * 1e-7 * *wn);

			    } else {

				if (*iunits == 3) {

/* ****                   convert from erg/s/cm**2/cm**-1 to W/m**2/nm */

				    *units = (real) (*wn * 1e-10 * *wn);

				} else {

/* ****                   convert from erg/s/cm**2/cm**-1 to */
/*                       photons/sec/m**2/micron */

				    *units = (real) (*wn * 1e-9 / hc);

				}

			    }

			}

		    } else {

/* ****             the input units are photons/sec/m**2/micron */

			if (*iunits == 1) {

/* ****              convert from photons/sec/m**2/micron to W/m**2/cm**-1 */

			    *units = (real) (hc * 1e6 / *wn);

			} else {

			    if (*iunits == 2) {

/* ****                 convert from photons/sec/m**2/micron to */
/*                     W/m**2/micron */

				*units = (real) (*wn * 100. * hc);

			    } else {

				if (*iunits == 3) {

/* ****                   convert from photons/sec/m**2/micron to */
/*                       W/m**2/nm */

				    *units = (real) (*wn * .1 * hc);

				} else {

/* ****                   convert from photons/sec/m**2/micron to */
/*                       erg/sec/cm**2/sr/cm-1 */

				    *units = (real) (hc * 1e9 / *wn);

				}

			    }

			}

		    }

		}

	    }

	}

    } else {

	*units = 1.f;

    }

    return 0;
} /* find_units__ */

