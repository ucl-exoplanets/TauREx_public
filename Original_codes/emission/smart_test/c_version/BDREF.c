/* BDREF.f -- translated by f2c (version 20100827).
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
static logical c_true = TRUE_;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* RCS version control information: */
/* $Header: BDREF.f,v 2.1 2000/03/27 21:40:51 laszlo Exp $ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
doublereal bdref_(real *wvnmlo, real *wvnmhi, real *mu, real *mup, real *dphi)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, "(//,7(1X,A,/))", 0 };


/*      Supplies surface bi-directional reflectivity. */

/*      This is only a "stub" version. The user must replace this */
/*      by his/her own BDREF function. */


/*      NOTE 1: Bidirectional reflectivity in DISORT is defined */
/*              by Eq. 39 in STWL. */
/*      NOTE 2: Both MU and MU0 (cosines of reflection and incidence */
/*              angles) are positive. */

/*  INPUT: */

/*    WVNMLO : Lower wavenumber (inv cm) of spectral interval */

/*    WVNMHI : Upper wavenumber (inv cm) of spectral interval */

/*    MU     : Cosine of angle of reflection (positive) */

/*    MUP    : Cosine of angle of incidence (positive) */

/*    DPHI   : Difference of azimuth angles of incidence and reflection */
/*                (radians) */


/*   Called by- DREF, SURFAC */
/* +-------------------------------------------------------------------+ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    s_wsfe(&io___1);
    do_fio(&c__1, "To use a bidirectionally reflecting lower boundary you mu"
	    "st", (ftnlen)59);
    do_fio(&c__1, "replace file BDREF.f with your own file. In that file, yo"
	    "u ", (ftnlen)59);
    do_fio(&c__1, "should supply the bidirectional reflectivity, as a functi"
	    "on ", (ftnlen)60);
    do_fio(&c__1, "of the cosine of angle of reflection, the cosine of angle "
	    , (ftnlen)58);
    do_fio(&c__1, "of incidence, and the difference of azimuth angles of ", (
	    ftnlen)54);
    do_fio(&c__1, "incidence and reflection. See DISORT.doc for more informa"
	    "tion", (ftnlen)61);
    do_fio(&c__1, "and subroutine BDREF in file DISOTEST.f for an example.", (
	    ftnlen)55);
    e_wsfe();
    errmsg_("BDREF--Please supply a surface BDRF model", &c_true, (ftnlen)41);
    ret_val = 0.f;
    return ret_val;
} /* bdref_ */

