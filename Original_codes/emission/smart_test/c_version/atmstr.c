/* atmstr.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int atmstr_(char *atmfile, integer *iuatm, integer *iopen, 
	integer *iform, char *formxy, integer *ioff, integer *nlmax, integer *
	icp, integer *ict, real *scp, real *ratm, real *radius, real *sgrav, 
	real *p, real *t, real *z__, real *alt, real *grav, integer *np_mx__, 
	integer *nlev, ftnlen atmfile_len, ftnlen formxy_len)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double log(doublereal);

    /* Local variables */
    extern /* Subroutine */ int xyinterp_(real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer k, l;
    static real ps[200], ts[200], zs[200];
    static integer nl0;
    static real add[20], scz[20];
    static integer icol[20], ncol, ierr, nptj, invz, ncdim, ncmax;
    static real array[20]	/* was [1][20] */;
    static integer nzdim, iskip, nxmax, nymax;
    static real p5ratm;
    extern /* Subroutine */ int readfil_(integer *, integer *, char *, char *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, real *, real *, real *, integer 
	    *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, "(/,/,5(1a,/),1a,i5,/,1a,2e13.5)", 0 };
    static cilist io___20 = { 0, 6, 0, "(1a,i5)", 0 };



/* ccccccccccccccccccccccccccc  a t m s t r  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    ********************   p u r p o s e   *********************    cc */
/* c                                                                    cc */
/* c    this subroutine reads in standard atmosphere pressures and      cc */
/* c    temperatures, and uses the hydrostatic equation to find the     cc */
/* c    altitudes and gravitational acceleration at each pressure level.cc */
/* c                                                                    cc */
/* c    ********************     i n p u t     **********************   cc */
/* c                                                                    cc */
/* c     iuatm - input unit number for atmospheric variables.           cc */
/* c     iopen - open unit flag - 0) don't open unit, 1) open unit.     cc */
/* c   atmfile - name of input file                                     cc */
/* c       icp - column index of p value in each record.                cc */
/* c       ict - column index of t value in each record.                cc */
/* c     iform - format of input file:                                  cc */
/* c             1) formatted, 2) unformatted, 3) list directed)        cc */
/* c    formxy - format for formatted data files (iform=1 only)         cc */
/* c       scp - multiplicative factor to convert pressure to pascals.  cc */
/* c      ioff - number of records to skip at top of file               cc */
/* c      nlev - number of levels in the input atmosphere:              cc */
/* c             (if nlev=0, the ps in adopted as the input p field.    cc */
/* c             otherwise, values are interpolated to input p-grid.    cc */
/* c     nlmax - maximum number of input levels to read from file.      cc */
/* c             if this variable is set to zero, the program reads     cc */
/* c             all levels in this file.                               cc */
/* c     ps(l) - pressure at level l of standard atmosphere (pascals)   cc */
/* c     ts(l) - temperature at level l of standard atmosphere (k)      cc */
/* c      ratm - gas constant of atmosphere (J/kg/K)                    cc */
/* c     sgrav - surface gravity (m/s**2)                               cc */
/* c    radius - radius of planet (km)                                  cc */
/* c                                                                    cc */
/* c    ********************   o u t p u t     ********************     cc */
/* c                                                                    cc */
/* c       alt(l) - altitude above planet's surface of kth level (km)   cc */
/* c         t(l) - temperature at model level l (k)                    cc */
/* c         p(l) - pressure at model level l (pascals)                 cc */
/* c         z(l) - log of pressure at each level, l                    cc */
/* c      grav(l) - gravitational acceleraton at each level l (m/s**2)  cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  a t m s t r cccccccccccccccccccccccccccccc */










/* ****   atmospheric structure variables */



/* *****    initialize i/o variables */

    /* Parameter adjustments */
    --grav;
    --alt;
    --z__;
    --t;
    --p;

    /* Function Body */
    ierr = 0;
    ncol = 2;
    ncdim = 20;
    nzdim = 1;
    iskip = 0;
    invz = 2;
    nptj = 1;
    icol[0] = *icp;
    icol[1] = *ict;
    ncmax = icol[1];
    if (icol[0] > icol[1]) {
	ncmax = icol[0];
    }
    scz[0] = *scp;
    scz[1] = 1.f;
    add[0] = 0.f;
    add[1] = 0.f;

    if (*iopen == 1) {

/* ****    o p e n   a t m o s p h e r i c   s t r u c t u r e    f i l e */

	readfil_(&c__0, iuatm, atmfile, formxy, iform, &nzdim, &ncdim, ioff, &
		iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &
		ierr, (ftnlen)132, (ftnlen)40);
	if (ierr != 0) {
	    s_wsle(&io___13);
	    do_lio(&c__9, &c__1, "atmstr error opening file: ", (ftnlen)27);
	    do_lio(&c__9, &c__1, atmfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
    }

/* ****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s . */

    if (*ioff >= 1) {
	readfil_(&c_n1, iuatm, atmfile, formxy, iform, &nzdim, &ncdim, ioff, &
		iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &
		ierr, (ftnlen)132, (ftnlen)40);

	if (ierr != 0) {
	    s_wsle(&io___14);
	    do_lio(&c__9, &c__1, "End-of-file Encountered while skipping rec"
		    "ords in file ", (ftnlen)55);
	    do_lio(&c__9, &c__1, atmfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);

	}
    }

/* *****   read each p-T record */

    *nlmax = *np_mx__;
    nl0 = 0;

L2401:
    nptj = 1;
    readfil_(&c__1, iuatm, atmfile, formxy, iform, &nzdim, &ncdim, ioff, &
	    iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &ierr, 
	    (ftnlen)132, (ftnlen)40);

    if (ierr != 0) {
	goto L2801;
    }
    ++nl0;

/* ****       determine if the number of vertical levels exceeds */
/*           the dimension bound:  Note, some of the partial */
/*           deriviative arrays use nlev + 1 levels, so the */
/*           maximum number of levels, nlev must be less than kp */

    if (nl0 >= *np_mx__) {
	s_wsfe(&io___16);
	do_fio(&c__1, " ****  E R R O R  ****    E R R O R  ****  E R R O R "
		" ****", (ftnlen)58);
	do_fio(&c__1, " The number of levels in the model atmosphere exceeds"
		" the", (ftnlen)57);
	do_fio(&c__1, " vertical level dimension, np_mx, in the subroutine A"
		"TMSTR.", (ftnlen)59);
	do_fio(&c__1, " Reduce the number of levels in the model atmosphere,"
		" or", (ftnlen)56);
	do_fio(&c__1, " change np_mx in the calling program and recompile.", (
		ftnlen)51);
	do_fio(&c__1, " current dimension bound, np_mx =", (ftnlen)33);
	do_fio(&c__1, (char *)&(*np_mx__), (ftnlen)sizeof(integer));
	do_fio(&c__1, " pressure and temperature at last level read =", (
		ftnlen)46);
	do_fio(&c__1, (char *)&array[0], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&array[1], (ftnlen)sizeof(real));
	e_wsfe();

	s_stop("", (ftnlen)0);
    }

/* ***      load p-t vectors */

    ps[nl0 - 1] = array[0];
    ts[nl0 - 1] = array[1];

/* ****      define the log of pressure */

    if (ps[nl0 - 1] > 0.f) {
	zs[nl0 - 1] = log(ps[nl0 - 1]);
    } else {
	zs[nl0 - 1] = -1e10f;
    }

    if (nl0 <= *nlmax || *nlmax == 0) {
	goto L2401;
    }

L2801:
    s_wsfe(&io___20);
    do_fio(&c__1, " Number of input levels in p-T file = ", (ftnlen)38);
    do_fio(&c__1, (char *)&nl0, (ftnlen)sizeof(integer));
    e_wsfe();

    if (*nlev == 0) {

/* ****     use the input grid as the default pressure grid */

	*nlev = nl0;
	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    p[k] = ps[k - 1];
	    t[k] = ts[k - 1];
	    z__[k] = zs[k - 1];
/* L3001: */
	}
    } else {

/* ****     interpolate the input temperature field to the input */
/*         pressure grid.  Assume temperatures vary linearly */
/*         with the log of pressure. */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    if (p[k] > 0.f) {
		z__[k] = log(p[k]);
	    } else {
		z__[k] = -1e10f;
	    }
/* L3201: */
	}

	nxmax = *np_mx__;
	nymax = 1;

	xyinterp_(zs, ts, &z__[1], &t[1], &nxmax, &nymax, &nl0, nlev, &c__1);

    }

/* ****   solve the hydrostatic equation to find the altitude and */
/*       gravitational acceleration at each pressure level. */

    p5ratm = *ratm * 5e-4f;
    alt[*nlev] = 0.f;
    grav[*nlev] = *sgrav;
    for (l = *nlev - 1; l >= 1; --l) {
	alt[l] = alt[l + 1] + p5ratm * (t[l + 1] + t[l]) * (z__[l + 1] - z__[
		l]) / grav[l + 1];
/* Computing 2nd power */
	r__1 = *radius / (*radius + alt[l]);
	grav[l] = *sgrav * (r__1 * r__1);
/* L4241: */
    }

/* ****  if top pressure is zero, set altitude there to a nominal value. */

    if (p[1] == 0.f) {
	alt[1] = alt[2] * 2.f;
/* Computing 2nd power */
	r__1 = *radius / (*radius + alt[1]);
	grav[1] = *sgrav * (r__1 * r__1);
    }

    return 0;
} /* atmstr_ */

