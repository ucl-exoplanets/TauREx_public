/* readsur.f -- translated by f2c (version 20100827).
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
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int readsur_(char *surfile, integer *iusur, integer *iopen, 
	integer *iform, char *formxy, integer *ioff, integer *nalbmx, integer 
	*iwn, integer *icwn, integer *nref, integer *icalb, real *scwn, real *
	scalb, doublereal *wnmin, doublereal *wnmax, integer *io_end__, 
	integer *io_err__, doublereal *wneof, ftnlen surfile_len, ftnlen 
	formxy_len)
{
    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_clos(cllist *), f_open(olist *), s_wsue(cilist *), do_uio(
	    integer *, char *, ftnlen), e_wsue(void), f_rew(alist *);

    /* Local variables */
    static integer k, nr;
    static real add[80], scz[80];
    static doublereal wn_0__[32766];
    static integer icol[80], ncol, ierr, nptj, invz, ncdim, ncmax, iskip;
    static real array[80]	/* was [1][80] */;
    static integer nzdim;
    static real surf_0__[131064]	/* was [32766][4] */;
    extern /* Subroutine */ int readfil_(integer *, integer *, char *, char *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, real *, real *, real *, integer 
	    *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___14 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___15 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___19 = { 0, 6, 0, "(/,1a,i5,/)", 0 };
    static cilist io___20 = { 0, 0, 0, 0, 0 };
    static cilist io___21 = { 0, 0, 0, 0, 0 };



/* cccccccccccccccccccccccccc  r e a d s u r  ccccccccccccccccccccccccccc */
/* c                                                                   cc */
/* c    ********************   p u r p o s e   *********************   cc */
/* c                                                                   cc */
/* c    this subroutine reads surface optical properties as a function cc */
/* c    ofwavelength or wavenumber.  These quantities are then         cc */
/* c    rewritten into a binary file that lists albedos monotonically  cc */
/* c    with wavenumber.                                               cc */
/* c                                                                   cc */
/* c    ********************     i n p u t     *********************   cc */
/* c                                                                   cc */
/* c    iusur - input unit number for surface albedos and              cc */
/* c             moments of surface phase functions.                   cc */
/* c     iusur - output unit number for surface albedos.               cc */
/* c     iopen - open unit flag - 0) don't open unit, 1) open unit.    cc */
/* c   surfile - name of input file                                    cc */
/* c       iwn - type of frequency coordinate for input albedos:       cc */
/* c              1) wavelength, 2) wavenumber.                        cc */
/* c      icwn - column index of spectral quantity in each record.     cc */
/* c     iclab - column index of surface optics in each record.        cc */
/* c     iform - format of input file:                                 cc */
/* c             1) formatted, 2) unformatted, 3) list directed)       cc */
/* c    formxy - format for formatted data files (iform=1 only)        cc */
/* c      scwn - multiplicative factor to convert wavelength - microns cc */
/* c     scalb - multiplicative factor to convert albedo to range(0-1).cc */
/* c      ioff - number of records to skip at top of file              cc */
/* c    nalbmx - maximum number of albedos in input file               cc */
/* c   wn_0(l) - wavelength or wavenumber of each input albedo.        cc */
/* c    surf_0 - input surface optical properties                      cc */
/* c     wnmin - minimum wavenumber where albedos are needed (cm**-1)  cc */
/* c     wnmax - maximum wavenumber where albedos are needed (cm**-1)  cc */
/* c                                                                   cc */
/* c    ********************   o u t p u t     *********************   cc */
/* c                                                                   cc */
/* c  wnalb(l) - wavenumber of each input albedo.                      cc */
/* c    surf_0 - input surface optical properties at each wavelength   cc */
/* c                                                                   cc */
/* cccccccccccccccccccccccccc  r e a d s u r cccccccccccccccccccccccccccc */




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








/* *****    initialize i/o variables */

    /* Parameter adjustments */
    --scalb;
    --icalb;
    --wneof;

    /* Function Body */
    ierr = 0;

    ncdim = 80;
    nzdim = 1;
    iskip = 0;
    invz = 2;
    icol[0] = *icwn;
    ncmax = icol[0];
    scz[0] = *scwn;
    add[0] = 0.f;
    ncol = *nref + 1;
    i__1 = *nref;
    for (nr = 1; nr <= i__1; ++nr) {
	icol[nr] = icalb[nr];
	if (icol[nr] > ncmax) {
	    ncmax = icol[nr];
	}
	scz[nr] = scalb[nr];
	add[nr] = 0.f;
/* L1001: */
    }

    if (*iopen == 1) {

/* ****    o p e n   s u r f a c e    a l b e d o    f i l e */

	readfil_(&c__0, iusur, surfile, formxy, iform, &nzdim, &ncdim, ioff, &
		iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &
		ierr, surfile_len, formxy_len);

	if (ierr != 0) {
	    s_wsfe(&io___14);
	    do_fio(&c__1, " readsur: error opening file: ", (ftnlen)30);
	    do_fio(&c__1, surfile, surfile_len);
	    e_wsfe();
	    s_stop("", (ftnlen)0);
	}
    }

/* ****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s . */

    if (*ioff >= 1) {

	readfil_(&c_n1, iusur, surfile, formxy, iform, &nzdim, &ncdim, ioff, &
		iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &
		ierr, surfile_len, formxy_len);

	if (ierr != 0) {
	    s_wsfe(&io___15);
	    do_fio(&c__1, "End-of-file Encountered while skipping records in"
		    " file ", (ftnlen)55);
	    do_fio(&c__1, surfile, surfile_len);
	    e_wsfe();
	    s_stop("", (ftnlen)0);

	}
    }

/* *****   read each surface albedo record */

    k = 0;

L2001:
    nptj = 1;

    readfil_(&c__1, iusur, surfile, formxy, iform, &nzdim, &ncdim, ioff, &
	    iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &ierr, 
	    surfile_len, formxy_len);

    if (ierr == 0) {

/* ***      load wavenumber vector */

	if (*iwn == 1) {
	    if (array[0] != 0.f) {
		wn_0__[k] = 1e4 / array[0];
	    } else {
		wn_0__[k] = 1e10;
	    }
	} else {
	    wn_0__[k] = array[0];
	}

/* ****     determine if this value is needed. */

	if (k == 0) {

	    ++k;

/* ****      initialize the very first point. */

	    i__1 = *nref;
	    for (nr = 1; nr <= i__1; ++nr) {
		surf_0__[k + nr * 32766 - 32767] = array[nr];
/* L2281: */
	    }

	} else {

/*  if(k .eq. 0) k is not equal to zero */
	    if (wn_0__[k] > wn_0__[k - 1]) {

/* ****        input data is ordered in order of increasing wavenumber. */

		if (wn_0__[k] < *wnmin) {

/* ****              wavenumber is too small - reset intial value */

		    k = 1;
		    wn_0__[k - 1] = wn_0__[k];
		    i__1 = *nref;
		    for (nr = 1; nr <= i__1; ++nr) {
			surf_0__[k + nr * 32766 - 32767] = array[nr];
/* L2201: */
		    }

		} else {

/* ****              wavenumber is too large - add this point and quit */

		    ++k;
		    i__1 = *nref;
		    for (nr = 1; nr <= i__1; ++nr) {
			surf_0__[k + nr * 32766 - 32767] = array[nr];
/* L2221: */
		    }

		    if (wn_0__[k - 1] > *wnmax) {
			goto L2402;
		    }

		}

/* if(wn_0(k+1) .lt. wnmin) */
	    } else {

/* ****        input data is ordered in order of increasing wavelength */

/* wn_0(k+1) < wn_0(k) */
		if (wn_0__[k] > *wnmax) {

/* ****          wavenumber is too large - reset intial value */

		    k = 1;
		    wn_0__[k - 1] = wn_0__[k];
		    i__1 = *nref;
		    for (nr = 1; nr <= i__1; ++nr) {
			surf_0__[k + nr * 32766 - 32767] = array[nr];
/* L2241: */
		    }

		} else {

/* ****              add this point */

/* if(wn_0(k+1) .gt. wnmax) */
		    ++k;
		    i__1 = *nref;
		    for (nr = 1; nr <= i__1; ++nr) {
			surf_0__[k + nr * 32766 - 32767] = array[nr];
/* L2261: */
		    }

		    if (wn_0__[k - 1] < *wnmin) {
			goto L2402;
		    }

		}

/* if(wn_0(k+1) .gt. wnmax) */
	    }

/*  if(wn_0(k+1) .gt. wn_0(k)) */
	}

/*  if(k .eq. 0) */
	if (k < *nalbmx || *nalbmx == 0) {
	    goto L2001;
	}

    }

/* ****    You are done.  reaorder arrays if necessary and close */

/* if(ierr .eq. 0) */
L2402:
    *nalbmx = k;
    if (ierr == 1) {
	*io_end__ = 1;
    } else {
	*io_err__ = 1;
    }

    s_wsfe(&io___19);
    do_fio(&c__1, " Number of input surface optical properties = ", (ftnlen)
	    46);
    do_fio(&c__1, (char *)&(*nalbmx), (ftnlen)sizeof(integer));
    e_wsfe();

    cl__1.cerr = 0;
    cl__1.cunit = *iusur;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* ****   open a scratch unit for albedos, monotonically ordered in */
/*       increasing wavenumber. */

    o__1.oerr = 0;
    o__1.ounit = *iusur;
    o__1.ofnm = 0;
    o__1.orl = 0;
    o__1.osta = "scratch";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);

    if (wn_0__[0] < wn_0__[*nalbmx - 1]) {
	i__1 = *nalbmx;
	for (k = 1; k <= i__1; ++k) {
	    io___20.ciunit = *iusur;
	    s_wsue(&io___20);
	    do_uio(&c__1, (char *)&wn_0__[k - 1], (ftnlen)sizeof(doublereal));
	    i__2 = *nref;
	    for (nr = 1; nr <= i__2; ++nr) {
		do_uio(&c__1, (char *)&surf_0__[k + nr * 32766 - 32767], (
			ftnlen)sizeof(real));
	    }
	    e_wsue();
/* L3001: */
	}
	wneof[1] = wn_0__[0];
	wneof[2] = wn_0__[*nalbmx - 1];
    } else {
	for (k = *nalbmx; k >= 1; --k) {
	    io___21.ciunit = *iusur;
	    s_wsue(&io___21);
	    do_uio(&c__1, (char *)&wn_0__[k - 1], (ftnlen)sizeof(doublereal));
	    i__1 = *nref;
	    for (nr = 1; nr <= i__1; ++nr) {
		do_uio(&c__1, (char *)&surf_0__[k + nr * 32766 - 32767], (
			ftnlen)sizeof(real));
	    }
	    e_wsue();
/* L3021: */
	}
	wneof[1] = wn_0__[*nalbmx - 1];
	wneof[2] = wn_0__[0];
    }

    al__1.aerr = 0;
    al__1.aunit = *iusur;
    f_rew(&al__1);

    return 0;
} /* readsur_ */

