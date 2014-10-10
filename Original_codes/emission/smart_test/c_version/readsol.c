/* readsol.f -- translated by f2c (version 20100827).
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
static integer c__3 = 3;

/* Subroutine */ int readsol_(char *solfile, integer *iusol1, integer *iusol2,
	 integer *iopen, integer *iform, char *formxy, integer *ioff, integer 
	*iwn, integer *icwn, integer *icsol, real *scwn, real *scsol, integer 
	*iuin, doublereal *wnmin, doublereal *wnmax, real *au, integer *
	io_end__, integer *io_err__, doublereal *wneof, ftnlen solfile_len, 
	ftnlen formxy_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     f_rew(alist *), s_wsue(cilist *), do_uio(integer *, char *, 
	    ftnlen), e_wsue(void), s_rsue(cilist *), e_rsue(void), f_clos(
	    cllist *);

    /* Local variables */
    static integer k, iu;
    static real au2;
    extern /* Subroutine */ int find_units__(integer *, integer *, doublereal 
	    *, real *);
    static doublereal wn0[32767];
    static real add[20], scz[20], sol0[32767];
    static integer ius0, icol[20], ncol, ierr, nptj, invz, nsol0, ncdim, 
	    ncmax, iskip;
    static real array[20]	/* was [1][20] */;
    static integer nzdim;
    static real units;
    static integer nstot, iunits;
    extern /* Subroutine */ int readfil_(integer *, integer *, char *, char *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, real *, real *, real *, integer 
	    *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, "(1x,1a,i10)", 0 };
    static cilist io___25 = { 0, 0, 0, 0, 0 };
    static cilist io___27 = { 0, 0, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 0, 1, 0, 0 };
    static cilist io___30 = { 0, 0, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };



/* ccccccccccccccccccccccccccc r e a d s o l cccccccccccccccccccccccccccc */
/* c                                                                   cc */
/* c    ********************   p u r p o s e   *********************   cc */
/* c                                                                   cc */
/* c    this subroutine reads solar fluxes as a function of wavelength cc */
/* c    or wavenumber.  These quantities are then rewritten into       cc */
/* c    a binary file that lists fluxes monotonically with wavenumber. cc */
/* c                                                                   cc */
/* c    ********************     i n p u t     *********************   cc */
/* c                                                                   cc */
/* c    iusol1 - input unit number for solar fluxes.                   cc */
/* c     iopen - open unit flag - 0) don't open unit, 1) open unit.    cc */
/* c   solfile - name of input file                                    cc */
/* c       iwn - type of frequency coordinate for input solar fluxes   cc */
/* c              1) wavelength, 2) wavenumber.                        cc */
/* c      iuin - units index for input solar fluxes:                   cc */
/* c              1) Watts/m**2/sr/cm**-1                              cc */
/* c              2) Watts/m**2/sr/micron                              cc */
/* c              3) Watts/m**2/sr/nanometer                           cc */
/* c              4) ergs/s/cm**2/sr/cm-1                              cc */
/* c              5) photons/s/m**2/sr/micron                          cc */
/* c      icwn - column index of spectral quantity in each record.     cc */
/* c     iclab - column index of solar flux in each record.            cc */
/* c     iform - format of input file:                                 cc */
/* c             1) formatted, 2) unformatted, 3) list directed)       cc */
/* c    formxy - format for formatted data files (iform=1 only)        cc */
/* c      scwn - multiplicative factor to convert wavelength - microns cc */
/* c      ioff - number of records to skip at top of file              cc */
/* c     nsol0 - maximum number of solar fluxe in input file           cc */
/* c    wn0(l) - wavelength or wavenumber of each input solar flux.    cc */
/* c   sol0(l) - input solar flux value                                cc */
/* c     wnmin - minimum wavenumber where fluxes are needed (cm**-1)   cc */
/* c     wnmax - maximum wavenumber where fluxes are needed (cm**-1)   cc */
/* c        au - relative distance from the sun in astronomical units  cc */
/* c                                                                   cc */
/* c    ********************   o u t p u t     *********************   cc */
/* c                                                                   cc */
/* c    wn0(l) - wavenumber of each input solar flux.                  cc */
/* c   sol0(l) - solar flux in units w/m/m/sr/cm-1 at distance, au     cc */
/* c                                                                   cc */
/* ccccccccccccccccccccccccccc r e a d s o l cccccccccccccccccccccccccccc */









    /* Parameter adjustments */
    --wneof;

    /* Function Body */
    ius0 = 0;
    nsol0 = 0;
    nstot = 0;

/* ****  set the default units index iunits = 1 to yield internal solar */
/*      flux units of W/m**2/cm-1 */

    iunits = 1;

/* ****   open a scratch unit for solar flux, monotonically ordered in */
/*       increasing wavenumber. */

    o__1.oerr = 0;
    o__1.ounit = *iusol2;
    o__1.ofnm = 0;
    o__1.orl = 0;
    o__1.osta = "scratch";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);

/* ****   find the square of the relative distance to the sun: */

    au2 = *au * *au;
    if (au2 == 0.f) {
	au2 = 1.f;
    }

/* *****    initialize i/o variables */

    ierr = 0;
    ncol = 2;
    ncdim = 20;
    nzdim = 1;
    iskip = 0;
    invz = 2;
    icol[0] = *icwn;
    icol[1] = *icsol;
    ncmax = icol[1];
    if (icol[0] > ncmax) {
	ncmax = icol[0];
    }
    if (*scwn != 0.f) {
	scz[0] = *scwn;
    } else {
	scz[0] = 1.f;
    }
    scz[1] = *scsol / au2;
    add[0] = 0.f;
    add[1] = 0.f;

    if (*iopen == 1) {

/* ****    o p e n   a t m o s p h e r i c   s t r u c t u r e    f i l e */

	readfil_(&c__0, iusol1, solfile, formxy, iform, &nzdim, &ncdim, ioff, 
		&iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &
		ierr, (ftnlen)132, (ftnlen)40);
	if (ierr != 0) {
	    s_wsle(&io___18);
	    do_lio(&c__9, &c__1, "readsol error opening file: ", (ftnlen)28);
	    do_lio(&c__9, &c__1, solfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
    }

/* ****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s . */

    if (*ioff >= 1) {
	readfil_(&c_n1, iusol1, solfile, formxy, iform, &nzdim, &ncdim, ioff, 
		&iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &
		ierr, (ftnlen)132, (ftnlen)40);

	if (ierr != 0) {
	    s_wsle(&io___19);
	    do_lio(&c__9, &c__1, "End-of-file Encountered while skipping rec"
		    "ords in file ", (ftnlen)55);
	    do_lio(&c__9, &c__1, solfile, (ftnlen)132);
	    e_wsle();
	    s_stop("", (ftnlen)0);

	}
    }

/* *****   read each solar flux record */

    k = 0;

L2403:
    nptj = 1;
    readfil_(&c__1, iusol1, solfile, formxy, iform, &nzdim, &ncdim, ioff, &
	    iskip, &invz, &ncol, icol, &nptj, scz, add, array, &ncmax, &ierr, 
	    (ftnlen)132, (ftnlen)40);

    if (ierr != 0) {
	goto L2803;
    }

/* ***      load wavenumber vector */

    if (*iwn == 1) {
	if (array[0] != 0.f) {
	    wn0[k] = 1e4 / array[0];
	} else {
	    wn0[k] = 1e10;
	}
    } else {
	wn0[k] = array[0];
	if (k > 0) {
	    if (wn0[k] == wn0[k - 1]) {
		goto L2403;
	    }
	}
    }

/* ****     determine if this value is needed. */

    if (k >= 1) {
	if (wn0[k] > wn0[k - 1]) {

/* ****        input data is ordered in order of increasing wavenumber. */

	    if (wn0[k] < *wnmin) {

/* ****            wavenumber is too small - discard */

		k = 1;
		wn0[k - 1] = wn0[k];
		sol0[k - 1] = array[1];
		goto L2403;
	    } else {
		if (wn0[k] > *wnmax) {

/* ****              wavenumber is too large - quit */

		    ++k;
		    sol0[k - 1] = array[1];
		    goto L2803;
		}
	    }
	} else {

/* ****      input data is ordered in order of increasing wavelength */

	    if (wn0[k] > *wnmax) {

/* ****            wavenumber is too large - discard */

		k = 1;
		wn0[k - 1] = wn0[k];
		sol0[k - 1] = array[1];
		goto L2403;
	    } else {
		if (wn0[k] < *wnmin) {

/* ****              wavenumber is too small - quit */

		    ++k;
		    sol0[k - 1] = array[1];
		    goto L2803;
		}
	    }
	}
    }

    ++k;

    sol0[k - 1] = array[1];

    if (k < 32766) {
	goto L2403;
    }

L2803:
    nsol0 = k;

    if (ierr == 1) {
	*io_end__ = 1;
    } else {
	*io_err__ = 1;
    }

    nstot += nsol0;
    s_wsfe(&io___23);
    do_fio(&c__1, "Number of input solar fluxes = ", (ftnlen)31);
    do_fio(&c__1, (char *)&nstot, (ftnlen)sizeof(integer));
    e_wsfe();
    if (wn0[0] < wn0[k - 1]) {
	if (wneof[1] < 0.f) {
	    wneof[1] = wn0[0];
	}
	wneof[2] = wn0[k - 1];
    } else {
	wneof[1] = wn0[k - 1];
	if (wneof[2] <= 0.f) {
	    wneof[2] = wn0[0];
	}
    }

/* ****    convert the input solar fluxes to w/m/m/sr/cm-1: */

    if (*iuin != 1) {

/* ****    find the conversoin factor needed to convert solar */
/*        fluxes from the input units to W/m/m/sr/cm-1 */

	i__1 = nsol0;
	for (k = 1; k <= i__1; ++k) {

	    find_units__(iuin, &iunits, &wn0[k - 1], &units);

	    sol0[k - 1] = units * sol0[k - 1];
/* L3001: */
	}
    }

    if (nsol0 == 0) {
	al__1.aerr = 0;
	al__1.aunit = *iusol2;
	f_rew(&al__1);
	return 0;
    }

    if (wn0[0] < wn0[nsol0 - 1]) {
	wneof[1] = wn0[0];
	wneof[2] = wn0[nsol0 - 1];
	i__1 = nsol0;
	for (k = 1; k <= i__1; ++k) {
	    io___25.ciunit = *iusol2;
	    s_wsue(&io___25);
	    do_uio(&c__1, (char *)&wn0[k - 1], (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&sol0[k - 1], (ftnlen)sizeof(real));
	    e_wsue();
/* L4001: */
	}

/* ****     read additional values if necessary */

	if (wn0[nsol0 - 1] < *wnmax && ierr == 0) {
	    k = 0;
	    nsol0 = 0;
	    goto L2403;
	}
    } else {
	++ius0;
	iu = *iusol2 + ius0;

	o__1.oerr = 0;
	o__1.ounit = iu;
	o__1.ofnm = 0;
	o__1.orl = 0;
	o__1.osta = "scratch";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);

	for (k = nsol0; k >= 1; --k) {
	    io___27.ciunit = iu;
	    s_wsue(&io___27);
	    do_uio(&c__1, (char *)&wn0[k - 1], (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&sol0[k - 1], (ftnlen)sizeof(real));
	    e_wsue();
/* L4021: */
	}

	if (wn0[nsol0 - 1] > *wnmin && ierr == 0) {

/* ****       get the next spectral segment */

	    k = 0;
	    nsol0 = 0;
	    goto L2403;
	} else {

/* ****     put all solar data segments in one file */

	    for (k = ius0; k >= 1; --k) {
		iu = *iusol2 + k;
		s_wsle(&io___28);
		do_lio(&c__9, &c__1, "readsol: rewinding and reading iu=", (
			ftnlen)34);
		do_lio(&c__3, &c__1, (char *)&iu, (ftnlen)sizeof(integer));
		e_wsle();
		al__1.aerr = 0;
		al__1.aunit = iu;
		f_rew(&al__1);
L5001:
		io___29.ciunit = iu;
		i__1 = s_rsue(&io___29);
		if (i__1 != 0) {
		    goto L5021;
		}
		i__1 = do_uio(&c__1, (char *)&wn0[0], (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L5021;
		}
		i__1 = do_uio(&c__1, (char *)&sol0[0], (ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L5021;
		}
		i__1 = e_rsue();
		if (i__1 != 0) {
		    goto L5021;
		}
		io___30.ciunit = *iusol2;
		s_wsue(&io___30);
		do_uio(&c__1, (char *)&wn0[0], (ftnlen)sizeof(doublereal));
		do_uio(&c__1, (char *)&sol0[0], (ftnlen)sizeof(real));
		e_wsue();
		if (wneof[1] <= 0.) {
		    wneof[1] = wn0[0];
		}
		goto L5001;
L5021:
		cl__1.cerr = 0;
		cl__1.cunit = iu;
		cl__1.csta = "delete";
		f_clos(&cl__1);
		s_wsle(&io___31);
		do_lio(&c__9, &c__1, "readsol: closing iu=", (ftnlen)20);
		do_lio(&c__3, &c__1, (char *)&iu, (ftnlen)sizeof(integer));
		e_wsle();
/* L5041: */
	    }
	    wneof[2] = wn0[0];
	}
    }

    al__1.aerr = 0;
    al__1.aunit = *iusol2;
    f_rew(&al__1);

/* ****   set file i/o termination status */

    if (ierr > 0) {
	*io_end__ = 1;
    }
    if (ierr < 0) {
	*io_err__ = 1;
    }

    return 0;
} /* readsol_ */

