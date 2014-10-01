/* readmie.f -- translated by f2c (version 20100827).
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

static integer c__132 = 132;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__200 = 200;

/* Subroutine */ int readmie_(integer *mode, integer *iumie, integer *iuaer0, 
	char *miefile, integer *iang, integer *ioffmie, integer *icwlq, 
	integer *icqext, integer *icqsca, integer *icg1, integer *icg2, 
	integer *icf, integer *ioffmom, integer *nstr, integer *io_end__, 
	integer *io_err__, doublereal *wnaer, doublereal *wneof, ftnlen 
	miefile_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , s_wsfe(cilist *), e_wsfe(void), f_open(olist *), f_clos(cllist *
	    ), s_rsle(cilist *), e_rsle(void), do_lio(integer *, integer *, 
	    char *, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double pow_ri(real *, integer *);
    integer s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void),
	     f_rew(alist *);

    /* Local variables */
    static real f[8192];
    static integer i__, k, l, m;
    static real g1[8192], g2[8192];
    static integer ia, ns0, nlb, ntb, mom;
    static doublereal wn0d[8192], dwn0, wl0p, wn0p;
    static real qsca[8192];
    static integer nmom;
    static real qext[8192];
    static integer nmom0;
    static real pmom0[1025];
    static integer icmax, iuaer;
    static doublereal dwnmn;
    static real qexts;
    static char miefil[132];
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer length;
    static char momfil[132];
    static real vector[100];
    static char fileout[1*132];
    static integer nmomstd[8192];
    static real pmomstd[1646592]	/* was [201][8192] */;

    /* Fortran I/O blocks */
    static icilist io___7 = { 0, miefil, 0, "(132a)", 132, 1 };
    static cilist io___9 = { 0, 6, 0, "(/,132a)", 0 };
    static cilist io___12 = { 0, 0, 1, 0, 0 };
    static cilist io___16 = { 1, 0, 1, 0, 0 };
    static icilist io___31 = { 0, momfil, 0, "(132a)", 132, 1 };
    static cilist io___32 = { 0, 6, 0, "(132a)", 0 };
    static cilist io___33 = { 0, 0, 0, 0, 0 };
    static cilist io___34 = { 0, 0, 1, 0, 0 };
    static cilist io___38 = { 0, 6, 0, "(/,1a,/,1a,1pe14.6,1a,1pe14.6)", 0 };
    static cilist io___40 = { 0, 6, 0, "(/,1a,i10,/,1a,i5,1a,1pe12.4)", 0 };
    static cilist io___41 = { 0, 0, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, "(1a,i5,1a,1pe14.6,1a,i5)", 0 };
    static cilist io___45 = { 0, 0, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, "(1a,i5,1a,1pe14.6,1a,i5)", 0 };
    static cilist io___47 = { 0, 0, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, "(/,1a,i5)", 0 };



/* ccccccccccccccccccccccccc  r e a d m i e  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    *******************   p u r p o s e   **********************    cc */
/* c                                                                    cc */
/* c    this subroutine reads  the aerosol extinction, absorption       cc */
/* c    and scattering efficiencies and the scattering asymmetry        cc */
/* c    parameters and phase function moments at all wavelengths.       cc */
/* c    these quantities are normalized to their values at a reference  cc */
/* c    wavelength and then written out into a binary file that is      cc */
/* c    monochromatically increasing in wavenumber.                     cc */
/* c                                                                    cc */
/* c    *******************     i n p u t     **********************    cc */
/* c                                                                    cc */
/* c          ns0 - number of wavelengths where scattering properties   cc */
/* c                are specified                                       cc */
/* c       wl0(i) - ith wavelength at which optical properties are      cc */
/* c                specified                                           cc */
/* c      qext(i) - extinction efficiency for the ith wavelength        cc */
/* c      qsca(i) - scattering efficiency for the ith wavelength        cc */
/* c        g1(i) - asymmetry parameter for forward scattering for      cc */
/* c                the  ith wavelength                                 cc */
/* c        g2(i) - asymmetry parameter for back scattering for         cc */
/* c                the  ith wavelength                                 cc */
/* c         f(i) - forward scattering fraction for two-term            cc */
/* c                Henyey-Greenstein phase function.                   cc */
/* c                                                                    cc */
/* c                      - calling arguments -                         cc */
/* c                                                                    cc */
/* c                                                                    cc */
/* c    *******************    o u t p u t    **********************    cc */
/* c                                                                    cc */
/* c       wl0(i) - ith wavelength at which optical properties are      cc */
/* c                specified                                           cc */
/* c      qext(i) - extinction efficiency for the ith wavelength        cc */
/* c      qsca(i) - scattering efficiency for the ith wavelength        cc */
/* c        g1(i) - asymmetry parameter for forward scattering for      cc */
/* c                the  ith wavelength                                 cc */
/* c        g2(i) - asymmetry parameter for back scattering for         cc */
/* c                the  ith wavelength                                 cc */
/* c         f(i) - forward scattering fraction for two-term            cc */
/* c                Henyey-Greenstein phase function.                   cc */
/* c          ns0 - number of wavelengths where scattering properties   cc */
/* c                are specified                                       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc  r e a d m i e  ccccccccccccccccccccccccccccc */




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





/* *****    a e r o s o l    p a r t i c l e    m o d e    l o o p */

    /* Parameter adjustments */
    --wneof;

    /* Function Body */
    m = *mode;

/* ****    open the input mie parameter file */

    charsp_(miefile, fileout, &length, &c__132, &nlb, &ntb, (ftnlen)132, (
	    ftnlen)1);

    s_copy(miefil, " ", (ftnlen)132, (ftnlen)1);
    s_wsfi(&io___7);
    i__1 = length;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, fileout + (i__ - 1), (ftnlen)1);
    }
    do_fio(&c__1, ".mie", (ftnlen)4);
    e_wsfi();
    s_wsfe(&io___9);
    do_fio(&c__1, " mie file: ", (ftnlen)11);
    i__1 = length;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, fileout + (i__ - 1), (ftnlen)1);
    }
    do_fio(&c__1, ".mie", (ftnlen)4);
    e_wsfe();

    o__1.oerr = 0;
    o__1.ounit = *iumie;
    o__1.ofnmlen = 132;
    o__1.ofnm = miefil;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = "formatted";
    o__1.oblnk = 0;
    f_open(&o__1);

/* ****   create a scatch unit for aerosol optical properties */

    iuaer = *iuaer0 + m;
    cl__1.cerr = 0;
    cl__1.cunit = iuaer;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = iuaer;
    o__1.ofnm = 0;
    o__1.orl = 0;
    o__1.osta = "scratch";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);

/* ****    skip unnecessary descriptive information at top of file */

    i__1 = *ioffmie;
    for (l = 1; l <= i__1; ++l) {
	io___12.ciunit = *iumie;
	i__2 = s_rsle(&io___12);
	if (i__2 != 0) {
	    goto L1241;
	}
	i__2 = e_rsle();
	if (i__2 != 0) {
	    goto L1241;
	}
/* L1201: */
    }
    goto L1301;
L1241:
    *io_end__ = 1;

/* ****    find the last column to read: */

L1301:
    icmax = 1;
    if (*icwlq > icmax) {
	icmax = *icwlq;
    }
    if (*icqext > icmax) {
	icmax = *icqext;
    }
    if (*icqsca > icmax) {
	icmax = *icqsca;
    }
    if (*icg1 > icmax) {
	icmax = *icg1;
    }
    if (*icg2 > icmax) {
	icmax = *icg2;
    }
    if (*icf > icmax) {
	icmax = *icf;
    }

/* ****   read the wavelength, and the extinction, and scattering */
/*       cross-sections or efficiencies and the scattering */
/*       assymetry parameter */

    dwnmn = 1e20;
    k = 0;

L1401:
    io___16.ciunit = *iumie;
    i__1 = s_rsle(&io___16);
    if (i__1 != 0) {
	goto L100001;
    }
    i__2 = icmax;
    for (ia = 1; ia <= i__2; ++ia) {
	i__1 = do_lio(&c__4, &c__1, (char *)&vector[ia - 1], (ftnlen)sizeof(
		real));
	if (i__1 != 0) {
	    goto L100001;
	}
    }
    i__1 = e_rsle();
L100001:
    if (i__1 < 0) {
	goto L1621;
    }
    if (i__1 > 0) {
	goto L1601;
    }

    ++k;
    wn0d[k - 1] = 1e4 / vector[*icwlq - 1];
    qext[k - 1] = vector[*icqext - 1];
    qsca[k - 1] = vector[*icqsca - 1];
    g1[k - 1] = vector[*icg1 - 1];
    if (*icg2 > 0) {
	g2[k - 1] = vector[*icg2 - 1];
	f[k - 1] = vector[*icf - 1];
    } else {
	g2[k - 1] = 0.f;
	f[k - 1] = 1.f;
    }

    dwn0 = (d__1 = *wnaer - wn0d[k - 1], abs(d__1));
    if (dwn0 < dwnmn) {
	dwnmn = dwn0;
	qexts = vector[*icqext - 1];
    }

    goto L1401;

/* ****   the following statement lets the program stumble over */
/*       unnecessary descriptive lines at the top of the file */

L1601:
    if (k == 0) {
	goto L1401;
    }
    *io_err__ = 1;

/* ****   set number of spectral intervals, and eof markers */

L1621:
    ns0 = k;
    if (*io_err__ == 0) {
	*io_end__ = 1;
    }
    if (wn0d[0] < wn0d[k - 1]) {
	wneof[1] = wn0d[0];
	wneof[2] = wn0d[k - 1];
    } else {
	wneof[1] = wn0d[k - 1];
	wneof[2] = wn0d[0];
    }

/* ****   normalize extinction and scattering efficiencies by qexts. */

    i__1 = ns0;
    for (k = 1; k <= i__1; ++k) {
	qext[k - 1] /= qexts;
	qsca[k - 1] /= qexts;
/* L2001: */
    }

/* ****    close efficiency factor file */

    cl__1.cerr = 0;
    cl__1.cunit = *iumie;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* ****   initialize phase function legendre coefficients */

    i__1 = ns0;
    for (k = 1; k <= i__1; ++k) {

	pmomstd[k * 201 - 201] = 1.f;
	for (mom = 1; mom <= 200; ++mom) {
	    pmomstd[mom + k * 201 - 201] = 0.f;
/* L2201: */
	}
/* L2221: */
    }

    if (*iang == 1) {

/* ****     create name of file with legendre-polynomial coefficients */
/*         for the scattering phase function */

	s_copy(momfil, " ", (ftnlen)132, (ftnlen)1);
	s_wsfi(&io___31);
	i__1 = length;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, fileout + (i__ - 1), (ftnlen)1);
	}
	do_fio(&c__1, ".mom", (ftnlen)4);
	e_wsfi();
	s_wsfe(&io___32);
	do_fio(&c__1, " pmom file: ", (ftnlen)12);
	i__1 = length;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, fileout + (i__ - 1), (ftnlen)1);
	}
	do_fio(&c__1, ".mom", (ftnlen)4);
	e_wsfe();

/* ****     open file with legendre coefficients of phase function */

	o__1.oerr = 0;
	o__1.ounit = *iumie;
	o__1.ofnmlen = 132;
	o__1.ofnm = momfil;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "formatted";
	o__1.oblnk = 0;
	f_open(&o__1);

/* ****     skip over unnecessary descriptive records at top of file */

	i__1 = *ioffmom;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___33.ciunit = *iumie;
	    s_rsle(&io___33);
	    e_rsle();
/* L3001: */
	}

	i__1 = ns0;
	for (k = 1; k <= i__1; ++k) {

/* ****         read wavelength and number of phase function coefficients */

	    io___34.ciunit = *iumie;
	    i__2 = s_rsle(&io___34);
	    if (i__2 != 0) {
		goto L3061;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&wl0p, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L3061;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nmom0, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L3061;
	    }
	    i__2 = e_rsle();
	    if (i__2 != 0) {
		goto L3061;
	    }

	    if (wl0p != 0.f) {
		wn0p = 1e4 / wl0p;
	    }
	    if ((d__1 = wn0p - wn0d[k - 1], abs(d__1)) > wn0d[k - 1] * 1e-4) {
		s_wsfe(&io___38);
		do_fio(&c__1, " Error reading phase function moments: ", (
			ftnlen)39);
		do_fio(&c__1, " wn0p =", (ftnlen)7);
		do_fio(&c__1, (char *)&wn0p, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, "  wn0d =", (ftnlen)8);
		do_fio(&c__1, (char *)&wn0d[k - 1], (ftnlen)sizeof(doublereal)
			);
		e_wsfe();
		s_stop("", (ftnlen)0);
	    }
	    nmomstd[k - 1] = nmom0;
	    if (nmomstd[k - 1] > 200) {
		s_wsfe(&io___40);
		do_fio(&c__1, " readmie warning: Number of phase function mo"
			"ments, ", (ftnlen)52);
		do_fio(&c__1, (char *)&nmomstd[k - 1], (ftnlen)sizeof(integer)
			);
		do_fio(&c__1, " exceeds dimension bounds: mxmom=", (ftnlen)33)
			;
		do_fio(&c__1, (char *)&c__200, (ftnlen)sizeof(integer));
		do_fio(&c__1, " at wavenumber, ", (ftnlen)16);
		do_fio(&c__1, (char *)&wn0p, (ftnlen)sizeof(doublereal));
		e_wsfe();
		nmomstd[k - 1] = 200;
	    }

/* ****         read the legendre polynomial coeficients: */

	    io___41.ciunit = *iumie;
	    s_rsle(&io___41);
	    i__2 = nmom0;
	    for (mom = 0; mom <= i__2; ++mom) {
		do_lio(&c__4, &c__1, (char *)&pmom0[mom], (ftnlen)sizeof(real)
			);
	    }
	    e_rsle();
	    i__2 = nmomstd[k - 1];
	    for (mom = 0; mom <= i__2; ++mom) {
		pmomstd[mom + k * 201 - 201] = pmom0[mom];
/* L3021: */
	    }

/* L3041: */
	}

/* ****         close phase function momemnt file */

L3061:
	cl__1.cerr = 0;
	cl__1.cunit = *iumie;
	cl__1.csta = 0;
	f_clos(&cl__1);

    } else {

/* ****     use a henyey-greenstein phase function */

	i__1 = ns0;
	for (k = 1; k <= i__1; ++k) {
	    if (*iang != 3) {
		i__2 = *nstr;
		for (mom = 1; mom <= i__2; ++mom) {
		    pmomstd[mom + k * 201 - 201] = pow_ri(&g1[k - 1], &mom);
/* L3201: */
		}
	    } else {

/* ****           use a double Henyey-Greenstein phase function */

		i__2 = *nstr;
		for (mom = 0; mom <= i__2; ++mom) {
		    pmomstd[mom + k * 201 - 201] = f[k - 1] * pow_ri(&g1[k - 
			    1], &mom) + (1.f - f[k - 1]) * pow_ri(&g2[k - 1], 
			    &mom);
/* L3221: */
		}
	    }
	    nmomstd[k - 1] = *nstr;
/* L3241: */
	}

    }

/* ****   write these values to iuaer unit - in order of increasing */
/*       wavenumber */

    if (wn0d[0] < wn0d[ns0 - 1]) {
	i__1 = ns0;
	for (k = 1; k <= i__1; ++k) {
	    if (nmomstd[k - 1] > 200) {
		s_wsfe(&io___43);
		do_fio(&c__1, "Warning, dimension mxmom=", (ftnlen)25);
		do_fio(&c__1, (char *)&c__200, (ftnlen)sizeof(integer));
		do_fio(&c__1, " is not adequate to resolve phase function at"
			" wn=", (ftnlen)49);
		do_fio(&c__1, (char *)&wn0d[k - 1], (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, " where nmomstd=", (ftnlen)15);
		do_fio(&c__1, (char *)&nmomstd[k - 1], (ftnlen)sizeof(integer)
			);
		e_wsfe();
		nmom = 200;
	    } else {
		nmom = nmomstd[k - 1];
	    }
	    io___45.ciunit = iuaer;
	    s_wsue(&io___45);
	    do_uio(&c__1, (char *)&nmom, (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&wn0d[k - 1], (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&qext[k - 1], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&qsca[k - 1], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&g1[k - 1], (ftnlen)sizeof(real));
	    i__2 = nmom;
	    for (mom = 0; mom <= i__2; ++mom) {
		do_uio(&c__1, (char *)&pmomstd[mom + k * 201 - 201], (ftnlen)
			sizeof(real));
	    }
	    e_wsue();
/* L3641: */
	}
    } else {
	for (k = ns0; k >= 1; --k) {
	    if (nmomstd[k - 1] > 200) {
		s_wsfe(&io___46);
		do_fio(&c__1, "Warning, dimension mxmom=", (ftnlen)25);
		do_fio(&c__1, (char *)&c__200, (ftnlen)sizeof(integer));
		do_fio(&c__1, " is not adequate to resolve phase function at"
			" wn=", (ftnlen)49);
		do_fio(&c__1, (char *)&wn0d[k - 1], (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, " where nmomstd=", (ftnlen)15);
		do_fio(&c__1, (char *)&nmomstd[k - 1], (ftnlen)sizeof(integer)
			);
		e_wsfe();
		nmom = 200;
	    } else {
		nmom = nmomstd[k - 1];
	    }
	    io___47.ciunit = iuaer;
	    s_wsue(&io___47);
	    do_uio(&c__1, (char *)&nmom, (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&wn0d[k - 1], (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&qext[k - 1], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&qsca[k - 1], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&g1[k - 1], (ftnlen)sizeof(real));
	    i__1 = nmom;
	    for (mom = 0; mom <= i__1; ++mom) {
		do_uio(&c__1, (char *)&pmomstd[mom + k * 201 - 201], (ftnlen)
			sizeof(real));
	    }
	    e_wsue();
/* L3661: */
	}
    }

    s_wsfe(&io___48);
    do_fio(&c__1, " Number of input aerosol wavelengths = ", (ftnlen)39);
    do_fio(&c__1, (char *)&ns0, (ftnlen)sizeof(integer));
    e_wsfe();

    al__1.aerr = 0;
    al__1.aunit = iuaer;
    f_rew(&al__1);

    return 0;
} /* readmie_ */

