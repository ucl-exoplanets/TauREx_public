/* smart_hdr.f -- translated by f2c (version 20100827).
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

static integer c__8 = 8;
static integer c__30 = 30;
static integer c__1 = 1;

/* Subroutine */ int smart_hdr__(integer *iuout, integer *iutrn, integer *
	iuflx, char *atmfile, char *aerfile, char *miefile, char *solfile, 
	char *surfile, char *mixfile, char *gasfile, integer *iunits, integer 
	*nmodes, integer *ncomp, integer *icomp, real *volmix, real *ts, real 
	*au, real *wgtatm, doublereal *wnmin, doublereal *wnmax, integer *
	isptype, integer *islit, doublereal *width, doublereal *dwn, real *
	tauerr, real *pi0err, real *phferr, real *surferr, integer *nstr, 
	integer *numu, real *umu, integer *nphi, real *phi, integer *nza, 
	real *umu0, real *phi0, integer *levout, integer *nlout, real *pout, 
	real *accur, logical *lamber, integer *isource, integer *irad, 
	integer *ifrmout, real *radius, real *sgrav, integer *igas, integer *
	ngases, integer *ngastyp, ftnlen atmfile_len, ftnlen aerfile_len, 
	ftnlen miefile_len, ftnlen solfile_len, ftnlen surfile_len, ftnlen 
	mixfile_len, ftnlen gasfile_len)
{
    /* System generated locals */
    integer i__1, i__2;
    alist al__1;

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_rew(alist *), s_rsfe(cilist *), e_rsfe(void), s_rsue(cilist *), 
	    do_uio(integer *, char *, ftnlen), e_rsue(void), s_wsfe(cilist *),
	     e_wsfe(void), s_wsue(cilist *), e_wsue(void);
    double acos(doublereal);

    /* Local variables */
    static integer i__, m, n, nl, nz;
    static char val[20];
    static integer iuf, ngt, iuo;
    static real sza0;
    static integer nrec;
    extern /* Subroutine */ int header_(integer *, integer *, char *, char *, 
	    char *, integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer lenrec;
    static char record[80];

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___7 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___9 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___11 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___12 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___13 = { 0, val, 0, "(1x,1pe12.4)", 20, 1 };
    static icilist io___14 = { 0, val, 0, "(1x,l7)", 20, 1 };
    static icilist io___15 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___16 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___17 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___18 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___19 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___20 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___21 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___22 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___23 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___24 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___25 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___26 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___27 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___28 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___29 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___30 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___31 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___32 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___33 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___34 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___35 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___36 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___37 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___39 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___40 = { 0, val, 0, "(1x,i5)", 20, 1 };
    static icilist io___41 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static cilist io___42 = { 0, 0, 1, "(1a)", 0 };
    static cilist io___44 = { 0, 0, 1, 0, 0 };
    static cilist io___46 = { 0, 0, 0, "(1a80)", 0 };
    static cilist io___47 = { 0, 0, 0, 0, 0 };
    static cilist io___49 = { 0, 0, 0, "(1a80)", 0 };
    static cilist io___50 = { 0, 0, 0, 0, 0 };
    static cilist io___52 = { 0, 0, 0, "(1a80)", 0 };
    static cilist io___53 = { 0, 0, 0, 0, 0 };
    static cilist io___54 = { 0, 0, 0, "(1a80)", 0 };
    static cilist io___55 = { 0, 0, 0, 0, 0 };
    static icilist io___57 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };
    static icilist io___58 = { 0, val, 0, "(1x,1pe15.6)", 20, 1 };



/* cccccccccccccccccccccccc  s m a r t _ h d r  cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This subroutine creates a header that documents the input       cc */
/* c    parameters used in the program smart                            cc */
/* c                                                                    cc */
/* c    i n p u t:                                                      cc */
/* c                                                                    cc */
/* c      iusol - unit number with solar radiances                      cc */
/* c     iusol1 - unit number for solar radiance scratch file           cc */
/* c     iuthrm - unit number for thermal radiances                     cc */
/* c      iuout - unit number of output radiance file                   cc */
/* c      iutrn - unit number of output transmission/pressure file      cc */
/* c    atmfile - name of input atmospheric structure file              cc */
/* c    aerfile - name of input aerosol vertical structure file         cc */
/* c    miefile - name of file with aerosol optical properties vs. wn   cc */
/* c    solfile - name of file with wn-dependent solar fluxes           cc */
/* c    surfile - name of file with wn-dependent surface optics         cc */
/* c    mixfile - name of file with gas mixing ratios                   cc */
/* c    gasfile - name fo file with gas absorption coeffiecients vs. wn cc */
/* c     nmodes - number of discrete aerosol partical modes             cc */
/* c      ncomp - number of rayleigh-scattering constituentes           cc */
/* c      icomp - index of each rayleigh scattering constituent         cc */
/* c     volmix - volume mixing ratio of each rayleigh scatterer        cc */
/* c         ts - sufrace temperature (K)                               cc */
/* c         au - distance to the sun (in AU's)                         cc */
/* c     tauerr - optical depth relative binning error (0. to ~0.8)     cc */
/* c     pi0err - co-single scattering albedo absolute binning error    cc */
/* c     phferr - asymmetry factor absolute binning error               cc */
/* c    surferr - surface optical property binning error                cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angles (degrees)                        cc */
/* c       pout - output pressure level (bars)                          cc */
/* c      accur - azimuth convergence accuracy for D/O routine          cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c              note: if lamber = F, a bi-directional reflection      cc */
/* c              function is required.                                 cc */
/* c    isource - index of source function type: (1) solar only         cc */
/* c              (2) thermal only, (3) both                            cc */
/* c       irad - index of output file type:                            cc */
/* c              1) fluxes, radiances, and heating rates,              cc */
/* c                 at computational azimuths and zenith angles,       cc */
/* c              2) fluxes, radiances, heating rates, and transmission cc */
/* c                functions at computational zenith angles,           cc */
/* c              3) fluxes, radiances, heating rates, and contribution cc */
/* c                 functions at computational zenith angles,          cc */
/* c              4) fluxes, radiances, heating rates, transmission     cc */
/* c                 functions and and contribution functions           cc */
/* c                 at computational zenith angles,                    cc */
/* c              5) fluxes, radiances, and heating rates,              cc */
/* c                 at computational azimuths and zenith angles,       cc */
/* c              6) fluxes, radiances and transmission functions       cc */
/* c                 at arbitrary zenith angles,                        cc */
/* c              7) fluxes, radiances, and contribution functions      cc */
/* c                 at arbitrary zenith angles,                        cc */
/* c              8) fluxes, radiances, transmission functions,and      cc */
/* c                 contribution functions at arbitrary zenith angles. cc */
/* c    ifmrout - index of output file format (1) ascii, (2) binary     cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       nlyr - number of computational model layers                  cc */
/* c        nza - number of solar zenith angles                         cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) Watts/m**2/sr/Angstrom                             cc */
/* c        phi - emission azimuth angles (degrees)                     cc */
/* c       umu0 - solar zenith angle cosines                            cc */
/* c       phi0 - solar azimuth angle cosines                           cc */
/* c      wnmin - minimum wavenumber of desired spectral window         cc */
/* c      wnmax - maximum wavenumber of desired spectral window         cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    header file with smart/dart header format                       cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccc  s m a r t _ h d r  cccccccccccccccccccccccccc */



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








    /* Parameter adjustments */
    --ngastyp;
    --igas;
    --pout;
    --phi0;
    --umu0;
    --phi;
    --umu;
    --volmix;
    --icomp;
    gasfile -= 528;
    mixfile -= 132;
    miefile -= 132;
    aerfile -= 132;

    /* Function Body */
    lenrec = 80;
    nrec = 0;

/* ****  input atmospheric structure file. */

    header_(iuout, ifrmout, "atmfile", atmfile, " Atmospheric structure file "
	    , &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)132, (ftnlen)28);
    ++nrec;

    i__1 = *nmodes;
    for (m = 1; m <= i__1; ++m) {
	header_(iuout, ifrmout, "aerfile", aerfile + m * 132, " aerosol vert"
		". str. file", &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)132, 
		(ftnlen)24);
	++nrec;
	header_(iuout, ifrmout, "miefile", miefile + m * 132, " aerosol mie "
		"scat. file", &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)132, (
		ftnlen)23);
	++nrec;
/* L1021: */
    }

/* ****    a t m o s p h e r i c    g a s e s */

    s_wsfi(&io___5);
    do_fio(&c__1, (char *)&(*ngases), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "ngases", val, " number of absorbing gases ", &
	    c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)27);
    ++nrec;

    i__1 = *ngases;
    for (n = 1; n <= i__1; ++n) {
	s_copy(val, " ", (ftnlen)20, (ftnlen)1);
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&igas[n], (ftnlen)sizeof(integer));
	e_wsfi();
	header_(iuout, ifrmout, "igas", val, " afgl gas code for gas", &c__8, 
		&c__30, &lenrec, (ftnlen)4, (ftnlen)20, (ftnlen)22);
	i__2 = ngastyp[n];
	for (ngt = 1; ngt <= i__2; ++ngt) {
	    ++nrec;
	    header_(iuout, ifrmout, "gasfile", gasfile + (ngt + n * 3) * 132, 
		    " Name of abs. coeff.", &c__8, &c__30, &lenrec, (ftnlen)7,
		     (ftnlen)132, (ftnlen)20);
/* L1041: */
	}
	++nrec;
	header_(iuout, ifrmout, "mixfile", mixfile + n * 132, " gas mixing r"
		"atio file", &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)132, (
		ftnlen)22);
	++nrec;
/* L1081: */
    }

/* ****    r a y l e i g h    s c a t t e r i n g */

    s_wsfi(&io___9);
    do_fio(&c__1, (char *)&(*ncomp), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "ncomp", val, " Number of Rayleigh scatterers", &
	    c__8, &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (ftnlen)30);
    ++nrec;
    i__1 = *ncomp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&icomp[i__], (ftnlen)sizeof(integer));
	e_wsfi();
	header_(iuout, ifrmout, "icomp", val, " Index: (1) air (2) co2 (3) n"
		"2 (4) o2", &c__8, &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (
		ftnlen)37);
	++nrec;
	s_wsfi(&io___12);
	do_fio(&c__1, (char *)&volmix[i__], (ftnlen)sizeof(real));
	e_wsfi();
	header_(iuout, ifrmout, "volmix", val, " Volume mixing ratio ", &c__8,
		 &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)21);
	++nrec;
/* L1661: */
    }

/* ****    s u r f a c e    c h a r a c t e r i s t i c s */

    s_wsfi(&io___13);
    do_fio(&c__1, (char *)&(*ts), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "ts", val, " Surface temperature ", &c__8, &c__30,
	     &lenrec, (ftnlen)2, (ftnlen)20, (ftnlen)21);
    ++nrec;

    s_wsfi(&io___14);
    do_fio(&c__1, (char *)&(*lamber), (ftnlen)sizeof(logical));
    e_wsfi();
    header_(iuout, ifrmout, "lamber", val, " Use Lambert surface?", &c__8, &
	    c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)21);
    ++nrec;

    header_(iuout, ifrmout, "surfile", surfile, " surface optical property f"
	    "ile ", &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)132, (ftnlen)31)
	    ;
    ++nrec;

    if (*isource == 1 || *isource == 3) {

/* ****     s o l a r    f l u x e s */

	header_(iuout, ifrmout, "solfile", solfile, " Name of solar flux fil"
		"e ", &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)132, (ftnlen)
		25);

    }

/* ****  read physical properties of planet and atmosphere. */

    s_wsfi(&io___15);
    do_fio(&c__1, (char *)&(*au), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "au", val, " Planets distance from sun (au)", &
	    c__8, &c__30, &lenrec, (ftnlen)2, (ftnlen)20, (ftnlen)31);
    ++nrec;

    s_wsfi(&io___16);
    do_fio(&c__1, (char *)&(*sgrav), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "sgrav", val, " surface gravity (m/s**2) ", &c__8,
	     &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (ftnlen)26);
    ++nrec;

    s_wsfi(&io___17);
    do_fio(&c__1, (char *)&(*radius), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "radius", val, " Radius of the planet (km) ", &
	    c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)27);
    ++nrec;

    s_wsfi(&io___18);
    do_fio(&c__1, (char *)&(*wgtatm), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "wgtatm", val, " Atmospheric molecular weight (k"
	    "g/kmole) ", &c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (
	    ftnlen)41);
    ++nrec;

/* ****    o u t p u t   s p e c t r a l   g r i d. */

    s_wsfi(&io___19);
    do_fio(&c__1, (char *)&(*wnmin), (ftnlen)sizeof(doublereal));
    e_wsfi();
    header_(iuout, ifrmout, "wnmin", val, "minimum output wavenumber", &c__8, 
	    &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (ftnlen)25);
    ++nrec;
    s_wsfi(&io___20);
    do_fio(&c__1, (char *)&(*wnmax), (ftnlen)sizeof(doublereal));
    e_wsfi();
    header_(iuout, ifrmout, "wnmax", val, "maximum output wavenumber", &c__8, 
	    &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (ftnlen)25);
    ++nrec;

    s_wsfi(&io___21);
    do_fio(&c__1, (char *)&(*isptype), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "isptype", val, " Index of output spectrum type ",
	     &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)20, (ftnlen)31);
    ++nrec;

    if (*isptype == 2) {

	s_wsfi(&io___22);
	do_fio(&c__1, (char *)&(*islit), (ftnlen)sizeof(integer));
	e_wsfi();
	header_(iuout, ifrmout, "islit", val, " Type of spectral response fu"
		"nction ", &c__8, &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (
		ftnlen)36);
	++nrec;

	s_wsfi(&io___23);
	do_fio(&c__1, (char *)&(*width), (ftnlen)sizeof(doublereal));
	e_wsfi();
	header_(iuout, ifrmout, "width", val, " Half-width-at-half-max ", &
		c__8, &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (ftnlen)24);
	++nrec;

	s_wsfi(&io___24);
	do_fio(&c__1, (char *)&(*dwn), (ftnlen)sizeof(doublereal));
	e_wsfi();
	header_(iuout, ifrmout, "dwn", val, " output sampling resolution (cm"
		"**-1)", &c__8, &c__30, &lenrec, (ftnlen)3, (ftnlen)20, (
		ftnlen)36);
	++nrec;
    }

    s_wsfi(&io___25);
    do_fio(&c__1, (char *)&(*iunits), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "iunits", val, " Index of output radiance units ",
	     &c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)32);
    ++nrec;

/* ****   set error limits for smt method */

    s_wsfi(&io___26);
    do_fio(&c__1, (char *)&(*tauerr), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "tauerr", val, " Fractional error for smt tau bi"
	    "nning ", &c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)
	    38);
    ++nrec;
    s_wsfi(&io___27);
    do_fio(&c__1, (char *)&(*pi0err), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "pi0err", val, " error for smt pi0 binning ", &
	    c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)27);
    ++nrec;
    s_wsfi(&io___28);
    do_fio(&c__1, (char *)&(*phferr), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "phferr", val, " error for smt <cos> binning ", &
	    c__8, &c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)29);
    ++nrec;
    s_wsfi(&io___29);
    do_fio(&c__1, (char *)&(*surferr), (ftnlen)sizeof(real));
    e_wsfi();
    header_(iuout, ifrmout, "surferr", val, " error for smt surface optical "
	    "property binning ", &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)20,
	     (ftnlen)48);
    ++nrec;

/* ****    d i s c r e t e   o r d i n a n t   m e t h o d */

    s_wsfi(&io___30);
    do_fio(&c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "nstr", val, " Number of streams for dom (> 2) ", 
	    &c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (ftnlen)33);
    ++nrec;

/* ****   o u t p u t    r a d i a n c e    a n g l e s */

    s_wsfi(&io___31);
    do_fio(&c__1, (char *)&(*irad), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "irad", val, " Index of output radiances/fluxes "
	    "dist.", &c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (ftnlen)39)
	    ;
    ++nrec;

    s_wsfi(&io___32);
    do_fio(&c__1, (char *)&(*numu), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "numu", val, " # of arbitrary emission zenith an"
	    "gles ", &c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (ftnlen)39)
	    ;
    ++nrec;

    i__1 = *numu;
    for (n = 1; n <= i__1; ++n) {
	s_wsfi(&io___33);
	do_fio(&c__1, (char *)&umu[n], (ftnlen)sizeof(real));
	e_wsfi();
	header_(iuout, ifrmout, "umu", val, " Emission zenith angle", &c__8, &
		c__30, &lenrec, (ftnlen)3, (ftnlen)20, (ftnlen)22);
	++nrec;
/* L3281: */
    }

    s_wsfi(&io___34);
    do_fio(&c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "nphi", val, " # of emission azimuth angles ", &
	    c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (ftnlen)30);
    ++nrec;

    i__1 = *nphi;
    for (n = 1; n <= i__1; ++n) {
	s_wsfi(&io___35);
	do_fio(&c__1, (char *)&phi[n], (ftnlen)sizeof(real));
	e_wsfi();
	header_(iuout, ifrmout, "phi", val, " Emission azimuth angle", &c__8, 
		&c__30, &lenrec, (ftnlen)3, (ftnlen)20, (ftnlen)23);
	++nrec;
/* L3641: */
    }

/* ****    o u t p u t    l e v e l s */

    s_wsfi(&io___36);
    do_fio(&c__1, (char *)&(*levout), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "levout", val, " output level index: ", &c__8, &
	    c__30, &lenrec, (ftnlen)6, (ftnlen)20, (ftnlen)21);
    ++nrec;

    s_wsfi(&io___37);
    do_fio(&c__1, (char *)&(*nlout), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "nlout", val, " number of output levels: ", &c__8,
	     &c__30, &lenrec, (ftnlen)5, (ftnlen)20, (ftnlen)26);
    ++nrec;

    i__1 = *nlout;
    for (nl = 1; nl <= i__1; ++nl) {
	s_wsfi(&io___39);
	do_fio(&c__1, (char *)&pout[nl], (ftnlen)sizeof(real));
	e_wsfi();
	header_(iuout, ifrmout, "pout", val, " output pressure level (bars): "
		, &c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (ftnlen)31);
	++nrec;
/* L3801: */
    }

/* ****   s o u r c e    f u n c t i o n s */

    s_wsfi(&io___40);
    do_fio(&c__1, (char *)&(*isource), (ftnlen)sizeof(integer));
    e_wsfi();
    header_(iuout, ifrmout, "isource", val, " Index of source function type ",
	     &c__8, &c__30, &lenrec, (ftnlen)7, (ftnlen)20, (ftnlen)31);
    ++nrec;

    if (*isource == 1 || *isource == 3) {

	s_wsfi(&io___41);
	do_fio(&c__1, (char *)&(*accur), (ftnlen)sizeof(real));
	e_wsfi();
	header_(iuout, ifrmout, "accur", val, " Convergence criteria for azi"
		"muthal series", &c__8, &c__30, &lenrec, (ftnlen)5, (ftnlen)20,
		 (ftnlen)42);
    }

/* ****   rewind output file and write header for optical depth file */
/*       and each solar zenith angle. */

    al__1.aerr = 0;
    al__1.aunit = *iuout;
    f_rew(&al__1);
L4001:
    if (*ifrmout == 1) {
	io___42.ciunit = *iuout;
	i__1 = s_rsfe(&io___42);
	if (i__1 != 0) {
	    goto L4201;
	}
	i__1 = do_fio(&c__1, record, (ftnlen)80);
	if (i__1 != 0) {
	    goto L4201;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L4201;
	}
    } else {
	io___44.ciunit = *iuout;
	i__1 = s_rsue(&io___44);
	if (i__1 != 0) {
	    goto L4201;
	}
	i__1 = do_uio(&c__1, record, (ftnlen)80);
	if (i__1 != 0) {
	    goto L4201;
	}
	i__1 = e_rsue();
	if (i__1 != 0) {
	    goto L4201;
	}
    }

/* ****   print out the headers for the files for the remainder of the */
/*       levels at the first zenith angle. */

    i__1 = *nlout;
    for (nl = 2; nl <= i__1; ++nl) {
	iuo = *iuout + nl - 1;
	if (*ifrmout == 1) {
	    io___46.ciunit = iuo;
	    s_wsfe(&io___46);
	    do_fio(&c__1, record, (ftnlen)80);
	    e_wsfe();
	} else {
	    io___47.ciunit = iuo;
	    s_wsue(&io___47);
	    do_uio(&c__1, record, (ftnlen)80);
	    e_wsue();
	}
/* L4101: */
    }

/* *****   if the number of zenith angles is greater than 1, write */
/*        out headers for remaining output zenith angles and levels */

    i__1 = *nza;
    for (nz = 2; nz <= i__1; ++nz) {
	i__2 = *nlout;
	for (nl = 1; nl <= i__2; ++nl) {
	    iuo = *iuout + (nz - 1) * *nlout + nl - 1;
	    if (*ifrmout == 1) {
		io___49.ciunit = iuo;
		s_wsfe(&io___49);
		do_fio(&c__1, record, (ftnlen)80);
		e_wsfe();
	    } else {
		io___50.ciunit = iuo;
		s_wsue(&io___50);
		do_uio(&c__1, record, (ftnlen)80);
		e_wsue();
	    }
/* L4121: */
	}
/* L4141: */
    }

    if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    iuf = *iuflx + nz - 1;
	    if (*ifrmout == 1) {
		io___52.ciunit = iuf;
		s_wsfe(&io___52);
		do_fio(&c__1, record, (ftnlen)80);
		e_wsfe();
	    } else {
		io___53.ciunit = iuf;
		s_wsue(&io___53);
		do_uio(&c__1, record, (ftnlen)80);
		e_wsue();
	    }
/* L4161: */
	}
    }

/* *****   write out to optical depth file */

    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	if (*ifrmout == 1) {
	    io___54.ciunit = *iutrn;
	    s_wsfe(&io___54);
	    do_fio(&c__1, record, (ftnlen)80);
	    e_wsfe();
	} else {
	    io___55.ciunit = *iutrn;
	    s_wsue(&io___55);
	    do_uio(&c__1, record, (ftnlen)80);
	    e_wsue();
	}
    }

    goto L4001;

L4201:

    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	i__2 = *nlout;
	for (nl = 1; nl <= i__2; ++nl) {
	    iuo = *iuout + (nz - 1) * *nlout + nl - 1;

	    if (*isource == 1 || *isource == 3) {

		sza0 = acos(umu0[nz]) * 180.f / 3.14159f;
/*                write(*,*) 'nz,umu0(nz),sza0',nz,umu0(nz),sza0 */
		s_wsfi(&io___57);
		do_fio(&c__1, (char *)&sza0, (ftnlen)sizeof(real));
		e_wsfi();
		header_(&iuo, ifrmout, "sza0", val, " solar zenith angle (de"
			"g)", &c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (
			ftnlen)25);
		s_wsfi(&io___58);
		do_fio(&c__1, (char *)&phi0[nz], (ftnlen)sizeof(real));
		e_wsfi();
		header_(&iuo, ifrmout, "phi0", val, " Solar azimuth angle (d"
			"eg)", &c__8, &c__30, &lenrec, (ftnlen)4, (ftnlen)20, (
			ftnlen)26);
	    }

/* ****            e n d   o f    h e a d e r */

	    header_(&iuo, ifrmout, "end", " 0 ", " End of file", &c__8, &
		    c__30, &lenrec, (ftnlen)3, (ftnlen)3, (ftnlen)12);

/* L5001: */
	}
/* L5021: */
    }

    return 0;
} /* smart_hdr__ */

