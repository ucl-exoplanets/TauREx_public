/* write_pd_flx.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int write_pd_flx__(integer *iuo, integer *npd0, integer *ipd,
	 integer *nlev, integer *l1, integer *l2, integer *ifrmout, integer *
	iunits, integer *isptype, integer *islit, integer *nfpd, integer *
	iquit_pd__, doublereal *wn, doublereal *wnmin, doublereal *wnmax, 
	doublereal *width, doublereal *dwn, real *points, real *umu0nz, real *
	var, real *p, real *x0, real *pd_pert__, real *up_flx__, real *
	dn_flx__, real *up_flx_trn__, real *dn_flx_trn__, real *trnflx1, real 
	*trnflx2, real *trnflx3, real *pd_trnflx1__, real *pd_trnflx2__, real 
	*pd_trnflx3__, real *dn_src__, real *up_src__, real *pd_dn_src__, 
	real *pd_up_src__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12, i__13, i__14;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int rad_slit__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, doublereal *, 
	    real *, integer *);
    static integer k, l, n, nn;
    static real wl;
    static integer l2t, iord;
    static real wnio, spect[5162];
    static integer iscwt;
    static real wnmins, wnmaxs;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, "(5i10)", 0 };
    static cilist io___3 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___5 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___8 = { 0, 0, 0, 0, 0 };
    static cilist io___9 = { 0, 0, 0, 0, 0 };
    static cilist io___10 = { 0, 0, 0, 0, 0 };
    static cilist io___14 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___16 = { 0, 0, 0, 0, 0 };



/* cccccccccccccccccccc  w r i t e _ p d _ f l x   ccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine prints partial derivatives for fluxes           cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        iuo - unit number for output level-dependent partials       cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) ergs/s/cm**2/sr/cm-1                               cc */
/* c              5) photons/s/m**2/sr/micron                           cc */
/* c    isptype - output spectrum type:                                 cc */
/* c              1) Full-resolution spectrum                           cc */
/* c              2) Sampled spectrum smoothed with a slit function     cc */
/* c      islit - index of slit function type:                          cc */
/* c              1) boxcar                                             cc */
/* c              2) triangular (approximate slit spectrometer)         cc */
/* c      width - half-width at half max of slit function  (cm**-1)     cc */
/* c        dwn - output sampling resolution (cm**-1)                   cc */
/* c         wn - current wavenumber (cm**-1)                           cc */
/* c      wnmin - minimum wavenumber in output spectral grid (cm**-1)   cc */
/* c      wnmax - maximum wavenumber in output spectral grid (cm**-1)   cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c         wl - current wavenumber (cm**-1)                           cc */
/* c       wnio - current wavenumber (cm**-1)                           cc */
/* c     pd_rad - radiance partial derivatives at wn                    cc */
/* c                                                                    cc */
/* cccccccccccccccccccc  w r i t e _ p d _ f l x   ccccccccccccccccccccccc */




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








/* ****    flux source functions and jacobians */


/* ****    flux transmission and absorption functions and jacobians */



    /* Parameter adjustments */
    pd_up_src__ -= 71;
    pd_dn_src__ -= 71;
    --up_src__;
    --dn_src__;
    pd_trnflx3__ -= 71;
    pd_trnflx2__ -= 71;
    pd_trnflx1__ -= 71;
    --trnflx3;
    --trnflx2;
    --trnflx1;
    --dn_flx_trn__;
    --up_flx_trn__;
    --dn_flx__;
    --up_flx__;
    --pd_pert__;
    x0 -= 71;
    --p;
    --points;
    --iquit_pd__;
    --l2;
    --l1;
    --ipd;

    /* Function Body */
    n = *npd0;

    if (*wn == *wnmin) {

/* ****    write a header, including number of model levels (nlev), */
/*        the range of levels used for this jacobian (l1,l2), */
/*        the jacobian type index (ipd), the units index (iunits) */

	if (*ifrmout == 1) {
	    io___2.ciunit = *iuo;
	    s_wsfe(&io___2);
	    do_fio(&c__1, (char *)&(*nlev), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&l1[n], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&l2[n], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ipd[n], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*iunits), (ftnlen)sizeof(integer));
	    e_wsfe();

/* ****       write the amplitude of the perturbation used for this */
/*           jacobian, pd_pert, and the nominal value of this state */
/*           vector element at each level where the jacobian is defined */

	    io___3.ciunit = *iuo;
	    s_wsfe(&io___3);
	    do_fio(&c__1, (char *)&(*umu0nz), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&pd_pert__[n], (ftnlen)sizeof(real));
	    i__1 = l2[n];
	    for (k = l1[n]; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&x0[k + n * 70], (ftnlen)sizeof(real));
	    }
	    e_wsfe();

/* ****       write the mininum and maximum wavenumber incldued in the */
/*           jacobian file, and the pressure at each model level */

	    io___5.ciunit = *iuo;
	    s_wsfe(&io___5);
	    do_fio(&c__1, (char *)&(*wnmin), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*wnmax), (ftnlen)sizeof(doublereal));
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&p[k], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	} else {
	    wnmins = (real) (*wnmin);
	    wnmaxs = (real) (*wnmax);
	    io___8.ciunit = *iuo;
	    s_wsue(&io___8);
	    do_uio(&c__1, (char *)&(*nlev), (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&l1[n], (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&l2[n], (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&ipd[n], (ftnlen)sizeof(integer));
	    do_uio(&c__1, (char *)&(*iunits), (ftnlen)sizeof(integer));
	    e_wsue();
	    io___9.ciunit = *iuo;
	    s_wsue(&io___9);
	    do_uio(&c__1, (char *)&(*umu0nz), (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&pd_pert__[n], (ftnlen)sizeof(real));
	    i__1 = l2[n];
	    for (k = l1[n]; k <= i__1; ++k) {
		do_uio(&c__1, (char *)&x0[k + n * 70], (ftnlen)sizeof(real));
	    }
	    e_wsue();
	    io___10.ciunit = *iuo;
	    s_wsue(&io___10);
	    do_uio(&c__1, (char *)&wnmins, (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&wnmaxs, (ftnlen)sizeof(real));
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_uio(&c__1, (char *)&p[k], (ftnlen)sizeof(real));
	    }
	    e_wsue();
	}

	if (l2[n] >= *nlev) {
	    l2t = *nlev - 1;
	} else {
	    l2t = l2[n];
	}

    }

/* *****  define the wavelength */

    wnio = (real) (*wn);
    if (*wn != 0.) {
	wl = (real) (1e4 / *wn);
    } else {
	wl = 1e6f;
    }

    if (*isptype == 1) {

/* ****     create output at full monochromatic resolution */

	if (*ifrmout == 1) {

/* ****        create ascii formatted output */

	    io___14.ciunit = *iuo;
	    s_wsfe(&io___14);
	    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*wn), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*var), (ftnlen)sizeof(real));
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&dn_src__[k], (ftnlen)sizeof(real));
	    }
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		do_fio(&c__1, (char *)&up_src__[k], (ftnlen)sizeof(real));
	    }
	    i__3 = l2[n];
	    for (l = l1[n]; l <= i__3; ++l) {
		do_fio(&c__1, (char *)&pd_dn_src__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__4 = l2[n];
	    for (l = l1[n]; l <= i__4; ++l) {
		do_fio(&c__1, (char *)&pd_up_src__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__5 = *nlev - 1;
	    for (k = 1; k <= i__5; ++k) {
		do_fio(&c__1, (char *)&dn_flx__[k], (ftnlen)sizeof(real));
	    }
	    i__6 = *nlev - 1;
	    for (k = 1; k <= i__6; ++k) {
		do_fio(&c__1, (char *)&up_flx__[k], (ftnlen)sizeof(real));
	    }
	    i__7 = *nlev - 1;
	    for (k = 1; k <= i__7; ++k) {
		do_fio(&c__1, (char *)&dn_flx_trn__[k], (ftnlen)sizeof(real));
	    }
	    i__8 = *nlev - 1;
	    for (k = 1; k <= i__8; ++k) {
		do_fio(&c__1, (char *)&up_flx_trn__[k], (ftnlen)sizeof(real));
	    }
	    i__9 = *nlev - 1;
	    for (k = 1; k <= i__9; ++k) {
		do_fio(&c__1, (char *)&trnflx1[k], (ftnlen)sizeof(real));
	    }
	    i__10 = *nlev - 1;
	    for (k = 1; k <= i__10; ++k) {
		do_fio(&c__1, (char *)&trnflx2[k], (ftnlen)sizeof(real));
	    }
	    i__11 = *nlev - 1;
	    for (k = 1; k <= i__11; ++k) {
		do_fio(&c__1, (char *)&trnflx3[k], (ftnlen)sizeof(real));
	    }
	    i__12 = l2t;
	    for (l = l1[n]; l <= i__12; ++l) {
		do_fio(&c__1, (char *)&pd_trnflx1__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__13 = l2t;
	    for (l = l1[n]; l <= i__13; ++l) {
		do_fio(&c__1, (char *)&pd_trnflx2__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__14 = l2t;
	    for (l = l1[n]; l <= i__14; ++l) {
		do_fio(&c__1, (char *)&pd_trnflx3__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    e_wsfe();

	} else {

/* ****        create binary output */

	    io___16.ciunit = *iuo;
	    s_wsue(&io___16);
	    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&(*var), (ftnlen)sizeof(real));
	    i__1 = *nlev;
	    for (k = 1; k <= i__1; ++k) {
		do_uio(&c__1, (char *)&dn_src__[k], (ftnlen)sizeof(real));
	    }
	    i__2 = *nlev;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&up_src__[k], (ftnlen)sizeof(real));
	    }
	    i__3 = l2[n];
	    for (l = l1[n]; l <= i__3; ++l) {
		do_uio(&c__1, (char *)&pd_dn_src__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__4 = l2[n];
	    for (l = l1[n]; l <= i__4; ++l) {
		do_uio(&c__1, (char *)&pd_up_src__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__5 = *nlev - 1;
	    for (k = 1; k <= i__5; ++k) {
		do_uio(&c__1, (char *)&dn_flx__[k], (ftnlen)sizeof(real));
	    }
	    i__6 = *nlev - 1;
	    for (k = 1; k <= i__6; ++k) {
		do_uio(&c__1, (char *)&up_flx__[k], (ftnlen)sizeof(real));
	    }
	    i__7 = *nlev - 1;
	    for (k = 1; k <= i__7; ++k) {
		do_uio(&c__1, (char *)&dn_flx_trn__[k], (ftnlen)sizeof(real));
	    }
	    i__8 = *nlev - 1;
	    for (k = 1; k <= i__8; ++k) {
		do_uio(&c__1, (char *)&up_flx_trn__[k], (ftnlen)sizeof(real));
	    }
	    i__9 = *nlev - 1;
	    for (k = 1; k <= i__9; ++k) {
		do_uio(&c__1, (char *)&trnflx1[k], (ftnlen)sizeof(real));
	    }
	    i__10 = *nlev - 1;
	    for (k = 1; k <= i__10; ++k) {
		do_uio(&c__1, (char *)&trnflx2[k], (ftnlen)sizeof(real));
	    }
	    i__11 = *nlev - 1;
	    for (k = 1; k <= i__11; ++k) {
		do_uio(&c__1, (char *)&trnflx3[k], (ftnlen)sizeof(real));
	    }
	    i__12 = l2t;
	    for (l = l1[n]; l <= i__12; ++l) {
		do_uio(&c__1, (char *)&pd_trnflx1__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__13 = l2t;
	    for (l = l1[n]; l <= i__13; ++l) {
		do_uio(&c__1, (char *)&pd_trnflx2__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    i__14 = l2t;
	    for (l = l1[n]; l <= i__14; ++l) {
		do_uio(&c__1, (char *)&pd_trnflx3__[l + n * 70], (ftnlen)
			sizeof(real));
	    }
	    e_wsue();

	}

    } else {

/* ****     increment counting variables */

	points[*nfpd] += 1.f;

	nn = 1;
	spect[nn - 1] = *var;

/* ****      downward and upward sources and their jacobians */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = dn_src__[k];
/* L2001: */
	}
	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = up_src__[k];
/* L2021: */
	}
	i__1 = l2[n];
	for (l = l1[n]; l <= i__1; ++l) {
	    ++nn;
	    spect[nn - 1] = pd_dn_src__[l + n * 70];
/* L2041: */
	}
	i__1 = l2[n];
	for (l = l1[n]; l <= i__1; ++l) {
	    ++nn;
	    spect[nn - 1] = pd_up_src__[l + n * 70];
/* L2061: */
	}

/* ***      product of flux and transmission */

	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = dn_flx__[k];
/* L2101: */
	}
	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = up_flx__[k];
/* L2121: */
	}
	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = dn_flx_trn__[k];
/* L2141: */
	}
	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = up_flx_trn__[k];
/* L2161: */
	}

/* ****     transmission variables and their jacobians */

	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = trnflx1[k];
/* L2221: */
	}
	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = trnflx2[k];
/* L2241: */
	}
	i__1 = *nlev - 1;
	for (k = 1; k <= i__1; ++k) {
	    ++nn;
	    spect[nn - 1] = trnflx3[k];
/* L2261: */
	}
	i__1 = l2t;
	for (l = l1[n]; l <= i__1; ++l) {
	    ++nn;
	    spect[nn - 1] = pd_trnflx1__[l + n * 70];
/* L2341: */
	}
	i__1 = l2t;
	for (l = l1[n]; l <= i__1; ++l) {
	    ++nn;
	    spect[nn - 1] = pd_trnflx2__[l + n * 70];
/* L2361: */
	}
	i__1 = l2t;
	for (l = l1[n]; l <= i__1; ++l) {
	    ++nn;
	    spect[nn - 1] = pd_trnflx3__[l + n * 70];
/* L2381: */
	}

	iord = 2;
	iscwt = 0;

	rad_slit__(iuo, nfpd, ifrmout, islit, &iscwt, &iord, &nn, wnmin, 
		wnmax, dwn, width, &points[1], wn, spect, &iquit_pd__[*nfpd]);

    }

    return 0;
} /* write_pd_flx__ */

