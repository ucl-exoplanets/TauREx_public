/* pd_out.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int pd_out__(integer *nz0, integer *npd, integer *npd1, 
	integer *ipd, integer *igs, integer *ifrmout, integer *levout, 
	integer *nlout, integer *nlyr, integer *nphi, integer *numu, integer *
	nzup, integer *iu_pd__, integer *iutpd, integer *iuspd, integer *
	modepd, integer *iupdrad, integer *iunits, integer *isptype, integer *
	islit, logical *lsolar, logical *lplanck, doublereal *wn, doublereal *
	wnmin, doublereal *wnmax, doublereal *width, doublereal *dwn, real *
	points, real *units, real *umu0nz, real *ts, real *t, real *p, real *
	rmix, real *alb, real *dtauaer, real *solflx, real *up_s_flx__, real *
	dn_s_flx__, real *dir_s_flx__, real *up_t_flx__, real *dn_t_flx__, 
	real *pd_pert__, real *dns_src__, real *pd_dns_src__, real *ups_src__,
	 real *pd_ups_src__, real *dnt_src__, real *pd_dnt_src__, real *
	upt_src__, real *pd_upt_src__, real *trn_dir__, real *trn_flx__, real 
	*ref_flx__, real *abs_flx__, real *pd_trndir__, real *pd_trnflx__, 
	real *pd_refflx__, real *pd_absflx__, real *rad, real *pd_rad__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int rad_slit__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, doublereal *, 
	    real *, integer *);
    static integer iquit_pd__[141], k, l, n, l1[10], l2[10];
    static real x0[700]	/* was [70][10] */;
    static integer nl, nn;
    static real wl, dn_flx_trn__[70], up_flx_trn__[70];
    static integer naz;
    static real var;
    static integer nze, iuo;
    static real rad0;
    static integer npd0, nfpd, iord, nlev;
    static real wnio;
    static integer nza_1__[3];
    extern /* Subroutine */ int write_pd_flx__(integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    static real spect[5162];
    static integer iscwt;
    static real dn_flx__[70], up_flx__[70], pd_rad0__[70], wnmins, wnmaxs;

    /* Fortran I/O blocks */
    static cilist io___23 = { 0, 0, 0, "(5i10)", 0 };
    static cilist io___24 = { 0, 0, 0, "(2(1pe14.6),10(1pe12.4),/,          "
	    "                                     10(12(1pe12.4),/))", 0 };
    static cilist io___25 = { 0, 0, 0, "(2(1pe14.6),10(1pe12.4),/,          "
	    "                                     10(12(1pe12.4),/))", 0 };
    static cilist io___28 = { 0, 0, 0, 0, 0 };
    static cilist io___29 = { 0, 0, 0, 0, 0 };
    static cilist io___30 = { 0, 0, 0, 0, 0 };
    static cilist io___33 = { 0, 0, 0, "(2(1pe14.6),10(1pe12.4),/,          "
	    "                                   10(12(1pe12.4),/))", 0 };
    static cilist io___34 = { 0, 0, 0, "(2(1pe14.6),10(1pe12.4),/,          "
	    "                                   10(12(1pe12.4),/))", 0 };
    static cilist io___35 = { 0, 0, 0, 0, 0 };
    static cilist io___36 = { 0, 0, 0, 0, 0 };



/* ccccccccccccccccccccccccccc   p d _ o u t   ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine prints monochromatic partial derivatives,       cc */
/* c    or values smoothed with a slit function.                        cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) ergs/s/cm**2/sr/cm-1                               cc */
/* c              5) photons/s/m**2/sr/micron                           cc */
/* c      iutpd - unit numbers for output level-dependent               cc */
/* c              thermal fluxes and thier partial derivatives          cc */
/* c      iuspd - unit numbers for output level-dependent               cc */
/* c              solar fluxes and thier partial derivatives            cc */
/* c    iupdrad - unit numbers for output level-dependent radiances     cc */
/* c              and thier partial derivatives                         cc */
/* c         nz - index of this solar zenith angle (1 to nsol)          cc */
/* c    ifrmout - output file format (1) ascii or (2) binary            cc */
/* c              (2) surface, (3) arbitrary level                      cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       numu - number of zenith angles in input file                 cc */
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
/* ccccccccccccccccccccccccccc   p d _ o u t   ccccccccccccccccccccccccccc */




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





/* ****    units for flux and radiance jacobians */





/* ****   slit function variables */


/* ****   layer-dependent source terms for solar and thermal fluxes */



/* ****    flux transmission and absorption functions at wavenuber, wn */


/* ****   output flux and radiance partial derivatives */


/* ****    output fluxes and radiances */



/* ****    flux transmission and absorption partial */
/*        derivatives for simplified adding method at wavenumber wn */



/* ****   specify the file index for the first output file */
/*       (required for slit function programs) */

    /* Parameter adjustments */
    pd_rad__ -= 54801;
    rad -= 273;
    pd_absflx__ -= 71;
    pd_refflx__ -= 71;
    pd_trnflx__ -= 71;
    pd_trndir__ -= 71;
    --abs_flx__;
    --ref_flx__;
    --trn_flx__;
    --trn_dir__;
    pd_upt_src__ -= 71;
    --upt_src__;
    pd_dnt_src__ -= 71;
    --dnt_src__;
    pd_ups_src__ -= 71;
    --ups_src__;
    pd_dns_src__ -= 71;
    --dns_src__;
    --pd_pert__;
    dn_t_flx__ -= 71;
    up_t_flx__ -= 71;
    dir_s_flx__ -= 771;
    dn_s_flx__ -= 771;
    up_s_flx__ -= 771;
    dtauaer -= 11;
    --alb;
    rmix -= 71;
    --p;
    --t;
    --points;
    iupdrad -= 10513;
    --modepd;
    iuspd -= 11;
    iutpd -= 11;
    --igs;
    --ipd;

    /* Function Body */
    nlev = *nlyr + 1;

/* ****   set io variables */

    if (*isptype == 2) {
	for (n = 1; n <= 10; ++n) {
	    if (*wn <= *wnmin) {
		iquit_pd__[n - 1] = 0;
	    } else {
		if (*wn >= *wnmax && iquit_pd__[n - 1] == 0) {
		    iquit_pd__[n - 1] = 1;
		}
	    }
/* L1001: */
	}
    }

/* *****  define the single precision wavenumber and wavelength */

    wnio = (real) (*wn);
    if (*wn != 0.f) {
	wl = (real) (1e4f / *wn);
    } else {
	wl = 1e6f;
    }

/* ****   initialize state vector varibles */

    i__1 = *npd;
    for (n = 1; n <= i__1; ++n) {

	if ((i__2 = ipd[n], abs(i__2)) == 1) {
	    l1[n - 1] = nlev;
	    l2[n - 1] = nlev;
	    x0[nlev + n * 70 - 71] = p[nlev];
	} else {
	    if ((i__2 = ipd[n], abs(i__2)) == 2) {
		l1[n - 1] = 1;
		l2[n - 1] = nlev;
		i__2 = l2[n - 1];
		for (l = l1[n - 1]; l <= i__2; ++l) {
		    x0[l + n * 70 - 71] = t[l];
/* L1021: */
		}
		x0[nlev + 1 + n * 70 - 71] = *ts;
	    } else {
		if ((i__2 = ipd[n], abs(i__2)) == 3) {
		    l1[n - 1] = 1;
		    l2[n - 1] = nlev;
		    i__2 = l2[n - 1];
		    for (l = l1[n - 1]; l <= i__2; ++l) {
			x0[l + n * 70 - 71] = rmix[l + igs[n - *npd1] * 70];
/* L1041: */
		    }
		} else {
		    if ((i__2 = ipd[n], abs(i__2)) == 4) {
			l1[n - 1] = 1;
			l2[n - 1] = nlev - 1;
			i__2 = l2[n - 1];
			for (l = l1[n - 1]; l <= i__2; ++l) {
			    x0[l + n * 70 - 71] = dtauaer[modepd[n - *npd1] + 
				    l * 10];
/* L1061: */
			}
		    } else {
			l1[n - 1] = nlev;
			l2[n - 1] = nlev;
			x0[nlev + n * 70 - 71] = alb[*nz0];
		    }
		}
	    }
	}

/* L1081: */
    }

/* ****   F l u x    J a c o b i a n s */

    i__1 = *npd;
    for (n = 1; n <= i__1; ++n) {

	npd0 = n;

	if (ipd[n] < 0) {

	    if (*lplanck && *nz0 == 1) {

/* ****          write out thermal flux partials at each level */

		nfpd = iutpd[n + *nz0 * 10] - *iu_pd__ + 5;
		var = alb[*nz0];

		i__2 = nlev;
		for (k = 1; k <= i__2; ++k) {
		    dnt_src__[k] = *units * dnt_src__[k];
		    upt_src__[k] = *units * upt_src__[k];
/* L2001: */
		}

		i__2 = l2[n - 1];
		for (k = l1[n - 1]; k <= i__2; ++k) {
		    pd_dnt_src__[k + n * 70] = *units * pd_dnt_src__[k + n * 
			    70];
		    pd_upt_src__[k + n * 70] = *units * pd_upt_src__[k + n * 
			    70];
/* L2021: */
		}

/* ***            corrections for transmission flux products */

		i__2 = nlev - 1;
		for (k = 1; k <= i__2; ++k) {
		    dn_flx__[k - 1] = dn_t_flx__[k + 140];
		    up_flx__[k - 1] = up_t_flx__[k + 141];
		    dn_flx_trn__[k - 1] = dn_t_flx__[k + 140] * trn_flx__[k];
		    up_flx_trn__[k - 1] = up_t_flx__[k + 141] * trn_flx__[k];
/* L2041: */
		}

/* ****            write thermal fluxes and their jacobians */

		write_pd_flx__(&iutpd[n + *nz0 * 10], &npd0, &ipd[1], &nlev, 
			l1, l2, ifrmout, iunits, isptype, islit, &nfpd, 
			iquit_pd__, wn, wnmin, wnmax, width, dwn, &points[1], 
			umu0nz, &var, &p[1], x0, &pd_pert__[1], up_flx__, 
			dn_flx__, up_flx_trn__, dn_flx_trn__, &trn_flx__[1], &
			ref_flx__[1], &abs_flx__[1], &pd_trnflx__[71], &
			pd_refflx__[71], &pd_absflx__[71], &dnt_src__[1], &
			upt_src__[1], &pd_dnt_src__[71], &pd_upt_src__[71]);

	    }

	    if (*lsolar) {

/* ****           write downward diffuse solar flux partials at each level */

		nfpd = iuspd[n + *nz0 * 10] - *iu_pd__ + 15;
		var = *units * *solflx;

		i__2 = nlev;
		for (k = 1; k <= i__2; ++k) {
		    dns_src__[k] = *units * dns_src__[k];
		    ups_src__[k] = *units * ups_src__[k];
/* L2101: */
		}

		i__2 = l2[n - 1];
		for (k = l1[n - 1]; k <= i__2; ++k) {
		    pd_dns_src__[k + n * 70] = *units * pd_dns_src__[k + n * 
			    70];
		    pd_ups_src__[k + n * 70] = *units * pd_ups_src__[k + n * 
			    70];
/* L2121: */
		}

		i__2 = nlev - 1;
		for (k = 1; k <= i__2; ++k) {
		    dn_flx__[k - 1] = dn_s_flx__[k + (*nz0 + 20) * 70] + 
			    dir_s_flx__[k + (*nz0 + 20) * 70];
		    up_flx__[k - 1] = up_s_flx__[k + 1 + (*nz0 + 20) * 70];
		    dn_flx_trn__[k - 1] = dn_s_flx__[k + (*nz0 + 20) * 70] * 
			    trn_flx__[k];
		    up_flx_trn__[k - 1] = up_s_flx__[k + 1 + (*nz0 + 20) * 70]
			     * trn_flx__[k];
/* L2141: */
		}

/* ****            write solar fluxes and their jacobians */

		write_pd_flx__(&iuspd[n + *nz0 * 10], &npd0, &ipd[1], &nlev, 
			l1, l2, ifrmout, iunits, isptype, islit, &nfpd, 
			iquit_pd__, wn, wnmin, wnmax, width, dwn, &points[1], 
			umu0nz, &var, &p[1], x0, &pd_pert__[1], up_flx__, 
			dn_flx__, up_flx_trn__, dn_flx_trn__, &trn_flx__[1], &
			ref_flx__[1], &trn_dir__[1], &pd_trnflx__[71], &
			pd_refflx__[71], &pd_trndir__[71], &dns_src__[1], &
			ups_src__[1], &pd_dns_src__[71], &pd_ups_src__[71]);
	    }

	}

/* L2201: */
    }

/* ****    R a d i a n c e    J a c o b i a n s */

/* ****    define the first output zenith angle for radiances */
/*        (only upward radiances are saved at the top fo the atmosphere) */

    i__1 = *nlout;
    for (nl = 1; nl <= i__1; ++nl) {
	if (*levout == 1 || *levout == 3 && nl == 1) {
	    nza_1__[nl - 1] = *nzup;
	} else {
	    nza_1__[nl - 1] = 1;
	}
/* L3001: */
    }

    i__1 = *npd;
    for (n = 1; n <= i__1; ++n) {

	if (ipd[n] > 0) {

	    i__2 = *nlout;
	    for (nl = 1; nl <= i__2; ++nl) {

		i__3 = *nphi;
		for (naz = 1; naz <= i__3; ++naz) {

		    i__4 = *numu;
		    for (nze = nza_1__[nl - 1]; nze <= i__4; ++nze) {

/* ****                     set the unit number */

			iuo = iupdrad[nze + (naz + (n + (nl + *nz0 * 3) * 10 
				<< 4) << 4)];

			if (*wn == *wnmin) {

/* *****                      print print out the pressure profile if */
/*                           this is the first point in the file */

			    if (*ifrmout == 1) {
				io___23.ciunit = iuo;
				s_wsfe(&io___23);
				do_fio(&c__1, (char *)&nlev, (ftnlen)sizeof(
					integer));
				do_fio(&c__1, (char *)&l1[n - 1], (ftnlen)
					sizeof(integer));
				do_fio(&c__1, (char *)&l2[n - 1], (ftnlen)
					sizeof(integer));
				do_fio(&c__1, (char *)&ipd[n], (ftnlen)sizeof(
					integer));
				do_fio(&c__1, (char *)&(*iunits), (ftnlen)
					sizeof(integer));
				e_wsfe();
				io___24.ciunit = iuo;
				s_wsfe(&io___24);
				do_fio(&c__1, (char *)&pd_pert__[n], (ftnlen)
					sizeof(real));
				i__5 = l2[n - 1];
				for (l = l1[n - 1]; l <= i__5; ++l) {
				    do_fio(&c__1, (char *)&x0[l + n * 70 - 71]
					    , (ftnlen)sizeof(real));
				}
				e_wsfe();
				io___25.ciunit = iuo;
				s_wsfe(&io___25);
				do_fio(&c__1, (char *)&(*wnmin), (ftnlen)
					sizeof(doublereal));
				do_fio(&c__1, (char *)&(*wnmax), (ftnlen)
					sizeof(doublereal));
				i__5 = nlev;
				for (k = 1; k <= i__5; ++k) {
				    do_fio(&c__1, (char *)&p[k], (ftnlen)
					    sizeof(real));
				}
				e_wsfe();
			    } else {
				wnmins = (real) (*wnmin);
				wnmaxs = (real) (*wnmax);
				io___28.ciunit = iuo;
				s_wsue(&io___28);
				do_uio(&c__1, (char *)&nlev, (ftnlen)sizeof(
					integer));
				do_uio(&c__1, (char *)&l1[n - 1], (ftnlen)
					sizeof(integer));
				do_uio(&c__1, (char *)&l2[n - 1], (ftnlen)
					sizeof(integer));
				do_uio(&c__1, (char *)&ipd[n], (ftnlen)sizeof(
					integer));
				do_uio(&c__1, (char *)&(*iunits), (ftnlen)
					sizeof(integer));
				e_wsue();
				io___29.ciunit = iuo;
				s_wsue(&io___29);
				do_uio(&c__1, (char *)&pd_pert__[n], (ftnlen)
					sizeof(real));
				i__5 = l2[n - 1];
				for (l = l1[n - 1]; l <= i__5; ++l) {
				    do_uio(&c__1, (char *)&x0[l + n * 70 - 71]
					    , (ftnlen)sizeof(real));
				}
				e_wsue();
				io___30.ciunit = iuo;
				s_wsue(&io___30);
				do_uio(&c__1, (char *)&wnmins, (ftnlen)sizeof(
					real));
				do_uio(&c__1, (char *)&wnmaxs, (ftnlen)sizeof(
					real));
				i__5 = nlev;
				for (k = 1; k <= i__5; ++k) {
				    do_uio(&c__1, (char *)&p[k], (ftnlen)
					    sizeof(real));
				}
				e_wsue();
			    }

			}

/* ****                     print level-dependent partials at wavelength */

			if (*isptype == 1) {

/* ****                       write monochromatic data at full resolution */

			    if (*ifrmout == 1) {

				if (*levout == 1 || *levout == 3 && nl == 1) {

				    rad0 = *units * rad[nze + (naz + 16 << 4)]
					    ;
				    i__5 = l2[n - 1];
				    for (l = l1[n - 1]; l <= i__5; ++l) {
					pd_rad0__[l - 1] = *units * pd_rad__[
						nze + (naz + (nl + (l + n * 
						70) * 3 << 4) << 4)];
/* L3101: */
				    }
				    io___33.ciunit = iuo;
				    s_wsfe(&io___33);
				    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(
					    real));
				    do_fio(&c__1, (char *)&(*wn), (ftnlen)
					    sizeof(doublereal));
				    do_fio(&c__1, (char *)&rad0, (ftnlen)
					    sizeof(real));
				    i__5 = l2[n - 1];
				    for (l = l1[n - 1]; l <= i__5; ++l) {
					do_fio(&c__1, (char *)&pd_rad0__[l - 
						1], (ftnlen)sizeof(real));
				    }
				    e_wsfe();
				} else {

				    if (*levout == 2 || *levout == 3 && nl == 
					    2) {
					rad0 = *units * rad[nze + (naz + (nl 
						<< 4) << 4)];
					i__5 = l2[n - 1];
					for (l = l1[n - 1]; l <= i__5; ++l) {
					    pd_rad0__[l - 1] = *units * 
						    pd_rad__[nze + (naz + (nl 
						    + (l + n * 70) * 3 << 4) 
						    << 4)];
/* L3121: */
					}

					io___34.ciunit = iuo;
					s_wsfe(&io___34);
					do_fio(&c__1, (char *)&wl, (ftnlen)
						sizeof(real));
					do_fio(&c__1, (char *)&(*wn), (ftnlen)
						sizeof(doublereal));
					do_fio(&c__1, (char *)&rad0, (ftnlen)
						sizeof(real));
					i__5 = l2[n - 1];
					for (l = l1[n - 1]; l <= i__5; ++l) {
					    do_fio(&c__1, (char *)&pd_rad0__[
						    l - 1], (ftnlen)sizeof(
						    real));
					}
					e_wsfe();
				    }
				}
			    } else {

				if (*levout == 1 || *levout == 3 && nl == 1) {

				    rad0 = *units * rad[nze + (naz + 16 << 4)]
					    ;
				    i__5 = l2[n - 1];
				    for (l = l1[n - 1]; l <= i__5; ++l) {
					pd_rad0__[l - 1] = *units * pd_rad__[
						nze + (naz + (nl + (l + n * 
						70) * 3 << 4) << 4)];
/* L3141: */
				    }
				    io___35.ciunit = iuo;
				    s_wsue(&io___35);
				    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(
					    real));
				    do_uio(&c__1, (char *)&wnio, (ftnlen)
					    sizeof(real));
				    do_uio(&c__1, (char *)&rad0, (ftnlen)
					    sizeof(real));
				    i__5 = l2[n - 1];
				    for (l = l1[n - 1]; l <= i__5; ++l) {
					do_uio(&c__1, (char *)&pd_rad0__[l - 
						1], (ftnlen)sizeof(real));
				    }
				    e_wsue();
				} else {
				    if (*levout == 2 || *levout == 3 && nl == 
					    2) {
					rad0 = *units * rad[nze + (naz + (nl 
						<< 4) << 4)];
					i__5 = l2[n - 1];
					for (l = l1[n - 1]; l <= i__5; ++l) {
					    pd_rad0__[l - 1] = *units * 
						    pd_rad__[nze + (naz + (nl 
						    + (l + n * 70) * 3 << 4) 
						    << 4)];
/* L3161: */
					}

					io___36.ciunit = iuo;
					s_wsue(&io___36);
					do_uio(&c__1, (char *)&wl, (ftnlen)
						sizeof(real));
					do_uio(&c__1, (char *)&wnio, (ftnlen)
						sizeof(real));
					do_uio(&c__1, (char *)&rad0, (ftnlen)
						sizeof(real));
					i__5 = l2[n - 1];
					for (l = l1[n - 1]; l <= i__5; ++l) {
					    do_uio(&c__1, (char *)&pd_rad0__[
						    l - 1], (ftnlen)sizeof(
						    real));
					}
					e_wsue();
				    }
				}

			    }

			} else {

/* ****                       convolve the data with a slit function. */
/*                           set slit function index */

			    nfpd = iuo - *iu_pd__ + 15;

/* ****                      increment counting variables */

			    points[nfpd] += 1.f;

			    nn = 1;
			    if (*levout == 1 || *levout == 3 && nl == 1) {
				spect[nn - 1] = *units * rad[nze + (naz + 16 
					<< 4)];
			    } else {
				if (*levout == 2 || *levout == 3 && nl == 2) {
				    spect[nn - 1] = *units * rad[nze + (naz + 
					    (nl << 4) << 4)];
				}
			    }

			    i__5 = l2[n - 1];
			    for (l = l1[n - 1]; l <= i__5; ++l) {
				++nn;
				spect[nn - 1] = *units * pd_rad__[nze + (naz 
					+ (nl + (l + n * 70) * 3 << 4) << 4)];
/* L3181: */
			    }

			    iord = 2;
			    iscwt = 0;

/* ****                        call the slit function program */

			    rad_slit__(&iuo, &nfpd, ifrmout, islit, &iscwt, &
				    iord, &nn, wnmin, wnmax, dwn, width, &
				    points[1], wn, spect, &iquit_pd__[nfpd - 
				    1]);

			}

/* L3201: */
		    }
/* L3221: */
		}
/* L3241: */
	    }

	}

/* L3261: */
    }

    return 0;
} /* pd_out__ */

