/* rad_out.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rad_out__(logical *lsolar, integer *iuout, integer *
	iutrn, integer *iu_flux__, integer *ifrmout, integer *irad, integer *
	levout, integer *nlout, integer *nlyr, integer *nz0, integer *nza, 
	integer *nphi, integer *numu, integer *nzup, integer *isptype, 
	integer *islit, doublereal *wn, doublereal *wnmin, doublereal *wnmax, 
	doublereal *width, doublereal *dwn, real *units, real *p, real *
	solflx, real *up_s_flx__, real *dn_s_flx__, real *dir_s_flx__, real *
	up_t_flx__, real *dn_t_flx__, real *up_flx__, real *dn_flx__, real *
	dir_flx__, real *rad, real *p_ray0__, real *p_gas0__, real *p_aer0__, 
	real *tau_ray_0__, real *tau_gas_0__, real *tau_aer_0__, real *
	p_ray_0__, real *p_gas_0__, real *p_aer_0__, real *points)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double exp(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void)
	    ;

    /* Local variables */
    extern /* Subroutine */ int rad_slit__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, doublereal *, 
	    real *, integer *);
    static real trn_aer_0__, trn_gas_0__;
    static integer k;
    static real trn_ray_0__;
    static integer iquit_rad__[10], iquit_flx__[10], iquit_trn__, nl, nn, no;
    static real wl;
    static integer nz;
    static real dn_dir_flx__[70];
    static integer nfl, naz, nze;
    static real dn_diff_flx__[70], up_diff_flx__[70], rad0[768]	/* was [16][
	    16][3] */;
    static integer iord, nlev;
    static real wnio;
    static integer ntrn;
    static real var_1__;
    static integer nze_1__[3], iurad;
    static real spect[5162];
    static integer iscwt;
    static real dn_flx0__[3], up_flx0__[3], dir_flx0__[3];

    /* Fortran I/O blocks */
    static cilist io___16 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___18 = { 0, 0, 0, 0, 0 };
    static cilist io___29 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___30 = { 0, 0, 0, 0, 0 };
    static cilist io___31 = { 0, 0, 0, "(6(1pe14.6))", 0 };
    static cilist io___32 = { 0, 0, 0, "(14x,9(1pe14.6))", 0 };
    static cilist io___33 = { 0, 0, 0, 0, 0 };
    static cilist io___34 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___35 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___36 = { 0, 0, 0, "(9(1pe14.6))", 0 };
    static cilist io___37 = { 0, 0, 0, 0, 0 };
    static cilist io___38 = { 0, 0, 0, 0, 0 };
    static cilist io___39 = { 0, 0, 0, 0, 0 };



/* ccccccccccccccccccccccccc   r a d _ o u t   ccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine prints backmapped radiances  and other          cc */
/* c    spectral quantities back to a high-resolution spectral grid.    cc */
/* c    It can print monochromatic variables, or smooth these variables cc */
/* c    with a slit function.                                           cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      iuout - unit number of output radiance file                   cc */
/* c      iutrn - unit number of output transmission file               cc */
/* c    iu_flux - unit number for output level-dependent fluxes         cc */
/* c         nz - index of this solar zenith angle (1 to nsol)          cc */
/* c        nza - number of solar zenith angles                         cc */
/* c    ifrmout - output file format (1) ascii or (2) binary            cc */
/* c     levout - output level index (1) top of atmosphere,             cc */
/* c              (2) surface, (3) arbitrary level                      cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       numu - number of zenith angles in input file                 cc */
/* c    isptype - output spectrum type:                                 cc */
/* c              1) Full-resolution spectrum                           cc */
/* c              2) Sampled spectrum smoothed with a slit function     cc */
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
/* c      islit - index of slit function type:                          cc */
/* c              1) boxcar                                             cc */
/* c              2) triangular (approximate slit spectrometer)         cc */
/* c      width - half-width at half max of slit function  (cm**-1)     cc */
/* c        dwn - output sampling resolution (cm**-1)                   cc */
/* c         wn - current wavenumber (cm**-1)                           cc */
/* c     solflx - solar flux at this wavenumber                         cc */
/* c      wnmin - minimum wavenumber in output spectral grid (cm**-1)   cc */
/* c      wnmax - maximum wavenumber in output spectral grid (cm**-1)   cc */
/* c    dir_flx - downward direct flux at wavenumber, wn                cc */
/* c     dn_flx - total downward flux (irradiance) at wavenumber wn     cc */
/* c     up_flx - upward flux (irradiance) at wavenumber wn             cc */
/* c        rad - radiance at each output zenith and azimuth angle      cc */
/* c  tau_ray_0 - normal-incidence rayleigh-scattering optical depth    cc */
/* c  tau_gas_0 - normal-incidence gas absorption optical depth         cc */
/* c  tau_ray_0 - normal-incidence aerosol extinction optical depth     cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c         wl - current wavenumber (cm**-1)                           cc */
/* c       wnio - current wavenumber (cm**-1)                           cc */
/* c    dir_flx - downward direct flux at wavenumber, wn                cc */
/* c     dn_flx - total downward flux (irradiance) at wavenumber wn     cc */
/* c     up_flx - upward flux (irradiance) at wavenumber wn             cc */
/* c        rad - radiance at each output zenith and azimuth angle      cc */
/* c  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc */
/* c  trn_gas_0 - normal-incidence gas transmission                     cc */
/* c  trn_ray_0 - normal-incidence aerosol transmission                 cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc   r a d _ o u t   ccccccccccccccccccccccccccc */




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








/* ****    output radiances */


/* ****    scaled output radiance (units converted to iunits) */



/* ****   slit function variables */


/* ****    number of levels */

    /* Parameter adjustments */
    --points;
    --p_aer0__;
    --p_gas0__;
    --p_ray0__;
    rad -= 273;
    --dir_flx__;
    --dn_flx__;
    --up_flx__;
    dn_t_flx__ -= 71;
    up_t_flx__ -= 71;
    dir_s_flx__ -= 771;
    dn_s_flx__ -= 771;
    up_s_flx__ -= 771;
    --p;

    /* Function Body */
    nlev = *nlyr + 1;

/* ****    index of solar zenith angle */

    nz = *nz0;

/* ****   specify the file index for the first flux output file */
/*       (required for slit function programs) */

    nfl = *nlout * *nza + nz;

/* ****    specify the file index for the transmission output file: */

    ntrn = (*nlout + 1) * *nza + 1;

/* ****   set io variables */

    if (*isptype == 2) {
	if (*wn <= *wnmin) {
	    iquit_rad__[nz - 1] = 0;
	    iquit_trn__ = 0;
	    iquit_flx__[nz - 1] = 0;
	} else {
	    if (*wn >= *wnmax) {
		if (iquit_rad__[nz - 1] == 0) {
		    iquit_rad__[nz - 1] = 1;
		}
		if (iquit_trn__ == 0) {
		    iquit_trn__ = 1;
		}
		if (iquit_flx__[nz - 1] == 0) {
		    iquit_flx__[nz - 1] = 1;
		}
	    }
	}
    }

/* *****  define the wavelength */

    wnio = (real) (*wn);
    if (*wn != 0.) {
	wl = (real) (1e4 / *wn);
    } else {
	wl = 1e9f;
    }

/* ****    define the first output zenith angle for radiances */
/*        (only upward radiances are saved at the top fo the atmosphere) */

    i__1 = *nlout;
    for (nl = 1; nl <= i__1; ++nl) {
	if (*levout == 1 || *levout == 3 && nl == 1) {
	    nze_1__[nl - 1] = *nzup;
	} else {
	    nze_1__[nl - 1] = 1;
	}
/* L1021: */
    }

    trn_ray_0__ = exp(-(*tau_ray_0__));
    trn_gas_0__ = exp(-(*tau_gas_0__));
    trn_aer_0__ = exp(-(*tau_aer_0__));

    if (*lsolar) {
	var_1__ = *units * *solflx;
    } else {
	var_1__ = *tau_ray_0__ + *tau_gas_0__ + *tau_aer_0__;
    }

    if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {
	if (*wn <= *wnmin) {
	    if (*ifrmout == 1) {
		io___16.ciunit = *iu_flux__;
		s_wsfe(&io___16);
		i__1 = nlev;
		for (k = 1; k <= i__1; ++k) {
		    do_fio(&c__1, (char *)&p[k], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    } else {
		io___18.ciunit = *iu_flux__;
		s_wsue(&io___18);
		do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		i__1 = nlev;
		for (k = 1; k <= i__1; ++k) {
		    do_uio(&c__1, (char *)&p[k], (ftnlen)sizeof(real));
		}
		e_wsue();
	    }
	}
    }

/* ****     find  total downward and upward diffuse fluxes at each level */
/*         and convert to W/m**2/cm-1 to iunits */

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	dn_dir_flx__[k - 1] = *units * dir_s_flx__[k + (nz + 20) * 70];
	dn_diff_flx__[k - 1] = *units * (dn_s_flx__[k + (nz + 20) * 70] + 
		dn_t_flx__[k + 140]);
	up_diff_flx__[k - 1] = *units * (up_s_flx__[k + (nz + 20) * 70] + 
		up_t_flx__[k + 140]);
/* L1041: */
    }

/* ***     convert radiances and fluxes form W/m**2/cm-1 to iunits */

    i__1 = *nlout;
    for (k = 1; k <= i__1; ++k) {
	dir_flx0__[k - 1] = *units * dir_flx__[k];
	dn_flx0__[k - 1] = *units * dn_flx__[k];
	up_flx0__[k - 1] = *units * up_flx__[k];
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *numu;
	    for (nze = nze_1__[k - 1]; nze <= i__3; ++nze) {
		rad0[nze + (naz + (k << 4) << 4) - 273] = *units * rad[nze + (
			naz + (k << 4) << 4)];
/* L1101: */
	    }
/* L1121: */
	}
/* L1141: */
    }

    if (*isptype == 1) {

/* ****     write monochromatic data at full resolution */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {

/* ****         compute the unit number */

	    iurad = *iuout + (nz - 1) * *nlout + k - 1;

/* ****         if only top-of-atmosphere radiances are to be saved, */
/*             print out only upward radiances and fluxes. */
/*             For output at other levels, print out both upward */
/*             and downward fluxes and radiances */

	    if (*levout == 1 || *levout == 3 && k == 1) {

/* ****          print out only upward streams at top of atmosphere */

		if (*ifrmout == 1) {

/* ****               print a formatted output file */

		    io___29.ciunit = iurad;
		    s_wsfe(&io___29);
		    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&var_1__, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&up_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__3 = *numu;
			for (nze = nze_1__[k - 1]; nze <= i__3; ++nze) {
			    do_fio(&c__1, (char *)&rad0[nze + (naz + (k << 4) 
				    << 4) - 273], (ftnlen)sizeof(real));
			}
		    }
		    e_wsfe();
		} else {

/* ****               print an unformatted output file */

		    io___30.ciunit = iurad;
		    s_wsue(&io___30);
		    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&var_1__, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&up_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    i__3 = *nphi;
		    for (naz = 1; naz <= i__3; ++naz) {
			i__2 = *numu;
			for (nze = nze_1__[k - 1]; nze <= i__2; ++nze) {
			    do_uio(&c__1, (char *)&rad0[nze + (naz + (k << 4) 
				    << 4) - 273], (ftnlen)sizeof(real));
			}
		    }
		    e_wsue();
		}

	    } else {

/* ****           print out both upward and downward streams. */

		if (*ifrmout == 1) {

/* ****               print a formatted output file */

		    io___31.ciunit = iurad;
		    s_wsfe(&io___31);
		    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&var_1__, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&dir_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&dn_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&up_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    e_wsfe();
		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			io___32.ciunit = iurad;
			s_wsfe(&io___32);
			i__3 = *numu;
			for (nze = nze_1__[k - 1]; nze <= i__3; ++nze) {
			    do_fio(&c__1, (char *)&rad0[nze + (naz + (k << 4) 
				    << 4) - 273], (ftnlen)sizeof(real));
			}
			e_wsfe();
/* L1201: */
		    }
		} else {

/* ****               print an unformatted output file */

		    io___33.ciunit = iurad;
		    s_wsue(&io___33);
		    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&var_1__, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&dir_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    do_uio(&c__1, (char *)&dn_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    do_uio(&c__1, (char *)&up_flx0__[k - 1], (ftnlen)sizeof(
			    real));
		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__3 = *numu;
			for (nze = nze_1__[k - 1]; nze <= i__3; ++nze) {
			    do_uio(&c__1, (char *)&rad0[nze + (naz + (k << 4) 
				    << 4) - 273], (ftnlen)sizeof(real));
			}
		    }
		    e_wsue();
		}
	    }

/* L1221: */
	}

/* ****       print other output radiance quantities: */

	if (*ifrmout == 1) {

/* ****       print formatted output files */

	    if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {


/* ****         print level-dependent fluxes at the wavelength */

		if (*lsolar) {
		    io___34.ciunit = *iu_flux__;
		    s_wsfe(&io___34);
		    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    i__1 = nlev;
		    for (k = 1; k <= i__1; ++k) {
			do_fio(&c__1, (char *)&dn_dir_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    i__3 = nlev;
		    for (k = 1; k <= i__3; ++k) {
			do_fio(&c__1, (char *)&dn_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    i__2 = nlev;
		    for (k = 1; k <= i__2; ++k) {
			do_fio(&c__1, (char *)&up_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    e_wsfe();
		} else {
		    io___35.ciunit = *iu_flux__;
		    s_wsfe(&io___35);
		    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    i__1 = nlev;
		    for (k = 1; k <= i__1; ++k) {
			do_fio(&c__1, (char *)&dn_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    i__3 = nlev;
		    for (k = 1; k <= i__3; ++k) {
			do_fio(&c__1, (char *)&up_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    e_wsfe();
		}
	    }

	    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {

/* ****         print the transmission to space and the p(tau=1) */

		if (nz == 1) {
		    io___36.ciunit = *iutrn;
		    s_wsfe(&io___36);
		    do_fio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&trn_ray_0__, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*p_ray_0__), (ftnlen)sizeof(real))
			    ;
		    i__1 = *numu;
		    for (nze = *nzup; nze <= i__1; ++nze) {
			do_fio(&c__1, (char *)&p_ray0__[nze], (ftnlen)sizeof(
				real));
		    }
		    do_fio(&c__1, (char *)&trn_gas_0__, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*p_gas_0__), (ftnlen)sizeof(real))
			    ;
		    i__3 = *numu;
		    for (nze = *nzup; nze <= i__3; ++nze) {
			do_fio(&c__1, (char *)&p_gas0__[nze], (ftnlen)sizeof(
				real));
		    }
		    do_fio(&c__1, (char *)&trn_aer_0__, (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&(*p_aer_0__), (ftnlen)sizeof(real))
			    ;
		    i__2 = *numu;
		    for (nze = *nzup; nze <= i__2; ++nze) {
			do_fio(&c__1, (char *)&p_aer0__[nze], (ftnlen)sizeof(
				real));
		    }
		    e_wsfe();
		}
	    }

	} else {

	    if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {

/* ****         print level-dependent fluxes at the wavelength */

		if (*lsolar) {
		    io___37.ciunit = *iu_flux__;
		    s_wsue(&io___37);
		    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    i__1 = nlev;
		    for (k = 1; k <= i__1; ++k) {
			do_uio(&c__1, (char *)&dn_dir_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    i__3 = nlev;
		    for (k = 1; k <= i__3; ++k) {
			do_uio(&c__1, (char *)&dn_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    i__2 = nlev;
		    for (k = 1; k <= i__2; ++k) {
			do_uio(&c__1, (char *)&up_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    e_wsue();
		} else {
		    io___38.ciunit = *iu_flux__;
		    s_wsue(&io___38);
		    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    i__1 = nlev;
		    for (k = 1; k <= i__1; ++k) {
			do_uio(&c__1, (char *)&dn_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    i__3 = nlev;
		    for (k = 1; k <= i__3; ++k) {
			do_uio(&c__1, (char *)&up_diff_flx__[k - 1], (ftnlen)
				sizeof(real));
		    }
		    e_wsue();
		}
	    }

	    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {

/* ****         print the transmission to space and the p(tau=1) */

		if (nz == 1) {
		    io___39.ciunit = *iutrn;
		    s_wsue(&io___39);
		    do_uio(&c__1, (char *)&wl, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&wnio, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&trn_ray_0__, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&(*p_ray_0__), (ftnlen)sizeof(real))
			    ;
		    i__1 = *numu;
		    for (nze = *nzup; nze <= i__1; ++nze) {
			do_uio(&c__1, (char *)&p_ray0__[nze], (ftnlen)sizeof(
				real));
		    }
		    do_uio(&c__1, (char *)&trn_gas_0__, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&(*p_gas_0__), (ftnlen)sizeof(real))
			    ;
		    i__3 = *numu;
		    for (nze = *nzup; nze <= i__3; ++nze) {
			do_uio(&c__1, (char *)&p_gas0__[nze], (ftnlen)sizeof(
				real));
		    }
		    do_uio(&c__1, (char *)&trn_aer_0__, (ftnlen)sizeof(real));
		    do_uio(&c__1, (char *)&(*p_aer_0__), (ftnlen)sizeof(real))
			    ;
		    i__2 = *numu;
		    for (nze = *nzup; nze <= i__2; ++nze) {
			do_uio(&c__1, (char *)&p_aer0__[nze], (ftnlen)sizeof(
				real));
		    }
		    e_wsue();
		}
	    }
	}

    } else {

/* ****      call slit function program */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {

/* ****         compute the unit number */

	    iurad = *iuout + (nz - 1) * *nlout + k - 1;
	    no = (nz - 1) * *nlout + k;
	    points[no] += 1.f;

/* ****         load variables into slit arrays */

	    nn = 1;
	    spect[nn - 1] = var_1__;
	    if (*levout == 2 || *levout >= 4 || *levout == 3 && k > 1) {
		++nn;
		spect[nn - 1] = dir_flx0__[k - 1];
		++nn;
		spect[nn - 1] = dn_flx0__[k - 1];
	    }
	    ++nn;
	    spect[nn - 1] = up_flx0__[k - 1];
	    i__3 = *nphi;
	    for (naz = 1; naz <= i__3; ++naz) {
		i__2 = *numu;
		for (nze = nze_1__[k - 1]; nze <= i__2; ++nze) {
		    ++nn;
		    spect[nn - 1] = rad0[nze + (naz + (k << 4) << 4) - 273];
/* L2001: */
		}
/* L2011: */
	    }

	    iord = 2;
	    iscwt = 0;

	    rad_slit__(&iurad, &no, ifrmout, islit, &iscwt, &iord, &nn, wnmin,
		     wnmax, dwn, width, &points[1], wn, spect, &iquit_rad__[
		    nz - 1]);
/* L2021: */
	}

	if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	    if (nz == 1) {

/* *****         add pressure levels of optical depth unity and column- */
/*              integrated transmission of gases and aerosols. */

		points[ntrn] += 1.f;
		nn = 1;
		spect[nn - 1] = trn_ray_0__;
		++nn;
		spect[nn - 1] = *p_ray_0__;
		i__1 = *numu;
		for (nze = *nzup; nze <= i__1; ++nze) {
		    ++nn;
		    spect[nn - 1] = p_ray0__[nze];
/* L3001: */
		}
		++nn;
		spect[nn - 1] = trn_gas_0__;
		++nn;
		spect[nn - 1] = *p_gas_0__;
		i__1 = *numu;
		for (nze = *nzup; nze <= i__1; ++nze) {
		    ++nn;
		    spect[nn - 1] = p_gas0__[nze];
/* L3021: */
		}
		++nn;
		spect[nn - 1] = trn_aer_0__;
		++nn;
		spect[nn - 1] = *p_aer_0__;
		i__1 = *numu;
		for (nze = *nzup; nze <= i__1; ++nze) {
		    ++nn;
		    spect[nn - 1] = p_aer0__[nze];
/* L3041: */
		}

		iord = 2;
		iscwt = 0;

		rad_slit__(iutrn, &ntrn, ifrmout, islit, &iscwt, &iord, &nn, 
			wnmin, wnmax, dwn, width, &points[1], wn, spect, &
			iquit_trn__);

	    }
	}

/* ****      call slit function program for level-dependent net fluxes */

	if (*irad == 3 || *irad == 4 || *irad == 7 || *irad == 8) {

/* ****        increment counting variables */

	    points[nfl] += 1.f;

	    if (*lsolar) {
		nn = nlev * 3;
		i__1 = nlev;
		for (k = 1; k <= i__1; ++k) {
		    spect[k - 1] = dn_dir_flx__[k - 1];
		    spect[k + nlev - 1] = dn_diff_flx__[k - 1];
		    spect[k + (nlev << 1) - 1] = up_diff_flx__[k - 1];
/* L4021: */
		}
	    } else {
		nn = nlev << 1;
		i__1 = nlev;
		for (k = 1; k <= i__1; ++k) {
		    spect[k - 1] = dn_diff_flx__[k - 1];
		    spect[k + nlev - 1] = up_diff_flx__[k - 1];
/* L4041: */
		}
	    }

	    iord = 2;
	    iscwt = 0;

	    rad_slit__(iu_flux__, &nfl, ifrmout, islit, &iscwt, &iord, &nn, 
		    wnmin, wnmax, dwn, width, &points[1], wn, spect, &
		    iquit_flx__[nz - 1]);

	}

    }

    return 0;
} /* rad_out__ */

