/* gas_tau.f -- translated by f2c (version 20100827).
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

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__4 = 4;

/* Subroutine */ int gas_tau__(integer *iuabc, integer *igtrn, integer *igpab,
	 integer *ngases, integer *ngastyp, integer *npg, integer *ntg, 
	integer *nbgas, integer *ne, integer *nlev, integer *nzup, integer *
	nt_pd__, integer *ntaug, integer *nstr, integer *iflext, integer *
	nsiext, integer *nlay, integer *ntau_pd0__, integer *istate, 
	doublereal *wn, doublereal *wnmin, doublereal *wnmax, real *ratm, 
	real *umu, real *p, real *t, real *alt, real *grav, real *z__, real *
	rmix, real *pd_frac__, doublereal *wnext, real *dtauex, real *
	tau_ext__, real *taugas, real *p_gas__, real *p_gas_0__, doublereal *
	wn_eof__, integer *io_end__, integer *io_err__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static real dtau_tst__;
    static integer k, l, m;
    extern /* Subroutine */ int gaslbltau_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, 
	    doublereal *, integer *, integer *);
    static integer ng, ni;
    static real dv;
    static integer ng0, nee, ngt, nze, ngt0, nlyr;
    static real dtau0;
    static integer ktau1;
    static real dtaug[94920]	/* was [70][2][113][2][3] */;
    static doublereal wneof[2];
    static real d_temp__[2], d_pres__[2], d_rmix__[100]	/* was [2][50] */, 
	    dtaugas[420]	/* was [70][2][3] */;
    extern /* Subroutine */ int gasetau_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, doublereal *, integer *, integer *);
    static real dtaugdv[47460]	/* was [70][113][2][3] */, ptaugas;
    extern /* Subroutine */ int gasptau_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 97, 0, 0, 0 };



/* cccccccccccccccccccccccccc  g a s _  t a u cccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine finds the gas absorption coefficients at each   cc */
/* c    model level, and computes the effective monochromatic gas       cc */
/* c    absorption optical depth for each atmospheric layer.            cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      iuabc - unit numbers for gas line parameters                  cc */
/* c       nlyr - number of layers in the model atmosphere (nlev-1).    cc */
/* c       nlev - number of levels in the model atmosphere (<kp).       cc */
/* c     ngases - number of absorbing gases included in model.          cc */
/* c         ne - extinction/source index.                              cc */
/* c     iflext - extinction flag: 0) do not read next spectral value,  cc */
/* c             1) read next spectral value.                           cc */
/* c      igtrn - type of energy transition: (1) line,                  cc */
/* c             (2) pressure-induced, (3) electronic                   cc */
/* c       ratm - mean atmospheric gas constant (j/kg/k).               cc */
/* c       nzup - index of first upward radiance stream.                cc */
/* c       nstr - number of computational zenith angles (streams).      cc */
/* c        umu - cosine of each upward radiance stream.                cc */
/* c         wn - current wavenumber (cm**-1).                          cc */
/* c      ntaug - number of optical depths needed for this gas:         cc */
/* c              =1 if this gas is not a variable part of state vector cc */
/* c              =3 if this gas is a variable part of state vector     cc */
/* c              (nominal, rmix perturbed at top of layer, rmix        cc */
/* c               perturbed at bottom of layer). */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c     io_end: end of file flag for each wavelength constituent:      cc */
/* c             (if io_end ne 0, the end of file has been encountered) cc */
/* c     io_err: i/o error flag  for each wavelength constituent:       cc */
/* c             (if io_err ne 0, an i/o error has been encountered)    cc */
/* c     wn_eof: last wavenumber in each input file.                    cc */
/* c     taugas - column-integrated gas optical depth.                  cc */
/* c      dtaug - monochromatic gas optical depth in each layer.        cc */
/* c     dtauex - monochromatic extinction optical depth in each layer. cc */
/* c    dtaugdv - derivative of dtaug with respect to wavenumber.       cc */
/* c      p_gas - pressure level of gas tau=1 for each upward stream    cc */
/* c              (bars).                                               cc */
/* c    p_gas_0 - pressure of gas tau=1 for a vertical stream.          cc */
/* c     d_pres - perturbed surface pressure used when pressure is a    cc */
/* c              variable part of the state vector                     cc */
/* c     d_temp - fractional temperature perturbation used when         cc */
/* c              temperature is a variable componen of the state       cc */
/* c              vector.                                               cc */
/* c     d_rmix - perturbed gas mixing ratio used when rmix is a        cc */
/* c              variable part of the state vector                     cc */
/* c                                                                    cc */
/* cccccccccccccccccccccccccc  g a s _  t a u cccccccccccccccccccccccccccc */




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
    --io_err__;
    --io_end__;
    wn_eof__ -= 3;
    --p_gas__;
    --taugas;
    tau_ext__ -= 281;
    dtauex -= 71;
    wnext -= 3;
    --pd_frac__;
    rmix -= 71;
    --z__;
    --grav;
    --alt;
    --t;
    --p;
    --umu;
    --istate;
    nsiext -= 3;
    --iflext;
    nbgas -= 4;
    ntg -= 4;
    npg -= 4;
    --ngastyp;
    igpab -= 4;
    igtrn -= 4;
    iuabc -= 4;

    /* Function Body */
    nlyr = *nlev - 1;

/* *****  set the wavenumber start counter */

    wneof[0] = wn_eof__[(*ne << 1) + 1];

/* *****   determine whether pressure is a variable part of */
/*        the state vector.  if it is, gas absorption must be */
/*        computed for the background and perterbed pressure profile */

    d_pres__[0] = 1.f;
    if (istate[1] != 0) {
	d_pres__[1] = pd_frac__[1] + 1.f;
    }

/* *****   determine whether temperature is a variable part of */
/*        the state vector.  if it is, gas absorption must be */
/*        computed for the background and perterbed temperature profile */

    d_temp__[0] = 1.f;
    if (istate[2] != 0) {
	d_temp__[1] = pd_frac__[2] + 1.f;
    }

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	taugas[k] = 0.f;
/* L1001: */
    }
    if (istate[1] != 0) {
	taugas[*nlev + 1] = 0.f;
    }

/* Danie Liang, 10/05/2007 */
/*      if(wn .eq. wnmin) then */
    if (*wn < *wnmin) {
	for (m = 1; m <= 3; ++m) {
	    for (ni = 1; ni <= 2; ++ni) {
		for (nee = 1; nee <= 113; ++nee) {
		    i__1 = *nlev;
		    for (k = 1; k <= i__1; ++k) {
			dtaugdv[k + (nee + (ni + (m << 1)) * 113) * 70 - 
				23801] = 0.f;
/* L1101: */
		    }
		    for (l = 1; l <= 2; ++l) {
			i__1 = *nlev;
			for (k = 1; k <= i__1; ++k) {
			    dtaug[k + (l + (nee + (ni + (m << 1)) * 113 << 1))
				     * 70 - 47671] = 0.f;
/* L1121: */
			}
/* L1141: */
		    }
/* L1161: */
		}
/* L1181: */
	    }
/* L1191: */
	}

    }

/* ****   enter the loop over absorbing gases */

    i__1 = *ngases;
    for (ng = 1; ng <= i__1; ++ng) {
	ng0 = ng;

/* *****         determine whether this gas is a variable element of */
/*              the state vector.  If it is, update extinction counter */
/*              used in calculating partial derivatives */

	*ntaug = 1;
	d_rmix__[(ng << 1) - 2] = 1.f;
	if ((i__2 = istate[ng + 2], abs(i__2)) == 3) {
	    ++(*ntau_pd0__);
	    *ntaug = 3;
	    d_rmix__[(ng << 1) - 1] = pd_frac__[ng + 2] + 1.f;

	    for (m = 1; m <= 3; ++m) {
		i__2 = *nlay;
		for (k = 1; k <= i__2; ++k) {
		    tau_ext__[k + (m + *ntau_pd0__ * 3) * 70] = 0.f;
/* L1201: */
		}
/* L1221: */
	    }

	}

/* ****        update extinction index and check extinction flag */

	i__2 = ngastyp[ng];
	for (ngt = 1; ngt <= i__2; ++ngt) {
	    ngt0 = ngt;
	    ++(*ne);
	    if (iflext[*ne] == 1) {

/* ****              determine the type of transition: (1) line, */
/*                  (2) pressure-induced, (3) electronic */

		if (igtrn[ngt + ng * 3] == 1) {

/* ****              read the next line absorption coefficient */

		    gaslbltau_(ne, &iuabc[4], &ng0, &ngt0, nlev, &npg[4], &
			    ntg[4], &nbgas[4], &nsiext[3], nlay, &istate[1], 
			    nt_pd__, ntaug, wn, wnmin, wnmax, &wnext[3], &p[1]
			    , &t[1], &z__[1], &grav[1], &rmix[71], d_pres__, 
			    d_temp__, d_rmix__, dtaug, dtaugdv, wneof, &
			    io_end__[1], &io_err__[1]);

		} else {

		    if (igtrn[ngt + ng * 3] == 2) {

/* ****                 read next p-induced transition */

			gasptau_(ne, &iuabc[4], &ng0, &ngt0, &igpab[4], nlev, 
				&nsiext[3], nlay, &istate[1], nt_pd__, ntaug, 
				wn, wnmin, wnmax, &wnext[3], &p[1], &t[1], &
				alt[1], &rmix[71], d_pres__, d_temp__, 
				d_rmix__, dtaug, dtaugdv, wneof, &io_end__[1],
				 &io_err__[1]);

		    } else {

/* ****                 read next electronic transition cross-section */

			gasetau_(ne, &iuabc[4], &ng0, &ngt0, nlev, nlay, &
				nsiext[3], nt_pd__, &istate[1], ntaug, wn, 
				wnmin, wnmax, &wnext[3], &p[1], &t[1], &alt[1]
				, &rmix[71], ratm, d_pres__, d_temp__, 
				d_rmix__, dtaug, dtaugdv, wneof, &io_end__[1],
				 &io_err__[1]);

		    }
		}

/* ****             turn off the extinction flag for this gas */

		iflext[*ne] = 0;

	    }
	    wn_eof__[(*ne << 1) + 1] = wneof[0];
	    wn_eof__[(*ne << 1) + 2] = wneof[1];

/* ****          interpolate gas optical depths to wavenumber wn */

	    dv = (real) (*wn - wnext[nsiext[(*ne << 1) + 1] + (*ne << 1)]);

	    i__3 = *ntaug;
	    for (m = 1; m <= i__3; ++m) {
		i__4 = *nt_pd__;
		for (l = 1; l <= i__4; ++l) {
		    i__5 = *nlay;
		    for (k = 1; k <= i__5; ++k) {
			dtau0 = dtaug[k + (nsiext[(*ne << 1) + 1] + (*ne + (l 
				+ (m << 1)) * 113 << 1)) * 70 - 47671] + 
				dtaugdv[k + (*ne + (l + (m << 1)) * 113) * 70 
				- 23801] * dv;
			if (dtau0 > 0.f) {
			    dtaugas[k + (l + (m << 1)) * 70 - 211] = dtau0;
			} else {
			    dtaugas[k + (l + (m << 1)) * 70 - 211] = 0.f;
			}
/* L1401: */
		    }
/* L1421: */
		}
/* L1441: */
	    }
	    i__3 = *nt_pd__;
	    for (l = 1; l <= i__3; ++l) {
		i__4 = *nlay;
		for (k = 1; k <= i__4; ++k) {
		    dtauex[k + l * 70] += dtaugas[k + (l + 2) * 70 - 211];
/* L1461: */
		}
/* L1481: */
	    }

	    i__3 = *nlay;
	    for (k = 1; k <= i__3; ++k) {
		taugas[k + 1] += dtaugas[k - 1];
/* L1501: */
	    }

/* ****                if this is a variable component of the state */
/*                    vector, save differential extinction optical depth */

	    if ((i__3 = istate[ng + 2], abs(i__3)) == 3) {

		i__3 = *ntaug;
		for (m = 1; m <= i__3; ++m) {
		    i__4 = *nlay;
		    for (k = 1; k <= i__4; ++k) {
			tau_ext__[k + (m + *ntau_pd0__ * 3) * 70] += dtaugas[
				k + ((m << 1) + 1) * 70 - 211];
/* L1601: */
		    }

/* L1621: */
		}

	    }

/* L1701: */
	}

/* L1801: */
    }

/* ****      find the layer of gas optical absorption optical depth */
/*          unity at each of the output emission angles */

    i__1 = *nstr;
    for (nze = *nzup; nze <= i__1; ++nze) {
	ptaugas = 0.f;
	ktau1 = 0;

	i__2 = nlyr;
	for (k = 1; k <= i__2; ++k) {
	    if (ptaugas == 0.f && taugas[k + 1] / umu[nze] >= 1.f) {
		ktau1 = k + 1;
		ptaugas = p[k + 1];
	    }
/* L2001: */
	}

/* ****       find the pressure of gas absorption optical depth unity: */

	if (ktau1 > 1) {

	    dtau_tst__ = (taugas[ktau1] - taugas[ktau1 - 1]) / umu[nze];
	    if (dtau_tst__ != 0.f) {
		p_gas__[nze] = (p[ktau1] - (taugas[ktau1] / umu[nze] - 1.f) * 
			(p[ktau1] - p[ktau1 - 1]) / dtau_tst__) * 1e-5f;
	    } else {
		p_gas__[nze] = p[ktau1] * 1e-5f;
	    }

	} else {

/* ****        tau < 1 in column.  set to surface value */

	    ktau1 = *nlev;
	    p_gas__[nze] = p[ktau1] * 1e-5f;

	}
/* L2041: */
    }

/* ****       find the pressure of tau = 1 for a vertical path */

    ptaugas = 0.f;
    ktau1 = 0;

    i__1 = nlyr;
    for (k = 1; k <= i__1; ++k) {
	if (ptaugas == 0.f && taugas[k + 1] >= 1.f) {
	    ktau1 = k + 1;
	    ptaugas = p[k + 1];
	}
/* L3001: */
    }

/* ****      find the pressure of gas absorption optical depth unity: */
/*          note: ktau can't equal to 1.0 because taugas(1) = 0.0 */

    if (ktau1 > 1) {

	dtau_tst__ = taugas[ktau1] - taugas[ktau1 - 1];
	if (dtau_tst__ != 0.f) {
	    *p_gas_0__ = (p[ktau1] - (taugas[ktau1] - 1.f) * (p[ktau1] - p[
		    ktau1 - 1]) / dtau_tst__) * 1e-5f;
	} else {
	    ktau1 = *nlev;
	    *p_gas_0__ = p[ktau1] * 1e-5f;
	}

    } else {

/* ****      tau < 1 in column.  set to surface value */

	*p_gas_0__ = p[*nlev] * 1e-5f;

    }

/* Danie Liang */
    i__1 = *nlay;
    for (k = 1; k <= i__1; ++k) {
	s_wsle(&io___24);
	do_lio(&c__5, &c__1, (char *)&(*wn), (ftnlen)sizeof(doublereal));
	do_lio(&c__4, &c__1, (char *)&p[k], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&taugas[k + 1], (ftnlen)sizeof(real));
	e_wsle();
    }
    return 0;
} /* gas_tau__ */

