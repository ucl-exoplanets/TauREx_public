/* init_interp.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int init_interp__(integer *nlev, integer *nlout, integer *
	numu, integer *nphi, integer *nz0, doublereal *wnmin, doublereal *
	wnout, real *upsflx, real *dnsflx, real *dirsflx, real *sol_rad__, 
	real *uptflx, real *dntflx, real *th_rad__, doublereal *up_s_src_i__, 
	doublereal *dn_s_src_i__, doublereal *dup_s_src__, doublereal *
	ddn_s_src__, doublereal *up_t_src_i__, doublereal *dn_t_src_i__, 
	doublereal *dup_t_src__, doublereal *ddn_t_src__, real *upsflx_i__, 
	real *dnsflx_i__, real *dirsflx_i__, real *dupsflx, real *ddnsflx, 
	real *ddirsflx, real *sol_rad0__, real *dsol_rad__, real *uptflx_i__, 
	real *dntflx_i__, real *duptflx, real *ddntflx, real *th_rad0__, real 
	*dth_rad__, real *ups, real *dns, real *dirs, real *upth, real *dnth, 
	real *rad_s__, real *rad_th__, real *pray0, real *pgas0, real *paer0, 
	real *pray_0__, real *pgas_0__, real *paer_0__, real *tau_ray__, real 
	*tau_gas__, real *tau_aer__, real *dpray, real *dpgas, real *dpaer, 
	real *tau_ray0__, real *tau_gas0__, real *tau_aer0__, real *pray_00__,
	 real *pgas_00__, real *paer_00__, real *dtau_ray__, real *dtau_gas__,
	 real *dtau_aer__, real *dpray0, real *dpgas0, real *dpaer0, real *
	pgas, real *pray, real *paer)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, ne, ni, nz, naz, nze;


/* ccccccccccccccccccccccccccc  init_interp  ccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine initializes extinction flags and counters and   cc */
/* c    variables that record the wavelength range for each wavelength  cc */
/* c    dependent input file.                                           cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    ngases : number of absorbing gases included in calculation.     cc */
/* c    nmodes : number of aerosol particle modes used in calculation.  cc */
/* c                                                                    cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  init_interp  ccccccccccccccccccccccccccccc */




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





/* ****     monochormatic fluxes at wavenumber wnout */


/* ****   downward and upward solar fluxes gradients at wavenumber wnout */


/* ****   monochormatic radiances and wavelength gradients at wnout */
/*       at output levels, k_out */


/* ****     monochormatic flux source interpolation variables */


/* ****   wavelength gradients in source terms */



/* ****    variables for the output wavenumber grid */





    /* Parameter adjustments */
    --paer;
    --pray;
    --pgas;
    --paer_00__;
    --pgas_00__;
    --pray_00__;
    --tau_aer0__;
    --tau_gas0__;
    --tau_ray0__;
    --dpaer;
    --dpgas;
    --dpray;
    paer0 -= 17;
    pgas0 -= 17;
    pray0 -= 17;
    rad_th__ -= 273;
    rad_s__ -= 273;
    --dnth;
    --upth;
    --dirs;
    --dns;
    --ups;
    dth_rad__ -= 273;
    th_rad0__ -= 1041;
    --ddntflx;
    --duptflx;
    dntflx_i__ -= 71;
    uptflx_i__ -= 71;
    dsol_rad__ -= 1041;
    sol_rad0__ -= 2577;
    ddirsflx -= 71;
    ddnsflx -= 71;
    dupsflx -= 71;
    dirsflx_i__ -= 211;
    dnsflx_i__ -= 211;
    upsflx_i__ -= 211;
    --ddn_t_src__;
    --dup_t_src__;
    dn_t_src_i__ -= 71;
    up_t_src_i__ -= 71;
    ddn_s_src__ -= 71;
    dup_s_src__ -= 71;
    dn_s_src_i__ -= 211;
    up_s_src_i__ -= 211;
    th_rad__ -= 273;
    --dntflx;
    --uptflx;
    sol_rad__ -= 273;
    --dirsflx;
    --dnsflx;
    --upsflx;
    wnout -= 229;

    /* Function Body */
    nz = *nz0;

/* ****   initialize spectral variables */

    for (ne = 1; ne <= 2; ++ne) {
	wnout[(ne + nz * 113 << 1) + 1] = -999.;
	wnout[(ne + nz * 113 << 1) + 2] = *wnmin;
/* L1001: */
    }

/* ****     initialize the thermal fluxes and radiances */

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	upsflx[k] = 0.f;
	dnsflx[k] = 0.f;
	dirsflx[k] = 0.f;
	uptflx[k] = 0.f;
	dntflx[k] = 0.f;
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *numu;
	    for (nze = 1; nze <= i__3; ++nze) {
		sol_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
		th_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L1021: */
	    }
/* L1041: */
	}
/* L1061: */
    }
    i__1 = *nlout;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *numu;
	    for (nze = 1; nze <= i__3; ++nze) {
		rad_s__[nze + (naz + (k << 4) << 4)] = 0.f;
		rad_th__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L1121: */
	    }
/* L1141: */
	}
/* L1161: */
    }

/* ****    initialize the optical depths and pressure of tau=1 */

    if (nz == 1) {
	*tau_gas__ = 0.f;
	*tau_ray__ = 0.f;
	*tau_aer__ = 0.f;
	*pgas_0__ = 0.f;
	*pray_0__ = 0.f;
	*paer_0__ = 0.f;
	*dtau_gas__ = 0.f;
	*dtau_ray__ = 0.f;
	*dtau_aer__ = 0.f;
	*dpgas0 = 0.f;
	*dpray0 = 0.f;
	*dpaer0 = 0.f;

	i__1 = *numu;
	for (nze = 1; nze <= i__1; ++nze) {
	    pgas[nze] = 0.f;
	    pray[nze] = 0.f;
	    paer[nze] = 0.f;
	    dpgas[nze] = 0.f;
	    dpray[nze] = 0.f;
	    dpaer[nze] = 0.f;
/* L1201: */
	}
    }

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	dup_s_src__[k + nz * 70] = 0.;
	ddn_s_src__[k + nz * 70] = 0.;
	dupsflx[k + nz * 70] = 0.f;
	ddnsflx[k + nz * 70] = 0.f;
	ddirsflx[k + nz * 70] = 0.f;
/* L1221: */
    }

    i__1 = *nlout;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *numu;
	    for (nze = 1; nze <= i__3; ++nze) {
		dsol_rad__[nze + (naz + (k + nz * 3 << 4) << 4)] = 0.f;
		dth_rad__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L1321: */
	    }
/* L1341: */
	}
/* L1361: */
    }

    for (ni = 1; ni <= 2; ++ni) {
	tau_gas0__[ni] = 0.f;
	tau_ray0__[ni] = 0.f;
	tau_aer0__[ni] = 0.f;
	pray_00__[ni] = 0.f;
	pgas_00__[ni] = 0.f;
	paer_00__[ni] = 0.f;
	i__1 = *numu;
	for (nze = 1; nze <= i__1; ++nze) {
	    pgas0[nze + (ni << 4)] = 0.f;
	    pray0[nze + (ni << 4)] = 0.f;
	    paer0[nze + (ni << 4)] = 0.f;
	    dpgas[nze] = 0.f;
	    dpray[nze] = 0.f;
	    dpaer[nze] = 0.f;
/* L1441: */
	}

/* ****       initialize the solar fluxes and radiances for this solar */
/*           zenith angle */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    up_s_src_i__[k + (ni + (nz << 1)) * 70] = 0.;
	    dn_s_src_i__[k + (ni + (nz << 1)) * 70] = 0.;
	    upsflx_i__[k + (ni + (nz << 1)) * 70] = 0.f;
	    dnsflx_i__[k + (ni + (nz << 1)) * 70] = 0.f;
	    dirsflx_i__[k + (ni + (nz << 1)) * 70] = 0.f;
/* L1541: */
	}

/* ****     initialize the spectral interpolation variables */

	i__1 = *nlev;
	for (k = 1; k <= i__1; ++k) {
	    uptflx_i__[k + ni * 70] = 0.f;
	    dntflx_i__[k + ni * 70] = 0.f;
/* L1641: */
	}

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {
	    ups[k] = 0.f;
	    dns[k] = 0.f;
	    dirs[k] = 0.f;
	    upth[k] = 0.f;
	    dnth[k] = 0.f;
	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = 1; nze <= i__3; ++nze) {
		    sol_rad0__[nze + (naz + (k + (ni + (nz << 1)) * 3 << 4) <<
			     4)] = 0.f;
		    th_rad0__[nze + (naz + (k + ni * 3 << 4) << 4)] = 0.f;
/* L1701: */
		}
/* L1721: */
	    }
/* L1741: */
	}

/* L1801: */
    }

    return 0;
} /* init_interp__ */

