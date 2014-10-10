/* int_rad_lev.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int int_rad_lev__(logical *lsolar, logical *lplanck, integer 
	*nphi, integer *numu, integer *nz, integer *nlyr, integer *nlout, 
	integer *levout, integer *k_out__, real *dp_dp__, real *upsflx, real *
	dnsflx, real *dirsflx, real *uptflx, real *dntflx, real *sol_rad__, 
	real *th_rad__, real *ups, real *dns, real *dirs, real *upth, real *
	dnth, real *rad_s__, real *rad_th__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, naz, nze;


/* ccccccccccccccccccccccc  i n t _ r a d _ l e v  ccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine interpolates radiances and fluxes to their      cc */
/* c    output levels, k_out.                                           cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       numu - number of zenith angles in input file                 cc */
/* c     uptflx - wn-dependent upward thermal flux                      cc */
/* c     dntflx - wn-dependent downward thermal flux                    cc */
/* c     th_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c     upsflx - wn-dependent upward solar flux                        cc */
/* c     dnsflx - wn-dependent downward diffuse + direct solar flux     cc */
/* c    dirsflx - wn-dependent downward direct solar flux               cc */
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    uptflx0 - wn-dependent upward thermal flux                      cc */
/* c    dntflx0 - wn-dependent downward thermal flux                    cc */
/* c    th_rad0 - wn-depndent, angle-dependent thermal radiances        cc */
/* c    upsflx0 - wn-dependent upward solar flux                        cc */
/* c    dnsflx0 - wn-dependent downward diffuse+direct solar flux       cc */
/* c   dirsflx0 - wn-dependent downward direct solar flux               cc */
/* c   sol_rad0 - wn-depndent, angle-dependent solar radiances          cc */
/* c              output stream                                         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccc  i n t _ r a d _ l e v  ccccccccccccccccccccccc */




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





/* *****   fluxes interpolated to the output levels */


/* ****    output fluxes and radiances */



    /* Parameter adjustments */
    rad_th__ -= 273;
    rad_s__ -= 273;
    --dnth;
    --upth;
    --dirs;
    --dns;
    --ups;
    th_rad__ -= 273;
    sol_rad__ -= 273;
    --dntflx;
    --uptflx;
    --dirsflx;
    --dnsflx;
    --upsflx;
    --dp_dp__;
    --k_out__;

    /* Function Body */
    if (*lsolar) {

/* ****        interpolate solar radiances and fluxes to output levels */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {
	    if (*levout <= 3) {

/* ****             load fluxes at the appropriate level */

		ups[k] = upsflx[k_out__[k]];
		dns[k] = dnsflx[k_out__[k]];
		dirs[k] = dirsflx[k_out__[k]];

/* ****             load solar radiances at the appropriate level */

		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__3 = *numu;
		    for (nze = 1; nze <= i__3; ++nze) {
			rad_s__[nze + (naz + (k << 4) << 4)] = sol_rad__[nze 
				+ (naz + (k_out__[k] << 4) << 4)];
/* L4201: */
		    }
/* L4211: */
		}

	    } else {

/* ****             interpolate fluxes to the output level */

		if (k_out__[k] != 1 && k_out__[k] != *nlyr + 1) {
		    ups[k] = upsflx[k_out__[k]] - dp_dp__[k] * (upsflx[
			    k_out__[k]] - upsflx[k_out__[k] - 1]);
		    dns[k] = dnsflx[k_out__[k]] - dp_dp__[k] * (dnsflx[
			    k_out__[k]] - dnsflx[k_out__[k] - 1]);
		    dirs[k] = dirsflx[k_out__[k]] - dp_dp__[k] * (dirsflx[
			    k_out__[k]] - dirsflx[k_out__[k] - 1]);

/* ****               interpolate solar radiances to the output level */

		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__3 = *numu;
			for (nze = 1; nze <= i__3; ++nze) {
			    rad_s__[nze + (naz + (k << 4) << 4)] = sol_rad__[
				    nze + (naz + (k_out__[k] << 4) << 4)] - 
				    dp_dp__[k] * (sol_rad__[nze + (naz + (
				    k_out__[k] << 4) << 4)] - sol_rad__[nze + 
				    (naz + (k_out__[k] - 1 << 4) << 4)]);
/* L4231: */
			}
/* L4241: */
		    }

		} else {

/* ****              load solar fluxes and radiances into the */
/*                  top and/or bottlm level */

		    ups[k] = upsflx[k_out__[k]];
		    dns[k] = dnsflx[k_out__[k]];
		    dirs[k] = dirsflx[k_out__[k]];

/* ****               load solar radiances at the appropriate level */

		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__3 = *numu;
			for (nze = 1; nze <= i__3; ++nze) {
			    rad_s__[nze + (naz + (k << 4) << 4)] = sol_rad__[
				    nze + (naz + (k_out__[k] << 4) << 4)];
/* L4251: */
			}
/* L4261: */
		    }
		}

	    }
/* L4281: */
	}

    }

/* ****       t h e r m a l    f l u x e s    a n d    r a d i a n c e s */
/*           (do only for nz = 1) */

    if (*lplanck && *nz == 1) {

/* ****        interpolate thermal radiances and fluxes to output levels */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {

	    if (*levout <= 3) {

/* ****             load thermal fluxes to the appropriate level */

		upth[k] = uptflx[k_out__[k]];
		dnth[k] = dntflx[k_out__[k]];

/* ****           interpolate thermal radiances */

		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__3 = *numu;
		    for (nze = 1; nze <= i__3; ++nze) {
			rad_th__[nze + (naz + (k << 4) << 4)] = th_rad__[nze 
				+ (naz + (k_out__[k] << 4) << 4)];
/* L4401: */
		    }
/* L4411: */
		}

	    } else {

/* ****         interpolate fluxes and radiances to the output level */

		if (k_out__[k] != 1 && k_out__[k] != *nlyr + 1) {

		    upth[k] = uptflx[k_out__[k]] - dp_dp__[k] * (uptflx[
			    k_out__[k]] - uptflx[k_out__[k] - 1]);
		    dnth[k] = dntflx[k_out__[k]] - dp_dp__[k] * (dntflx[
			    k_out__[k]] - dntflx[k_out__[k] - 1]);

/* ****               interpolate thermal radiances */

		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__3 = *numu;
			for (nze = 1; nze <= i__3; ++nze) {
			    rad_th__[nze + (naz + (k << 4) << 4)] = th_rad__[
				    nze + (naz + (k_out__[k] << 4) << 4)] - 
				    dp_dp__[k] * (th_rad__[nze + (naz + (
				    k_out__[k] << 4) << 4)] - th_rad__[nze + (
				    naz + (k_out__[k] - 1 << 4) << 4)]);
/* L4421: */
			}
/* L4431: */
		    }

		} else {

/* ****              load thermal fluxes and radiances into the */
/*                  top and/or bottom level */

		    upth[k] = uptflx[k_out__[k]];
		    dnth[k] = dntflx[k_out__[k]];

		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__3 = *numu;
			for (nze = 1; nze <= i__3; ++nze) {
			    rad_th__[nze + (naz + (k << 4) << 4)] = th_rad__[
				    nze + (naz + (k_out__[k] << 4) << 4)];
/* L4441: */
			}
/* L4451: */
		    }

		}
	    }
/* L4461: */
	}

    }

    return 0;
} /* int_rad_lev__ */

