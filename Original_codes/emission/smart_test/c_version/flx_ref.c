/* flx_ref.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int flx_ref__(integer *l0, integer *nlyr, integer *ng0, real 
	*alb0, doublereal *trnflx_b__, doublereal *refflx_b__, doublereal *
	absflx_b__, real *dx_b_i__, doublereal *flx_dd__, doublereal *
	flx_du__, doublereal *flx_urt__, doublereal *flx_drb__, doublereal *
	flx_rb__, doublereal *flx_rt__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k, l, ibn, nlev;


/* ccccccccccccccccccccccccccc  f l x _ r e f   cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine computes the effective layer radiance and flux  cc */
/* c    transmittances, reflectances and and absorptances for each      cc */
/* c    spectral bin.  These quantities are used in the spectral        cc */
/* c    mapping and jacobian calculations.                              cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c        nlyr: number of layers on the atmosphere.                   cc */
/* c         ng0: spectral bin number                                   cc */
/* c      dtau_b: layer optical depth for bin, ibn0.                    cc */
/* c     copi0_b: layer single scattering albedo for bin, ibn0.         cc */
/* c      pmom_b: layer particle phase function for bin, ibn0.          cc */
/* c         umu: zenith angle of each stream.                          cc */
/* c         gwt: gaussian weight for each stream.                      cc */
/* c    trn_rad: layer transmittance for each radiance stream.          cc */
/* c    abs_rad: layer absorptance for each rediance stream             cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c      flx_rt(k): reflectance of an inhomogeneous layer extending    cc */
/* c             from level 1 to level k+1 (layers 1 - k)               cc */
/* c      flx_rb(k): reflectance of an inhomogeneous layer extending    cc */
/* c             from the surface (nlev) to level k                     cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  f l x _ r e f   cccccccccccccccccccccccccc */




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





/*      to DISORT for the flux transmission calculation. */


/* ****   binned layer flux transmittances and absorptances */



/* ****  find layer-integrated reflectivities.  the quantities are: */
/*      flx_rt(k): reflectance at the bottom of an inhomogeneous layer */
/*             extending from level 1 to level k+1 (layers 1 - k) */
/*      flx_dd(k): denominator of flx_rt(k) */
/*      flx_rb(k): reflectance of the top of an inhomogeneous layer */
/*             extending from the surface (nlev) to level k */
/*      flx_du(k): denominator of flx_rb(k) */

    /* Parameter adjustments */
    flx_rt__ -= 71;
    flx_rb__ -= 71;
    flx_drb__ -= 71;
    flx_urt__ -= 71;
    flx_du__ -= 71;
    flx_dd__ -= 71;
    dx_b_i__ -= 421;
    absflx_b__ -= 421;
    refflx_b__ -= 421;
    trnflx_b__ -= 421;

    /* Function Body */
    l = *l0;
    ibn = *ng0;
    nlev = *nlyr + 1;

    flx_rt__[l * 70 + 1] = 0.f;
    flx_rb__[nlev + l * 70] = refflx_b__[nlev + (l + ibn * 5) * 70];

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {

/* ****       find the layer reflectances adding from top */

	flx_dd__[k + l * 70] = 1. / (1. - flx_rt__[k + l * 70] * refflx_b__[k 
		+ (l + ibn * 5) * 70]);
	flx_urt__[k + l * 70] = trnflx_b__[k + (l + ibn * 5) * 70] * (
		refflx_b__[k + (l + ibn * 5) * 70] * flx_rt__[k + l * 70] * 
		flx_dd__[k + l * 70] + 1.f);
/* Computing 2nd power */
	d__1 = trnflx_b__[k + (l + ibn * 5) * 70];
	flx_rt__[k + 1 + l * 70] = refflx_b__[k + (l + ibn * 5) * 70] + d__1 *
		 d__1 * flx_rt__[k + l * 70] * flx_dd__[k + l * 70];

/* ***        find layer reflectances, adding up from the surface */

	flx_du__[nlev - k + l * 70] = 1. / (1. - refflx_b__[nlev - k + (l + 
		ibn * 5) * 70] * flx_rb__[nlev - k + 1 + l * 70]);
	flx_drb__[nlev - k + 1 + l * 70] = trnflx_b__[nlev - k + (l + ibn * 5)
		 * 70] * (refflx_b__[nlev - k + (l + ibn * 5) * 70] * 
		flx_rb__[nlev - k + 1 + l * 70] * flx_du__[nlev - k + l * 70] 
		+ 1.f);
/* Computing 2nd power */
	d__1 = trnflx_b__[nlev - k + (l + ibn * 5) * 70];
	flx_rb__[nlev - k + l * 70] = refflx_b__[nlev - k + (l + ibn * 5) * 
		70] + d__1 * d__1 * flx_rb__[nlev - k + 1 + l * 70] * 
		flx_du__[nlev - k + l * 70];
/* L1001: */
    }

/*      if(ibn .eq. 1) then */
/*      write(*,'(/,1a,2i5,1pe12.4)') 'in flx_ref: ibn,l,alb0 ',ibn,l,alb0 */
/*      write(*,'(1a,16(1pe12.4))') 'refflx: ', */
/*     -                             (refflx_b(k,l,ibn),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'trnflx: ', */
/*     -                             (trnflx_b(k,l,ibn),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_dd: ',(flx_dd(k,l),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_urt:',(flx_urt(k,l),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_du: ',(flx_du(k,l),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_drb:',(flx_drb(k,l),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_rt: ',(flx_rt(k,l),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_rb: ',(flx_rb(k,l),k=1,nlev) */
/*      endif */

    return 0;
} /* flx_ref__ */

