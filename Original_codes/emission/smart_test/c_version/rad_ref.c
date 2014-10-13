/* rad_ref.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int rad_ref__(integer *l0, integer *nlyr, integer *ng0, 
	integer *numu, integer *nzdn, integer *nzup, integer *nphi, 
	doublereal *trnrad_b__, doublereal *refrad_b__, doublereal *brdf_b__, 
	doublereal *flx_urt__, doublereal *flx_drb__, doublereal *rad_rb__, 
	doublereal *rad_rt__)
{
    /* Initialized data */

    static doublereal pi = 3.141592654;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, l, ibn, naz, nze;
    static doublereal dflx;
    static integer nlev;
    static doublereal uflx;


/* ccccccccccccccccccccccccccc  r a d _ r e f   cccccccccccccccccccccccccc */
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
/* c      rad_rt(k): reflectance of an inhomogeneous layer extendind    cc */
/* c             from level 1 to level k+1 (layers 1 - k)               cc */
/* c      rad_rb(k): reflectance of an inhomogeneous layer extendind    cc */
/* c             from the surface (nlev) to level k                     cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc  r a d _ r e f   cccccccccccccccccccccccccc */




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




/* **** flux transmittances, reflectances, and fluxes */


/* ****   binned layer radiance transmittances and absorptances */




    /* Parameter adjustments */
    rad_rt__ -= 18193;
    rad_rb__ -= 18193;
    flx_drb__ -= 71;
    flx_urt__ -= 71;
    brdf_b__ -= 1553;
    refrad_b__ -= 6737;
    trnrad_b__ -= 6737;

    /* Function Body */

/* ****  find layer-integrated reflectivities.  the quantities are: */
/*      rad_rt(k): reflectance of an inhomogeneous layer extendind from */
/*                 level 1 to level k+1 (layers 1 - k) */
/*      rad_dd(k): denominator of rad_rt(k) */
/*      rad_rb(k): reflectance of an inhomogeneous layer extendind from */
/*                 the surface (nlev) to level k */
/*      rad_du(k): denominator of rad_rb(k) */

    l = *l0;
    ibn = *ng0;
    nlev = *nlyr + 1;

/* ****    set reflectance at the top of the atmosphere for adding down */

    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {
	i__2 = *nzdn;
	for (nze = 1; nze <= i__2; ++nze) {
	    rad_rt__[nze + (naz + (l * 70 + 1 << 4) << 4)] = refrad_b__[nze + 
		    ((l + ibn * 5) * 70 + 1 << 4)];
/* L1001: */
	}
/* L1021: */
    }

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	uflx = flx_urt__[k + l * 70] / pi;
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *nzdn;
	    for (nze = 1; nze <= i__3; ++nze) {
		rad_rt__[nze + (naz + (k + 1 + l * 70 << 4) << 4)] = 
			refrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] + 
			trnrad_b__[nze + (k + (l + ibn * 5) * 70 << 4)] * 
			rad_rt__[nze + (naz + (k + l * 70 << 4) << 4)] * uflx;
/* L1201: */
	    }
/* L1221: */
	}
/* L1241: */
    }
/*      if(ng0 .eq. 1) then */
/*      do 1 nze=1,nzdn */
/*      write(*,'(/,1a,3i5)') 'in rad_ref: nze,naz,l',nze,l */
/*      write(*,'(1a,16(1pe12.4))') 'trnrad_b:  ', */
/*     -         (trnrad_b(nze,k,l,ibn),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'refrad_b:  ', */
/*     -         (refrad_b(nze,k,l,ibn),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_urt/pi:', */
/*     -         (flx_urt(k,l)/pi,k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'rad_rt:    ', */
/*     -         (rad_rt(nze,1,k,l),k=1,nlev) */
/* 1     continue */
/*      endif */

    i__1 = *nphi;
    for (naz = 1; naz <= i__1; ++naz) {
	i__2 = *numu;
	for (nze = *nzup; nze <= i__2; ++nze) {
	    rad_rb__[nze + (naz + (nlev + l * 70 << 4) << 4)] = brdf_b__[nze 
		    + (naz + (l + ibn * 5 << 4) << 4)];
/* L2001: */
	}
/* L2021: */
    }

    i__1 = *nlyr;
    for (k = 1; k <= i__1; ++k) {
	dflx = flx_drb__[nlev - k + 1 + l * 70] / pi;
	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__3 = *numu;
	    for (nze = *nzup; nze <= i__3; ++nze) {
		rad_rb__[nze + (naz + (nlev - k + l * 70 << 4) << 4)] = 
			refrad_b__[nze + (nlev - k + (l + ibn * 5) * 70 << 4)]
			 + trnrad_b__[nze + (nlev - k + (l + ibn * 5) * 70 << 
			4)] * rad_rb__[nze + (naz + (nlev - k + 1 + l * 70 << 
			4) << 4)] * dflx;
/* L2201: */
	    }
/* L2221: */
	}
/* L2241: */
    }

/* ****  rad_rb(nze,nlev,l) is defined in rad_src, */
/*      pi*uu(numu-nze+1,k,naz)/rfldn(k) */

/*      if(ng0 .eq. 1) then */
/*      do 2 nze=nzup,numu */
/*      write(*,'(/,1a,3i5)') 'in rad_ref: nze, l:',nze,l */
/*      write(*,'(1a,16(1pe12.4))') 'trnrad_b:  ', */
/*     -         (trnrad_b(nze,k,l,ibn),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'refrad_b:  ', */
/*     -         (refrad_b(nze,k,l,ibn),k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'flx_drb/pi:', */
/*     -         (flx_drb(k,l)/pi,k=1,nlev) */
/*      write(*,'(1a,16(1pe12.4))') 'rad_rb:    ', */
/*     -         (rad_rb(nze,1,k,l),k=1,nlev) */
/* 2     continue */
/*      endif */

    return 0;
} /* rad_ref__ */

