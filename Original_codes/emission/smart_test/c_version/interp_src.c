/* interp_src.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int interp_src__(logical *lsolar, logical *lplanck, integer *
	ng0, integer *nlev, integer *l_1__, integer *l_2__, integer *lev1, 
	integer *l2tau, integer *nzup, integer *nzdn, integer *numu, integer *
	nphi, doublereal *dnsflxsrc_b__, doublereal *upsflxsrc_b__, 
	doublereal *sradsrc_b__, doublereal *ddnsflxdx, doublereal *dupsflxdx,
	 doublereal *dsraddx, doublereal *dntflxsrc_b__, doublereal *
	uptflxsrc_b__, doublereal *tradsrc_b__, doublereal *ddntflxdx, 
	doublereal *duptflxdx, doublereal *dtraddx, real *dalb, real *deltau, 
	real *delpi0, real *delg, real *alb0, doublereal *wn_io__, doublereal 
	*up_s_src__, doublereal *dn_s_src__, doublereal *rad_s_src__, 
	doublereal *up_t_src__, doublereal *dn_t_src__, doublereal *
	rad_t_src__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer k, ibn, naz, nze;


/* cccccccccccccccccccccc   i n t e r p _ s r c   cccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine uses linear interpolation to scale the solar    cc */
/* c    and thermal radiative sources to this spetral point.            cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c        ng0 - spectral mapping bin index for this spectral point    cc */
/* c       nlev - number of model vertical levels                       cc */
/* c         l1 - starting vertical index for source calculation        cc */
/* c         l2 - ending vertical index for source calculation          cc */
/* c       lev1 - first vertical level for source interpolation         cc */
/* c      l2tau - last vertical level for source function interpolation cc */
/* c       nzup - number of upward steams for radiances                 cc */
/* c       nzdn - number of downward streams for radiances              cc */
/* c       numu - total number of streams for radiances                 cc */
/* c       nphi - total number of azimuths for radiances                cc */
/* c dnsflxsrc_b - downward solar flux source at each level for bin     cc */
/* c upsflxsrc_b - downward solar flux source at each level for bin     cc */
/* c   sradsrc_b - angle-dependent solar radiance source at each level  cc */
/* c              for this bin                                          cc */
/* c  ddnsflxdx - downward solar flux jacobian at each level for  bin   cc */
/* c  dupsflxdx - upward solar flux jacobian at each level for bin      cc */
/* c    dsraddx - solar radiance jacobian at each level for bin         cc */
/* c dntflxsrc_b - downward thermal flux source at each level for bin   cc */
/* c uptflxsrc_b - downward thermal flux source at each level for bin   cc */
/* c   tradsrc_b - angle-dependent thermal radiance source at each      cc */
/* c              level for this bin                                    cc */
/* c  ddntflxdx - downward thermal flux jacobian at each level for bin  cc */
/* c  duptflxdx - upward thermal flux jacobian at each level for bin    cc */
/* c    dtraddx - thermal radiance jacobian at each level for bin       cc */
/* c       wnio - current wavenumber (cm**-1)                           cc */
/* c       dalb - difference between albedo at wnio and binned albedo   cc */
/* c     deltau - difference between optical depth in each layer and    cc */
/* c              binned value                                          cc */
/* c     delpi0 - difference between single scattering albedo in each   cc */
/* c              layer and binned value                                cc */
/* c       delg - difference between scattering asymmetry parameter in  cc */
/* c              each layer and binned value                           cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c   up_s_src - interpolated/scaled upward solar flux source for      cc */
/* c              each level at this spectral point                     cc */
/* c   dn_s_src - interpolated/scaled downward solar flux source for    cc */
/* c              each level at this spectral point                     cc */
/* c  rad_s_src - interpolated/scaled solar radinace source for each    cc */
/* c              level at this spectral point                          cc */
/* c   up_s_src - interpolated/scaled upward thermal flux source for    cc */
/* c              each level at this spectral point                     cc */
/* c   dn_s_src - interpolated/scaled downward thermal flux source for  cc */
/* c              each level at this spectral point                     cc */
/* c  rad_s_src - interpolated/scaled thermal radinace source for each  cc */
/* c              level at this spectral point                          cc */
/* c                                                                    cc */
/* cccccccccccccccccccccc   i n t e r p _ s r c   cccccccccccccccccccccccc */




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



/*      integer kk */

/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */



/* ****   interpolated solar and thermal source functions */




/* ****   specify the bin index */

    /* Parameter adjustments */
    rad_t_src__ -= 273;
    --dn_t_src__;
    --up_t_src__;
    rad_s_src__ -= 273;
    --dn_s_src__;
    --up_s_src__;
    --delg;
    --delpi0;
    --deltau;
    dtraddx -= 89873;
    duptflxdx -= 351;
    ddntflxdx -= 351;
    tradsrc_b__ -= 18193;
    uptflxsrc_b__ -= 71;
    dntflxsrc_b__ -= 71;
    dsraddx -= 89873;
    dupsflxdx -= 351;
    ddnsflxdx -= 351;
    sradsrc_b__ -= 18193;
    upsflxsrc_b__ -= 71;
    dnsflxsrc_b__ -= 71;

    /* Function Body */
    ibn = *ng0;

    if (*lsolar) {

/* ********************************************************************** */

/* *****         f i n d   s o l a r    s o u r c e s */

/* ********************************************************************** */

	if (*l_1__ == 1) {

/* ****        set downward flux and radiance source functions */
/*            at the top boundary */

	    k = 1;
	    dn_s_src__[k] = 0.f;
	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {
		i__2 = *nzdn;
		for (nze = 1; nze <= i__2; ++nze) {
		    rad_s_src__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L1001: */
		}
/* L1021: */
	    }

	}

	if (*l_2__ == *nlev) {

/* ****      set upward radiances at surface */

	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {
		i__2 = *numu;
		for (nze = *nzup; nze <= i__2; ++nze) {
		    rad_s_src__[nze + (naz + (*nlev << 4) << 4)] = 
			    sradsrc_b__[nze + (naz + (*nlev + ibn * 70 << 4) 
			    << 4)] + dsraddx[nze + (naz + (*nlev + ((ibn << 2)
			     + 4) * 70 << 4) << 4)] * *dalb;
/* L1041: */
		}
/* L1061: */
	    }

	}

/* ****    find the downward radiance sources for all other levels */

	i__1 = *l2tau + 1;
	for (k = *lev1; k <= i__1; ++k) {

/* ****              interpolate wavenumber-dependent downward fluxes */

	    dn_s_src__[k] = dnsflxsrc_b__[k + ibn * 70] + ddnsflxdx[k + ((ibn 
		    << 2) + 1) * 70] * deltau[k - 1] + ddnsflxdx[k + ((ibn << 
		    2) + 2) * 70] * delpi0[k - 1] + ddnsflxdx[k + ((ibn << 2) 
		    + 3) * 70] * delg[k - 1] + ddnsflxdx[k + ((ibn << 2) + 4) 
		    * 70] * *dalb;

/* ****              interpolate wavenumber-dependent downward radiances */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
		    rad_s_src__[nze + (naz + (k << 4) << 4)] = sradsrc_b__[
			    nze + (naz + (k + ibn * 70 << 4) << 4)] + dsraddx[
			    nze + (naz + (k + ((ibn << 2) + 1) * 70 << 4) << 
			    4)] * deltau[k - 1] + dsraddx[nze + (naz + (k + ((
			    ibn << 2) + 2) * 70 << 4) << 4)] * delpi0[k - 1] 
			    + dsraddx[nze + (naz + (k + ((ibn << 2) + 3) * 70 
			    << 4) << 4)] * delg[k - 1] + dsraddx[nze + (naz + 
			    (k + ((ibn << 2) + 4) * 70 << 4) << 4)] * *dalb;
/* L1201: */
		}
/* L1221: */
	    }

/* ****        interpolate upward radiances and fluxes */

	    up_s_src__[k - 1] = upsflxsrc_b__[k - 1 + ibn * 70] + dupsflxdx[k 
		    - 1 + ((ibn << 2) + 1) * 70] * deltau[k - 1] + dupsflxdx[
		    k - 1 + ((ibn << 2) + 2) * 70] * delpi0[k - 1] + 
		    dupsflxdx[k - 1 + ((ibn << 2) + 3) * 70] * delg[k - 1] + 
		    dupsflxdx[k - 1 + ((ibn << 2) + 4) * 70] * *dalb;

/* ****           interpolate wavenumber-dependent upward radiances */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
		    rad_s_src__[nze + (naz + (k - 1 << 4) << 4)] = 
			    sradsrc_b__[nze + (naz + (k - 1 + ibn * 70 << 4) 
			    << 4)] + dsraddx[nze + (naz + (k - 1 + ((ibn << 2)
			     + 1) * 70 << 4) << 4)] * deltau[k - 1] + dsraddx[
			    nze + (naz + (k - 1 + ((ibn << 2) + 2) * 70 << 4) 
			    << 4)] * delpi0[k - 1] + dsraddx[nze + (naz + (k 
			    - 1 + ((ibn << 2) + 3) * 70 << 4) << 4)] * delg[k 
			    - 1] + dsraddx[nze + (naz + (k - 1 + ((ibn << 2) 
			    + 4) * 70 << 4) << 4)] * *dalb;
/* L1241: */
		}
/* L1261: */
	    }

/* L1281: */
	}

    }

    if (*lplanck) {

/* ********************************************************************** */

/* *****         f i n d   t h e r m a l    s o u r c e s */

/* ********************************************************************** */

	if (*l_1__ == 1) {

/* ****        set downward flux and radiance source functions */
/*            at the top boundary */

	    k = 1;
	    dn_t_src__[k] = 0.f;
	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {
		i__2 = *nzdn;
		for (nze = 1; nze <= i__2; ++nze) {
		    rad_t_src__[nze + (naz + (k << 4) << 4)] = 0.f;
/* L2201: */
		}
/* L2221: */
	    }

	}

	if (*l_2__ == *nlev) {

/* ****      set upward fluxes and radiances at surface */

	    k = *nlev;

	    up_t_src__[k] = uptflxsrc_b__[k + ibn * 70] + duptflxdx[k + ((ibn 
		    << 2) + 4) * 70] * *dalb;

	    i__1 = *nphi;
	    for (naz = 1; naz <= i__1; ++naz) {
		i__2 = *numu;
		for (nze = *nzup; nze <= i__2; ++nze) {
		    rad_t_src__[nze + (naz + (k << 4) << 4)] = tradsrc_b__[
			    nze + (naz + (k + ibn * 70 << 4) << 4)] + dtraddx[
			    nze + (naz + (k + ((ibn << 2) + 4) * 70 << 4) << 
			    4)] * *dalb;
/* L2441: */
		}
/* L2461: */
	    }

	}

/* ****    find the downward radiance sources for all other levels */

	i__1 = *l2tau + 1;
	for (k = *lev1; k <= i__1; ++k) {

/* ****              interpolate wavenumber-dependent downward fluxes */

	    dn_t_src__[k] = dntflxsrc_b__[k + ibn * 70] + ddntflxdx[k + ((ibn 
		    << 2) + 1) * 70] * deltau[k - 1] + ddntflxdx[k + ((ibn << 
		    2) + 2) * 70] * delpi0[k - 1] + ddntflxdx[k + ((ibn << 2) 
		    + 3) * 70] * delg[k - 1] + ddntflxdx[k + ((ibn << 2) + 4) 
		    * 70] * *dalb;

/* ****              interpolate wavenumber-dependent downward radiances */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *nzdn;
		for (nze = 1; nze <= i__3; ++nze) {
		    rad_t_src__[nze + (naz + (k << 4) << 4)] = tradsrc_b__[
			    nze + (naz + (k + ibn * 70 << 4) << 4)] + dtraddx[
			    nze + (naz + (k + ((ibn << 2) + 1) * 70 << 4) << 
			    4)] * deltau[k - 1] + dtraddx[nze + (naz + (k + ((
			    ibn << 2) + 2) * 70 << 4) << 4)] * delpi0[k - 1] 
			    + dtraddx[nze + (naz + (k + ((ibn << 2) + 3) * 70 
			    << 4) << 4)] * delg[k - 1] + dtraddx[nze + (naz + 
			    (k + ((ibn << 2) + 4) * 70 << 4) << 4)] * *dalb;
/* L2601: */
		}
/* L2621: */
	    }

/* ****        interpolate upward radiances and flux sourceses */

	    up_t_src__[k - 1] = uptflxsrc_b__[k - 1 + ibn * 70] + duptflxdx[k 
		    - 1 + ((ibn << 2) + 1) * 70] * deltau[k - 1] + duptflxdx[
		    k - 1 + ((ibn << 2) + 2) * 70] * delpi0[k - 1] + 
		    duptflxdx[k - 1 + ((ibn << 2) + 3) * 70] * delg[k - 1] + 
		    duptflxdx[k - 1 + ((ibn << 2) + 4) * 70] * *dalb;

/* ****           interpolate wavenumber-dependent upward radiances */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__3 = *numu;
		for (nze = *nzup; nze <= i__3; ++nze) {
		    rad_t_src__[nze + (naz + (k - 1 << 4) << 4)] = 
			    tradsrc_b__[nze + (naz + (k - 1 + ibn * 70 << 4) 
			    << 4)] + dtraddx[nze + (naz + (k - 1 + ((ibn << 2)
			     + 1) * 70 << 4) << 4)] * deltau[k - 1] + dtraddx[
			    nze + (naz + (k - 1 + ((ibn << 2) + 2) * 70 << 4) 
			    << 4)] * delpi0[k - 1] + dtraddx[nze + (naz + (k 
			    - 1 + ((ibn << 2) + 3) * 70 << 4) << 4)] * delg[k 
			    - 1] + dtraddx[nze + (naz + (k - 1 + ((ibn << 2) 
			    + 4) * 70 << 4) << 4)] * *dalb;
/* L2641: */
		}
/* L2661: */
	    }

/* L2681: */
	}
/*       write(*,'(/,1a,2i5,1pe14.6)') 'interp_src: ibn,l_1, wn_io', */
/*     -  ibn,l_1,wn_io */
/*       write(*,'(1a,16(1pe12.4))') 'dn_t_src  ', */
/*     - (dn_t_src(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'dntflxsrc_b: ', */
/*     - (dntflxsrc_b(k,ibn),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'ddntf*dtau', */
/*     - ((ddntflxdx(k-1,1,ibn)*deltau(k-1)),k=2,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'ddntf*dpi0', */
/*     - ((ddntflxdx(k-1,2,ibn)*delpi0(k-1)),k=2,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'ddntf*dg  ', */
/*     - ((ddntflxdx(k-1,3,ibn)*delg(k-1)),k=2,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'ddntf*dalb', */
/*     - ((ddntflxdx(k-1,4,ibn)*dalb),k=2,nlev) */

/*       write(*,'(/,1a,16(1pe12.4))') 'up_t_src: ', */
/*     - (up_t_src(k),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'uptflxsrc_b: ', */
/*     - (uptflxsrc_b(k,ibn),k=1,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'duptflx*dtau', */
/*     - ((duptflxdx(k-1,1,ibn)*deltau(k-1)),k=2,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'duptflx*dpi0', */
/*     - ((duptflxdx(k-1,2,ibn)*delpi0(k-1)),k=2,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'duptflx*dg  ', */
/*     - ((duptflxdx(k-1,3,ibn)*delg(k-1)),k=2,nlev) */
/*       write(*,'(1a,16(1pe12.4))') 'duptflx*dalb', */
/*     - ((duptflxdx(k-1,4,ibn)*dalb),k=2,nlev) */

    }

    return 0;
} /* interp_src__ */

