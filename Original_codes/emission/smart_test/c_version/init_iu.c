/* init_iu.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int init_iu__(integer *iuheat, integer *iustat, integer *
	iuatm, integer *iumie, integer *iumix, integer *iusur, integer *iugrp,
	 integer *iuthrm, integer *iutrn, integer *iuaer, integer *iuflx, 
	integer *iuout, integer *iugas, integer *iusol0, integer *iusol1, 
	integer *iusol2, integer *iu_pd__)
{

/* ccccccccccccccccccccccccc   i n i t _ i u    cccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This subroutine initializes the i/o unit numbers for SMART      cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c    nsol - maximum number of solar zenith angles                    cc */
/* c    nmode - maximum number of particle modes                        cc */
/* c    mxlout - maximum number of output radiance/flux levels          cc */
/* c    ngas - maximum number of absorbing gases                        cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    iuheat - output unit for heating and cooling rates.             cc */
/* c    iustat - output unit for binning statistics                     cc */
/* c    iuatm - input unit for input atmospheric thermal structure      cc */
/* c    iumie - first input unit for aerosol optical properties         cc */
/* c    iumix - first input unit for gas mixing ratios                  cc */
/* c    iusur - input unit for surface optical properties               cc */
/* c    iugrp - scratch i/o unit for binned optcal properties           cc */
/* c    iuthrm - scratch i/o unit for binned thermal radiances          cc */
/* c    iutrn - output unit with transmission/pressure file             cc */
/* c    iuaer - i/o unit with combined aerosol optical properties       cc */
/* c    iuflx - output unit with flx vs pressure                        cc */
/* c    iuout - output unit for wavenumber-dependent radiances/fluxes   cc */
/* c    iugas - first input unit with gas gas optical properties        cc */
/* c    iusol0 - solar flux scratch file used by map_back               cc */
/* c    iusol1 - first of two solar flux input units                    cc */
/* c    iusol2 - second of two solar flux input units                   cc */
/* c    iu_pd - output unit for flux and radiance jacobians             cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc   i n i t _ i u    cccccccccccccccccccccccccc */



/* ****   define the output unit numbers for heating rates */
/*        (iuheat), and binning statistics (iustat) */


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



    *iuheat = 9;
    *iustat = 10;

/* ****  define the unit numbers for the input atmospheric thermal */
/*      structure, aerosol optical properties, gas mixing ratios, */
/*      and surface albedos. */

    *iuatm = 11;
    *iumie = 12;
    *iumix = 13;
    *iusur = 14;

/* ****  define a scratch unit for bin numbers at each wavenumber */

    *iugrp = 15;

/* ****  define a scratch unit for thermal fluxes, iuthrm */

    *iuthrm = 18;

/* ****   define an output unit for optical depths, iutrn */

    *iutrn = 19;

/* ****   define input unit for aerosol optical properties, iuaer */
/*       note: nmode units are needed. */

    *iuaer = 20;

/* ****   define the first output unit for the level-dependent fluxes. */
/*       note: there are up to nsol of these units */

    *iuflx = *iuaer + 10;

/* ****   define units for output flux/radiance files */
/*       note: up to nsol units are needed */

    *iuout = *iuflx + 10;

/* ****   define input units for gas absorption coefficeints */
/*       note: ngas units are needed (this can be a large number). */

    *iugas = *iuout + 30;

/* ****  note: for solar flux files ordered in increasing wavenumber, */
/*      only 2 unit numbers are needed.  for large files ordered in */
/*      increasing wavelength, a large number of files are needed */
/*      to reverse the order of the solar flux files.  The unit iusol1 */
/*      is released after the subroutine 'readsol' completes, and is */
/*      then reused in 'backmap' */

    *iusol0 = *iugas + 51;
    *iusol1 = *iusol0 + 1;
    *iusol2 = *iusol1 + 1;

/* ****    define unit number for partial derivitives */

    *iu_pd__ = *iusol2 + 1;

    return 0;
} /* init_iu__ */

