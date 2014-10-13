/* smart_spectra.f -- translated by f2c (version 20100827).
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
static integer c__132 = 132;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double acos(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int smart_in__(integer *, char *, integer *, char 
	    *, integer *, integer *, integer *, real *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    char *, integer *, char *, integer *, char *, integer *, integer *
	    , integer *, integer *, integer *, real *, real *, integer *, 
	    integer *, char *, integer *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, real *, real *, logical *, integer *, real *, integer *, char *
	    , integer *, char *, integer *, real *, real *, integer *, 
	    integer *, integer *, real *, real *, integer *, char *, integer *
	    , char *, integer *, integer *, integer *, integer *, real *, 
	    real *, integer *, integer *, integer *, real *, real *, real *, 
	    real *, real *, real *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, logical *, 
	    integer *, real *, integer *, real *, integer *, real *, real *, 
	    integer *, integer *, real *, real *, logical *, logical *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static real dtau_tot__;
    static doublereal dirsoflx[700]	/* was [70][10] */;
    static integer k, m, n;
    static real p[70], t[70], z__[70];
    extern /* Subroutine */ int smart_hdr__(integer *, integer *, integer *, 
	    char *, char *, char *, char *, char *, char *, char *, integer *,
	     integer *, integer *, integer *, real *, real *, real *, real *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, integer *, integer *
	    , real *, integer *, real *, integer *, real *, real *, integer *,
	     integer *, real *, real *, logical *, integer *, integer *, 
	    integer *, real *, real *, integer *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ii, ne, ng;
    static real au, pi, ts, ws;
    static integer icf[10];
    static real ang;
    static integer icp, nlb, len, ict;
    static real phi[16], alt[70];
    static integer igs[10], ntb;
    static real scp;
    static integer ngt, nza;
    static doublereal dwn;
    static integer nze;
    static real umu[16];
    static integer icg1[10], icg2[10];
    static real phi0[10];
    static integer ixs0;
    static real umu0[10];
    static integer iang[10], irad;
    static char name__[1*132];
    static integer igas[50], mode, iref;
    static real pbar;
    static logical lcsh[10];
    static integer nref, nphi;
    static real grav[70], ratm;
    static integer nlev, imix[50];
    static real phiw;
    static integer next;
    static real rmix[3500]	/* was [70][50] */;
    static integer nlyr, numu, nstr;
    static real pout[3], umu_1__[16];
    static char frms0[40];
    static integer icalb[4], igpab[150]	/* was [3][50] */, iuabc[150]	/* 
	    was [3][50] */;
    static real scalb[4], dp_dp__[3];
    static integer iccsh[10];
    static real accur;
    static integer iu_pd__;
    static real sccsh[10];
    static integer iuaer, ictau[10], icomp[6], iumie, iugas, icsol, ncomp, 
	    iopen, icmix[50], iuatm, icwlq[10], k_out__[3];
    static doublereal width;
    static integer igtrn[150]	/* was [3][50] */, islit, iugrp, iuflx, iutpd[
	    100]	/* was [10][10] */, iuspd[100]	/* was [10][10] */;
    static doublereal wnmin;
    static integer iumix, iztau[10];
    static doublereal wnmax, wneof[2];
    static integer izmix[50], iutrn, nlout;
    static doublereal wnaer[10];
    static real sctau[10];
    static integer iuout;
    static real scsol;
    static integer iusur;
    static real scmix[50], umu_o__[16], gwt_o__[16], sgrav;
    extern /* Subroutine */ int smart_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, char *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, real *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, logical *, 
	    integer *, integer *, real *, integer *, real *, integer *, real *
	    , real *, integer *, integer *, real *, real *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, real *, real *, char *, real *, real *, 
	    real *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, ftnlen, ftnlen);
    static integer ioffs0, ifrms0;
    static real pi0err;
    static integer icwns0, iuins0, iusol0, iusol1, iusol2;
    static real scwns0;
    static doublereal aid_lr__[70];
    static integer io_end__[113];
    static logical lamber;
    static integer icqsca[10];
    static real scwalb;
    static doublereal wn_eof__[226]	/* was [2][113] */, thheat[70];
    static integer iuheat, ngases, nalbmx, io_err__[113];
    static doublereal soheat[700]	/* was [70][10] */;
    extern /* Subroutine */ int charsp_(char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer nmodes;
    static char frmatm[40];
    static integer icptau[10], iskatm, icpmix[50], istate[113];
    static real phferr;
    static logical lsolar;
    static integer new_pt__;
    static logical usrang;
    static integer icqext[10], nstate;
    static char frmmix[40*50];
    static doublereal wn_tol__;
    static real scptau[10], scpmix[50], wgtgas[50], wgtatm;
    static integer iuthrm, iustat;
    static real tauerr;
    static integer iunits;
    static real radius;
    extern /* Subroutine */ int cldstr_(integer *, integer *, integer *, 
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    real *, real *, logical *, integer *, real *, real *, real *, 
	    real *, real *, ftnlen);
    static char frmsur[40];
    static integer levout;
    static real volmix[6];
    extern /* Subroutine */ int qgausn_(integer *, real *, real *);
    static integer iwnsur;
    static real pd_frac__[113];
    static char io_file__[132*113], aerfile[132*10];
    extern /* Subroutine */ int readmie_(integer *, integer *, integer *, 
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, ftnlen);
    static char gasfile[132*3*50], miefile[132*10];
    static integer ioffmie[10];
    extern /* Subroutine */ int mod_atm__(integer *, integer *, integer *, 
	    integer *, integer *, char *, integer *, char *, integer *, 
	    integer *, integer *, real *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, integer *, integer *
	    , char *, integer *, char *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static char atmfile[132];
    static logical lplanck;
    static real dtauaer[700]	/* was [10][70] */;
    static integer iupdrad[76800]	/* was [16][16][10][3][10] */;
    extern /* Subroutine */ int readsol_(char *, integer *, integer *, 
	    integer *, integer *, char *, integer *, integer *, integer *, 
	    integer *, real *, real *, integer *, doublereal *, doublereal *, 
	    real *, integer *, integer *, doublereal *, ftnlen, ftnlen);
    static integer ntau_pd__, ioffmom[10];
    static char mixfile[132*50], solfile[132];
    static integer ifrmatm, iofftau[10], ioffmix[50];
    extern /* Subroutine */ int init_iu__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), readsur_(char *, integer *, integer *, 
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, ftnlen, ftnlen);
    static integer new_mix__;
    static doublereal dnthflx[70];
    static char surfile[132];
    static integer isource, ifrmmix[50], ioffsur;
    static doublereal dnsoflx[700]	/* was [70][10] */;
    static real tau_tot__;
    static integer ngastyp[50], ifrmout, ifrmsur;
    static real surferr;
    static integer icwnsur;
    static doublereal upthflx[70];
    static integer isptype;
    static doublereal upsoflx[700]	/* was [70][10] */;

    /* Fortran I/O blocks */
    static cilist io___161 = { 0, 6, 0, "(/,1a,/,/,1a)", 0 };
    static cilist io___165 = { 0, 6, 0, "(5(1pe12.4))", 0 };
    static cilist io___180 = { 0, 6, 0, "(/,/,1a,/,/,1a)", 0 };
    static cilist io___184 = { 0, 6, 0, "(1x,i3,2(1pe13.5))", 0 };
    static cilist io___185 = { 0, 6, 0, "(/,/,1a,/,/,2a)", 0 };
    static cilist io___190 = { 0, 6, 0, "(3i8,2x,2(1pe14.6),5x,132a)", 0 };
    static cilist io___192 = { 0, 6, 0, "(3i8,2x,2(1pe14.6),5x,132a)", 0 };
    static cilist io___193 = { 0, 6, 0, "(3i8,2x,2(1pe14.6),5x,132a)", 0 };



/* ccccccccccccccccccccc    s m a r t _ s p e c t r a  ccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This program combines spectral mapping methods and the          cc */
/* c    discrete ordinate model (Stamnes et al. 1988) to generate       cc */
/* c    high-resolution synthetic monochromatic radiances in            cc */
/* c    vertically-inhomogeneous scattering absorbing, emitting,        cc */
/* c    planetary atmospheres.                                          cc */
/* c                                                                    cc */
/* c    note:  This version of the program uses albedo gradients.       cc */
/* c           It also outputs wavelength-dependent, level-dependent    cc */
/* c           contribution functions.                                  cc */
/* c                                                                    cc */
/* c   9/28/97:   this version prints the direct solar flux as well as  cc */
/* c              the total downward and upward fluxes, and allows      cc */
/* c              wavelength-dependent radiances to be printed          cc */
/* c              at mxlout levels.                                     cc */
/* c                                                                    cc */
/* c              as part of this modification, the slit function       cc */
/* c              routine, slit, was modified such that it only allows  cc */
/* c              boxcar or triangular slits (other functions required  cc */
/* c              far too much RAM.                                     cc */
/* c                                                                    cc */
/* c    5/3/98:   the aerosol, gas, and rayliegh scattering optical     cc */
/* c              depth sections have been extracted from the main      cc */
/* c              program and are called as subroutines.                cc */
/* c                                                                    cc */
/* c    6/3/98:   smart.f was modified to include a 0th order           cc */
/* c              correction for the solar zenith angle to account for  cc */
/* c              the sphericity of the planet.                         cc */
/* c                                                                    cc */
/* c    5/28/01   smart modified to evaluate radiances for the          cc */
/* c              minimum and maximum optical depth within each         cc */
/* c              bin, and interpolate linearly in tau                  cc */
/* c                                                                    cc */
/* c    5/29/03:  smart modified to use DISORT 2.0                      cc */
/* c                                                                    cc */
/* c    12/26/03: smart modified to provide an option to elminate       cc */
/* c              headers from binary out files (ifrmout>2)             cc */
/* c    11/14/04: bdrf mode (iref) propagated through program           cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      iuabc - unit numbers for gas line parameters                  cc */
/* c     iusol0 - unit number as a scratch file for solar fluxes        cc */
/* c     iusol1 - unit number for solar fluxes                          cc */
/* c     iusol2 - unit number for solar flux scratch file               cc */
/* c     iuthrm - unit number for thermal radiances                     cc */
/* c      iuout - unit number of output radiance file                   cc */
/* c      iutrn - unit number of output transmission/pressure file      cc */
/* c      iutpd - unit numbers for output level-dependent               cc */
/* c              thermal fluxes and thier partial derivatives          cc */
/* c      iuspd - unit numbers for output level-dependent               cc */
/* c              solar fluxes and thier partial derivatives            cc */
/* c    iupdrad - unit numbers for output level-dependent radiances     cc */
/* c              and thier partial derivatives                         cc */
/* c    atmfile - name of input atmospheric structure file              cc */
/* c    aerfile - name of input aerosol vertical structure file         cc */
/* c    miefile - name of file with aerosol optical properties vs. wn   cc */
/* c    solfile - name of file with wn-dependent solar fluxes           cc */
/* c    surfile - name of file with wn-dependent surface albedos        cc */
/* c    mixfile - name of file with gas mixing ratios                   cc */
/* c    gasfile - name fo file with gas absorption coeffiecients vs. wn cc */
/* c     nmodes - number of discrete aerosol partical modes             cc */
/* c      ncomp - number of rayleigh-scattering constituentes           cc */
/* c      icomp - index of each rayleigh scattering constituent         cc */
/* c     volmix - volume mixing ratio of each rayleigh scatterer        cc */
/* c         ts - surface temperature (K)                               cc */
/* c         au - distance to the sun (in AU's)                         cc */
/* c     tauerr - optical depth relative binning error (0. to ~0.8)     cc */
/* c     pi0err - co-single scattering albedo absolute binning error    cc */
/* c     phferr - asymmetry factor absolute binning error               cc */
/* c    surferr - surface albedo absolute binning error                 cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angles (degrees)                        cc */
/* c       pout - output pressure level (bars)                          cc */
/* c      accur - azimuth convergence accuracy for D/O routine          cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c              note: if lamber = F, use a BRDF is used.              cc */
/* c    isource - index of source function type: (1) solar only         cc */
/* c              (2) thermal only, (3) both                            cc */
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
/* c    ifmrout - index of output file format (1) ascii, (2) binary,    cc */
/* c                                          (3) binary, no header     cc */
/* c       iref - bidirectional reflectance options                     cc */
/* c              0 - Lambert                                           cc */
/* c              1 - Hapke's BDR model                                 cc */
/* c              2 - Breon's BDR model; combination of Li + Roujean    cc */
/* c              3 - Roujean's BDR model                               cc */
/* c              4 - Cox and Munk glint model                          cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) Watts/m**2/sr/Angstrom                             cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       nlev - number of output levels                               cc */
/* c       nlyr - number of computational model layers                  cc */
/* c     levout - output level index (1) top of atmosphere,             cc */
/* c              (2) surface, (3) arbitrary level                      cc */
/* c        nza - number of solar zenith angles                         cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c       nref - number of surface optical properties specified at     cc */
/* c              each wavelength.                                      cc */
/* c  SURF_PR : Wavelength dependent surface properties array           cc */
/* c            IREF= 0 - Lambert albedo                                cc */
/* c            IREF= 1 - Hapke : HH, W                                 cc */
/* c            IREF= 2 - Breon's BDR model: k0, k1, k2                 cc */
/* c            IREF= 3 - Roujean's BDR model: k0, k1, k2               cc */
/* c            IREF= 4 - Cox and Munk glint model: n, k, ws, phiw      cc */
/* c     usrang - output radiances at user angles? (logical: T/F)       cc */
/* c        phi - emission azimuth angles (degrees)                     cc */
/* c       umu0 - solar zenith angle cosines                            cc */
/* c       phi0 - solar azimuth angle cosines                           cc */
/* c      wnmin - minimum wavenumber of desired spectral window         cc */
/* c      wnmax - maximum wavenumber of desired spectral window         cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c         wl - current wavenumber (cm**-1)                           cc */
/* c         wn - current wavenumber (cm**-1)                           cc */
/* c    dir_flx - direct downward solar flux at wavenumber, wn          cc */
/* c     dn_flx - total downward flux (irradiance) at wavenumber wn     cc */
/* c     up_flx - upward flux (irradiance) at wavenumber wn             cc */
/* c        rad - radiance at each output zenith and azimuth angle      cc */
/* c  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc */
/* c  trn_gas_0 - normal-incidence gas transmission                     cc */
/* c  trn_ray_0 - normal-incidence aerosol transmission                 cc */
/* c                                                                    cc */
/* ccccccccccccccccccccc    s m a r t _ s p e c t r a  ccccccccccccccccccc */




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



/* ****   define variables used by discrete ordinate method */








/* ****    units for flux and radiance jacobians */



/* ****   atmospheric structure variables */


/* ****   double precision variables for heating rates */


/* ****   define the minimum spectral interval width: */

    wn_tol__ = 5e-7;

/* ****   define the input and output unit numbers for all i/o */

    init_iu__(&iuheat, &iustat, &iuatm, &iumie, &iumix, &iusur, &iugrp, &
	    iuthrm, &iutrn, &iuaer, &iuflx, &iuout, &iugas, &iusol0, &iusol1, 
	    &iusol2, &iu_pd__);

/* ****   read in all user-supplied input */

    smart_in__(&iuatm, atmfile, &ifrmatm, frmatm, &iskatm, &icp, &ict, &scp, &
	    ngases, ngastyp, igas, &iugas, iuabc, igtrn, igpab, wgtgas, 
	    gasfile, &iumix, mixfile, ifrmmix, frmmix, ioffmix, izmix, imix, 
	    icpmix, icmix, scpmix, scmix, &nmodes, &iuaer, aerfile, &iumie, 
	    miefile, ioffmie, ioffmom, iang, icwlq, icqext, icqsca, icg1, 
	    icg2, icf, wnaer, iofftau, iztau, icptau, ictau, scptau, sctau, 
	    lcsh, iccsh, sccsh, &iusur, surfile, &ifrmsur, frmsur, &ioffsur, &
	    ws, &phiw, &iwnsur, &icwnsur, icalb, &scwalb, scalb, &iusol1, 
	    solfile, &ifrms0, frms0, &ioffs0, &ixs0, &icwns0, &icsol, &scwns0,
	     &scsol, &iuins0, &ncomp, icomp, volmix, &ts, &au, &radius, &
	    sgrav, &wgtatm, &wnmin, &wnmax, &isptype, &islit, &width, &dwn, &
	    nstr, &usrang, &numu, umu, &nphi, phi, &nza, umu0, phi0, &levout, 
	    &nlout, pout, &accur, &lamber, &lplanck, &lsolar, &isource, &irad,
	     &iref, &nref, &iuout, &iunits, &ifrmout, &iuheat, &iustat, &
	    iutrn, &iuflx, &tauerr, &pi0err, &phferr, &surferr, pd_frac__, &
	    nstate, istate, &iu_pd__, iutpd, iuspd, iupdrad, (ftnlen)132, (
	    ftnlen)40, (ftnlen)132, (ftnlen)132, (ftnlen)40, (ftnlen)132, (
	    ftnlen)132, (ftnlen)132, (ftnlen)40, (ftnlen)132, (ftnlen)40);

/* ****  read the atmospheric structure file and the gas mixing ratio */
/*      files and set up the model atmosphere: */

    nlev = 0;
    new_pt__ = 1;
    new_mix__ = 1;
    for (k = 1; k <= 70; ++k) {
	p[k - 1] = 0.f;
	t[k - 1] = 0.f;
	alt[k - 1] = 0.f;
	z__[k - 1] = 0.f;
	grav[k - 1] = 0.f;
/* L1001: */
    }

    mod_atm__(&iuatm, &iumix, &new_pt__, &new_mix__, istate, atmfile, &
	    ifrmatm, frmatm, &iskatm, &icp, &ict, &scp, &nlev, &nlyr, &levout,
	     &nlout, k_out__, pout, &wgtatm, wgtgas, &ngases, igas, mixfile, 
	    ifrmmix, frmmix, ioffmix, izmix, imix, icpmix, icmix, scpmix, 
	    scmix, pd_frac__, &radius, &sgrav, &ratm, p, t, alt, z__, grav, 
	    rmix, dp_dp__, (ftnlen)132, (ftnlen)40, (ftnlen)132, (ftnlen)40);

/* ****   define optical depth counters for */
/*       rayleigh scattering (ne = 1) gas absorption */

    ne = 1;
    wn_eof__[(ne << 1) - 2] = wnmin;
    wn_eof__[(ne << 1) - 1] = wnmax;

/* ****  initialize the counter for perturbations in optical depth */

    ntau_pd__ = 0;

/* ****   increment optical depth and optical depth perturbation */
/*       counter, ntau_pd for each radiatively active gas that is a */
/*       variable part of the state vector */

    i__1 = ngases;
    for (ng = 1; ng <= i__1; ++ng) {
	if ((i__2 = istate[ng + 1], abs(i__2)) == 3) {
	    ++ntau_pd__;
	    igs[ntau_pd__ - 1] = ng;
	}
	i__2 = ngastyp[ng - 1];
	for (ngt = 1; ngt <= i__2; ++ngt) {
	    ++ne;
	    wn_eof__[(ne << 1) - 2] = -999.;
	    wn_eof__[(ne << 1) - 1] = -999.;
	    s_copy(io_file__ + (ne - 1) * 132, gasfile + (ngt + ng * 3 - 4) * 
		    132, (ftnlen)132, (ftnlen)132);
/* L1801: */
	}
/* L1821: */
    }

/* ****   define aerosol optical properties and vertical distribution */

    i__1 = nmodes;
    for (m = 1; m <= i__1; ++m) {
	++ne;
	s_copy(io_file__ + (ne - 1) * 132, miefile + (m - 1) * 132, (ftnlen)
		132, (ftnlen)132);
	wn_eof__[(ne << 1) - 2] = -999.;
	wn_eof__[(ne << 1) - 1] = -999.;
	mode = m;
	iopen = 1;

/* *****       if this particle mode is a variable part of the state */
/*            vector, update the state vector counter and optical depth */
/*            counter */

	if ((i__2 = istate[m + ngases + 1], abs(i__2)) == 4) {
	    ++ntau_pd__;
	}

/* ****      read wavelength-dependent optical properties of aerosols. */

	readmie_(&mode, &iumie, &iuaer, miefile + (m - 1) * 132, &iang[m - 1],
		 &ioffmie[m - 1], &icwlq[m - 1], &icqext[m - 1], &icqsca[m - 
		1], &icg1[m - 1], &icg2[m - 1], &icf[m - 1], &ioffmom[m - 1], 
		&nstr, &io_end__[ne - 1], &io_err__[ne - 1], &wnaer[m - 1], 
		wneof, (ftnlen)132);

	wn_eof__[(ne << 1) - 2] = wneof[0];
	wn_eof__[(ne << 1) - 1] = wneof[1];

/* ****      read the aerosol vertical structure */

	cldstr_(&nlev, &iopen, &mode, &iuaer, aerfile + (m - 1) * 132, &
		iofftau[m - 1], &iztau[m - 1], &icptau[m - 1], &ictau[m - 1], 
		&scptau[m - 1], &sctau[m - 1], &lcsh[m - 1], &iccsh[m - 1], &
		sccsh[m - 1], p, alt, z__, dtauaer, (ftnlen)132);

/* L2221: */
    }

/* ****   print the integrated cloud optical depth */

    s_wsfe(&io___161);
    do_fio(&c__1, " Total Cloud and Aerosol Optical Depth at reference wavel"
	    "ength:", (ftnlen)63);
    do_fio(&c__1, "    alt(km)     p (bar)      T(K)       dtau        tau", (
	    ftnlen)55);
    e_wsfe();
    tau_tot__ = 0.f;
    i__1 = nlyr;
    for (k = 1; k <= i__1; ++k) {
	dtau_tot__ = 0.f;
	i__2 = nmodes;
	for (m = 1; m <= i__2; ++m) {
	    dtau_tot__ += dtauaer[m + k * 10 - 11];
/* L2441: */
	}
	pbar = p[k] * 1e-5f;
	tau_tot__ += dtau_tot__;
	s_wsfe(&io___165);
	do_fio(&c__1, (char *)&alt[k], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&pbar, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&t[k], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dtau_tot__, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&tau_tot__, (ftnlen)sizeof(real));
	e_wsfe();
/* L2461: */
    }

/* ****   read surface albedos and bi-directional reflection functions */

    nalbmx = 32766;
    iopen = 1;
    ++ne;
    wn_eof__[(ne << 1) - 2] = -999.;
    wn_eof__[(ne << 1) - 1] = -999.;
    s_copy(io_file__ + (ne - 1) * 132, surfile, (ftnlen)132, (ftnlen)132);

    readsur_(surfile, &iusur, &iopen, &ifrmsur, frmsur, &ioffsur, &nalbmx, &
	    iwnsur, &icwnsur, &nref, icalb, &scwalb, scalb, &wnmin, &wnmax, &
	    io_end__[ne - 1], &io_err__[ne - 1], wneof, (ftnlen)132, (ftnlen)
	    40);

    wn_eof__[(ne << 1) - 2] = wneof[0];
    wn_eof__[(ne << 1) - 1] = wneof[1];

/* ****   read solar flux input file. */

    if (lsolar) {
	++ne;
	wn_eof__[(ne << 1) - 2] = -999.f;
	wn_eof__[(ne << 1) - 1] = -999.f;
	iopen = 1;
	s_copy(io_file__ + (ne - 1) * 132, solfile, (ftnlen)132, (ftnlen)132);

	readsol_(solfile, &iusol1, &iusol2, &iopen, &ifrms0, frms0, &ioffs0, &
		ixs0, &icwns0, &icsol, &scwns0, &scsol, &iuins0, &wnmin, &
		wnmax, &au, &io_end__[ne - 1], &io_err__[ne - 1], wneof, (
		ftnlen)132, (ftnlen)40);

	wn_eof__[(ne << 1) - 2] = wneof[0];
	wn_eof__[(ne << 1) - 1] = wneof[1];

    } else {
	nza = 1;
    }

/* ****  save input file names and other data in a header file. */

    if (ifrmout <= 2) {

/* ****      find gaussian points and weights. These quantities are */
/*          needed to define the number of upward streams, nzup, */
/*          and the number of downward streams, nzdn,and */
/*          for flux integration for grey calculations */

	i__1 = nstr / 2;
	qgausn_(&i__1, umu_o__, gwt_o__);

/* ****      restack the values into order used by disort */

	i__1 = nstr / 2;
	for (nze = 1; nze <= i__1; ++nze) {
	    umu_1__[nze + nstr / 2 - 1] = umu_o__[nze - 1];
	    umu_1__[nze - 1] = -umu_o__[nstr / 2 - nze];
/* L2601: */
	}

	smart_hdr__(&iuout, &iutrn, &iuflx, atmfile, aerfile, miefile, 
		solfile, surfile, mixfile, gasfile, &iunits, &nmodes, &ncomp, 
		icomp, volmix, &ts, &au, &wgtatm, &wnmin, &wnmax, &isptype, &
		islit, &width, &dwn, &tauerr, &pi0err, &phferr, &surferr, &
		nstr, &numu, umu_1__, &nphi, phi, &nza, umu0, phi0, &levout, &
		nlout, pout, &accur, &lamber, &isource, &irad, &ifrmout, &
		radius, &sgrav, igas, &ngases, ngastyp, (ftnlen)132, (ftnlen)
		132, (ftnlen)132, (ftnlen)132, (ftnlen)132, (ftnlen)132, (
		ftnlen)132);

    }

/* ****   initialze th number of sources of extinction */

    next = 0;

/* ****  Call the spectral mapping Algorithm */

    smart_(&iugrp, &iuthrm, iuabc, &iuaer, &iusur, &iusol0, &iusol2, &iuout, &
	    iuheat, &iustat, &iutrn, &iuflx, &nmodes, &ngases, ngastyp, igtrn,
	     igpab, gasfile, &iref, &nref, &nlev, &nlyr, &next, k_out__, 
	    dp_dp__, &ncomp, icomp, volmix, &radius, &wgtatm, &ts, grav, p, t,
	     alt, z__, rmix, dtauaer, &wnmin, &wnmax, &wn_tol__, &isptype, &
	    islit, &width, &dwn, &usrang, &nstr, &numu, umu, &nphi, phi, &nza,
	     umu0, phi0, &levout, &nlout, &ws, &phiw, &lamber, &lplanck, &
	    lsolar, &irad, &iunits, &ifrmout, io_end__, io_err__, wn_eof__, &
	    accur, &ratm, io_file__, &tauerr, &pi0err, &phferr, &surferr, 
	    aid_lr__, dirsoflx, dnsoflx, upsoflx, dnthflx, upthflx, soheat, 
	    thheat, pd_frac__, &nstate, istate, &ntau_pd__, igs, &iu_pd__, 
	    iutpd, iuspd, iupdrad, (ftnlen)132, (ftnlen)132);

/* ****      Print final output */

    s_wsfe(&io___180);
    do_fio(&c__1, " O u t p u t    R a d i a n c e    A n g l e s    U s e d"
	    ": ", (ftnlen)59);
    do_fio(&c__1, "   i   ang(deg)     cos(ang)", (ftnlen)28);
    e_wsfe();
    pi = acos(-1.f);
    i__1 = numu;
    for (n = 1; n <= i__1; ++n) {
	ang = acos(umu[n - 1]) * 180.f / pi;
	s_wsfe(&io___184);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ang, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&umu[n - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L7001: */
    }

/* ****    spectral ranges used for each wavelength-dependent file: */

    s_wsfe(&io___185);
    do_fio(&c__1, " I/O Status and Final Wavenumber for Each Constituent:", (
	    ftnlen)54);
    do_fio(&c__1, "       n     io_end  io_err", (ftnlen)27);
    do_fio(&c__1, " wavenumber extent (cm**-1)     filename", (ftnlen)40);
    e_wsfe();
    i__1 = ngases + 1;
    for (ne = 1; ne <= i__1; ++ne) {

	charsp_(io_file__ + (ne - 1) * 132, name__, &len, &c__132, &nlb, &ntb,
		 (ftnlen)132, (ftnlen)1);

	s_wsfe(&io___190);
	do_fio(&c__1, (char *)&ne, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&io_end__[ne - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&io_err__[ne - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&wn_eof__[(ne << 1) - 2], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&wn_eof__[(ne << 1) - 1], (ftnlen)sizeof(
		doublereal));
	i__2 = len;
	for (ii = 1; ii <= i__2; ++ii) {
	    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
	}
	e_wsfe();
/* L7101: */
    }
    i__1 = next;
    for (ne = ngases + 2; ne <= i__1; ++ne) {

	charsp_(io_file__ + (ne - 1) * 132, name__, &len, &c__132, &nlb, &ntb,
		 (ftnlen)132, (ftnlen)1);

	s_wsfe(&io___192);
	do_fio(&c__1, (char *)&ne, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&io_end__[ne - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&io_err__[ne - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&wn_eof__[(ne << 1) - 2], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&wn_eof__[(ne << 1) - 1], (ftnlen)sizeof(
		doublereal));
	i__2 = len;
	for (ii = 1; ii <= i__2; ++ii) {
	    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
	}
	e_wsfe();
/* L7121: */
    }

    if (lsolar) {
	ne = next + 1;

	charsp_(io_file__ + (ne - 1) * 132, name__, &len, &c__132, &nlb, &ntb,
		 (ftnlen)132, (ftnlen)1);

	s_wsfe(&io___193);
	do_fio(&c__1, (char *)&ne, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&io_end__[ne - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&io_err__[ne - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&wn_eof__[(ne << 1) - 2], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&wn_eof__[(ne << 1) - 1], (ftnlen)sizeof(
		doublereal));
	i__1 = len;
	for (ii = 1; ii <= i__1; ++ii) {
	    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
	}
	e_wsfe();
    }

    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int smart_spectra__ () { MAIN__ (); return 0; }
