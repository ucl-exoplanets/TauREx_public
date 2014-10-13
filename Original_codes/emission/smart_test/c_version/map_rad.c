/* map_rad.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int map_rad__(logical *lsolar, logical *lplanck, logical *
	lamber, logical *usrang, integer *iugrp, integer *iuthrm, integer *
	ifl, integer *iend_ne__, integer *ne0, integer *nz0, integer *nphi, 
	integer *nstr, integer *numu, integer *nzup, integer *nzdn, integer *
	nza, integer *nlyr, integer *nlay, integer *nlout, integer *levout, 
	integer *k_out__, integer *iref, integer *ibin, integer *nsiout, 
	integer *igs, integer *npd, integer *npd1, integer *ipd, integer *
	ntau_pd__, integer *nt_pd__, integer *modepd, integer *irad, 
	doublereal *wnmin, doublereal *wnmax, doublereal *wn, doublereal *
	wn_io__, doublereal *wnout, doublereal *wn_tol__, doublereal *dvi, 
	real *umu0nz, real *phi0nz, real *solflx, real *umu_f__, real *
	gwt_f__, real *umu, real *phi, real *ts, real *t, real *p, real *
	dp_dp__, real *rmix, real *alb, real *dtauaer, real *tau_ext__, real *
	pd_pert__, real *alb_b__, real *dtau_b__, real *copi0_b__, real *
	g_b__, real *p_ray_0__, real *p_gas_0__, real *p_aer_0__, real *
	p_ray0__, real *p_gas0__, real *p_aer0__, real *tau_ray_0__, real *
	tau_gas_0__, real *tau_aer_0__, doublereal *trnflx_b__, doublereal *
	dtrnflxdx, doublereal *refflx_b__, doublereal *drefflxdx, doublereal *
	absflx_b__, doublereal *dabsflxdx, doublereal *refrad_b__, doublereal 
	*drefraddx, doublereal *absrad_b__, doublereal *dabsraddx, doublereal 
	*brdf_b__, doublereal *dbrdfdx, doublereal *dnsflxsrc_b__, doublereal 
	*ddnsflxdx, doublereal *upsflxsrc_b__, doublereal *dupsflxdx, 
	doublereal *dntflxsrc_b__, doublereal *ddntflxdx, doublereal *
	uptflxsrc_b__, doublereal *duptflxdx, doublereal *tradsrc_b__, 
	doublereal *dtraddx, doublereal *sradsrc_b__, doublereal *dsraddx, 
	real *up_s_flx__, real *dn_s_flx__, real *dir_s_flx__, real *
	up_t_flx__, real *dn_t_flx__, real *dns_src__, real *ups_src__, real *
	pd_dns_src__, real *pd_ups_src__, real *dnt_src__, real *upt_src__, 
	real *pd_dnt_src__, real *pd_upt_src__, real *pd_rad__, real *
	pd_trndir__, real *pd_trnflx__, real *pd_refflx__, real *pd_absflx__, 
	real *up_flx__, real *dn_flx__, real *dir_flx__, real *rad, real *
	trn_dir__, real *trn_flx__, real *ref_flx__, real *abs_flx__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12, i__13;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void)
	    ;
    double exp(doublereal);
    integer s_wsue(cilist *), e_wsue(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real tau_ray0__[2], dtau_aer__, dsol_rad__[7680]	/* was [16][
	    16][3][10] */;
    static doublereal absflx_i__[140]	/* was [70][2] */, dn_s_src__[70], 
	    dn_t_src__[70], refflx_i__[140]	/* was [70][2] */;
    static real dtau_gas__, dnsflx_i__[1400]	/* was [70][2][10] */, 
	    dntflx_i__[140]	/* was [70][2] */, dtau_ray__;
    static doublereal trndir_i__[1400]	/* was [70][2][10] */, wnio_old__, 
	    up_s_src__[70];
    static real ddirsflx[700]	/* was [70][10] */;
    static doublereal up_t_src__[70], trnflx_i__[140]	/* was [70][2] */;
    static real upsflx_i__[1400]	/* was [70][2][10] */, uptflx_i__[140]
	    	/* was [70][2] */, surf_opt__[4];
    static doublereal absflx1_i__[1400]	/* was [70][10][2] */, dsol_rad1__[
	    5376000]	/* was [16][16][3][70][10][10] */, dn_s_src0__[70], 
	    dn_t_src0__[70], dn_s_src1__[14000]	/* was [70][10][2][10] */, 
	    dn_t_src1__[1400]	/* was [70][10][2] */, refflx1_i__[1400]	
	    /* was [70][10][2] */;
    static real g[70];
    static doublereal trndir1_i__[14000]	/* was [70][10][2][10] */;
    static integer k, l, m, n;
    static doublereal up_s_src0__[70], up_t_src0__[70], up_s_src1__[14000]	
	    /* was [70][10][2][10] */, up_t_src1__[1400]	/* was [70][
	    10][2] */, trnflx1_i__[1400]	/* was [70][10][2] */, 
	    bb_flx_dn__[70], dabsflx_i__[70], ddn_s_src__[700]	/* was [70][
	    10] */, rad_s_src__[17920]	/* was [16][16][70] */, bb_flx_up__[
	    70], rad_t_src__[17920]	/* was [16][16][70] */, ddn_t_src__[
	    70], drefflx_i__[70];
    static real g0[70];
    static doublereal dtrndir_i__[700]	/* was [70][10] */, dup_s_src__[700]	
	    /* was [70][10] */, dup_t_src__[70];
    static real dirsflx_i__[1400]	/* was [70][2][10] */;
    static doublereal dtrnflx_i__[70], bb_flx_dn0__[70], dabsflx1_i__[700]	
	    /* was [70][10] */, rad_s_src0__[17920]	/* was [16][16][70] */
	    , bb_flx_up0__[70], rad_t_src0__[17920]	/* was [16][16][70] */
	    , ddn_s_src1__[7000]	/* was [70][10][10] */, bb[70], 
	    ddn_t_src1__[700]	/* was [70][10] */, drefflx1_i__[700]	/* 
	    was [70][10] */, dtrndir1_i__[7000]	/* was [70][10][10] */, 
	    dup_s_src1__[7000]	/* was [70][10][10] */;
    static integer ii, ne;
    static doublereal dup_t_src1__[700]	/* was [70][10] */;
    static integer ni;
    static doublereal dtrnflx1_i__[700]	/* was [70][10] */, dv;
    static integer nn;
    static real dx[700]	/* was [70][10] */;
    static integer nz;
    static doublereal dn_s_src_i__[1400]	/* was [70][2][10] */, 
	    dn_t_src_i__[140]	/* was [70][2] */, bb0[70];
    static integer l_1__, l_2__;
    static doublereal up_s_src_i__[1400]	/* was [70][2][10] */;
    static integer ng0;
    static doublereal up_t_src_i__[140]	/* was [70][2] */;
    static integer ni0;
    extern /* Subroutine */ int interp_rad__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, real *, real *, real *, real *
	    , real *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static real dns[3];
    static integer naz, nze;
    static real ups[3], alb0;
    extern /* Subroutine */ int int_rad_lev__(logical *, logical *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static real dns0[60]	/* was [3][2][10] */, dnt0[6]	/* was [3][2] 
	    */;
    extern /* Subroutine */ int perturb_rad__(logical *, logical *, logical *,
	     logical *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), grey_eq_trn__(logical *, logical *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, doublereal *, doublereal *, doublereal *, doublereal *), 
	    init_interp__(integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, real *, real *, real *, real *, 
	    real *, real *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static real ups0[60]	/* was [3][2][10] */, upt0[6]	/* was [3][2] 
	    */, dalb, delg[70], paer[16], ddns[30]	/* was [3][10] */, 
	    pgas[16], ddnt[3];
    static doublereal d_wn__;
    static real dtau[70], dnth[3], dirs[3];
    static integer nlev;
    static real tatm[70], pray[16], dups[30]	/* was [3][10] */, dupt[3], 
	    upth[3], dalb0, delg0[70], paer0[32]	/* was [16][2] */, 
	    copi0[70], pgas0[32]	/* was [16][2] */, dirs0[60]	/* 
	    was [3][2][10] */, pray0[32]	/* was [16][2] */, fbeam, 
	    g_sca__[700]	/* was [70][10] */, rad_s__[768]	/* 
	    was [16][16][3] */, dpaer[16], dpgas[16], ddirs[30]	/* was [3][10]
	     */, btemp, dvi_s__, dwn_s__, dpray[16], temis;
    static doublereal wnavg;
    static real ttemp, tsurf, paer_0__, co_pi0__[140]	/* was [70][2] */, 
	    pgas_0__, delpi0[70], dpgas0, dpaer0, pray_0__, dpray0;
    static doublereal bb_rad__[17920]	/* was [16][16][70] */;
    static real paer_00__[2], pgas_00__[2];
    static doublereal absrad[1120]	/* was [16][70] */;
    static real delpi00[70], th_rad__[17920]	/* was [16][16][70] */;
    static doublereal refrad[1120]	/* was [16][70] */;
    static real rad_th__[768]	/* was [16][16][3] */, pray_00__[2];
    static doublereal bb_sur__;
    static real deltau[70];
    static doublereal absflx[70], refflx[70], bb_rad0__[17920]	/* was [16][
	    16][70] */;
    static real dtauex[140]	/* was [70][2] */;
    static doublereal trnrad[1120]	/* was [16][70] */;
    static real dnsflx[70], dntflx[70];
    static doublereal trndir[70];
    static integer istart;
    static doublereal absrad0[1120]	/* was [16][70] */, trnflx[70];
    static real upsflx[70], uptflx[70], th_rad0__[1536]	/* was [16][16][3][2] 
	    */;
    static doublereal refrad0[1120]	/* was [16][70] */, th_rad1__[1075200]
	    	/* was [16][16][3][70][10][2] */, bb_sur0__;
    static real deltau0[70];
    static doublereal absflx0[70], refflx0[70], trnrad0[1120]	/* was [16][
	    70] */, trndir0[70], trnflx0[70];
    extern /* Subroutine */ int find_pd__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static real dth_rad__[17920]	/* was [16][16][70] */;
    static doublereal abs_cmu__[1120]	/* was [16][70] */;
    static real tau_sca__[700]	/* was [70][10] */, tau_aer__, tau_gas__, 
	    sol_rad__[17920]	/* was [16][16][70] */;
    extern /* Subroutine */ int rad_trn__(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static real ddnsflx[700]	/* was [70][10] */, ddntflx[70], tau_ray__;
    static doublereal trn_cmu__[1120]	/* was [16][70] */;
    static real dirsflx[70], dupsflx[700]	/* was [70][10] */, duptflx[
	    70];
    static doublereal dth_rad1__[537600]	/* was [16][16][3][70][10] */;
    static real tau_aer0__[2], sol_rad0__[15360]	/* was [16][16][3][2][
	    10] */;
    static doublereal sol_rad1__[10752000]	/* was [16][16][3][70][10][2][
	    10] */;
    static real tau_gas0__[2];

    /* Fortran I/O blocks */
    static cilist io___70 = { 1, 0, 1, 0, 0 };
    static cilist io___83 = { 1, 0, 1, 0, 0 };
    static cilist io___85 = { 0, 6, 0, "(1a,5(1pe15.8),4i5)", 0 };
    static cilist io___191 = { 0, 0, 0, 0, 0 };
    static cilist io___192 = { 0, 0, 0, 0, 0 };
    static cilist io___193 = { 0, 6, 0, "(/,1a,1pe14.6)", 0 };



/* ccccccccccccccccccccccccc   m a p _ r a d  cccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    this subroutine reads the wavenumber, interval width, bin       cc */
/* c    bin index, and albedo of the next spectral point.               cc */
/* c                                                                    cc */
/* c    NOTE: this version has been modified to produce layer dependent cc */
/* c    flux source functions and their partial derivatives, as well    cc */
/* c    as the corresponding flux values.                               cc */
/* c                                                                    cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c      iugrp - unit number of file with binned radiances             cc */
/* c         ne - index of extinction source                            cc */
/* c        nz0 - index of this solar zenith angle (1 to nsol)          cc */
/* c    iend_ne - end of file flag (0 - okay, 1 - end of file)          cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c       nlyr - number of computational model layers                  cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       numu - number of zenith angles in input file                 cc */
/* c      wnmin - minimum wavenumber of desired spectral window         cc */
/* c      wnmax - maximum wavenumber of desired spectral window         cc */
/* c         wn - current wavenumber                                    cc */
/* c     uptflx - wn-dependent upward thermal flux                      cc */
/* c     dntflx - wn-dependent downward thermal flux                    cc */
/* c     th_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c     upsflx - wn-dependent upward solar flux                        cc */
/* c     dnsflx - wn-dependent downward diffuse + direct solar flux     cc */
/* c    dirsflx - wn-dependent downward direct solar flux               cc */
/* c    sol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c    tau_ray - wn-dependent rayleigh-scattering optical depth        cc */
/* c    tau_gas - wn-dependent gas optical depth                        cc */
/* c    tau_aer - wn-dependent aerosol extinction optical depth         cc */
/* c     pray_0 - pressure of normial-incident rayleigh optical depth   cc */
/* c              unity                                                 cc */
/* c     pgas_0 - pressure of normial-incident gas optical depth unity  cc */
/* c     paer_0 - pressure of normial-incident aerosol optical depth    cc */
/* c              unity                                                 cc */
/* c       pray - pressure of rayleigh optical depth unity along each   cc */
/* c              output stream                                         cc */
/* c       pgas - pressure of gas optical depth unity along each stream cc */
/* c       paer - pressure of aerosol optical depth unity along each    cc */
/* c              output stream                                         cc */
/* c                                                                    cc */
/* c    o u t p u t :                                                   cc */
/* c                                                                    cc */
/* c    uptflx_i - wn-dependent upward thermal flux                      cc */
/* c    dntflx_i - wn-dependent downward thermal flux                    cc */
/* c    th_rad0 - wn-depndent, angle-dependent thermal radiances        cc */
/* c    upsflx_i - wn-dependent upward solar flux                        cc */
/* c    dnsflx_i - wn-dependent downward diffuse+direct solar flux       cc */
/* c   dirsflx_i - wn-dependent downward direct solar flux               cc */
/* c   sol_rad0 - wn-depndent, angle-dependent solar radiances          cc */
/* c   tau_ray0 - wn-dependent rayleigh-scattering optical depth        cc */
/* c   tau_gas0 - wn-dependent gas optical depth                        cc */
/* c   tau_aer0 - wn-dependent aerosol extinction optical depth         cc */
/* c    pray_00 - pressure of normial-incident rayleigh optical depth   cc */
/* c              unity                                                 cc */
/* c    pgas_00 - pressure of normial-incident gas optical depth unity  cc */
/* c    paer_00 - pressure of normial-incident aerosol optical depth    cc */
/* c              unity                                                 cc */
/* c      pray0 - pressure of rayleigh optical depth unity along each   cc */
/* c              output stream                                         cc */
/* c      pgas0 - pressure of gas optical depth unity along each stream cc */
/* c      paer0 - pressure of aerosol optical depth unity along each    cc */
/* c              output stream                                         cc */
/* c    duptflx - wn-dependent upward thermal flux                      cc */
/* c    ddntflx - wn-dependent downward thermal flux                    cc */
/* c    dth_rad - wn-depndent, angle-dependent thermal radiances        cc */
/* c    dupsflx - wn-dependent upward solar flux                        cc */
/* c    ddnsflx - wn-dependent downward diffuse + direct solar flux     cc */
/* c   ddirsflx - wn-dependent downward direct solar flux               cc */
/* c   dsol_rad - wn-depndent, angle-dependent solar radiances          cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccc   m a p _ r a d  cccccccccccccccccccccccccccc */




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





/* ****   counters used in map_back */


/* ****   internal spectral flag for starting wavenumber loop */


/* ****   spectral binning parameters */


/* ****   spectral binning parameters */





/* ****   gaussian angles for flux integration */



/* ****   perturbations for finite derivative jacobians */


/* ****    aerosol optical depth for each particle mode */


/* ****    stored monochromatic optical properties */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   binned layer flux transmittances and absorptances */


/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */



/* ****  pressures and optical depths of tau = 1 level */


/* ****   temporary variables used in interp_rad */


/* ****   spectrally-dependent output flux and radiance from interp_rad */


/* ****     monochormatic fluxes at wavenumber wn_io */


/* ****   downward and upward solar fluxes gradients at wavenumber wn_io */



/* ****   monochormatic radiances and wavelength gradients at wn_io */
/*       at output levels, k_out */



/* ****    optical property differences for this wn */


/* ****    planck functions used in interp_rad */


/* ****    layer radiance transmission values for simplified adding method */


/* ****    layer flux transmission values for simplified adding method */


/* ****   interpolated solar and thermal source functions */


/* ****    baseline optical property differences for this wn */


/* ****    baseline planck functions used in interp_rad */


/* ****    baseline layer transmission values for simplified adding method */


/* ****   baseline flux and radiance source functions */


/* ****   variables for interpolating source functions to wn */



/* ****   wavelength gradients in source terms */



/* ****   flux source functions at wavenumber wn_io for perturbed state */



/* ****   flux source functions wavenumber derivatives */



/* ****    flux transmission and absorption needed for simplified */
/*        adding method at wavenumber wn_io */



/* ****    flux transmission and absorption needed for simplified */
/*        adding method at wavenumber wn_io */



/* ****   fluxes and radiances at output levels k_out */



/* ****    flux transmission and absorption functions at wavenuber, wn */


/* ****     perturrbed monochormatic radiance and wn interpolation values */



/* ****    flux transmission and absorption partial */
/*        derivatives for simplified adding method at wavenumber wn */


/* ****    optical depth variables for the output wavenumber grid */







/* *****   variables passed to int_rad_lev */


/* ****    output fluxes and radiances */



/* ****   layer-dependent source terms for solar and thermal fluxes */


/* ****   output flux and radiance partial derivatives */




/* ****    specify the solar zenith angle and gas index */

    /* Parameter adjustments */
    --abs_flx__;
    --ref_flx__;
    --trn_flx__;
    --trn_dir__;
    rad -= 273;
    --dir_flx__;
    --dn_flx__;
    --up_flx__;
    pd_absflx__ -= 71;
    pd_refflx__ -= 71;
    pd_trnflx__ -= 71;
    pd_trndir__ -= 71;
    pd_rad__ -= 54801;
    pd_upt_src__ -= 71;
    pd_dnt_src__ -= 71;
    --upt_src__;
    --dnt_src__;
    pd_ups_src__ -= 71;
    pd_dns_src__ -= 71;
    --ups_src__;
    --dns_src__;
    dn_t_flx__ -= 71;
    up_t_flx__ -= 71;
    dir_s_flx__ -= 771;
    dn_s_flx__ -= 771;
    up_s_flx__ -= 771;
    dsraddx -= 89873;
    sradsrc_b__ -= 18193;
    dtraddx -= 89873;
    tradsrc_b__ -= 18193;
    duptflxdx -= 351;
    uptflxsrc_b__ -= 71;
    ddntflxdx -= 351;
    dntflxsrc_b__ -= 71;
    dupsflxdx -= 351;
    upsflxsrc_b__ -= 71;
    ddnsflxdx -= 351;
    dnsflxsrc_b__ -= 71;
    dbrdfdx -= 273;
    brdf_b__ -= 1553;
    dabsraddx -= 4497;
    absrad_b__ -= 6737;
    drefraddx -= 4497;
    refrad_b__ -= 6737;
    dabsflxdx -= 281;
    absflx_b__ -= 421;
    drefflxdx -= 281;
    refflx_b__ -= 421;
    dtrnflxdx -= 281;
    trnflx_b__ -= 421;
    --p_aer0__;
    --p_gas0__;
    --p_ray0__;
    g_b__ -= 421;
    copi0_b__ -= 421;
    dtau_b__ -= 421;
    alb_b__ -= 61;
    --pd_pert__;
    tau_ext__ -= 281;
    dtauaer -= 11;
    --alb;
    rmix -= 71;
    --dp_dp__;
    --p;
    --t;
    --phi;
    --umu;
    --gwt_f__;
    --umu_f__;
    wnout -= 229;
    --modepd;
    --ipd;
    --igs;
    nsiout -= 229;
    --k_out__;

    /* Function Body */
    nz = *nz0;

    ne = *ne0;
    nlev = *nlyr + 1;

    if (nsiout[(ne + nz * 113 << 1) + 2] == 0) {

/* ****    initialize spectral interval counters */

	nsiout[(ne + nz * 113 << 1) + 1] = 1;
	nsiout[(ne + nz * 113 << 1) + 2] = 2;

/* ****   initialize spectral radiance interpolation variables */

	*ibin = -1;
	*iend_ne__ = 0;
	istart = 0;
	*wn_io__ = -9999.;

	init_interp__(&nlev, nlout, numu, nphi, nz0, wnmin, &wnout[229], 
		upsflx, dnsflx, dirsflx, sol_rad__, uptflx, dntflx, th_rad__, 
		up_s_src_i__, dn_s_src_i__, dup_s_src__, ddn_s_src__, 
		up_t_src_i__, dn_t_src_i__, dup_t_src__, ddn_t_src__, 
		upsflx_i__, dnsflx_i__, dirsflx_i__, dupsflx, ddnsflx, 
		ddirsflx, sol_rad0__, dsol_rad__, uptflx_i__, dntflx_i__, 
		duptflx, ddntflx, th_rad0__, dth_rad__, ups, dns, dirs, upth, 
		dnth, rad_s__, rad_th__, pray0, pgas0, paer0, &pray_0__, &
		pgas_0__, &paer_0__, &tau_ray__, &tau_gas__, &tau_aer__, 
		dpray, dpgas, dpaer, tau_ray0__, tau_gas0__, tau_aer0__, 
		pray_00__, pgas_00__, paer_00__, &dtau_ray__, &dtau_gas__, &
		dtau_aer__, &dpray0, &dpgas0, &dpaer0, pgas, pray, paer);

    }

/* ****  end of initialization block */

/* ****  determine if new spectral radiances and fluxes are needed */
/*      (ifl = 1) or if the existing ones span the range of including */
/*      wavenumber wn, and can be used to interpolate the values */

    if (*ifl == 1) {

/* ****   read the binned optical properties at the next wavelength */

/* ****    s p e c t r a l    i n t e r v a l    l o o p */

L2002:
	ni = nsiout[(ne + nz * 113 << 1) + 1];
	nsiout[(ne + nz * 113 << 1) + 1] = nsiout[(ne + nz * 113 << 1) + 2];
	nsiout[(ne + nz * 113 << 1) + 2] = ni;

/* ****     read the next wavenumber, group index, and surface albedo. */

	wnio_old__ = *wn_io__;
	io___70.ciunit = *iugrp;
	i__1 = s_rsue(&io___70);
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = do_uio(&c__1, (char *)&(*wn_io__), (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = do_uio(&c__1, (char *)&(*ibin), (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L100001;
	}
	for (ii = 1; ii <= 4; ++ii) {
	    i__1 = do_uio(&c__1, (char *)&surf_opt__[ii - 1], (ftnlen)sizeof(
		    real));
	    if (i__1 != 0) {
		goto L100001;
	    }
	}
	i__2 = *nza;
	for (nn = 1; nn <= i__2; ++nn) {
	    i__1 = do_uio(&c__1, (char *)&alb[nn], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100001;
	    }
	}
	i__3 = *nt_pd__;
	for (l = 1; l <= i__3; ++l) {
	    i__4 = *nlay;
	    for (k = 1; k <= i__4; ++k) {
		i__1 = do_uio(&c__1, (char *)&dtauex[k + l * 70 - 71], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L100001;
		}
	    }
	}
	i__5 = *nt_pd__;
	for (l = 1; l <= i__5; ++l) {
	    i__6 = *nlay;
	    for (k = 1; k <= i__6; ++k) {
		i__1 = do_uio(&c__1, (char *)&co_pi0__[k + l * 70 - 71], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L100001;
		}
	    }
	}
	i__7 = *nlyr;
	for (k = 1; k <= i__7; ++k) {
	    i__1 = do_uio(&c__1, (char *)&g[k - 1], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100001;
	    }
	}
	i__8 = *ntau_pd__;
	for (n = 1; n <= i__8; ++n) {
	    for (m = 1; m <= 3; ++m) {
		i__9 = nlev;
		for (k = 1; k <= i__9; ++k) {
		    i__1 = do_uio(&c__1, (char *)&tau_ext__[k + (m + n * 3) * 
			    70], (ftnlen)sizeof(real));
		    if (i__1 != 0) {
			goto L100001;
		    }
		}
	    }
	}
	i__10 = *ntau_pd__;
	for (n = 1; n <= i__10; ++n) {
	    i__11 = nlev;
	    for (k = 1; k <= i__11; ++k) {
		i__1 = do_uio(&c__1, (char *)&tau_sca__[k + n * 70 - 71], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L100001;
		}
	    }
	}
	i__12 = *ntau_pd__;
	for (n = 1; n <= i__12; ++n) {
	    i__13 = nlev;
	    for (k = 1; k <= i__13; ++k) {
		i__1 = do_uio(&c__1, (char *)&g_sca__[k + n * 70 - 71], (
			ftnlen)sizeof(real));
		if (i__1 != 0) {
		    goto L100001;
		}
	    }
	}
	i__1 = e_rsue();
L100001:
	if (i__1 < 0) {
	    goto L2601;
	}
	if (i__1 > 0) {
	    goto L6001;
	}
	if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
	    io___83.ciunit = *iugrp;
	    i__1 = s_rsue(&io___83);
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__1 = do_uio(&c__1, (char *)&pray_0__, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__2 = *numu;
	    for (nze = *nzup; nze <= i__2; ++nze) {
		i__1 = do_uio(&c__1, (char *)&pray[nze - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L100002;
		}
	    }
	    i__1 = do_uio(&c__1, (char *)&tau_ray__, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__1 = do_uio(&c__1, (char *)&pgas_0__, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__4 = *numu;
	    for (nze = *nzup; nze <= i__4; ++nze) {
		i__1 = do_uio(&c__1, (char *)&pgas[nze - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L100002;
		}
	    }
	    i__1 = do_uio(&c__1, (char *)&tau_gas__, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__1 = do_uio(&c__1, (char *)&paer_0__, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__3 = *numu;
	    for (nze = *nzup; nze <= i__3; ++nze) {
		i__1 = do_uio(&c__1, (char *)&paer[nze - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		    goto L100002;
		}
	    }
	    i__1 = do_uio(&c__1, (char *)&tau_aer__, (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L100002;
	    }
	    i__1 = e_rsue();
L100002:
	    if (i__1 < 0) {
		goto L2601;
	    }
	    if (i__1 > 0) {
		goto L6001;
	    }
	}

	goto L2621;

L2601:
	*iend_ne__ = 1;
	*wn_io__ = *wnmax;
	ni = nsiout[(ne + nz * 113 << 1) + 2];
	s_wsfe(&io___85);
	do_fio(&c__1, "End of iugrp file encountered:", (ftnlen)30);
	do_fio(&c__1, (char *)&(*wn), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&wnio_old__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*wnmax), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&wnout[(ne + nz * 113 << 1) + 1], (ftnlen)
		sizeof(doublereal));
	do_fio(&c__1, (char *)&wnout[(ne + nz * 113 << 1) + 2], (ftnlen)
		sizeof(doublereal));
	do_fio(&c__1, (char *)&nz, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) + 1], (ftnlen)
		sizeof(integer));
	do_fio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) + 2], (ftnlen)
		sizeof(integer));
	do_fio(&c__1, (char *)&ni, (ftnlen)sizeof(integer));
	e_wsfe();

	return 0;

/* *****    set the radiance wavelength marker to wn_io (which is */
/*         presumably beyond the current wavelength, wn) */

L2621:
	wnout[ni + (ne + nz * 113 << 1)] = *wn_io__;

/* ***      define wavenumber offset of wn */

	d_wn__ = *wn - wnout[(ne + nz * 113 << 1) + 1];
	dwn_s__ = (real) d_wn__;

/* ****     find the spectral interval width */

	dv = wnout[(ne + nz * 113 << 1) + 2] - wnout[(ne + nz * 113 << 1) + 1]
		;
	wnavg = (wnout[(ne + nz * 113 << 1) + 2] + wnout[(ne + nz * 113 << 1) 
		+ 1]) * .5;

	if (abs(dv) >= *wn_tol__ * wnavg) {
	    *dvi = 1. / dv;
	} else {
	    *dvi = 0.;
	}
	dvi_s__ = (real) (*dvi);

/* ****     load tau=1 arrays for spectral interpolation */

	if (nz == 1) {
	    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
		tau_gas0__[ni - 1] = tau_gas__;
		tau_ray0__[ni - 1] = tau_ray__;
		tau_aer0__[ni - 1] = tau_aer__;
		pgas_00__[ni - 1] = pgas_0__;
		pray_00__[ni - 1] = pray_0__;
		paer_00__[ni - 1] = paer_0__;
		i__1 = *numu;
		for (nze = 1; nze <= i__1; ++nze) {
		    pgas0[nze + (ni << 4) - 17] = pgas[nze - 1];
		    pray0[nze + (ni << 4) - 17] = pray[nze - 1];
		    paer0[nze + (ni << 4) - 17] = paer[nze - 1];
/* L2801: */
		}

/* ****        find optical depth and p(tau=1) spectral gradients */

		dtau_gas__ = (tau_gas0__[1] - tau_gas0__[0]) * dvi_s__;
		dtau_ray__ = (tau_ray0__[1] - tau_ray0__[0]) * dvi_s__;
		dtau_aer__ = (tau_aer0__[1] - tau_aer0__[0]) * dvi_s__;
		dpgas0 = (pgas_00__[1] - pgas_00__[0]) * dvi_s__;
		dpray0 = (pray_00__[1] - pray_00__[0]) * dvi_s__;
		dpaer0 = (paer_00__[1] - paer_00__[0]) * dvi_s__;
		i__1 = *numu;
		for (nze = 1; nze <= i__1; ++nze) {
		    dpgas[nze - 1] = (pgas0[nze + 15] - pgas0[nze - 1]) * 
			    dvi_s__;
		    dpray[nze - 1] = (pray0[nze + 15] - pray0[nze - 1]) * 
			    dvi_s__;
		    dpaer[nze - 1] = (paer0[nze + 15] - paer0[nze - 1]) * 
			    dvi_s__;
/* L2821: */
		}

	    }

	}

/* ****         evaluate the radiances and partial derivatives */
/*             at this wavenumber for each optical property */
/*             that can vary */

	i__1 = *nlyr;
	for (k = 1; k <= i__1; ++k) {
	    tatm[k - 1] = t[k];
	    dtau[k - 1] = dtauex[k - 1];
	    copi0[k - 1] = co_pi0__[k - 1];
	    g0[k - 1] = g[k - 1];
/* L3001: */
	}

	tatm[nlev - 1] = t[nlev];
	tsurf = *ts;
	alb0 = alb[nz];

/* ****       set layer index range, l_1 and L_2 */

	l_1__ = 1;
	l_2__ = nlev;
	ng0 = *ibin;
	if (*lsolar) {

/* ****       find the direct beam transmission for sun */

	    i__1 = *nlyr;
	    for (k = 1; k <= i__1; ++k) {
		trndir[k - 1] = exp(-dtau[k - 1] / *umu0nz);
/* L3011: */
	    }

	}

/* ****  set the surface pressure interpolation index */

	if (*ibin > 0) {

/* ****       calculate radiances for wavenumber wn_io */

	    interp_rad__(lsolar, lplanck, &ng0, &l_1__, &l_2__, &nlev, nz0, 
		    nphi, numu, nzup, nzdn, wn_io__, umu0nz, &umu[1], &
		    alb_b__[61], &dtau_b__[421], &copi0_b__[421], &g_b__[421],
		     &trnflx_b__[421], &dtrnflxdx[281], &refflx_b__[421], &
		    drefflxdx[281], &absflx_b__[421], &dabsflxdx[281], &
		    refrad_b__[6737], &drefraddx[4497], &absrad_b__[6737], &
		    dabsraddx[4497], &brdf_b__[1553], &dbrdfdx[273], &tsurf, 
		    tatm, &alb0, dtau, copi0, g0, &dalb, deltau, delpi0, delg,
		     &dnsflxsrc_b__[71], &ddnsflxdx[351], &upsflxsrc_b__[71], 
		    &dupsflxdx[351], &dntflxsrc_b__[71], &ddntflxdx[351], &
		    uptflxsrc_b__[71], &duptflxdx[351], &sradsrc_b__[18193], &
		    dsraddx[89873], &tradsrc_b__[18193], &dtraddx[89873], bb, 
		    &bb_sur__, bb_flx_up__, bb_flx_dn__, bb_rad__, dn_s_src__,
		     up_s_src__, rad_s_src__, dn_t_src__, up_t_src__, 
		    rad_t_src__, trndir, trnflx, refflx, absflx, trnrad, 
		    refrad, absrad, upsflx, dnsflx, dirsflx, sol_rad__, 
		    uptflx, dntflx, th_rad__);

	} else {

/* ****       use the pure absorption model to calculate the radiances */
/*           at wn_io */

	    fbeam = 1.f;
	    temis = 1.f;
	    btemp = tsurf;
	    ttemp = 0.f;
	    l_1__ = 1;
	    l_2__ = nlev;

	    rad_trn__(usrang, lplanck, &l_1__, &l_2__, &ng0, nstr, numu, &umu[
		    1], &umu_f__[1], &gwt_f__[1], dtau, copi0, g0, trnrad, 
		    absrad, trn_cmu__, abs_cmu__);

	    grey_eq_trn__(lsolar, lplanck, usrang, lamber, &nlev, nz0, nstr, 
		    numu, nphi, nzup, nzdn, iref, wn_io__, dtau, copi0, tatm, 
		    &umu[1], &umu_f__[1], &gwt_f__[1], &phi[1], umu0nz, 
		    phi0nz, &fbeam, &alb0, &btemp, &ttemp, &temis, trndir, 
		    trnflx, absflx, trn_cmu__, abs_cmu__, trnrad, absrad, 
		    upsflx, dnsflx, dirsflx, sol_rad__, uptflx, dntflx, 
		    th_rad__, dn_s_src__, up_s_src__, dn_t_src__, up_t_src__);

	}

/* ****  load arrays for wavelength interpolation */

	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {
	    trnflx_i__[k + ni * 70 - 71] = trnflx[k - 1];
	    refflx_i__[k + ni * 70 - 71] = refflx[k - 1];
	    absflx_i__[k + ni * 70 - 71] = absflx[k - 1];
/* L3021: */
	}

	if (*lsolar) {
	    trndir_i__[(ni + (nz << 1)) * 70 - 210] = trndir[0];
	    i__1 = nlev;
	    for (k = 2; k <= i__1; ++k) {
		trndir_i__[k + (ni + (nz << 1)) * 70 - 211] = trndir[k - 1] * 
			trndir_i__[k - 1 + (ni + (nz << 1)) * 70 - 211];
/* L3031: */
	    }
	}

/* ****     find the spectral gradients in the flux transmission and */
/*         absorption values used for the simplified adding method */

	i__1 = nlev;
	for (k = 1; k <= i__1; ++k) {

/* ****         define solar flux wavelength derivatives */

	    dtrndir_i__[k + nz * 70 - 71] = (trndir_i__[k + ((nz << 1) + 2) * 
		    70 - 211] - trndir_i__[k + ((nz << 1) + 1) * 70 - 211]) * 
		    *dvi;
	    dtrnflx_i__[k - 1] = (trnflx_i__[k + 69] - trnflx_i__[k - 1]) * *
		    dvi;
	    drefflx_i__[k - 1] = (refflx_i__[k + 69] - refflx_i__[k - 1]) * *
		    dvi;
	    dabsflx_i__[k - 1] = (absflx_i__[k + 69] - absflx_i__[k - 1]) * *
		    dvi;
/* L3101: */
	}

/* ****       interpolate flux and radiance values to their output levels */

	int_rad_lev__(lsolar, lplanck, nphi, numu, &nz, nlyr, nlout, levout, &
		k_out__[1], &dp_dp__[1], upsflx, dnsflx, dirsflx, uptflx, 
		dntflx, sol_rad__, th_rad__, ups, dns, dirs, upth, dnth, 
		rad_s__, rad_th__);

/* ****       find the spectral gradients in the fluxes and */
/*           flux source functions */


	if (*lsolar) {
	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {

/* ****         define layer dependent solar source terms at wn_io */

		up_s_src_i__[k + (ni + (nz << 1)) * 70 - 211] = up_s_src__[k 
			- 1];
		dn_s_src_i__[k + (ni + (nz << 1)) * 70 - 211] = dn_s_src__[k 
			- 1];

/* ****         define solar source term wavelength derivatives */

		dup_s_src__[k + nz * 70 - 71] = (up_s_src_i__[k + ((nz << 1) 
			+ 2) * 70 - 211] - up_s_src_i__[k + ((nz << 1) + 1) * 
			70 - 211]) * dvi_s__;
		ddn_s_src__[k + nz * 70 - 71] = (dn_s_src_i__[k + ((nz << 1) 
			+ 2) * 70 - 211] - dn_s_src_i__[k + ((nz << 1) + 1) * 
			70 - 211]) * dvi_s__;
/* L3121: */
	    }

	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {

/* ****           define solar fluxes at wn_io */

		upsflx_i__[k + (ni + (nz << 1)) * 70 - 211] = upsflx[k - 1];
		dnsflx_i__[k + (ni + (nz << 1)) * 70 - 211] = dnsflx[k - 1];
		dirsflx_i__[k + (ni + (nz << 1)) * 70 - 211] = dirsflx[k - 1];

/* ****           define solar flux wavelength derivatives */

		dupsflx[k + nz * 70 - 71] = (upsflx_i__[k + ((nz << 1) + 2) * 
			70 - 211] - upsflx_i__[k + ((nz << 1) + 1) * 70 - 211]
			) * dvi_s__;
		ddnsflx[k + nz * 70 - 71] = (dnsflx_i__[k + ((nz << 1) + 2) * 
			70 - 211] - dnsflx_i__[k + ((nz << 1) + 1) * 70 - 211]
			) * dvi_s__;
		ddirsflx[k + nz * 70 - 71] = (dirsflx_i__[k + ((nz << 1) + 2) 
			* 70 - 211] - dirsflx_i__[k + ((nz << 1) + 1) * 70 - 
			211]) * dvi_s__;
/* L3201: */
	    }

	    i__1 = *nlout;
	    for (k = 1; k <= i__1; ++k) {

/* ****           load fluxes at the output levels, k_out, and */
/*               define solar flux wavelength derivatives */

		ups0[k + (ni + (nz << 1)) * 3 - 10] = ups[k - 1];
		dups[k + nz * 3 - 4] = (ups0[k + ((nz << 1) + 2) * 3 - 10] - 
			ups0[k + ((nz << 1) + 1) * 3 - 10]) * dvi_s__;
		dns0[k + (ni + (nz << 1)) * 3 - 10] = dns[k - 1];
		ddns[k + nz * 3 - 4] = (dns0[k + ((nz << 1) + 2) * 3 - 10] - 
			dns0[k + ((nz << 1) + 1) * 3 - 10]) * dvi_s__;
		dirs0[k + (ni + (nz << 1)) * 3 - 10] = dirs[k - 1];
		ddirs[k + nz * 3 - 4] = (dirs0[k + ((nz << 1) + 2) * 3 - 10] 
			- dirs0[k + ((nz << 1) + 1) * 3 - 10]) * dvi_s__;
		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__4 = *numu;
		    for (nze = 1; nze <= i__4; ++nze) {

/* ****                   define solar radiances and wavelength */
/*                       derivatives at output levels, k_out */

			sol_rad0__[nze + (naz + (k + (ni + (nz << 1)) * 3 << 
				4) << 4) - 2577] = rad_s__[nze + (naz + (k << 
				4) << 4) - 273];
			dsol_rad__[nze + (naz + (k + nz * 3 << 4) << 4) - 
				1041] = (sol_rad0__[nze + (naz + (k + ((nz << 
				1) + 2) * 3 << 4) << 4) - 2577] - sol_rad0__[
				nze + (naz + (k + ((nz << 1) + 1) * 3 << 4) <<
				 4) - 2577]) * dvi_s__;
/* L3221: */
		    }
/* L3241: */
		}
/* L3261: */
	    }
	}

	if (*lplanck && nz == 1) {

	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {

/* ****         define layer dependent solar source terms at wn_io */

		up_t_src_i__[k + ni * 70 - 71] = up_t_src__[k - 1];
		dn_t_src_i__[k + ni * 70 - 71] = dn_t_src__[k - 1];

/* ****         define solar source term wavelength derivatives */

		dup_t_src__[k - 1] = (up_t_src_i__[k + 69] - up_t_src_i__[k - 
			1]) * dvi_s__;
		ddn_t_src__[k - 1] = (dn_t_src_i__[k + 69] - dn_t_src_i__[k - 
			1]) * dvi_s__;
/* L3321: */
	    }

	    i__1 = nlev;
	    for (k = 1; k <= i__1; ++k) {

/* ****           define thermal fluxes at wn_io */

		uptflx_i__[k + ni * 70 - 71] = uptflx[k - 1];
		dntflx_i__[k + ni * 70 - 71] = dntflx[k - 1];

/* ****           define thermal flux wavelength derivatives */

		duptflx[k - 1] = (uptflx_i__[k + 69] - uptflx_i__[k - 1]) * 
			dvi_s__;
		ddntflx[k - 1] = (dntflx_i__[k + 69] - dntflx_i__[k - 1]) * 
			dvi_s__;
/* L3401: */
	    }

	    i__1 = *nlout;
	    for (k = 1; k <= i__1; ++k) {

/* ****           load fluxes at the output levels, k_out, and */
/*               define solar flux wavelength derivatives */

		upt0[k + ni * 3 - 4] = upth[k - 1];
		dupt[k - 1] = (upt0[k + 2] - upt0[k - 1]) * dvi_s__;
		dnt0[k + ni * 3 - 4] = dnth[k - 1];
		ddnt[k - 1] = (dnt0[k + 2] - dnt0[k - 1]) * dvi_s__;

/* ****           define thermal radiance wavelength derivatives */

		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__4 = *numu;
		    for (nze = 1; nze <= i__4; ++nze) {

/* ****                   define thermal radiances and wavelength */
/*                       derivatives at output levels, k_out */

			th_rad0__[nze + (naz + (k + ni * 3 << 4) << 4) - 1041]
				 = rad_th__[nze + (naz + (k << 4) << 4) - 273]
				;
			dth_rad__[nze + (naz + (k << 4) << 4) - 273] = (
				th_rad0__[nze + (naz + (k + 6 << 4) << 4) - 
				1041] - th_rad0__[nze + (naz + (k + 3 << 4) <<
				 4) - 1041]) * dvi_s__;
/* L3421: */
		    }
/* L3441: */
		}
/* L3461: */
	    }

	}

	if (*npd > 0) {

/* ****       initialze optical properties, flux and radiance transmission */
/*           and absorption functions in each layer */

	    dalb0 = dalb;
	    i__1 = nlev - 1;
	    for (k = 1; k <= i__1; ++k) {
		deltau0[k - 1] = deltau[k - 1];
		delpi00[k - 1] = delpi0[k - 1];
		delg0[k - 1] = delg[k - 1];
		trndir0[k - 1] = trndir[k - 1];
		trnflx0[k - 1] = trnflx[k - 1];
		refflx0[k - 1] = refflx[k - 1];
		absflx0[k - 1] = absflx[k - 1];
		i__2 = *numu;
		for (nze = 1; nze <= i__2; ++nze) {
		    trnrad0[nze + (k << 4) - 17] = trnrad[nze + (k << 4) - 17]
			    ;
		    refrad0[nze + (k << 4) - 17] = refrad[nze + (k << 4) - 17]
			    ;
		    absrad0[nze + (k << 4) - 17] = absrad[nze + (k << 4) - 17]
			    ;
/* L3601: */
		}
/* L3631: */
	    }

	    if (*lsolar) {
		i__1 = nlev;
		for (k = 1; k <= i__1; ++k) {
		    dn_s_src0__[k - 1] = dn_s_src__[k - 1];
		    up_s_src0__[k - 1] = up_s_src__[k - 1];
		    i__2 = *nphi;
		    for (naz = 1; naz <= i__2; ++naz) {
			i__4 = *numu;
			for (nze = 1; nze <= i__4; ++nze) {
			    rad_s_src0__[nze + (naz + (k << 4) << 4) - 273] = 
				    rad_s_src__[nze + (naz + (k << 4) << 4) - 
				    273];
/* L3641: */
			}
/* L3661: */
		    }
/* L3681: */
		}
	    }

	    if (*lplanck) {
		bb_sur0__ = bb_sur__;
		i__1 = nlev;
		for (k = 1; k <= i__1; ++k) {
		    bb0[k - 1] = bb[k - 1];
		    bb_flx_up0__[k - 1] = bb_flx_up__[k - 1];
		    bb_flx_dn0__[k - 1] = bb_flx_dn__[k - 1];
		    dn_t_src0__[k - 1] = dn_t_src__[k - 1];
		    up_t_src0__[k - 1] = up_t_src__[k - 1];
		    i__2 = *numu;
		    for (nze = 1; nze <= i__2; ++nze) {
			i__4 = *nphi;
			for (naz = 1; naz <= i__4; ++naz) {
			    bb_rad0__[nze + (naz + (k << 4) << 4) - 273] = 
				    bb_rad__[nze + (naz + (k << 4) << 4) - 
				    273];
			    rad_t_src0__[nze + (naz + (k << 4) << 4) - 273] = 
				    rad_t_src__[nze + (naz + (k << 4) << 4) - 
				    273];
/* L3701: */
			}
/* L3721: */
		    }
/* L3781: */
		}
	    }

/* ****        find flux and radiance jacobians at wavenumer wn_io */

	    ni0 = ni;

	    perturb_rad__(lsolar, lplanck, lamber, usrang, &ni0, nz0, nlyr, 
		    iref, nlout, levout, &k_out__[1], ibin, nphi, nstr, numu, 
		    nzup, nzdn, npd, npd1, &ipd[1], &igs[1], &modepd[1], 
		    wn_io__, dvi, &umu_f__[1], &gwt_f__[1], &umu[1], &phi[1], 
		    umu0nz, phi0nz, &pd_pert__[1], dx, &alb[1], ts, &t[1], &p[
		    1], &dp_dp__[1], &rmix[71], &dtauaer[11], &tau_ext__[281],
		     tau_sca__, g_sca__, dtauex, co_pi0__, g, &alb_b__[61], &
		    dtau_b__[421], &copi0_b__[421], &g_b__[421], &trnflx_b__[
		    421], &dtrnflxdx[281], &refflx_b__[421], &drefflxdx[281], 
		    &absflx_b__[421], &dabsflxdx[281], &refrad_b__[6737], &
		    drefraddx[4497], &absrad_b__[6737], &dabsraddx[4497], &
		    brdf_b__[1553], &dbrdfdx[273], &dnsflxsrc_b__[71], &
		    upsflxsrc_b__[71], &sradsrc_b__[18193], &dntflxsrc_b__[71]
		    , &uptflxsrc_b__[71], &tradsrc_b__[18193], &ddntflxdx[351]
		    , &duptflxdx[351], &dtraddx[89873], &ddnsflxdx[351], &
		    dupsflxdx[351], &dsraddx[89873], &dalb0, deltau0, delpi00,
		     delg0, bb0, &bb_sur0__, bb_flx_up0__, bb_flx_dn0__, 
		    bb_rad0__, dn_s_src0__, up_s_src0__, rad_s_src0__, 
		    dn_t_src0__, up_t_src0__, rad_t_src0__, trndir0, trnflx0, 
		    refflx0, absflx0, trnrad0, refrad0, absrad0, trndir1_i__, 
		    trnflx1_i__, refflx1_i__, absflx1_i__, dtrndir1_i__, 
		    dtrnflx1_i__, drefflx1_i__, dabsflx1_i__, dn_s_src1__, 
		    up_s_src1__, dup_s_src1__, ddn_s_src1__, dn_t_src1__, 
		    up_t_src1__, dup_t_src1__, ddn_t_src1__, sol_rad1__, 
		    dsol_rad1__, th_rad1__, dth_rad1__);

	}

/* ******  if this is is the first point, get the second one for the */
/*        spectral interpolation */

	if (istart == 0) {
	    istart = 1;
	    goto L2002;
	}

/* ****     ensure that at least one valid radiance point beyond wn has */
/*         been found.  if not, read the next point. */

	if (wnout[ni + (ne + nz * 113 << 1)] < *wn && wnout[ni + (ne + nz * 
		113 << 1)] < *wnmax && *iend_ne__ != 1) {
	    goto L2002;
	}

    }

    if (nz == 1) {
	if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {

/* ****       save optical depths and pressure of tau = 1 */

	    *tau_gas_0__ = tau_gas0__[0] + dtau_gas__ * dwn_s__;
	    *tau_ray_0__ = tau_ray0__[0] + dtau_ray__ * dwn_s__;
	    *tau_aer_0__ = tau_aer0__[0] + dtau_aer__ * dwn_s__;
	    *p_gas_0__ = pgas_00__[0] + dpgas0 * dwn_s__;
	    *p_ray_0__ = pray_00__[0] + dpray0 * dwn_s__;
	    *p_aer_0__ = paer_00__[0] + dpaer0 * dwn_s__;
	    i__1 = *numu;
	    for (nze = 1; nze <= i__1; ++nze) {
		p_gas0__[nze] = pgas0[nze - 1] + dpgas[nze - 1] * dwn_s__;
		p_ray0__[nze] = pray0[nze - 1] + dpray[nze - 1] * dwn_s__;
		p_aer0__[nze] = paer0[nze - 1] + dpaer[nze - 1] * dwn_s__;
/* L4001: */
	    }

	}

    }

/* ****           interpolate flux transmision and absorption functions */
/*               to wavenumber, wn */

    i__1 = nlev;
    for (k = 1; k <= i__1; ++k) {
	if (*lsolar) {
	    trn_dir__[k] = (real) (trndir_i__[k + ((nz << 1) + 1) * 70 - 211] 
		    + dtrndir_i__[k + nz * 70 - 71] * d_wn__);
	}
	trn_flx__[k] = (real) (trnflx_i__[k - 1] + dtrnflx_i__[k - 1] * 
		d_wn__);
	ref_flx__[k] = (real) (refflx_i__[k - 1] + drefflx_i__[k - 1] * 
		d_wn__);
	abs_flx__[k] = (real) (absflx_i__[k - 1] + dabsflx_i__[k - 1] * 
		d_wn__);
/* L4021: */
    }

    if (*lsolar) {

/* ****     Interpolate solar fluxes, flux source functions, and */
/*          radiances to wavenumber, wn */

	i__1 = *nlyr + 1;
	for (k = 1; k <= i__1; ++k) {

/* ****        interpolate solar source functions to this wn */

	    ups_src__[k] = (real) (up_s_src_i__[k + ((nz << 1) + 1) * 70 - 
		    211] + dup_s_src__[k + nz * 70 - 71] * d_wn__);
	    dns_src__[k] = (real) (dn_s_src_i__[k + ((nz << 1) + 1) * 70 - 
		    211] + ddn_s_src__[k + nz * 70 - 71] * d_wn__);

/* ****        interpolate solar fluxes and radiances to this wn */
/*            Note: fluxes must be retained at wn adjacent wavenumbers */
/*            for the trapezoid integration over wavenumber in the */
/*            heating rate calculation */

	    up_s_flx__[k + (nz + 20) * 70] = *solflx * (upsflx_i__[k + ((nz <<
		     1) + 1) * 70 - 211] + dupsflx[k + nz * 70 - 71] * 
		    dwn_s__);
	    dn_s_flx__[k + (nz + 20) * 70] = *solflx * (dnsflx_i__[k + ((nz <<
		     1) + 1) * 70 - 211] + ddnsflx[k + nz * 70 - 71] * 
		    dwn_s__);
	    dir_s_flx__[k + (nz + 20) * 70] = *solflx * (dirsflx_i__[k + ((nz 
		    << 1) + 1) * 70 - 211] + ddirsflx[k + nz * 70 - 71] * 
		    dwn_s__);

	    if (up_s_flx__[k + (nz + 20) * 70] < 0.f) {
		up_s_flx__[k + (nz + 20) * 70] = 0.f;
	    }
	    if (dn_s_flx__[k + (nz + 20) * 70] < 0.f) {
		dn_s_flx__[k + (nz + 20) * 70] = 0.f;
	    }
	    if (dir_s_flx__[k + (nz + 20) * 70] < 0.f) {
		dir_s_flx__[k + (nz + 20) * 70] = 0.f;
	    }
/* L4101: */
	}

/* ****        interpolate solar fluxes and radiances at levels nlout */
/*            to this wn */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {

/* ****         find the solar fluxes at levels k_out */

	    ups[k - 1] = *solflx * (ups0[k + ((nz << 1) + 1) * 3 - 10] + dups[
		    k + nz * 3 - 4] * dwn_s__);
	    dns[k - 1] = *solflx * (dns0[k + ((nz << 1) + 1) * 3 - 10] + ddns[
		    k + nz * 3 - 4] * dwn_s__);
	    dirs[k - 1] = *solflx * (dirs0[k + ((nz << 1) + 1) * 3 - 10] + 
		    ddirs[k + nz * 3 - 4] * dwn_s__);

/* ****         find radiances at levels k_out */

	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__4 = *numu;
		for (nze = 1; nze <= i__4; ++nze) {
		    rad_s__[nze + (naz + (k << 4) << 4) - 273] = *solflx * (
			    sol_rad0__[nze + (naz + (k + ((nz << 1) + 1) * 3 
			    << 4) << 4) - 2577] + dsol_rad__[nze + (naz + (k 
			    + nz * 3 << 4) << 4) - 1041] * dwn_s__);

		    if (rad_s__[nze + (naz + (k << 4) << 4) - 273] < 0.f) {
			rad_s__[nze + (naz + (k << 4) << 4) - 273] = 0.f;
		    }
/* L4121: */
		}
/* L4141: */
	    }
/* L4161: */
	}

    }

/* ****       t h e r m a l    f l u x e s    a n d    r a d i a n c e s */
/*           (do only for nz = 1) */

    if (*lplanck && nz == 1) {

/* ****     interpolate thermal fluxes and radiances to this wn */

	i__1 = *nlyr + 1;
	for (k = 1; k <= i__1; ++k) {

/* ****        interpolate solar source functions to this wn */

	    upt_src__[k] = (real) (up_t_src_i__[k - 1] + dup_t_src__[k - 1] * 
		    d_wn__);
	    dnt_src__[k] = (real) (dn_t_src_i__[k - 1] + ddn_t_src__[k - 1] * 
		    d_wn__);

/* ****     interpolate thermal fluxes to this wn */

	    up_t_flx__[k + 140] = uptflx_i__[k - 1] + duptflx[k - 1] * 
		    dwn_s__;
	    dn_t_flx__[k + 140] = dntflx_i__[k - 1] + ddntflx[k - 1] * 
		    dwn_s__;

	    if (up_t_flx__[k + 140] < 0.f) {
		up_t_flx__[k + 140] = 0.f;
	    }
	    if (dn_t_flx__[k + 140] < 0.f) {
		dn_t_flx__[k + 140] = 0.f;
	    }

/* L4201: */
	}

/* ****     interpolate thermal fluxes and radiances at levels */
/*         k_out to this wn */

	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {

/* ****         interpolate thermal fluxes at levels k_out to this wn */

	    upth[k - 1] = upt0[k - 1] + dupt[k - 1] * dwn_s__;
	    dnth[k - 1] = dnt0[k - 1] + ddnt[k - 1] * dwn_s__;
	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__4 = *numu;
		for (nze = 1; nze <= i__4; ++nze) {
		    rad_th__[nze + (naz + (k << 4) << 4) - 273] = th_rad0__[
			    nze + (naz + (k + 3 << 4) << 4) - 1041] + 
			    dth_rad__[nze + (naz + (k << 4) << 4) - 273] * 
			    dwn_s__;
		    if (rad_th__[nze + (naz + (k << 4) << 4) - 273] < 0.f) {
			rad_th__[nze + (naz + (k << 4) << 4) - 273] = 0.f;
		    }
/* L4221: */
		}
/* L4241: */
	    }
/* L4261: */
	}

	io___191.ciunit = *iuthrm;
	s_wsue(&io___191);
	do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) + 1], (ftnlen)
		sizeof(integer));
	do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) + 2], (ftnlen)
		sizeof(integer));
	do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) + 1] + (ne + 
		nz * 113 << 1)], (ftnlen)sizeof(doublereal));
	do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) + 2] + (ne + 
		nz * 113 << 1)], (ftnlen)sizeof(doublereal));
	do_uio(&c__1, (char *)&(*wn_io__), (ftnlen)sizeof(doublereal));
	i__1 = *nlout;
	for (k = 1; k <= i__1; ++k) {
	    do_uio(&c__1, (char *)&upth[k - 1], (ftnlen)sizeof(real));
	    do_uio(&c__1, (char *)&dnth[k - 1], (ftnlen)sizeof(real));
	    i__2 = *nphi;
	    for (naz = 1; naz <= i__2; ++naz) {
		i__4 = *numu;
		for (nze = 1; nze <= i__4; ++nze) {
		    do_uio(&c__1, (char *)&rad_th__[nze + (naz + (k << 4) << 
			    4) - 273], (ftnlen)sizeof(real));
		}
	    }
	}
	e_wsue();

    } else {

/* ****      nz ne 1.  Read stored values */

	if (*lplanck) {
	    io___192.ciunit = *iuthrm;
	    s_rsue(&io___192);
	    do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) + 1], (ftnlen)
		    sizeof(integer));
	    do_uio(&c__1, (char *)&nsiout[(ne + nz * 113 << 1) + 2], (ftnlen)
		    sizeof(integer));
	    do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) + 1] + (
		    ne + nz * 113 << 1)], (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&wnout[nsiout[(ne + nz * 113 << 1) + 2] + (
		    ne + nz * 113 << 1)], (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&(*wn_io__), (ftnlen)sizeof(doublereal));
	    i__4 = *nlout;
	    for (k = 1; k <= i__4; ++k) {
		do_uio(&c__1, (char *)&upth[k - 1], (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&dnth[k - 1], (ftnlen)sizeof(real));
		i__2 = *nphi;
		for (naz = 1; naz <= i__2; ++naz) {
		    i__1 = *numu;
		    for (nze = 1; nze <= i__1; ++nze) {
			do_uio(&c__1, (char *)&rad_th__[nze + (naz + (k << 4) 
				<< 4) - 273], (ftnlen)sizeof(real));
		    }
		}
	    }
	    e_rsue();
	}

    }

/* ****    c o m b i n e d   f l u x e s   a n d   r a d i a n c e s */

/* ****     define the combined solar+thermal upward and downward */
/*         fluxes at the output levels */

    i__1 = *nlout;
    for (k = 1; k <= i__1; ++k) {

/* ****               store fluxes at top or bottom of atmosphere: */

	up_flx__[k] = upth[k - 1] + ups[k - 1];
	dn_flx__[k] = dnth[k - 1] + dns[k - 1];
	dir_flx__[k] = dirs[k - 1];

/* ****          define the combined solar+thermal radiances at the */
/*              output level, */

	i__2 = *nphi;
	for (naz = 1; naz <= i__2; ++naz) {
	    i__4 = *numu;
	    for (nze = 1; nze <= i__4; ++nze) {
		rad[nze + (naz + (k << 4) << 4)] = rad_s__[nze + (naz + (k << 
			4) << 4) - 273] + rad_th__[nze + (naz + (k << 4) << 4)
			 - 273];
/* L4601: */
	    }
/* L4621: */
	}
/* L4641: */
    }

    if (*npd > 0) {

/* ****      find partial derivatives */

	ni0 = ni;

	find_pd__(lsolar, lplanck, nz0, nlyr, nlout, nphi, numu, npd, &ipd[1],
		 wn, &d_wn__, dx, solflx, &rad[273], &trn_dir__[1], &
		trn_flx__[1], &ref_flx__[1], &abs_flx__[1], &dns_src__[1], &
		ups_src__[1], &dnt_src__[1], &upt_src__[1], trndir1_i__, 
		trnflx1_i__, refflx1_i__, absflx1_i__, dtrndir1_i__, 
		dtrnflx1_i__, drefflx1_i__, dabsflx1_i__, dn_s_src1__, 
		up_s_src1__, ddn_s_src1__, dup_s_src1__, dn_t_src1__, 
		up_t_src1__, ddn_t_src1__, dup_t_src1__, sol_rad1__, 
		dsol_rad1__, th_rad1__, dth_rad1__, &pd_trndir__[71], &
		pd_trnflx__[71], &pd_refflx__[71], &pd_absflx__[71], &
		pd_dns_src__[71], &pd_ups_src__[71], &pd_dnt_src__[71], &
		pd_upt_src__[71], &pd_rad__[54801]);

    }

    return 0;
L6001:
    s_wsfe(&io___193);
    do_fio(&c__1, "Input Error on unit iugrp in subroutine map_rad after: wn"
	    " =", (ftnlen)59);
    do_fio(&c__1, (char *)&wnio_old__, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_stop("", (ftnlen)0);

    return 0;
} /* map_rad__ */

