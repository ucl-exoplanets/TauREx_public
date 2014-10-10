/* smart.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int smart_(integer *iugrp, integer *iuthrm, integer *iuabc, 
	integer *iuaer, integer *iusur, integer *iusol0, integer *iusol2, 
	integer *iuout, integer *iuheat, integer *iustat, integer *iutrn, 
	integer *iuflx, integer *nmodes, integer *ngases, integer *ngastyp, 
	integer *igtrn, integer *igpab, char *gasfile, integer *iref, integer 
	*nref, integer *nlev, integer *nlyr, integer *next, integer *k_out__, 
	real *dp_dp__, integer *ncomp, integer *icomp, real *volmix, real *
	radius, real *wgtatm, real *ts, real *grav, real *p, real *t, real *
	alt, real *z__, real *rmix, real *dtauaer, doublereal *wnmin, 
	doublereal *wnmax, doublereal *wn_tol__, integer *isptype, integer *
	islit, doublereal *width, doublereal *dwn, logical *usrang, integer *
	nstr, integer *numu, real *umu, integer *nphi, real *phi, integer *
	nza, real *umu0, real *phi0, integer *levout, integer *nlout, real *
	ws, real *phiw, logical *lamber, logical *lplanck, logical *lsolar, 
	integer *irad, integer *iunits, integer *ifrmout, integer *io_end__, 
	integer *io_err__, doublereal *wn_eof__, real *accur, real *ratm, 
	char *io_file__, real *tauerr, real *pi0err, real *phferr, real *
	surferr, doublereal *aid_lr__, doublereal *dirsoflx, doublereal *
	dnsoflx, doublereal *upsoflx, doublereal *dnthflx, doublereal *
	upthflx, doublereal *soheat, doublereal *thheat, real *pd_frac__, 
	integer *nstate, integer *istate, integer *ntau_pd__, integer *igs, 
	integer *iu_pd__, integer *iutpd, integer *iuspd, integer *iupdrad, 
	ftnlen gasfile_len, ftnlen io_file_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12;
    real r__1, r__2;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double acos(doublereal), sqrt(doublereal);
    integer f_clos(cllist *), f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_wsue(cilist *), e_wsue(void), i_nint(real *), f_rew(alist *);

    /* Local variables */
    static integer itau1cnt;
    static doublereal absrad_b__[2867200]	/* was [16][70][5][512] */;
    extern /* Subroutine */ int map_back__(logical *, logical *, logical *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal refrad_b__[2867200]	/* was [16][70][5][512] */, 
	    absflx_b__[179200]	/* was [70][5][512] */, refflx_b__[179200]	
	    /* was [70][5][512] */;
    static real tau_scai__[2100]	/* was [70][10][3] */, tau_aeri__[3];
    static doublereal trnrad_b__[2867200]	/* was [16][70][5][512] */;
    static real tau_gasi__[3];
    extern /* Subroutine */ int init_bin__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *), 
	    cpu_time__(real *);
    static real bcputime, tau_rayi__[3];
    static doublereal trnflx_b__[179200]	/* was [70][5][512] */;
    static real tau_exti__[6300]	/* was [70][3][10][3] */;
    static doublereal dwn_step__;
    static real rcputime, surf_opt__[4], small_pi0__, d__, g[70];
    static integer k, l, m, n;
    static real x, small_alb__;
    static doublereal dabsraddx[1720320]	/* was [16][70][3][512] */, 
	    sradsrc_b__[9175040]	/* was [16][16][70][512] */, 
	    tradsrc_b__[9175040]	/* was [16][16][70][512] */, 
	    drefraddx[1720320]	/* was [16][70][3][512] */;
    static real h0;
    extern /* Subroutine */ int map_spect__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, doublereal *, doublereal *
	    , real *, real *, real *, real *, real *, real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer l0;
    static doublereal dabsflxdx[107520]	/* was [70][3][512] */;
    static real small_tau__, t0;
    static doublereal drefflxdx[107520]	/* was [70][3][512] */;
    extern /* Subroutine */ int sm_eq_trn__(logical *, logical *, logical *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal dtrnraddx[1720320]	/* was [16][70][3][512] */, 
	    ddnsflxdx[143360]	/* was [70][4][512] */, ddntflxdx[143360]	
	    /* was [70][4][512] */;
    extern /* Subroutine */ int layer_trn__(logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static real smcputime, surf_opti__[12]	/* was [4][3] */;
    static doublereal dtrnflxdx[107520]	/* was [70][3][512] */, dupsflxdx[
	    143360]	/* was [70][4][512] */, duptflxdx[143360]	/* 
	    was [70][4][512] */;
    static real gi[210]	/* was [70][3] */;
    static integer ii, ne, nf, ng, ij;
    static real pi;
    static doublereal wn;
    static integer nz, nc0, ij0, ng0;
    static real pi0;
    static doublereal wn0, wn1;
    static integer nz0;
    static real g_b__[179200]	/* was [70][5][512] */, alb[10], eff;
    extern /* Subroutine */ int norm_pi0_ph__(integer *, integer *, integer *,
	     real *, real *, real *, real *, real *, real *);
    static integer nlb, ipd[10], len, ntb, npg[150]	/* was [3][50] */, 
	    npd, ngt, nze, mom, ntg[150]	/* was [3][50] */;
    static doublereal wni[3];
    static real rtd;
    extern /* Subroutine */ int load_optics__(logical *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, real *, real *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static integer ipi0;
    static doublereal dnsflxsrc_b__[35840]	/* was [70][512] */, 
	    dntflxsrc_b__[35840]	/* was [70][512] */;
    static real sza0;
    static doublereal upsflxsrc_b__[35840]	/* was [70][512] */, 
	    uptflxsrc_b__[35840]	/* was [70][512] */;
    static real albi[30]	/* was [10][3] */;
    static integer ieff;
    static char name__[1*132];
    static integer icpu;
    static real ggrp[107520]	/* was [70][512][3] */;
    static integer nlay, nmom, nzdn;
    static doublereal dist[113];
    static real tcpu;
    static integer nzup;
    static real aerg0[20]	/* was [2][10] */, dg0dv[10], time0, time1, 
	    time2, time3;
    static doublereal dist0;
    static real alb_b__[25600]	/* was [10][5][512] */;
    static integer nlyr0, numu0;
    static real g_sca__[700]	/* was [70][10] */;
    static integer n_rad__[2560]	/* was [5][512] */;
    static real p_aer__[16], p_gas__[16];
    static integer nbgas[150]	/* was [3][50] */;
    static real btime;
    static integer igcnt, nt_pd__, nscat;
    static doublereal delnu;
    static real gwt_f__[16], umu_f__[16];
    static integer ntaug;
    static real rtime, gwt_o__[16];
    static integer ngrey, igtyp[150]	/* was [3][50] */;
    static doublereal dwnmx;
    static real taumn, umu_o__[16], umu0_1__[10], units, p_ray__[16], 
	    co_pi0__[140]	/* was [70][2] */;
    static doublereal wntot;
    static real phmom[14070]	/* was [201][70] */, sur_b__[10240]	/* 
	    was [4][5][512] */;
    static doublereal wnext[226]	/* was [2][113] */, wngrp[1536]	/* 
	    was [512][3] */;
    static real surf_0__[8]	/* was [2][4] */, pi0grp[107520]	/* 
	    was [70][512][3] */, phi0nz;
    static doublereal brdf_b__[655360]	/* was [16][16][5][512] */;
    static real dx_b_i__[179200]	/* was [70][5][512] */, g_scai__[2100]
	    	/* was [70][10][3] */;
    extern /* Subroutine */ int albedo_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, logical *, real *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, real *, 
	    real *, real *, real *, real *, real *, doublereal *, integer *, 
	    integer *);
    static real umu0nz, p_aeri__[48]	/* was [16][3] */, dtau_b__[179200]	
	    /* was [70][5][512] */;
    extern doublereal cp_atm__(integer *, integer *, real *, real *);
    static real p_gasi__[48]	/* was [16][3] */, albgrp[15360]	/* 
	    was [512][10][3] */;
    static integer modepd[10], ibtime;
    static real pmom_b__[36019200]	/* was [201][70][5][512] */, tauabs, 
	    tausca;
    static integer iwnflg, irtime, iflext[113], nspect;
    static real dtausc[70], p_aer_0__;
    static integer iwngrp;
    static real taugas[70];
    static integer nsiext[226]	/* was [2][113] */, ngroup, ngrtot;
    static doublereal wnsmin, wnsmax;
    static real p_gas_0__, tauray[70], tauaer[70], dtauex[140]	/* was [70][2]
	     */, co_pi0i__[420]	/* was [70][2][3] */, phmomi[42210]	/* 
	    was [201][70][3] */, p_ray_0__, p_rayi__[48]	/* was [16][3]
	     */;
    static doublereal wnskip;
    static real copi0_b__[179200]	/* was [70][5][512] */, ss_min__, 
	    points[141], taugrp[107520]	/* was [70][512][3] */;
    static doublereal dnugrp[512];
    extern /* Subroutine */ int qgausn_(integer *, real *, real *), charsp_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer ipi0cnt;
    static logical usrang0;
    static integer levtau1[1024]	/* was [512][2] */;
    static doublereal dbrdfdx[131072]	/* was [16][16][512] */;
    static real small_g__, aerqsca[20]	/* was [2][10] */, tau_aer__, 
	    tau_sca__[700]	/* was [70][10] */;
    extern /* Subroutine */ int aer_tau__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, real *, real *, real *, real *, real *
	    , real *, real *, real *, doublereal *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *);
    static real tau_gas__, dqscadv[10], dtausca;
    static doublereal dsraddx[36700160]	/* was [16][16][70][4][512] */;
    static real ang_umu__, theight, dtausci[210]	/* was [70][3] */;
    static integer nmomaer[20]	/* was [2][10] */;
    static doublereal dtraddx[36700160]	/* was [16][16][70][4][512] */;
    static real aerpmom[4020]	/* was [2][201][10] */, pd_pert__[10];
    static integer iwnmode;
    static real dtauexi[420]	/* was [70][2][3] */, tau_ray__;
    extern /* Subroutine */ int ray_tau__(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, integer *, integer *, 
	    integer *, integer *, integer *, real *, doublereal *, doublereal 
	    *, real *, real *, real *, real *, real *, real *, real *);
    static real dpmomdv[2010]	/* was [201][10] */;
    static integer ismtime, itaucnt, nmom_mx__[10], iu_flux__;
    static doublereal distmin, distmax;
    static integer nmomgrp[512];
    static real tau_tot__, dsurfdv[4], aerqext[20]	/* was [2][10] */, 
	    dqextdv[10], tau_ext__[2100]	/* was [70][3][10] */;
    static integer ismterr, isurcnt;
    static real p_gas_0i__[3], p_aer_0i__[3], pmomgrp[21611520]	/* was [201][
	    70][512][3] */, surfgrp[6144]	/* was [512][4][3] */;
    static doublereal wn_step__;
    static real p_ray_0i__[3];
    extern /* Subroutine */ int gas_tau__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, doublereal *, real *, real *, 
	    real *, real *, real *, doublereal *, integer *, integer *), 
	    skip_wn__(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, doublereal *, 
	    real *, real *, real *, real *, real *, real *, real *, 
	    doublereal *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *), heating_(integer *, integer *, integer *, real *, 
	    doublereal *, doublereal *, real *, real *, real *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer ntau_pd0__;
    static real dalb_b_i__[5120]	/* was [10][512] */;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, "(/,1a,/,1a,/1a)", 0 };
    static cilist io___19 = { 0, 6, 0, "(2(1pe12.4))", 0 };
    static cilist io___30 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___33 = { 0, 6, 0, "(3(1pe12.4))", 0 };
    static cilist io___34 = { 0, 6, 0, "(/,1a,/,1a)", 0 };
    static cilist io___36 = { 0, 6, 0, "(3(1pe12.4))", 0 };
    static cilist io___60 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___61 = { 0, 6, 0, "(/,1a)", 0 };
    static cilist io___68 = { 0, 0, 0, 0, 0 };
    static cilist io___73 = { 0, 6, 0, "(5i5,2x,132a)", 0 };
    static cilist io___75 = { 0, 6, 0, "(3i5,2x,132a)", 0 };
    static cilist io___178 = { 0, 0, 0, 0, 0 };
    static cilist io___180 = { 0, 0, 0, 0, 0 };
    static cilist io___227 = { 0, 0, 0, "(/,/1a,/,/,4a,/,3a/)", 0 };
    static cilist io___232 = { 0, 0, 0, "(2(1pe11.4),i5,2i8,i5,2i10,3i9,5i8,"
	    "i10)", 0 };
    static cilist io___233 = { 0, 0, 0, "(/,1a,1pe15.5)", 0 };
    static cilist io___234 = { 0, 0, 0, "(1a,1pe15.5)", 0 };
    static cilist io___235 = { 0, 0, 0, "(1a,i10)", 0 };
    static cilist io___237 = { 0, 0, 0, "(1a,1pe12.4)", 0 };
    static cilist io___238 = { 0, 0, 0, "(/,1a,1pe12.4,1a)", 0 };



/* ccccccccccccccccccccccccccc   s m a r t  cccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c    p u r p o s e :                                                 cc */
/* c                                                                    cc */
/* c    This program combines spectral mapping methods and the          cc */
/* c    discrete ordinate model (Stamnes et al. 1988) to generate       cc */
/* c    high-resolution synthetic monochromatic radiances in            cc */
/* c    vertically-inhomogeneous scattering absorbing, emitting,        cc */
/* c    planetary atmospheres.                                          cc */
/* c                                                                    cc */
/* c    note:  This version of the program uses a full brdf and         cc */
/* c           surface optical propertiy gradients for interpolation.   cc */
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
/* c 12/04 -4/05: radiance partial derivatives found for optical depth  cc */
/* c              single scattering albedos, atmospheric temperatures,  cc */
/* c              surface pressures and surface albedos,                cc */
/* c    05/06/05: radiance partial derivatives for scattering phase     cc */
/* c              functions added.                                      cc */
/* c    05/08/05: redundant spectral grid points are skipped (skip_wn)  cc */
/* c    02/11/07: layer transmittance and absorptance routines changed  cc */
/* c              from a full-range interpolation to a bin-by-bin       cc */
/* c              interpolation based on disort trnmed, albmed          cc */
/* c    02/11/07: calls to flush and timing routines swapped out to     cc */
/* c              make the routines more compatible with f95.           cc */
/* c                                                                    cc */
/* c    i n p u t :                                                     cc */
/* c                                                                    cc */
/* c     iusol0 - unit number as a scratch file for solar fluxes        cc */
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
/* c     nstate - number of elements in the state vecror                cc */
/* c     istate - state vector flag indicating which state variables    cc */
/* c              are variable components of the state vector.          cc */
/* c              1 - surface pressure                                  cc */
/* c              2 - surface/atmospheric temperature                   cc */
/* c              3 - gas absorption coeffient                          cc */
/* c              4 - cloud/aerosol optical depth                       cc */
/* c              5 - surface albedo                                    cc */
/* c    ntau_pd - number of gases and aerosols that are variable        cc */
/* c              components of thestate vector that affect the optical cc */
/* c              depth.                                                cc */
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
/* c    surferr - surface optical property binning error                cc */
/* c       umu0 - cosine of solar zenith angles                         cc */
/* c       phi0 - solar azimuth angles (degrees)                        cc */
/* c      accur - azimuth convergence accuracy for D/O routine          cc */
/* c     lamber - Include a lambertian surface? (Logical: T/F)          cc */
/* c              note: if lamber = F, use a BRDF is used.              cc */
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
/* c       nref - number of surface optical properties specified at     cc */
/* c              each wavelength.                                      cc */
/* c       nstr - number of gaussian zenith angles used in D/O code     cc */
/* c       nphi - number of output azimuth angles                       cc */
/* c       nlev - number of output levels                               cc */
/* c       nlyr - number of computational model layers                  cc */
/* c     levout - output level index (1) top of atmosphere,             cc */
/* c              (2) surface, (3) arbitrary level                      cc */
/* c        nza - number of solar zenith angles                         cc */
/* c     lsolar - include solar fluxes? (logical: T/F)                  cc */
/* c    lplanck - include thermal fluxes? (logical: T/F)                cc */
/* c     usrang - output radiances at user angles? (logical: T/F)       cc */
/* c     iunits - index of output radiance units:                       cc */
/* c              1) Watts/m**2/sr/cm**-1                               cc */
/* c              2) Watts/m**2/sr/micron                               cc */
/* c              3) Watts/m**2/sr/nanometer                            cc */
/* c              4) Watts/m**2/sr/Angstrom                             cc */
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
/* c        rad - radiance at each output zenith and azimuth angle      cc */
/* c  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc */
/* c  trn_gas_0 - normal-incidence gas transmission                     cc */
/* c  trn_ray_0 - normal-incidence aerosol transmission                 cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccc   s m a r t  cccccccccccccccccccccccccccccc */




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




/* ****   input aerosol wavelength-dependent variables set in aer_tau */


/* ****   variables used for lbl gas absorption coefficients */


/* ****   index of optical depth unity for each bin */


/* ****    units for flux and radiance jacobians */


/* ****  spectral grid parameters */






/* ***    set up direction cosines for flux integration */



/* ****    pressure of tau = 1 and column optical depths */


/* ****   input aerosol wavelength-dependent variables set in aer_tau */


/* ****   atmospheric structure variables */


/* ****    gas mixing ratios */


/* ****   atmospheric optical properties */


/* ****    optical depths of each variable component of state vector */


/* ****   variables used in skip_wn */






/* ****   binned solar source function variables and jacobians */



/* ****   binned thermal source function variables and jacobians */



/* ****     monochormatic radiance interpolation variables */


/* ****   variables for single/multiple scattering test */


/* ****   partial derivative fractional change */


/* ****   number of points in slit function algorithm */


/* ****    variables for spectral grouping algorithm. */


/* ****   spectral grid monitoring varibles */



/* ****   double precision variables for heating rates */


/* ****   binned layer radiance transmittances and absorptances */


/* ****   binned layer flux transmittances and absorptances */


/* ****    initialize timing variables */

    /* Parameter adjustments */
    iupdrad -= 10513;
    iuspd -= 11;
    iutpd -= 11;
    --igs;
    --istate;
    --pd_frac__;
    --thheat;
    soheat -= 71;
    --upthflx;
    --dnthflx;
    upsoflx -= 71;
    dnsoflx -= 71;
    dirsoflx -= 71;
    --aid_lr__;
    io_file__ -= 132;
    wn_eof__ -= 3;
    --io_err__;
    --io_end__;
    --phi0;
    --umu0;
    --phi;
    --umu;
    dtauaer -= 11;
    rmix -= 71;
    --z__;
    --alt;
    --t;
    --p;
    --grav;
    --volmix;
    --icomp;
    --dp_dp__;
    --k_out__;
    gasfile -= 528;
    igpab -= 4;
    igtrn -= 4;
    --ngastyp;
    iuabc -= 4;

    /* Function Body */
    tcpu = 0.f;

/* ****   define a small value that is added to optical depths.  this */
/*       quantity is needed to prevent an instability that develops */
/*       in the thermal radiance calculation if there are layers */
/*       with  near-zero optical depth. */

    small_tau__ = 1e-5f;
/* Danie Liang */
    small_tau__ = 1e-15f;

/* ****    set small values of single scattering parameter and */
/*        asymmetry parameter for use in the spectral skipping */
/*        algorithm. */

    small_g__ = 1e-5f;
    small_pi0__ = 1e-5f;
    small_alb__ = 1e-5f;
/* Danie Liang */
    small_g__ = 0.f;
    small_pi0__ = 0.f;
    small_alb__ = 0.f;

/* ****   Define a minimum single scattering albedo for scattering */
/*       calculations */

    ss_min__ = *pi0err * .25f;

/* ****    define a maximum spectral interval between output points. */

    if (*isptype == 1) {
	dwnmx = (*wnmax - *wnmin) * .001;
    } else {
	dwnmx = *dwn;
    }

/* ****   set the output interval for the spectral grid monitor */

    dwn_step__ = (*wnmax - *wnmin) * .01;
    wn_step__ = *wnmin + dwn_step__;

/* ****   modify solar zenith angles to approximate the sphericity of the */
/*       planet.  To do this, assume that the planet is spherical with */
/*       a 3-scale-height thick atmosphere. Define the ray tangent height */
/*                     theight = sin(umu0). */
/*       The distance between the top of the atmosphere and the */
/*       axis of rotation, x, is then given by the pythagorean theorm: */
/*           d = sqrt((radius + H)**2 - rt**2) */
/*       and the solar path length between the top of the atmosphere and */
/*       the surface is given by: */
/*           x = d - radius*cos(umu0) */
/*       The effective solar zenith angle cosine is therefore */
/*       given by H/x. */

    if (*lsolar) {

/* ****    find the effective scale height: */

	s_wsfe(&io___10);
	do_fio(&c__1, "Effective Solar Zenith Angle: ", (ftnlen)30);
	do_fio(&c__1, "    sza0       sza_eff", (ftnlen)22);
	do_fio(&c__1, "  (degrees)   (degrees)", (ftnlen)23);
	e_wsfe();

	h0 = *ratm * .003f * t[*nlev] / grav[*nlev];
	rtd = 180.f / acos(-1.f);

/* ****    find the effective solar zenith angles: */

	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    sza0 = rtd * acos(umu0[nz]);
/* Computing 2nd power */
	    r__1 = umu0[nz];
	    theight = *radius * sqrt(1.f - r__1 * r__1);
/* Computing 2nd power */
	    r__1 = *radius + h0;
/* Computing 2nd power */
	    r__2 = theight;
	    d__ = sqrt(r__1 * r__1 - r__2 * r__2);
	    x = d__ - *radius * umu0[nz];
	    umu0_1__[nz - 1] = umu0[nz];
	    umu0[nz] = h0 / x;
	    if (umu0[nz] > 1.f) {
		umu0[nz] = 1.f;
	    }
	    s_wsfe(&io___19);
	    do_fio(&c__1, (char *)&sza0, (ftnlen)sizeof(real));
	    r__1 = rtd * acos(umu0[nz]);
	    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	    e_wsfe();
/* L1141: */
	}

    }

/* ****    compute the adiabatic lapse rate, g/cp, at each level. */
/*        this quantity is needed for heating rates. */

    nc0 = *ncomp;
    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {

/* ****    define the adiabatic lapse rate (k/m) */

	t0 = t[k];
	aid_lr__[k] = grav[k] / cp_atm__(&nc0, &icomp[1], &volmix[1], &t0);

/* L1201: */
    }

/* ****   initialize point counter in slit function routine */

    if (*isptype != 1) {
	for (nf = 1; nf <= 141; ++nf) {
	    points[nf - 1] = 0.f;
/* L1241: */
	}
    }

/* ****  find gaussian points and weights for flux integration */

    i__1 = *nstr / 2;
    qgausn_(&i__1, umu_o__, gwt_o__);

/* ****  restack the values into order used by disort */

    i__1 = *nstr / 2;
    for (nze = 1; nze <= i__1; ++nze) {
	umu_f__[nze + *nstr / 2 - 1] = umu_o__[nze - 1];
	gwt_f__[nze + *nstr / 2 - 1] = gwt_o__[nze - 1];
/* L1301: */
    }
    i__1 = *nstr / 2;
    for (nze = 1; nze <= i__1; ++nze) {
	umu_f__[nze - 1] = -umu_o__[*nstr / 2 - nze];
	gwt_f__[nze - 1] = -gwt_o__[*nstr / 2 - nze];
/* L1321: */
    }
    s_wsfe(&io___30);
    do_fio(&c__1, "Computational Gaussian Points and Weights", (ftnlen)41);
    do_fio(&c__1, "  Angle (deg)    umu       gwt", (ftnlen)30);
    e_wsfe();
    pi = acos(-1.f);
    i__1 = *nstr;
    for (nze = 1; nze <= i__1; ++nze) {
	ang_umu__ = acos(umu_f__[nze - 1]) * 180.f / pi;
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&ang_umu__, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&umu_f__[nze - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&gwt_f__[nze - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L1341: */
    }

/* ****   define the last downward and first upward looking angles */

    if (*usrang) {
	s_wsfe(&io___34);
	do_fio(&c__1, "User-Defined Emission Angles:", (ftnlen)29);
	do_fio(&c__1, "  Angle (deg)    umu ", (ftnlen)21);
	e_wsfe();
	nzdn = 0;
	i__1 = *numu;
	for (nze = 1; nze <= i__1; ++nze) {
	    if (umu[nze] < 0.f) {
		nzdn = nze;
	    }
	    ang_umu__ = acos(umu[nze]) * 180.f / pi;
	    s_wsfe(&io___36);
	    do_fio(&c__1, (char *)&ang_umu__, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&umu[nze], (ftnlen)sizeof(real));
	    e_wsfe();
/* L1361: */
	}
	nzup = nzdn + 1;
    } else {
	nzdn = *nstr / 2;
	nzup = nzdn + 1;

/* ****     initialize umu values for use in tau=1 calculation */

	i__1 = *numu;
	for (nze = 1; nze <= i__1; ++nze) {
	    umu[nze] = umu_f__[nze - 1];
/* L1381: */
	}
    }

/* ****   initialze timing routines: */

    cpu_time__(&time0);
    tcpu = 0.f;

/* ****   initialize wavenumbers */

    wn = *wnmin;
    wn0 = *wnmin;
    wn1 = *wnmin;
    wntot = 1.;
    wnsmin = *wnmin;
    wnsmax = *wnmin;
    wnskip = 0.;

/* ****    set a maximum step size */

    distmax = (*wnmax - *wnmin) * .01;
    if (distmax < *wnmin * .01) {
	distmax = *wnmin * .01;
    }

/* ****    initialize flux summation arrays for heating rate calculation */

    i__1 = *nlev;
    for (k = 1; k <= i__1; ++k) {
	upthflx[k] = 0.f;
	dnthflx[k] = 0.f;
	i__2 = *nza;
	for (nz = 1; nz <= i__2; ++nz) {
	    upsoflx[k + nz * 70] = 0.f;
	    dnsoflx[k + nz * 70] = 0.f;
	    dirsoflx[k + nz * 70] = 0.f;
/* L1501: */
	}
/* L1521: */
    }

/* ****    initialize spectral mapping variables */

    ngrtot = 0;
    ngroup = 0;

/* ****    define a small number as a minimum fractional tau error */

/* Computing 2nd power */
    r__1 = *tauerr;
    taumn = (r__1 * r__1 + small_tau__) * 1e-5f;

/* ****    determine whether pressure is a variable part of the state */
/*        vector.  If it is, gas optical depths are computed for the */
/*        usual nlyr layers, and for a a "perturnbed" lower layer, */
/*        which has 1% higher pressure. */

    if (abs(istate[1]) == 1) {
	nlay = *nlyr + 1;
    } else {
	nlay = *nlyr;
    }

/* *****   determine whether temperature is a variable part of */
/*        the state vector.  if it is, extinction optical depth must be */
/*        computed for the background and perterbed temperature profile */

    if (istate[2] == 0) {
	nt_pd__ = 1;
    } else {
	nt_pd__ = 2;
    }

/* ****     set jacobian indexes and fractional changes for each variable */
/*         component of the state vector */

    npd = 0;
    i__1 = *nstate;
    for (n = 1; n <= i__1; ++n) {
	if (istate[n] != 0) {
	    ++npd;
	    ipd[npd - 1] = istate[n];
	    pd_pert__[npd - 1] = pd_frac__[n];
	}
/* L1601: */
    }

/* ****   initialize absorption coefficient read flag, iflext */

    for (n = 1; n <= 113; ++n) {
	iflext[n - 1] = 1;
	nsiext[(n << 1) - 2] = 0;
	nsiext[(n << 1) - 1] = 0;
	wnext[(n << 1) - 1] = -999.;
	io_err__[n] = 0;
	io_end__[n] = 0;
/* L1701: */
    }

/* ****    open the scratch unit for group numbers at each wavenumber */

    cl__1.cerr = 0;
    cl__1.cunit = *iugrp;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = *iugrp;
    o__1.ofnm = 0;
    o__1.orl = 0;
    o__1.osta = "scratch";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);

/* ****    open the scratch unit for thermal fluxes at each wavenumber */

    cl__1.cerr = 0;
    cl__1.cunit = *iuthrm;
    cl__1.csta = 0;
    f_clos(&cl__1);
    o__1.oerr = 0;
    o__1.ounit = *iuthrm;
    o__1.ofnm = 0;
    o__1.orl = 0;
    o__1.osta = "scratch";
    o__1.oacc = 0;
    o__1.ofm = "unformatted";
    o__1.oblnk = 0;
    f_open(&o__1);

/* ****   set io file properties for rayleigh scattering */

    ne = 1;
    s_copy(io_file__ + ne * 132, "Rayleigh Scattering", (ftnlen)132, (ftnlen)
	    19);
    wn_eof__[(ne << 1) + 1] = *wnmin;
    wn_eof__[(ne << 1) + 2] = *wnmax;
    io_err__[ne] = 0;
    io_end__[ne] = 0;

/* ****      read header of each gas absorption line file: */

    if (*ngases > 0) {
	s_wsfe(&io___60);
	do_fio(&c__1, " Gas Optical Property Files:", (ftnlen)28);
	e_wsfe();
    }
    s_wsfe(&io___61);
    do_fio(&c__1, "  ngt   ng iuabc gasfile", (ftnlen)24);
    e_wsfe();
    i__1 = *ngases;
    for (ng = 1; ng <= i__1; ++ng) {
	i__2 = ngastyp[ng];
	for (ngt = 1; ngt <= i__2; ++ngt) {
	    ++ne;

	    charsp_(gasfile + (ngt + ng * 3) * 132, name__, &len, &c__132, &
		    nlb, &ntb, (ftnlen)132, (ftnlen)1);

	    if (igtrn[ngt + ng * 3] == 1) {

		io___68.ciunit = iuabc[ngt + ng * 3];
		s_rsue(&io___68);
		do_uio(&c__1, (char *)&npg[ngt + ng * 3 - 4], (ftnlen)sizeof(
			integer));
		do_uio(&c__1, (char *)&ntg[ngt + ng * 3 - 4], (ftnlen)sizeof(
			integer));
		do_uio(&c__1, (char *)&igtyp[ngt + ng * 3 - 4], (ftnlen)
			sizeof(integer));
		do_uio(&c__1, (char *)&nbgas[ngt + ng * 3 - 4], (ftnlen)
			sizeof(integer));
		e_rsue();

		s_wsfe(&io___73);
		do_fio(&c__1, (char *)&ngt, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ng, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iuabc[ngt + ng * 3], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&npg[ngt + ng * 3 - 4], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&ntg[ngt + ng * 3 - 4], (ftnlen)sizeof(
			integer));
		i__3 = len;
		for (ii = 1; ii <= i__3; ++ii) {
		    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
		}
		e_wsfe();
	    } else {

		s_wsfe(&io___75);
		do_fio(&c__1, (char *)&ngt, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ng, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iuabc[ngt + ng * 3], (ftnlen)sizeof(
			integer));
		i__3 = len;
		for (ii = 1; ii <= i__3; ++ii) {
		    do_fio(&c__1, name__ + (ii - 1), (ftnlen)1);
		}
		e_wsfe();

	    }

/* L2401: */
	}
/* L2421: */
    }

/* ****   initialze spectral optical property buffers used to determine */
/*       if a given spectral point contributes to the spectrum when */
/*       compared to its surrounding points */

    ij0 = 1;
    iwnflg = 1;
    wni[ij0 - 1] = wn;

    for (ij = 1; ij <= 3; ++ij) {
	wni[ij - 1] = 0.;
	i__1 = nlay;
	for (k = 1; k <= i__1; ++k) {
	    dtauexi[k + ((ij << 1) + 1) * 70 - 211] = small_tau__;
	    dtauexi[k + ((ij << 1) + 2) * 70 - 211] = small_tau__;
	    dtausci[k + ij * 70 - 71] = 0.f;
	    co_pi0i__[k + ((ij << 1) + 1) * 70 - 211] = 0.f;
	    co_pi0i__[k + ((ij << 1) + 2) * 70 - 211] = 0.f;
	    gi[k + ij * 70 - 71] = 0.f;
	    for (mom = 0; mom <= 200; ++mom) {
		phmomi[mom + (k + ij * 70) * 201 - 14271] = 0.f;
/* L2801: */
	    }
	    for (n = 1; n <= 10; ++n) {
		tau_exti__[k + ((n + ij * 10) * 3 + 1) * 70 - 2381] = 0.f;
		tau_exti__[k + ((n + ij * 10) * 3 + 2) * 70 - 2381] = 0.f;
		tau_exti__[k + ((n + ij * 10) * 3 + 3) * 70 - 2381] = 0.f;
		tau_scai__[k + (n + ij * 10) * 70 - 771] = 0.f;
		g_scai__[k + (n + ij * 10) * 70 - 771] = 0.f;
/* L2811: */
	    }
/* L2821: */
	}
	for (l = 1; l <= 4; ++l) {
	    surf_opti__[l + (ij << 2) - 5] = 0.f;
/* L2841: */
	}
	i__1 = *nza;
	for (nz = 1; nz <= i__1; ++nz) {
	    albi[nz + ij * 10 - 11] = 0.f;
/* L2851: */
	}

	if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {

	    p_ray_0i__[ij - 1] = 0.f;
	    p_gas_0i__[ij - 1] = 0.f;
	    p_aer_0i__[ij - 1] = 0.f;
	    tau_rayi__[ij - 1] = 0.f;
	    tau_gasi__[ij - 1] = 0.f;
	    tau_aeri__[ij - 1] = 0.f;
	    i__1 = *numu;
	    for (nze = 1; nze <= i__1; ++nze) {
		p_rayi__[nze + (ij << 4) - 17] = 0.f;
		p_gasi__[nze + (ij << 4) - 17] = 0.f;
		p_aeri__[nze + (ij << 4) - 17] = 0.f;
/* L2861: */
	    }
	}
/* L2881: */
    }

/*       write(*,'(/,1a,/,1a)') ' Progress Report: ', */
/*     -                       '     v(cm**-1)    # intervals' */

/* ****        s p e c t r a l    m a p p i n g    l o o p           ***** */

L3001:
    nspect = 0;

/* ****         Initialize counter for number of grey and */
/*             number of multiple-scattering calculations */

    nscat = 0;
    ngrey = 0;

/* ****   initialize spectral binning parameters */

    init_bin__(nlyr, nref, nza, &ngroup, &iwngrp, &ismterr, &itau1cnt, &
	    itaucnt, &ipi0cnt, &igcnt, &isurcnt, levtau1, nmomgrp, dnugrp, 
	    wngrp, surfgrp, albgrp, taugrp, pi0grp, ggrp, pmomgrp, dtau_b__, 
	    copi0_b__, g_b__, pmom_b__, sur_b__, alb_b__, dx_b_i__, 
	    dalb_b_i__);

/* ****           s t a r t    s p e c t r a l    l o o p            ***** */

L3302:
    ++nspect;

/* ****      set extinction counter */

    ne = 0;

/* ****       initialize extinction and scattering optical depths */
/*           and moments of the scattering phase function at wn */
/*           (note: values go to nlay = nlyr+1 to accommodate pressure */
/*            partial derivatives). */

    i__1 = nlay;
    for (k = 1; k <= i__1; ++k) {
	dtauex[k - 1] = small_tau__;
	dtauex[k + 69] = small_tau__;
	dtausc[k - 1] = 0.f;
	for (n = 1; n <= 10; ++n) {
	    tau_ext__[k + (n * 3 + 1) * 70 - 281] = 0.f;
	    tau_ext__[k + (n * 3 + 2) * 70 - 281] = 0.f;
	    tau_ext__[k + (n * 3 + 3) * 70 - 281] = 0.f;
	    tau_sca__[k + n * 70 - 71] = 0.f;
	    g_sca__[k + n * 70 - 71] = 0.f;
/* L3401: */
	}
	g[k - 1] = 0.f;
	for (mom = 0; mom <= 200; ++mom) {
	    phmom[mom + k * 201 - 201] = 0.f;
/* L3421: */
	}
/* L3441: */
    }
    for (l = 1; l <= 4; ++l) {
	surf_opt__[l - 1] = 0.f;
/* L3461: */
    }
    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	alb[nz - 1] = 0.f;
/* L3481: */
    }

/* ****        initialize the counter for variable optical depth */
/*            components of the state vector */

    ntau_pd0__ = 0;

/* ****       initialize the number of phase function moments to 2, */
/*           to account for Raleigh scattering.  More modes can be */
/*           added for aerosols in aer_tau. */

    nmom = *nstr;

/* ****            r a y l e i g h    s c a t t e r i n g            ***** */


    ray_tau__(&ne, nlyr, ncomp, &icomp[1], &volmix[1], &p[1], &grav[1], 
	    wgtatm, &nzup, nstr, iflext, nsiext, &istate[1], &pd_frac__[1], 
	    wnext, &wn, &umu[1], dtauex, dtausc, phmom, tauray, p_ray__, &
	    p_ray_0__);

    tau_ray__ = tauray[*nlev - 1];

/* ****            g a s    a b s o r p t i o n                      ***** */


    gas_tau__(&iuabc[4], &igtrn[4], &igpab[4], ngases, &ngastyp[1], npg, ntg, 
	    nbgas, &ne, nlev, &nzup, &nt_pd__, &ntaug, nstr, iflext, nsiext, &
	    nlay, &ntau_pd0__, &istate[1], &wn, wnmin, wnmax, ratm, &umu[1], &
	    p[1], &t[1], &alt[1], &grav[1], &z__[1], &rmix[71], &pd_frac__[1],
	     wnext, dtauex, tau_ext__, taugas, p_gas__, &p_gas_0__, &wn_eof__[
	    3], &io_end__[1], &io_err__[1]);


    tau_gas__ = taugas[*nlev - 1];

/* ****       a e r o s o l    o p t i c a l    p r o p e r t i e s */


    aer_tau__(&ne, nlyr, iuaer, nmodes, &nzup, nstr, &nmom, &wn, wnmin, wnmax,
	     &umu[1], &p[1], dtauex, g, phmom, tau_ext__, tau_sca__, g_sca__, 
	    wnext, aerqext, aerqsca, aerg0, aerpmom, dqextdv, dqscadv, 
	    dpmomdv, dg0dv, dtausc, &dtauaer[11], p_aer__, &p_aer_0__, tauaer,
	     iflext, nsiext, nmomaer, nmom_mx__, &ntau_pd0__, &istate[1], 
	    ngases, modepd, &wn_eof__[3], &io_end__[1], &io_err__[1]);


    tau_aer__ = tauaer[*nlev - 1];

/* ****       normalize single scattering albedos and phase functions */


    norm_pi0_ph__(&nmom, &nlay, &nt_pd__, &small_tau__, dtauex, dtausc, 
	    co_pi0__, g, phmom);


/* *****      s u r f a c e    o p t i c a l    p r o p e r t i e s */


    albedo_(&ne, iusur, iflext, nsiext, iref, nref, nza, lamber, &umu0[1], &
	    wn, wnmin, wnmax, wnext, surf_0__, dsurfdv, ws, phiw, surf_opt__, 
	    alb, &wn_eof__[3], &io_end__[1], &io_err__[1]);


/* *****     find the next wavenumber where monochromatic properties */
/*          are specified.  Check all input data sets and select the */
/*          wavenumber that is closest to the current wavenumber. */

/*          note: the minimum distance must exceed the round-off */
/*                tolerance of the computer, wn_tol */

    distmin = 1e30;
    dist0 = *wn_tol__ * wn;

    *next = ne;
    i__1 = *next;
    for (n = 1; n <= i__1; ++n) {
	dist[n - 1] = wnext[nsiext[(n << 1) - 1] + (n << 1) - 3] - wn;
	if (dist[n - 1] < distmin) {
	    if (dist[n - 1] >= dist0) {
		distmin = dist[n - 1];
	    } else {
		distmin = dist0;
	    }
	}
/* L4001: */
    }

/* ****       set the next spectral point */

    if (distmin > distmax) {
	distmin = distmax;
    }
    wn1 = wn + distmin;

/* ****       for iwnmode = 1, determine if the optical properties are */
/*           changing rapidly enough to justify an RT calculation for */
/*           this point */

    iwnmode = 1;


    skip_wn__(&iwnmode, &ij0, &nlay, &iwnflg, ntau_pd__, numu, irad, nza, &
	    dwnmx, &small_tau__, &small_pi0__, &small_g__, &small_alb__, 
	    tauerr, pi0err, phferr, surferr, &wn, dtauex, dtausc, co_pi0__, g,
	     phmom, surf_opt__, alb, wni, dtauexi, dtausci, co_pi0i__, gi, 
	    phmomi, surf_opti__, albi, tau_ext__, tau_sca__, g_sca__, 
	    tau_exti__, tau_scai__, g_scai__, &p_ray_0__, p_ray__, &tau_ray__,
	     p_ray_0i__, p_rayi__, tau_rayi__, &p_gas_0__, p_gas__, &
	    tau_gas__, p_gas_0i__, p_gasi__, tau_gasi__, &p_aer_0__, p_aer__, 
	    &tau_aer__, p_aer_0i__, p_aeri__, tau_aeri__);

/* *****     check to see all three spectral intervals ar loaded */

    if (ij0 < 3) {
	++ij0;
	goto L3302;
    }

/* ****       check to see if this spectral point is needed (iwnflg = 1) */

/* ****        the following flag should shut down skip_wn */

/*           iwnflg =1 */

    if (iwnflg == 1) {

/* *****       All 3 points are good. put x(1) into working arrays */
/*            and process that point, and then pack x(2) into x(1) */
/*            and x(3) into x(2) */

	iwnmode = 2;

/*            write(*,*) 'above 2rd call skip_wn' */
/*            call flush(0) */

	skip_wn__(&iwnmode, &ij0, &nlay, &iwnflg, ntau_pd__, numu, irad, nza, 
		&dwnmx, &small_tau__, &small_pi0__, &small_g__, &small_alb__, 
		tauerr, pi0err, phferr, surferr, &wn, dtauex, dtausc, 
		co_pi0__, g, phmom, surf_opt__, alb, wni, dtauexi, dtausci, 
		co_pi0i__, gi, phmomi, surf_opti__, albi, tau_ext__, 
		tau_sca__, g_sca__, tau_exti__, tau_scai__, g_scai__, &
		p_ray_0__, p_ray__, &tau_ray__, p_ray_0i__, p_rayi__, 
		tau_rayi__, &p_gas_0__, p_gas__, &tau_gas__, p_gas_0i__, 
		p_gasi__, tau_gasi__, &p_aer_0__, p_aer__, &tau_aer__, 
		p_aer_0i__, p_aeri__, tau_aeri__);


/* ****         set the next wavenumber and the effective spectral */
/*             interval width, delnu.  The spectral interval is assumed */
/*             to extend between the current wavenumber, wn, and half-way */
/*             to the previous value, wn0, and the next value, wni(2). */

	delnu = (real) (wni[0] - wn0) * .5f;
	wn0 = wn;

/* ****         s p e c t r a l   b i n n i n g  s e c t i o n */

/* ****       define the column-integrated optical depth and */
/*           determine if scattering is important in this interval */

	ipi0 = 0;
	tauabs = 0.f;
	tausca = 0.f;
	tau_tot__ = 0.f;
	i__1 = nlay;
	for (k = 1; k <= i__1; ++k) {
	    tau_tot__ += dtauex[k - 1];
	    pi0 = 1.f - co_pi0__[k - 1];
	    tauabs += co_pi0__[k - 1] * dtauex[k - 1];
	    dtausca = pi0 * dtauex[k - 1];
	    tausca += dtausca;
	    if (dtausca > small_tau__ && pi0 > ss_min__) {
		ipi0 = 1;
	    }
/* L4201: */
	}

/* ****           determine whether a full multiple scattering */
/*               calcluation must be performed (map_spect) or whether */
/*               a simple extinction calculation will work. */

/* ****           NOTE: setting ipi0 = 1 turns off the grey calculation */

	if (ipi0 == 1 || ! (*lamber)) {

/* ****         c r e a t e    s p e c t r a l    m a p */

	    map_spect__(nlyr, &nlay, nstr, &nmom, levout, levtau1, nref, nza, 
		    &wn, &delnu, dtauex, co_pi0__, g, phmom, surf_opt__, alb, 
		    &taumn, tauerr, pi0err, phferr, surferr, dnugrp, wngrp, 
		    surfgrp, albgrp, taugrp, pi0grp, ggrp, pmomgrp, nmomgrp, &
		    iwngrp, &ngroup, &ismterr, &itau1cnt, &itaucnt, &ipi0cnt, 
		    &igcnt, &isurcnt);

	} else {

/* ****            scattering is not important.  Set the bin number */
/*                to zero and perform a single-scattering calculation */

	    iwngrp = 0;
	    ++ngrey;
	    ismterr = 0;

	}

/* ****          write out the group number at each each wavenumber */

	if (ismterr == 0) {
	    wnsmax = wn;

/* ****            write group parameters for this wn */


	    io___178.ciunit = *iugrp;
	    s_wsue(&io___178);
	    do_uio(&c__1, (char *)&wn, (ftnlen)sizeof(doublereal));
	    do_uio(&c__1, (char *)&iwngrp, (ftnlen)sizeof(integer));
	    for (ii = 1; ii <= 4; ++ii) {
		do_uio(&c__1, (char *)&surf_opt__[ii - 1], (ftnlen)sizeof(
			real));
	    }
	    i__1 = *nza;
	    for (nz = 1; nz <= i__1; ++nz) {
		do_uio(&c__1, (char *)&alb[nz - 1], (ftnlen)sizeof(real));
	    }
	    i__2 = nt_pd__;
	    for (l = 1; l <= i__2; ++l) {
		i__3 = nlay;
		for (k = 1; k <= i__3; ++k) {
		    do_uio(&c__1, (char *)&dtauex[k + l * 70 - 71], (ftnlen)
			    sizeof(real));
		}
	    }
	    i__4 = nt_pd__;
	    for (l = 1; l <= i__4; ++l) {
		i__5 = nlay;
		for (k = 1; k <= i__5; ++k) {
		    do_uio(&c__1, (char *)&co_pi0__[k + l * 70 - 71], (ftnlen)
			    sizeof(real));
		}
	    }
	    i__6 = *nlyr;
	    for (k = 1; k <= i__6; ++k) {
		do_uio(&c__1, (char *)&g[k - 1], (ftnlen)sizeof(real));
	    }
	    i__7 = *ntau_pd__;
	    for (n = 1; n <= i__7; ++n) {
		for (m = 1; m <= 3; ++m) {
		    i__8 = *nlev;
		    for (k = 1; k <= i__8; ++k) {
			do_uio(&c__1, (char *)&tau_ext__[k + (m + n * 3) * 70 
				- 281], (ftnlen)sizeof(real));
		    }
		}
	    }
	    i__9 = *ntau_pd__;
	    for (n = 1; n <= i__9; ++n) {
		i__10 = *nlev;
		for (k = 1; k <= i__10; ++k) {
		    do_uio(&c__1, (char *)&tau_sca__[k + n * 70 - 71], (
			    ftnlen)sizeof(real));
		}
	    }
	    i__11 = *ntau_pd__;
	    for (n = 1; n <= i__11; ++n) {
		i__12 = *nlev;
		for (k = 1; k <= i__12; ++k) {
		    do_uio(&c__1, (char *)&g_sca__[k + n * 70 - 71], (ftnlen)
			    sizeof(real));
		}
	    }
	    e_wsue();

	    if (*irad == 2 || *irad == 4 || *irad == 6 || *irad == 8) {
		io___180.ciunit = *iugrp;
		s_wsue(&io___180);
		do_uio(&c__1, (char *)&p_ray_0__, (ftnlen)sizeof(real));
		i__1 = *numu;
		for (nze = nzup; nze <= i__1; ++nze) {
		    do_uio(&c__1, (char *)&p_ray__[nze - 1], (ftnlen)sizeof(
			    real));
		}
		do_uio(&c__1, (char *)&tau_ray__, (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&p_gas_0__, (ftnlen)sizeof(real));
		i__3 = *numu;
		for (nze = nzup; nze <= i__3; ++nze) {
		    do_uio(&c__1, (char *)&p_gas__[nze - 1], (ftnlen)sizeof(
			    real));
		}
		do_uio(&c__1, (char *)&tau_gas__, (ftnlen)sizeof(real));
		do_uio(&c__1, (char *)&p_aer_0__, (ftnlen)sizeof(real));
		i__2 = *numu;
		for (nze = nzup; nze <= i__2; ++nze) {
		    do_uio(&c__1, (char *)&p_aer__[nze - 1], (ftnlen)sizeof(
			    real));
		}
		do_uio(&c__1, (char *)&tau_aer__, (ftnlen)sizeof(real));
		e_wsue();
	    }


	}

/* ****         update wavenumber counter */

	wntot += 1.f;

/* ****        determine if this wavenumber is beyond end of desired */
/*            spectral range.  if so, jump out of wavenumber loop. */

	if (ismterr != 0 || wn >= *wnmax) {
	    goto L4601;
	}

    } else {

/* ****       skip point x(2).  Load x(3) into x(2) */

	wnskip += 1.;
	iwnmode = 3;

	skip_wn__(&iwnmode, &ij0, &nlay, &iwnflg, ntau_pd__, numu, irad, nza, 
		&dwnmx, &small_tau__, &small_pi0__, &small_g__, &small_alb__, 
		tauerr, pi0err, phferr, surferr, &wn, dtauex, dtausc, 
		co_pi0__, g, phmom, surf_opt__, alb, wni, dtauexi, dtausci, 
		co_pi0i__, gi, phmomi, surf_opti__, albi, tau_ext__, 
		tau_sca__, g_sca__, tau_exti__, tau_scai__, g_scai__, &
		p_ray_0__, p_ray__, &tau_ray__, p_ray_0i__, p_rayi__, 
		tau_rayi__, &p_gas_0__, p_gas__, &tau_gas__, p_gas_0i__, 
		p_gasi__, tau_gasi__, &p_aer_0__, p_aer__, &tau_aer__, 
		p_aer_0i__, p_aeri__, tau_aeri__);

    }

/* *****        specify the next wavenumber grid point and reset */
/*             wavenumber grid flags. */

    wn = wn1;

/* ****         set extinction flags for each constituent */

    i__1 = *next;
    for (n = 1; n <= i__1; ++n) {

/* ****             find the distance between the last location where */
/*                 the extintion from this agent was specified and the */
/*                 current wavelength. */

	dist[n - 1] = wnext[nsiext[(n << 1) - 1] + (n << 1) - 3] - wn;
	if (dist[n - 1] <= 0.) {

/* ****              wn is beyond latest input point for this constituent. */
/*                  set flag to read next point. */

	    iflext[n - 1] = 1;
	} else {
	    iflext[n - 1] = 0;
	}
/* L4401: */
    }

    if (wn > wn_step__) {

/* ****         update spectral grid monitor */

	wn_step__ += dwn_step__;

    }

/* ****        bin monochromatic optical properties for next wavenumber. */

    goto L3302;

/* ****          Record time required for binning */

L4601:
    cpu_time__(&time1);
    smcputime = time1 - time0;
    ismtime = i_nint(&smcputime);
    tcpu += smcputime;

/* ****        Update group counters */

    ngrtot += ngroup;

    i__1 = ngroup;
    for (n = 1; n <= i__1; ++n) {

/* ****      load optical properties into binned radiance arrays */

	ng0 = n;

	load_optics__(lamber, &ng0, nlyr, nza, nstate, &istate[1], n_rad__, 
		nref, iref, nmomgrp, &umu0[1], tauerr, pi0err, phferr, 
		surferr, &taumn, surfgrp, albgrp, taugrp, pi0grp, ggrp, 
		pmomgrp, alb_b__, sur_b__, phiw, ws, dtau_b__, copi0_b__, 
		g_b__, pmom_b__, dx_b_i__, dalb_b_i__);

	usrang0 = *usrang;
	numu0 = *numu;

	for (l = 1; l <= 5; ++l) {

/* ****               find the layer transmittances, reflectances and */
/*                   absorptances for each profile */

	    nlyr0 = *nlyr;
	    l0 = l;
	    nmom = nmomgrp[n - 1];
	    numu0 = *numu;

	    layer_trn__(&usrang0, &l0, &ng0, &nlyr0, nstr, &numu0, &nmom, 
		    nphi, n_rad__, dtau_b__, copi0_b__, pmom_b__, dx_b_i__, &
		    umu[1], &phi[1], umu_f__, gwt_f__, trnflx_b__, dtrnflxdx, 
		    refflx_b__, drefflxdx, absflx_b__, dabsflxdx, trnrad_b__, 
		    dtrnraddx, refrad_b__, drefraddx, absrad_b__, dabsraddx);

/* L4801: */
	}

/* L4881: */
    }

/* ****        s o l a r     z e n i t h    a n g l e    l o o p */

    rtime = 0.f;
    btime = 0.f;
    i__1 = *nza;
    for (nz = 1; nz <= i__1; ++nz) {
	nz0 = nz;
	iu_flux__ = *iuflx + nz - 1;
	umu0nz = umu0[nz];
	phi0nz = phi0[nz];

/* ****             r a d i a n c e s   f o r    e a c h   b i n */

	i__3 = ngroup;
	for (n = 1; n <= i__3; ++n) {
	    ng0 = n;

	    sm_eq_trn__(usrang, lamber, lplanck, lsolar, &ng0, &nz0, nlyr, 
		    nstr, numu, &nzdn, &nzup, nphi, nmomgrp, nstate, &istate[
		    1], n_rad__, iref, nref, &nscat, tauerr, pi0err, phferr, 
		    surferr, ws, phiw, &taumn, &umu[1], umu_f__, gwt_f__, &
		    phi[1], &umu0nz, &phi0nz, ts, &t[1], accur, wngrp, 
		    surfgrp, taugrp, pi0grp, ggrp, pmomgrp, alb_b__, sur_b__, 
		    dtau_b__, copi0_b__, g_b__, pmom_b__, dx_b_i__, 
		    dalb_b_i__, trnflx_b__, refflx_b__, absflx_b__, dtrnflxdx,
		     drefflxdx, dabsflxdx, trnrad_b__, refrad_b__, absrad_b__,
		     brdf_b__, dtrnraddx, drefraddx, dabsraddx, dbrdfdx, 
		    dnsflxsrc_b__, upsflxsrc_b__, sradsrc_b__, ddnsflxdx, 
		    dupsflxdx, dsraddx, dntflxsrc_b__, uptflxsrc_b__, 
		    tradsrc_b__, ddntflxdx, duptflxdx, dtraddx);

/* L5221: */
	}

/* ****            Record time needed for radiance calculations */

	cpu_time__(&time2);
	rcputime = time2 - time1;
	rtime += rcputime;
	tcpu += rcputime;

/* ****             Rewind unit with group index data for backmapping */

	al__1.aerr = 0;
	al__1.aunit = *iugrp;
	f_rew(&al__1);

/* ****             s m t    b a c k    m a p p i n g */

	map_back__(lsolar, lplanck, lamber, usrang, iugrp, iuthrm, iusol0, 
		iusol2, iutrn, &iu_flux__, iuout, &npd, ipd, iunits, irad, 
		iu_pd__, &iutpd[11], &iuspd[11], &iupdrad[10513], ifrmout, 
		nza, &nz0, nstr, numu, &nzup, &nzdn, nphi, nlyr, &nlay, 
		levout, nlout, &k_out__[1], &igs[1], iref, modepd, &nt_pd__, 
		nstate, &istate[1], ntau_pd__, isptype, islit, width, dwn, 
		points, wnmin, wnmax, &wnsmin, &wnsmax, wn_tol__, &units, &
		umu0nz, &phi0nz, umu_f__, gwt_f__, &umu[1], &phi[1], 
		pd_pert__, &p[1], &t[1], ts, &dp_dp__[1], &rmix[71], &dtauaer[
		11], alb_b__, dtau_b__, copi0_b__, g_b__, trnflx_b__, 
		dtrnflxdx, refflx_b__, drefflxdx, absflx_b__, dabsflxdx, 
		refrad_b__, drefraddx, absrad_b__, dabsraddx, brdf_b__, 
		dbrdfdx, dnsflxsrc_b__, ddnsflxdx, upsflxsrc_b__, dupsflxdx, 
		dntflxsrc_b__, ddntflxdx, uptflxsrc_b__, duptflxdx, 
		sradsrc_b__, dsraddx, tradsrc_b__, dtraddx, &dirsoflx[71], &
		dnsoflx[71], &upsoflx[71], &dnthflx[1], &upthflx[1]);

/* ****             Record time needed for back-mapping */

	cpu_time__(&time3);
	bcputime = time3 - time2;
	btime += bcputime;
	tcpu += bcputime;

/* L5301: */
    }

/* ****         print spectral grouping statistics */

    if (wnsmin == *wnmin) {
	io___227.ciunit = *iustat;
	s_wsfe(&io___227);
	do_fio(&c__1, "   Spectral Binning Statistics: ", (ftnlen)32);
	do_fio(&c__1, "   wnmin      wnmax   # bins  total  # mono   bin", (
		ftnlen)49);
	do_fio(&c__1, "    bin rejection criteria", (ftnlen)26);
	do_fio(&c__1, "                       calc. type", (ftnlen)33);
	do_fio(&c__1, "        timing(seconds)     total", (ftnlen)33);
	do_fio(&c__1, "   cm**-1     cm**-1  in int  bins segments  ratio", (
		ftnlen)50);
	do_fio(&c__1, "   tau=1    tau <>1     pi0        g   albedo", (
		ftnlen)45);
	do_fio(&c__1, "    #scat   #grey     bin     rad bk_map       cpu", (
		ftnlen)50);
	e_wsfe();

    }

    irtime = i_nint(&rtime);
    ibtime = i_nint(&btime);
    icpu = i_nint(&tcpu);

    ieff = 1;
    if (ngroup > 0) {
	ieff = nspect / ngroup;
    }

    io___232.ciunit = *iustat;
    s_wsfe(&io___232);
    do_fio(&c__1, (char *)&wnsmin, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&wnsmax, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ngroup, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ngrtot, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nspect, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ieff, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&itau1cnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&itaucnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ipi0cnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&igcnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&isurcnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nscat, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ngrey, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ismtime, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&irtime, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ibtime, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&icpu, (ftnlen)sizeof(integer));
    e_wsfe();

/* ****        find the solar and thermal heating rates */

    heating_(iuheat, nlyr, nza, umu0_1__, wnmin, &wnsmax, &p[1], &t[1], &alt[
	    1], &aid_lr__[1], &dirsoflx[71], &dnsoflx[71], &upsoflx[71], &
	    dnthflx[1], &upthflx[1], &soheat[71], &thheat[1]);

    if (wn < *wnmax) {

/* ****         rewind units for spectral bins and thermal fluxes */
/*             and empty all output buffers the spectral output units. */

	al__1.aerr = 0;
	al__1.aunit = *iugrp;
	f_rew(&al__1);
	al__1.aerr = 0;
	al__1.aunit = *iuthrm;
	f_rew(&al__1);

/* ****           finish remainder of spectrum */

	wnsmin = wn1;
	wn = wn1;

	goto L3001;

    } else {

	cl__1.cerr = 0;
	cl__1.cunit = *iugrp;
	cl__1.csta = "delete";
	f_clos(&cl__1);
	cl__1.cerr = 0;
	cl__1.cunit = *iuthrm;
	cl__1.csta = "delete";
	f_clos(&cl__1);

    }

/* ****          b i n n i n g   e f f i c i e n c y */

    io___233.ciunit = *iustat;
    s_wsfe(&io___233);
    do_fio(&c__1, " Total number of spectral intervals skipped =   ", (ftnlen)
	    48);
    do_fio(&c__1, (char *)&wnskip, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___234.ciunit = *iustat;
    s_wsfe(&io___234);
    do_fio(&c__1, " Total number of monochromatic intervals used = ", (ftnlen)
	    48);
    do_fio(&c__1, (char *)&wntot, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___235.ciunit = *iustat;
    s_wsfe(&io___235);
    do_fio(&c__1, " Total number of spectral mapping bins =        ", (ftnlen)
	    48);
    do_fio(&c__1, (char *)&ngrtot, (ftnlen)sizeof(integer));
    e_wsfe();
    eff = 1.f;
    if (ngrtot > 0) {
	eff = (real) wntot / (real) ngrtot;
    }
    io___237.ciunit = *iustat;
    s_wsfe(&io___237);
    do_fio(&c__1, " Net efficiency =", (ftnlen)17);
    do_fio(&c__1, (char *)&eff, (ftnlen)sizeof(real));
    e_wsfe();
    io___238.ciunit = *iustat;
    s_wsfe(&io___238);
    do_fio(&c__1, " Total CPU time =", (ftnlen)17);
    do_fio(&c__1, (char *)&tcpu, (ftnlen)sizeof(real));
    do_fio(&c__1, " seconds", (ftnlen)8);
    e_wsfe();

/* ********        E n d    o f    s p e c t r a l    l o o p    ********* */

    return 0;
} /* smart_ */

