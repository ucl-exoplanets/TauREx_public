/* DISORT.f -- translated by f2c (version 20100827).
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

static integer c__4 = 4;
static integer c__71 = 71;
static integer c__70 = 70;
static integer c__16 = 16;
static integer c__1000 = 1000;
static integer c__17 = 17;
static integer c__256 = 256;
static integer c__1190 = 1190;
static integer c__272 = 272;
static integer c__1120 = 1120;
static integer c__17920 = 17920;
static integer c__8 = 8;
static integer c__64 = 64;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__25 = 25;
static integer c__9 = 9;
static integer c__2 = 2;
static integer c__3 = 3;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* rcs version control information: */
/*  disort.f,v 2.1 2000/04/04 18:21:55 laszlo exp $ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Subroutine */ int disort_(integer *nlyr, real *dtauc, real *ssalb, integer 
	*nmom, real *pmom, real *temper, real *wvnm, logical *usrtau, integer 
	*ntau, real *utau, integer *nstr, integer *iunits, logical *usrang, 
	integer *numu, real *umu, integer *nphi, real *phi, integer *ibcnd, 
	real *fbeam, real *umu0, real *phi0, real *fisot, logical *lamber, 
	integer *iref, real *surf_pr__, real *albedo, real *btemp, real *
	ttemp, real *temis, logical *plank, logical *onlyfl, real *accur, 
	logical *prnt, integer *maxcly, integer *maxulv, integer *maxumu, 
	integer *maxphi, integer *maxmom, real *rfldir, real *rfldn, real *
	flup, real *dfdt, real *uavg, real *uu, real *albmed, real *trnmed)
{
    /* Initialized data */

    static logical pass1 = TRUE_;
    static logical prntu0[2] = { FALSE_,FALSE_ };

    /* System generated locals */
    integer pmom_dim1, pmom_offset, uu_dim1, uu_dim2, uu_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    double asin(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    static real b[1120];
    static integer j, l;
    static real z__[1120], z0[16], z1[16];
    static doublereal bb[70];
    static real cc[256]	/* was [16][16] */, gc[17920]	/* was [16][16][70] */
	    ;
    static integer lc;
    static real gl[1190]	/* was [17][70] */, kk[1120]	/* was [16][
	    70] */, ll[1120]	/* was [16][70] */, pi;
    static integer iq, nn;
    static real gu[17920]	/* was [16][16][70] */;
    static integer iu, lu, ns;
    static real wk[16], zj[16], zz[1120]	/* was [16][70] */, u0c[1120]	
	    /* was [16][70] */, u0u[1120]	/* was [16][70] */, xr0[70], 
	    xr1[70], z0u[1120]	/* was [16][70] */, z1u[1120]	/* was [16][
	    70] */;
    static doublereal aad[64]	/* was [8][8] */;
    static real amb[64]	/* was [8][8] */, apb[64]	/* was [8][8] */, bem[
	    8], bdr[72]	/* was [8][9] */, cmu[16], dum;
    static integer lev;
    static real rpd;
    static integer naz;
    static real sgn, emu[16];
    static doublereal wkd[16];
    static real cwt[16], rmu[144]	/* was [16][9] */, sqt[1000], uum[
	    1120]	/* was [16][70] */, psi0[16], psi1[16], ylm0[17], 
	    pkag[71], fldn[70], eval[8];
    static integer ncol;
    static real tauc[71];
    static integer ncos;
    static real ylmc[272]	/* was [17][16] */;
    static integer ncut;
    static real flyr[70];
    static integer ipvt[1120];
    static real ylmu[272]	/* was [17][16] */, delm0, zplk0[1120]	/* 
	    was [16][70] */, zplk1[1120]	/* was [16][70] */, cband[
	    78400]	/* was [70][1120] */, evecc[256]	/* was [16][
	    16] */;
    static doublereal evald[8];
    static real phasa[70], zbeam[1120]	/* was [16][70] */, fldir[70], phasm[
	    70];
    static integer mazim;
    static real array[256]	/* was [16][16] */, phast[70];
    static integer kconv;
    extern doublereal ratio_(real *, real *);
    static real azerr, oprim[70];
    static integer layru[70];
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int solve0_(real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, integer *, logical *, 
	    real *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static doublereal eveccd[64]	/* was [8][8] */;
    extern /* Subroutine */ int chekin_(integer *, real *, real *, integer *, 
	    real *, real *, real *, logical *, integer *, integer *, real *, 
	    integer *, logical *, integer *, real *, integer *, real *, 
	    integer *, real *, real *, real *, real *, logical *, real *, 
	    real *, real *, real *, real *, logical *, logical *, logical *, 
	    logical *, real *, real *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static real expbea[71];
    static logical deltam;
    static real bplank, phirad[16], tempbb[70], angcos[1];
    extern /* Subroutine */ int upbeam_(real *, real *, real *, real *, real *
	    , real *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *), planck_(
	    integer *, integer *, integer *, integer *, real *, real *, 
	    doublereal *);
    static real dither, dtaucp[70];
    static logical compar;
    extern /* Subroutine */ int albtrn_(real *, real *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *, real *, real *, integer *, 
	    integer *, integer *, integer *, logical *, real *, real *, real *
	    , real *, real *, real *, real *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, real *, real *, real *), soleig_(
	    real *, real *, real *, real *, real *, real *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static real cosphi;
    extern /* Subroutine */ int surfac_(real *, real *, real *, real *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, real *, real *, real *, logical *, real *, 
	    real *, real *, real *, real *);
    static real tplank;
    extern /* Subroutine */ int cmpint_(real *, real *, real *, integer *, 
	    real *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, logical *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static logical corint;
    static real taucpr[71], azterm;
    extern /* Subroutine */ int intcor_(real *, real *, real *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *), fluxes_(real *, 
	    real *, real *, real *, real *, integer *, real *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, real *, logical *, logical *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *), lepoly_(
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    ), pravin_(real *, integer *, integer *, real *, integer *, real *
	    ), setdis_(real *, real *, logical *, real *, real *, real *, 
	    real *, real *, real *, integer *, integer *, logical *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, logical *, integer *, logical *, logical *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    logical *, logical *), terpev_(real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *), prtinp_(integer *, real *, real *, real *
	    , integer *, real *, real *, real *, integer *, real *, integer *,
	     integer *, real *, integer *, real *, integer *, real *, real *, 
	    real *, real *, logical *, real *, real *, real *, real *, 
	    logical *, logical *, logical *, logical *, real *, real *, 
	    logical *, real *, real *, real *, integer *, logical *);
    static real utaupr[70];
    extern /* Subroutine */ int prtint_(real *, real *, integer *, real *, 
	    integer *, real *, integer *, integer *, integer *);
    static logical lyrcut;
    extern /* Subroutine */ int setmtx_(real *, real *, real *, real *, real *
	    , real *, real *, real *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, real *, real *), slftst_(logical *, real *, real *, 
	    real *, logical *, real *, real *, real *, integer *, logical *, 
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    logical *, real *, real *, integer *, real *, logical *, logical *
	    , real *, real *, real *, real *, real *, logical *, logical *, 
	    real *, real *, real *, logical *, real *, real *, real *, real *)
	    , terpso_(real *, real *, real *, real *, integer *, integer *, 
	    logical *, integer *, integer *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *), upisot_(real *, real *, real *, integer *
	    , integer *, integer *, integer *, real *, real *, real *, real *,
	     real *, real *, real *, real *), usrint_(real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, logical *, integer *, real *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, integer *, integer *, real *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *), zeroal_(integer *, real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, integer 
	    *, real *, real *, real *, real *, real *, real *, real *, real *,
	     integer *, real *, integer *, real *, real *, real *, integer *, 
	    real *, integer *, real *, integer *, real *, integer *, real *, 
	    real *, real *, real *, real *, integer *, real *, integer *, 
	    integer *, real *, integer *, real *, integer *, real *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *, 
	    integer *, real *, integer *, real *, real *, real *, real *, 
	    real *, integer *, real *, real *, integer *, real *, integer *, 
	    real *), zeroit_(real *, integer *);

/* ******************************************************************* */
/*       plane-parallel discrete ordinates radiative transfer program */
/*             ( see disort.doc for complete documentation ) */
/* ******************************************************************* */

/* ***** this version passes the surface reflection properties */
/*      through to bdref */

/* +------------------------------------------------------------------+ */
/*  calling tree (omitting calls to errmsg): */
/*  (routines in parentheses are not in this file) */

/*  disort-+-(r1mach) */
/*         +-slftst-+-(tstbad) */
/*         +-zeroit */
/*         +-chekin-+-(wrtbad) */
/*         |        +-(wrtdim) */
/*         |        +-dref */
/*         +-zeroal */
/*         +-setdis-+-qgausn-+-(d1mach) */
/*         +-prtinp */
/*         +-albtrn-+-lepoly */
/*         |        +-zeroit */
/*         |        +-soleig-+-asymtx-+-(d1mach) */
/*         |        +-terpev */
/*         |        +-setmtx-+-zeroit */
/*         |        +-(sgbco) */
/*         |        +-solve1-+-zeroit */
/*         |        |        +-(sgbsl) */
/*         |        +-altrin */
/*         |        +-spaltr */
/*         |        +-praltr */
/*         +-planck-+-(r1mach) */
/*         +-lepoly */
/*         +-surfac-+-qgausn-+-(d1mach) */
/*         |        +-bdref */
/*         |        +-zeroit */
/*         +-soleig-+-asymtx-+-(d1mach) */
/*         +-upbeam-+-(sgeco) */
/*         |        +-(sgesl) */
/*         +-upisot-+-(sgeco) */
/*         |        +-(sgesl) */
/*         +-terpev */
/*         +-terpso */
/*         +-setmtx-+-zeroit */
/*         +-solve0-+-zeroit */
/*         |        +-(sgbco) */
/*         |        +-(sgbsl) */
/*         +-fluxes--zeroit */
/*         +-zeroit */
/*         +-usrint */
/*         +-cmpint */
/*         +-pravin */
/*         +-zeroit */
/*         +-ratio--(r1mach) */
/*         +-intcor-+-sinsca */
/*         |        +-secsca-+-xifunc */
/*         +-prtint */

/* *** intrinsic functions used in disort package which take */
/*     non-negligible amount of time: */

/*    exp :  called by- albtrn, altrin, cmpint, fluxes, setdis, */
/*                      setmtx, spaltr, usrint, planck */

/*    sqrt : called by- asymtx, soleig */

/* +-------------------------------------------------------------------+ */

/*  index conventions (for all do-loops and all variable descriptions): */

/*     iu     :  for user polar angles */

/*  iq,jq,kq  :  for computational polar angles ('quadrature angles') */

/*   iq/2     :  for half the computational polar angles (just the ones */
/*               in either 0-90 degrees, or 90-180 degrees) */

/*     j      :  for user azimuthal angles */

/*     k,l    :  for legendre expansion coefficients or, alternatively, */
/*               subscripts of associated legendre polynomials */

/*     lu     :  for user levels */

/*     lc     :  for computational layers (each having a different */
/*               single-scatter albedo and/or phase function) */

/*    lev     :  for computational levels */

/*    mazim   :  for azimuthal components in fourier cosine expansion */
/*               of intensity and phase function */
/* +-------------------------------------------------------------------+ */

/*   variables added to the subroutine call list */

/*     iref   :  surface reflectance mode */
/*             0 - Lambert */
/*             1 - Hapke's BDR model */
/*             2 - Breon's BDR model; combination of Li + Roujean */
/*             3 - Roujean's BDR model */
/*             4 - Cox and Munk glint model */
/*   SURF_PR : Wavelength dependent surface properties array */
/*             IREF= 0 - Lambert albedo */
/*             IREF= 1 - Hapke : HH, W */
/*             IREF= 2 - Breon's BDR model: k0, k1, k2 */
/*             IREF= 3 - Roujean's BDR model: k0, k1, k2 */
/*             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw */

/* +------------------------------------------------------------------+ */

/*               i n t e r n a l    v a r i a b l e s */

/*   amb(iq/2,iq/2)    first matrix factor in reduced eigenvalue problem */
/*                     of eqs. ss(12), stwj(8e), stwl(23f) */
/*                     (used only in soleig) */

/*   apb(iq/2,iq/2)    second matrix factor in reduced eigenvalue problem */
/*                     of eqs. ss(12), stwj(8e), stwl(23f) */
/*                     (used only in soleig) */

/*   array(iq,iq)      scratch matrix for soleig, upbeam and upisot */
/*                     (see each subroutine for definition) */

/*   b()               right-hand side vector of eq. sc(5) going into */
/*                     solve0,1;  returns as solution vector */
/*                     vector  l, the constants of integration */

/*   bdr(iq/2,0:iq/2)  bottom-boundary bidirectional reflectivity for a */
/*                     given azimuthal component.  first index always */
/*                     refers to a computational angle.  second index: */
/*                     if zero, refers to incident beam angle umu0; */
/*                     if non-zero, refers to a computational angle. */

/*   bem(iq/2)         bottom-boundary directional emissivity at compu- */
/*                     tational angles. */

/*   bplank            intensity emitted from bottom boundary */

/*   cband()           matrix of left-hand side of the linear system */
/*                     eq. sc(5), scaled by eq. sc(12);  in banded */
/*                     form required by linpack solution routines */

/*   cc(iq,iq)         c-sub-ij in eq. ss(5) */

/*   cmu(iq)           computational polar angles (gaussian) */

/*   cwt(iq)           quadrature weights corresponding to cmu */

/*   corint            when set true, correct intensities for */
/*                     delta-scaling effects (see nakajima and tanaka, */
/*                     1988). when false, intensities are not corrected. */
/*                     in general, corint should be set true when beam */
/*                     source is present (fbeam is not zero) asurf_prnd deltam */
/*                     is true in a problem including scattering. */
/*                     however, execution is faster when corint is false, */
/*                     and intensities outside the aureole may still be */
/*                     accurate enough.  when corint is true, it is */
/*                     important to have a sufficiently high order of */
/*                     legendre approximation of the phase function. this */
/*                     is because the intensities are corrected by */
/*                     calculating the single-scattered radiation, for */
/*                     which an adequate representation of the phase */
/*                     function is crucial.  in case of a low order */
/*                     legendre approximation of an otherwise highly */
/*                     anisotropic phase function, the intensities might */
/*                     actually be more accurate when corint is false. */
/*                     when only fluxes are calculated (onlyfl is true), */
/*                     or there is no beam source (fbeam=0.0),surf_pr or there */
/*                     is no scattering (ssalb=0.0 for all layers) corint */
/*                     is set false by the code. */

/*   delm0             kronecker delta, delta-sub-m0, where m = mazim */
/*                     is the number of the fourier component in the */
/*                     azimuth cosine expansion */

/*   deltam            true,  use delta-m method ( see wiscombe, 1977 ); */
/*                     false, do not use delta-m method. in general, for */
/*                     a given number of streams, intensities and */
/*                     fluxes will be more accurate for phase functions */
/*                     with a large forward peak if deltam is set true. */
/*                     intensities close to the forward scattering */
/*                     direction are often less accurate, however, when */
/*                     the delta-m method is applied. the intensity */
/*                     correction of nakajima and tanaka is used to */
/*                     improve the accuracy of the intensities. */

/*   dither            small quantity subtracted from single-scattering */
/*                     albedos of unity, in order to avoid usurf_prsing special */
/*                     case formulas;  prevents an eigenvalue of exactly */
/*                     zero from occurring, which would cause an */
/*                     immediate overflow */

/*   dtaucp(lc)        computational-layer optical depths (delta-m-scaled */
/*                     if deltam = true, otherwise equal to dtauc) */

/*   emu(iu)           bottom-boundary directional emissivity at user */
/*                     angles. */

/*   eval(iq)          temporary storage for eigenvalues of eq. ss(12) */

/*   evecc(iq,iq)      complete eigenvectors of ss(7) on return from */
/*                     soleig; stored permanently in  gc */

/*   expbea(lc)        transmission of direct beam in delta-m optical */
/*                     depth coordinates */

/*   flyr(lc)          separated fraction in delta-m method */
/* surf_pr */
/*   gl(k,lc)          phase function legendre polynomial expansion */
/*                     coefficients, calculated from pmom by */
/*                     including single-scattering albedo, factor */
/*                     2k+1, and (if deltam=true) the delta-m */
/*                     scaling */

/*   gc(iq,iq,lc)      eigenvectors at polar quadrature angles, */
/*                     g  in eq. sc(1) */

/*   gu(iu,iq,lc)      eigenvectors interpolated to user polar angles */
/*                     ( g  in eqs. sc(3) and s1(8-9), i.e. */
/*                       g without the l factor ) */

/*   ipvt(lc*iq)       integer vector of pivot indices for linpack */
/*                     routines */
/* surf_pr */
/*   kk(iq,lc)         eigenvalues of coeff. matrix in eq. ss(7) */

/*   kconv             counter in azimuth convergence test */

/*   layru(lu)         computational layer in which user output level */
/*                     utau(lu) is located */

/*   ll(iq,lc)         constants of integration l in eq. sc(1), */
/*                     obtained by solving scaled version of eq. sc(5) */

/*   lyrcut            true, radiation is assumed zero below layer */
/*                     ncut because of almost complete absorption */

/*   naz               number of azimuthal components considered */

/*   ncut              computational layer number in which absorption */
/*                     optical depth first exceeds abscut */

/*   oprim(lc)         single scattering albedo after delta-m scaling */

/*   pass1             true on first entry, false thereafter */

/*   pkag(0:lc)        integrated planck function for internal emission */

/*   prntu0(l)         logical flag to trigger printing of azimuthally- */
/*                     averaged intensities: */
/*                       l    quantities printed */
/*                      --    ------------------ */
/*                       1    azimuthally-averaged intensities at user */
/*                               levels and computational polar angles */
/*                       2    azimuthally-averaged intensities at user */
/*                               levels and user polar angles */

/*   psi0(iq)          sum just after square bracket in  eq. sd(9) */

/*   psi1(iq)          sum in  eq. stwl(31d) */

/*   rmu(iu,0:iq)      bottom-boundary bidirectional reflectivity for a */
/*                     given azimuthal component.  first index always */
/*                     refers to a user angle.  second index: */
/*                     if zero, refers to incident beam angle umu0; */
/*                     if non-zero, refers to a computational angle. */

/*   sqt(k)            square root of k (used only in lepoly for */
/*                     computing associated legendre polynomials) */

/*   tauc(0:lc)        cumulative optical depth (un-delta-m-scaled) */

/*   taucpr(0:lc)      cumulative optical depth (delta-m-scaled if */
/*                     deltam = true, otherwise equal to tauc) */

/*   tplank            intensity emitted from top boundary */

/*   uum(iu,lu)        expansion coefficients when the intensity */
/*                     (u-super-m) is expanded in fourier cosine series */
/*                     in azimuth angle */

/*   u0c(iq,lu)        azimuthally-averaged intensity at quadrature */
/*                     angle */

/*   u0u(iu,lu)        if onlyfl = false, azimuthally-averaged intensity */
/*                     at user angles and user levels */

/*                     if onlyfl = true and maxumu.ge.nstr, */
/*                     azimuthally-averaged intensity at computational */
/*                     (gaussian quadrature) angles and user levels; */
/*                     the corresponding quadrature angle cosines are */
/*                     returned in umu.  if maxumu.lt.nstr, u0u will be */
/*                     zeroed, and umu, numu will not be set. */

/*   utaupr(lu)        optical depths of user output levels in delta-m */
/*                     coordinates;  equal to  utau(lu) if no delta-m */

/*   wk()              scratch array */

/*   xr0(lc)           x-sub-zero in expansion of thermal source func- */
/*                     tion preceding eq. ss(14)(has no mu-dependence); */
/*                     b-sub-zero in eq. stwl(24d) */

/*   xr1(lc)           x-sub-one in expansion of thermal source func- */
/*                     tion; see  eqs. ss(14-16); b-sub-one in stwl(24d) */

/*   ylm0(l)           normalized associated legendre polynomial */
/*                     of subscript l at the beam angle (not saved */
/*                     as function of superscipt m) */

/*   ylmc(l,iq)        normalized associated legendre polynomial */
/*                     of subscript l at the computational angles */
/*                     (not saved as function of superscipt m) */

/*   ylmu(l,iu)        normalized associated legendre polynomial */
/*                     of subscript l at the user angles */
/*                     (not saved as function of superscipt m) */

/*   z()               scratch array used in solve0, albtrn to solve */
/*                     a linear system for the constants of integration */

/*   z0(iq)            solution vectors z-sub-zero of eq. ss(16) */

/*   z0u(iu,lc)        z-sub-zero in eq. ss(16) interpolated to user */
/*                     angles from an equation derived from ss(16) */

/*   z1(iq)            solution vectors z-sub-one  of eq. ss(16) */

/*   z1u(iu,lc)        z-sub-one in eq. ss(16) interpolated to user */
/*                     angles from an equation derived from ss(16) */

/*   zbeam(iu,lc)      particular solution for beam source */

/*   zj(iq)            right-hand side vector  x-sub-zero in */
/*                     eq. ss(19), also the solution vector */
/*                     z-sub-zero after solving that system */

/*   zz(iq,lc)         permanent storage for the beam source vectors zj */

/*   zplk0(iq,lc)      permanent storage for the thermal source */
/*                     vectors  z0  obtained by solving  eq. ss(16) */

/*   zplk1(iq,lc)      permanent storage for the thermal source */
/*                     vectors  z1  obtained by solving  eq. ss(16) */

/* +-------------------------------------------------------------------+ */

/*  local symbolic dimensions (have big effect on storage requirements): */

/*       mxcly  = max no. of computational layers */
/*       mxulv  = max no. of output levels */
/*       mxcmu  = max no. of computation polar angles */
/*       mxumu  = max no. of output polar angles */
/*       mxphi  = max no. of output azimuthal angles */
/*       mxsqt  = max no. of square roots of integers (for lepoly) */
/* +-------------------------------------------------------------------+ */
/*     .. parameters .. */
/*     .. */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. local arrays .. */

/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --surf_pr__;
    --prnt;
    --ssalb;
    --dtauc;
    --uavg;
    --dfdt;
    --flup;
    --rfldn;
    --rfldir;
    --utau;
    --trnmed;
    --albmed;
    --umu;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2);
    uu -= uu_offset;
    --phi;
    pmom_dim1 = *maxmom - 0 + 1;
    pmom_offset = 0 + pmom_dim1;
    pmom -= pmom_offset;

    /* Function Body */
    deltam = TRUE_;
    corint = TRUE_;
    if (pass1) {
	pi = asin(1.f) * 2.f;
	dither = r1mach_(&c__4) * 10.f;
/*                            ** must dither more on high (14-digit) */
/*                            ** precision machine */
	if (dither < 1e-10f) {
	    dither *= 10.f;
	}
	rpd = pi / 180.f;
	for (ns = 1; ns <= 1000; ++ns) {
	    sqt[ns - 1] = sqrt((real) ns);
/* L10: */
	}
/*                            ** set input values for self-test. */
/*                            ** be sure slftst sets all print flags off. */
	compar = FALSE_;
	slftst_(&corint, accur, albedo, btemp, &deltam, &dtauc[1], fbeam, 
		fisot, ibcnd, lamber, nlyr, plank, nphi, numu, nstr, ntau, 
		onlyfl, &phi[1], phi0, nmom, &pmom[pmom_dim1], &prnt[1], 
		prntu0, &ssalb[1], temis, temper, ttemp, &umu[1], usrang, 
		usrtau, &utau[1], umu0, wvnm, &compar, &dum, &dum, &dum, &dum)
		;
    }
L20:

/*                                  ** calculate cumulative optical depth */
/*                                  ** and dither single-scatter albedo */
/*                                  ** to improve numerical behavior of */
/*                                  ** eigenvalue/vector computation */
    zeroit_(tauc, &c__71);
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	if (ssalb[lc] == 1.f) {
	    ssalb[lc] = 1.f - dither;
	}
	tauc[lc] = tauc[lc - 1] + dtauc[lc];
/* L30: */
    }
/*                                ** check input dimensions and variables */
    chekin_(nlyr, &dtauc[1], &ssalb[1], nmom, &pmom[pmom_offset], temper, 
	    wvnm, usrtau, ntau, iref, &utau[1], nstr, usrang, numu, &umu[1], 
	    nphi, &phi[1], ibcnd, fbeam, umu0, phi0, fisot, lamber, albedo, 
	    btemp, ttemp, temis, &surf_pr__[1], plank, onlyfl, &deltam, &
	    corint, accur, tauc, maxcly, maxulv, maxumu, maxphi, maxmom, &
	    c__70, &c__70, &c__16, &c__16, &c__16, &c__1000);
/*                                 ** zero internal and output arrays */
    i__1 = *maxumu * *maxulv * *maxphi;
    zeroal_(&c__70, &expbea[1], flyr, oprim, phasa, phast, phasm, &taucpr[1], 
	    xr0, xr1, &c__16, cmu, cwt, psi0, psi1, wk, z0, z1, zj, &c__17, 
	    ylm0, &c__256, array, cc, evecc, &c__1190, gl, &c__272, ylmc, &
	    c__272, ylmu, &c__1120, kk, ll, zz, zplk0, zplk1, &c__17920, gc, &
	    c__70, layru, utaupr, &c__17920, gu, &c__1120, z0u, z1u, zbeam, &
	    c__8, eval, &c__64, amb, apb, &c__1120, ipvt, z__, maxulv, &
	    rfldir[1], &rfldn[1], &flup[1], &uavg[1], &dfdt[1], maxumu, &
	    albmed[1], &trnmed[1], &c__1120, u0u, &i__1, &uu[uu_offset]);
/*                                 ** perform various setup operations */
    setdis_(cmu, cwt, &deltam, &dtauc[1], dtaucp, expbea, fbeam, flyr, gl, 
	    ibcnd, layru, &lyrcut, maxmom, maxumu, &c__16, &ncut, nlyr, ntau, 
	    &nn, nstr, plank, numu, onlyfl, &corint, oprim, &pmom[pmom_offset]
	    , &ssalb[1], tauc, taucpr, &utau[1], utaupr, &umu[1], umu0, 
	    usrtau, usrang);
/*                                 ** print input information */
    if (prnt[1]) {
	prtinp_(nlyr, &dtauc[1], dtaucp, &ssalb[1], nmom, &pmom[pmom_offset], 
		temper, wvnm, ntau, &utau[1], nstr, numu, &umu[1], nphi, &phi[
		1], ibcnd, fbeam, umu0, phi0, fisot, lamber, albedo, btemp, 
		ttemp, temis, &deltam, plank, onlyfl, &corint, accur, flyr, &
		lyrcut, oprim, tauc, taucpr, maxmom, &prnt[5]);
    }
/*                              ** handle special case for getting albedo */
/*                              ** and transmissivity of medium for many */
/*                              ** beam angles at once */
    if (*ibcnd == 1) {
	albtrn_(albedo, amb, apb, array, b, bdr, cband, cc, cmu, cwt, dtaucp, 
		eval, evecc, gl, gc, gu, ipvt, kk, ll, nlyr, &nn, nstr, numu, 
		&prnt[1], taucpr, &umu[1], u0u, wk, ylmc, ylmu, z__, aad, 
		evald, eveccd, wkd, &c__8, &c__70, maxumu, &c__16, &c__16, &
		c__1120, sqt, &albmed[1], &trnmed[1]);
	return 0;
    }
/*                                   ** calculate planck functions */
    if (! (*plank)) {
	bplank = 0.f;
	tplank = 0.f;
	zeroit_(pkag, &c__71);
    } else {
	tempbb[0] = *ttemp;
	planck_(&c__1, &c__1, &c__1, iunits, wvnm, tempbb, bb);
	tplank = (real) (*temis * bb[0]);
	tempbb[0] = *btemp;
	planck_(&c__1, &c__1, &c__1, &c__1, wvnm, tempbb, bb);
	bplank = (real) bb[0];
	i__1 = *nlyr;
	for (lev = 0; lev <= i__1; ++lev) {
	    tempbb[lev] = temper[lev];
/* L39: */
	}

	i__1 = *nlyr + 1;
	planck_(&c__1, &i__1, &c__1, iunits, wvnm, tempbb, bb);

	i__1 = *nlyr;
	for (lev = 0; lev <= i__1; ++lev) {
	    pkag[lev] = (real) bb[lev];
/* L40: */
	}
    }

/* ========  begin loop to sum azimuthal components of intensity  ======= */
/*           (eq stwj 5, stwl 6) */
    kconv = 0;
    naz = *nstr - 1;
/*                                    ** azimuth-independent case */
    if (*fbeam == 0.f || (r__1 = 1.f - *umu0, dabs(r__1)) < 1e-5f || *onlyfl 
	    || *numu == 1 && (r__2 = 1.f - umu[1], dabs(r__2)) < 1e-5f || *
	    numu == 1 && (r__3 = umu[1] + 1.f, dabs(r__3)) < 1e-5f || *numu ==
	     2 && (r__4 = umu[1] + 1.f, dabs(r__4)) < 1e-5f && (r__5 = 1.f - 
	    umu[2], dabs(r__5)) < 1e-5f) {
	naz = 0;
    }
    i__1 = naz;
    for (mazim = 0; mazim <= i__1; ++mazim) {
	if (mazim == 0) {
	    delm0 = 1.f;
	}
	if (mazim > 0) {
	    delm0 = 0.f;
	}
/*                             ** get normalized associated legendre */
/*                             ** polynomials for */
/*                             ** (a) incident beam angle cosine */
/*                             ** (b) computational and user polar angle */
/*                             **     cosines */
	if (*fbeam > 0.f) {
	    ncos = 1;
	    angcos[0] = -(*umu0);
	    i__2 = *nstr - 1;
	    lepoly_(&ncos, &mazim, &c__16, &i__2, angcos, sqt, ylm0);
	}
	if (! (*onlyfl) && *usrang) {
	    i__2 = *nstr - 1;
	    lepoly_(numu, &mazim, &c__16, &i__2, &umu[1], sqt, ylmu);
	}
	i__2 = *nstr - 1;
	lepoly_(&nn, &mazim, &c__16, &i__2, cmu, sqt, ylmc);
/*                       ** get normalized associated legendre polys. */
/*                       ** with negative arguments from those with */
/*                       ** positive arguments; dave/armstrong eq. (15), */
/*                       ** stwl(59) */
	sgn = -1.f;
	i__2 = *nstr - 1;
	for (l = mazim; l <= i__2; ++l) {
	    sgn = -sgn;
	    i__3 = *nstr;
	    for (iq = nn + 1; iq <= i__3; ++iq) {
		ylmc[l + iq * 17 - 17] = sgn * ylmc[l + (iq - nn) * 17 - 17];
/* L60: */
	    }
/* L70: */
	}
/*                                 ** specify users bottom reflectivity */
/*                                 ** and emissivity properties */
	if (! lyrcut) {
	    surfac_(albedo, &delm0, cmu, fbeam, lamber, iref, &c__8, &mazim, &
		    c__16, &nn, numu, onlyfl, &pi, &umu[1], umu0, usrang, &
		    surf_pr__[1], bdr, emu, bem, rmu);
	}
/* ===================  begin loop on computational layers  ============= */
	i__2 = ncut;
	for (lc = 1; lc <= i__2; ++lc) {
/*                      ** solve eigenfunction problem in eq. stwj(8b), */
/*                      ** stwl(23f); return eigenvalues and eigenvectors */
	    soleig_(amb, apb, array, cmu, cwt, &gl[lc * 17 - 17], &c__8, &
		    mazim, &c__16, &nn, nstr, ylmc, cc, evecc, eval, &kk[(lc 
		    << 4) - 16], &gc[((lc << 4) + 1 << 4) - 272], aad, eveccd,
		     evald, wkd);
/*                                  ** calculate particular solutions of */
/*                                  ** eq. ss(18), stwl(24a) for incident */
/*                                  ** beam source */
	    if (*fbeam > 0.f) {
		upbeam_(array, cc, cmu, &delm0, fbeam, &gl[lc * 17 - 17], 
			ipvt, &mazim, &c__16, &nn, nstr, &pi, umu0, wk, ylm0, 
			ylmc, zj, &zz[(lc << 4) - 16]);
	    }
/*                              ** calculate particular solutions of eq. */
/*                              ** ss(15), stwl(25) for thermal emission */
/*                              ** source */

	    if (*plank && mazim == 0) {
		xr1[lc - 1] = 0.f;
		if (dtaucp[lc - 1] > 0.f) {
		    xr1[lc - 1] = (pkag[lc] - pkag[lc - 1]) / dtaucp[lc - 1];
		}
		xr0[lc - 1] = pkag[lc - 1] - xr1[lc - 1] * taucpr[lc - 1];
		upisot_(array, cc, cmu, ipvt, &c__16, &nn, nstr, &oprim[lc - 
			1], wk, &xr0[lc - 1], &xr1[lc - 1], z0, z1, &zplk0[(
			lc << 4) - 16], &zplk1[(lc << 4) - 16]);
	    }
	    if (! (*onlyfl) && *usrang) {
/*                                            ** interpolate eigenvectors */
/*                                            ** to user angles */
		terpev_(cwt, evecc, &gl[lc * 17 - 17], &gu[((lc << 4) + 1 << 
			4) - 272], &mazim, &c__16, &c__16, &nn, nstr, numu, 
			wk, ylmc, ylmu);
/*                                            ** interpolate source terms */
/*                                            ** to user angles */
		terpso_(cwt, &delm0, fbeam, &gl[lc * 17 - 17], &mazim, &c__16,
			 plank, numu, nstr, &oprim[lc - 1], &pi, ylm0, ylmc, 
			ylmu, psi0, psi1, &xr0[lc - 1], &xr1[lc - 1], z0, z1, 
			zj, &zbeam[(lc << 4) - 16], &z0u[(lc << 4) - 16], &
			z1u[(lc << 4) - 16]);
	    }
/* L80: */
	}
/* ===================  end loop on computational layers  =============== */
/*                      ** set coefficient matrix of equations combining */
/*                      ** boundary and layer interface conditions */
	setmtx_(bdr, cband, cmu, cwt, &delm0, dtaucp, gc, kk, lamber, &lyrcut,
		 &c__8, &c__70, &c__16, &ncol, &ncut, &c__1120, &nn, nstr, 
		taucpr, wk);
/*                      ** solve for constants of integration in homo- */
/*                      ** geneous solution (general boundary conditions) */
	solve0_(b, bdr, bem, &bplank, cband, cmu, cwt, expbea, fbeam, fisot, 
		ipvt, lamber, ll, &lyrcut, &mazim, &c__8, &c__70, &c__16, &
		ncol, &ncut, &nn, nstr, &c__1120, &pi, &tplank, taucpr, umu0, 
		z__, zz, zplk0, zplk1);
/*                                  ** compute upward and downward fluxes */
	if (mazim == 0) {
	    fluxes_(cmu, cwt, fbeam, gc, kk, layru, ll, &lyrcut, maxulv, &
		    c__16, &c__70, &ncut, &nn, nstr, ntau, &pi, &prnt[1], 
		    prntu0, &ssalb[1], taucpr, umu0, &utau[1], utaupr, xr0, 
		    xr1, zz, zplk0, zplk1, &dfdt[1], &flup[1], fldn, fldir, &
		    rfldir[1], &rfldn[1], &uavg[1], u0c);
	}
	if (*onlyfl) {
	    if (*maxumu >= *nstr) {
/*                                     ** save azimuthal-avg intensities */
/*                                     ** at quadrature angles */
		i__2 = *ntau;
		for (lu = 1; lu <= i__2; ++lu) {
		    i__3 = *nstr;
		    for (iq = 1; iq <= i__3; ++iq) {
			u0u[iq + (lu << 4) - 17] = u0c[iq + (lu << 4) - 17];
/* L90: */
		    }
/* L100: */
		}
	    }
	    goto L190;
	}
	zeroit_(uum, &c__1120);
	if (*usrang) {
/*                                     ** compute azimuthal intensity */
/*                                     ** components at user angles */
	    usrint_(&bplank, cmu, cwt, &delm0, dtaucp, emu, expbea, fbeam, 
		    fisot, gc, gu, kk, lamber, layru, ll, &lyrcut, &mazim, &
		    c__16, &c__70, &c__16, &ncut, nlyr, &nn, nstr, plank, 
		    numu, ntau, &pi, rmu, taucpr, &tplank, &umu[1], umu0, 
		    utaupr, wk, zbeam, z0u, z1u, zz, zplk0, zplk1, uum);
	} else {
/*                                     ** compute azimuthal intensity */
/*                                     ** components at quadrature angles */
	    cmpint_(fbeam, gc, kk, layru, ll, &lyrcut, &mazim, &c__16, &c__70,
		     &c__16, &ncut, &nn, nstr, plank, ntau, taucpr, umu0, 
		    utaupr, zz, zplk0, zplk1, uum);
	}
	if (mazim == 0) {
/*                               ** save azimuthally averaged intensities */
	    i__2 = *ntau;
	    for (lu = 1; lu <= i__2; ++lu) {
		i__3 = *numu;
		for (iu = 1; iu <= i__3; ++iu) {
		    u0u[iu + (lu << 4) - 17] = uum[iu + (lu << 4) - 17];
		    i__4 = *nphi;
		    for (j = 1; j <= i__4; ++j) {
			uu[iu + (lu + j * uu_dim2) * uu_dim1] = uum[iu + (lu 
				<< 4) - 17];
/* L110: */
		    }
/* L120: */
		}
/* L130: */
	    }
/*                              ** print azimuthally averaged intensities */
/*                              ** at user angles */
	    if (prntu0[1]) {
		pravin_(&umu[1], numu, &c__16, &utau[1], ntau, u0u);
	    }
	    if (naz > 0) {
		zeroit_(phirad, &c__16);
		i__2 = *nphi;
		for (j = 1; j <= i__2; ++j) {
		    phirad[j - 1] = rpd * (phi[j] - *phi0);
/* L140: */
		}
	    }
	} else {
/*                                ** increment intensity by current */
/*                                ** azimuthal component (fourier */
/*                                ** cosine series);  eq sd(2), stwl(6) */
	    azerr = 0.f;
	    i__2 = *nphi;
	    for (j = 1; j <= i__2; ++j) {
		cosphi = cos(mazim * phirad[j - 1]);
		i__3 = *ntau;
		for (lu = 1; lu <= i__3; ++lu) {
		    i__4 = *numu;
		    for (iu = 1; iu <= i__4; ++iu) {
			azterm = uum[iu + (lu << 4) - 17] * cosphi;
			uu[iu + (lu + j * uu_dim2) * uu_dim1] += azterm;
/* Computing MAX */
			r__4 = dabs(azterm);
			r__5 = (r__1 = uu[iu + (lu + j * uu_dim2) * uu_dim1], 
				dabs(r__1));
			r__2 = azerr, r__3 = ratio_(&r__4, &r__5);
			azerr = dmax(r__2,r__3);
/* L150: */
		    }
/* L160: */
		}
/* L170: */
	    }
	    if (azerr <= *accur) {
		++kconv;
	    }
	    if (kconv >= 2) {
		goto L190;
	    }
	}
/* L180: */
    }
/* ===================  end loop on azimuthal components  =============== */
L190:
/*                                    ** apply nakajima/tanaka intensity */
/*                                    ** corrections */
    if (corint) {
	intcor_(&dither, fbeam, flyr, layru, &lyrcut, maxmom, maxulv, maxumu, 
		nmom, &ncut, nphi, nstr, ntau, numu, oprim, phasa, phast, 
		phasm, phirad, &pi, &rpd, &pmom[pmom_offset], &ssalb[1], &
		dtauc[1], tauc, taucpr, &umu[1], umu0, &utau[1], utaupr, &uu[
		uu_offset]);
    }
/*                                          ** print intensities */
    if (prnt[3] && ! (*onlyfl)) {
	prtint_(&uu[uu_offset], &utau[1], ntau, &umu[1], numu, &phi[1], nphi, 
		maxulv, maxumu);
    }
    if (pass1) {
/*                                    ** compare test case results with */
/*                                    ** correct answers and abort if bad */
	compar = TRUE_;
	slftst_(&corint, accur, albedo, btemp, &deltam, &dtauc[1], fbeam, 
		fisot, ibcnd, lamber, nlyr, plank, nphi, numu, nstr, ntau, 
		onlyfl, &phi[1], phi0, nmom, &pmom[pmom_dim1], &prnt[1], 
		prntu0, &ssalb[1], temis, temper, ttemp, &umu[1], usrang, 
		usrtau, &utau[1], umu0, wvnm, &compar, &flup[1], &rfldir[1], &
		rfldn[1], &uu[(uu_dim2 + 1) * uu_dim1 + 1]);
	pass1 = FALSE_;
	goto L20;
    }
    return 0;
} /* disort_ */

/* Subroutine */ int asymtx_(real *aa, real *evec, real *eval, integer *m, 
	integer *ia, integer *ievec, integer *ier, doublereal *wkd, 
	doublereal *aad, doublereal *evecd, doublereal *evald)
{
    /* Initialized data */

    static doublereal c4 = .95;
    static doublereal c5 = 16.;
    static doublereal c6 = 256.;
    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal c1 = .4375;
    static doublereal c2 = .5;
    static doublereal c3 = .75;

    /* System generated locals */
    integer aa_dim1, aa_offset, evec_dim1, evec_offset, aad_dim1, aad_offset, 
	    evecd_dim1, evecd_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__, j, k, l, n;
    static doublereal p, q, r__, s, t, w, x, y, z__;
    static integer n1, n2, ka, lb, ii, in;
    static logical tf;
    static doublereal uu, vv, col;
    static integer kkk, lll;
    static doublereal sgn, tol, row, repl, scale, rnorm;
    extern doublereal d1mach_(integer *);
    static doublereal discri;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static logical notlas, noconv;

/*    =======  d o u b l e    p r e c i s i o n    v e r s i o n  ====== */
/*       solves eigenfunction problem for real asymmetric matrix */
/*       for which it is known a priori that the eigenvalues are real. */
/*       this is an adaptation of a subroutine eigrf in the imsl */
/*       library to use real instead of complex arithmetic, accounting */
/*       for the known fact that the eigenvalues and eigenvectors in */
/*       the discrete ordinate solution are real.  other changes include */
/*       putting all the called subroutines in-line, deleting the */
/*       performance index calculation, updating many do-loops */
/*       to fortran77, and in calculating the machine precision */
/*       tol instead of specifying it in a data statement. */
/*       eigrf is based primarily on eispack routines.  the matrix is */
/*       first balanced using the parlett-reinsch algorithm.  then */
/*       the martin-wilkinson algorithm is applied. */
/*       there is a statement 'j  = wkd( i )' that converts a double */
/*       precision variable to an integer variable, that seems dangerous */
/*       to us in principle, but seems to work fine in practice. */
/*       references: */
/*          dongarra, j. and c. moler, eispack -- a package for solving */
/*             matrix eigenvalue problems, in cowell, ed., 1984: */
/*             sources and development of mathematical software, */
/*             prentice-hall, englewood cliffs, nj */
/*         parlett and reinsch, 1969: balancing a matrix for calculation */
/*             of eigenvalues and eigenvectors, num. math. 13, 293-304 */
/*         wilkinson, j., 1965: the algebraic eigenvalue problem, */
/*             clarendon press, oxford */
/*   i n p u t    v a r i a b l e s: */

/*       aa    :  input asymmetric matrix, destroyed after solved */

/*        m    :  order of  aa */

/*       ia    :  first dimension of  aa */

/*    ievec    :  first dimension of  evec */


/*   o u t p u t    v a r i a b l e s: */

/*       evec  :  (unnormalized) eigenvectors of  aa */
/*                   ( column j corresponds to eval(j) ) */

/*       eval  :  (unordered) eigenvalues of aa ( dimension at least m ) */

/*       ier   :  if .ne. 0, signals that eval(ier) failed to converge; */
/*                   in that case eigenvalues ier+1,ier+2,...,m  are */
/*                   correct but eigenvalues 1,...,ier are set to zero. */


/*   s c r a t c h   v a r i a b l e s: */

/*       wkd   :  work area ( dimension at least 2*m ) */
/*       aad   :  double precision stand-in for aa */
/*       evecd :  double precision stand-in for evec */
/*       evald :  double precision stand-in for eval */

/*   called by- soleig */
/*   calls- d1mach, errmsg */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --evald;
    --eval;
    evecd_dim1 = *ia;
    evecd_offset = 1 + evecd_dim1;
    evecd -= evecd_offset;
    aad_dim1 = *ia;
    aad_offset = 1 + aad_dim1;
    aad -= aad_offset;
    aa_dim1 = *ia;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    evec_dim1 = *ievec;
    evec_offset = 1 + evec_dim1;
    evec -= evec_offset;
    --wkd;

    /* Function Body */
    *ier = 0;
    tol = d1mach_(&c__4);
    if (*m < 1 || *ia < *m || *ievec < *m) {
	tf = TRUE_;
	errmsg_("asymtx--bad input variable(s)", &tf, (ftnlen)29);
    }
/*                           ** handle 1x1 and 2x2 special cases */
    if (*m == 1) {
	eval[1] = aa[aa_dim1 + 1];
	evec[evec_dim1 + 1] = 1.f;
	return 0;
    } else if (*m == 2) {
/* Computing 2nd power */
	r__1 = aa[aa_dim1 + 1] - aa[(aa_dim1 << 1) + 2];
	discri = r__1 * r__1 + aa[(aa_dim1 << 1) + 1] * 4.f * aa[aa_dim1 + 2];
	tf = TRUE_;
	if (discri < 0.f) {
	    errmsg_("asymtx--complex evals in 2x2 case", &tf, (ftnlen)33);
	}
	sgn = one;
	if (aa[aa_dim1 + 1] < aa[(aa_dim1 << 1) + 2]) {
	    sgn = -one;
	}
	eval[1] = (aa[aa_dim1 + 1] + aa[(aa_dim1 << 1) + 2] + sgn * sqrt(
		discri)) * .5f;
	eval[2] = (aa[aa_dim1 + 1] + aa[(aa_dim1 << 1) + 2] - sgn * sqrt(
		discri)) * .5f;
	evec[evec_dim1 + 1] = 1.f;
	evec[(evec_dim1 << 1) + 2] = 1.f;
	if (aa[aa_dim1 + 1] == aa[(aa_dim1 << 1) + 2] && (aa[aa_dim1 + 2] == 
		0.f || aa[(aa_dim1 << 1) + 1] == 0.f)) {
	    rnorm = (r__1 = aa[aa_dim1 + 1], dabs(r__1)) + (r__2 = aa[(
		    aa_dim1 << 1) + 1], dabs(r__2)) + (r__3 = aa[aa_dim1 + 2],
		     dabs(r__3)) + (r__4 = aa[(aa_dim1 << 1) + 2], dabs(r__4))
		    ;
	    w = tol * rnorm;
	    evec[evec_dim1 + 2] = aa[aa_dim1 + 2] / w;
	    evec[(evec_dim1 << 1) + 1] = -aa[(aa_dim1 << 1) + 1] / w;
	} else {
	    evec[evec_dim1 + 2] = aa[aa_dim1 + 2] / (eval[1] - aa[(aa_dim1 << 
		    1) + 2]);
	    evec[(evec_dim1 << 1) + 1] = aa[(aa_dim1 << 1) + 1] / (eval[2] - 
		    aa[aa_dim1 + 1]);
	}
	return 0;
    }
/*                               ** convert single-prec. matrix to double */
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
	    aad[j + k * aad_dim1] = aa[j + k * aa_dim1];
/* L10: */
	}
/* L20: */
    }
/*                                ** initialize output variables */
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	evald[i__] = zero;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    evecd[i__ + j * evecd_dim1] = zero;
/* L30: */
	}
	evecd[i__ + i__ * evecd_dim1] = one;
/* L40: */
    }
/*                  ** balance the input matrix and reduce its norm by */
/*                  ** diagonal similarity transformation stored in wk; */
/*                  ** then search for rows isolating an eigenvalue */
/*                  ** and push them down */
    rnorm = zero;
    l = 1;
    k = *m;
L50:
    kkk = k;
    for (j = kkk; j >= 1; --j) {
	row = zero;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ != j) {
		row += (d__1 = aad[j + i__ * aad_dim1], abs(d__1));
	    }
/* L60: */
	}
	if (row == zero) {
	    wkd[k] = (doublereal) j;
	    if (j != k) {
		i__1 = k;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    repl = aad[i__ + j * aad_dim1];
		    aad[i__ + j * aad_dim1] = aad[i__ + k * aad_dim1];
		    aad[i__ + k * aad_dim1] = repl;
/* L70: */
		}
		i__1 = *m;
		for (i__ = l; i__ <= i__1; ++i__) {
		    repl = aad[j + i__ * aad_dim1];
		    aad[j + i__ * aad_dim1] = aad[k + i__ * aad_dim1];
		    aad[k + i__ * aad_dim1] = repl;
/* L80: */
		}
	    }
	    --k;
	    goto L50;
	}
/* L90: */
    }
/*                                ** search for columns isolating an */
/*                                ** eigenvalue and push them left */
L100:
    lll = l;
    i__1 = k;
    for (j = lll; j <= i__1; ++j) {
	col = zero;
	i__2 = k;
	for (i__ = l; i__ <= i__2; ++i__) {
	    if (i__ != j) {
		col += (d__1 = aad[i__ + j * aad_dim1], abs(d__1));
	    }
/* L110: */
	}
	if (col == zero) {
	    wkd[l] = (doublereal) j;
	    if (j != l) {
		i__2 = k;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    repl = aad[i__ + j * aad_dim1];
		    aad[i__ + j * aad_dim1] = aad[i__ + l * aad_dim1];
		    aad[i__ + l * aad_dim1] = repl;
/* L120: */
		}
		i__2 = *m;
		for (i__ = l; i__ <= i__2; ++i__) {
		    repl = aad[j + i__ * aad_dim1];
		    aad[j + i__ * aad_dim1] = aad[l + i__ * aad_dim1];
		    aad[l + i__ * aad_dim1] = repl;
/* L130: */
		}
	    }
	    ++l;
	    goto L100;
	}
/* L140: */
    }
/*                           ** balance the submatrix in rows l through k */
    i__1 = k;
    for (i__ = l; i__ <= i__1; ++i__) {
	wkd[i__] = one;
/* L150: */
    }
L160:
    noconv = FALSE_;
    i__1 = k;
    for (i__ = l; i__ <= i__1; ++i__) {
	col = zero;
	row = zero;
	i__2 = k;
	for (j = l; j <= i__2; ++j) {
	    if (j != i__) {
		col += (d__1 = aad[j + i__ * aad_dim1], abs(d__1));
		row += (d__1 = aad[i__ + j * aad_dim1], abs(d__1));
	    }
/* L170: */
	}
	f = one;
	g = row / c5;
	h__ = col + row;
L180:
	if (col < g) {
	    f *= c5;
	    col *= c6;
	    goto L180;
	}
	g = row * c5;
L190:
	if (col >= g) {
	    f /= c5;
	    col /= c6;
	    goto L190;
	}
/*                                                ** now balance */
	if ((col + row) / f < c4 * h__) {
	    wkd[i__] *= f;
	    noconv = TRUE_;
	    i__2 = *m;
	    for (j = l; j <= i__2; ++j) {
		aad[i__ + j * aad_dim1] /= f;
/* L200: */
	    }
	    i__2 = k;
	    for (j = 1; j <= i__2; ++j) {
		aad[j + i__ * aad_dim1] *= f;
/* L210: */
	    }
	}
/* L220: */
    }
    if (noconv) {
	goto L160;
    }
/*                                   ** is a already in hessenberg form? */
    if (k - 1 < l + 1) {
	goto L370;
    }
/*                                   ** transfer a to a hessenberg form */
    i__1 = k - 1;
    for (n = l + 1; n <= i__1; ++n) {
	h__ = zero;
	wkd[n + *m] = zero;
	scale = zero;
/*                                                 ** scale column */
	i__2 = k;
	for (i__ = n; i__ <= i__2; ++i__) {
	    scale += (d__1 = aad[i__ + (n - 1) * aad_dim1], abs(d__1));
/* L230: */
	}
	if (scale != zero) {
	    i__2 = n;
	    for (i__ = k; i__ >= i__2; --i__) {
		wkd[i__ + *m] = aad[i__ + (n - 1) * aad_dim1] / scale;
/* Computing 2nd power */
		d__1 = wkd[i__ + *m];
		h__ += d__1 * d__1;
/* L240: */
	    }
	    d__1 = sqrt(h__);
	    g = -d_sign(&d__1, &wkd[n + *m]);
	    h__ -= wkd[n + *m] * g;
	    wkd[n + *m] -= g;
/*                                            ** form (i-(u*ut)/h)*a */
	    i__2 = *m;
	    for (j = n; j <= i__2; ++j) {
		f = zero;
		i__3 = n;
		for (i__ = k; i__ >= i__3; --i__) {
		    f += wkd[i__ + *m] * aad[i__ + j * aad_dim1];
/* L250: */
		}
		i__3 = k;
		for (i__ = n; i__ <= i__3; ++i__) {
		    aad[i__ + j * aad_dim1] -= wkd[i__ + *m] * f / h__;
/* L260: */
		}
/* L270: */
	    }
/*                                    ** form (i-(u*ut)/h)*a*(i-(u*ut)/h) */
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		f = zero;
		i__3 = n;
		for (j = k; j >= i__3; --j) {
		    f += wkd[j + *m] * aad[i__ + j * aad_dim1];
/* L280: */
		}
		i__3 = k;
		for (j = n; j <= i__3; ++j) {
		    aad[i__ + j * aad_dim1] -= wkd[j + *m] * f / h__;
/* L290: */
		}
/* L300: */
	    }
	    wkd[n + *m] = scale * wkd[n + *m];
	    aad[n + (n - 1) * aad_dim1] = scale * g;
	}
/* L310: */
    }
    i__1 = l;
    for (n = k - 2; n >= i__1; --n) {
	n1 = n + 1;
	n2 = n + 2;
	f = aad[n + 1 + n * aad_dim1];
	if (f != zero) {
	    f *= wkd[n + 1 + *m];
	    i__2 = k;
	    for (i__ = n + 2; i__ <= i__2; ++i__) {
		wkd[i__ + *m] = aad[i__ + n * aad_dim1];
/* L320: */
	    }
	    if (n + 1 <= k) {
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
		    g = zero;
		    i__3 = k;
		    for (i__ = n + 1; i__ <= i__3; ++i__) {
			g += wkd[i__ + *m] * evecd[i__ + j * evecd_dim1];
/* L330: */
		    }
		    g /= f;
		    i__3 = k;
		    for (i__ = n + 1; i__ <= i__3; ++i__) {
			evecd[i__ + j * evecd_dim1] += g * wkd[i__ + *m];
/* L340: */
		    }
/* L350: */
		}
	    }
	}
/* L360: */
    }
L370:
    n = 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = n; j <= i__2; ++j) {
	    rnorm += (d__1 = aad[i__ + j * aad_dim1], abs(d__1));
/* L380: */
	}
	n = i__;
	if (i__ < l || i__ > k) {
	    evald[i__] = aad[i__ + i__ * aad_dim1];
	}
/* L390: */
    }
    n = k;
    t = zero;
/*                                      ** search for next eigenvalues */
L400:
    if (n < l) {
	goto L550;
    }
    in = 0;
    n1 = n - 1;
    n2 = n - 2;
/*                          ** look for single small sub-diagonal element */
L410:
    i__1 = n;
    for (i__ = l; i__ <= i__1; ++i__) {
	lb = n + l - i__;
	if (lb == l) {
	    goto L430;
	}
	s = (d__1 = aad[lb - 1 + (lb - 1) * aad_dim1], abs(d__1)) + (d__2 = 
		aad[lb + lb * aad_dim1], abs(d__2));
	if (s == zero) {
	    s = rnorm;
	}
	if ((d__1 = aad[lb + (lb - 1) * aad_dim1], abs(d__1)) <= tol * s) {
	    goto L430;
	}
/* L420: */
    }
L430:
    x = aad[n + n * aad_dim1];
    if (lb == n) {
/*                                        ** one eigenvalue found */
	aad[n + n * aad_dim1] = x + t;
	evald[n] = aad[n + n * aad_dim1];
	n = n1;
	goto L400;
    }
    y = aad[n1 + n1 * aad_dim1];
    w = aad[n + n1 * aad_dim1] * aad[n1 + n * aad_dim1];
    if (lb == n1) {
/*                                        ** two eigenvalues found */
	p = (y - x) * c2;
/* Computing 2nd power */
	d__1 = p;
	q = d__1 * d__1 + w;
	z__ = sqrt((abs(q)));
	aad[n + n * aad_dim1] = x + t;
	x = aad[n + n * aad_dim1];
	aad[n1 + n1 * aad_dim1] = y + t;
/*                                        ** real pair */
	z__ = p + d_sign(&z__, &p);
	evald[n1] = x + z__;
	evald[n] = evald[n1];
	if (z__ != zero) {
	    evald[n] = x - w / z__;
	}
	x = aad[n + n1 * aad_dim1];
/*                                  ** employ scale factor in case */
/*                                  ** x and z are very small */
	r__ = sqrt(x * x + z__ * z__);
	p = x / r__;
	q = z__ / r__;
/*                                             ** row modification */
	i__1 = *m;
	for (j = n1; j <= i__1; ++j) {
	    z__ = aad[n1 + j * aad_dim1];
	    aad[n1 + j * aad_dim1] = q * z__ + p * aad[n + j * aad_dim1];
	    aad[n + j * aad_dim1] = q * aad[n + j * aad_dim1] - p * z__;
/* L440: */
	}
/*                                             ** column modification */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__ = aad[i__ + n1 * aad_dim1];
	    aad[i__ + n1 * aad_dim1] = q * z__ + p * aad[i__ + n * aad_dim1];
	    aad[i__ + n * aad_dim1] = q * aad[i__ + n * aad_dim1] - p * z__;
/* L450: */
	}
/*                                          ** accumulate transformations */
	i__1 = k;
	for (i__ = l; i__ <= i__1; ++i__) {
	    z__ = evecd[i__ + n1 * evecd_dim1];
	    evecd[i__ + n1 * evecd_dim1] = q * z__ + p * evecd[i__ + n * 
		    evecd_dim1];
	    evecd[i__ + n * evecd_dim1] = q * evecd[i__ + n * evecd_dim1] - p 
		    * z__;
/* L460: */
	}
	n = n2;
	goto L400;
    }
    if (in == 30) {
/*                    ** no convergence after 30 iterations; set error */
/*                    ** indicator to the index of the current eigenvalue */
	*ier = n;
	goto L700;
    }
/*                                                  ** form shift */
    if (in == 10 || in == 20) {
	t += x;
	i__1 = n;
	for (i__ = l; i__ <= i__1; ++i__) {
	    aad[i__ + i__ * aad_dim1] -= x;
/* L470: */
	}
	s = (d__1 = aad[n + n1 * aad_dim1], abs(d__1)) + (d__2 = aad[n1 + n2 *
		 aad_dim1], abs(d__2));
	x = c3 * s;
	y = x;
/* Computing 2nd power */
	d__1 = s;
	w = -c1 * (d__1 * d__1);
    }
    ++in;
/*                ** look for two consecutive small sub-diagonal elements */
    i__1 = n2;
    for (j = lb; j <= i__1; ++j) {
	i__ = n2 + lb - j;
	z__ = aad[i__ + i__ * aad_dim1];
	r__ = x - z__;
	s = y - z__;
	p = (r__ * s - w) / aad[i__ + 1 + i__ * aad_dim1] + aad[i__ + (i__ + 
		1) * aad_dim1];
	q = aad[i__ + 1 + (i__ + 1) * aad_dim1] - z__ - r__ - s;
	r__ = aad[i__ + 2 + (i__ + 1) * aad_dim1];
	s = abs(p) + abs(q) + abs(r__);
	p /= s;
	q /= s;
	r__ /= s;
	if (i__ == lb) {
	    goto L490;
	}
	uu = (d__1 = aad[i__ + (i__ - 1) * aad_dim1], abs(d__1)) * (abs(q) + 
		abs(r__));
	vv = abs(p) * ((d__1 = aad[i__ - 1 + (i__ - 1) * aad_dim1], abs(d__1))
		 + abs(z__) + (d__2 = aad[i__ + 1 + (i__ + 1) * aad_dim1], 
		abs(d__2)));
	if (uu <= tol * vv) {
	    goto L490;
	}
/* L480: */
    }
L490:
    aad[i__ + 2 + i__ * aad_dim1] = zero;
    i__1 = n;
    for (j = i__ + 3; j <= i__1; ++j) {
	aad[j + (j - 2) * aad_dim1] = zero;
	aad[j + (j - 3) * aad_dim1] = zero;
/* L500: */
    }
/*             ** double qr step involving rows k to n and columns m to n */
    i__1 = n1;
    for (ka = i__; ka <= i__1; ++ka) {
	notlas = ka != n1;
	if (ka == i__) {
	    d__1 = sqrt(p * p + q * q + r__ * r__);
	    s = d_sign(&d__1, &p);
	    if (lb != i__) {
		aad[ka + (ka - 1) * aad_dim1] = -aad[ka + (ka - 1) * aad_dim1]
			;
	    }
	} else {
	    p = aad[ka + (ka - 1) * aad_dim1];
	    q = aad[ka + 1 + (ka - 1) * aad_dim1];
	    r__ = zero;
	    if (notlas) {
		r__ = aad[ka + 2 + (ka - 1) * aad_dim1];
	    }
	    x = abs(p) + abs(q) + abs(r__);
	    if (x == zero) {
		goto L540;
	    }
	    p /= x;
	    q /= x;
	    r__ /= x;
	    d__1 = sqrt(p * p + q * q + r__ * r__);
	    s = d_sign(&d__1, &p);
	    aad[ka + (ka - 1) * aad_dim1] = -s * x;
	}
	p += s;
	x = p / s;
	y = q / s;
	z__ = r__ / s;
	q /= p;
	r__ /= p;
/*                                              ** row modification */
	i__2 = *m;
	for (j = ka; j <= i__2; ++j) {
	    p = aad[ka + j * aad_dim1] + q * aad[ka + 1 + j * aad_dim1];
	    if (notlas) {
		p += r__ * aad[ka + 2 + j * aad_dim1];
		aad[ka + 2 + j * aad_dim1] -= p * z__;
	    }
	    aad[ka + 1 + j * aad_dim1] -= p * y;
	    aad[ka + j * aad_dim1] -= p * x;
/* L510: */
	}
/*                                                 ** column modification */
/* Computing MIN */
	i__3 = n, i__4 = ka + 3;
	i__2 = min(i__3,i__4);
	for (ii = 1; ii <= i__2; ++ii) {
	    p = x * aad[ii + ka * aad_dim1] + y * aad[ii + (ka + 1) * 
		    aad_dim1];
	    if (notlas) {
		p += z__ * aad[ii + (ka + 2) * aad_dim1];
		aad[ii + (ka + 2) * aad_dim1] -= p * r__;
	    }
	    aad[ii + (ka + 1) * aad_dim1] -= p * q;
	    aad[ii + ka * aad_dim1] -= p;
/* L520: */
	}
/*                                          ** accumulate transformations */
	i__2 = k;
	for (ii = l; ii <= i__2; ++ii) {
	    p = x * evecd[ii + ka * evecd_dim1] + y * evecd[ii + (ka + 1) * 
		    evecd_dim1];
	    if (notlas) {
		p += z__ * evecd[ii + (ka + 2) * evecd_dim1];
		evecd[ii + (ka + 2) * evecd_dim1] -= p * r__;
	    }
	    evecd[ii + (ka + 1) * evecd_dim1] -= p * q;
	    evecd[ii + ka * evecd_dim1] -= p;
/* L530: */
	}
L540:
	;
    }
    goto L410;
/*                     ** all evals found, now backsubstitute real vector */
L550:
    if (rnorm != zero) {
	for (n = *m; n >= 1; --n) {
	    n2 = n;
	    aad[n + n * aad_dim1] = one;
	    for (i__ = n - 1; i__ >= 1; --i__) {
		w = aad[i__ + i__ * aad_dim1] - evald[n];
		if (w == zero) {
		    w = tol * rnorm;
		}
		r__ = aad[i__ + n * aad_dim1];
		i__1 = n - 1;
		for (j = n2; j <= i__1; ++j) {
		    r__ += aad[i__ + j * aad_dim1] * aad[j + n * aad_dim1];
/* L560: */
		}
		aad[i__ + n * aad_dim1] = -r__ / w;
		n2 = i__;
/* L570: */
	    }
/* L580: */
	}
/*                      ** end backsubstitution vectors of isolated evals */
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ < l || i__ > k) {
		i__2 = *m;
		for (j = i__; j <= i__2; ++j) {
		    evecd[i__ + j * evecd_dim1] = aad[i__ + j * aad_dim1];
/* L590: */
		}
	    }
/* L600: */
	}
/*                                   ** multiply by transformation matrix */
	if (k != 0) {
	    i__1 = l;
	    for (j = *m; j >= i__1; --j) {
		i__2 = k;
		for (i__ = l; i__ <= i__2; ++i__) {
		    z__ = zero;
		    i__3 = min(j,k);
		    for (n = l; n <= i__3; ++n) {
			z__ += evecd[i__ + n * evecd_dim1] * aad[n + j * 
				aad_dim1];
/* L610: */
		    }
		    evecd[i__ + j * evecd_dim1] = z__;
/* L620: */
		}
/* L630: */
	    }
	}
    }
    i__1 = k;
    for (i__ = l; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    evecd[i__ + j * evecd_dim1] *= wkd[i__];
/* L640: */
	}
/* L650: */
    }
/*                           ** interchange rows if permutations occurred */
    for (i__ = l - 1; i__ >= 1; --i__) {
	j = (integer) wkd[i__];
	if (i__ != j) {
	    i__1 = *m;
	    for (n = 1; n <= i__1; ++n) {
		repl = evecd[i__ + n * evecd_dim1];
		evecd[i__ + n * evecd_dim1] = evecd[j + n * evecd_dim1];
		evecd[j + n * evecd_dim1] = repl;
/* L660: */
	    }
	}
/* L670: */
    }
    i__1 = *m;
    for (i__ = k + 1; i__ <= i__1; ++i__) {
	j = (integer) wkd[i__];
	if (i__ != j) {
	    i__2 = *m;
	    for (n = 1; n <= i__2; ++n) {
		repl = evecd[i__ + n * evecd_dim1];
		evecd[i__ + n * evecd_dim1] = evecd[j + n * evecd_dim1];
		evecd[j + n * evecd_dim1] = repl;
/* L680: */
	    }
	}
/* L690: */
    }
/*                         ** put results into output arrays */
L700:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	eval[j] = evald[j];
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
	    evec[j + k * evec_dim1] = evecd[j + k * evecd_dim1];
/* L710: */
	}
/* L720: */
    }
    return 0;
} /* asymtx_ */

/* Subroutine */ int cmpint_(real *fbeam, real *gc, real *kk, integer *layru, 
	real *ll, logical *lyrcut, integer *mazim, integer *mxcmu, integer *
	mxulv, integer *mxumu, integer *ncut, integer *nn, integer *nstr, 
	logical *plank, integer *ntau, real *taucpr, real *umu0, real *utaupr,
	 real *zz, real *zplk0, real *zplk1, real *uum)
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, uum_dim1, uum_offset, zplk0_dim1, zplk0_offset, 
	    zplk1_dim1, zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer iq, jq, lu, lyu;
    static real zint;

/*          calculates the fourier intensity components at the quadrature */
/*          angles for azimuthal expansion terms (mazim) in eq. sd(2), */
/*          stwl(6) */


/*    i n p u t    v a r i a b l e s: */

/*       kk      :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b) */

/*       gc      :  eigenvectors at polar quadrature angles in eq. sc(1) */

/*       ll      :  constants of integration in eq. sc(1), obtained */
/*                  by solving scaled version of eq. sc(5); */
/*                  exponential term of eq. sc(12) not included */

/*       lyrcut  :  logical flag for truncation of computational layer */

/*       mazim   :  order of azimuthal component */

/*       ncut    :  number of computational layer where absorption */
/*                  optical depth exceeds abscut */

/*       nn      :  order of double-gauss quadrature (nstr/2) */

/*       taucpr  :  cumulative optical depth (delta-m-scaled) */

/*       utaupr  :  optical depths of user output levels in delta-m */
/*                  coordinates;  equal to utau if no delta-m */

/*       zz      :  beam source vectors in eq. ss(19), stwl(24b) */

/*       zplk0   :  thermal source vectors z0, by solving eq. ss(16), */
/*                  y-sub-zero in stwl(26ab) */

/*       zplk1   :  thermal source vectors z1, by solving eq. ss(16), */
/*                  y-sub-one in stwl(26ab) */

/*       (remainder are 'disort' input variables) */


/*    o u t p u t   v a r i a b l e s: */

/*       uum     :  fourier components of the intensity in eq. sd(12) */
/*                    (at polar quadrature angles) */


/*    i n t e r n a l   v a r i a b l e s: */

/*       fact    :  exp( - utaupr / umu0 ) */
/*       zint    :  intensity of m=0 case, in eq. sc(1) */

/*   called by- disort */
/* +-------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
/*                                       ** loop over user levels */
    /* Parameter adjustments */
    --layru;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1;
    zz -= zz_offset;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2);
    gc -= gc_offset;
    --utaupr;
    uum_dim1 = *mxumu;
    uum_offset = 1 + uum_dim1;
    uum -= uum_offset;

    /* Function Body */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	lyu = layru[lu];
	if (*lyrcut && lyu > *ncut) {
	    goto L40;
	}
	i__2 = *nstr;
	for (iq = 1; iq <= i__2; ++iq) {
	    zint = 0.f;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu]));
/* L10: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu - 1]));
/* L20: */
	    }
	    uum[iq + lu * uum_dim1] = zint;
	    if (*fbeam > 0.f) {
		uum[iq + lu * uum_dim1] = zint + zz[iq + lyu * zz_dim1] * exp(
			-utaupr[lu] / *umu0);
	    }
	    if (*plank && *mazim == 0) {
		uum[iq + lu * uum_dim1] = uum[iq + lu * uum_dim1] + zplk0[iq 
			+ lyu * zplk0_dim1] + zplk1[iq + lyu * zplk1_dim1] * 
			utaupr[lu];
	    }
/* L30: */
	}
L40:
	;
    }
    return 0;
} /* cmpint_ */

/* Subroutine */ int fluxes_(real *cmu, real *cwt, real *fbeam, real *gc, 
	real *kk, integer *layru, real *ll, logical *lyrcut, integer *maxulv, 
	integer *mxcmu, integer *mxulv, integer *ncut, integer *nn, integer *
	nstr, integer *ntau, real *pi, logical *prnt, logical *prntu0, real *
	ssalb, real *taucpr, real *umu0, real *utau, real *utaupr, real *xr0, 
	real *xr1, real *zz, real *zplk0, real *zplk1, real *dfdt, real *flup,
	 real *fldn, real *fldir, real *rfldir, real *rfldn, real *uavg, real 
	*u0c)
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, u0c_dim1, u0c_offset, zplk0_dim1, zplk0_offset, 
	    zplk1_dim1, zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double exp(doublereal), acos(doublereal);

    /* Local variables */
    static integer iq, jq, lu, lyu;
    static real ang1, ang2, fact, fnet, zint, dirint, fdntot, plsorc;
    extern /* Subroutine */ int zeroit_(real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___148 = { 0, 6, 0, "(//,21x,a,/,2a,/,2a,/)", 0 };
    static cilist io___159 = { 0, 6, 0, "(f10.4,i7,1p,7e12.3,e14.3)", 0 };
    static cilist io___160 = { 0, 6, 0, "(//,2a)", 0 };
    static cilist io___161 = { 0, 6, 0, "(/,a,f10.4,//,2a)", 0 };
    static cilist io___164 = { 0, 6, 0, "(2(0p,f16.4,f13.5,1p,e14.3))", 0 };


/*       calculates the radiative fluxes, mean intensity, and flux */
/*       derivative with respect to optical depth from the m=0 intensity */
/*       components (the azimuthally-averaged intensity) */


/*    i n p u t     v a r i a b l e s: */

/*       cmu      :  abscissae for gauss quadrature over angle cosine */

/*       cwt      :  weights for gauss quadrature over angle cosine */

/*       gc       :  eigenvectors at polar quadrature angles, sc(1) */

/*       kk       :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b) */

/*       layru    :  layer number of user level utau */

/*       ll       :  constants of integration in eq. sc(1), obtained */
/*                   by solving scaled version of eq. sc(5); */
/*                   exponential term of eq. sc(12) not included */

/*       lyrcut   :  logical flag for truncation of comput. layer */

/*       nn       :  order of double-gauss quadrature (nstr/2) */

/*       ncut     :  number of computational layer where absorption */
/*                   optical depth exceeds abscut */

/*       prntu0   :  true, print azimuthally-averaged intensity at */
/*                   quadrature angles */

/*       taucpr   :  cumulative optical depth (delta-m-scaled) */

/*       utaupr   :  optical depths of user output levels in delta-m */
/*                   coordinates;  equal to utau if no delta-m */

/*       xr0      :  expansion of thermal source function in eq. ss(14), */
/*                   stwl(24c) */

/*       xr1      :  expansion of thermal source function eq. ss(16), */
/*                   stwl(24c) */

/*       zz       :  beam source vectors in eq. ss(19), stwl(24b) */

/*       zplk0    :  thermal source vectors z0, by solving eq. ss(16), */
/*                   y0 in stwl(26b) */

/*       zplk1    :  thermal source vectors z1, by solving eq. ss(16), */
/*                   y1 in stwl(26a) */

/*       (remainder are disort input variables) */


/*    o u t p u t     v a r i a b l e s: */

/*       u0c      :  azimuthally averaged intensities */
/*                   ( at polar quadrature angles ) */

/*       (rfldir, rfldn, flup, dfdt, uavg are disort output variables) */


/*    i n t e r n a l       v a r i a b l e s: */

/*       dirint   :  direct intensity attenuated */
/*       fdntot   :  total downward flux (direct + diffuse) */
/*       fldir    :  direct-beam flux (delta-m scaled) */
/*       fldn     :  diffuse down-flux (delta-m scaled) */
/*       fnet     :  net flux (total-down - diffuse-up) */
/*       fact     :  exp( - utaupr / umu0 ) */
/*       plsorc   :  planck source function (thermal) */
/*       zint     :  intensity of m = 0 case, in eq. sc(1) */

/*   called by- disort */
/*   calls- zeroit */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --uavg;
    --rfldn;
    --rfldir;
    --flup;
    --dfdt;
    --utau;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1;
    zz -= zz_offset;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2);
    gc -= gc_offset;
    --cwt;
    --cmu;
    u0c_dim1 = *mxcmu;
    u0c_offset = 1 + u0c_dim1;
    u0c -= u0c_offset;
    --fldir;
    --fldn;
    --utaupr;
    --layru;
    --prnt;
    --ssalb;
    --xr0;
    --xr1;

    /* Function Body */
    if (prnt[2]) {
	s_wsfe(&io___148);
	do_fio(&c__1, "<----------------------- fluxes ---------------------"
		"-->", (ftnlen)56);
	do_fio(&c__1, "   optical  compu    downward    downward    downward"
		"     ", (ftnlen)58);
	do_fio(&c__1, " upward                    mean      planck   d(net f"
		"lux)", (ftnlen)57);
	do_fio(&c__1, "     depth  layer      direct     diffuse       total"
		"     ", (ftnlen)58);
	do_fio(&c__1, "diffuse         net   intensity      source   / d(op "
		"dep)", (ftnlen)57);
	e_wsfe();
    }
/*                                        ** zero disort output arrays */
    i__1 = *mxulv * *mxcmu;
    zeroit_(&u0c[u0c_offset], &i__1);
    zeroit_(&fldir[1], mxulv);
    zeroit_(&fldn[1], mxulv);
/*                                        ** loop over user levels */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	lyu = layru[lu];
	if (*lyrcut && lyu > *ncut) {
/*                                                ** no radiation reaches */
/*                                                ** this level */
	    fdntot = 0.f;
	    fnet = 0.f;
	    plsorc = 0.f;
	    goto L70;
	}
	if (*fbeam > 0.f) {
	    fact = exp(-utaupr[lu] / *umu0);
	    dirint = *fbeam * fact;
	    fldir[lu] = *umu0 * (*fbeam * fact);
	    rfldir[lu] = *umu0 * *fbeam * exp(-utau[lu] / *umu0);
	} else {
	    dirint = 0.f;
	    fldir[lu] = 0.f;
	    rfldir[lu] = 0.f;
	}
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    zint = 0.f;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu]));
/* L10: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu - 1]));
/* L20: */
	    }
	    u0c[iq + lu * u0c_dim1] = zint;
	    if (*fbeam > 0.f) {
		u0c[iq + lu * u0c_dim1] = zint + zz[iq + lyu * zz_dim1] * 
			fact;
	    }
	    u0c[iq + lu * u0c_dim1] = u0c[iq + lu * u0c_dim1] + zplk0[iq + 
		    lyu * zplk0_dim1] + zplk1[iq + lyu * zplk1_dim1] * utaupr[
		    lu];
	    uavg[lu] += cwt[*nn + 1 - iq] * u0c[iq + lu * u0c_dim1];
	    fldn[lu] += cwt[*nn + 1 - iq] * cmu[*nn + 1 - iq] * u0c[iq + lu * 
		    u0c_dim1];
/* L30: */
	}
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    zint = 0.f;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu]));
/* L40: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu - 1]));
/* L50: */
	    }
	    u0c[iq + lu * u0c_dim1] = zint;
	    if (*fbeam > 0.f) {
		u0c[iq + lu * u0c_dim1] = zint + zz[iq + lyu * zz_dim1] * 
			fact;
	    }
	    u0c[iq + lu * u0c_dim1] = u0c[iq + lu * u0c_dim1] + zplk0[iq + 
		    lyu * zplk0_dim1] + zplk1[iq + lyu * zplk1_dim1] * utaupr[
		    lu];
	    uavg[lu] += cwt[iq - *nn] * u0c[iq + lu * u0c_dim1];
	    flup[lu] += cwt[iq - *nn] * cmu[iq - *nn] * u0c[iq + lu * 
		    u0c_dim1];
/* L60: */
	}
	flup[lu] = *pi * 2.f * flup[lu];
	fldn[lu] = *pi * 2.f * fldn[lu];
	fdntot = fldn[lu] + fldir[lu];
	fnet = fdntot - flup[lu];
	rfldn[lu] = fdntot - rfldir[lu];
	uavg[lu] = (*pi * 2.f * uavg[lu] + dirint) / (*pi * 4.f);
	plsorc = xr0[lyu] + xr1[lyu] * utaupr[lu];
	dfdt[lu] = (1.f - ssalb[lyu]) * 4.f * *pi * (uavg[lu] - plsorc);
L70:
	if (prnt[2]) {
	    s_wsfe(&io___159);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&lyu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rfldir[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rfldn[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&fdntot, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&flup[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&fnet, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&uavg[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&plsorc, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dfdt[lu], (ftnlen)sizeof(real));
	    e_wsfe();
	}
/* L80: */
    }
    if (*prntu0) {
	s_wsfe(&io___160);
	do_fio(&c__1, " ******** azimuthally averaged ", (ftnlen)31);
	do_fio(&c__1, "intensities ( at polar quadrature angles ) *******", (
		ftnlen)50);
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    s_wsfe(&io___161);
	    do_fio(&c__1, " optical depth =", (ftnlen)16);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, "     angle (deg)   cos(angle)     intensity", (
		    ftnlen)43);
	    do_fio(&c__1, "     angle (deg)   cos(angle)     intensity", (
		    ftnlen)43);
	    e_wsfe();
	    i__2 = *nn;
	    for (iq = 1; iq <= i__2; ++iq) {
		ang1 = 180.f / *pi * acos(cmu[(*nn << 1) - iq + 1]);
		ang2 = 180.f / *pi * acos(cmu[iq]);
		s_wsfe(&io___164);
		do_fio(&c__1, (char *)&ang1, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&cmu[(*nn << 1) - iq + 1], (ftnlen)
			sizeof(real));
		do_fio(&c__1, (char *)&u0c[iq + lu * u0c_dim1], (ftnlen)
			sizeof(real));
		do_fio(&c__1, (char *)&ang2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&cmu[iq], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&u0c[iq + *nn + lu * u0c_dim1], (ftnlen)
			sizeof(real));
		e_wsfe();
/* L90: */
	    }
/* L100: */
	}
    }
    return 0;
} /* fluxes_ */

/* Subroutine */ int intcor_(real *dither, real *fbeam, real *flyr, integer *
	layru, logical *lyrcut, integer *maxmom, integer *maxulv, integer *
	maxumu, integer *nmom, integer *ncut, integer *nphi, integer *nstr, 
	integer *ntau, integer *numu, real *oprim, real *phasa, real *phast, 
	real *phasm, real *phirad, real *pi, real *rpd, real *pmom, real *
	ssalb, real *dtauc, real *tauc, real *taucpr, real *umu, real *umu0, 
	real *utau, real *utaupr, real *uu)
{
    /* System generated locals */
    integer pmom_dim1, pmom_offset, uu_dim1, uu_dim2, uu_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2;

    /* Builtin functions */
    double acos(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    static integer k, lc, jp;
    static real pl;
    static integer iu, lu;
    static real plm1, plm2;
    static integer ltau;
    static real ussp, duims, theta0;
    extern doublereal secsca_(real *, real *, integer *, integer *, integer *,
	     integer *, real *, real *, real *, real *, real *, real *, real *
	    , real *, real *);
    static real ctheta, dtheta;
    extern doublereal sinsca_(real *, integer *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *);
    static real thetap, ussndm;

/*       corrects intensity field by using nakajima-tanaka algorithm */
/*       (1988). for more details, see section 3.6 of stwl nasa report. */
/*                i n p u t   v a r i a b l e s */

/*       dither  10 times machine precision */

/*       dtauc   computational-layer optical depths */

/*       fbeam   incident beam radiation at top */

/*       flyr    separated fraction in delta-m method */

/*       layru   index of utau in multi-layered system */

/*       lyrcut  logical flag for truncation of computational layer */

/*       nmom    number of phase function legendre coefficients supplied */

/*       ncut    total number of computational layers considered */

/*       nphi    number of user azimuthal angles */

/*       nstr    number of polar quadrature angles */

/*       ntau    number of user-defined optical depths */

/*       numu    number of user polar angles */

/*       oprim   delta-m-scaled single-scatter albedo */

/*       phirad  azimuthal angles in radians */

/*       pmom    phase function legendre coefficients (k, lc) */
/*                   k = 0 to nmom, lc = 1 to nlyr with pmom(0,lc)=1 */

/*       rpd     pi/180 */

/*       ssalb   single scattering albedo at computational layers */

/*       tauc    optical thickness at computational levels */

/*       taucpr  delta-m-scaled optical thickness */

/*       umu     cosine of emergent angle */

/*       umu0    cosine of incident zenith angle */

/*       utau    user defined optical depths */

/*       utaupr  delta-m-scaled version of utau */

/*                o u t p u t   v a r i a b l e s */

/*       uu      corrected intensity field; uu(iu,lu,j) */
/*                         iu=1,numu; lu=1,ntau; j=1,nphi */

/*                i n t e r n a l   v a r i a b l e s */

/*       ctheta  cosine of scattering angle */
/*       dtheta  angle (degrees) to define aureole region as */
/*                    direction of beam source +/- dtheta */
/*       phasa   actual (exact) phase function */
/*       phasm   delta-m-scaled phase function */
/*       phast   phase function used in tms correction; actual phase */
/*                    function divided by (1-flyr*ssalb) */
/*       pl      ordinary legendre polynomial of degree l, p-sub-l */
/*       plm1    ordinary legendre polynomial of degree l-1, p-sub-(l-1) */
/*       plm2    ordinary legendre polynomial of degree l-2, p-sub-(l-2) */
/*       theta0  incident zenith angle (degrees) */
/*       thetap  emergent angle (degrees) */
/*       ussndm  single-scattered intensity computed by using exact */
/*                   phase function and scaled optical depth */
/*                   (first term in stwl(68a)) */
/*       ussp    single-scattered intensity from delta-m method */
/*                   (second term in stwl(68a)) */
/*       duims   intensity correction term from ims method */
/*                   (delta-i-sub-ims in stwl(a.19)) */

/*   called by- disort */
/*   calls- sinsca, secsca */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --flyr;
    --layru;
    pmom_dim1 = *maxmom - 0 + 1;
    pmom_offset = 0 + pmom_dim1;
    pmom -= pmom_offset;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2);
    uu -= uu_offset;
    --oprim;
    --phasa;
    --phast;
    --phasm;
    --phirad;
    --ssalb;
    --dtauc;
    --umu;
    --utau;
    --utaupr;

    /* Function Body */
    dtheta = 10.f;
/*                                ** start loop over zenith angles */
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	if (umu[iu] < 0.f) {
/*                                ** calculate zenith angles of icident */
/*                                ** and emerging directions */
	    theta0 = acos(-(*umu0)) / *rpd;
	    thetap = acos(umu[iu]) / *rpd;
	}
/*                                ** start loop over azimuth angles */
	i__2 = *nphi;
	for (jp = 1; jp <= i__2; ++jp) {
/*                                ** calculate cosine of scattering */
/*                                ** angle, eq. stwl(4) */
/* Computing 2nd power */
	    r__1 = *umu0;
/* Computing 2nd power */
	    r__2 = umu[iu];
	    ctheta = -(*umu0) * umu[iu] + sqrt((1.f - r__1 * r__1) * (1.f - 
		    r__2 * r__2)) * cos(phirad[jp]);
/*                                ** initialize phase function */
	    i__3 = *ncut;
	    for (lc = 1; lc <= i__3; ++lc) {
		phasa[lc] = 1.f;
		phasm[lc] = 1.f;
/* L10: */
	    }
/*                                ** initialize legendre poly. recurrence */
	    plm1 = 1.f;
	    plm2 = 0.f;
	    i__3 = *nmom;
	    for (k = 1; k <= i__3; ++k) {
/*                                ** calculate legendre polynomial of */
/*                                ** p-sub-l by upward recurrence */
		pl = (((k << 1) - 1) * ctheta * plm1 - (k - 1) * plm2) / k;
		plm2 = plm1;
		plm1 = pl;
/*                                ** calculate actual phase function */
		i__4 = *ncut;
		for (lc = 1; lc <= i__4; ++lc) {
		    phasa[lc] += ((k << 1) + 1) * pl * pmom[k + lc * 
			    pmom_dim1];
/* L20: */
		}
/*                                ** calculate delta-m transformed */
/*                                ** phase function */
		if (k <= *nstr - 1) {
		    i__4 = *ncut;
		    for (lc = 1; lc <= i__4; ++lc) {
			phasm[lc] += ((k << 1) + 1) * pl * (pmom[k + lc * 
				pmom_dim1] - flyr[lc]) / (1.f - flyr[lc]);
/* L30: */
		    }
		}
/* L40: */
	    }
/*                                ** apply tms method, eq. stwl(68) */
	    i__3 = *ncut;
	    for (lc = 1; lc <= i__3; ++lc) {
		phast[lc] = phasa[lc] / (1.f - flyr[lc] * ssalb[lc]);
/* L70: */
	    }
	    i__3 = *ntau;
	    for (lu = 1; lu <= i__3; ++lu) {
		if (! (*lyrcut) || layru[lu] < *ncut) {
		    ussndm = sinsca_(dither, &layru[lu], ncut, &phast[1], &
			    ssalb[1], taucpr, &umu[iu], umu0, &utaupr[lu], 
			    fbeam, pi);
		    ussp = sinsca_(dither, &layru[lu], ncut, &phasm[1], &
			    oprim[1], taucpr, &umu[iu], umu0, &utaupr[lu], 
			    fbeam, pi);
		    uu[iu + (lu + jp * uu_dim2) * uu_dim1] = uu[iu + (lu + jp 
			    * uu_dim2) * uu_dim1] + ussndm - ussp;
		}
/* L80: */
	    }
	    if (umu[iu] < 0.f && (r__1 = theta0 - thetap, dabs(r__1)) <= 
		    dtheta) {
/*                                ** emerging direction is in the aureole */
/*                                ** (theta0 +/- dtheta). apply ims */
/*                                ** method for correction of secondary */
/*                                ** scattering below top level. */
		ltau = 1;
		if (utau[1] <= *dither) {
		    ltau = 2;
		}
		i__3 = *ntau;
		for (lu = ltau; lu <= i__3; ++lu) {
		    if (! (*lyrcut) || layru[lu] < *ncut) {
			duims = secsca_(&ctheta, &flyr[1], &layru[lu], maxmom,
				 nmom, nstr, &pmom[pmom_offset], &ssalb[1], &
				dtauc[1], tauc, &umu[iu], umu0, &utau[lu], 
				fbeam, pi);
			uu[iu + (lu + jp * uu_dim2) * uu_dim1] -= duims;
		    }
/* L90: */
		}
	    }
/*                                ** end loop over azimuth angles */
/* L100: */
	}
/*                                ** end loop over zenith angles */
/* L110: */
    }
    return 0;
} /* intcor_ */

doublereal secsca_(real *ctheta, real *flyr, integer *layru, integer *maxmom, 
	integer *nmom, integer *nstr, real *pmom, real *ssalb, real *dtauc, 
	real *tauc, real *umu, real *umu0, real *utau, real *fbeam, real *pi)
{
    /* System generated locals */
    integer pmom_dim1, pmom_offset, i__1, i__2;
    real ret_val, r__1, r__2;

    /* Local variables */
    static integer k;
    static real pl;
    static integer lyr;
    static real plm1, plm2, fbar, gbar, wbar, dtau, stau, zero, umu0p, pspike;
    extern doublereal xifunc_(real *, real *, real *, real *);

/*          calculates secondary scattered intensity of eq. stwl (a7) */
/*                i n p u t   v a r i a b l e s */
/*        ctheta  cosine of scattering angle */

/*        dtauc   computational-layer optical depths */

/*        flyr    separated fraction f in delta-m method */

/*        layru   index of utau in multi-layered system */

/*        maxmom  maximum number of phase function moment coefficients */

/*        nmom    number of phase function legendre coefficients supplied */

/*        nstr    number of polar quadrature angles */

/*        pmom    phase function legendre coefficients (k, lc) */
/*                k = 0 to nmom, lc = 1 to nlyr, with pmom(0,lc)=1 */

/*        ssalb   single scattering albedo of computational layers */

/*        tauc    cumulative optical depth at computational layers */

/*        umu     cosine of emergent angle */

/*        umu0    cosine of incident zenith angle */

/*        utau    user defined optical depth for output intensity */

/*        fbeam   incident beam radiation at top */

/*        pi       3.1415... */

/*   local variables */

/*        pspike  2*p"-p"**2, where p" is the residual phase function */
/*        wbar    mean value of single scattering albedo */
/*        fbar    mean value of separated fraction f */
/*        dtau    layer optical depth */
/*        stau    sum of layer optical depths between top of atmopshere */
/*                and layer layru */

/*   called by- intcor */
/*   calls- xifunc */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
    /* Parameter adjustments */
    --flyr;
    pmom_dim1 = *maxmom - 0 + 1;
    pmom_offset = 0 + pmom_dim1;
    pmom -= pmom_offset;
    --ssalb;
    --dtauc;

    /* Function Body */
    zero = 1e-4f;
/*                          ** calculate vertically averaged value of */
/*                          ** single scattering albedo and separated */
/*                          ** fraction f, eq. stwl (a.15) */
    dtau = *utau - tauc[*layru - 1];
    wbar = ssalb[*layru] * dtau;
    fbar = flyr[*layru] * wbar;
    stau = dtau;
    i__1 = *layru - 1;
    for (lyr = 1; lyr <= i__1; ++lyr) {
	wbar += ssalb[lyr] * dtauc[lyr];
	fbar += ssalb[lyr] * dtauc[lyr] * flyr[lyr];
	stau += dtauc[lyr];
/* L10: */
    }
    if (wbar <= zero || fbar <= zero || stau <= zero || *fbeam <= zero) {
	ret_val = 0.f;
	return ret_val;
    }
    fbar /= wbar;
    wbar /= stau;
/*                          ** calculate pspike=(2p"-p"**2) */
    pspike = 1.f;
    gbar = 1.f;
    plm1 = 1.f;
    plm2 = 0.f;
/*                                   ** pspike for l<=2n-1 */
    i__1 = *nstr - 1;
    for (k = 1; k <= i__1; ++k) {
	pl = (((k << 1) - 1) * *ctheta * plm1 - (k - 1) * plm2) / k;
	plm2 = plm1;
	plm1 = pl;
/* Computing 2nd power */
	r__1 = gbar;
	pspike += (gbar * 2.f - r__1 * r__1) * ((k << 1) + 1) * pl;
/* L20: */
    }
/*                                   ** pspike for l>2n-1 */
    i__1 = *nmom;
    for (k = *nstr; k <= i__1; ++k) {
	pl = (((k << 1) - 1) * *ctheta * plm1 - (k - 1) * plm2) / k;
	plm2 = plm1;
	plm1 = pl;
	dtau = *utau - tauc[*layru - 1];
	gbar = pmom[k + *layru * pmom_dim1] * ssalb[*layru] * dtau;
	i__2 = *layru - 1;
	for (lyr = 1; lyr <= i__2; ++lyr) {
	    gbar += pmom[k + lyr * pmom_dim1] * ssalb[lyr] * dtauc[lyr];
/* L30: */
	}
	if (fbar * wbar * stau <= zero) {
	    gbar = 0.f;
	} else {
	    gbar /= fbar * wbar * stau;
	}
/* Computing 2nd power */
	r__1 = gbar;
	pspike += (gbar * 2.f - r__1 * r__1) * ((k << 1) + 1) * pl;
/* L40: */
    }
    umu0p = *umu0 / (1.f - fbar * wbar);
/*                              ** calculate ims correction term, */
/*                              ** eq. stwl (a.13) */
/* Computing 2nd power */
    r__1 = fbar * wbar;
    r__2 = -(*umu);
    ret_val = *fbeam / (*pi * 4.f) * (r__1 * r__1) / (1.f - fbar * wbar) * 
	    pspike * xifunc_(&r__2, &umu0p, &umu0p, utau);
    return ret_val;
} /* secsca_ */

/* Subroutine */ int setdis_(real *cmu, real *cwt, logical *deltam, real *
	dtauc, real *dtaucp, real *expbea, real *fbeam, real *flyr, real *gl, 
	integer *ibcnd, integer *layru, logical *lyrcut, integer *maxmom, 
	integer *maxumu, integer *mxcmu, integer *ncut, integer *nlyr, 
	integer *ntau, integer *nn, integer *nstr, logical *plank, integer *
	numu, logical *onlyfl, logical *corint, real *oprim, real *pmom, real 
	*ssalb, real *tauc, real *taucpr, real *utau, real *utaupr, real *umu,
	 real *umu0, logical *usrtau, logical *usrang)
{
    /* Initialized data */

    static real abscut = 10.f;

    /* System generated locals */
    integer gl_dim1, gl_offset, pmom_dim1, pmom_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static real f;
    static integer k, lc;
    static logical tf;
    static integer iq, iu, lu;
    static real abstau;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen), qgausn_(
	    integer *, real *, real *);
    static real yessct;

/*          perform miscellaneous setting-up operations */

/*    input :  all are disort input variables (see doc file) */


/*    o u t p u t     v a r i a b l e s: */

/*       ntau,utau   if usrtau = false (defined in disort.doc) */
/*       numu,umu    if usrang = false (defined in disort.doc) */

/*       cmu,cwt     computational polar angles and */
/*                   corresponding quadrature weights */

/*       expbea      transmission of direct beam */

/*       flyr        separated fraction in delta-m method */

/*       gl          phase function legendre coefficients multiplied */
/*                   by (2l+1) and single-scatter albedo */

/*       layru       computational layer in which utau falls */

/*       lyrcut      flag as to whether radiation will be zeroed */
/*                   below layer ncut */

/*       ncut        computational layer where absorption */
/*                   optical depth first exceeds  abscut */

/*       nn          nstr / 2 */

/*       oprim       delta-m-scaled single-scatter albedo */

/*       taucpr      delta-m-scaled optical depth */

/*       utaupr      delta-m-scaled version of  utau */

/*   called by- disort */
/*   calls- qgausn, errmsg */
/* --------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --dtauc;
    --dtaucp;
    --flyr;
    --layru;
    pmom_dim1 = *maxmom - 0 + 1;
    pmom_offset = 0 + pmom_dim1;
    pmom -= pmom_offset;
    --umu;
    gl_dim1 = *mxcmu - 0 + 1;
    gl_offset = 0 + gl_dim1;
    gl -= gl_offset;
    --cwt;
    --cmu;
    --oprim;
    --ssalb;
    --utau;
    --utaupr;

    /* Function Body */
    if (! (*usrtau)) {
/*                              ** set output levels at computational */
/*                              ** layer boundaries */
	*ntau = *nlyr + 1;
	i__1 = *ntau - 1;
	for (lc = 0; lc <= i__1; ++lc) {
	    utau[lc + 1] = tauc[lc];
/* L10: */
	}
    }
/*                        ** apply delta-m scaling and move description */
/*                        ** of computational layers to local variables */
    expbea[0] = 1.f;
    taucpr[0] = 0.f;
    abstau = 0.f;
    yessct = 0.f;
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	yessct += ssalb[lc];
	pmom[lc * pmom_dim1] = 1.f;
	if (abstau < abscut) {
	    *ncut = lc;
	}
	abstau += (1.f - ssalb[lc]) * dtauc[lc];
	if (! (*deltam)) {
	    oprim[lc] = ssalb[lc];
	    dtaucp[lc] = dtauc[lc];
	    taucpr[lc] = tauc[lc];
	    i__2 = *nstr - 1;
	    for (k = 0; k <= i__2; ++k) {
		gl[k + lc * gl_dim1] = ((k << 1) + 1) * oprim[lc] * pmom[k + 
			lc * pmom_dim1];
/* L20: */
	    }
	    f = 0.f;
	} else {
/*                                    ** do delta-m transformation */
	    f = pmom[*nstr + lc * pmom_dim1];
	    oprim[lc] = ssalb[lc] * (1.f - f) / (1.f - f * ssalb[lc]);
	    dtaucp[lc] = (1.f - f * ssalb[lc]) * dtauc[lc];
	    taucpr[lc] = taucpr[lc - 1] + dtaucp[lc];
	    i__2 = *nstr - 1;
	    for (k = 0; k <= i__2; ++k) {
		gl[k + lc * gl_dim1] = ((k << 1) + 1) * oprim[lc] * (pmom[k + 
			lc * pmom_dim1] - f) / (1.f - f);
/* L30: */
	    }
	}
	flyr[lc] = f;
	expbea[lc] = 0.f;
	if (*fbeam > 0.f) {
	    expbea[lc] = exp(-taucpr[lc] / *umu0);
	}
/* L40: */
    }
/*                      ** if no thermal emission, cut off medium below */
/*                      ** absorption optical depth = abscut ( note that */
/*                      ** delta-m transformation leaves absorption */
/*                      ** optical depth invariant ).  not worth the */
/*                      ** trouble for one-layer problems, though. */
    *lyrcut = FALSE_;
    if (abstau >= abscut && ! (*plank) && *ibcnd != 1 && *nlyr > 1) {
	*lyrcut = TRUE_;
    }
    if (! (*lyrcut)) {
	*ncut = *nlyr;
    }
/*                             ** set arrays defining location of user */
/*                             ** output levels within delta-m-scaled */
/*                             ** computational mesh */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	i__2 = *nlyr;
	for (lc = 1; lc <= i__2; ++lc) {
	    if (utau[lu] >= tauc[lc - 1] && utau[lu] <= tauc[lc]) {
		goto L60;
	    }
/* L50: */
	}
	lc = *nlyr;
L60:
	utaupr[lu] = utau[lu];
	if (*deltam) {
	    utaupr[lu] = taucpr[lc - 1] + (1.f - ssalb[lc] * flyr[lc]) * (
		    utau[lu] - tauc[lc - 1]);
	}
	layru[lu] = lc;
/* L70: */
    }
/*                      ** calculate computational polar angle cosines */
/*                      ** and associated quadrature weights for gaussian */
/*                      ** quadrature on the interval (0,1) (upward) */
    *nn = *nstr / 2;
    qgausn_(nn, &cmu[1], &cwt[1]);
/*                                  ** downward (neg) angles and weights */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	cmu[iq + *nn] = -cmu[iq];
	cwt[iq + *nn] = cwt[iq];
/* L80: */
    }
    if (*fbeam > 0.f) {
/*                     ** compare beam angle to comput. angles */
	tf = TRUE_;
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    if ((r__1 = *umu0 - cmu[iq], dabs(r__1)) / *umu0 < 1e-4f) {
		errmsg_("setdis--beam angle=computational angle; change nstr",
			 &tf, (ftnlen)51);
	    }
/* L90: */
	}
    }
    if (! (*usrang) || *onlyfl && *maxumu >= *nstr) {
/*                                   ** set output polar angles to */
/*                                   ** computational polar angles */
	*numu = *nstr;
	i__1 = *nn;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[iu] = -cmu[*nn + 1 - iu];
/* L100: */
	}
	i__1 = *nstr;
	for (iu = *nn + 1; iu <= i__1; ++iu) {
	    umu[iu] = cmu[iu - *nn];
/* L110: */
	}
    }
    if (*usrang && *ibcnd == 1) {
/*                               ** shift positive user angle cosines to */
/*                               ** upper locations and put negatives */
/*                               ** in lower locations */
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[iu + *numu] = umu[iu];
/* L120: */
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[iu] = -umu[(*numu << 1) + 1 - iu];
/* L130: */
	}
	*numu <<= 1;
    }
/*                               ** turn off intensity correction when */
/*                               ** only fluxes are calculated, there */
/*                               ** is no beam source, no scattering, */
/*                               ** or delta-m transformation is not */
/*                               ** applied */

    if (*onlyfl || *fbeam == 0.f || yessct == 0.f || ! (*deltam)) {
	*corint = FALSE_;
    }
    return 0;
} /* setdis_ */

/* Subroutine */ int setmtx_(real *bdr, real *cband, real *cmu, real *cwt, 
	real *delm0, real *dtaucp, real *gc, real *kk, logical *lamber, 
	logical *lyrcut, integer *mi, integer *mi9m2, integer *mxcmu, integer 
	*ncol, integer *ncut, integer *nnlyri, integer *nn, integer *nstr, 
	real *taucpr, real *wk)
{
    /* System generated locals */
    integer bdr_dim1, bdr_offset, cband_dim1, cband_offset, gc_dim1, gc_dim2, 
	    gc_offset, kk_dim1, kk_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer k, lc, iq, jq, lda, ncd;
    static real sum;
    static integer jcol;
    static real expa;
    static integer irow, nncol, nshift;
    extern /* Subroutine */ int zeroit_(real *, integer *);

/*        calculate coefficient matrix for the set of equations */
/*        obtained from the boundary conditions and the continuity- */
/*        of-intensity-at-layer-interface equations;  store in the */
/*        special banded-matrix format required by linpack routines */


/*    i n p u t      v a r i a b l e s: */

/*       bdr      :  surface bidirectional reflectivity */

/*       cmu,cwt     abscissae, weights for gauss quadrature */
/*                   over angle cosine */

/*       delm0    :  kronecker delta, delta-sub-m0 */

/*       gc       :  eigenvectors at polar quadrature angles, sc(1) */

/*       kk       :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b) */

/*       lyrcut   :  logical flag for truncation of computational layers */

/*       nn       :  number of streams in a hemisphere (nstr/2) */

/*       ncut     :  total number of computational layers considered */

/*       taucpr   :  cumulative optical depth (delta-m-scaled) */

/*       (remainder are disort input variables) */


/*   o u t p u t     v a r i a b l e s: */

/*       cband    :  left-hand side matrix of linear system eq. sc(5), */
/*                   scaled by eq. sc(12); in banded form required */
/*                   by linpack solution routines */

/*       ncol     :  number of columns in cband */


/*   i n t e r n a l    v a r i a b l e s: */

/*       irow     :  points to row in cband */
/*       jcol     :  points to position in layer block */
/*       lda      :  row dimension of cband */
/*       ncd      :  number of diagonals below or above main diagonal */
/*       nshift   :  for positioning number of rows in band storage */
/*       wk       :  temporary storage for exp evaluations */


/*   band storage */

/*      linpack requires band matrices to be input in a special */
/*      form where the elements of each diagonal are moved up or */
/*      down (in their column) so that each diagonal becomes a row. */
/*      (the column locations of diagonal elements are unchanged.) */

/*      example:  if the original matrix is */

/*          11 12 13  0  0  0 */
/*          21 22 23 24  0  0 */
/*           0 32 33 34 35  0 */
/*           0  0 43 44 45 46 */
/*           0  0  0 54 55 56 */
/*           0  0  0  0 65 66 */

/*      then its linpack input form would be: */

/*           *  *  *  +  +  +  , * = not used */
/*           *  * 13 24 35 46  , + = used for pivoting */
/*           * 12 23 34 45 56 */
/*          11 22 33 44 55 66 */
/*          21 32 43 54 65  * */

/*      if a is a band matrix, the following program segment */
/*      will convert it to the form (abd) required by linpack */
/*      band-matrix routines: */

/*               n  = (column dimension of a, abd) */
/*               ml = (band width below the diagonal) */
/*               mu = (band width above the diagonal) */
/*               m = ml + mu + 1 */
/*               do j = 1, n */
/*                  i1 = max(1, j-mu) */
/*                  i2 = min(n, j+ml) */
/*                  do i = i1, i2 */
/*                     k = i - j + m */
/*                     abd(k,j) = a(i,j) */
/*                  end do */
/*               end do */

/*      this uses rows  ml+1  through  2*ml+mu+1  of abd. */
/*      the total number of rows needed in abd is  2*ml+mu+1 . */
/*      in the example above, n = 6, ml = 1, mu = 2, and the */
/*      row dimension of abd must be >= 5. */


/*   called by- disort, albtrn */
/*   calls- zeroit */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --dtaucp;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    --wk;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2);
    gc -= gc_offset;
    --cwt;
    --cmu;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1;
    cband -= cband_offset;

    /* Function Body */
    i__1 = *mi9m2 * *nnlyri;
    zeroit_(&cband[cband_offset], &i__1);
    ncd = *nn * 3 - 1;
    lda = ncd * 3 + 1;
    nshift = lda - (*nstr << 1) + 1;
    *ncol = 0;
/*                         ** use continuity conditions of eq. stwj(17) */
/*                         ** to form coefficient matrix in stwj(20); */
/*                         ** employ scaling transformation stwj(22) */
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    wk[iq] = exp(kk[iq + lc * kk_dim1] * dtaucp[lc]);
/* L10: */
	}
	jcol = 0;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ++(*ncol);
	    irow = nshift - jcol;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		cband[irow + *nstr + *ncol * cband_dim1] = gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1];
		cband[irow + *ncol * cband_dim1] = -gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1] * wk[iq];
		++irow;
/* L20: */
	    }
	    ++jcol;
/* L30: */
	}
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    ++(*ncol);
	    irow = nshift - jcol;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		cband[irow + *nstr + *ncol * cband_dim1] = gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1] * wk[*nstr + 1 - iq];
		cband[irow + *ncol * cband_dim1] = -gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1];
		++irow;
/* L40: */
	    }
	    ++jcol;
/* L50: */
	}
/* L60: */
    }
/*                  ** use top boundary condition of stwj(20a) for */
/*                  ** first layer */
    jcol = 0;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	expa = exp(kk[iq + kk_dim1] * taucpr[1]);
	irow = nshift - jcol + *nn;
	for (jq = *nn; jq >= 1; --jq) {
	    cband[irow + (jcol + 1) * cband_dim1] = gc[jq + (iq + gc_dim2) * 
		    gc_dim1] * expa;
	    ++irow;
/* L70: */
	}
	++jcol;
/* L80: */
    }
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	irow = nshift - jcol + *nn;
	for (jq = *nn; jq >= 1; --jq) {
	    cband[irow + (jcol + 1) * cband_dim1] = gc[jq + (iq + gc_dim2) * 
		    gc_dim1];
	    ++irow;
/* L90: */
	}
	++jcol;
/* L100: */
    }
/*                           ** use bottom boundary condition of */
/*                           ** stwj(20c) for last layer */
    nncol = *ncol - *nstr;
    jcol = 0;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	++nncol;
	irow = nshift - jcol + *nstr;
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    if (*lyrcut || *lamber && *delm0 == 0.f) {
/*                          ** no azimuthal-dependent intensity if lam- */
/*                          ** bert surface; no intensity component if */
/*                          ** truncated bottom layer */
		cband[irow + nncol * cband_dim1] = gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1];
	    } else {
		sum = 0.f;
		i__3 = *nn;
		for (k = 1; k <= i__3; ++k) {
		    sum += cwt[k] * cmu[k] * bdr[jq - *nn + k * bdr_dim1] * 
			    gc[*nn + 1 - k + (iq + *ncut * gc_dim2) * gc_dim1]
			    ;
/* L110: */
		}
		cband[irow + nncol * cband_dim1] = gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1] - (*delm0 + 1.f) * sum;
	    }
	    ++irow;
/* L120: */
	}
	++jcol;
/* L130: */
    }
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	++nncol;
	irow = nshift - jcol + *nstr;
	expa = wk[*nstr + 1 - iq];
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    if (*lyrcut || *lamber && *delm0 == 0.f) {
		cband[irow + nncol * cband_dim1] = gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1] * expa;
	    } else {
		sum = 0.f;
		i__3 = *nn;
		for (k = 1; k <= i__3; ++k) {
		    sum += cwt[k] * cmu[k] * bdr[jq - *nn + k * bdr_dim1] * 
			    gc[*nn + 1 - k + (iq + *ncut * gc_dim2) * gc_dim1]
			    ;
/* L140: */
		}
		cband[irow + nncol * cband_dim1] = (gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1] - (*delm0 + 1.f) * sum) * expa;
	    }
	    ++irow;
/* L150: */
	}
	++jcol;
/* L160: */
    }
    return 0;
} /* setmtx_ */

doublereal sinsca_(real *dither, integer *layru, integer *nlyr, real *phase, 
	real *omega, real *tau, real *umu, real *umu0, real *utau, real *
	fbeam, real *pi)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer lyr;
    static real exp0, exp1;

/*        calculates single-scattered intensity from eqs. stwl (65b,d,e) */
/*                i n p u t   v a r i a b l e s */
/*        dither   10 times machine precision */

/*        layru    index of utau in multi-layered system */

/*        nlyr     number of sublayers */

/*        phase    phase functions of sublayers */

/*        omega    single scattering albedos of sublayers */

/*        tau      optical thicknesses of sublayers */

/*        umu      cosine of emergent angle */

/*        umu0     cosine of incident zenith angle */

/*        utau     user defined optical depth for output intensity */

/*        fbeam   incident beam radiation at top */

/*        pi       3.1415... */

/*   called by- intcor */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --omega;
    --phase;

    /* Function Body */
    ret_val = 0.f;
    exp0 = exp(-(*utau) / *umu0);
    if ((r__1 = *umu + *umu0, dabs(r__1)) <= *dither) {
/*                                 ** calculate downward intensity when */
/*                                 ** umu=umu0, eq. stwl (65e) */
	i__1 = *layru - 1;
	for (lyr = 1; lyr <= i__1; ++lyr) {
	    ret_val += omega[lyr] * phase[lyr] * (tau[lyr] - tau[lyr - 1]);
/* L10: */
	}
	ret_val = *fbeam / (*pi * 4.f * *umu0) * exp0 * (ret_val + omega[*
		layru] * phase[*layru] * (*utau - tau[*layru - 1]));
	return ret_val;
    }
    if (*umu > 0.f) {
/*                                 ** upward intensity, eq. stwl (65b) */
	i__1 = *nlyr;
	for (lyr = *layru; lyr <= i__1; ++lyr) {
	    exp1 = exp(-((tau[lyr] - *utau) / *umu + tau[lyr] / *umu0));
	    ret_val += omega[lyr] * phase[lyr] * (exp0 - exp1);
	    exp0 = exp1;
/* L20: */
	}
    } else {
/*                                 ** downward intensity, eq. stwl (65d) */
	for (lyr = *layru; lyr >= 1; --lyr) {
	    exp1 = exp(-((tau[lyr - 1] - *utau) / *umu + tau[lyr - 1] / *umu0)
		    );
	    ret_val += omega[lyr] * phase[lyr] * (exp0 - exp1);
	    exp0 = exp1;
/* L30: */
	}
    }
    ret_val = *fbeam / (*pi * 4.f * (*umu / *umu0 + 1.f)) * ret_val;
    return ret_val;
} /* sinsca_ */

/* Subroutine */ int soleig_(real *amb, real *apb, real *array, real *cmu, 
	real *cwt, real *gl, integer *mi, integer *mazim, integer *mxcmu, 
	integer *nn, integer *nstr, real *ylmc, real *cc, real *evecc, real *
	eval, real *kk, real *gc, doublereal *aad, doublereal *eveccd, 
	doublereal *evald, doublereal *wkd)
{
    /* System generated locals */
    integer amb_dim1, amb_offset, apb_dim1, apb_offset, array_dim1, 
	    array_offset, cc_dim1, cc_offset, evecc_dim1, evecc_offset, 
	    gc_dim1, gc_offset, ylmc_dim1, ylmc_offset, aad_dim1, aad_offset, 
	    eveccd_dim1, eveccd_offset, i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static integer l;
    static logical tf;
    static integer iq, jq, kq, ier;
    static real sum, beta, alpha, gpmigm, gpplgm;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen), asymtx_(
	    real *, real *, real *, integer *, integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___228 = { 0, 6, 0, "(//,a,i4,a)", 0 };


/*         solves eigenvalue/vector problem necessary to construct */
/*         homogeneous part of discrete ordinate solution; stwj(8b), */
/*         stwl(23f) */
/*         ** note ** eigenvalue problem is degenerate when single */
/*                    scattering albedo = 1;  present way of doing it */
/*                    seems numerically more stable than alternative */
/*                    methods that we tried */


/*   i n p u t     v a r i a b l e s: */

/*       gl     :  delta-m scaled legendre coefficients of phase function */
/*                 (including factors 2l+1 and single-scatter albedo) */

/*       cmu    :  computational polar angle cosines */

/*       cwt    :  weights for quadrature over polar angle cosine */

/*       mazim  :  order of azimuthal component */

/*       nn     :  half the total number of streams */

/*       ylmc   :  normalized associated legendre polynomial */
/*                 at the quadrature angles cmu */

/*       (remainder are disort input variables) */


/*   o u t p u t    v a r i a b l e s: */

/*       cc     :  c-sub-ij in eq. ss(5); needed in ss(15&18) */

/*       eval   :  nn eigenvalues of eq. ss(12), stwl(23f) on return */
/*                 from asymtx but then square roots taken */

/*       evecc  :  nn eigenvectors  (g+) - (g-)  on return */
/*                 from asymtx ( column j corresponds to eval(j) ) */
/*                 but then  (g+) + (g-)  is calculated from ss(10), */
/*                 g+  and  g-  are separated, and  g+  is stacked on */
/*                 top of  g-  to form nstr eigenvectors of ss(7) */

/*       gc     :  permanent storage for all nstr eigenvectors, but */
/*                 in an order corresponding to kk */

/*       kk     :  permanent storage for all nstr eigenvalues of ss(7), */
/*                 but re-ordered with negative values first ( square */
/*                 roots of eval taken and negatives added ) */


/*   i n t e r n a l   v a r i a b l e s: */

/*       amb,apb :  matrices (alpha-beta), (alpha+beta) in reduced */
/*                    eigenvalue problem */
/*       array   :  complete coefficient matrix of reduced eigenvalue */
/*                    problem: (alfa+beta)*(alfa-beta) */
/*       gpplgm  :  (g+) + (g-) (cf. eqs. ss(10-11)) */
/*       gpmigm  :  (g+) - (g-) (cf. eqs. ss(10-11)) */
/*       wkd     :  scratch array required by asymtx */

/*   called by- disort, albtrn */
/*   calls- asymtx, errmsg */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
/*                             ** calculate quantities in eqs. ss(5-6), */
/*                             ** stwl(8b,15,23f) */
    /* Parameter adjustments */
    --evald;
    eveccd_dim1 = *mi;
    eveccd_offset = 1 + eveccd_dim1;
    eveccd -= eveccd_offset;
    aad_dim1 = *mi;
    aad_offset = 1 + aad_dim1;
    aad -= aad_offset;
    --eval;
    array_dim1 = *mi;
    array_offset = 1 + array_dim1;
    array -= array_offset;
    apb_dim1 = *mi;
    apb_offset = 1 + apb_dim1;
    apb -= apb_offset;
    amb_dim1 = *mi;
    amb_offset = 1 + amb_dim1;
    amb -= amb_offset;
    --wkd;
    gc_dim1 = *mxcmu;
    gc_offset = 1 + gc_dim1;
    gc -= gc_offset;
    --kk;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1;
    evecc -= evecc_offset;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1;
    cc -= cc_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1;
    ylmc -= ylmc_offset;
    --cwt;
    --cmu;

    /* Function Body */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    sum = 0.f;
	    i__3 = *nstr - 1;
	    for (l = *mazim; l <= i__3; ++l) {
		sum += gl[l] * ylmc[l + iq * ylmc_dim1] * ylmc[l + jq * 
			ylmc_dim1];
/* L10: */
	    }
	    cc[iq + jq * cc_dim1] = sum * .5f * cwt[jq];
/* L20: */
	}
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
/*                             ** fill remainder of array using symmetry */
/*                             ** relations  c(-mui,muj) = c(mui,-muj) */
/*                             ** and        c(-mui,-muj) = c(mui,muj) */
	    cc[iq + *nn + jq * cc_dim1] = cc[iq + (jq + *nn) * cc_dim1];
	    cc[iq + *nn + (jq + *nn) * cc_dim1] = cc[iq + jq * cc_dim1];
/*                                       ** get factors of coeff. matrix */
/*                                       ** of reduced eigenvalue problem */
	    alpha = cc[iq + jq * cc_dim1] / cmu[iq];
	    beta = cc[iq + (jq + *nn) * cc_dim1] / cmu[iq];
	    amb[iq + jq * amb_dim1] = alpha - beta;
	    apb[iq + jq * apb_dim1] = alpha + beta;
/* L30: */
	}
	amb[iq + iq * amb_dim1] -= 1.f / cmu[iq];
	apb[iq + iq * apb_dim1] -= 1.f / cmu[iq];
/* L40: */
    }
/*                      ** finish calculation of coefficient matrix of */
/*                      ** reduced eigenvalue problem:  get matrix */
/*                      ** product (alfa+beta)*(alfa-beta); ss(12), */
/*                      ** stwl(23f) */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    sum = 0.f;
	    i__3 = *nn;
	    for (kq = 1; kq <= i__3; ++kq) {
		sum += apb[iq + kq * apb_dim1] * amb[kq + jq * amb_dim1];
/* L50: */
	    }
	    array[iq + jq * array_dim1] = sum;
/* L60: */
	}
/* L70: */
    }
/*                      ** find (real) eigenvalues and eigenvectors */
    asymtx_(&array[array_offset], &evecc[evecc_offset], &eval[1], nn, mi, 
	    mxcmu, &ier, &wkd[1], &aad[aad_offset], &eveccd[eveccd_offset], &
	    evald[1]);
    tf = TRUE_;
    if (ier > 0) {
	s_wsfe(&io___228);
	do_fio(&c__1, " asymtx--eigenvalue no. ", (ftnlen)24);
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	do_fio(&c__1, "  didnt converge.  lower-numbered eigenvalues wrong.", 
		(ftnlen)52);
	e_wsfe();
	errmsg_("asymtx--convergence problems", &tf, (ftnlen)28);
    }
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	eval[iq] = sqrt((r__1 = eval[iq], dabs(r__1)));
	kk[iq + *nn] = eval[iq];
/*                                      ** add negative eigenvalue */
	kk[*nn + 1 - iq] = -eval[iq];
/* L80: */
    }
/*                          ** find eigenvectors (g+) + (g-) from ss(10) */
/*                          ** and store temporarily in apb array */
    i__1 = *nn;
    for (jq = 1; jq <= i__1; ++jq) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    sum = 0.f;
	    i__3 = *nn;
	    for (kq = 1; kq <= i__3; ++kq) {
		sum += amb[iq + kq * amb_dim1] * evecc[kq + jq * evecc_dim1];
/* L90: */
	    }
	    apb[iq + jq * apb_dim1] = sum / eval[jq];
/* L100: */
	}
/* L110: */
    }
    i__1 = *nn;
    for (jq = 1; jq <= i__1; ++jq) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    gpplgm = apb[iq + jq * apb_dim1];
	    gpmigm = evecc[iq + jq * evecc_dim1];
/*                                ** recover eigenvectors g+,g- from */
/*                                ** their sum and difference; stack them */
/*                                ** to get eigenvectors of full system */
/*                                ** ss(7) (jq = eigenvector number) */
	    evecc[iq + jq * evecc_dim1] = (gpplgm + gpmigm) * .5f;
	    evecc[iq + *nn + jq * evecc_dim1] = (gpplgm - gpmigm) * .5f;
/*                                ** eigenvectors corresponding to */
/*                                ** negative eigenvalues (corresp. to */
/*                                ** reversing sign of 'k' in ss(10) ) */
	    gpplgm = -gpplgm;
	    evecc[iq + (jq + *nn) * evecc_dim1] = (gpplgm + gpmigm) * .5f;
	    evecc[iq + *nn + (jq + *nn) * evecc_dim1] = (gpplgm - gpmigm) * 
		    .5f;
	    gc[iq + *nn + (jq + *nn) * gc_dim1] = evecc[iq + jq * evecc_dim1];
	    gc[*nn + 1 - iq + (jq + *nn) * gc_dim1] = evecc[iq + *nn + jq * 
		    evecc_dim1];
	    gc[iq + *nn + (*nn + 1 - jq) * gc_dim1] = evecc[iq + (jq + *nn) * 
		    evecc_dim1];
	    gc[*nn + 1 - iq + (*nn + 1 - jq) * gc_dim1] = evecc[iq + *nn + (
		    jq + *nn) * evecc_dim1];
/* L120: */
	}
/* L130: */
    }
    return 0;
} /* soleig_ */

/* Subroutine */ int solve0_(real *b, real *bdr, real *bem, real *bplank, 
	real *cband, real *cmu, real *cwt, real *expbea, real *fbeam, real *
	fisot, integer *ipvt, logical *lamber, real *ll, logical *lyrcut, 
	integer *mazim, integer *mi, integer *mi9m2, integer *mxcmu, integer *
	ncol, integer *ncut, integer *nn, integer *nstr, integer *nnlyri, 
	real *pi, real *tplank, real *taucpr, real *umu0, real *z__, real *zz,
	 real *zplk0, real *zplk1)
{
    /* System generated locals */
    integer bdr_dim1, bdr_offset, cband_dim1, cband_offset, ll_dim1, 
	    ll_offset, zplk0_dim1, zplk0_offset, zplk1_dim1, zplk1_offset, 
	    zz_dim1, zz_offset, i__1, i__2;

    /* Local variables */
    static integer lc;
    static logical tf;
    static integer iq, it, jq, ncd;
    static real sum;
    static integer ipnt;
    extern /* Subroutine */ int sgbco_(real *, integer *, integer *, integer *
	    , integer *, integer *, real *, real *);
    static real rcond;
    extern /* Subroutine */ int sgbsl_(real *, integer *, integer *, integer *
	    , integer *, integer *, real *, integer *), errmsg_(char *, 
	    logical *, ftnlen), zeroit_(real *, integer *);

/*        construct right-hand side vector b for general boundary */
/*        conditions stwj(17) and solve system of equations obtained */
/*        from the boundary conditions and the continuity-of- */
/*        intensity-at-layer-interface equations. */
/*        thermal emission contributes only in azimuthal independence. */


/*    i n p u t      v a r i a b l e s: */

/*       bdr      :  surface bidirectional reflectivity */

/*       bem      :  surface bidirectional emissivity */

/*       bplank   :  bottom boundary thermal emission */

/*       cband    :  left-hand side matrix of linear system eq. sc(5), */
/*                   scaled by eq. sc(12); in banded form required */
/*                   by linpack solution routines */

/*       cmu,cwt  :  abscissae, weights for gauss quadrature */
/*                   over angle cosine */

/*       expbea   :  transmission of incident beam, exp(-taucpr/umu0) */

/*       lyrcut   :  logical flag for truncation of computational layers */

/*       mazim    :  order of azimuthal component */

/*       ncol     :  number of columns in cband */

/*       nn       :  order of double-gauss quadrature (nstr/2) */

/*       ncut     :  total number of computational layers considered */

/*       tplank   :  top boundary thermal emission */

/*       taucpr   :  cumulative optical depth (delta-m-scaled) */

/*       zz       :  beam source vectors in eq. ss(19), stwl(24b) */

/*       zplk0    :  thermal source vectors z0, by solving eq. ss(16), */
/*                   y0 in stwl(26b) */

/*       zplk1    :  thermal source vectors z1, by solving eq. ss(16), */
/*                   y1 in stwl(26a) */

/*       (remainder are disort input variables) */


/*    o u t p u t     v a r i a b l e s: */

/*       b        :  right-hand side vector of eq. sc(5) going into */
/*                   sgbsl; returns as solution vector of eq. sc(12), */
/*                   constants of integration without exponential term */

/*      ll        :  permanent storage for b, but re-ordered */


/*   i n t e r n a l    v a r i a b l e s: */

/*       ipvt     :  integer vector of pivot indices */
/*       it       :  pointer for position in  b */
/*       ncd      :  number of diagonals below or above main diagonal */
/*       rcond    :  indicator of singularity for cband */
/*       z        :  scratch array required by sgbco */

/*   called by- disort */
/*   calls- zeroit, sgbco, errmsg, sgbsl */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
    /* Parameter adjustments */
    --ipvt;
    --bem;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1;
    zz -= zz_offset;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    --cwt;
    --cmu;
    --z__;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1;
    cband -= cband_offset;
    --b;

    /* Function Body */
    zeroit_(&b[1], nnlyri);
/*                              ** construct b,  stwj(20a,c) for */
/*                              ** parallel beam + bottom reflection + */
/*                              ** thermal emission at top and/or bottom */
    if (*mazim > 0 && *fbeam > 0.f) {
/*                                         ** azimuth-dependent case */
/*                                         ** (never called if fbeam = 0) */
	if (*lyrcut || *lamber) {
/*               ** no azimuthal-dependent intensity for lambert surface; */
/*               ** no intensity component for truncated bottom layer */
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
/*                                                  ** top boundary */
		b[iq] = -zz[*nn + 1 - iq + zz_dim1];
/*                                                  ** bottom boundary */
		b[*ncol - *nn + iq] = -zz[iq + *nn + *ncut * zz_dim1] * 
			expbea[*ncut];
/* L10: */
	    }
	} else {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		b[iq] = -zz[*nn + 1 - iq + zz_dim1];
		sum = 0.f;
		i__2 = *nn;
		for (jq = 1; jq <= i__2; ++jq) {
		    sum += cwt[jq] * cmu[jq] * bdr[iq + jq * bdr_dim1] * zz[*
			    nn + 1 - jq + *ncut * zz_dim1] * expbea[*ncut];
/* L20: */
		}
		b[*ncol - *nn + iq] = sum;
		if (*fbeam > 0.f) {
		    b[*ncol - *nn + iq] = sum + (bdr[iq] * *umu0 * *fbeam / *
			    pi - zz[iq + *nn + *ncut * zz_dim1]) * expbea[*
			    ncut];
		}
/* L30: */
	    }
	}
/*                             ** continuity condition for layer */
/*                             ** interfaces of eq. stwj(20b) */
	it = *nn;
	i__1 = *ncut - 1;
	for (lc = 1; lc <= i__1; ++lc) {
	    i__2 = *nstr;
	    for (iq = 1; iq <= i__2; ++iq) {
		++it;
		b[it] = (zz[iq + (lc + 1) * zz_dim1] - zz[iq + lc * zz_dim1]) 
			* expbea[lc];
/* L40: */
	    }
/* L50: */
	}
    } else {
/*                                   ** azimuth-independent case */
	if (*fbeam == 0.f) {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
/*                                      ** top boundary */
		b[iq] = -zplk0[*nn + 1 - iq + zplk0_dim1] + *fisot + *tplank;
/* L60: */
	    }
	    if (*lyrcut) {
/*                               ** no intensity component for truncated */
/*                               ** bottom layer */
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
/*                                      ** bottom boundary */
		    b[*ncol - *nn + iq] = -zplk0[iq + *nn + *ncut * 
			    zplk0_dim1] - zplk1[iq + *nn + *ncut * zplk1_dim1]
			     * taucpr[*ncut];
/* L70: */
		}
	    } else {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    sum = 0.f;
		    i__2 = *nn;
		    for (jq = 1; jq <= i__2; ++jq) {
			sum += cwt[jq] * cmu[jq] * bdr[iq + jq * bdr_dim1] * (
				zplk0[*nn + 1 - jq + *ncut * zplk0_dim1] + 
				zplk1[*nn + 1 - jq + *ncut * zplk1_dim1] * 
				taucpr[*ncut]);
/* L80: */
		    }
		    b[*ncol - *nn + iq] = sum * 2.f + bem[iq] * *bplank - 
			    zplk0[iq + *nn + *ncut * zplk0_dim1] - zplk1[iq + 
			    *nn + *ncut * zplk1_dim1] * taucpr[*ncut];
/* L90: */
		}
	    }
/*                             ** continuity condition for layer */
/*                             ** interfaces, stwj(20b) */
	    it = *nn;
	    i__1 = *ncut - 1;
	    for (lc = 1; lc <= i__1; ++lc) {
		i__2 = *nstr;
		for (iq = 1; iq <= i__2; ++iq) {
		    ++it;
		    b[it] = zplk0[iq + (lc + 1) * zplk0_dim1] - zplk0[iq + lc 
			    * zplk0_dim1] + (zplk1[iq + (lc + 1) * zplk1_dim1]
			     - zplk1[iq + lc * zplk1_dim1]) * taucpr[lc];
/* L100: */
		}
/* L110: */
	    }
	} else {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		b[iq] = -zz[*nn + 1 - iq + zz_dim1] - zplk0[*nn + 1 - iq + 
			zplk0_dim1] + *fisot + *tplank;
/* L120: */
	    }
	    if (*lyrcut) {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    b[*ncol - *nn + iq] = -zz[iq + *nn + *ncut * zz_dim1] * 
			    expbea[*ncut] - zplk0[iq + *nn + *ncut * 
			    zplk0_dim1] - zplk1[iq + *nn + *ncut * zplk1_dim1]
			     * taucpr[*ncut];
/* L130: */
		}
	    } else {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    sum = 0.f;
		    i__2 = *nn;
		    for (jq = 1; jq <= i__2; ++jq) {
			sum += cwt[jq] * cmu[jq] * bdr[iq + jq * bdr_dim1] * (
				zz[*nn + 1 - jq + *ncut * zz_dim1] * expbea[*
				ncut] + zplk0[*nn + 1 - jq + *ncut * 
				zplk0_dim1] + zplk1[*nn + 1 - jq + *ncut * 
				zplk1_dim1] * taucpr[*ncut]);
/* L140: */
		    }
		    b[*ncol - *nn + iq] = sum * 2.f + (bdr[iq] * *umu0 * *
			    fbeam / *pi - zz[iq + *nn + *ncut * zz_dim1]) * 
			    expbea[*ncut] + bem[iq] * *bplank - zplk0[iq + *
			    nn + *ncut * zplk0_dim1] - zplk1[iq + *nn + *ncut 
			    * zplk1_dim1] * taucpr[*ncut];
/* L150: */
		}
	    }
	    it = *nn;
	    i__1 = *ncut - 1;
	    for (lc = 1; lc <= i__1; ++lc) {
		i__2 = *nstr;
		for (iq = 1; iq <= i__2; ++iq) {
		    ++it;
		    b[it] = (zz[iq + (lc + 1) * zz_dim1] - zz[iq + lc * 
			    zz_dim1]) * expbea[lc] + zplk0[iq + (lc + 1) * 
			    zplk0_dim1] - zplk0[iq + lc * zplk0_dim1] + (
			    zplk1[iq + (lc + 1) * zplk1_dim1] - zplk1[iq + lc 
			    * zplk1_dim1]) * taucpr[lc];
/* L160: */
		}
/* L170: */
	    }
	}
    }
/*                     ** find l-u (lower/upper triangular) decomposition */
/*                     ** of band matrix cband and test if it is nearly */
/*                     ** singular (note: cband is destroyed) */
/*                     ** (cband is in linpack packed format) */
    rcond = 0.f;
    ncd = *nn * 3 - 1;
    sgbco_(&cband[cband_offset], mi9m2, ncol, &ncd, &ncd, &ipvt[1], &rcond, &
	    z__[1]);
    tf = FALSE_;
    if (rcond + 1.f == 1.f) {
	errmsg_("solve0--sgbco says matrix near singular", &tf, (ftnlen)39);
    }
/*                   ** solve linear system with coeff matrix cband */
/*                   ** and r.h. side(s) b after cband has been l-u */
/*                   ** decomposed.  solution is returned in b. */
    sgbsl_(&cband[cband_offset], mi9m2, ncol, &ncd, &ncd, &ipvt[1], &b[1], &
	    c__0);
/*                   ** zero cband (it may contain 'foreign' */
/*                   ** elements upon returning from linpack); */
/*                   ** necessary to prevent errors */
    i__1 = *mi9m2 * *nnlyri;
    zeroit_(&cband[cband_offset], &i__1);
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	ipnt = lc * *nstr - *nn;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ll[*nn + 1 - iq + lc * ll_dim1] = b[ipnt + 1 - iq];
	    ll[iq + *nn + lc * ll_dim1] = b[iq + ipnt];
/* L180: */
	}
/* L190: */
    }
    return 0;
} /* solve0_ */

/* Subroutine */ int surfac_(real *albedo, real *delm0, real *cmu, real *
	fbeam, logical *lamber, integer *iref, integer *mi, integer *mazim, 
	integer *mxumu, integer *nn, integer *numu, logical *onlyfl, real *pi,
	 real *umu, real *umu0, logical *usrang, real *surf_pr__, real *bdr, 
	real *emu, real *bem, real *rmu)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer bdr_dim1, bdr_offset, rmu_dim1, rmu_offset, i__1, i__2;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer k, jg, iq, jq, iu;
    static real gmu[50], gwt[50], sum, dref;
    extern doublereal bdref_(real *, real *, real *, real *, integer *);
    static real pigmu;
    extern /* Subroutine */ int qgausn_(integer *, real *, real *), zeroit_(
	    real *, integer *);

/*       computes user's surface bidirectional properties, stwl(41) */

/*   i n p u t     v a r i a b l e s: */

/*       cmu    :  computational polar angle cosines (gaussian) */

/*       delm0  :  kronecker delta, delta-sub-m0 */

/*       mazim  :  order of azimuthal component */

/*       nn     :  order of double-gauss quadrature (nstr/2) */

/*       iref   : bidirectional reflectance options */
/*             0 - Lambert */
/*             1 - Hapke's BDR model */
/*             2 - Breon's BDR model; combination of Li + Roujean */
/*             3 - Roujean's BDR model */
/*             4 - Cox and Munk glint model */

/*   SURF_PR : Wavelength dependent surface properties array */
/*             IREF= 0 - Lambert albedo */
/*             IREF= 1 - Hapke : HH, W */
/*             IREF= 2 - Breon's BDR model: k0, k1, k2 */
/*             IREF= 3 - Roujean's BDR model: k0, k1, k2 */
/*             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw */

/*       (remainder are 'disort' input variables) */

/*    o u t p u t     v a r i a b l e s: */

/*       bdr :  fourier expansion coefficient of surface bidirectional */
/*                 reflectivity (computational angles) */

/*       rmu :  surface bidirectional reflectivity (user angles) */

/*       bem :  surface directional emissivity (computational angles) */

/*       emu :  surface directional emissivity (user angles) */

/*    i n t e r n a l     v a r i a b l e s: */
/*       dref   :  directional reflectivity */

/*       nmug   :  number of angle cosine quadrature points on (-1,1) */
/*                 for integrating bidirectional reflectivity to get */
/*                 directional emissivity (it is necessary to use a */
/*                 quadrature set distinct from the computational angles, */
/*                 because the computational angles may not be dense */
/*                 enough -- i.e. 'nstr' may be too small-- to give an */
/*                 accurate approximation for the integration). */

/*       gmu    :  the 'nmug' angle cosine quadrature points on (0,1) */

/*       gwt    :  the 'nmug' angle cosine quadrature weights on (0,1) */

/*   called by- disort */
/*   calls- qgausn, bdref, zeroit */
/* +---------------------------------------------------------------------+ */
/*     .. parameters .. */
/*     .. */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */

/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. local arrays .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --cmu;
    --bem;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    rmu_dim1 = *mxumu;
    rmu_offset = 1 + rmu_dim1 * 0;
    rmu -= rmu_offset;
    --emu;
    --umu;
    --surf_pr__;

    /* Function Body */
    if (pass1) {
	pass1 = FALSE_;
	qgausn_(&c__25, gmu, gwt);
	for (k = 1; k <= 25; ++k) {
	    gmu[k + 24] = -gmu[k - 1];
	    gwt[k + 24] = gwt[k - 1];
/* L10: */
	}
    }
    i__1 = *mi * (*mi + 1);
    zeroit_(&bdr[bdr_offset], &i__1);
    zeroit_(&bem[1], mi);
/*                             ** compute fourier expansion coefficient */
/*                             ** of surface bidirectional reflectance */
/*                             ** at computational angles eq. stwl (41) */
    if (*lamber && *mazim == 0) {
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    bem[iq] = 1.f - *albedo;
	    i__2 = *nn;
	    for (jq = 0; jq <= i__2; ++jq) {
		bdr[iq + jq * bdr_dim1] = *albedo;
/* L20: */
	    }
/* L30: */
	}
    } else if (! (*lamber)) {
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    i__2 = *nn;
	    for (jq = 1; jq <= i__2; ++jq) {
		sum = 0.f;
		for (k = 1; k <= 50; ++k) {
		    pigmu = *pi * gmu[k - 1];
		    sum += gwt[k - 1] * bdref_(&cmu[iq], &cmu[jq], &pigmu, &
			    surf_pr__[1], iref) * cos(*mazim * pigmu);
/* L40: */
		}
		bdr[iq + jq * bdr_dim1] = (2.f - *delm0) * .5f * sum;
/* L50: */
	    }
	    if (*fbeam > 0.f) {
		sum = 0.f;
		for (k = 1; k <= 50; ++k) {
		    pigmu = *pi * gmu[k - 1];
		    sum += gwt[k - 1] * bdref_(&cmu[iq], umu0, &pigmu, &
			    surf_pr__[1], iref) * cos(*mazim * pigmu);
/* L60: */
		}
		bdr[iq] = (2.f - *delm0) * .5f * sum;
	    }
/* L70: */
	}
	if (*mazim == 0) {
/*                             ** integrate bidirectional reflectivity */
/*                             ** at reflection polar angle cosines -cmu- */
/*                             ** and incident angle cosines -gmu- to get */
/*                             ** directional emissivity at computational */
/*                             ** angle cosines -cmu-. */
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		dref = 0.f;
		for (jg = 1; jg <= 50; ++jg) {
		    pigmu = *pi * gmu[jg - 1];
		    sum = 0.f;
		    for (k = 1; k <= 25; ++k) {
			sum += gwt[k - 1] * gmu[k - 1] * bdref_(&cmu[iq], &
				gmu[k - 1], &pigmu, &surf_pr__[1], iref);
/* L80: */
		    }
		    dref += gwt[jg - 1] * sum;
/* L90: */
		}
		bem[iq] = 1.f - dref;
/* L100: */
	    }
	}
    }
/*                             ** compute fourier expansion coefficient */
/*                             ** of surface bidirectional reflectance */
/*                             ** at user angles eq. stwl (41) */
    if (! (*onlyfl) && *usrang) {
	zeroit_(&emu[1], mxumu);
	i__1 = *mxumu * (*mi + 1);
	zeroit_(&rmu[rmu_offset], &i__1);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    if (umu[iu] > 0.f) {
		if (*lamber && *mazim == 0) {
		    i__2 = *nn;
		    for (iq = 0; iq <= i__2; ++iq) {
			rmu[iu + iq * rmu_dim1] = *albedo;
/* L110: */
		    }
		    emu[iu] = 1.f - *albedo;
		} else if (! (*lamber)) {
		    i__2 = *nn;
		    for (iq = 1; iq <= i__2; ++iq) {
			sum = 0.f;
			for (k = 1; k <= 50; ++k) {
			    pigmu = *pi * gmu[k - 1];
			    sum += gwt[k - 1] * bdref_(&umu[iu], &cmu[iq], &
				    pigmu, &surf_pr__[1], iref) * cos(*mazim *
				     pigmu);
/* L120: */
			}
			rmu[iu + iq * rmu_dim1] = (2.f - *delm0) * .5f * sum;
/* L130: */
		    }
		    if (*fbeam > 0.f) {
			sum = 0.f;
			for (k = 1; k <= 50; ++k) {
			    pigmu = *pi * gmu[k - 1];
			    sum += gwt[k - 1] * bdref_(&umu[iu], umu0, &pigmu,
				     &surf_pr__[1], iref) * cos(*mazim * 
				    pigmu);
/* L140: */
			}
			rmu[iu] = (2.f - *delm0) * .5f * sum;
		    }
		    if (*mazim == 0) {
/*                               ** integrate bidirectional reflectivity */
/*                               ** at reflection angle cosines -umu- and */
/*                               ** incident angle cosines -gmu- to get */
/*                               ** directional emissivity at */
/*                               ** user angle cosines -umu-. */
			dref = 0.f;
			for (jg = 1; jg <= 50; ++jg) {
			    pigmu = *pi * gmu[jg - 1];
			    sum = 0.f;
			    for (k = 1; k <= 25; ++k) {
				sum += gwt[k - 1] * gmu[k - 1] * bdref_(&umu[
					iu], &gmu[k - 1], &pigmu, &surf_pr__[
					1], iref);
/* L150: */
			    }
			    dref += gwt[jg - 1] * sum;
/* L160: */
			}
			emu[iu] = 1.f - dref;
		    }
		}
	    }
/* L170: */
	}
    }
    return 0;
} /* surfac_ */

/* Subroutine */ int terpev_(real *cwt, real *evecc, real *gl, real *gu, 
	integer *mazim, integer *mxcmu, integer *mxumu, integer *nn, integer *
	nstr, integer *numu, real *wk, real *ylmc, real *ylmu)
{
    /* System generated locals */
    integer evecc_dim1, evecc_offset, gu_dim1, gu_offset, ylmc_dim1, 
	    ylmc_offset, ylmu_dim1, ylmu_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer l, iq, jq, iu;
    static real sum;

/*         interpolate eigenvectors to user angles; eq sd(8) */
/*   called by- disort, albtrn */
/* --------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
    /* Parameter adjustments */
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1;
    ylmc -= ylmc_offset;
    --wk;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1;
    evecc -= evecc_offset;
    --cwt;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1;
    ylmu -= ylmu_offset;
    gu_dim1 = *mxumu;
    gu_offset = 1 + gu_dim1;
    gu -= gu_offset;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr - 1;
	for (l = *mazim; l <= i__2; ++l) {
/*                                   ** inner sum in sd(8) times all */
/*                                   ** factors in outer sum but plm(mu) */
	    sum = 0.f;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		sum += cwt[jq] * ylmc[l + jq * ylmc_dim1] * evecc[jq + iq * 
			evecc_dim1];
/* L10: */
	    }
	    wk[l + 1] = gl[l] * .5f * sum;
/* L20: */
	}
/*                                    ** finish outer sum in sd(8) */
/*                                    ** and store eigenvectors */
	i__2 = *numu;
	for (iu = 1; iu <= i__2; ++iu) {
	    sum = 0.f;
	    i__3 = *nstr - 1;
	    for (l = *mazim; l <= i__3; ++l) {
		sum += wk[l + 1] * ylmu[l + iu * ylmu_dim1];
/* L30: */
	    }
	    if (iq <= *nn) {
		gu[iu + (iq + *nn) * gu_dim1] = sum;
	    }
	    if (iq > *nn) {
		gu[iu + (*nstr + 1 - iq) * gu_dim1] = sum;
	    }
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* terpev_ */

/* Subroutine */ int terpso_(real *cwt, real *delm0, real *fbeam, real *gl, 
	integer *mazim, integer *mxcmu, logical *plank, integer *numu, 
	integer *nstr, real *oprim, real *pi, real *ylm0, real *ylmc, real *
	ylmu, real *psi0, real *psi1, real *xr0, real *xr1, real *z0, real *
	z1, real *zj, real *zbeam, real *z0u, real *z1u)
{
    /* System generated locals */
    integer ylmc_dim1, ylmc_offset, ylmu_dim1, ylmu_offset, i__1, i__2;

    /* Local variables */
    static integer iq, jq, iu;
    static real sum, sum0, sum1, fact, psum, psum0, psum1;

/*         interpolates source functions to user angles, eq. stwl(30) */


/*    i n p u t      v a r i a b l e s: */

/*       cwt    :  weights for gauss quadrature over angle cosine */

/*       delm0  :  kronecker delta, delta-sub-m0 */

/*       gl     :  delta-m scaled legendre coefficients of phase function */
/*                 (including factors 2l+1 and single-scatter albedo) */

/*       mazim  :  order of azimuthal component */

/*       oprim  :  single scattering albedo */

/*       xr0    :  expansion of thermal source function, eq. stwl(24d) */

/*       xr1    :  expansion of thermal source function eq. stwl(24d) */

/*       ylm0   :  normalized associated legendre polynomial */
/*                 at the beam angle */

/*       ylmc   :  normalized associated legendre polynomial */
/*                 at the quadrature angles */

/*       ylmu   :  normalized associated legendre polynomial */
/*                 at the user angles */

/*       z0     :  solution vectors z-sub-zero of eq. ss(16), stwl(26a) */

/*       z1     :  solution vectors z-sub-one  of eq. ss(16), stwl(26b) */

/*       zj     :  solution vector z-sub-zero after solving eq. ss(19), */
/*                 stwl(24b) */

/*       (remainder are disort input variables) */


/*    o u t p u t     v a r i a b l e s: */

/*       zbeam  :  incident-beam source function at user angles */

/*       z0u,z1u:  components of a linear-in-optical-depth-dependent */
/*                 source (approximating the planck emission source) */


/*   i n t e r n a l    v a r i a b l e s: */

/*       psi0  :  sum just after square bracket in  eq. sd(9) */
/*       psi1  :  sum in eq. stwl(31d) */

/*   called by- disort */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
    /* Parameter adjustments */
    --zj;
    --z1;
    --z0;
    --psi1;
    --psi0;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1;
    ylmu -= ylmu_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1;
    ylmc -= ylmc_offset;
    --cwt;
    --zbeam;
    --z0u;
    --z1u;

    /* Function Body */
    if (*fbeam > 0.f) {
/*                                  ** beam source terms; eq. sd(9) */
	i__1 = *nstr - 1;
	for (iq = *mazim; iq <= i__1; ++iq) {
	    psum = 0.f;
	    i__2 = *nstr;
	    for (jq = 1; jq <= i__2; ++jq) {
		psum += cwt[jq] * ylmc[iq + jq * ylmc_dim1] * zj[jq];
/* L10: */
	    }
	    psi0[iq + 1] = gl[iq] * .5f * psum;
/* L20: */
	}
	fact = (2.f - *delm0) * *fbeam / (*pi * 4.f);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    sum = 0.f;
	    i__2 = *nstr - 1;
	    for (iq = *mazim; iq <= i__2; ++iq) {
		sum += ylmu[iq + iu * ylmu_dim1] * (psi0[iq + 1] + fact * gl[
			iq] * ylm0[iq]);
/* L30: */
	    }
	    zbeam[iu] = sum;
/* L40: */
	}
    }
    if (*plank && *mazim == 0) {
/*                          ** thermal source terms, stwj(27c), stwl(31c) */

	i__1 = *nstr - 1;
	for (iq = *mazim; iq <= i__1; ++iq) {
	    psum0 = 0.f;
	    psum1 = 0.f;
	    i__2 = *nstr;
	    for (jq = 1; jq <= i__2; ++jq) {
		psum0 += cwt[jq] * ylmc[iq + jq * ylmc_dim1] * z0[jq];
		psum1 += cwt[jq] * ylmc[iq + jq * ylmc_dim1] * z1[jq];
/* L50: */
	    }
	    psi0[iq + 1] = gl[iq] * .5f * psum0;
	    psi1[iq + 1] = gl[iq] * .5f * psum1;
/* L60: */
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    sum0 = 0.f;
	    sum1 = 0.f;
	    i__2 = *nstr - 1;
	    for (iq = *mazim; iq <= i__2; ++iq) {
		sum0 += ylmu[iq + iu * ylmu_dim1] * psi0[iq + 1];
		sum1 += ylmu[iq + iu * ylmu_dim1] * psi1[iq + 1];
/* L70: */
	    }
	    z0u[iu] = sum0 + (1.f - *oprim) * *xr0;
	    z1u[iu] = sum1 + (1.f - *oprim) * *xr1;
/* L80: */
	}
    }
    return 0;
} /* terpso_ */

/* Subroutine */ int upbeam_(real *array, real *cc, real *cmu, real *delm0, 
	real *fbeam, real *gl, integer *ipvt, integer *mazim, integer *mxcmu, 
	integer *nn, integer *nstr, real *pi, real *umu0, real *wk, real *
	ylm0, real *ylmc, real *zj, real *zz)
{
    /* System generated locals */
    integer array_dim1, array_offset, cc_dim1, cc_offset, ylmc_dim1, 
	    ylmc_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    static logical tf;
    static integer iq, jq, job;
    static real sum;
    extern /* Subroutine */ int sgeco_(real *, integer *, integer *, integer *
	    , real *, real *);
    static real rcond;
    extern /* Subroutine */ int sgesl_(real *, integer *, integer *, integer *
	    , real *, integer *), errmsg_(char *, logical *, ftnlen);

/*         finds the incident-beam particular solution of ss(18), */
/*         stwl(24a) */

/*   i n p u t    v a r i a b l e s: */

/*       cc     :  c-sub-ij in eq. ss(5) */

/*       cmu    :  abscissae for gauss quadrature over angle cosine */

/*       delm0  :  kronecker delta, delta-sub-m0 */

/*       gl     :  delta-m scaled legendre coefficients of phase function */
/*                 (including factors 2l+1 and single-scatter albedo) */

/*       mazim  :  order of azimuthal component */

/*       ylm0   :  normalized associated legendre polynomial */
/*                 at the beam angle */

/*       ylmc   :  normalized associated legendre polynomial */
/*                 at the quadrature angles */

/*       (remainder are disort input variables) */


/*   o u t p u t    v a r i a b l e s: */

/*       zj     :  right-hand side vector x-sub-zero in ss(19),stwl(24b); */
/*                 also the solution vector z-sub-zero after solving */
/*                 that system */

/*       zz     :  permanent storage for zj, but re-ordered */


/*   i n t e r n a l    v a r i a b l e s: */

/*       array  :  coefficient matrix in left-hand side of eq. ss(19), */
/*                   stwl(24b) */
/*       ipvt   :  integer vector of pivot indices required by linpack */
/*       wk     :  scratch array required by linpack */

/*   called by- disort */
/*   calls- sgeco, errmsg, sgesl */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
    /* Parameter adjustments */
    --ipvt;
    --zz;
    --zj;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1;
    ylmc -= ylmc_offset;
    --wk;
    --cmu;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1;
    cc -= cc_offset;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1;
    array -= array_offset;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    array[iq + jq * array_dim1] = -cc[iq + jq * cc_dim1];
/* L10: */
	}
	array[iq + iq * array_dim1] = cmu[iq] / *umu0 + 1.f + array[iq + iq * 
		array_dim1];
	sum = 0.f;
	i__2 = *nstr - 1;
	for (k = *mazim; k <= i__2; ++k) {
	    sum += gl[k] * ylmc[k + iq * ylmc_dim1] * ylm0[k];
/* L20: */
	}
	zj[iq] = (2.f - *delm0) * *fbeam * sum / (*pi * 4.f);
/* L30: */
    }
/*                  ** find l-u (lower/upper triangular) decomposition */
/*                  ** of array and see if it is nearly singular */
/*                  ** (note:  array is altered) */
    rcond = 0.f;
    sgeco_(&array[array_offset], mxcmu, nstr, &ipvt[1], &rcond, &wk[1]);
    tf = FALSE_;
    if (rcond + 1.f == 1.f) {
	errmsg_("upbeam--sgeco says matrix near singular", &tf, (ftnlen)39);
    }
/*                ** solve linear system with coeff matrix array */
/*                ** (assumed already l-u decomposed) and r.h. side(s) */
/*                ** zj;  return solution(s) in zj */
    job = 0;
    sgesl_(&array[array_offset], mxcmu, nstr, &ipvt[1], &zj[1], &job);
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zz[iq + *nn] = zj[iq];
	zz[*nn + 1 - iq] = zj[iq + *nn];
/* L40: */
    }
    return 0;
} /* upbeam_ */

/* Subroutine */ int upisot_(real *array, real *cc, real *cmu, integer *ipvt, 
	integer *mxcmu, integer *nn, integer *nstr, real *oprim, real *wk, 
	real *xr0, real *xr1, real *z0, real *z1, real *zplk0, real *zplk1)
{
    /* System generated locals */
    integer array_dim1, array_offset, cc_dim1, cc_offset, i__1, i__2;

    /* Local variables */
    static logical tf;
    static integer iq, jq;
    extern /* Subroutine */ int sgeco_(real *, integer *, integer *, integer *
	    , real *, real *);
    static real rcond;
    extern /* Subroutine */ int sgesl_(real *, integer *, integer *, integer *
	    , real *, integer *), errmsg_(char *, logical *, ftnlen);

/*       finds the particular solution of thermal radiation of stwl(25) */



/*    i n p u t     v a r i a b l e s: */

/*       cc     :  c-sub-ij in eq. ss(5), stwl(8b) */

/*       cmu    :  abscissae for gauss quadrature over angle cosine */

/*       oprim  :  delta-m scaled single scattering albedo */

/*       xr0    :  expansion coefficient b-sub-zero of thermal source */
/*                   function, eq. stwl(24c) */

/*       xr1    :  expansion coefficient b-sub-one of thermal source */
/*                   function eq. stwl(24c) */

/*       (remainder are disort input variables) */


/*    o u t p u t    v a r i a b l e s: */

/*       z0     :  solution vectors z-sub-zero of eq. ss(16), stwl(26a) */

/*       z1     :  solution vectors z-sub-one  of eq. ss(16), stwl(26b) */

/*       zplk0, :  permanent storage for z0,z1, but re-ordered */
/*        zplk1 */


/*   i n t e r n a l    v a r i a b l e s: */

/*       array  :  coefficient matrix in left-hand side of eq. ss(16) */
/*       ipvt   :  integer vector of pivot indices required by linpack */
/*       wk     :  scratch array required by linpack */

/*   called by- disort */
/*   calls- sgeco, errmsg, sgesl */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
    /* Parameter adjustments */
    --ipvt;
    --zplk1;
    --zplk0;
    --z1;
    --z0;
    --wk;
    --cmu;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1;
    cc -= cc_offset;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1;
    array -= array_offset;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    array[iq + jq * array_dim1] = -cc[iq + jq * cc_dim1];
/* L10: */
	}
	array[iq + iq * array_dim1] += 1.f;
	z1[iq] = (1.f - *oprim) * *xr1;
/* L20: */
    }
/*                       ** solve linear equations: same as in upbeam, */
/*                       ** except zj replaced by z1 and z0 */
    rcond = 0.f;
    sgeco_(&array[array_offset], mxcmu, nstr, &ipvt[1], &rcond, &wk[1]);
    tf = FALSE_;
    if (rcond + 1.f == 1.f) {
	errmsg_("upisot--sgeco says matrix near singular", &tf, (ftnlen)39);
    }
    sgesl_(&array[array_offset], mxcmu, nstr, &ipvt[1], &z1[1], &c__0);
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	z0[iq] = (1.f - *oprim) * *xr0 + cmu[iq] * z1[iq];
/* L30: */
    }
    sgesl_(&array[array_offset], mxcmu, nstr, &ipvt[1], &z0[1], &c__0);
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zplk0[iq + *nn] = z0[iq];
	zplk1[iq + *nn] = z1[iq];
	zplk0[*nn + 1 - iq] = z0[iq + *nn];
	zplk1[*nn + 1 - iq] = z1[iq + *nn];
/* L40: */
    }
    return 0;
} /* upisot_ */

/* Subroutine */ int usrint_(real *bplank, real *cmu, real *cwt, real *delm0, 
	real *dtaucp, real *emu, real *expbea, real *fbeam, real *fisot, real 
	*gc, real *gu, real *kk, logical *lamber, integer *layru, real *ll, 
	logical *lyrcut, integer *mazim, integer *mxcmu, integer *mxulv, 
	integer *mxumu, integer *ncut, integer *nlyr, integer *nn, integer *
	nstr, logical *plank, integer *numu, integer *ntau, real *pi, real *
	rmu, real *taucpr, real *tplank, real *umu, real *umu0, real *utaupr, 
	real *wk, real *zbeam, real *z0u, real *z1u, real *zz, real *zplk0, 
	real *zplk1, real *uum)
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, gu_dim1, gu_dim2, gu_offset, kk_dim1,
	     kk_offset, ll_dim1, ll_offset, rmu_dim1, rmu_offset, uum_dim1, 
	    uum_offset, z0u_dim1, z0u_offset, z1u_dim1, z1u_offset, 
	    zbeam_dim1, zbeam_offset, zplk0_dim1, zplk0_offset, zplk1_dim1, 
	    zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer lc, iq, jq, iu, lu;
    static real f0n, f1n, sgn;
    static integer lyu;
    static real exp0, exp1, exp2, fact, dtau, expn, dtau1, dtau2, denom, 
	    bnddfu, bnddir, bndint, palint, dfuint;
    static integer lyrend;
    static logical negumu;
    static real plkint;
    static integer lyrstr;

/*       computes intensity components at user output angles */
/*       for azimuthal expansion terms in eq. sd(2), stwl(6) */


/*   i n p u t    v a r i a b l e s: */

/*       bplank :  integrated planck function for emission from */
/*                 bottom boundary */

/*       cmu    :  abscissae for gauss quadrature over angle cosine */

/*       cwt    :  weights for gauss quadrature over angle cosine */

/*       delm0  :  kronecker delta, delta-sub-m0 */

/*       emu    :  surface directional emissivity (user angles) */

/*       expbea :  transmission of incident beam, exp(-taucpr/umu0) */

/*       gc     :  eigenvectors at polar quadrature angles, sc(1) */

/*       gu     :  eigenvectors interpolated to user polar angles */
/*                    (i.e., g in eq. sc(1) ) */

/*       kk     :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b) */

/*       layru  :  layer number of user level utau */

/*       ll     :  constants of integration in eq. sc(1), obtained */
/*                 by solving scaled version of eq. sc(5); */
/*                 exponential term of eq. sc(12) not included */

/*       lyrcut :  logical flag for truncation of computational layer */

/*       mazim  :  order of azimuthal component */

/*       ncut   :  total number of computational layers considered */

/*       nn     :  order of double-gauss quadrature (nstr/2) */

/*       rmu    :  surface bidirectional reflectivity (user angles) */

/*       taucpr :  cumulative optical depth (delta-m-scaled) */

/*       tplank :  integrated planck function for emission from */
/*                 top boundary */

/*       utaupr :  optical depths of user output levels in delta-m */
/*                 coordinates;  equal to utau if no delta-m */

/*       z0u    :  z-sub-zero in eq. ss(16) interpolated to user */
/*                 angles from an equation derived from ss(16), */
/*                 y-sub-zero on stwl(26b) */

/*       z1u    :  z-sub-one in eq. ss(16) interpolated to user */
/*                 angles from an equation derived from ss(16), */
/*                 y-sub-one in stwl(26a) */

/*       zz     :  beam source vectors in eq. ss(19), stwl(24b) */

/*       zplk0  :  thermal source vectors z0, by solving eq. ss(16), */
/*                 y-sub-zero in stwl(26) */

/*       zplk1  :  thermal source vectors z1, by solving eq. ss(16), */
/*                 y-sub-one in stwl(26) */

/*       zbeam  :  incident-beam source vectors */

/*       (remainder are disort input variables) */


/*    o u t p u t    v a r i a b l e s: */

/*       uum    :  azimuthal components of the intensity in eq. stwj(5), */
/*                 stwl(6) */


/*    i n t e r n a l    v a r i a b l e s: */

/*       bnddir :  direct intensity down at the bottom boundary */
/*       bnddfu :  diffuse intensity down at the bottom boundary */
/*       bndint :  intensity attenuated at both boundaries, stwj(25-6) */
/*       dtau   :  optical depth of a computational layer */
/*       lyrend :  end layer of integration */
/*       lyrstr :  start layer of integration */
/*       palint :  intensity component from parallel beam */
/*       plkint :  intensity component from planck source */
/*       wk     :  scratch vector for saving exp evaluations */

/*       all the exponential factors ( exp1, expn,... etc.) */
/*       come from the substitution of constants of integration in */
/*       eq. sc(12) into eqs. s1(8-9).  they all have negative */
/*       arguments so there should never be overflow problems. */

/*   called by- disort */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
/*                          ** incorporate constants of integration into */
/*                          ** interpolated eigenvectors */
    /* Parameter adjustments */
    --dtaucp;
    --layru;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1;
    zz -= zz_offset;
    --wk;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2);
    gc -= gc_offset;
    --cwt;
    --cmu;
    --utaupr;
    uum_dim1 = *mxumu;
    uum_offset = 1 + uum_dim1;
    uum -= uum_offset;
    z1u_dim1 = *mxumu;
    z1u_offset = 1 + z1u_dim1;
    z1u -= z1u_offset;
    z0u_dim1 = *mxumu;
    z0u_offset = 1 + z0u_dim1;
    z0u -= z0u_offset;
    zbeam_dim1 = *mxumu;
    zbeam_offset = 1 + zbeam_dim1;
    zbeam -= zbeam_offset;
    rmu_dim1 = *mxumu;
    rmu_offset = 1 + rmu_dim1 * 0;
    rmu -= rmu_offset;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2);
    gu -= gu_offset;
    --emu;
    --umu;

    /* Function Body */
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	i__2 = *nstr;
	for (iq = 1; iq <= i__2; ++iq) {
	    i__3 = *numu;
	    for (iu = 1; iu <= i__3; ++iu) {
		gu[iu + (iq + lc * gu_dim2) * gu_dim1] *= ll[iq + lc * 
			ll_dim1];
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }
/*                           ** loop over levels at which intensities */
/*                           ** are desired ('user output levels') */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	if (*fbeam > 0.f) {
	    exp0 = exp(-utaupr[lu] / *umu0);
	}
	lyu = layru[lu];
/*                              ** loop over polar angles at which */
/*                              ** intensities are desired */
	i__2 = *numu;
	for (iu = 1; iu <= i__2; ++iu) {
	    if (*lyrcut && lyu > *ncut) {
		goto L150;
	    }
	    negumu = umu[iu] < 0.f;
	    if (negumu) {
		lyrstr = 1;
		lyrend = lyu - 1;
		sgn = -1.f;
	    } else {
		lyrstr = lyu + 1;
		lyrend = *ncut;
		sgn = 1.f;
	    }
/*                          ** for downward intensity, integrate from top */
/*                          ** to lyu-1 in eq. s1(8); for upward, */
/*                          ** integrate from bottom to lyu+1 in s1(9) */
	    palint = 0.f;
	    plkint = 0.f;
	    i__3 = lyrend;
	    for (lc = lyrstr; lc <= i__3; ++lc) {
		dtau = dtaucp[lc];
		exp1 = exp((utaupr[lu] - taucpr[lc - 1]) / umu[iu]);
		exp2 = exp((utaupr[lu] - taucpr[lc]) / umu[iu]);
		if (*plank && *mazim == 0) {
/*                          ** eqs. stwl(36b,c, 37b,c) */

		    f0n = sgn * (exp1 - exp2);
		    f1n = sgn * ((taucpr[lc - 1] + umu[iu]) * exp1 - (taucpr[
			    lc] + umu[iu]) * exp2);
		    plkint = plkint + z0u[iu + lc * z0u_dim1] * f0n + z1u[iu 
			    + lc * z1u_dim1] * f1n;
		}
		if (*fbeam > 0.f) {
		    denom = umu[iu] / *umu0 + 1.f;
		    if (dabs(denom) < 1e-4f) {
/*                                                   ** l'hospital limit */
			expn = dtau / *umu0 * exp0;
		    } else {
			expn = (exp1 * expbea[lc - 1] - exp2 * expbea[lc]) * 
				sgn / denom;
		    }
		    palint += zbeam[iu + lc * zbeam_dim1] * expn;
		}
/*                                                   ** kk is negative */
		i__4 = *nn;
		for (iq = 1; iq <= i__4; ++iq) {
		    wk[iq] = exp(kk[iq + lc * kk_dim1] * dtau);
		    denom = umu[iu] * kk[iq + lc * kk_dim1] + 1.f;
		    if (dabs(denom) < 1e-4f) {
/*                                                   ** l'hospital limit */
			expn = dtau / umu[iu] * exp2;
		    } else {
			expn = sgn * (exp1 * wk[iq] - exp2) / denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * expn;
/* L40: */
		}
/*                                                   ** kk is positive */
		i__4 = *nstr;
		for (iq = *nn + 1; iq <= i__4; ++iq) {
		    denom = umu[iu] * kk[iq + lc * kk_dim1] + 1.f;
		    if (dabs(denom) < 1e-4f) {
/*                                                   ** l'hospital limit */
			expn = -dtau / umu[iu] * exp1;
		    } else {
			expn = sgn * (exp1 - exp2 * wk[*nstr + 1 - iq]) / 
				denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * expn;
/* L50: */
		}
/* L60: */
	    }
/*                           ** calculate contribution from user */
/*                           ** output level to next computational level */
	    dtau1 = utaupr[lu] - taucpr[lyu - 1];
	    dtau2 = utaupr[lu] - taucpr[lyu];
	    if (dabs(dtau1) < 1e-6f && negumu) {
		goto L90;
	    }
	    if (dabs(dtau2) < 1e-6f && ! negumu) {
		goto L90;
	    }
	    if (negumu) {
		exp1 = exp(dtau1 / umu[iu]);
	    }
	    if (! negumu) {
		exp2 = exp(dtau2 / umu[iu]);
	    }
	    if (*fbeam > 0.f) {
		denom = umu[iu] / *umu0 + 1.f;
		if (dabs(denom) < 1e-4f) {
		    expn = dtau1 / *umu0 * exp0;
		} else if (negumu) {
		    expn = (exp0 - expbea[lyu - 1] * exp1) / denom;
		} else {
		    expn = (exp0 - expbea[lyu] * exp2) / denom;
		}
		palint += zbeam[iu + lyu * zbeam_dim1] * expn;
	    }
/*                                                   ** kk is negative */
	    dtau = dtaucp[lyu];
	    i__3 = *nn;
	    for (iq = 1; iq <= i__3; ++iq) {
		denom = umu[iu] * kk[iq + lyu * kk_dim1] + 1.f;
		if (dabs(denom) < 1e-4f) {
		    expn = -dtau2 / umu[iu] * exp2;
		} else if (negumu) {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau2) - exp(kk[iq 
			    + lyu * kk_dim1] * dtau) * exp1) / denom;
		} else {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau2) - exp2) / 
			    denom;
		}
		palint += gu[iu + (iq + lyu * gu_dim2) * gu_dim1] * expn;
/* L70: */
	    }
/*                                                   ** kk is positive */
	    i__3 = *nstr;
	    for (iq = *nn + 1; iq <= i__3; ++iq) {
		denom = umu[iu] * kk[iq + lyu * kk_dim1] + 1.f;
		if (dabs(denom) < 1e-4f) {
		    expn = -dtau1 / umu[iu] * exp1;
		} else if (negumu) {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau1) - exp1) / 
			    denom;
		} else {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau1) - exp(-kk[iq 
			    + lyu * kk_dim1] * dtau) * exp2) / denom;
		}
		palint += gu[iu + (iq + lyu * gu_dim2) * gu_dim1] * expn;
/* L80: */
	    }
	    if (*plank && *mazim == 0) {
/*                            ** eqs. stwl (35-37) with tau-sub-n-1 */
/*                            ** replaced by tau for upward, and */
/*                            ** tau-sub-n replaced by tau for downward */
/*                            ** directions */
		if (negumu) {
		    expn = exp1;
		    fact = taucpr[lyu - 1] + umu[iu];
		} else {
		    expn = exp2;
		    fact = taucpr[lyu] + umu[iu];
		}
		f0n = 1.f - expn;
		f1n = utaupr[lu] + umu[iu] - fact * expn;
		plkint = plkint + z0u[iu + lyu * z0u_dim1] * f0n + z1u[iu + 
			lyu * z1u_dim1] * f1n;
	    }
/*                            ** calculate intensity components */
/*                            ** attenuated at both boundaries. */
/*                            ** note: no azimuthal intensity */
/*                            ** component for isotropic surface */
L90:
	    bndint = 0.f;
	    if (negumu && *mazim == 0) {
		bndint = (*fisot + *tplank) * exp(utaupr[lu] / umu[iu]);
	    } else if (! negumu) {
		if (*lyrcut || *lamber && *mazim > 0) {
		    goto L140;
		}
		i__3 = *nstr;
		for (jq = *nn + 1; jq <= i__3; ++jq) {
		    wk[jq] = exp(-kk[jq + *nlyr * kk_dim1] * dtaucp[*nlyr]);
/* L100: */
		}
		bnddfu = 0.f;
		for (iq = *nn; iq >= 1; --iq) {
		    dfuint = 0.f;
		    i__3 = *nn;
		    for (jq = 1; jq <= i__3; ++jq) {
			dfuint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * 
				ll[jq + *nlyr * ll_dim1];
/* L110: */
		    }
		    i__3 = *nstr;
		    for (jq = *nn + 1; jq <= i__3; ++jq) {
			dfuint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * 
				ll[jq + *nlyr * ll_dim1] * wk[jq];
/* L120: */
		    }
		    if (*fbeam > 0.f) {
			dfuint += zz[iq + *nlyr * zz_dim1] * expbea[*nlyr];
		    }
		    dfuint += *delm0 * (zplk0[iq + *nlyr * zplk0_dim1] + 
			    zplk1[iq + *nlyr * zplk1_dim1] * taucpr[*nlyr]);
		    bnddfu += (*delm0 + 1.f) * rmu[iu + (*nn + 1 - iq) * 
			    rmu_dim1] * cmu[*nn + 1 - iq] * cwt[*nn + 1 - iq] 
			    * dfuint;
/* L130: */
		}
		bnddir = 0.f;
		if (*fbeam > 0.f) {
		    bnddir = *umu0 * *fbeam / *pi * rmu[iu] * expbea[*nlyr];
		}
		bndint = (bnddfu + bnddir + *delm0 * emu[iu] * *bplank) * exp(
			(utaupr[lu] - taucpr[*nlyr]) / umu[iu]);
	    }
L140:
	    uum[iu + lu * uum_dim1] = palint + plkint + bndint;
L150:
	    ;
	}
/* L160: */
    }
    return 0;
} /* usrint_ */

doublereal xifunc_(real *umu1, real *umu2, real *umu3, real *tau)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static real x1, x2, exp1;

/*          calculates xi function of eq. stwl (72) */
/*                    i n p u t   v a r i a b l e s */
/*        tau         optical thickness of the layer */

/*        umu1,2,3    cosine of zenith angle_1, _2, _3 */

/*   called by- secsca */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    x1 = 1.f / *umu1 - 1.f / *umu2;
    x2 = 1.f / *umu1 - 1.f / *umu3;
    exp1 = exp(-(*tau) / *umu1);
    if (*umu2 == *umu3 && *umu1 == *umu2) {
	ret_val = *tau * *tau * exp1 / (*umu1 * 2.f * *umu2);
    } else if (*umu2 == *umu3 && *umu1 != *umu2) {
	ret_val = ((*tau - 1.f / x1) * exp(-(*tau) / *umu2) + exp1 / x1) / (
		x1 * *umu1 * *umu2);
    } else if (*umu2 != *umu3 && *umu1 == *umu2) {
	ret_val = ((exp(-(*tau) / *umu3) - exp1) / x2 - *tau * exp1) / (x2 * *
		umu1 * *umu2);
    } else if (*umu2 != *umu3 && *umu1 == *umu3) {
	ret_val = ((exp(-(*tau) / *umu2) - exp1) / x1 - *tau * exp1) / (x1 * *
		umu1 * *umu2);
    } else {
	ret_val = ((exp(-(*tau) / *umu3) - exp1) / x2 - (exp(-(*tau) / *umu2) 
		- exp1) / x1) / (x2 * *umu1 * *umu2);
    }
    return ret_val;
} /* xifunc_ */

/* ****************************************************************** */
/* ********** disort service routines ************************ */
/* ****************************************************************** */
/* Subroutine */ int chekin_(integer *nlyr, real *dtauc, real *ssalb, integer 
	*nmom, real *pmom, real *temper, real *wvnm, logical *usrtau, integer 
	*ntau, integer *iref, real *utau, integer *nstr, logical *usrang, 
	integer *numu, real *umu, integer *nphi, real *phi, integer *ibcnd, 
	real *fbeam, real *umu0, real *phi0, real *fisot, logical *lamber, 
	real *albedo, real *btemp, real *ttemp, real *temis, real *surf_pr__, 
	logical *plank, logical *onlyfl, logical *deltam, logical *corint, 
	real *accur, real *tauc, integer *maxcly, integer *maxulv, integer *
	maxumu, integer *maxphi, integer *maxmom, integer *mxcly, integer *
	mxulv, integer *mxumu, integer *mxcmu, integer *mxphi, integer *mxsqt)
{
    /* System generated locals */
    integer pmom_dim1, pmom_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer j, k, lc;
    static logical tf;
    static integer iu, lu;
    static real rmu;
    extern doublereal dref_(real *, real *, integer *);
    static integer irmu;
    static real flxalb;
    extern logical wrtbad_(char *, ftnlen);
    static logical inperr;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    extern logical wrtdim_(char *, integer *, ftnlen);
    static real yessct;
    static integer numsqt;

    /* Fortran I/O blocks */
    static cilist io___318 = { 0, 6, 0, 0, 0 };
    static cilist io___319 = { 0, 6, 0, 0, 0 };
    static cilist io___320 = { 0, 6, 0, 0, 0 };
    static cilist io___321 = { 0, 6, 0, 0, 0 };
    static cilist io___322 = { 0, 6, 0, 0, 0 };
    static cilist io___323 = { 0, 6, 0, 0, 0 };
    static cilist io___324 = { 0, 6, 0, 0, 0 };
    static cilist io___325 = { 0, 6, 0, 0, 0 };
    static cilist io___326 = { 0, 6, 0, 0, 0 };
    static cilist io___327 = { 0, 6, 0, 0, 0 };
    static cilist io___328 = { 0, 6, 0, 0, 0 };
    static cilist io___330 = { 0, 6, 0, 0, 0 };


/*           checks the input dimensions and variables */
/*   calls- wrtbad, wrtdim, dref, errmsg */
/*   called by- disort */
/* -------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*      integer   irmu */
/*      real      flxalb, rmu */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --surf_pr__;
    --ssalb;
    --dtauc;
    --utau;
    --umu;
    --phi;
    pmom_dim1 = *maxmom - 0 + 1;
    pmom_offset = 0 + pmom_dim1;
    pmom -= pmom_offset;

    /* Function Body */
    inperr = FALSE_;
    if (*nstr < 2 || *nstr % 2 != 0) {
	inperr = wrtbad_("nstr", (ftnlen)4);
    }
    tf = TRUE_;
    if (*nstr == 2) {
	errmsg_("chekin--2 streams not recommended; use specialized 2-stream"
		" code twostr instead", &tf, (ftnlen)79);
    }
    if (*nlyr < 1) {
	inperr = wrtbad_("nlyr", (ftnlen)4);
    }
    if (*nlyr > *maxcly) {
	inperr = wrtbad_("maxcly", (ftnlen)6);
    }
    yessct = 0.f;
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	if (dtauc[lc] < 0.f) {
	    inperr = wrtbad_("dtauc", (ftnlen)5);
	}
	if (ssalb[lc] < 0.f || ssalb[lc] > 1.f) {
	    inperr = wrtbad_("ssalb", (ftnlen)5);
	}
	yessct += ssalb[lc];
	if (*plank && *ibcnd != 1) {
	    if (lc == 1 && temper[0] < 0.f) {
		inperr = wrtbad_("temper", (ftnlen)6);
	    }
	    if (temper[lc] < 0.f) {
		inperr = wrtbad_("temper", (ftnlen)6);
	    }
	}
/* L10: */
    }
    if (*nmom < 0 || yessct > 0.f && *nmom < *nstr) {
	inperr = wrtbad_("nmom", (ftnlen)4);
    }
    if (*maxmom < *nmom) {
	inperr = wrtbad_("maxmom", (ftnlen)6);
    }
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	i__2 = *nmom;
	for (k = 0; k <= i__2; ++k) {
	    if (pmom[k + lc * pmom_dim1] < -1.f || pmom[k + lc * pmom_dim1] > 
		    1.f) {
		inperr = wrtbad_("pmom", (ftnlen)4);
	    }
/* L20: */
	}
/* L30: */
    }
    if (*ibcnd == 1) {
	if (*maxulv < 2) {
	    inperr = wrtbad_("maxulv", (ftnlen)6);
	}
    } else if (*usrtau) {
	if (*ntau < 1) {
	    inperr = wrtbad_("ntau", (ftnlen)4);
	}
	if (*maxulv < *ntau) {
	    inperr = wrtbad_("maxulv", (ftnlen)6);
	}
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    if ((r__1 = utau[lu] - tauc[*nlyr], dabs(r__1)) <= 1e-4f) {
		utau[lu] = tauc[*nlyr];
	    }
	    if (utau[lu] < 0.f || utau[lu] > tauc[*nlyr]) {
		inperr = wrtbad_("utau", (ftnlen)4);
	    }
/* L40: */
	}
    } else {
	if (*maxulv < *nlyr + 1) {
	    inperr = wrtbad_("maxulv", (ftnlen)6);
	}
    }
    if (*usrang) {
	if (*numu < 0) {
	    inperr = wrtbad_("numu", (ftnlen)4);
	}
	if (! (*onlyfl) && *numu == 0) {
	    inperr = wrtbad_("numu", (ftnlen)4);
	}
	if (*numu > *maxumu) {
	    inperr = wrtbad_("maxumu", (ftnlen)6);
	}
	if (*ibcnd == 1 && *numu << 1 > *maxumu) {
	    inperr = wrtbad_("maxumu", (ftnlen)6);
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    if (umu[iu] < -1.f || umu[iu] > 1.f || umu[iu] == 0.f) {
		inperr = wrtbad_("umu", (ftnlen)3);
	    }
	    if (*ibcnd == 1 && umu[iu] < 0.f) {
		inperr = wrtbad_("umu", (ftnlen)3);
	    }
	    if (iu > 1) {
		if (umu[iu] < umu[iu - 1]) {
		    inperr = wrtbad_("umu", (ftnlen)3);
		}
	    }
/* L50: */
	}
    } else {
	if (*maxumu < *nstr) {
	    inperr = wrtbad_("maxumu", (ftnlen)6);
	}
    }
    if (! (*onlyfl) && *ibcnd != 1) {
	if (*nphi <= 0) {
	    inperr = wrtbad_("nphi", (ftnlen)4);
	}
	if (*nphi > *maxphi) {
	    inperr = wrtbad_("maxphi", (ftnlen)6);
	}
	i__1 = *nphi;
	for (j = 1; j <= i__1; ++j) {
	    if (phi[j] < 0.f || phi[j] > 360.f) {
		inperr = wrtbad_("phi", (ftnlen)3);
	    }
/* L60: */
	}
    }
    if (*ibcnd < 0 || *ibcnd > 1) {
	inperr = wrtbad_("ibcnd", (ftnlen)5);
    }
    if (*ibcnd == 0) {
	if (*fbeam < 0.f) {
	    inperr = wrtbad_("fbeam", (ftnlen)5);
	}
	if (*fbeam > 0.f && (*umu0 <= 0.f || *umu0 > 1.f)) {
	    inperr = wrtbad_("umu0", (ftnlen)4);
	}
	if (*fbeam > 0.f && (*phi0 < 0.f || *phi0 > 360.f)) {
	    inperr = wrtbad_("phi0", (ftnlen)4);
	}
	if (*fisot < 0.f) {
	    inperr = wrtbad_("fisot", (ftnlen)5);
	}
	if (*lamber) {
	    if (*albedo < 0.f || *albedo > 1.f) {
		inperr = wrtbad_("albedo", (ftnlen)6);
	    }
	} else {
/*                    ** make sure flux albedo at dense mesh of incident */
/*                    ** angles does not assume unphysical values */

/*  ********note, the following check was commented out because it */
/*                dramatically increases the amount of time needed */
/*                in cases with a non-lambertian surface */
/*            do 70 irmu = 3, 100 */
	    for (irmu = 1; irmu <= 10; ++irmu) {

		rmu = irmu * .1f;
		flxalb = dref_(&rmu, &surf_pr__[1], iref);

		if (flxalb < 0.f || flxalb > 1.f) {
		    inperr = wrtbad_("function bdref", (ftnlen)14);
		}
		if (inperr) {
		    s_wsle(&io___318);
		    do_lio(&c__9, &c__1, "at rmu=", (ftnlen)7);
		    do_lio(&c__4, &c__1, (char *)&rmu, (ftnlen)sizeof(real));
		    do_lio(&c__9, &c__1, "function bdref gives flxalb=", (
			    ftnlen)28);
		    do_lio(&c__4, &c__1, (char *)&flxalb, (ftnlen)sizeof(real)
			    );
		    e_wsle();
		}
/* L70: */
	    }
	}
    } else if (*ibcnd == 1) {
	if (*albedo < 0.f || *albedo > 1.f) {
	    inperr = wrtbad_("albedo", (ftnlen)6);
	}
	if (inperr) {
	    s_wsle(&io___319);
	    do_lio(&c__9, &c__1, "albedo=", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&(*albedo), (ftnlen)sizeof(real));
	    e_wsle();
	}
    }
    if (*plank && *ibcnd != 1) {
	if (*wvnm < 0.f) {
	    inperr = wrtbad_("wvnm", (ftnlen)4);
	}
	if (inperr) {
	    s_wsle(&io___320);
	    do_lio(&c__9, &c__1, "wvnm", (ftnlen)4);
	    do_lio(&c__4, &c__1, (char *)&(*wvnm), (ftnlen)sizeof(real));
	    e_wsle();
	}
	if (*temis < 0.f || *temis > 1.f) {
	    inperr = wrtbad_("temis", (ftnlen)5);
	}
	if (inperr) {
	    s_wsle(&io___321);
	    do_lio(&c__9, &c__1, "temis =", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&(*temis), (ftnlen)sizeof(real));
	    e_wsle();
	}
	if (*btemp < 0.f) {
	    inperr = wrtbad_("btemp", (ftnlen)5);
	}
	if (inperr) {
	    s_wsle(&io___322);
	    do_lio(&c__9, &c__1, "btemp =", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&(*btemp), (ftnlen)sizeof(real));
	    e_wsle();
	}
	if (*ttemp < 0.f) {
	    inperr = wrtbad_("ttemp", (ftnlen)5);
	}
	if (inperr) {
	    s_wsle(&io___323);
	    do_lio(&c__9, &c__1, "ttemp =", (ftnlen)7);
	    do_lio(&c__4, &c__1, (char *)&(*ttemp), (ftnlen)sizeof(real));
	    e_wsle();
	}
    }
    if (*accur < 0.f || *accur > .01f) {
	inperr = wrtbad_("accur", (ftnlen)5);
    }
    if (inperr) {
	s_wsle(&io___324);
	do_lio(&c__9, &c__1, "accur =", (ftnlen)7);
	do_lio(&c__4, &c__1, (char *)&(*accur), (ftnlen)sizeof(real));
	e_wsle();
    }
    if (*mxcly < *nlyr) {
	inperr = wrtdim_("mxcly", nlyr, (ftnlen)5);
    }
    if (*ibcnd != 1) {
	if (*usrtau && *mxulv < *ntau) {
	    inperr = wrtdim_("mxulv", ntau, (ftnlen)5);
	}
	if (! (*usrtau) && *mxulv < *nlyr + 1) {
	    i__1 = *nlyr + 1;
	    inperr = wrtdim_("mxulv", &i__1, (ftnlen)5);
	}
    } else {
	if (*mxulv < 2) {
	    inperr = wrtdim_("mxulv", &c__2, (ftnlen)5);
	}
    }
    if (*mxcmu < *nstr) {
	inperr = wrtdim_("mxcmu", nstr, (ftnlen)5);
    }
    if (inperr) {
	s_wsle(&io___325);
	do_lio(&c__9, &c__1, "mxcmu,nstr", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&(*mxcmu), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (*usrang && *mxumu < *numu) {
	inperr = wrtdim_("mxumu", numu, (ftnlen)5);
    }
    if (inperr) {
	s_wsle(&io___326);
	do_lio(&c__9, &c__1, "mxumu, numu", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*mxumu), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*numu), (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (*usrang && *ibcnd == 1 && *mxumu < *numu << 1) {
	i__1 = *numu << 1;
	inperr = wrtdim_("mxumu", &i__1, (ftnlen)5);
    }
    if (! (*usrang) && *mxumu < *nstr) {
	inperr = wrtdim_("mxumu", nstr, (ftnlen)5);
    }
    if (inperr) {
	s_wsle(&io___327);
	do_lio(&c__9, &c__1, "mxumu, nstr", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*mxumu), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (! (*onlyfl) && *ibcnd != 1 && *mxphi < *nphi) {
	inperr = wrtdim_("mxphi", nphi, (ftnlen)5);
    }
    if (inperr) {
	s_wsle(&io___328);
	do_lio(&c__9, &c__1, "mxphi, nphi", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*mxphi), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
	e_wsle();
    }
    numsqt = max(100,*nstr) << 1;
    if (*mxsqt < numsqt) {
	inperr = wrtdim_("mxsqt", &numsqt, (ftnlen)5);
    }
    if (inperr) {
	s_wsle(&io___330);
	do_lio(&c__9, &c__1, "mxsqt,numsqt", (ftnlen)12);
	do_lio(&c__3, &c__1, (char *)&(*mxsqt), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&numsqt, (ftnlen)sizeof(integer));
	e_wsle();
    }
    tf = TRUE_;
    if (inperr) {
	errmsg_("disort--input and/or dimension errors", &tf, (ftnlen)37);
    }
    if (*plank) {
	tf = FALSE_;
	i__1 = *nlyr;
	for (lc = 1; lc <= i__1; ++lc) {
/*            if( abs( temper( lc )-temper( lc-1 ) ).gt. 10.0 ) */
	    if ((r__1 = temper[lc] - temper[lc - 1], dabs(r__1)) > 50.f) {
		errmsg_("chekin--vertical temperature step may be too large "
			"for good accuracy", &tf, (ftnlen)68);
	    }
/* L80: */
	}
    }
    tf = FALSE_;
    if (! (*corint) && ! (*onlyfl) && *fbeam > 0.f && yessct > 0.f && *deltam)
	     {
	errmsg_("chekin--intensity correction is off; intensities may be les"
		"s accurate", &tf, (ftnlen)69);
    }
    return 0;
} /* chekin_ */

doublereal dref_(real *mu, real *surf_pr__, integer *iref)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double asin(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer k, jg;
    static real pi;
    static logical tf;
    static real gmu[50], gwt[50], sum;
    extern doublereal bdref_(real *, real *, real *, real *, integer *);
    static real pigmu;
    extern /* Subroutine */ int qgausn_(integer *, real *, real *), errmsg_(
	    char *, logical *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___340 = { 0, 6, 0, "(/,1a,1pe12.4)", 0 };


/*        flux albedo for given angle of incidence, given */
/*        a bidirectional reflectivity. */

/*  input :   mu      cosine of incidence angle */

/*            wvnm  wavenumber (inv-cm) of spectral interval */

/*            iref  : bidirectional reflectance options */
/*             0 - Lambert */
/*             1 - Hapke's BDR model */
/*             2 - Breon's BDR model; combination of Li + Roujean */
/*             3 - Roujean's BDR model */
/*             4 - Cox and Munk glint model */

/*   SURF_PR : Wavelength dependent surface properties array */
/*             IREF= 0 - Lambert albedo */
/*             IREF= 1 - Hapke : HH, W */
/*             IREF= 2 - Breon's BDR model: k0, k1, k2 */
/*             IREF= 3 - Roujean's BDR model: k0, k1, k2 */
/*             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw */

/*  internal variables : */

/*       nmug   :  number of angle cosine quadrature points on (-1,1) */
/*                 for integrating bidirectional reflectivity to get */
/*                 directional emissivity (it is necessary to use a */
/*                 quadrature set distinct from the computational angles, */
/*                 because the computational angles may not be dense */
/*                 enough -- i.e. 'nstr' may be too small -- to give an */
/*                 accurate approximation for the integration). */

/*       gmu    :  the 'nmug' angle cosine quadrature points on (0,1) */

/*       gwt    :  the 'nmug' angle cosine quadrature weights on (0,1) */

/*   called by- chekin */
/*   calls- qgausn, errmsg, bdref */
/* +--------------------------------------------------------------------+ */
/*     .. parameters .. */
/*     .. */
/*     .. scalar arguments .. */

/*     .. array arguments */

/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. local arrays .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --surf_pr__;

    /* Function Body */
    if (pass1) {
	pass1 = FALSE_;
	pi = asin(1.f) * 2.f;
	qgausn_(&c__25, gmu, gwt);
	for (k = 1; k <= 25; ++k) {
	    gmu[k + 24] = -gmu[k - 1];
	    gwt[k + 24] = gwt[k - 1];
/* L10: */
	}
    }
    tf = TRUE_;
    if (dabs(*mu) > 1.f) {
	errmsg_("dref--input argument error(s)", &tf, (ftnlen)29);
    }
    ret_val = 0.f;
/*                       ** loop over azimuth angle difference */
    for (jg = 1; jg <= 50; ++jg) {
	sum = 0.f;
/*                       ** loop over angle of reflection */
	for (k = 1; k <= 25; ++k) {
	    pigmu = pi * gmu[jg - 1];
	    sum += gwt[k - 1] * gmu[k - 1] * bdref_(&gmu[k - 1], mu, &pigmu, &
		    surf_pr__[1], iref);
/* L20: */
	}
	ret_val += gwt[jg - 1] * sum;
/* L30: */
    }
    if (ret_val < 0.f || ret_val > 1.f) {
	tf = FALSE_;
	s_wsfe(&io___340);
	do_fio(&c__1, " dref = ", (ftnlen)8);
	do_fio(&c__1, (char *)&ret_val, (ftnlen)sizeof(real));
	e_wsfe();
	errmsg_("dref--albedo value not in (0,1)", &tf, (ftnlen)31);
    }
    return ret_val;
} /* dref_ */

/* Subroutine */ int lepoly_(integer *nmu, integer *m, integer *maxmu, 
	integer *twonm1, real *mu, real *sqt, real *ylm)
{
    /* System generated locals */
    integer ylm_dim1, ylm_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, l;
    static real tmp1, tmp2;

/*       computes the normalized associated legendre polynomial, */
/*       defined in terms of the associated legendre polynomial */
/*       plm = p-sub-l-super-m as */

/*             ylm(mu) = sqrt( (l-m)!/(l+m)! ) * plm(mu) */

/*       for fixed order m and all degrees from l = m to twonm1. */
/*       when m.gt.0, assumes that y-sub(m-1)-super(m-1) is available */
/*       from a prior call to the routine. */

/*       reference: dave, j.v. and b.h. armstrong, computations of */
/*                  high-order associated legendre polynomials, */
/*                  j. quant. spectrosc. radiat. transfer 10, */
/*                  557-562, 1970.  (hereafter d/a) */

/*       method: varying degree recurrence relationship. */

/*       notes: */
/*       (1) the d/a formulas are transformed by setting m=n-1; l=k-1. */
/*       (2) assumes that routine is called first with  m = 0, then with */
/*           m = 1, etc. up to  m = twonm1. */


/*  i n p u t     v a r i a b l e s: */

/*       nmu    :  number of arguments of ylm */

/*       m      :  order of ylm */

/*       maxmu  :  first dimension of ylm */

/*       twonm1 :  max degree of ylm */

/*       mu(i)  :  arguments of ylm (i = 1 to nmu) */

/*       sqt(k) :  square root of k */

/*       if m.gt.0, ylm(m-1,i) for i = 1 to nmu is assumed to exist */
/*       from a prior call. */


/*  o u t p u t     v a r i a b l e: */

/*       ylm(l,i) :  l = m to twonm1, normalized associated legendre */
/*                   polynomials evaluated at argument mu(i) */

/*   called by- disort, albtrn */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
    /* Parameter adjustments */
    ylm_dim1 = *maxmu - 0 + 1;
    ylm_offset = 0 + ylm_dim1;
    ylm -= ylm_offset;
    --mu;
    --sqt;

    /* Function Body */
    if (*m == 0) {
/*                             ** upward recurrence for ordinary */
/*                             ** legendre polynomials */
	i__1 = *nmu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ylm[i__ * ylm_dim1] = 1.f;
	    ylm[i__ * ylm_dim1 + 1] = mu[i__];
/* L20: */
	}
	i__1 = *twonm1;
	for (l = 2; l <= i__1; ++l) {
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1] = (((l << 1) - 1) * mu[i__] * ylm[l - 
			1 + i__ * ylm_dim1] - (l - 1) * ylm[l - 2 + i__ * 
			ylm_dim1]) / l;
/* L30: */
	    }
/* L40: */
	}
    } else {
	i__1 = *nmu;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*                               ** y-sub-m-super-m; derived from */
/*                               ** d/a eqs. (11,12), stwl(58c) */
/* Computing 2nd power */
	    r__1 = mu[i__];
	    ylm[*m + i__ * ylm_dim1] = -sqt[(*m << 1) - 1] / sqt[*m * 2] * 
		    sqrt(1.f - r__1 * r__1) * ylm[*m - 1 + i__ * ylm_dim1];
/*                              ** y-sub-(m+1)-super-m; derived from */
/*                              ** d/a eqs.(13,14) using eqs.(11,12), */
/*                              ** stwl(58f) */
	    ylm[*m + 1 + i__ * ylm_dim1] = sqt[(*m << 1) + 1] * mu[i__] * ylm[
		    *m + i__ * ylm_dim1];
/* L50: */
	}
/*                                   ** upward recurrence; d/a eq.(10), */
/*                                   ** stwl(58a) */
	i__1 = *twonm1;
	for (l = *m + 2; l <= i__1; ++l) {
	    tmp1 = sqt[l - *m] * sqt[l + *m];
	    tmp2 = sqt[l - *m - 1] * sqt[l + *m - 1];
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1] = (((l << 1) - 1) * mu[i__] * ylm[l - 
			1 + i__ * ylm_dim1] - tmp2 * ylm[l - 2 + i__ * 
			ylm_dim1]) / tmp1;
/* L60: */
	    }
/* L70: */
	}
    }
    return 0;
} /* lepoly_ */

/* Subroutine */ int pravin_(real *umu, integer *numu, integer *mxumu, real *
	utau, integer *ntau, real *u0u)
{
    /* System generated locals */
    integer u0u_dim1, u0u_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer iu, np, lu, iumin, iumax, npass, lenfmt;

    /* Fortran I/O blocks */
    static cilist io___345 = { 0, 6, 0, "(//,a)", 0 };
    static cilist io___348 = { 0, 6, 0, "(/,a,/,a)", 0 };
    static cilist io___352 = { 0, 6, 0, "(/,10x,8f14.5)", 0 };
    static cilist io___355 = { 0, 6, 0, "(0p,f10.4,1p,8e14.4)", 0 };


/*        print azimuthally averaged intensities at user angles */
/*   called by- disort */
/*     lenfmt   max number of polar angle cosines umu that can be */
/*              printed on one line, as set in format statement */
/* -------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --umu;
    u0u_dim1 = *mxumu;
    u0u_offset = 1 + u0u_dim1;
    u0u -= u0u_offset;
    --utau;

    /* Function Body */
    if (*numu < 1) {
	return 0;
    }
    s_wsfe(&io___345);
    do_fio(&c__1, " *******  azimuthally averaged intensities (at user polar"
	    " angles)  ********", (ftnlen)75);
    e_wsfe();
    lenfmt = 8;
    npass = (*numu - 1) / lenfmt + 1;
    s_wsfe(&io___348);
    do_fio(&c__1, "   optical   polar angle cosines", (ftnlen)32);
    do_fio(&c__1, "     depth", (ftnlen)10);
    e_wsfe();
    i__1 = npass;
    for (np = 1; np <= i__1; ++np) {
	iumin = lenfmt * (np - 1) + 1;
/* Computing MIN */
	i__2 = lenfmt * np;
	iumax = min(i__2,*numu);
	s_wsfe(&io___352);
	i__2 = iumax;
	for (iu = iumin; iu <= i__2; ++iu) {
	    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	}
	e_wsfe();
	i__2 = *ntau;
	for (lu = 1; lu <= i__2; ++lu) {
	    s_wsfe(&io___355);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    i__3 = iumax;
	    for (iu = iumin; iu <= i__3; ++iu) {
		do_fio(&c__1, (char *)&u0u[iu + lu * u0u_dim1], (ftnlen)
			sizeof(real));
	    }
	    e_wsfe();
/* L10: */
	}
/* L20: */
    }
    return 0;
} /* pravin_ */

/* Subroutine */ int prtinp_(integer *nlyr, real *dtauc, real *dtaucp, real *
	ssalb, integer *nmom, real *pmom, real *temper, real *wvnm, integer *
	ntau, real *utau, integer *nstr, integer *numu, real *umu, integer *
	nphi, real *phi, integer *ibcnd, real *fbeam, real *umu0, real *phi0, 
	real *fisot, logical *lamber, real *albedo, real *btemp, real *ttemp, 
	real *temis, logical *deltam, logical *plank, logical *onlyfl, 
	logical *corint, real *accur, real *flyr, logical *lyrcut, real *
	oprim, real *tauc, real *taucpr, integer *maxmom, logical *prtmom)
{
    /* System generated locals */
    integer pmom_dim1, pmom_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer j, k, lc, iu, lu;
    static real yessct;

    /* Fortran I/O blocks */
    static cilist io___356 = { 0, 6, 0, "(/,a,i4,a,i4)", 0 };
    static cilist io___357 = { 0, 6, 0, "(i4,a,10f10.4,/,(26x,10f10.4))", 0 };
    static cilist io___359 = { 0, 6, 0, "(i4,a,10f9.5,/,(31x,10f9.5))", 0 };
    static cilist io___361 = { 0, 6, 0, "(i4,a,10f9.2,/,(28x,10f9.2))", 0 };
    static cilist io___363 = { 0, 6, 0, "(a)", 0 };
    static cilist io___364 = { 0, 6, 0, "(a,i2)", 0 };
    static cilist io___365 = { 0, 6, 0, "(a,1p,e11.3,a,0p,f8.5,a,f7.2,/,a,1p"
	    ",e11.3)", 0 };
    static cilist io___366 = { 0, 6, 0, "(a,0p,f8.4)", 0 };
    static cilist io___367 = { 0, 6, 0, "(a)", 0 };
    static cilist io___368 = { 0, 6, 0, "(a,f14.4,/,a,f10.2,a,f10.2,a,f8.4)", 
	    0 };
    static cilist io___369 = { 0, 6, 0, "(a)", 0 };
    static cilist io___370 = { 0, 6, 0, "(a,0p,f8.4)", 0 };
    static cilist io___371 = { 0, 6, 0, "(a)", 0 };
    static cilist io___372 = { 0, 6, 0, "(a)", 0 };
    static cilist io___373 = { 0, 6, 0, "(a)", 0 };
    static cilist io___374 = { 0, 6, 0, "(a)", 0 };
    static cilist io___375 = { 0, 6, 0, "(a)", 0 };
    static cilist io___376 = { 0, 6, 0, "(a)", 0 };
    static cilist io___377 = { 0, 6, 0, "(a)", 0 };
    static cilist io___378 = { 0, 6, 0, "(a,1p,e11.2)", 0 };
    static cilist io___379 = { 0, 6, 0, "(a)", 0 };
    static cilist io___380 = { 0, 6, 0, "(/,37x,a,3(/,2a))", 0 };
    static cilist io___381 = { 0, 6, 0, "(/,37x,a,3(/,2a))", 0 };
    static cilist io___384 = { 0, 6, 0, "(i4,2f10.4,f10.5,f12.5,2f10.4,f10.5"
	    ",f9.4,f14.3)", 0 };
    static cilist io___385 = { 0, 6, 0, "(i4,2f10.4,f10.5,f12.5,2f10.4,f10.5"
	    ",f9.4)", 0 };
    static cilist io___386 = { 0, 6, 0, "(85x,f14.3)", 0 };
    static cilist io___387 = { 0, 6, 0, "(/,a,i5)", 0 };
    static cilist io___388 = { 0, 6, 0, "(a)", 0 };
    static cilist io___389 = { 0, 6, 0, "(i6,10f11.6,/,(6x,10f11.6))", 0 };


/*        print values of input variables */
/*   called by- disort */
/* -------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
    /* Parameter adjustments */
    --dtauc;
    --dtaucp;
    --ssalb;
    --utau;
    --umu;
    --phi;
    --flyr;
    --oprim;
    pmom_dim1 = *maxmom - 0 + 1;
    pmom_offset = 0 + pmom_dim1;
    pmom -= pmom_offset;

    /* Function Body */
    s_wsfe(&io___356);
    do_fio(&c__1, " no. streams =", (ftnlen)14);
    do_fio(&c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
    do_fio(&c__1, "     no. computational layers =", (ftnlen)31);
    do_fio(&c__1, (char *)&(*nlyr), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ibcnd != 1) {
	s_wsfe(&io___357);
	do_fio(&c__1, (char *)&(*ntau), (ftnlen)sizeof(integer));
	do_fio(&c__1, " user optical depths :", (ftnlen)22);
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    if (! (*onlyfl)) {
	s_wsfe(&io___359);
	do_fio(&c__1, (char *)&(*numu), (ftnlen)sizeof(integer));
	do_fio(&c__1, " user polar angle cosines :", (ftnlen)27);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    if (! (*onlyfl) && *ibcnd != 1) {
	s_wsfe(&io___361);
	do_fio(&c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
	do_fio(&c__1, " user azimuthal angles :", (ftnlen)24);
	i__1 = *nphi;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&phi[j], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    if (! (*plank) || *ibcnd == 1) {
	s_wsfe(&io___363);
	do_fio(&c__1, " no thermal emission", (ftnlen)20);
	e_wsfe();
    }
    s_wsfe(&io___364);
    do_fio(&c__1, " boundary condition flag: ibcnd =", (ftnlen)33);
    do_fio(&c__1, (char *)&(*ibcnd), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ibcnd == 0) {
	s_wsfe(&io___365);
	do_fio(&c__1, "    incident beam with intensity =", (ftnlen)34);
	do_fio(&c__1, (char *)&(*fbeam), (ftnlen)sizeof(real));
	do_fio(&c__1, " and polar angle cosine = ", (ftnlen)26);
	do_fio(&c__1, (char *)&(*umu0), (ftnlen)sizeof(real));
	do_fio(&c__1, "  and azimuth angle =", (ftnlen)21);
	do_fio(&c__1, (char *)&(*phi0), (ftnlen)sizeof(real));
	do_fio(&c__1, "    plus isotropic incident intensity =", (ftnlen)39);
	do_fio(&c__1, (char *)&(*fisot), (ftnlen)sizeof(real));
	e_wsfe();
	if (*lamber) {
	    s_wsfe(&io___366);
	    do_fio(&c__1, "    bottom albedo (lambertian) =", (ftnlen)32);
	    do_fio(&c__1, (char *)&(*albedo), (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (! (*lamber)) {
	    s_wsfe(&io___367);
	    do_fio(&c__1, "    bidirectional reflectivity at bottom", (ftnlen)
		    40);
	    e_wsfe();
	}
	if (*plank) {
	    s_wsfe(&io___368);
	    do_fio(&c__1, "    thermal emission at wavenumber: ", (ftnlen)36);
	    do_fio(&c__1, (char *)&(*wvnm), (ftnlen)sizeof(real));
	    do_fio(&c__1, "    bottom temperature =", (ftnlen)24);
	    do_fio(&c__1, (char *)&(*btemp), (ftnlen)sizeof(real));
	    do_fio(&c__1, "    top temperature =", (ftnlen)21);
	    do_fio(&c__1, (char *)&(*ttemp), (ftnlen)sizeof(real));
	    do_fio(&c__1, "    top emissivity =", (ftnlen)20);
	    do_fio(&c__1, (char *)&(*temis), (ftnlen)sizeof(real));
	    e_wsfe();
	}
    } else if (*ibcnd == 1) {
	s_wsfe(&io___369);
	do_fio(&c__1, "    isotropic illumination from top and bottom", (
		ftnlen)46);
	e_wsfe();
	s_wsfe(&io___370);
	do_fio(&c__1, "    bottom albedo (lambertian) =", (ftnlen)32);
	do_fio(&c__1, (char *)&(*albedo), (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*deltam) {
	s_wsfe(&io___371);
	do_fio(&c__1, " uses delta-m method", (ftnlen)20);
	e_wsfe();
    }
    if (! (*deltam)) {
	s_wsfe(&io___372);
	do_fio(&c__1, " does not use delta-m method", (ftnlen)28);
	e_wsfe();
    }
    if (*corint) {
	s_wsfe(&io___373);
	do_fio(&c__1, " uses tms/ims method", (ftnlen)20);
	e_wsfe();
    }
    if (! (*corint)) {
	s_wsfe(&io___374);
	do_fio(&c__1, " does not use tms/ims method", (ftnlen)28);
	e_wsfe();
    }
    if (*ibcnd == 1) {
	s_wsfe(&io___375);
	do_fio(&c__1, " calculate albedo and transmissivity of medium vs. in"
		"cident beam angle", (ftnlen)70);
	e_wsfe();
    } else if (*onlyfl) {
	s_wsfe(&io___376);
	do_fio(&c__1, " calculate fluxes only", (ftnlen)22);
	e_wsfe();
    } else {
	s_wsfe(&io___377);
	do_fio(&c__1, " calculate fluxes and intensities", (ftnlen)33);
	e_wsfe();
    }
    s_wsfe(&io___378);
    do_fio(&c__1, " relative convergence criterion for azimuth series =", (
	    ftnlen)52);
    do_fio(&c__1, (char *)&(*accur), (ftnlen)sizeof(real));
    e_wsfe();
    if (*lyrcut) {
	s_wsfe(&io___379);
	do_fio(&c__1, " sets radiation = 0 below absorption optical depth 10",
		 (ftnlen)53);
	e_wsfe();
    }
/*                                    ** print layer variables */
/*                                    ** (to read, skip every other line) */
    if (*plank) {
	s_wsfe(&io___380);
	do_fio(&c__1, "<------------- delta-m --------------->", (ftnlen)39);
	do_fio(&c__1, "                   total    single                   "
		"        ", (ftnlen)61);
	do_fio(&c__1, "total    single", (ftnlen)15);
	do_fio(&c__1, "       optical   optical   scatter   separated   ", (
		ftnlen)49);
	do_fio(&c__1, "optical   optical   scatter    asymm", (ftnlen)36);
	do_fio(&c__1, "         depth     depth    albedo    fraction     ", (
		ftnlen)51);
	do_fio(&c__1, "depth     depth    albedo   factor   temperature", (
		ftnlen)48);
	e_wsfe();
    }
    if (! (*plank)) {
	s_wsfe(&io___381);
	do_fio(&c__1, "<------------- delta-m --------------->", (ftnlen)39);
	do_fio(&c__1, "                   total    single                   "
		"        ", (ftnlen)61);
	do_fio(&c__1, "total    single", (ftnlen)15);
	do_fio(&c__1, "       optical   optical   scatter   separated   ", (
		ftnlen)49);
	do_fio(&c__1, "optical   optical   scatter    asymm", (ftnlen)36);
	do_fio(&c__1, "         depth     depth    albedo    fraction     ", (
		ftnlen)51);
	do_fio(&c__1, "depth     depth    albedo   factor", (ftnlen)34);
	e_wsfe();
    }
    yessct = 0.f;
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	yessct += ssalb[lc];
/*                                       ** f90 nonadvancing i/o would */
/*                                       ** simplify this a lot (also the */
/*                                       ** two writes above) */
	if (*plank) {
	    s_wsfe(&io___384);
	    do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dtauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&tauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ssalb[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&flyr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dtaucp[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&taucpr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&oprim[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&pmom[lc * pmom_dim1 + 1], (ftnlen)sizeof(
		    real));
	    do_fio(&c__1, (char *)&temper[lc - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (! (*plank)) {
	    s_wsfe(&io___385);
	    do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dtauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&tauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ssalb[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&flyr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dtaucp[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&taucpr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&oprim[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&pmom[lc * pmom_dim1 + 1], (ftnlen)sizeof(
		    real));
	    e_wsfe();
	}
/* L10: */
    }
    if (*plank) {
	s_wsfe(&io___386);
	do_fio(&c__1, (char *)&temper[*nlyr], (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*prtmom && yessct > 0.f) {
	s_wsfe(&io___387);
	do_fio(&c__1, " number of phase function moments = ", (ftnlen)36);
	i__1 = *nmom + 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___388);
	do_fio(&c__1, " layer   phase function moments", (ftnlen)31);
	e_wsfe();
	i__1 = *nlyr;
	for (lc = 1; lc <= i__1; ++lc) {
	    if (ssalb[lc] > 0.f) {
		s_wsfe(&io___389);
		do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
		i__2 = *nmom;
		for (k = 0; k <= i__2; ++k) {
		    do_fio(&c__1, (char *)&pmom[k + lc * pmom_dim1], (ftnlen)
			    sizeof(real));
		}
		e_wsfe();
	    }
/* L20: */
	}
    }
    return 0;
} /* prtinp_ */

/* Subroutine */ int prtint_(real *uu, real *utau, integer *ntau, real *umu, 
	integer *numu, real *phi, integer *nphi, integer *maxulv, integer *
	maxumu)
{
    /* System generated locals */
    integer uu_dim1, uu_dim2, uu_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer j, iu, np, lu, jmin, jmax, npass, lenfmt;

    /* Fortran I/O blocks */
    static cilist io___391 = { 0, 6, 0, "(//,a)", 0 };
    static cilist io___394 = { 0, 6, 0, "(/,a,/,a,/,a)", 0 };
    static cilist io___399 = { 0, 6, 0, "(/,18x,10f11.2)", 0 };
    static cilist io___401 = { 0, 6, 0, "(f10.4,f8.4,1p,10e11.3)", 0 };
    static cilist io___402 = { 0, 6, 0, "(10x,f8.4,1p,10e11.3)", 0 };
    static cilist io___404 = { 0, 6, 0, "(10x,f8.4,1p,10e11.3)", 0 };


/*         prints the intensity at user polar and azimuthal angles */
/*     all arguments are disort input or output variables */
/*   called by- disort */
/*     lenfmt   max number of azimuth angles phi that can be printed */
/*                on one line, as set in format statement */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --utau;
    --umu;
    --phi;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2);
    uu -= uu_offset;

    /* Function Body */
    if (*nphi < 1) {
	return 0;
    }
    s_wsfe(&io___391);
    do_fio(&c__1, " *********  i n t e n s i t i e s  *********", (ftnlen)44);
    e_wsfe();
    lenfmt = 10;
    npass = (*nphi - 1) / lenfmt + 1;
    s_wsfe(&io___394);
    do_fio(&c__1, "             polar   azimuth angles (degrees)", (ftnlen)45)
	    ;
    do_fio(&c__1, "   optical   angle", (ftnlen)18);
    do_fio(&c__1, "    depth   cosine", (ftnlen)18);
    e_wsfe();
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	i__2 = npass;
	for (np = 1; np <= i__2; ++np) {
	    jmin = lenfmt * (np - 1) + 1;
/* Computing MIN */
	    i__3 = lenfmt * np;
	    jmax = min(i__3,*nphi);
	    s_wsfe(&io___399);
	    i__3 = jmax;
	    for (j = jmin; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&phi[j], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    if (np == 1) {
		s_wsfe(&io___401);
		do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&umu[1], (ftnlen)sizeof(real));
		i__3 = jmax;
		for (j = jmin; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&uu[(lu + j * uu_dim2) * uu_dim1 + 
			    1], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    }
	    if (np > 1) {
		s_wsfe(&io___402);
		do_fio(&c__1, (char *)&umu[1], (ftnlen)sizeof(real));
		i__3 = jmax;
		for (j = jmin; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&uu[(lu + j * uu_dim2) * uu_dim1 + 
			    1], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    }
	    i__3 = *numu;
	    for (iu = 2; iu <= i__3; ++iu) {
		s_wsfe(&io___404);
		do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
		i__4 = jmax;
		for (j = jmin; j <= i__4; ++j) {
		    do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
			    uu_dim1], (ftnlen)sizeof(real));
		}
		e_wsfe();
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }
    return 0;
} /* prtint_ */

/* Subroutine */ int qgausn_(integer *m, real *gmu, real *gwt)
{
    /* Initialized data */

    static real pi = 0.f;
    static integer maxit = 1000;
    static doublereal one = 1.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double asin(doublereal), tan(doublereal), cos(doublereal);

    /* Local variables */
    static integer k;
    static doublereal p;
    static real t;
    static doublereal x, en;
    static logical tf;
    static integer nn;
    static doublereal xi, pm1;
    static integer np1;
    static doublereal pm2;
    static integer lim;
    static doublereal tol, tmp, ppr, nnp1;
    static real cona;
    static integer iter;
    static doublereal prod, p2pri;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*       compute weights and abscissae for ordinary gaussian quadrature */
/*       on the interval (0,1);  that is, such that */
/*           sum(i=1 to m) ( gwt(i) f(gmu(i)) ) */
/*       is a good approximation to */
/*           integral(0 to 1) ( f(x) dx ) */
/*   input :    m       order of quadrature rule */
/*   output :  gmu(i)   array of abscissae (i = 1 to m) */
/*             gwt(i)   array of weights (i = 1 to m) */
/*   reference:  davis, p.j. and p. rabinowitz, methods of numerical */
/*                   integration, academic press, new york, pp. 87, 1975 */
/*   method:  compute the abscissae as roots of the legendre */
/*            polynomial p-sub-m using a cubically convergent */
/*            refinement of newton's method.  compute the */
/*            weights from eq. 2.7.3.8 of davis/rabinowitz.  note */
/*            that newton's method can very easily diverge; only a */
/*            very good initial guess can guarantee convergence. */
/*            the initial guess used here has never led to divergence */
/*            even for m up to 1000. */
/*   accuracy:  relative error no better than tol or computer */
/*              precision (machine epsilon), whichever is larger */
/*   internal variables: */
/*    iter      : number of newton method iterations */
/*    maxit     : maximum allowed iterations of newton method */
/*    pm2,pm1,p : 3 successive legendre polynomials */
/*    ppr       : derivative of legendre polynomial */
/*    p2pri     : 2nd derivative of legendre polynomial */
/*    tol       : convergence criterion for legendre poly root iteration */
/*    x,xi      : successive iterates in cubically-convergent version */
/*                of newtons method (seeking roots of legendre poly.) */
/*   called by- dref, setdis, surfac */
/*   calls- d1mach, errmsg */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --gwt;
    --gmu;

    /* Function Body */
    if (pi == 0.f) {
	pi = asin(1.f) * 2.f;
	tol = d1mach_(&c__4) * 10.f;
    }
    tf = TRUE_;
    if (*m < 1) {
	errmsg_("qgausn--bad value of m", &tf, (ftnlen)22);
    }
    if (*m == 1) {
	gmu[1] = .5f;
	gwt[1] = 1.f;
	return 0;
    }
    en = (doublereal) (*m);
    np1 = *m + 1;
    nnp1 = (doublereal) (*m * np1);
/* Computing 3rd power */
    i__1 = *m;
    cona = (real) (*m - 1) / (i__1 * (i__1 * i__1) << 3);
    lim = *m / 2;
    i__1 = lim;
    for (k = 1; k <= i__1; ++k) {
/*                                        ** initial guess for k-th root */
/*                                        ** of legendre polynomial, from */
/*                                        ** davis/rabinowitz (2.7.3.3a) */
	t = ((k << 2) - 1) * pi / ((*m << 2) + 2);
	x = cos(t + cona / tan(t));
	iter = 0;
/*                                        ** upward recurrence for */
/*                                        ** legendre polynomials */
L10:
	++iter;
	pm2 = one;
	pm1 = x;
	i__2 = *m;
	for (nn = 2; nn <= i__2; ++nn) {
	    p = (((nn << 1) - 1) * x * pm1 - (nn - 1) * pm2) / nn;
	    pm2 = pm1;
	    pm1 = p;
/* L20: */
	}
/*                                              ** newton method */
/* Computing 2nd power */
	d__1 = x;
	tmp = one / (one - d__1 * d__1);
	ppr = en * (pm2 - x * p) * tmp;
	p2pri = (two * x * ppr - nnp1 * p) * tmp;
	xi = x - p / ppr * (one + p / ppr * p2pri / (two * ppr));
/*                                              ** check for convergence */
	if ((d__1 = xi - x, abs(d__1)) > tol) {
	    tf = TRUE_;
	    if (iter > maxit) {
		errmsg_("qgausn--max iteration count", &tf, (ftnlen)27);
	    }
	    x = xi;
	    goto L10;
	}
/*                             ** iteration finished--calculate weights, */
/*                             ** abscissae for (-1,1) */
	gmu[k] = -x;
/* Computing 2nd power */
	d__1 = en * pm2;
	gwt[k] = two / (tmp * (d__1 * d__1));
	gmu[np1 - k] = -gmu[k];
	gwt[np1 - k] = gwt[k];
/* L30: */
    }
/*                                    ** set middle abscissa and weight */
/*                                    ** for rules of odd order */
    if (*m % 2 != 0) {
	gmu[lim + 1] = 0.f;
	prod = one;
	i__1 = *m;
	for (k = 3; k <= i__1; k += 2) {
	    prod = prod * k / (k - 1);
/* L40: */
	}
/* Computing 2nd power */
	d__1 = prod;
	gwt[lim + 1] = two / (d__1 * d__1);
    }
/*                                        ** convert from (-1,1) to (0,1) */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	gmu[k] = gmu[k] * .5f + .5f;
	gwt[k] *= .5f;
/* L50: */
    }
    return 0;
} /* qgausn_ */

doublereal ratio_(real *a, real *b)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double r_lg10(real *), r_sign(real *, real *);

    /* Local variables */
    static real absa, absb, huge__, powa, powb, tiny;
    extern doublereal r1mach_(integer *);
    static real powmin, powmax;

/*        calculate ratio  a/b  with over- and under-flow protection */
/*        (thanks to prof. jeff dozier for some suggestions here). */
/*        since this routine takes two logs, it is no speed demon, */
/*        but it is invaluable for comparing results from two runs */
/*        of a program under development. */

/*        note:  in fortran90, built-in functions tiny and huge */
/*               can replace the r1mach calls. */

/*   called by- disort */
/*   calls- r1mach */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
/*     .. save statement .. */
/*     .. */
/*     .. data statements .. */
/*     .. */
    if (pass1) {
	tiny = r1mach_(&c__1);
	huge__ = r1mach_(&c__2);
	powmax = r_lg10(&huge__);
	powmin = r_lg10(&tiny);
	pass1 = FALSE_;
    }
    if (*a == 0.f) {
	if (*b == 0.f) {
	    ret_val = 1.f;
	} else {
	    ret_val = 0.f;
	}
    } else if (*b == 0.f) {
	ret_val = r_sign(&huge__, a);
    } else {
	absa = dabs(*a);
	absb = dabs(*b);
	powa = r_lg10(&absa);
	powb = r_lg10(&absb);
	if (absa < tiny && absb < tiny) {
	    ret_val = 1.f;
	} else if (powa - powb >= powmax) {
	    ret_val = huge__;
	} else if (powa - powb <= powmin) {
	    ret_val = tiny;
	} else {
	    ret_val = absa / absb;
	}
/*                      ** dont use old trick of determining sign */
/*                      ** from a*b because a*b may (over/under)flow */
	if (*a > 0.f && *b < 0.f || *a < 0.f && *b > 0.f) {
	    ret_val = -ret_val;
	}
    }
    return ret_val;
} /* ratio_ */

/* Subroutine */ int slftst_(logical *corint, real *accur, real *albedo, real 
	*btemp, logical *deltam, real *dtauc, real *fbeam, real *fisot, 
	integer *ibcnd, logical *lamber, integer *nlyr, logical *plank, 
	integer *nphi, integer *numu, integer *nstr, integer *ntau, logical *
	onlyfl, real *phi, real *phi0, integer *nmom, real *pmom, logical *
	prnt, logical *prntu0, real *ssalb, real *temis, real *temper, real *
	ttemp, real *umu, logical *usrang, logical *usrtau, real *utau, real *
	umu0, real *wvnm, logical *compar, real *flup, real *rfldir, real *
	rfldn, real *uu)
{
    /* Initialized data */

    static real acc = 1e-4f;

    static integer i__, n;
    static logical tf, ok;
    static real phis, umus, phi0s, umu0s;
    static integer nphis, nmoms, ntaus;
    static real pmoms[5], utaus;
    static logical prnts[5];
    static integer nlyrs, numus, nstrs;
    static real wvnms, error1, error2, error3, error4;
    static logical prnu0s[2];
    static real albeds, fbeams;
    static logical lambes;
    static integer ibcnds;
    static logical deltas;
    static real accurs, dtaucs;
    extern logical tstbad_(char *, real *, ftnlen);
    static real ssalbs;
    static logical planks;
    static real btemps;
    static logical corins;
    static real tempes[2];
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static real temiss, fisots;
    static logical onlyfs, usrans;
    static real ttemps;
    static logical usrtas;

/*       if  compar = false, save user input values that would otherwise */
/*       be destroyed and replace them with input values for self-test. */
/*       if  compar = true, compare self-test case results with correct */
/*       answers and restore user input values if test is passed. */

/*       (see file 'disort.doc' for variable definitions.) */


/*     i n t e r n a l    v a r i a b l e s: */

/*         acc     relative accuracy required for passing self-test */

/*         errorn  relative errors in disort output variables */

/*         ok      logical variable for determining failure of self-test */

/*         all variables ending in 's' are temporary 's'torage for input */

/*   called by- disort */
/*   calls- tstbad, errmsg */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. local arrays .. */
/*     .. */
/*     .. external functions .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --prntu0;
    --prnt;

    /* Function Body */
    if (! (*compar)) {
/*                                     ** save user input values */
	nlyrs = *nlyr;
	dtaucs = *dtauc;
	ssalbs = *ssalb;
	for (n = 0; n <= 4; ++n) {
	    pmoms[n] = pmom[n];
/* L10: */
	}
	nstrs = *nstr;
	nmoms = *nmom;
	usrans = *usrang;
	numus = *numu;
	umus = *umu;
	usrtas = *usrtau;
	ntaus = *ntau;
	utaus = *utau;
	nphis = *nphi;
	phis = *phi;
	ibcnds = *ibcnd;
	fbeams = *fbeam;
	umu0s = *umu0;
	phi0s = *phi0;
	fisots = *fisot;
	lambes = *lamber;
	albeds = *albedo;
	deltas = *deltam;
	onlyfs = *onlyfl;
	corins = *corint;
	accurs = *accur;
	planks = *plank;
	wvnms = *wvnm;
	btemps = *btemp;
	ttemps = *ttemp;
	temiss = *temis;
	tempes[0] = temper[0];
	tempes[1] = temper[1];
	for (i__ = 1; i__ <= 5; ++i__) {
	    prnts[i__ - 1] = prnt[i__];
/* L20: */
	}
	for (i__ = 1; i__ <= 2; ++i__) {
	    prnu0s[i__ - 1] = prntu0[i__];
/* L30: */
	}
/*                                     ** set input values for self-test */
	*nstr = 4;
	*nlyr = 1;
	*dtauc = 1.f;
	*ssalb = .9f;
	*nmom = 4;
/*                          ** haze l moments */
	pmom[0] = 1.f;
	pmom[1] = .8042f;
	pmom[2] = .646094f;
	pmom[3] = .481851f;
	pmom[4] = .359056f;
	*usrang = TRUE_;
	*numu = 1;
	*umu = .5f;
	*usrtau = TRUE_;
	*ntau = 1;
	*utau = .5f;
	*nphi = 1;
	*phi = 90.f;
	*ibcnd = 0;
	*fbeam = 3.14159265f;
	*umu0 = .866f;
	*phi0 = 0.f;
	*fisot = 1.f;
	*lamber = TRUE_;
	*albedo = .7f;
	*deltam = TRUE_;
	*onlyfl = FALSE_;
	*corint = TRUE_;
	*accur = 1e-4f;
	*plank = TRUE_;
	*wvnm = 0.f;
	*btemp = 300.f;
	*ttemp = 100.f;
	*temis = .8f;
	temper[0] = 210.f;
	temper[1] = 200.f;
	for (i__ = 1; i__ <= 5; ++i__) {
	    prnt[i__] = FALSE_;
/* L40: */
	}
	for (i__ = 1; i__ <= 2; ++i__) {
	    prntu0[i__] = FALSE_;
/* L50: */
	}
    } else {
/*                                    ** compare test case results with */
/*                                    ** correct answers and abort if bad */
	ok = TRUE_;
	error1 = (*uu - 47.865571f) / 47.865571f;
	error2 = (*rfldir - 1.527286f) / 1.527286f;
	error3 = (*rfldn - 28.372225f) / 28.372225f;
	error4 = (*flup - 152.585284f) / 152.585284f;
	if (dabs(error1) > acc) {
	    ok = tstbad_("uu", &error1, (ftnlen)2);
	}
	if (dabs(error2) > acc) {
	    ok = tstbad_("rfldir", &error2, (ftnlen)6);
	}
	if (dabs(error3) > acc) {
	    ok = tstbad_("rfldn", &error3, (ftnlen)5);
	}
	if (dabs(error4) > acc) {
	    ok = tstbad_("flup", &error4, (ftnlen)4);
	}
	tf = TRUE_;
	if (! ok) {
	    errmsg_("disort--self-test failed", &tf, (ftnlen)24);
	}
/*                                      ** restore user input values */
	*nlyr = nlyrs;
	*dtauc = dtaucs;
	*ssalb = ssalbs;
	for (n = 0; n <= 4; ++n) {
	    pmom[n] = pmoms[n];
/* L60: */
	}
	*nstr = nstrs;
	*nmom = nmoms;
	*usrang = usrans;
	*numu = numus;
	*umu = umus;
	*usrtau = usrtas;
	*ntau = ntaus;
	*utau = utaus;
	*nphi = nphis;
	*phi = phis;
	*ibcnd = ibcnds;
	*fbeam = fbeams;
	*umu0 = umu0s;
	*phi0 = phi0s;
	*fisot = fisots;
	*lamber = lambes;
	*albedo = albeds;
	*deltam = deltas;
	*onlyfl = onlyfs;
	*corint = corins;
	*accur = accurs;
	*plank = planks;
	*wvnm = wvnms;
	*btemp = btemps;
	*ttemp = ttemps;
	*temis = temiss;
	temper[0] = tempes[0];
	temper[1] = tempes[1];
	for (i__ = 1; i__ <= 5; ++i__) {
	    prnt[i__] = prnts[i__ - 1];
/* L70: */
	}
	for (i__ = 1; i__ <= 2; ++i__) {
	    prntu0[i__] = prnu0s[i__ - 1];
/* L80: */
	}
    }
    return 0;
} /* slftst_ */

/* Subroutine */ int zeroal_(integer *nd1, real *expbea, real *flyr, real *
	oprim, real *phasa, real *phast, real *phasm, real *taucpr, real *xr0,
	 real *xr1, integer *nd2, real *cmu, real *cwt, real *psi0, real *
	psi1, real *wk, real *z0, real *z1, real *zj, integer *nd3, real *
	ylm0, integer *nd4, real *array, real *cc, real *evecc, integer *nd5, 
	real *gl, integer *nd6, real *ylmc, integer *nd7, real *ylmu, integer 
	*nd8, real *kk, real *ll, real *zz, real *zplk0, real *zplk1, integer 
	*nd9, real *gc, integer *nd10, integer *layru, real *utaupr, integer *
	nd11, real *gu, integer *nd12, real *z0u, real *z1u, real *zbeam, 
	integer *nd13, real *eval, integer *nd14, real *amb, real *apb, 
	integer *nd15, integer *ipvt, real *z__, integer *nd16, real *rfldir, 
	real *rfldn, real *flup, real *uavg, real *dfdt, integer *nd17, real *
	albmed, real *trnmed, integer *nd18, real *u0u, integer *nd19, real *
	uu)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n;

/*         zero arrays; ndn is dimension of all arrays following */
/*         it in the argument list */
/*   called by- disort */
/* -------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
    /* Parameter adjustments */
    --uu;
    --u0u;
    --trnmed;
    --albmed;
    --dfdt;
    --uavg;
    --flup;
    --rfldn;
    --rfldir;
    --z__;
    --ipvt;
    --apb;
    --amb;
    --eval;
    --zbeam;
    --z1u;
    --z0u;
    --gu;
    --utaupr;
    --layru;
    --gc;
    --zplk1;
    --zplk0;
    --zz;
    --ll;
    --kk;
    --ylmu;
    --ylmc;
    --gl;
    --evecc;
    --cc;
    --array;
    --ylm0;
    --zj;
    --z1;
    --z0;
    --wk;
    --psi1;
    --psi0;
    --cwt;
    --cmu;
    --xr1;
    --xr0;
    --taucpr;
    --phasm;
    --phast;
    --phasa;
    --oprim;
    --flyr;
    --expbea;

    /* Function Body */
    i__1 = *nd1;
    for (n = 1; n <= i__1; ++n) {
	expbea[n] = 0.f;
	flyr[n] = 0.f;
	oprim[n] = 0.f;
	phasa[n] = 0.f;
	phast[n] = 0.f;
	phasm[n] = 0.f;
	taucpr[n] = 0.f;
	xr0[n] = 0.f;
	xr1[n] = 0.f;
/* L10: */
    }
    i__1 = *nd2;
    for (n = 1; n <= i__1; ++n) {
	cmu[n] = 0.f;
	cwt[n] = 0.f;
	psi0[n] = 0.f;
	psi1[n] = 0.f;
	wk[n] = 0.f;
	z0[n] = 0.f;
	z1[n] = 0.f;
	zj[n] = 0.f;
/* L20: */
    }
    i__1 = *nd3;
    for (n = 1; n <= i__1; ++n) {
	ylm0[n] = 0.f;
/* L30: */
    }
    i__1 = *nd4;
    for (n = 1; n <= i__1; ++n) {
	array[n] = 0.f;
	cc[n] = 0.f;
	evecc[n] = 0.f;
/* L40: */
    }
    i__1 = *nd5;
    for (n = 1; n <= i__1; ++n) {
	gl[n] = 0.f;
/* L50: */
    }
    i__1 = *nd6;
    for (n = 1; n <= i__1; ++n) {
	ylmc[n] = 0.f;
/* L60: */
    }
    i__1 = *nd7;
    for (n = 1; n <= i__1; ++n) {
	ylmu[n] = 0.f;
/* L70: */
    }
    i__1 = *nd8;
    for (n = 1; n <= i__1; ++n) {
	kk[n] = 0.f;
	ll[n] = 0.f;
	zz[n] = 0.f;
	zplk0[n] = 0.f;
	zplk1[n] = 0.f;
/* L80: */
    }
    i__1 = *nd9;
    for (n = 1; n <= i__1; ++n) {
	gc[n] = 0.f;
/* L90: */
    }
    i__1 = *nd10;
    for (n = 1; n <= i__1; ++n) {
	layru[n] = 0;
	utaupr[n] = 0.f;
/* L100: */
    }
    i__1 = *nd11;
    for (n = 1; n <= i__1; ++n) {
	gu[n] = 0.f;
/* L110: */
    }
    i__1 = *nd12;
    for (n = 1; n <= i__1; ++n) {
	z0u[n] = 0.f;
	z1u[n] = 0.f;
	zbeam[n] = 0.f;
/* L120: */
    }
    i__1 = *nd13;
    for (n = 1; n <= i__1; ++n) {
	eval[n] = 0.f;
/* L130: */
    }
    i__1 = *nd14;
    for (n = 1; n <= i__1; ++n) {
	amb[n] = 0.f;
	apb[n] = 0.f;
/* L140: */
    }
    i__1 = *nd15;
    for (n = 1; n <= i__1; ++n) {
	ipvt[n] = 0;
	z__[n] = 0.f;
/* L150: */
    }
    i__1 = *nd16;
    for (n = 1; n <= i__1; ++n) {
	rfldir[n] = 0.f;
	rfldn[n] = 0.f;
	flup[n] = 0.f;
	uavg[n] = 0.f;
	dfdt[n] = 0.f;
/* L160: */
    }
    i__1 = *nd17;
    for (n = 1; n <= i__1; ++n) {
	albmed[n] = 0.f;
	trnmed[n] = 0.f;
/* L170: */
    }
    i__1 = *nd18;
    for (n = 1; n <= i__1; ++n) {
	u0u[n] = 0.f;
/* L180: */
    }
    i__1 = *nd19;
    for (n = 1; n <= i__1; ++n) {
	uu[n] = 0.f;
/* L190: */
    }
    return 0;
} /* zeroal_ */

/* Subroutine */ int zeroit_(real *a, integer *length)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l;

/*         zeros a real array a having length elements */

/*   called by- disort, albtrn, solve1, surfac, setmtx, solve0, fluxes */
/* -------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
    /* Parameter adjustments */
    --a;

    /* Function Body */
    i__1 = *length;
    for (l = 1; l <= i__1; ++l) {
	a[l] = 0.f;
/* L10: */
    }
    return 0;
} /* zeroit_ */

/* ****************************************************************** */
/* ********** end of disort service routines ************************ */
/* ****************************************************************** */
/* ****************************************************************** */
/* ********** ibcnd=1 special case routines ************************* */
/* ****************************************************************** */
/* Subroutine */ int albtrn_(real *albedo, real *amb, real *apb, real *array, 
	real *b, real *bdr, real *cband, real *cc, real *cmu, real *cwt, real 
	*dtaucp, real *eval, real *evecc, real *gl, real *gc, real *gu, 
	integer *ipvt, real *kk, real *ll, integer *nlyr, integer *nn, 
	integer *nstr, integer *numu, logical *prnt, real *taucpr, real *umu, 
	real *u0u, real *wk, real *ylmc, real *ylmu, real *z__, doublereal *
	aad, doublereal *evald, doublereal *eveccd, doublereal *wkd, integer *
	mi, integer *mi9m2, integer *maxumu, integer *mxcmu, integer *mxumu, 
	integer *nnlyri, real *sqt, real *albmed, real *trnmed)
{
    /* System generated locals */
    integer amb_dim1, amb_offset, apb_dim1, apb_offset, array_dim1, 
	    array_offset, bdr_dim1, bdr_offset, cband_dim1, cband_offset, 
	    cc_dim1, cc_offset, evecc_dim1, evecc_offset, gc_dim1, gc_dim2, 
	    gc_offset, gl_dim1, gl_offset, gu_dim1, gu_dim2, gu_offset, 
	    kk_dim1, kk_offset, ll_dim1, ll_offset, u0u_dim1, u0u_offset, 
	    ylmc_dim1, ylmc_offset, ylmu_dim1, ylmu_offset, aad_dim1, 
	    aad_offset, eveccd_dim1, eveccd_offset, i__1, i__2;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer l, lc;
    static logical tf;
    static integer iq, iu, ncd;
    static real sgn;
    static integer ncol, ncut;
    static real delm0;
    extern /* Subroutine */ int sgbco_(real *, integer *, integer *, integer *
	    , integer *, integer *, real *, real *);
    static real rcond;
    static integer mazim;
    static real fisot;
    extern /* Subroutine */ int solve1_(real *, real *, real *, integer *, 
	    integer *, real *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static logical lamber;
    static real sphalb;
    extern /* Subroutine */ int soleig_(real *, real *, real *, real *, real *
	    , real *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), altrin_(real *, real *,
	     real *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *), errmsg_(
	    char *, logical *, ftnlen), lepoly_(integer *, integer *, integer 
	    *, integer *, real *, real *, real *), praltr_(real *, integer *, 
	    real *, real *), spaltr_(real *, real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    ), terpev_(real *, real *, real *, real *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    ), zeroit_(real *, integer *);
    static real sphtrn;
    static logical lyrcut;
    extern /* Subroutine */ int setmtx_(real *, real *, real *, real *, real *
	    , real *, real *, real *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, real *, real *);

/*    disort special case to get only albedo and transmissivity */
/*    of entire medium as a function of incident beam angle */
/*    (many simplifications because boundary condition is just */
/*    isotropic illumination, there are no thermal sources, and */
/*    particular solutions do not need to be computed).  see */
/*    ref. s2 and references therein for details. */
/*    the basic idea is as follows.  the reciprocity principle leads to */
/*    the following relationships for a plane-parallel, vertically */
/*    inhomogeneous medium lacking thermal (or other internal) sources: */

/*       albedo(theta) = u_0(theta) for unit-intensity isotropic */
/*                       illumination at *top* boundary */

/*       trans(theta) =  u_0(theta) for unit-intensity isotropic */
/*                       illumination at *bottom* boundary */

/*    where */

/*       albedo(theta) = albedo for beam incidence at angle theta */
/*       trans(theta) = transmissivity for beam incidence at angle theta */
/*       u_0(theta) = upward azim-avg intensity at top boundary */
/*                    at angle theta */
/*   o u t p u t    v a r i a b l e s: */

/*       albmed(iu)   albedo of the medium as a function of incident */
/*                    beam angle cosine umu(iu) */

/*       trnmed(iu)   transmissivity of the medium as a function of */
/*                    incident beam angle cosine umu(iu) */
/*    i n t e r n a l   v a r i a b l e s: */
/*       ncd         number of diagonals below/above main diagonal */
/*       rcond       estimate of the reciprocal condition of matrix */
/*                   cband; for system  cband*x = b, relative */
/*                   perturbations in cband and b of size epsilon may */
/*                   cause relative perturbations in x of size */
/*                   epsilon/rcond.  if rcond is so small that */
/*                          1.0 + rcond .eq. 1.0 */
/*                   is true, then cband may be singular to working */
/*                   precision. */
/*       cband       left-hand side matrix of linear system eq. sc(5), */
/*                   scaled by eq. sc(12); in banded form required */
/*                   by linpack solution routines */
/*       ncol        number of columns in cband matrix */
/*       ipvt        integer vector of pivot indices */
/*       (most others documented in disort) */
/*   called by- disort */
/*   calls- lepoly, zeroit, sgbco, soleig, terpev, setmtx, solve1, */
/*          altrin, spaltr, praltr */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --dtaucp;
    --ipvt;
    --prnt;
    eveccd_dim1 = *mi;
    eveccd_offset = 1 + eveccd_dim1;
    eveccd -= eveccd_offset;
    --evald;
    aad_dim1 = *mi;
    aad_offset = 1 + aad_dim1;
    aad -= aad_offset;
    --eval;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    apb_dim1 = *mi;
    apb_offset = 1 + apb_dim1;
    apb -= apb_offset;
    amb_dim1 = *mi;
    amb_offset = 1 + amb_dim1;
    amb -= amb_offset;
    --trnmed;
    --albmed;
    --umu;
    --wkd;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1;
    ylmu -= ylmu_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1;
    ylmc -= ylmc_offset;
    --wk;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2);
    gc -= gc_offset;
    gl_dim1 = *mxcmu - 0 + 1;
    gl_offset = 0 + gl_dim1;
    gl -= gl_offset;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1;
    evecc -= evecc_offset;
    --cwt;
    --cmu;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1;
    cc -= cc_offset;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1;
    array -= array_offset;
    u0u_dim1 = *mxumu;
    u0u_offset = 1 + u0u_dim1;
    u0u -= u0u_offset;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2);
    gu -= gu_offset;
    --z__;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1;
    cband -= cband_offset;
    --b;
    --sqt;

    /* Function Body */
    mazim = 0;
    delm0 = 1.f;
/*                    ** set disort variables that are ignored in this */
/*                    ** special case but are needed below in argument */
/*                    ** lists of subroutines shared with general case */
    ncut = *nlyr;
    lyrcut = FALSE_;
    fisot = 1.f;
    lamber = TRUE_;
/*                          ** get legendre polynomials for computational */
/*                          ** and user polar angle cosines */
    i__1 = *nstr - 1;
    lepoly_(numu, &mazim, mxcmu, &i__1, &umu[1], &sqt[1], &ylmu[ylmu_offset]);
    i__1 = *nstr - 1;
    lepoly_(nn, &mazim, mxcmu, &i__1, &cmu[1], &sqt[1], &ylmc[ylmc_offset]);
/*                       ** evaluate legendre polynomials with negative */
/*                       ** arguments from those with positive arguments; */
/*                       ** dave/armstrong eq. (15), stwl(59) */
    sgn = -1.f;
    i__1 = *nstr - 1;
    for (l = mazim; l <= i__1; ++l) {
	sgn = -sgn;
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    ylmc[l + iq * ylmc_dim1] = sgn * ylmc[l + (iq - *nn) * ylmc_dim1];
/* L10: */
	}
/* L20: */
    }
/*                                  ** zero out bottom reflectivity */
/*                                  ** (albedo is used only in analytic */
/*                                  ** formulae involving albedo = 0 */
/*                                  ** solutions; eqs 16-17 of ref s2) */
    i__1 = *mi * (*mi + 1);
    zeroit_(&bdr[bdr_offset], &i__1);
/* ===================  begin loop on computational layers  ============= */
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
/*                                       ** solve eigenfunction problem */
/*                                       ** in eq. stwj(8b), stwl(23f) */
	soleig_(&amb[amb_offset], &apb[apb_offset], &array[array_offset], &
		cmu[1], &cwt[1], &gl[lc * gl_dim1], mi, &mazim, mxcmu, nn, 
		nstr, &ylmc[ylmc_offset], &cc[cc_offset], &evecc[evecc_offset]
		, &eval[1], &kk[lc * kk_dim1 + 1], &gc[(lc * gc_dim2 + 1) * 
		gc_dim1 + 1], &aad[aad_offset], &eveccd[eveccd_offset], &
		evald[1], &wkd[1]);
/*                          ** interpolate eigenvectors to user angles */
	terpev_(&cwt[1], &evecc[evecc_offset], &gl[lc * gl_dim1], &gu[(lc * 
		gu_dim2 + 1) * gu_dim1 + 1], &mazim, mxcmu, mxumu, nn, nstr, 
		numu, &wk[1], &ylmc[ylmc_offset], &ylmu[ylmu_offset]);
/* L30: */
    }
/* ===================  end loop on computational layers  =============== */
/*                      ** set coefficient matrix (cband) of equations */
/*                      ** combining boundary and layer interface */
/*                      ** conditions (in band-storage mode required by */
/*                      ** linpack routines) */
    setmtx_(&bdr[bdr_offset], &cband[cband_offset], &cmu[1], &cwt[1], &delm0, 
	    &dtaucp[1], &gc[gc_offset], &kk[kk_offset], &lamber, &lyrcut, mi, 
	    mi9m2, mxcmu, &ncol, &ncut, nnlyri, nn, nstr, taucpr, &wk[1]);
/*                      ** lu-decompose the coeff. matrix (linpack) */
    ncd = *nn * 3 - 1;
    sgbco_(&cband[cband_offset], mi9m2, &ncol, &ncd, &ncd, &ipvt[1], &rcond, &
	    z__[1]);
    tf = FALSE_;
    if (rcond + 1.f == 1.f) {
	errmsg_("albtrn--sgbco says matrix near singular", &tf, (ftnlen)39);
    }
/*                             ** first, illuminate from top; if only */
/*                             ** one layer, this will give us everything */
/*                             ** solve for constants of integration in */
/*                             ** homogeneous solution */
    solve1_(&b[1], &cband[cband_offset], &fisot, &c__1, &ipvt[1], &ll[
	    ll_offset], mi9m2, mxcmu, &ncol, nlyr, nn, nnlyri, nstr);
/*                             ** compute azimuthally-averaged intensity */
/*                             ** at user angles; gives albedo if multi- */
/*                             ** layer (eq. 9 of ref s2); gives both */
/*                             ** albedo and transmissivity if single */
/*                             ** layer (eqs. 3-4 of ref s2) */
    altrin_(&gu[gu_offset], &kk[kk_offset], &ll[ll_offset], mxcmu, mxumu, 
	    maxumu, nlyr, nn, nstr, numu, taucpr, &umu[1], &u0u[u0u_offset], &
	    wk[1]);
/*                               ** get beam-incidence albedos from */
/*                               ** reciprocity principle */
    i__1 = *numu / 2;
    for (iu = 1; iu <= i__1; ++iu) {
	albmed[iu] = u0u[iu + *numu / 2 + u0u_dim1];
/* L40: */
    }
    if (*nlyr == 1) {
	i__1 = *numu / 2;
	for (iu = 1; iu <= i__1; ++iu) {
/*                               ** get beam-incidence transmissivities */
/*                               ** from reciprocity principle (1 layer); */
/*                               ** flip them end over end to correspond */
/*                               ** to positive umu instead of negative */
	    trnmed[iu] = u0u[*numu / 2 + 1 - iu + (u0u_dim1 << 1)] + exp(
		    -taucpr[*nlyr] / umu[iu + *numu / 2]);
/* L50: */
	}
    } else {
/*                             ** second, illuminate from bottom */
/*                             ** (if multiple layers) */
	solve1_(&b[1], &cband[cband_offset], &fisot, &c__2, &ipvt[1], &ll[
		ll_offset], mi9m2, mxcmu, &ncol, nlyr, nn, nnlyri, nstr);
	altrin_(&gu[gu_offset], &kk[kk_offset], &ll[ll_offset], mxcmu, mxumu, 
		maxumu, nlyr, nn, nstr, numu, taucpr, &umu[1], &u0u[
		u0u_offset], &wk[1]);
/*                               ** get beam-incidence transmissivities */
/*                               ** from reciprocity principle */
	i__1 = *numu / 2;
	for (iu = 1; iu <= i__1; ++iu) {
	    trnmed[iu] = u0u[iu + *numu / 2 + u0u_dim1] + exp(-taucpr[*nlyr] /
		     umu[iu + *numu / 2]);
/* L60: */
	}
    }
    if (*albedo > 0.f) {
/*                             ** get spherical albedo and transmissivity */
	if (*nlyr == 1) {
	    spaltr_(&cmu[1], &cwt[1], &gc[gc_offset], &kk[kk_offset], &ll[
		    ll_offset], mxcmu, nlyr, nn, nstr, taucpr, &sphalb, &
		    sphtrn);
	} else {
	    spaltr_(&cmu[1], &cwt[1], &gc[gc_offset], &kk[kk_offset], &ll[
		    ll_offset], mxcmu, nlyr, nn, nstr, taucpr, &sphtrn, &
		    sphalb);
	}
/*                                ** ref. s2, eqs. 16-17 (these eqs. have */
/*                                ** a simple physical interpretation */
/*                                ** like that of adding-doubling eqs.) */
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    albmed[iu] += *albedo / (1.f - *albedo * sphalb) * sphtrn * 
		    trnmed[iu];
	    trnmed[iu] += *albedo / (1.f - *albedo * sphalb) * sphalb * 
		    trnmed[iu];
/* L70: */
	}
    }
/*                          ** return umu to all positive values, to */
/*                          ** agree with ordering in albmed, trnmed */
    *numu /= 2;
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	umu[iu] = umu[iu + *numu];
/* L80: */
    }
    if (prnt[4]) {
	praltr_(&umu[1], numu, &albmed[1], &trnmed[1]);
    }
    return 0;
} /* albtrn_ */

/* Subroutine */ int altrin_(real *gu, real *kk, real *ll, integer *mxcmu, 
	integer *mxumu, integer *maxumu, integer *nlyr, integer *nn, integer *
	nstr, integer *numu, real *taucpr, real *umu, real *u0u, real *wk)
{
    /* System generated locals */
    integer gu_dim1, gu_dim2, gu_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, u0u_dim1, u0u_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer lc, iq, iu, lu;
    static real mu, sgn, exp1, exp2, dtau, expn, denom;
    static integer iumin, iumax;
    static real palint, utaupr[2];

/*       computes azimuthally-averaged intensity at top and bottom */
/*       of medium (related to albedo and transmission of medium by */
/*       reciprocity principles; see ref s2).  user polar angles are */
/*       used as incident beam angles. (this is a very specialized */
/*       version of usrint) */

/*       ** note **  user input values of umu (assumed positive) are */
/*                   temporarily in upper locations of  umu  and */
/*                   corresponding negatives are in lower locations */
/*                   (this makes gu come out right).  i.e. the contents */
/*                   of the temporary umu array are: */

/*                     -umu(numu),..., -umu(1), umu(1),..., umu(numu) */


/*   i n p u t    v a r i a b l e s: */

/*       gu     :  eigenvectors interpolated to user polar angles */
/*                   (i.e., g in eq. sc(1), stwl(31ab)) */

/*       kk     :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b) */

/*       ll     :  constants of integration in eq. sc(1), obtained */
/*                   by solving scaled version of eq. sc(5); */
/*                   exponential term of eq. sc(12) not included */

/*       nn     :  order of double-gauss quadrature (nstr/2) */

/*       taucpr :  cumulative optical depth (delta-m-scaled) */

/*       (remainder are disort input variables) */


/*   o u t p u t    v a r i a b l e: */

/*       u0u  :    diffuse azimuthally-averaged intensity at top and */
/*                 bottom of medium (directly transmitted component, */
/*                 corresponding to bndint in usrint, is omitted). */


/*   i n t e r n a l    v a r i a b l e s: */

/*       dtau   :  optical depth of a computational layer */
/*       palint :  non-boundary-forced intensity component */
/*       utaupr :  optical depths of user output levels (delta-m scaled) */
/*       wk     :  scratch vector for saving 'exp' evaluations */
/*       all the exponential factors (i.e., exp1, expn,... etc.) */
/*       come from the substitution of constants of integration in */
/*       eq. sc(12) into eqs. s1(8-9).  all have negative arguments. */

/*   called by- albtrn */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. local arrays .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --wk;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    u0u_dim1 = *mxumu;
    u0u_offset = 1 + u0u_dim1;
    u0u -= u0u_offset;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2);
    gu -= gu_offset;
    --umu;

    /* Function Body */
    utaupr[0] = 0.f;
    utaupr[1] = taucpr[*nlyr];
    for (lu = 1; lu <= 2; ++lu) {
	if (lu == 1) {
	    iumin = *numu / 2 + 1;
	    iumax = *numu;
	    sgn = 1.f;
	} else {
	    iumin = 1;
	    iumax = *numu / 2;
	    sgn = -1.f;
	}
/*                                   ** loop over polar angles at which */
/*                                   ** albedos/transmissivities desired */
/*                                   ** ( upward angles at top boundary, */
/*                                   ** downward angles at bottom ) */
	i__1 = iumax;
	for (iu = iumin; iu <= i__1; ++iu) {
	    mu = umu[iu];
/*                                     ** integrate from top to bottom */
/*                                     ** computational layer */
	    palint = 0.f;
	    i__2 = *nlyr;
	    for (lc = 1; lc <= i__2; ++lc) {
		dtau = taucpr[lc] - taucpr[lc - 1];
		exp1 = exp((utaupr[lu - 1] - taucpr[lc - 1]) / mu);
		exp2 = exp((utaupr[lu - 1] - taucpr[lc]) / mu);
/*                                      ** kk is negative */
		i__3 = *nn;
		for (iq = 1; iq <= i__3; ++iq) {
		    wk[iq] = exp(kk[iq + lc * kk_dim1] * dtau);
		    denom = mu * kk[iq + lc * kk_dim1] + 1.f;
		    if (dabs(denom) < 1e-4f) {
/*                                                   ** l'hospital limit */
			expn = dtau / mu * exp2;
		    } else {
			expn = (exp1 * wk[iq] - exp2) * sgn / denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * ll[iq 
			    + lc * ll_dim1] * expn;
/* L10: */
		}
/*                                        ** kk is positive */
		i__3 = *nstr;
		for (iq = *nn + 1; iq <= i__3; ++iq) {
		    denom = mu * kk[iq + lc * kk_dim1] + 1.f;
		    if (dabs(denom) < 1e-4f) {
			expn = -dtau / mu * exp1;
		    } else {
			expn = (exp1 - exp2 * wk[*nstr + 1 - iq]) * sgn / 
				denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * ll[iq 
			    + lc * ll_dim1] * expn;
/* L20: */
		}
/* L30: */
	    }
	    u0u[iu + lu * u0u_dim1] = palint;
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* altrin_ */

/* Subroutine */ int praltr_(real *umu, integer *numu, real *albmed, real *
	trnmed)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double acos(doublereal);

    /* Local variables */
    static integer iu;

    /* Fortran I/O blocks */
    static cilist io___514 = { 0, 6, 0, "(///,a,//,a)", 0 };
    static cilist io___516 = { 0, 6, 0, "(0p,f13.4,f20.6,f12.5,1p,e17.4)", 0 }
	    ;


/*        print planar albedo and transmissivity of medium */
/*        as a function of incident beam angle */
/*   called by- albtrn */
/* -------------------------------------------------------------------- */
/*     .. parameters .. */
/*     .. */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    --trnmed;
    --albmed;
    --umu;

    /* Function Body */
    s_wsfe(&io___514);
    do_fio(&c__1, " *******  flux albedo and/or transmissivity of entire med"
	    "ium  ********", (ftnlen)70);
    do_fio(&c__1, " beam zen ang   cos(beam zen ang)      albedo   transmiss"
	    "ivity", (ftnlen)62);
    e_wsfe();
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	s_wsfe(&io___516);
	r__1 = acos(umu[iu]) * 57.295779578552292f;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&albmed[iu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&trnmed[iu], (ftnlen)sizeof(real));
	e_wsfe();
/* L10: */
    }
    return 0;
} /* praltr_ */

/* Subroutine */ int solve1_(real *b, real *cband, real *fisot, integer *ihom,
	 integer *ipvt, real *ll, integer *mi9m2, integer *mxcmu, integer *
	ncol, integer *ncut, integer *nn, integer *nnlyri, integer *nstr)
{
    /* System generated locals */
    integer cband_dim1, cband_offset, ll_dim1, ll_offset, i__1, i__2;

    /* Local variables */
    static integer i__, lc, iq, ncd, ipnt;
    extern /* Subroutine */ int sgbsl_(real *, integer *, integer *, integer *
	    , integer *, integer *, real *, integer *), zeroit_(real *, 
	    integer *);

/*        construct right-hand side vector b for isotropic incidence */
/*        (only) on either top or bottom boundary and solve system */
/*        of equations obtained from the boundary conditions and the */
/*        continuity-of-intensity-at-layer-interface equations */


/*     i n p u t      v a r i a b l e s: */

/*       cband    :  left-hand side matrix of banded linear system */
/*                   eq. sc(5), scaled by eq. sc(12); assumed already */
/*                   in lu-decomposed form, ready for linpack solver */

/*       ihom     :  direction of illumination flag (1, top; 2, bottom) */

/*       ncol     :  number of columns in cband */

/*       nn       :  order of double-gauss quadrature (nstr/2) */

/*       (remainder are disort input variables) */


/*    o u t p u t     v a r i a b l e s: */

/*       b        :  right-hand side vector of eq. sc(5) going into */
/*                   sgbsl; returns as solution vector of eq. */
/*                   sc(12), constants of integration without */
/*                   exponential term */

/*       ll      :   permanent storage for b, but re-ordered */


/*    i n t e r n a l    v a r i a b l e s: */

/*       ipvt     :  integer vector of pivot indices */
/*       ncd      :  number of diagonals below or above main diagonal */

/*   called by- albtrn */
/*   calls- zeroit, sgbsl */
/* +-------------------------------------------------------------------+ */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. external subroutines .. */
/*     .. */
    /* Parameter adjustments */
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    --ipvt;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1;
    cband -= cband_offset;
    --b;

    /* Function Body */
    zeroit_(&b[1], nnlyri);
    if (*ihom == 1) {
/*                             ** because there are no beam or emission */
/*                             ** sources, remainder of b array is zero */
	i__1 = *nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__] = *fisot;
	    b[*ncol - *nn + i__] = 0.f;
/* L10: */
	}
    } else if (*ihom == 2) {
	i__1 = *nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__] = 0.f;
	    b[*ncol - *nn + i__] = *fisot;
/* L20: */
	}
    }
    ncd = *nn * 3 - 1;
    sgbsl_(&cband[cband_offset], mi9m2, ncol, &ncd, &ncd, &ipvt[1], &b[1], &
	    c__0);
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	ipnt = lc * *nstr - *nn;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ll[*nn + 1 - iq + lc * ll_dim1] = b[ipnt + 1 - iq];
	    ll[iq + *nn + lc * ll_dim1] = b[iq + ipnt];
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* solve1_ */

/* Subroutine */ int spaltr_(real *cmu, real *cwt, real *gc, real *kk, real *
	ll, integer *mxcmu, integer *nlyr, integer *nn, integer *nstr, real *
	taucpr, real *sflup, real *sfldn)
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, i__1, i__2;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer iq, jq;
    static real zint;

/*       calculates spherical albedo and transmissivity for the entire */
/*       medium from the m=0 intensity components */
/*       (this is a very specialized version of fluxes) */


/*    i n p u t    v a r i a b l e s: */

/*       cmu,cwt    abscissae, weights for gauss quadrature */
/*                  over angle cosine */

/*       kk      :  eigenvalues of coeff. matrix in eq. ss(7) */

/*       gc      :  eigenvectors at polar quadrature angles, sc(1) */

/*       ll      :  constants of integration in eq. sc(1), obtained */
/*                  by solving scaled version of eq. sc(5); */
/*                  exponential term of eq. sc(12) not included */

/*       nn      :  order of double-gauss quadrature (nstr/2) */

/*       (remainder are disort input variables) */


/*    o u t p u t   v a r i a b l e s: */

/*       sflup   :  up-flux at top (equivalent to spherical albedo due to */
/*                  reciprocity).  for illumination from below it gives */
/*                  spherical transmissivity */

/*       sfldn   :  down-flux at bottom (for single layer, equivalent to */
/*                  spherical transmissivity due to reciprocity) */


/*    i n t e r n a l   v a r i a b l e s: */

/*       zint    :  intensity of m=0 case, in eq. sc(1) */

/*   called by- albtrn */
/* +-------------------------------------------------------------------- */
/*     .. scalar arguments .. */
/*     .. */
/*     .. array arguments .. */
/*     .. */
/*     .. local scalars .. */
/*     .. */
/*     .. intrinsic functions .. */
/*     .. */
    /* Parameter adjustments */
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2);
    gc -= gc_offset;
    --cwt;
    --cmu;

    /* Function Body */
    *sflup = 0.f;
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	zint = 0.f;
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + gc_dim2) * gc_dim1] * ll[jq + ll_dim1] * 
		    exp(kk[jq + kk_dim1] * taucpr[1]);
/* L10: */
	}
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + gc_dim2) * gc_dim1] * ll[jq + ll_dim1];
/* L20: */
	}
	*sflup += cwt[iq - *nn] * cmu[iq - *nn] * zint;
/* L30: */
    }
    *sfldn = 0.f;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zint = 0.f;
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * ll[jq + *nlyr 
		    * ll_dim1];
/* L40: */
	}
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * ll[jq + *nlyr 
		    * ll_dim1] * exp(-kk[jq + *nlyr * kk_dim1] * (taucpr[*
		    nlyr] - taucpr[*nlyr - 1]));
/* L50: */
	}
	*sfldn += cwt[*nn + 1 - iq] * cmu[*nn + 1 - iq] * zint;
/* L60: */
    }
    *sflup *= 2.f;
    *sfldn *= 2.f;
    return 0;
} /* spaltr_ */

