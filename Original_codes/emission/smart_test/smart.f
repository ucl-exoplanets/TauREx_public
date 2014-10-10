      subroutine smart(iugrp,iuthrm,iuabc,iuaer,iusur,iusol0,
     -     iusol2,iuout,iuheat,iustat,iutrn,iuflx,nmodes,
     -     ngases,ngastyp,igtrn,igpab,gasfile,iref,nref,
     -     nlev,nlyr,next,k_out,dp_dp,ncomp,icomp,volmix,
     -     radius,wgtatm,ts,grav,p,t,alt,z,rmix,dtauaer,
     -     wnmin,wnmax,wn_tol,isptype,islit,width,dwn,
     -     usrang,nstr,numu,umu,nphi,phi,nza,umu0,phi0,
     -     levout,nlout,ws,phiw,lamber,lplanck,lsolar,irad,
     -     iunits,ifrmout,io_end,io_err,wn_eof,accur,ratm,
     -     io_file,tauerr,pi0err,phferr,surferr,aid_lr,
     -     dirsoflx,dnsoflx,upsoflx,dnthflx,upthflx,
     -     soheat,thheat,pd_frac,nstate,istate,ntau_pd,
     -     igs,iu_pd,iutpd,iuspd,iupdrad)
c
cccccccccccccccccccccccccccc   s m a r t  cccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This program combines spectral mapping methods and the          cc
cc    discrete ordinate model (Stamnes et al. 1988) to generate       cc
cc    high-resolution synthetic monochromatic radiances in            cc
cc    vertically-inhomogeneous scattering absorbing, emitting,        cc
cc    planetary atmospheres.                                          cc
cc                                                                    cc
cc    note:  This version of the program uses a full brdf and         cc
cc           surface optical propertiy gradients for interpolation.   cc
cc           It also outputs wavelength-dependent, level-dependent    cc
cc           contribution functions.                                  cc
cc                                                                    cc
cc   9/28/97:   this version prints the direct solar flux as well as  cc
cc              the total downward and upward fluxes, and allows      cc
cc              wavelength-dependent radiances to be printed          cc
cc              at mxlout levels.                                     cc
cc                                                                    cc
cc              as part of this modification, the slit function       cc
cc              routine, slit, was modified such that it only allows  cc
cc              boxcar or triangular slits (other functions required  cc
cc              far too much RAM.                                     cc
cc                                                                    cc
cc    5/3/98:   the aerosol, gas, and rayliegh scattering optical     cc
cc              depth sections have been extracted from the main      cc
cc              program and are called as subroutines.                cc
cc                                                                    cc
cc    6/3/98:   smart.f was modified to include a 0th order           cc
cc              correction for the solar zenith angle to account for  cc
cc              the sphericity of the planet.                         cc
cc                                                                    cc
cc    5/28/01   smart modified to evaluate radiances for the          cc
cc              minimum and maximum optical depth within each         cc
cc              bin, and interpolate linearly in tau                  cc
cc                                                                    cc
cc    5/29/03:  smart modified to use DISORT 2.0                      cc
cc                                                                    cc
cc    12/26/03: smart modified to provide an option to elminate       cc
cc              headers from binary out files (ifrmout>2)             cc
cc    11/14/04: bdrf mode (iref) propagated through program           cc
cc 12/04 -4/05: radiance partial derivatives found for optical depth  cc
cc              single scattering albedos, atmospheric temperatures,  cc
cc              surface pressures and surface albedos,                cc
cc    05/06/05: radiance partial derivatives for scattering phase     cc
cc              functions added.                                      cc
cc    05/08/05: redundant spectral grid points are skipped (skip_wn)  cc
cc    02/11/07: layer transmittance and absorptance routines changed  cc
cc              from a full-range interpolation to a bin-by-bin       cc
cc              interpolation based on disort trnmed, albmed          cc
cc    02/11/07: calls to flush and timing routines swapped out to     cc
cc              make the routines more compatible with f95.           cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc     iusol0 - unit number as a scratch file for solar fluxes        cc
cc     iusol2 - unit number for solar flux scratch file               cc
cc     iuthrm - unit number for thermal radiances                     cc
cc      iuout - unit number of output radiance file                   cc
cc      iutrn - unit number of output transmission/pressure file      cc
cc      iutpd - unit numbers for output level-dependent               cc
cc              thermal fluxes and thier partial derivatives          cc
cc      iuspd - unit numbers for output level-dependent               cc
cc              solar fluxes and thier partial derivatives            cc
cc    iupdrad - unit numbers for output level-dependent radiances     cc
cc              and thier partial derivatives                         cc
cc     nstate - number of elements in the state vecror                cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc              1 - surface pressure                                  cc
cc              2 - surface/atmospheric temperature                   cc
cc              3 - gas absorption coeffient                          cc
cc              4 - cloud/aerosol optical depth                       cc
cc              5 - surface albedo                                    cc
cc    ntau_pd - number of gases and aerosols that are variable        cc
cc              components of thestate vector that affect the optical cc
cc              depth.                                                cc
cc    gasfile - name fo file with gas absorption coeffiecients vs. wn cc
cc     nmodes - number of discrete aerosol partical modes             cc
cc      ncomp - number of rayleigh-scattering constituentes           cc
cc      icomp - index of each rayleigh scattering constituent         cc
cc     volmix - volume mixing ratio of each rayleigh scatterer        cc
cc         ts - surface temperature (K)                               cc
cc         au - distance to the sun (in AU's)                         cc
cc     tauerr - optical depth relative binning error (0. to ~0.8)     cc
cc     pi0err - co-single scattering albedo absolute binning error    cc
cc     phferr - asymmetry factor absolute binning error               cc
cc    surferr - surface optical property binning error                cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, use a BRDF is used.              cc
cc       irad - index of output file type:                            cc
cc              1) fluxes, radiances, and heating rates,              cc
cc                 at computational azimuths and zenith angles,       cc
cc              2) fluxes, radiances, heating rates, and transmission cc
cc                functions at computational zenith angles,           cc
cc              3) fluxes, radiances, heating rates, and contribution cc
cc                 functions at computational zenith angles,          cc
cc              4) fluxes, radiances, heating rates, transmission     cc
cc                 functions and and contribution functions           cc
cc                 at computational zenith angles,                    cc
cc              5) fluxes, radiances, and heating rates,              cc
cc                 at computational azimuths and zenith angles,       cc
cc              6) fluxes, radiances and transmission functions       cc
cc                 at arbitrary zenith angles,                        cc
cc              7) fluxes, radiances, and contribution functions      cc
cc                 at arbitrary zenith angles,                        cc
cc              8) fluxes, radiances, transmission functions,and      cc
cc                 contribution functions at arbitrary zenith angles. cc
cc    ifmrout - index of output file format (1) ascii, (2) binary,    cc
cc                                          (3) binary, no header     cc
cc       iref - bidirectional reflectance options                     cc
cc              0 - Lambert                                           cc
cc              1 - Hapke's BDR model                                 cc
cc              2 - Breon's BDR model; combination of Li + Roujean    cc
cc              3 - Roujean's BDR model                               cc
cc              4 - Cox and Munk glint model                          cc
cc       nref - number of surface optical properties specified at     cc
cc              each wavelength.                                      cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       nphi - number of output azimuth angles                       cc
cc       nlev - number of output levels                               cc
cc       nlyr - number of computational model layers                  cc
cc     levout - output level index (1) top of atmosphere,             cc
cc              (2) surface, (3) arbitrary level                      cc
cc        nza - number of solar zenith angles                         cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) Watts/m**2/sr/Angstrom                             cc
cc        phi - emission azimuth angles (degrees)                     cc
cc       umu0 - solar zenith angle cosines                            cc
cc       phi0 - solar azimuth angle cosines                           cc
cc      wnmin - minimum wavenumber of desired spectral window         cc
cc      wnmax - maximum wavenumber of desired spectral window         cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc         wl - current wavenumber (cm**-1)                           cc
cc         wn - current wavenumber (cm**-1)                           cc
cc        rad - radiance at each output zenith and azimuth angle      cc
cc  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc
cc  trn_gas_0 - normal-incidence gas transmission                     cc
cc  trn_ray_0 - normal-incidence aerosol transmission                 cc
cc                                                                    cc
cccccccccccccccccccccccccccc   s m a r t  cccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nfmx
      parameter (nfmx = (mxlout+mxpd+1)*nsol+1)
c
      character*132 io_file(nex),gasfile(ngtmax,ngas)
      character*1  name(132)
c
c****   define variables used by discrete ordinate method
c
      logical usrang,lamber,lplanck,lsolar
      logical usrang0
c
      integer ncomp,nlev,nlyr,k,l,ng,ngt,nf
      integer next,n,nze,nc0,nz,ngrtot
      integer ipi0
      integer iref,nref
      integer nspect,mom,nmom
      integer ismtime,irtime,icpu,ibtime
      integer ismterr,itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt,ieff
      integer numu,nstr,nphi,nzup,nzdn
      integer ne,nscat,ngrey,l0,nlyr0,numu0,ng0
c
      integer icomp(6),nmodes
      integer ngases,ngastyp(ngas)
      integer iuabc(ngtmax,ngas),igtrn(ngtmax,ngas),igpab(ngtmax,ngas)
      integer io_err(nex),io_end(nex)
      integer iustat,iuaer,iuheat,ii,len,nlb,ntb
      integer iugrp,iuthrm,iusol0,iusol2,iuout,
     -        iutrn,iuflx,iusur,iu_flux
      integer ifrmout,iunits,levout,nlout,irad,
     -        nza,nz0,isptype,islit
      integer k_out(mxlout)
      integer iflext(nex),nsiext(2,nex)
c
c****   input aerosol wavelength-dependent variables set in aer_tau
c
      integer nmomaer(2,nmode),nmom_mx(nmode)
c
c****   variables used for lbl gas absorption coefficients
c
      integer npg(ngtmax,ngas),ntg(ngtmax,ngas),
     -          igtyp(ngtmax,ngas),nbgas(ngtmax,ngas)
c
c****   index of optical depth unity for each bin
c
      integer levtau1(ngrp,2),nmomgrp(ngrp),iwngrp,ngroup
      integer ntau_pd,ntau_pd0,nlay,nt_pd,modepd(mxpd),igs(mxpd),ntaug
      integer n_rad(5,ngrp)
      integer ij,ij0,iwnflg,iwnmode,m
c
c****    units for flux and radiance jacobians
c 
      integer nstate,istate(nex),npd,ipd(mxpd),
     -        iu_pd,iutpd(mxpd,nsol),iuspd(mxpd,nsol),
     -        iupdrad(mxumu,mxphi,mxpd,mxlout,nsol)
c
c****  spectral grid parameters
c
      double precision wnmin,wnmax,width,dwn,
     -                 wnsmin,wnsmax,wn_tol,delnu
c
      double precision wn,wn0,wn1,distmin,distmax,dist(nex),dist0
      double precision wni(3),dwnmx,wn_eof(2,nex),wntot
c
      real ratm,tau_tot,tauabs,tausca
      real h0,theight,d,x,t0
      real taumn,pi0
      real umu0nz,phi0nz
      real tcpu,bcputime,smcputime,rcputime,eff
      real rtime,btime
      real time0,time1,time2,time3
c      
      real cp_atm
c
      real pi,ang_umu,umu_o(mxumu), gwt_o(mxumu)
c
c***    set up direction cosines for flux integration
c
      real umu_f(mxumu),gwt_f(mxumu)
c      
      real phi(mxphi),umu(mxumu)
      real phi0(nsol),umu0(nsol),umu0_1(nsol)
      real dtauaer(nmode,kp)
      real volmix(6),ts,wgtatm
      real dtausc(kp)
      real surf_0(2,4),dsurfdv(4)
      real tauerr,pi0err,phferr,surferr
      real dp_dp(mxlout),accur,ws,phiw
      real units,sza0,rtd
c
c****    pressure of tau = 1 and column optical depths
c
      real p_ray_0,p_gas_0,p_aer_0
      real tau_ray,tau_gas,tau_aer
      real taugas(kp),tauray(kp),tauaer(kp)
      real p_ray(mxumu),p_gas(mxumu),p_aer(mxumu)
c
c****   input aerosol wavelength-dependent variables set in aer_tau
c
      real aerqext(2,nmode),aerqsca(2,nmode),aerg0(2,nmode),
     -     aerpmom(2,0:mxmom,nmode),dpmomdv(0:mxmom,nmode),
     -     dqextdv(nmode),dqscadv(nmode),dg0dv(nmode)
c
c****   atmospheric structure variables
c
      real p(kp),t(kp),alt(kp),grav(kp),z(kp),radius
c
c****    gas mixing ratios 
c
      real rmix(kp,ngas)
c
c****   atmospheric optical properties
c
      real dtauex(kp,2),co_pi0(kp,2),g(kp),phmom(0:mxmom,kp),
     -     surf_opt(4),alb(nsol)
c
c****    optical depths of each variable component of state vector
c
      real tau_ext(kp,3,mxpd),tau_sca(kp,mxpd),g_sca(kp,mxpd)
c
c****   variables used in skip_wn
c
      real dtauexi(kp,2,3),dtausci(kp,3),co_pi0i(kp,2,3),gi(kp,3),
     -     phmomi(0:mxmom,kp,3),surf_opti(4,3),albi(nsol,3),
     -     small_tau,small_pi0,small_g,small_alb
c
      real tau_exti(kp,3,mxpd,3),tau_scai(kp,mxpd,3),g_scai(kp,mxpd,3)
c
      real p_ray_0i(3),p_rayi(mxumu,3),tau_rayi(3),
     -     p_gas_0i(3),p_gasi(mxumu,3),tau_gasi(3),
     -     p_aer_0i(3),p_aeri(mxumu,3),tau_aeri(3)
c
      double precision wnskip
c
      real sur_b(4,5,ngrp),alb_b(nsol,5,ngrp),
     -     dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp),pmom_b(0:mxmom,kp,5,ngrp),
     -     dx_b_i(kp,5,ngrp),dalb_b_i(nsol,ngrp)
c
c****   binned solar source function variables and jacobians
c
      double precision upsflxsrc_b(kp,ngrp),dnsflxsrc_b(kp,ngrp),
     -                 sradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddnsflxdx(kp,4,ngrp),dupsflxdx(kp,4,ngrp),
     -                 dsraddx(mxumu,mxphi,kp,4,ngrp)
c
c****   binned thermal source function variables and jacobians
c
      double precision uptflxsrc_b(kp,ngrp),dntflxsrc_b(kp,ngrp),
     -                 tradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddntflxdx(kp,4,ngrp),duptflxdx(kp,4,ngrp),
     -                 dtraddx(mxumu,mxphi,kp,4,ngrp)
c
c****     monochormatic radiance interpolation variables
c 
      double precision wnext(2,nex)
c
c****   variables for single/multiple scattering test
c
      real ss_min,dtausca
c
c****   partial derivative fractional change
c
      real pd_frac(nex),pd_pert(mxpd)
c
c****   number of points in slit function algorithm
c
      real points(nfmx)
c
c****    variables for spectral grouping algorithm.
c
      real taugrp(kp,ngrp,3),pi0grp(kp,ngrp,3),
     -     ggrp(kp,ngrp,3),pmomgrp(0:mxmom,kp,ngrp,3),
     -     surfgrp(ngrp,4,3),albgrp(ngrp,nsol,3)
c
c****   spectral grid monitoring varibles
c
      double precision wn_step,dwn_step
c
      double precision wngrp(ngrp,3),dnugrp(ngrp)
c
c****   double precision variables for heating rates
c
      double precision aid_lr(mxulv),dirsoflx(mxulv,nsol),
     -       dnsoflx(mxulv,nsol),upsoflx(mxulv,nsol),dnthflx(mxulv),
     -       upthflx(mxulv),soheat(mxulv,nsol),thheat(mxulv)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp),
     -                 dtrnraddx(mxumu,kp,3,ngrp),
     -                 refrad_b(mxumu,kp,5,ngrp),
     -                 drefraddx(mxumu,kp,3,ngrp),
     -                 absrad_b(mxumu,kp,5,ngrp),
     -                 dabsraddx(mxumu,kp,3,ngrp),
     -                 brdf_b(mxumu,mxphi,5,ngrp),
     -                 dbrdfdx(mxumu,mxphi,ngrp)
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),dtrnflxdx(kp,3,ngrp),
     -                 refflx_b(kp,5,ngrp),drefflxdx(kp,3,ngrp),
     -                 absflx_b(kp,5,ngrp),dabsflxdx(kp,3,ngrp)
c
c****    initialize timing variables
c
      tcpu = 0.0
c
c****   define a small value that is added to optical depths.  this
c       quantity is needed to prevent an instability that develops
c       in the thermal radiance calculation if there are layers 
c       with  near-zero optical depth.
c
      small_tau = 1.e-5
c Danie Liang
      small_tau = 1.e-15
c
c****    set small values of single scattering parameter and 
c        asymmetry parameter for use in the spectral skipping
c        algorithm.
c
      small_g = 1.e-5
      small_pi0 = 1.e-5
      small_alb = 1.e-5
c Danie Liang
      small_g = 0.e-15
      small_pi0 = 0.e-15
      small_alb = 0.e-15
c
c****   Define a minimum single scattering albedo for scattering
c       calculations
c
      ss_min = 0.25*pi0err
c
c****    define a maximum spectral interval between output points.
c
      if(isptype .eq. 1) then
        dwnmx = 1.0d-3*(wnmax - wnmin)
      else
        dwnmx = dwn
      endif
c
c****   set the output interval for the spectral grid monitor
c
      dwn_step = 1.0d-2*(wnmax - wnmin)
      wn_step = wnmin + dwn_step
c
c****   modify solar zenith angles to approximate the sphericity of the
c       planet.  To do this, assume that the planet is spherical with
c       a 3-scale-height thick atmosphere. Define the ray tangent height
c                     theight = sin(umu0).  
c       The distance between the top of the atmosphere and the 
c       axis of rotation, x, is then given by the pythagorean theorm:
c           d = sqrt((radius + H)**2 - rt**2)
c       and the solar path length between the top of the atmosphere and
c       the surface is given by:
c           x = d - radius*cos(umu0)
c       The effective solar zenith angle cosine is therefore 
c       given by H/x.
c
      if(lsolar) then
c
c****    find the effective scale height:
c 
        write(*,'(/,1a,/,1a,/1a)') 'Effective Solar Zenith Angle: ',
     -              '    sza0       sza_eff',
     -              '  (degrees)   (degrees)'
c
      h0 = 0.003*ratm*t(nlev)/grav(nlev)
      rtd = 180./acos(-1.0)
c
c****    find the effective solar zenith angles:
c
        do 1141 nz=1,nza
            sza0 = rtd*acos(umu0(nz))
            theight = radius*sqrt(1.0 - umu0(nz)**2)
            d = sqrt((radius + h0)**2 - theight**2)
            x = d - radius*umu0(nz)
            umu0_1(nz) = umu0(nz)
            umu0(nz) = h0/x
            if(umu0(nz) .gt. 1.0) umu0(nz) = 1.0
            write(*,'(2(1pe12.4))') sza0,rtd*acos(umu0(nz))
1141    continue
c
      endif
c
c****    compute the adiabatic lapse rate, g/cp, at each level.
c        this quantity is needed for heating rates.
c
      nc0 = ncomp
      do 1201 k=1,nlev
c
c****    define the adiabatic lapse rate (k/m)
c
          t0 = t(k)
          aid_lr(k) = grav(k)/cp_atm(nc0,icomp,volmix,t0)
c
1201  continue
c
c****   initialize point counter in slit function routine
c
      if(isptype .ne. 1) then
        do 1241 nf=1,nfmx
            points(nf) = 0.0d0
1241    continue
      endif
c
c****  find gaussian points and weights for flux integration
c
      call qgausn(nstr/2,umu_o,gwt_o)
c
c****  restack the values into order used by disort
c
      do 1301 nze=1,nstr/2
           umu_f(nze+nstr/2) = umu_o(nze)
           gwt_f(nze+nstr/2) = gwt_o(nze)
1301  continue
      do 1321 nze=1,nstr/2
          umu_f(nze) = -umu_o(nstr/2-nze+1)
          gwt_f(nze) = -gwt_o(nstr/2-nze+1)
1321  continue
      write(*,'(/,1a,/,1a)') 
     - 'Computational Gaussian Points and Weights',
     - '  Angle (deg)    umu       gwt'
      pi = acos(-1.0)
      do 1341 nze=1,nstr
          ang_umu = 180.*acos(umu_f(nze))/pi
          write(*,'(3(1pe12.4))') ang_umu,umu_f(nze),gwt_f(nze)
1341  continue
c
c****   define the last downward and first upward looking angles
c
      if(usrang) then
      write(*,'(/,1a,/,1a)') 
     - 'User-Defined Emission Angles:',
     - '  Angle (deg)    umu '
        nzdn = 0
        do 1361 nze=1,numu
            if(umu(nze) .lt. 0.0) nzdn = nze
            ang_umu = 180.*acos(umu(nze))/pi
            write(*,'(3(1pe12.4))') ang_umu,umu(nze)
1361    continue
        nzup = nzdn + 1
      else
        nzdn = nstr/2
        nzup = nzdn+1
c
c****     initialize umu values for use in tau=1 calculation
c
        do 1381 nze=1,numu
             umu(nze) = umu_f(nze)
1381    continue
      endif
c
c****   initialze timing routines:
c
      call cpu_time(time0)
      tcpu = 0.0
c
c****   initialize wavenumbers
c
      wn = wnmin
      wn0 = wnmin
      wn1 = wnmin
      wntot = 1.0d0
      wnsmin = wnmin
      wnsmax = wnmin
      wnskip = 0.0d0
c
c****    set a maximum step size
c
      distmax = 0.01d0*(wnmax - wnmin)
      if(distmax .lt. 0.01d0*wnmin) distmax = 0.01d0*wnmin
c
c****    initialize flux summation arrays for heating rate calculation
c
      do 1521 k=1,nlev
          upthflx(k) = 0.0
          dnthflx(k) = 0.0
          do 1501 nz=1,nza
              upsoflx(k,nz) = 0.0
              dnsoflx(k,nz) = 0.0
              dirsoflx(k,nz) = 0.0
1501      continue
1521  continue
c
c****    initialize spectral mapping variables
c
      ngrtot = 0
      ngroup = 0
c
c****    define a small number as a minimum fractional tau error
c
      taumn = 1.e-5*(tauerr**2 + small_tau)
c
c****    determine whether pressure is a variable part of the state 
c        vector.  If it is, gas optical depths are computed for the
c        usual nlyr layers, and for a a "perturnbed" lower layer, 
c        which has 1% higher pressure.
c
      if(abs(istate(1)) .eq. 1) then
        nlay = nlyr + 1
      else
        nlay = nlyr
      endif 
c
c*****   determine whether temperature is a variable part of
c        the state vector.  if it is, extinction optical depth must be
c        computed for the background and perterbed temperature profile
c
      if(istate(2) .eq. 0) then
        nt_pd = 1
      else
        nt_pd = 2
      endif
c
c****     set jacobian indexes and fractional changes for each variable
c         component of the state vector 
c
      npd = 0
      do 1601 n=1,nstate
          if(istate(n) .ne. 0) then
            npd = npd + 1
            ipd(npd) = istate(n)
            pd_pert(npd) = pd_frac(n)
          endif
1601  continue          
c
c****   initialize absorption coefficient read flag, iflext
c
      do 1701 n=1,nex
          iflext(n) = 1
          nsiext(1,n) = 0
          nsiext(2,n) = 0
          wnext(2,n) = -999.0d0
          io_err(n) = 0
          io_end(n) = 0
1701  continue
c
c****    open the scratch unit for group numbers at each wavenumber
c
      close(iugrp)
      open(iugrp,form='unformatted',status='scratch')
c
c****    open the scratch unit for thermal fluxes at each wavenumber
c
      close(iuthrm)
      open(iuthrm,form='unformatted',status='scratch')
c
c****   set io file properties for rayleigh scattering
c
      ne = 1  
      io_file(ne) = 'Rayleigh Scattering'
      wn_eof(1,ne) = wnmin
      wn_eof(2,ne) = wnmax
      io_err(ne) = 0
      io_end(ne) = 0
c
c****      read header of each gas absorption line file:
c
      if(ngases .gt. 0) write(*,'(/,1a)') ' Gas Optical Property Files:'
      write(*,'(/,1a)') '  ngt   ng iuabc gasfile'
      do 2421 ng=1,ngases
          do 2401 ngt = 1,ngastyp(ng)
              ne = ne + 1
c
              call charsp(gasfile(ngt,ng),name,len,132,nlb,ntb) 
c
              if(igtrn(ngt,ng) .eq. 1) then
c    
                read(iuabc(ngt,ng)) npg(ngt,ng),ntg(ngt,ng),
     -                              igtyp(ngt,ng),nbgas(ngt,ng)
c
                write(*,'(5i5,2x,132a)') 
     -                    ngt,ng,iuabc(ngt,ng),npg(ngt,ng),ntg(ngt,ng),
     -                   (name(ii),ii=1,len)
              else
c
                write(*,'(3i5,2x,132a)') 
     -                    ngt,ng,iuabc(ngt,ng),(name(ii),ii=1,len)
c
              endif
c
2401      continue
2421  continue 
c
c****   initialze spectral optical property buffers used to determine
c       if a given spectral point contributes to the spectrum when
c       compared to its surrounding points
c
      ij0 = 1
      iwnflg = 1
      wni(ij0) = wn
c        
      do 2881 ij=1,3
          wni(ij) = 0.0d0
          do 2821 k=1,nlay
              dtauexi(k,1,ij) = small_tau
              dtauexi(k,2,ij) = small_tau
              dtausci(k,ij) = 0.0
              co_pi0i(k,1,ij) = 0.0
              co_pi0i(k,2,ij) = 0.0
              gi(k,ij) = 0.0
              do 2801 mom=0,mxmom
                  phmomi(mom,k,ij) = 0.0
2801          continue
              do 2811 n=1,mxpd
                  tau_exti(k,1,n,ij) = 0.0
                  tau_exti(k,2,n,ij) = 0.0
                  tau_exti(k,3,n,ij) = 0.0
                  tau_scai(k,n,ij) = 0.0
                  g_scai(k,n,ij) = 0.0
2811          continue
2821      continue
          do 2841 l=1,4
              surf_opti(l,ij) = 0.0
2841      continue
          do 2851 nz=1,nza
              albi(nz,ij) = 0.0
2851      continue
c
          if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -          .or. irad .eq. 8) then
c
            p_ray_0i(ij) = 0.0
            p_gas_0i(ij) = 0.0
            p_aer_0i(ij) = 0.0
            tau_rayi(ij) = 0.0
            tau_gasi(ij) = 0.0
            tau_aeri(ij) = 0.0
            do 2861 nze=1,numu
                p_rayi(nze,ij) = 0.0
                p_gasi(nze,ij) = 0.0
                p_aeri(nze,ij) = 0.0
2861        continue
          endif
2881  continue
c
c       write(*,'(/,1a,/,1a)') ' Progress Report: ',
c     -                       '     v(cm**-1)    # intervals'
c
c****        s p e c t r a l    m a p p i n g    l o o p           *****
c
3001  nspect = 0
c
c****         Initialize counter for number of grey and 
c             number of multiple-scattering calculations
c
      nscat = 0
      ngrey = 0
c
c****   initialize spectral binning parameters
c
      call init_bin(nlyr,nref,nza,ngroup,iwngrp,ismterr,
     -                     itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt,
     -                     levtau1,nmomgrp,dnugrp,wngrp,surfgrp,
     -                     albgrp,taugrp,pi0grp,ggrp,pmomgrp,
     -                     dtau_b,copi0_b,g_b,pmom_b,sur_b,alb_b,
     -                     dx_b_i,dalb_b_i)
c
c****           s t a r t    s p e c t r a l    l o o p            *****
c
3302  nspect = nspect + 1
c
c****      set extinction counter
c
          ne = 0
c
c****       initialize extinction and scattering optical depths  
c           and moments of the scattering phase function at wn
c           (note: values go to nlay = nlyr+1 to accommodate pressure 
c            partial derivatives).
c
          do 3441 k=1,nlay
              dtauex(k,1) = small_tau
              dtauex(k,2) = small_tau
              dtausc(k) = 0.0
              do 3401 n=1,mxpd
                  tau_ext(k,1,n) = 0.0
                  tau_ext(k,2,n) = 0.0
                  tau_ext(k,3,n) = 0.0
                  tau_sca(k,n) = 0.0
                  g_sca(k,n) = 0.0
3401          continue
              g(k) = 0.0
              do 3421 mom=0,mxmom
                  phmom(mom,k) = 0.0
3421          continue
3441      continue
          do 3461 l=1,4
              surf_opt(l) = 0.0
3461      continue
          do 3481 nz=1,nza
              alb(nz) = 0.0
3481      continue
c
c****        initialize the counter for variable optical depth 
c            components of the state vector
c
          ntau_pd0 = 0
c
c****       initialize the number of phase function moments to 2,
c           to account for Raleigh scattering.  More modes can be
c           added for aerosols in aer_tau.
c 
          nmom = nstr
c
c****            r a y l e i g h    s c a t t e r i n g            *****
c
c
          call ray_tau(ne,nlyr,ncomp,icomp,volmix,p,grav,
     -                 wgtatm,nzup,nstr,iflext,nsiext,istate,
     -                 pd_frac,wnext,wn,umu,dtauex,dtausc,phmom,
     -                 tauray,p_ray,p_ray_0)

c
          tau_ray = tauray(nlev)
c
c****            g a s    a b s o r p t i o n                      *****
c
c
          call gas_tau(iuabc,igtrn,igpab,ngases,ngastyp,
     -                  npg,ntg,nbgas,ne,nlev,nzup,nt_pd,ntaug,
     -                  nstr,iflext,nsiext,nlay,ntau_pd0,istate,
     -                  wn,wnmin,wnmax,ratm,umu,p,t,alt,grav,
     -                  z,rmix,pd_frac,wnext,dtauex,tau_ext,    
     -                  taugas,p_gas,p_gas_0,wn_eof,io_end,io_err)
c
c
          tau_gas = taugas(nlev)
c
c****       a e r o s o l    o p t i c a l    p r o p e r t i e s
c
c
          call aer_tau(ne,nlyr,iuaer,nmodes,nzup,nstr,nmom,
     -                   wn,wnmin,wnmax,umu,p,dtauex,g,phmom,
     -                   tau_ext,tau_sca,g_sca,wnext,
     -                   aerqext,aerqsca,aerg0,aerpmom,
     -                   dqextdv,dqscadv,dpmomdv,dg0dv,
     -                   dtausc,dtauaer,p_aer,p_aer_0,tauaer,
     -                   iflext,nsiext,nmomaer,nmom_mx,ntau_pd0,
     -                   istate,ngases,modepd,wn_eof,io_end,io_err)
c
c
          tau_aer = tauaer(nlev)
c
c****       normalize single scattering albedos and phase functions
c
c
          call norm_pi0_ph(nmom,nlay,nt_pd,small_tau,
     -                     dtauex,dtausc,co_pi0,g,phmom)
c
c
c*****      s u r f a c e    o p t i c a l    p r o p e r t i e s
c
c
          call albedo(ne,iusur,iflext,nsiext,iref,nref,nza,
     -                lamber,umu0,wn,wnmin,wnmax,wnext,
     -                surf_0,dsurfdv,ws,phiw,surf_opt,alb,
     -                wn_eof,io_end,io_err)
c
c
c*****     find the next wavenumber where monochromatic properties
c          are specified.  Check all input data sets and select the
c          wavenumber that is closest to the current wavenumber.
c
c          note: the minimum distance must exceed the round-off
c                tolerance of the computer, wn_tol
c
          distmin = 1.0d30
          dist0 = wn_tol*wn
c
          next = ne
          do 4001 n=1,next
              dist(n) = wnext(nsiext(2,n),n) - wn
              if(dist(n) .lt. distmin) then
                if(dist(n) .ge. dist0) then
                  distmin = dist(n)
                else
                  distmin = dist0
                endif
              endif          
4001      continue
c
c****       set the next spectral point 
c
          if(distmin .gt. distmax) distmin = distmax
          wn1 = wn + distmin
c
c****       for iwnmode = 1, determine if the optical properties are 
c           changing rapidly enough to justify an RT calculation for 
c           this point
c
          iwnmode = 1
c
c
          call skip_wn(iwnmode,ij0,nlay,iwnflg,ntau_pd,numu,irad,nza,
     -                   dwnmx,small_tau,small_pi0,small_g,small_alb,
     -                   tauerr,pi0err,phferr,surferr,
     -                   wn,dtauex,dtausc,co_pi0,g,phmom,surf_opt,alb,
     -                   wni,dtauexi,dtausci,co_pi0i,gi,phmomi,
     -                   surf_opti,albi,tau_ext,tau_sca,g_sca,
     -                   tau_exti,tau_scai,g_scai,
     -                   p_ray_0,p_ray,tau_ray,p_ray_0i,p_rayi,tau_rayi,
     -                   p_gas_0,p_gas,tau_gas,p_gas_0i,p_gasi,tau_gasi,
     -                   p_aer_0,p_aer,tau_aer,p_aer_0i,p_aeri,tau_aeri)
c
c*****     check to see all three spectral intervals ar loaded
c
          if(ij0 .lt. 3) then
            ij0 = ij0 + 1
            go to 3302
          endif   
c
c****       check to see if this spectral point is needed (iwnflg = 1)
c
c****        the following flag should shut down skip_wn
c
c           iwnflg =1      
c
          if(iwnflg .eq. 1) then
c
c*****       All 3 points are good. put x(1) into working arrays
c            and process that point, and then pack x(2) into x(1)
c            and x(3) into x(2)
c
            iwnmode = 2
c
c            write(*,*) 'above 2rd call skip_wn'
c            call flush(0)
c
            call skip_wn(iwnmode,ij0,nlay,iwnflg,ntau_pd,numu,irad,nza,
     -                   dwnmx,small_tau,small_pi0,small_g,small_alb,
     -                   tauerr,pi0err,phferr,surferr,
     -                   wn,dtauex,dtausc,co_pi0,g,phmom,surf_opt,alb,
     -                   wni,dtauexi,dtausci,co_pi0i,gi,phmomi,
     -                   surf_opti,albi,tau_ext,tau_sca,g_sca,
     -                   tau_exti,tau_scai,g_scai,
     -                   p_ray_0,p_ray,tau_ray,p_ray_0i,p_rayi,tau_rayi,
     -                   p_gas_0,p_gas,tau_gas,p_gas_0i,p_gasi,tau_gasi,
     -                   p_aer_0,p_aer,tau_aer,p_aer_0i,p_aeri,tau_aeri)
c
c
c****         set the next wavenumber and the effective spectral 
c             interval width, delnu.  The spectral interval is assumed
c             to extend between the current wavenumber, wn, and half-way
c             to the previous value, wn0, and the next value, wni(2).
c
            delnu = 0.5*real(wni(1) - wn0)
            wn0 = wn
c
c****         s p e c t r a l   b i n n i n g  s e c t i o n
c
c****       define the column-integrated optical depth and 
c           determine if scattering is important in this interval
c
            ipi0 = 0
            tauabs = 0.0
            tausca = 0.0
            tau_tot = 0.0
            do 4201 k=1,nlay
                tau_tot = tau_tot + dtauex(k,1)
                pi0 = (1.0 - co_pi0(k,1))
                tauabs = tauabs + co_pi0(k,1)*dtauex(k,1)
                dtausca = pi0*dtauex(k,1)
                tausca = tausca + dtausca
                if(dtausca .gt. small_tau .and. pi0 .gt. ss_min) 
     -             ipi0 = 1
4201        continue
c
c****           determine whether a full multiple scattering 
c               calcluation must be performed (map_spect) or whether
c               a simple extinction calculation will work.
c
c****           NOTE: setting ipi0 = 1 turns off the grey calculation
c
            if(ipi0 .eq. 1 .or. .not. lamber) then
c
c****         c r e a t e    s p e c t r a l    m a p
c
              call map_spect(nlyr,nlay,nstr,nmom,
     -                     levout,levtau1,nref,nza,
     -                     wn,delnu,dtauex,co_pi0,g,
     -                     phmom,surf_opt,alb,taumn,
     -                     tauerr,pi0err,phferr,surferr,
     -                     dnugrp,wngrp,surfgrp,albgrp,
     -                     taugrp,pi0grp,ggrp,pmomgrp,
     -                     nmomgrp,iwngrp,ngroup,ismterr,
     -                     itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt)
c
            else
c
c****            scattering is not important.  Set the bin number
c                to zero and perform a single-scattering calculation
c
              iwngrp = 0
              ngrey = ngrey + 1
              ismterr = 0
c
            endif
c
c****          write out the group number at each each wavenumber
c
            if(ismterr .eq. 0) then
              wnsmax = wn
c
c****            write group parameters for this wn
c
c
              write(iugrp) wn,iwngrp,
     -           (surf_opt(ii),ii=1,4),(alb(nz),nz=1,nza),
     -           ((dtauex(k,l),k=1,nlay),l=1,nt_pd),
     -           ((co_pi0(k,l),k=1,nlay),l=1,nt_pd),
     -           (g(k),k=1,nlyr),
     -           (((tau_ext(k,m,n),k=1,nlev),m=1,3),n=1,ntau_pd),
     -           ((tau_sca(k,n),k=1,nlev),n=1,ntau_pd),
     -           ((g_sca(k,n),k=1,nlev),n=1,ntau_pd)
c
              if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -           .or. irad .eq. 8) write(iugrp)
     -           p_ray_0,(p_ray(nze),nze=nzup,numu),tau_ray,
     -           p_gas_0,(p_gas(nze),nze=nzup,numu),tau_gas,
     -           p_aer_0,(p_aer(nze),nze=nzup,numu),tau_aer
c
c
            endif
c
c****         update wavenumber counter
c
            wntot = wntot + 1.0
c
c****        determine if this wavenumber is beyond end of desired
c            spectral range.  if so, jump out of wavenumber loop.
c
            if(ismterr .ne. 0 .or. wn .ge. wnmax) go to 4601
c
          else
c
c****       skip point x(2).  Load x(3) into x(2)
c
            wnskip = wnskip + 1.0d0
            iwnmode = 3
c
            call skip_wn(iwnmode,ij0,nlay,iwnflg,ntau_pd,numu,irad,nza,
     -                   dwnmx,small_tau,small_pi0,small_g,small_alb,
     -                   tauerr,pi0err,phferr,surferr,
     -                   wn,dtauex,dtausc,co_pi0,g,phmom,surf_opt,alb,
     -                   wni,dtauexi,dtausci,co_pi0i,gi,phmomi,
     -                   surf_opti,albi,tau_ext,tau_sca,g_sca,
     -                   tau_exti,tau_scai,g_scai,
     -                   p_ray_0,p_ray,tau_ray,p_ray_0i,p_rayi,tau_rayi,
     -                   p_gas_0,p_gas,tau_gas,p_gas_0i,p_gasi,tau_gasi,
     -                   p_aer_0,p_aer,tau_aer,p_aer_0i,p_aeri,tau_aeri)
c        
          endif
c
c*****        specify the next wavenumber grid point and reset
c             wavenumber grid flags.
c
          wn = wn1
c
c****         set extinction flags for each constituent
c
          do 4401 n=1,next
c
c****             find the distance between the last location where
c                 the extintion from this agent was specified and the
c                 current wavelength.  
c
              dist(n) = wnext(nsiext(2,n),n) - wn
              if(dist(n) .le. 0.0d0) then
c
c****              wn is beyond latest input point for this constituent.
c                  set flag to read next point.
c
                iflext(n) = 1
              else
                iflext(n) = 0
              endif
4401      continue
c
          if(wn .gt. wn_step) then
c
c****         update spectral grid monitor
c
            wn_step = wn_step + dwn_step
c
          endif
c
c****        bin monochromatic optical properties for next wavenumber.
c
          go to 3302
c
c****          Record time required for binning
c
4601      call cpu_time(time1)
          smcputime = time1 - time0
          ismtime = nint(smcputime)
          tcpu = tcpu + smcputime
c
c****        Update group counters
c
          ngrtot = ngrtot + ngroup
c
          do 4881 n=1,ngroup
c
c****      load optical properties into binned radiance arrays
c
              ng0 = n          
c
              call load_optics(lamber,ng0,nlyr,nza,nstate,istate,
     -                         n_rad,nref,iref,nmomgrp,
     -                         umu0,tauerr,pi0err,phferr,surferr,
     -                         taumn,surfgrp,albgrp,
     -                         taugrp,pi0grp,ggrp,pmomgrp,
     -                         alb_b,sur_b,phiw,ws,
     -                         dtau_b,copi0_b,g_b,pmom_b,
     -                         dx_b_i,dalb_b_i)
c
              usrang0 = usrang
              numu0 = numu
c
              do 4801 l=1,5
c
c****               find the layer transmittances, reflectances and 
c                   absorptances for each profile
c
                  nlyr0 = nlyr
                  l0 = l
                  nmom = nmomgrp(n)
                  numu0 = numu
c
                  call layer_trn(usrang0,l0,ng0,nlyr0,
     -                     nstr,numu0,nmom,nphi,n_rad,
     -                     dtau_b,copi0_b,pmom_b,dx_b_i,
     -                     umu,phi,umu_f,gwt_f,
     -                     trnflx_b,dtrnflxdx,refflx_b,drefflxdx,
     -                     absflx_b,dabsflxdx,trnrad_b,dtrnraddx,
     -                     refrad_b,drefraddx,absrad_b,dabsraddx)
c
4801          continue
c
4881      continue
c
c****        s o l a r     z e n i t h    a n g l e    l o o p
c
          rtime = 0.0
          btime = 0.0
          do 5301 nz=1,nza
              nz0 = nz
              iu_flux = iuflx + nz - 1
              umu0nz = umu0(nz)
              phi0nz = phi0(nz)
c
c****             r a d i a n c e s   f o r    e a c h   b i n 
c
              do 5221 n=1,ngroup
                  ng0 = n
c
                  call sm_eq_trn(usrang,lamber,lplanck,lsolar,ng0,nz0,
     -                     nlyr,nstr,numu,nzdn,nzup,nphi,nmomgrp,
     -                     nstate,istate,n_rad,iref,nref,nscat,
     -                     tauerr,pi0err,phferr,surferr,ws,phiw,taumn,
     -                     umu,umu_f,gwt_f,phi,umu0nz,phi0nz,ts,t,accur,
     -                     wngrp,surfgrp,taugrp,pi0grp,ggrp,pmomgrp,
     -                     alb_b,sur_b,dtau_b,copi0_b,g_b,pmom_b,
     -                     dx_b_i,dalb_b_i,trnflx_b,refflx_b,absflx_b,
     -                     dtrnflxdx,drefflxdx,dabsflxdx,
     -                     trnrad_b,refrad_b,absrad_b,brdf_b,
     -                     dtrnraddx,drefraddx,dabsraddx,dbrdfdx,
     -                     dnsflxsrc_b,upsflxsrc_b,sradsrc_b,
     -                     ddnsflxdx,dupsflxdx,dsraddx,
     -                     dntflxsrc_b,uptflxsrc_b,tradsrc_b,
     -                     ddntflxdx,duptflxdx,dtraddx)
c
5221          continue
c
c****            Record time needed for radiance calculations
c
              call cpu_time(time2)
              rcputime = time2 - time1
              rtime = rtime + rcputime
              tcpu = tcpu + rcputime
c
c****             Rewind unit with group index data for backmapping
c
              rewind(iugrp)
c
c****             s m t    b a c k    m a p p i n g
c
              call map_back(lsolar,lplanck,lamber,usrang,iugrp,iuthrm,
     -          iusol0,iusol2,iutrn,iu_flux,iuout,npd,ipd,iunits,irad,
     -          iu_pd,iutpd,iuspd,iupdrad,ifrmout,nza,nz0,nstr,numu,
     -          nzup,nzdn,nphi,nlyr,nlay,levout,nlout,k_out,igs,iref,
     -          modepd,nt_pd,nstate,istate,ntau_pd,isptype,islit,
     -          width,dwn,points,wnmin,wnmax,wnsmin,wnsmax,wn_tol,
     -          units,umu0nz,phi0nz,umu_f,gwt_f,umu,phi,pd_pert,
     -          p,t,ts,dp_dp,rmix,dtauaer,alb_b,dtau_b,copi0_b,g_b,
     -          trnflx_b,dtrnflxdx,refflx_b,drefflxdx,
     -          absflx_b,dabsflxdx,refrad_b,drefraddx,
     -          absrad_b,dabsraddx,brdf_b,dbrdfdx,
     -          dnsflxsrc_b,ddnsflxdx,upsflxsrc_b,dupsflxdx,
     -          dntflxsrc_b,ddntflxdx,uptflxsrc_b,duptflxdx,
     -          sradsrc_b,dsraddx,tradsrc_b,dtraddx,
     -          dirsoflx,dnsoflx,upsoflx,dnthflx,upthflx)
c
c****             Record time needed for back-mapping
c
              call cpu_time(time3)
              bcputime = time3 - time2
              btime = btime + bcputime
              tcpu = tcpu + bcputime
c
5301      continue
c
c****         print spectral grouping statistics
c
          if(wnsmin .eq. wnmin) then
            write(iustat,'(/,/1a,/,/,4a,/,3a/)')
     -        '   Spectral Binning Statistics: ',
     -        '   wnmin      wnmax   # bins  total  # mono   bin',
     -        '    bin rejection criteria',
     -        '                       calc. type',
     -        '        timing(seconds)     total',
     -        '   cm**-1     cm**-1  in int  bins segments  ratio',
     -        '   tau=1    tau <>1     pi0        g   albedo',
     -        '    #scat   #grey     bin     rad bk_map       cpu'
c
          endif
c
          irtime = nint(rtime)
          ibtime = nint(btime)
          icpu = nint(tcpu)
c
          ieff = 1
          if(ngroup .gt. 0) ieff = nspect/ngroup
c
          write(iustat,'(2(1pe11.4),i5,2i8,i5,2i10,3i9,5i8,i10)') 
     -           wnsmin,wnsmax,ngroup,ngrtot,nspect,ieff,
     -           itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt,nscat,
     -           ngrey,ismtime,irtime,ibtime,icpu
c
c****        find the solar and thermal heating rates
c
          call heating(iuheat,nlyr,nza,umu0_1,wnmin,wnsmax,p,t,alt,
     -                    aid_lr,dirsoflx,dnsoflx,upsoflx,
     -                    dnthflx,upthflx,soheat,thheat)
c
          if(wn .lt.wnmax) then
c
c****         rewind units for spectral bins and thermal fluxes
c             and empty all output buffers the spectral output units.
c
            rewind(iugrp)
            rewind(iuthrm)
c
c****           finish remainder of spectrum
c
             wnsmin = wn1
             wn = wn1
c 
             go to 3001
c    
          else
c
            close(iugrp,status='delete')
            close(iuthrm,status='delete')
c
          endif
c
c****          b i n n i n g   e f f i c i e n c y
c
          write(iustat,'(/,1a,1pe15.5)') 
     -      ' Total number of spectral intervals skipped =   ',wnskip
          write(iustat,'(1a,1pe15.5)') 
     -       ' Total number of monochromatic intervals used = ',wntot
          write(iustat,'(1a,i10)') 
     -      ' Total number of spectral mapping bins =        ',ngrtot
          eff = 1.0
          if(ngrtot .gt. 0) eff = real(wntot)/float(ngrtot)
          write(iustat,'(1a,1pe12.4)') ' Net efficiency =',eff 
          write(iustat,'(/,1a,1pe12.4,1a)') 
     -       ' Total CPU time =',tcpu,' seconds'
c
c********        E n d    o f    s p e c t r a l    l o o p    *********
c
      return
      end
     
