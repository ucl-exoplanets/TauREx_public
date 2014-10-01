      subroutine sm_eq_trn(usrang,lamber,lplanck,lsolar,ng0,nz0,
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
ccccccccccccccccccccccccc  s m _ e q _ t r n   ccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this routine uses a discrete ordinate method to find radiances  cc
cc    and fluxes for spectrally mapped optical properties.            cc
cc                                                                    cc
cc 3/18/07: This version of sm_eq_trn uses flux layer adding  for     cc
cc          both fluxes and radiances.  The same source function are  cc
cc          used for both solar and thermal regions.                  cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        ng0 - index  of current spectral bin                        cc
cc       nlyr - number of computational model layers                  cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       numu - number of output zenith angles used in D/O code       cc
cc       nphi - number of output azimuth angles                       cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc        umu - emission zenith angle cosines                         cc
cc        phi - emission azimuth angles (degrees)                     cc
cc        nz0 - index of current solar zenith angle (1 to nza)        cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, use a BRDF is used.              cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc      icbnd - disort variable (0- radiances, 1- albedos only)       cc
cc      fbeam - intensity of collimated flux at top of atmosphere     cc
cc      fisot - thermal flux at top of atmosphere                     cc
cc        wng - wavenumber of bin                                     cc
cc         ts - surface temperature (K)                               cc
cc     nstate - number of variable elements in the state vecror       cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc              1 - surface pressure                                  cc
cc              2 - surface/atmospheric temperature                   cc
cc              3 - gas absorption coeffient                          cc
cc              4 - cloud/aerosol optical depth                       cc
cc              5 - surface albedo                                    cc
cc       iref - bidirectional reflectance options                     cc
cc            0 - Lambert                                             cc
cc            1 - Hapke's BDR model                                   cc
cc            2 - Breon's BDR model; combination of Li + Roujean      cc
cc            3 - Roujean's BDR model                                 cc
cc            4 - Cox and Munk glint model                            cc
cc         ws - wind speed (m/s) for Cox/Munk model                   cc
cc       phiw - wind azimuth (deg) for Cox/Munk model                 cc
cc      n_rad - radiance calculations flag for each bin               cc
cc             1 : radiances for nominal state structure              cc
cc             2 : radiances for perturbed surface pressures          cc
cc             3 : radiances for perturbed optical depths             cc
cc             4 : radiances for perturbed single scattering albedos  cc
cc             5 : radiances for perturbed surface reflectance        cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      temis - emissivity of uppermost model layer (space)           cc
cc      ttemp - temperature of uppermost model layer (space)          cc
cc      sur_b - surface albedo                                        cc
cc     dtau_b - differential optical depth in each layer              cc
cc    copi0_b - single scattering co-albedo in each layer             cc 
cc      btemp - surface temperature                                   cc
cc     rfldir - downward direct flux at each model level              cc
cc      rfldn - downward diffuse flux at each model level             cc
cc         uu - upward radiance at each level, zenith angle, azimuth  cc
cc     albmed - albedo of system                                      cc
cc     trnmed - transmissivity of system                              cc
cc      nscat - counter for number of scattering calculations         cc
cc                                                                    cc
ccccccccccccccccccccccccc  s m _ e q _ t r n   ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
c****   define variables used by discrete ordinate method
c
      logical usrang,lamber,lplanck,lsolar,source
c
      integer ng0,nz0,nlyr,numu,nstr,nphi
      integer nzup,nzdn,nscat,nmom
      integer iref,nref,nstate,istate(nex)
      integer n_rad(5,ngrp)
      integer nlev,nz,k,l,mom,ibn,nr,l0,naz,nze
c      integer nlyr0
c      integer nze,naz
c
c****   spectral binning parameters
c
      integer nmomgrp(ngrp)
c
c****   spectral binning parameters
c
      double precision wngrp(ngrp,3)
c
      real taugrp(kp,ngrp,3),pi0grp(kp,ngrp,3),
     -     ggrp(kp,ngrp,3),pmomgrp(0:mxmom,kp,ngrp,3),
     -     surfgrp(ngrp,4,3)
      real tauerr,pi0err,phferr,surferr,taumn
c
c****   emission angle variables
c
      real phi(mxphi),umu(mxumu),umu0nz,phi0nz
      real umu_f(mxumu),gwt_f(mxumu)
c
c****   DISORT input variables
c
      real umu0,phi0,alb0,fbeam
      real surf_pr(4),ws,phiw,dtauc(kp),ssalb(kp),pmom(0:mxmom,kp),
     -     ts,t(kp),accur
c
c****   DISORT output
c
      real rfldir( mxulv ), rfldn( mxulv ), uu( mxumu, mxulv, mxphi )
c
c*****   state vector variables.
c
      real sur_b(4,5,ngrp),alb_b(nsol,5,ngrp),
     -     dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp),pmom_b(0:mxmom,kp,5,ngrp),
     -     dx_b_i(kp,5,ngrp),dalb_b_i(nsol,ngrp)
c
      real wng0
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),dtrnflxdx(kp,3,ngrp),
     -                 refflx_b(kp,5,ngrp),drefflxdx(kp,3,ngrp),
     -                 absflx_b(kp,5,ngrp),dabsflxdx(kp,3,ngrp)
c
      double precision flx_rb(kp,5),flx_rt(kp,5),
     -                 flx_dd(kp,5),flx_du(kp,5),
     -                 flx_urt(kp,5),flx_drb(kp,5)
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
c****   variables for radiance layer adding method
c
      double precision rad_rb(mxumu,mxphi,kp,5),rad_rt(mxumu,mxphi,kp,5)
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
c****   internal binned source function variables
c
      double precision upsflxsrc(kp,5),dnsflxsrc(kp,5),
     -                 sradsrc(mxumu,mxphi,kp,5)
c
c****   internal binned source function variables
c
      double precision uptflxsrc(kp,5),dntflxsrc(kp,5),
     -                 tradsrc(mxumu,mxphi,kp,5)
c
c****   source_fcn output variables
c
      double precision upflxsrc(kp,5),dnflxsrc(kp,5),
     -                 radsrc(mxumu,mxphi,kp,5)
c
c****  quantities for radiance calculation: downward and upward fluxes
c      at layer interfaces
c
      double precision flx_ft(kp,5),flx_fb(kp,5)
      double precision flx_dfb(kp,5),flx_uft(kp,5)
      double precision flx_dft(kp,5),flx_ufb(kp,5)
c
c****   specify number of levels, and solar zenith angle
c
      nlev = nlyr+1
      nz = nz0
c
      if(lsolar) then
        umu0 = umu0nz
        phi0 = phi0nz
      else
        umu0 = 0.0
        phi0 = 0.0
      endif
c
c****    set group counters
c
      ibn = ng0
      nmom = nmomgrp(ibn)
c
      wng0 = real(wngrp(ibn,1))
c
c****   find the layer reflectance quantities for fluxes and radiances
c       for the simplified adding method
c
      do 1401 l=1,5
c
          l0 = l
c
          if(l .eq. 1 .or. n_rad(l,ibn) .ne. 0) then
c
c****       define the surface optical properties
c
            alb0 = alb_b(nz,l,ibn)
            trnflx_b(nlev,l,ibn) = 0.0
            refflx_b(nlev,l,ibn) = alb0
c
            if(lamber) then
c
              do 1241 naz=1,nphi              
                  do 1221 nze=1,numu
                      brdf_b(nze,naz,l,ibn) = alb0
                      if(l .eq. 5) dbrdfdx(nze,naz,ibn) = (alb0 - 
     -                               brdf_b(nze,naz,1,ibn))*
     -                               dalb_b_i(nz,ibn)
1221              continue                  
1241          continue
c              
            else
c
c****          use disort to find brdf of surface
c
              call surf_brdf(usrang,lamber,
     -                     l0,ng0,nz0,nstr,numu,nphi,
     -                     nmom,n_rad,nref,umu,phi,dalb_b_i,
     -                     alb_b,sur_b,ws,phiw,umu0,phi0,
     -                     brdf_b,dbrdfdx)
c
            endif
c
c
c****         find the layer transmittances, reflectances and 
c             absorptances for combined layers
c
c
            call flx_ref(l0,nlyr,ng0,alb0,trnflx_b,refflx_b,absflx_b,
     -                    dx_b_i,flx_dd,flx_du,flx_urt,flx_drb,
     -                    flx_rb,flx_rt)
c
c****          find the radiance sources
c
            call rad_ref(l0,nlyr,ng0,numu,nzdn,nzup,nphi,
     -                    trnrad_b,refrad_b,brdf_b,flx_urt,flx_drb,
     -                    rad_rb,rad_rt)
c
          endif
c
1401  continue     
c
c******  o p t i c a l    d e p t h   p r o f i l e   l o o p   *******
c
c****   perform radiance calculations for the nomical case and 
c       for each variable component of the optical depth structure
c       and set the radiance calculation flag, n_rad, where
c       n_rad = 1 : radiances for nominal state structure
c             = 2 : radiances for perturbed optical depths
c             = 3 : radiances for perturbed single scattering albedos
c             = 4 : radiances for perturbed phase functions
c             = 5 : radiances for perturbed surface albedos
c
      do 2601 l=1,5
          l0 = l
c
          if(l .eq. 1 .or. n_rad(l,ibn) .ne. 0) then
c
c****        define the surface optical properties
c
            alb0 = alb_b(nz,l,ibn)
c
c****        set the surface optical properties for disort
c
            do 2021 nr=1,nref
                surf_pr(nr) = sur_b(nr,l,ibn)
2021        continue
c
            if(.not. lamber) then 
c
              if(iref .eq. 4) then
c
c****            set wind speed and direction for Cox/Munk model
c
                surf_pr(3) = ws
                surf_pr(4) = phiw
              endif
c
            endif
c
c****            set number of l-p phase function moments
c
            nmom = nmomgrp(ibn)
c
            do 2141 k=1,nlyr
c
c****            load the optical depth and single scattering albedos 
c                for use in the discrete ordinate algorithm
c
                dtauc(k) = dtau_b(k,l,ibn)
                ssalb(k) = 1.0 - copi0_b(k,l,ibn)
                if(ssalb(k) .lt. 0.0) ssalb(k) = 0.0
c
c****            load the phase function moments at each level
c
                pmom(0,k) = 1.0
                do 2121 mom=1,nmomgrp(ibn)
                    pmom(mom,k) =  pmom_b(mom,k,l,ibn)
2121            continue
2141        continue
c
          endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          if(lsolar) then
c
c****         find solar fluxes, radiances and source functions
c             in each layer
c 
            if(l .eq. 1 .or. n_rad(l,ibn) .ne. 0) then
c
c****             turn off thermal sources and set a unit direct flux
c                 at the top of the atmosphere
c
                source = .false.
                fbeam = 1.0
c
c****             find the solar flux and radince sources in each layer
c                 note: source_fcn output variables dnsflxsrc, 
c                 upsflxsrc, and sradsrc are specific to this call
c
                call  source_fcn(ng0,l0,nlyr,nzdn,nzup,nstr,nmom,
     -                      numu,nphi,iref,usrang,lamber,source,
     -                      umu,phi,umu0,phi0,accur,alb0,surf_pr,
     -                      dtauc,ssalb,pmom,ts,t,fbeam,wng0,
     -                      trnflx_b,refflx_b,trnrad_b,
     -                      flx_rt,flx_rb,flx_du,flx_dd,
     -                      flx_dft,flx_uft,flx_dfb,flx_ufb,
     -                      rad_rt,rad_rb,dnflxsrc,upflxsrc,
     -                      radsrc,rfldir,rfldn,
     -                      flx_ft,flx_fb,uu,nscat)
c
            endif
c
            call solar_src(ng0,l0,nz0,nlyr,nzdn,nzup,numu,nphi,n_rad,
     -                     alb0,umu0,rfldir,rfldn,dx_b_i,dalb_b_i,
     -                     trnrad_b,rad_rt,rad_rb,flx_ft,flx_fb,
     -                     trnflx_b,flx_rt,flx_rb,flx_du,flx_dd,
     -                     flx_dft,flx_uft,flx_dfb,flx_ufb,
     -                     dnflxsrc,upflxsrc,radsrc,
     -                     dnsflxsrc,upsflxsrc,sradsrc,
     -                     dnsflxsrc_b,upsflxsrc_b,sradsrc_b,
     -                     ddnsflxdx,dupsflxdx,dsraddx)
c
          endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          if(lplanck) then
c
c****         find thermal fluxes, radiances and source functions
c             in each layer
c 
            if(l .eq. 1 .or. n_rad(l,ibn) .ne. 0) then
c
c****             turn on thermal sources and set the direct flux
c                 at the top of the atmosphere to zero
c
                source = .true.
                fbeam = 0.0
c
c                 note: source_fcn output variables dntflxsrc, 
c                 uptflxsrc, and tradsrc are specific to this call
c
                call  source_fcn(ng0,l0,nlyr,nzdn,nzup,nstr,nmom,
     -                      numu,nphi,iref,usrang,lamber,source,
     -                      umu,phi,umu0,phi0,accur,alb0,surf_pr,
     -                      dtauc,ssalb,pmom,ts,t,fbeam,wng0,
     -                      trnflx_b,refflx_b,trnrad_b,
     -                      flx_rt,flx_rb,flx_du,flx_dd,
     -                      flx_dft,flx_uft,flx_dfb,flx_ufb,
     -                      rad_rt,rad_rb,dnflxsrc,upflxsrc,
     -                      radsrc,rfldir,rfldn,
     -                      flx_ft,flx_fb,uu,nscat)
c
            endif
c                
c****             find thermal flux and radince sources in each layer
c
            call thermal_src(ng0,l0,nz0,nlyr,nzdn,nzup,numu,nphi,n_rad,
     -                       wng0,alb0,brdf_b,ts,t,dx_b_i,dalb_b_i,
     -                       trnrad_b,rad_rt,rad_rb,flx_uft,flx_dfb,
     -                       dnflxsrc,upflxsrc,radsrc,
     -                       dntflxsrc,uptflxsrc,tradsrc,
     -                       copi0_b,absflx_b,absrad_b,
     -                       dntflxsrc_b,uptflxsrc_b,tradsrc_b,
     -                       ddntflxdx,duptflxdx,dtraddx)
c
          endif
c
2601  continue
c
      return
      end
