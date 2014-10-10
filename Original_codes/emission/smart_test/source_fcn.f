      subroutine source_fcn(ng0,l0,nlyr,nzdn,nzup,nstr,nmom,
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
ccccccccccccccccccccccccc  s o u r c e _ f c n   ccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the radiance and flux source functions cc
cc    for a given bin for smart                                       cc
cc                                                                    cc
cc    Note: in the current formulation, if usrang = .true., each      cc
cc          downward / upward stream much be accompanied by its       cc
cc          conjugate, such that ang(nze) = 180 - ang(numu-nze+1)=    cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        ng0 - index  of current spectral bin                        cc
cc         l0 - index of state vector sounding                        cc
cc       nlyr - number of computational model layers                  cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       nmom - number of phase function moments used in D/O code     cc
cc       numu - number of output zenith angles used in D/O code       cc
cc       nphi - number of output azimuth angles                       cc
cc      ibcnd - boundy condition flag (0- flx/rad, 1-albedos only)    cc
cc     iunits - output units desired for planck function:             cc
cc           0: unit flux: b(k) = 1.0                                 cc
cc           1: Watts/m**2/cm**-1                                     cc
cc           2: Watts/m**2/micron                                     cc
cc           3: Watts/m**2/nanometer                                  cc
cc           4: Watts/m**2/Angstrom                                   cc
cc           5: Watts/m**2/Hz                                         cc
cc       iref - bidirectional reflectance options                     cc
cc              0 - lambert                                           cc
cc              1 - Hapke's BDR model                                 cc
cc              2 - Breon's BDR model; combination of Li + Roujean    cc
cc              3 - Roujean's BDR model                               cc
cc              4 - Cox and Munk glint model                          cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, use a BRDF is used.              cc
cc     source - include thermal fluxes? (logical: T/F)                cc
cc        umu - emission zenith angle cosines                         cc
cc        phi - emission azimuth angles (degrees)                     cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc       alb0 - surface albedo                                        cc
cc    surf_pr - surface proerties for non-lambertian BRDF             cc
cc      dtauc - layer optical depth                                   cc
cc      ssalb - layer single scattering albedo                        cc
cc       pmom - layer particle phase function                         cc
cc          t - temperature in each atmospheric                       cc
cc      ttemp - temperature of uppermost model layer (space)          cc
cc      btemp - surface temperature                                   cc
cc      temis - emissivity of uppermost model layer (space)           cc
cc      fbeam - intensity of collimated flux at top of atmosphere     cc
cc      fisot - thermal flux at top of atmosphere                     cc
cc       wng0 - wavenumber of bin                                     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc     rfldir - downward direct solar flux at each level for this bin cc
cc   dnflxsrc - downward flux source at each level for this bin       cc
cc   upflxsrc - downward flux source at each level for this bin       cc
cc     radsrc - angle-dependent radiance source at each level for     cc
cc              this bin                                              cc
cc                                                                    cc
ccccccccccccccccccccccccc  s o u r c e _ f c n   ccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang,lamber,source
c
      integer l0,ng0,nlyr,nzdn,nzup,nstr,nmom,numu,nphi,iref,nscat
      integer ibcnd,iunits
      integer ibn,l
c      integer nze,naz
c
      real accur, alb0, fbeam,  phi0, umu0, t(kp), ts
c
      real umu(mxumu),phi(mxphi),surf_pr(4),dtauc(kp),ssalb(kp),
     -       pmom(0:mxmom,kp)
c
      real wng0
c
      real btemp, fisot, temis, ttemp
c
c***    output variables
c
      real rfldir( mxulv ), rfldn( mxulv ), flup( mxulv ),
     -     uu( mxumu, mxulv, mxphi ),albmed( mxumu ),trnmed( mxumu )
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),refflx_b(kp,5,ngrp)
c
c****   variables for flux layer adding method
c
      double precision flx_rb(kp,5),flx_rt(kp,5),
     -                 flx_dd(kp,5),flx_du(kp,5)
c
c****  quantities for radiance calculation: downward and upward fluxes
c      at layer interfaces
c
      double precision flx_ft(kp,5),flx_fb(kp,5)
      double precision flx_dfb(kp,5),flx_uft(kp,5)
      double precision flx_dft(kp,5),flx_ufb(kp,5)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp)
c
c****   variables for radiance layer adding method
c
      double precision rad_rb(mxumu,mxphi,kp,5),rad_rt(mxumu,mxphi,kp,5)
c
c****   internal binned source function variables
c
      double precision upflxsrc(kp,5),dnflxsrc(kp,5),
     -                 radsrc(mxumu,mxphi,kp,5)

c      integer k
c
      ibn = ng0
      l = l0
c
c****     define thermal source function units
c
      iunits = 1
      ibcnd = 0
c
c****            set surface temperature and albedo
c
      btemp = ts
c
c****             set top boundary condition -
c                 no emission or reflection
c
      fisot = 0.0
      ttemp = 0.0
      temis = 1.0
c
c****     call discrete ordinate driver
c
      call do_eq_trn(nlyr,nstr,nmom,numu,nphi,ibcnd,
     -                     iunits,iref,usrang,lamber,source,
     -                     umu,phi,umu0,phi0,accur,
     -                     alb0,surf_pr,dtauc,ssalb,pmom,t,
     -                     ttemp,btemp,temis,fbeam,fisot,wng0,
     -                     flup,rfldir,rfldn,uu,albmed,trnmed,nscat)

c      write(*,'(/,1a,2(1pe14.6),l2))') 'source_fcn: wng0,alb0',
c     -                        wng0,alb0
c      write(*,*) nlyr,nstr,nmom,numu,nphi,ibcnd,
c     -           iunits,iref,usrang,lamber,umu0,phi0,accur,
c     -           ttemp,btemp,temis,fbeam,fisot,(surf_pr(k),k=1,4)
c     - 
c      write(*,'(1a,16(1pe12.4))') 'dtauc ',(dtauc(k),k=1,nlyr)
c      write(*,'(1a,16(1pe12.4))') 'ssalb ',(ssalb(k),k=1,nlyr)
c      write(*,'(1a,16(1pe12.4))') 'pmom  ',(pmom(1,k),k=1,nlyr)
c      write(*,'(1a,16(1pe12.4))') 't     ',(t(k),k=1,nlyr)
c      write(*,'(1a,16(1pe12.4))') 'rfldir',(rfldir(k),k=1,nlyr+1)
c      write(*,'(1a,16(1pe12.4))') 'rfldn ',(rfldn(k),k=1,nlyr+1)
c      write(*,'(1a,16(1pe12.4))') 'flup  ',(flup(k),k=1,nlyr+1)
c
c****   use the flux adding method to solve for the 
c       upward flux at the top of each layer, f+, and the
c       downward flux at the base of each layer, f-.
c
c****   find layer source functions for fluxes
c
      call flx_src(ng0,l0,nlyr,trnflx_b,refflx_b,flx_rt,flx_rb,
     -                   flx_du,flx_dd,rfldn,flup,
     -                   flx_fb,flx_ft,flx_dft,flx_uft,flx_dfb,flx_ufb,
     -                   upflxsrc,dnflxsrc)
c
c****   find radiance source terms
c

      call rad_src(ng0,l0,nlyr,nphi,nzdn,numu,
     -                    trnrad_b,rad_rt,rad_rb,
     -                    flx_ft,flx_fb,flx_dfb,flx_uft,rfldn,flup,uu,
     -                    radsrc)
c
      return
      end
