      subroutine flx_ref(l0,nlyr,ng0,alb0,trnflx_b,refflx_b,absflx_b,
     -                    dx_b_i,flx_dd,flx_du,flx_urt,flx_drb,
     -                    flx_rb,flx_rt)
c
cccccccccccccccccccccccccccc  f l x _ r e f   cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the effective layer radiance and flux  cc
cc    transmittances, reflectances and and absorptances for each      cc
cc    spectral bin.  These quantities are used in the spectral        cc
cc    mapping and jacobian calculations.                              cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        nlyr: number of layers on the atmosphere.                   cc
cc         ng0: spectral bin number                                   cc
cc      dtau_b: layer optical depth for bin, ibn0.                    cc
cc     copi0_b: layer single scattering albedo for bin, ibn0.         cc
cc      pmom_b: layer particle phase function for bin, ibn0.          cc
cc         umu: zenith angle of each stream.                          cc
cc         gwt: gaussian weight for each stream.                      cc
cc    trn_rad: layer transmittance for each radiance stream.          cc
cc    abs_rad: layer absorptance for each rediance stream             cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      flx_rt(k): reflectance of an inhomogeneous layer extending    cc
cc             from level 1 to level k+1 (layers 1 - k)               cc
cc      flx_rb(k): reflectance of an inhomogeneous layer extending    cc
cc             from the surface (nlev) to level k                     cc
cc                                                                    cc
cccccccccccccccccccccccccccc  f l x _ r e f   cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer l0,nlyr,ng0
c
      integer l,k,ibn
c
      integer nlev
c
c      to DISORT for the flux transmission calculation.
c
      real alb0,dx_b_i(kp,5,ngrp)
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),
     -                 refflx_b(kp,5,ngrp),
     -                 absflx_b(kp,5,ngrp)
c
      double precision flx_rb(kp,5),flx_rt(kp,5),
     -                 flx_dd(kp,5),flx_du(kp,5),
     -                 flx_urt(kp,5),flx_drb(kp,5)
c
c****  find layer-integrated reflectivities.  the quantities are:
c      flx_rt(k): reflectance at the bottom of an inhomogeneous layer
c             extending from level 1 to level k+1 (layers 1 - k)
c      flx_dd(k): denominator of flx_rt(k)
c      flx_rb(k): reflectance of the top of an inhomogeneous layer 
c             extending from the surface (nlev) to level k
c      flx_du(k): denominator of flx_rb(k)
c
      l = l0
      ibn = ng0
      nlev = nlyr + 1
c
      flx_rt(1,l) = 0.0
      flx_rb(nlev,l) = refflx_b(nlev,l,ibn)
c
      do 1001 k=1,nlyr
c
c****       find the layer reflectances adding from top
c
          flx_dd(k,l) = 1.0d0/(1.0d0 - flx_rt(k,l)*refflx_b(k,l,ibn))
          flx_urt(k,l) = trnflx_b(k,l,ibn)*(1.0 + 
     -                   refflx_b(k,l,ibn)*flx_rt(k,l)*flx_dd(k,l))
          flx_rt(k+1,l) = refflx_b(k,l,ibn) + 
     -                    trnflx_b(k,l,ibn)**2*flx_rt(k,l)*flx_dd(k,l)
c
c***        find layer reflectances, adding up from the surface
c
          flx_du(nlev-k,l) = 1.0d0/(1.0d0 - 
     -                       refflx_b(nlev-k,l,ibn)*flx_rb(nlev-k+1,l))
          flx_drb(nlev-k+1,l) = trnflx_b(nlev-k,l,ibn)*(1.0 + 
     -                          refflx_b(nlev-k,l,ibn)*
     -                          flx_rb(nlev-k+1,l)*flx_du(nlev-k,l))
          flx_rb(nlev-k,l) = refflx_b(nlev-k,l,ibn) + 
     -                           trnflx_b(nlev-k,l,ibn)**2*
     -                           flx_rb(nlev-k+1,l)*
     -                           flx_du(nlev-k,l)
1001  continue
c
c      if(ibn .eq. 1) then
c      write(*,'(/,1a,2i5,1pe12.4)') 'in flx_ref: ibn,l,alb0 ',ibn,l,alb0
c      write(*,'(1a,16(1pe12.4))') 'refflx: ',
c     -                             (refflx_b(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'trnflx: ',
c     -                             (trnflx_b(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_dd: ',(flx_dd(k,l),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_urt:',(flx_urt(k,l),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_du: ',(flx_du(k,l),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_drb:',(flx_drb(k,l),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_rt: ',(flx_rt(k,l),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_rb: ',(flx_rb(k,l),k=1,nlev)
c      endif
c
      return
      end
