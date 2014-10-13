      subroutine layer_ref(l0,nlyr,ng0,n_rad,alb0,
     -                    trnflx_b,refflx_b,absflx_b,
     -                    dd,du,rb,rt)
c
cccccccccccccccccccccccccc  l a y e r _ r e f   cccccccccccccccccccccccc
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
cc      rt(k): reflectance of an inhomogeneous layer extendind from   cc
cc             level 1 to level k+1 (layers 1 - k)                    cc
cc      rb(k): reflectance of an inhomogeneous layer extendind from   cc
cc             the surface (nlev) to level k                          cc
cc                                                                    cc
cccccccccccccccccccccccccc  l a y e r _ r e f   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer l0,nlyr,ng0
      integer n_rad(5,ngrp)
c
      integer l,k,ibn
c
      integer nlev
c
c      to DISORT for the flux transmission calculation.
c
      real alb0
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),
     -                 refflx_b(kp,5,ngrp),
     -                 absflx_b(kp,5,ngrp)
c
      double precision rb(kp,5,ngrp),rt(kp,5,ngrp),
     -                 dd(kp,5,ngrp),du(kp,5,ngrp)
c
c****  find layer-integrated reflectivities.  the quantities are:
c      rt(k): reflectance of an inhomogeneous layer extendind from 
c             level 1 to level k+1 (layers 1 - k)
c      dd(k): denominator of rt(k)
c      rb(k): reflectance of an inhomogeneous layer extendind from 
c             the surface (nlev) to level k
c      du(k): denominator of rb(k)
c
      nlev = nlyr+1
      refflx_b(nlev,l,ibn) = alb0
      trnflx_b(nlev,l,ibn) = 0.0d0
      absflx_b(nlev,l,ibn) = 1.0d0 - alb0
      rt(1,l,ibn) = refflx_b(1,l,ibn)
      rb(nlev,l,ibn) = refflx_b(nlev,l,ibn)
c
      do 1001 k=1,nlyr
          dd(k,l,ibn) = 1.0d0/(1.0d0 - rt(k,l,ibn)*
     -                      refflx_b(k+1,l,ibn))
          rt(k+1,l,ibn) = refflx_b(k+1,l,ibn) + 
     -                        trnflx_b(k+1,l,ibn)**2*
     -                        rt(k,l,ibn)*dd(k,l,ibn)
          du(nlev-k,l,ibn) = 1.0d0/(1.0d0 - 
     -                           refflx_b(nlev-k,l,ibn)*
     -                           rb(nlev+1-k,l,ibn))
          rb(nlev-k,l,ibn) = refflx_b(nlev-k,l,ibn) + 
     -                           trnflx_b(nlev-k,l,ibn)**2*
     -                           rb(nlev+1-k,l,ibn)*du(nlev-k,l,ibn)
1001  continue
c      write(*,'(/,1a,2i5)') 'in layer_trn: ibn,l ',ibn,l
c      write(*,'(1a,16(1pe12.4))') 'refflx:',(refflx_b(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'trnflx:',(trnflx_b(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'dd:    ',(dd(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'du:    ',(du(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'rt:    ',(rt(k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'rb:    ',(rb(k,l,ibn),k=1,nlev)
c
      return
      end
