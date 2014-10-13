      subroutine rad_ref(l0,nlyr,ng0,numu,nzdn,nzup,nphi,
     -                    trnrad_b,refrad_b,brdf_b,flx_urt,flx_drb,
     -                    rad_rb,rad_rt)
c
cccccccccccccccccccccccccccc  r a d _ r e f   cccccccccccccccccccccccccc
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
cc      rad_rt(k): reflectance of an inhomogeneous layer extendind    cc
cc             from level 1 to level k+1 (layers 1 - k)               cc
cc      rad_rb(k): reflectance of an inhomogeneous layer extendind    cc
cc             from the surface (nlev) to level k                     cc
cc                                                                    cc
cccccccccccccccccccccccccccc  r a d _ r e f   cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer l0,nlyr,ng0,numu,nzdn,nzup,nphi
c
      integer l,k,ibn,nze,naz,nlev
c
c**** flux transmittances, reflectances, and fluxes
c
      double precision flx_urt(kp,5),flx_drb(kp,5)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp),
     -                 refrad_b(mxumu,kp,5,ngrp),
     -                 brdf_b(mxumu,mxphi,5,ngrp)
c
      double precision rad_rb(mxumu,mxphi,kp,5),rad_rt(mxumu,mxphi,kp,5)
c
      double precision pi,dflx,uflx
c
      data pi /3.141592654d0 /
c
c****  find layer-integrated reflectivities.  the quantities are:
c      rad_rt(k): reflectance of an inhomogeneous layer extendind from 
c                 level 1 to level k+1 (layers 1 - k)
c      rad_dd(k): denominator of rad_rt(k)
c      rad_rb(k): reflectance of an inhomogeneous layer extendind from 
c                 the surface (nlev) to level k
c      rad_du(k): denominator of rad_rb(k)
c
      l = l0
      ibn = ng0
      nlev = nlyr+1
c
c****    set reflectance at the top of the atmosphere for adding down
c
      do 1021 naz=1,nphi
          do 1001 nze=1,nzdn
              rad_rt(nze,naz,1,l) = refrad_b(nze,1,l,ibn)
1001      continue
1021  continue
c
      do 1241 k=1,nlyr
          uflx = flx_urt(k,l)/pi
          do 1221 naz=1,nphi
              do 1201 nze=1,nzdn
                  rad_rt(nze,naz,k+1,l) = refrad_b(nze,k,l,ibn) + 
     -                                    trnrad_b(nze,k,l,ibn)*
     -                                    rad_rt(nze,naz,k,l)*uflx
1201          continue
1221      continue
1241  continue
c      if(ng0 .eq. 1) then
c      do 1 nze=1,nzdn
c      write(*,'(/,1a,3i5)') 'in rad_ref: nze,naz,l',nze,l
c      write(*,'(1a,16(1pe12.4))') 'trnrad_b:  ',
c     -         (trnrad_b(nze,k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'refrad_b:  ',
c     -         (refrad_b(nze,k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_urt/pi:',
c     -         (flx_urt(k,l)/pi,k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'rad_rt:    ',
c     -         (rad_rt(nze,1,k,l),k=1,nlev)
c1     continue
c      endif
c
      do 2021 naz=1,nphi
          do 2001 nze=nzup,numu
              rad_rb(nze,naz,nlev,l) = brdf_b(nze,naz,l,ibn)
2001      continue
2021  continue
c
      do 2241 k=1,nlyr
          dflx = flx_drb(nlev-k+1,l)/pi
          do 2221 naz=1,nphi
              do 2201 nze=nzup,numu
                  rad_rb(nze,naz,nlev-k,l) = 
     -                                  refrad_b(nze,nlev-k,l,ibn) + 
     -                                  trnrad_b(nze,nlev-k,l,ibn)*
     -                                  rad_rb(nze,naz,nlev-k+1,l)*dflx
2201          continue
2221      continue
2241  continue
c
c****  rad_rb(nze,nlev,l) is defined in rad_src,
c      pi*uu(numu-nze+1,k,naz)/rfldn(k)
c      
c      if(ng0 .eq. 1) then
c      do 2 nze=nzup,numu
c      write(*,'(/,1a,3i5)') 'in rad_ref: nze, l:',nze,l
c      write(*,'(1a,16(1pe12.4))') 'trnrad_b:  ',
c     -         (trnrad_b(nze,k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'refrad_b:  ',
c     -         (refrad_b(nze,k,l,ibn),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'flx_drb/pi:',
c     -         (flx_drb(k,l)/pi,k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'rad_rb:    ',
c     -         (rad_rb(nze,1,k,l),k=1,nlev)
c2     continue
c      endif
c
      return
      end
