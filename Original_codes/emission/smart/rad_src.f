      subroutine rad_src(ng0,l0,nlyr,nphi,nzdn,numu,
     -                    trnrad_b,rad_rt,rad_rb,
     -                    flx_ft,flx_fb,flx_dfb,flx_uft,rfldn,flup,uu,
     -                    radsrc)
c
cccccccccccccccccccccccccccc  r a d _ s r c   cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine used a simplified flux adding method to find    cc
cc    the effective radiance source term in each layer, given the     cc
cc    radiance and flux field from DISORT.                            cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        ng0 - index  of current spectral bin                        cc
cc       nlyr - number of computational model layers                  cc
cc   trnrad_b - scattering transmittance of an isolate layer for bin  cc
cc     rad_rt - radiance reflectance of layer over this level         cc
cc     rad_rb - radiance reflectance of layer below this level        cc
cc         uu - angle dependent radiance at interface of 2 layers     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       radsrc: radiance source for an isolated layer                cc
cc                                                                    cc
cccccccccccccccccccccccccccc  r a d _ s r c   cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,ng0,l0,nphi,nzdn,numu
      integer k,ibn,l,nze,naz,nlev
c
c****   DISORT output
c
      real rfldn( mxulv ), flup( mxulv ),uu( mxumu, mxulv, mxphi )
c
c****   binned layer transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp)
c
c****   variables for radiance layer adding method
c
      double precision rad_rb(mxumu,mxphi,kp,5),rad_rt(mxumu,mxphi,kp,5)
c
c****  quantities for radiance calculation: downward and upward fluxes
c      at layer interfaces
c
      double precision flx_ft(kp,5),flx_fb(kp,5)
      double precision flx_dfb(kp,5),flx_uft(kp,5)
c
c****   internal binned source function variables
c
      double precision radsrc(mxumu,mxphi,kp,5)
      double precision rad_ib(kp),rad_it(kp)
c
      double precision pi
c
c****   define pi
c
      data pi /3.141592654d0 /
c
      ibn = ng0
      l = l0
      nlev = nlyr+1
c
c****    use the modified flux adding method to find radiance sources
c
      do 1241 naz=1,nphi
c
          do 1221 nze=1,nzdn
c
              do 1001 k=1,nlev
                  rad_it(k) = uu(nze,k,naz) - 
     -                        rad_rt(nze,naz,k,l)*flx_fb(k,l)/pi
                  rad_ib(k) = uu(numu-nze+1,k,naz) - 
     -                        rad_rb(numu-nze+1,naz,k,l)*flx_ft(k,l)/pi
1001          continue
c
c****            find downward radiance sources, starting at top
c
              radsrc(nze,naz,1,l) = rad_it(1)
c
              do 1041 k=1,nlyr
c
c****               find downward radiances 
c
                  radsrc(nze,naz,k+1,l) = rad_it(k+1) - 
     -                   trnrad_b(nze,k,l,ibn)*
     -                   (rad_it(k) + rad_rt(nze,naz,k,l)*
     -                    flx_uft(k,l)/pi)
1041          continue
c
              do 1061 k=1,nlyr
c
c****               find the upward radiances
c
                  radsrc(numu-nze+1,naz,k,l) = rad_ib(k) - 
     -                        trnrad_b(numu-nze+1,k,l,ibn)*
     -                        (rad_ib(k+1) + 
     -                         rad_rb(numu-nze+1,naz,k+1,l)*
     -                         flx_dfb(k+1,l)/pi)
1061          continue
c
              radsrc(numu-nze+1,naz,nlev,l) = rad_ib(nlev)
c
1221      continue
1241  continue
c
      return
      end
