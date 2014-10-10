      subroutine flx_src(ng0,l0,nlyr,trnflx_b,refflx_b,flx_rt,flx_rb,
     -                   flx_du,flx_dd,rfldn,flup,
     -                   flx_fb,flx_ft,flx_dft,flx_uft,flx_dfb,flx_ufb,
     -                   upflxsrc,dnflxsrc)
c
ccccccccccccccccccccccccccc  f l x _ s r c   ccccccccccccccccccccccccccc
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
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       upflxsrc: f- (Crisp 1986 eq 7) upward flux at base           cc
cc                     of an isolated layer.                          cc
cc       dnflxsrc: f+ (Crisp 1986 eq 8) downward flux at top          cc
cc                    of isolated  layer.                             cc
cc                                                                    cc
ccccccccccccccccccccccccccc  f l x _ s r c   ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,ng0,l0
      integer nlev,k,l,ibn
c
c****   DISORT output
c
      real rfldn( mxulv ), flup( mxulv )
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),refflx_b(kp,5,ngrp)
c
c****   variables for radiance layer adding method
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
c****   output variables
c
      double precision upflxsrc(kp,5),dnflxsrc(kp,5)
c
      nlev = nlyr+1
      l = l0
      ibn = ng0
c
c****     use the flux adding method to solve for flx_ft, the
c         downward diffuse flux at the base of each 
c         inhomogeneous layer extending from layers 1 to k
c         and flx_fb, the upward flux at the top of each
c         inhomogeneous layer extending from the surface to k 
c
c****   find the downward flux at the base of layer k 
c                  (level k+1) (dflux) eq 22 Crisp 1986)
c
      flx_ft(1,l) = 0.0d0
      do 1001 k=2,nlev
          flx_ft(k,l) = rfldn(k) - flx_rt(k,l)*flup(k)
1001  continue
c
c****   find the upward flux at the top of an inhomogeneous
c       layer extending from level k to the surface (uflux)
c       (simplifed version of Eq 23 of Crisp 1986)
c                  
      do 1201 k=1,nlev
          flx_fb(k,l) = flup(k) - flx_rb(k,l)*rfldn(k)
1201  continue
c
c****   find the fluxes at the top of each homogeneous layer 
c       (upflxsrc) and the bottom of each layer
c
      dnflxsrc(1,l) = 0.0
      do 2001 k=1,nlyr
c
c            solve for the upward flux at the top of each
c            homogenous layer,upflxsrc (ufl in delted).
c            To derive this equation, I started with Eqs. 18 and 
c            21 of Crisp 1986, solved for f+ and f-, respectively,
c            and then substituted the value of f+ into the expression
c            for f-
c
c            upflxsrc = f+
c
          upflxsrc(k,l) = (flx_fb(k,l) - trnflx_b(k,l,ibn)*
     -                    (flx_fb(k+1,l) + flx_rb(k+1,l)*
     -                     (flx_ft(k+1,l) - trnflx_b(k,l,ibn)*
     -                     flx_ft(k,l)*flx_dd(k,l)))*flx_du(k,l))/
     -                     (1.0d0 - trnflx_b(k,l,ibn)**2*
     -                     flx_rb(k+1,l)*flx_rt(k,l)*
     -                     flx_dd(k,l)*flx_du(k,l))
c
c             solve for the downward flux at the base of each
c                   homogenous layer,f- = dnflxsrc
c
          dnflxsrc(k+1,l) = flx_ft(k+1,l) - trnflx_b(k,l,ibn)*
     -                      (flx_ft(k,l) + flx_rt(k,l)*upflxsrc(k,l))*
     -                      flx_dd(k,l)
c
2001  continue
      upflxsrc(nlev,l) = flx_fb(nlev,l)
c
c****   find the variables flx_dft and flx_ufb for radiance adding
c
      do 3001 k=1,nlev
          flx_uft(k,l) = (upflxsrc(k,l) + 
     -                    refflx_b(k,l,ibn)*flx_ft(k,l))*flx_dd(k,l)
          flx_dft(k,l) = (flx_ft(k,l) + flx_rt(k,l)*upflxsrc(k,l))*
     -                    flx_dd(k,l)
3001  continue
c
      flx_dfb(1,l) = 0.0
c
      do 3021 k=1,nlyr
          flx_dfb(k+1,l) = (dnflxsrc(k+1,l) + refflx_b(k,l,ibn)*
     -                      flx_fb(k+1,l))*flx_du(k,l)
          flx_ufb(k,l) = flx_fb(k,l) + flx_rb(k,l)*flx_dfb(k,l)
3021  continue
      flx_ufb(nlev,l) = flx_fb(nlev,l) + flx_rb(nlev,l)*flx_dfb(nlev,l)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       write(*,'(1a,16(1pe12.4))') 'rfldn   ',
c     - (rfldn(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'rfldn   ',(rfldn(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flup    ',(flup(k),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flx_ft  ',(flx_ft(k,l),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flx_fb  ',(flx_fb(k,l),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flx_rt  ',(flx_rt(k,l),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flx_rb  ',(flx_rb(k,l),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flx_rb  ',(flx_dd(k,l),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'flx_rb  ',(flx_du(k,l),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'refflx_b',
c     -                             (refflx_b(k,l,ibn),k=1,nlev)
c       write(*,'(1a,16(1pe12.4))') 'trnflx_b',
c     -                             (trnflx_b(k,l,ibn),k=1,nlev)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
