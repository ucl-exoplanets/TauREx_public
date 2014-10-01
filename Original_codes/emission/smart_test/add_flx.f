      subroutine add_flx(nlyr,tl,rl,fup,fdn,urt,drb,uft,dfb,ft,fb,
     -                   flxu,flxd)
c
ccccccccccccccccccccccccccc  a d d _ f l x   ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this program uses the layer adding method for finding fluxes    cc
cc    in a vertically inhomogeneous scattering, absorbing atmosphere. cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       nlyr - number of computational model layers                  cc
cc         rl - diffuse flux reflectance of each layer                cc
cc         tl - diffuse flux transmittance of each homogeneous layer  cc
cc        fup - f+, upward flux at the top of each homogeneous layer  cc
cc        fdn - f-, downward flux at the base of each layer           cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      flxu - upward flux at nlyr+1 layer boundries.                 cc
cc             (flxu(l) refers to the upward flux at the top          cc
cc              of layer l)                                           cc
cc      flxd - downward flux at nlyr+1 layer boundries.               cc
cc             (flxd(l) refers to the downward flux at the bottom     cc
cc              of layer l-1)                                         cc
cc                                                                    cc
ccccccccccccccccccccccccccc  a d d _ f l x   ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr
      integer nlev,k,l
c
c****   layer flux transmittances, absorptances, and fluxes
c
      double precision tl(kp),rl(kp),fup(kp),fdn(kp)
c
c****   variables for flux layer adding method
c
      double precision rb(kp),rt(kp),dd(kp),du(kp)
c
c****    quantitities used on radiance adding
c
      double precision urt(kp),drb(kp),uft(kp),dfb(kp),ft(kp),fb(kp)
c
c****   output
c
      double precision flxu(kp),flxd(kp)
c
      nlev = nlyr + 1
c
c****  find layer-integrated reflectivities.  the quantities are:
c      rt(k): reflectance of an inhomogeneous layer extendind from 
c             level 1 to level k
c      rb(k): reflectance of an inhomogeneous layer extendind from 
c             the surface (nlev) to level k
c
      rt(1) = 0.0
      rb(nlev) = rl(nlev)
      urt(nlev) = 0.0d0
      drb(1) = 0.0d0
      do 1201 k=1,nlyr
          dd(k) = 1.0d0/(1.0d0 - rt(k)*rl(k))
          urt(k) = tl(k)*dd(k)
          rt(k+1) = rl(k) + tl(k)**2*rt(k)*dd(k)
          du(nlev-k) = 1.0d0/(1.0d0 - rl(nlev-k)*rb(nlev-k+1))
          drb(nlev-k+1) = tl(nlev-k)*du(nlev-k)
          rb(nlev-k) = rl(nlev-k) + tl(nlev-k)**2*
     -                       rb(nlev-k+1)*du(nlev-k)
1201  continue
c
c****   use adding method to find upward and downward fluxes
c       for combined layers.  start at top, adding homogeneous 
c       layers to the base of the existing inhomogeneous layer
c
      ft(1) = 0.0d0
      do 1401 k=1,nlyr
          ft(k+1) = fdn(k+1) + tl(k)*(ft(k) + rt(k)*fup(k))*dd(k) 
1401  continue
c
c****   use adding method to find upward and downward fluxes
c       for combined layers.  start at bottom.
c
      fb(nlev) = fup(nlev)
      do 1601 l=1,nlyr
          k = nlev - l
          fb(k) = fup(k) + tl(k)*(fb(k+1) + rb(k+1)*fdn(k+1))*du(k)
1601  continue
c
c****  find the total upward and downward fluxes at interfaces
c      between inhomogeneous layers.
c
      flxd(1) = fdn(1)
      do 2001 k=2,nlev
          flxd(k) = (ft(k) + rt(k)*fb(k))/(1.0d0 - rt(k)*rb(k))
2001  continue 
c
      do 2021 k=1,nlev
          flxu(k) = fb(k) + rb(k)*flxd(k)
2021  continue
c
c****   find the variables dft and ufb for radiance adding
c
      do 3001 k=1,nlyr
          uft(k) = (fup(k) + rl(k)*ft(k))*dd(k)
          dfb(k+1) = (fdn(k+1) + rl(k)*fb(k+1))*du(k)
3001  continue
      write(*,'(/,1a)') 'in add_flx: '
      write(*,'(1a,16(1pe12.4))') 'rl:    ',(rl(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'tl:    ',(tl(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'dd:    ',(dd(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'du:    ',(du(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'rt:    ',(rt(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'rb:    ',(rb(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'ft:    ',(ft(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'fb:    ',(fb(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'fdn:   ',(fdn(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'fup:   ',(fup(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'uft:   ',(uft(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'dfb:   ',(dfb(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'flxd:  ',(flxd(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4))') 'flxu:  ',(flxu(k),k=1,nlev)
c
      return
      end
     
