      subroutine add_rad(nlyr,nzdn,nzup,numu,nphi,trnrad,refrad,brdf,
     -                    rad_src0,ft,fb,urt,drb,uft,dfb,rad0)
c
ccccccccccccccccccccccccccc  a d d _ r a d   ccccccccccccccccccccccccccc
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
cc        fup - f-, upward flux at the top of each homogeneous layer  cc
cc        fdn - f+, downward flux at the base of each layer           cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      radu - upward radiance at nlyr+1 layer boundries.             cc
cc             (flxu(l) refers to the upward flux at the top          cc
cc              of layer l)                                           cc
cc      radd - downward radiance at nlyr+1 layer boundries.           cc
cc             (flxd(l) refers to the downward flux at the bottom     cc
cc              of layer l-1)                                         cc
cc                                                                    cc
ccccccccccccccccccccccccccc  a d d _ r a d   ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,nzdn,nzup,numu,nphi
      integer nlev,k,nze,naz
c
c****   layer flux transmittances, absorptances, and fluxes
c
      double precision trnrad(mxumu,kp),refrad(mxumu,kp)
c
c****   surface bi-directional reflection function
c
      double precision brdf(mxumu,mxphi)
c
c****    radiance source terms
c
      double precision rad_src0(mxumu,mxphi,kp)
c
c****    flux quantitities used on radiance adding
c
      double precision urt(kp),drb(kp),uft(kp),dfb(kp),ft(kp),fb(kp)
c
      double precision fup(mxumu,mxphi,kp),fdn(mxumu,mxphi,kp)
c
c****   variables for flux layer adding method
c
      double precision rb(mxumu,mxphi,kp),rt(mxumu,mxphi,kp)
c
      double precision rad_it(mxumu,mxphi,kp),rad_ib(mxumu,mxphi,kp)
c
c****   output upward and downward radiances
c
      double precision rad0(mxumu,mxphi,kp)
c
      double precision pi
c
      data pi /3.141592654d0 /
c
      nlev = nlyr + 1
c
c****  find layer-integrated reflectivities.  the quantities are:
c      rt(k): reflectance of an inhomogeneous layer extendind from 
c             level 1 to level k
c      rb(k): reflectance of an inhomogeneous layer extendind from 
c             the surface (nlev) to level k
c
      do 1281 naz=1,nphi
c
c****       reflectances of layers extending from the top of the
c           atmosphere to each level for downward streams 
c
          do 1221 nze=1,nzdn
              rt(nze,naz,1) = refrad(nze,1)
              do 1201 k=1,nlyr
                  rt(nze,naz,k+1) = refrad(nze,k) + 
     -                              trnrad(nze,k)*
     -                              rt(nze,naz,k)*urt(k)/pi
1201          continue
c
1221      continue
c
c****       reflectances of layers extending from the surface to 
c           each level of the atmosphere
c
          do 1261 nze=nzup,numu
              rb(nze,naz,nlev) = brdf(nze,naz)
              do 1241 k=1,nlyr
                  rb(nze,naz,nlev-k) = refrad(nze,nlev-k) + 
     -                                 trnrad(nze,nlev-k)*
     -                                 rb(nze,naz,nlev-k+1)*
     -                                 drb(nlev-k+1)/pi
1241          continue
1261      continue
1281  continue
c
c****   rescale source functions with current values
c       note: the upward radiances at the surface (k=nlev) are
c             not scaled since tl(nlev) = 0.0
c
      do 1481 naz=1,nphi
c
c****       scale downward radiance source terms
c
          do 1421 nze=1,nzdn 
              fdn(nze,naz,1) = rad_src0(nze,naz,1)
              do 1401 k=1,nlyr
                  fdn(nze,naz,k+1) = rad_src0(nze,naz,k+1) - 
     -                                trnrad(nze,k)*rt(nze,naz,k)*
     -                                uft(k)/pi
c                  fdn(nze,naz,k+1) = rad_src0(nze,naz,k+1)
1401          continue
1421      continue
c
c****       scale upward radiance source terms
c
          do 1461 nze=nzup,numu
              do 1441 k=1,nlyr        
                  fup(nze,naz,k) = rad_src0(nze,naz,k) - 
     -                              trnrad(nze,k)*rb(nze,naz,k+1)*
     -                              dfb(k+1)/pi
c                  fup(nze,naz,k) = rad_src0(nze,naz,k)
1441          continue
              fup(nze,naz,nlev) = rad_src0(nze,naz,nlev)
1461      continue
1481  continue
c
      do 2081 naz=1,nphi
c
c****       use adding method to find downward radinaces
c           for combined layers.  
c
          do 2021 nze=1,nzdn
c
c              start at top, adding homogeneous layers to the base 
c               of the existing inhomogeneous layer
c
              rad_it(nze,naz,1) = fdn(nze,naz,1)
              do 2001 k=1,nlyr
                rad_it(nze,naz,k+1) = fdn(nze,naz,k+1) + 
     -                                trnrad(nze,k)*
     -                                (rad_it(nze,naz,k) + 
     -                                 rt(nze,naz,k)*uft(k)/pi)
2001          continue
2021      continue
c
          do 2061 nze=nzup,numu
c
c****            use adding method to find upward and downward radiances
c                for combined layers starting at bottom of atmosphere.  
c
              rad_ib(nze,naz,nlev) = fup(nze,naz,nlev)
              do 2041 k=1,nlyr
                  rad_ib(nze,naz,nlev-k) = fup(nze,naz,nlev-k) + 
     -                                     trnrad(nze,nlev-k)*
     -                                     (rad_ib(nze,naz,nlev-k+1) +
     -                                      rb(nze,naz,nlev-k+1)*
     -                                      dfb(nlev-k+1)/pi)
2041          continue
2061      continue
2081  continue
c
c****  find the total upward and downward fluxes at interfaces
c      between inhomogeneous layers.
c
      do 3061 k=1,nlev
          do 3041 naz=1,nphi
c
c****           use adding method to find downward radinaces
c               at layer interfaces  
c
              do 3001 nze=1,nzdn
                  rad0(nze,naz,k) = rad_it(nze,naz,k) + 
     -                              rt(nze,naz,k)*fb(k)/pi
3001          continue
c
c****           use adding method to find downward radinaces
c               at layer interfaces  
c
              do 3021 nze=nzup,numu
                  rad0(nze,naz,k) = rad_ib(nze,naz,k) + 
     -                              rb(nze,naz,k)*ft(k)/pi
3021          continue 
3041      continue
3061  continue
c
      return
      end
     
