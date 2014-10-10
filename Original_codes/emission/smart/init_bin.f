      subroutine init_bin(nlyr,nref,nza,ngroup,iwngrp,ismterr,
     -                     itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt,
     -                     levtau1,nmomgrp,dnugrp,wngrp,surfgrp,
     -                     albgrp,taugrp,pi0grp,ggrp,pmomgrp,
     -                     dtau_b,copi0_b,g_b,pmom_b,sur_b,alb_b,
     -                     dx_b_i,dalb_b_i)
c
ccccccccccccccccccccccc    i n i t _ b i n     ccccccccccccccccccccccccc
cc                                                                    cc
cc     p u r p o s e :                                                cc
cc                                                                    cc
cc    this subroutine initializes values set in the binning routine.  cc
cc                                                                    cc
cc     i n p u t                                                      cc
cc                                                                    cc
cc       nlyr - number of layers in nominal atmosphere                cc
cc       nref - number of surface optical properties specified at     cc
cc              each wavelength.                                      cc
cc     ngroup - number of spectral bins (0 to ngrp)                   cc
cc     iwngrp - index of spectral bin                                 cc
cc    ismterr - error flag: ngroup > ngrp                             cc
cc   itau1cnt - number of bins rejected by poor fit at tau=1          cc
cc    itaucnt - number of bins rejected by poor fit away from tau=1   cc
cc    ipi0cnt - number of bins rejected by poor pi0 fit               cc
cc      igcnt - number of bins rejected by poor g fit                 cc
cc    isurcnt - number of bins rejected by poor albedo fit            cc
cc    levtau1 - index of tau=1 level                                  cc
cc    nmomgrp - maximum number of phase function moments in bin       cc
cc     taugrp - mean, min, and max optical depth in each bin          cc
cc     piogrp - mean, min, and max single scattering albedo in  bin   cc
cc    pmomgrp - mean, min, and max phase function in each bin         cc
cc       ggrp - mean, min, and max asymmetry parmeter in each bin     cc
cc     dtau_b - perturbed differential optical depth in each bin      cc
cc      pi0_b - perturbed single scattering albedo in each bin        cc
cc     pmom_b - perturbed scattering phase function in each bin       cc
cc        g_b - perturbed scattering asymmetery parametr in each bin  cc
cc      sur_b - perturbed surface albedo in each bin                  cc
cc                                                                    cc
cc     o u t p u t                                                    cc
cc                                                                    cc
cc    zeroed variables.                                               cc
cc                                                                    cc
ccccccccccccccccccccccc    i n i t _ b i n     ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,ismterr,itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt
      integer levtau1(ngrp,2),nmomgrp(ngrp),iwngrp,ngroup,nref,nza
      integer k,l,m,n,nr,nz
c
      real taugrp(kp,ngrp,3),pi0grp(kp,ngrp,3),
     -     ggrp(kp,ngrp,3),pmomgrp(0:mxmom,kp,ngrp,3),
     -     surfgrp(ngrp,4,3),albgrp(ngrp,nsol,3)
c
      real sur_b(4,5,ngrp),alb_b(nsol,5,ngrp),
     -     dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp),pmom_b(0:mxmom,kp,5,ngrp),
     -     dx_b_i(kp,5,ngrp),dalb_b_i(nsol,ngrp)
c
      double precision wngrp(ngrp,3),dnugrp(ngrp)
c
c****   initialize spectral binning parameters
c
      ismterr = 0
      itau1cnt = 0
      itaucnt = 0
      ipi0cnt = 0
      igcnt = 0
      isurcnt = 0
      iwngrp = 0
      ngroup = 0
c
      do 3281 n=1,ngrp
          levtau1(n,1) = nlyr
          levtau1(n,2) = nlyr
          nmomgrp(n) = 0
          dnugrp(n) = 0.0d0
c
          do 2801 l=1,3
              wngrp(n,l) = 0.0d0
              do 2001 nr=1,nref
                 surfgrp(n,nr,l) = 0.0
2001          continue
              do 2201 nz=1,nza
                 albgrp(n,nz,l) = 0.0
2201          continue
              do 2601 k=1,kp
                  taugrp(k,n,l) = 0.0
                  pi0grp(k,n,l) = 0.0
                  ggrp(k,n,l) = 0.0
                  pmomgrp(0,k,n,l) = 1.0
                  do 2401 m=1,mxmom
                      pmomgrp(m,k,n,l) = 0.0
2401              continue
2601          continue
2801      continue
c
          do 3061 l=1,5
              do 3021 k=1,kp
                  dtau_b(k,l,n) = 0.0
                  copi0_b(k,l,n) = 0.0
                  g_b(k,l,n) = 0.0
                  pmom_b(0,k,l,n) = 1.0
                  dx_b_i(k,l,n) = 0.0
                  do 3001 m=1,mxmom
                      pmom_b(m,k,l,n) = 0.0
3001              continue
3021          continue
3061      continue
c
          do 3101 nz=1,nza
              dalb_b_i(nz,n) = 0.0
3101      continue
c
          do 3241 l=1,5
              do 3201 nr=1,nref
                  sur_b(nr,l,n) = 0.0
3201          continue
              do 3221 nz=1,nza
                  alb_b(nz,l,n) = 0.0
3221          continue
3241      continue
3281  continue
c
      return
      end
