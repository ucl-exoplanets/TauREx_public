      subroutine load_optics(lamber,ng0,nlyr,nza,nstate,istate,
     -                         n_rad,nref,iref,nmomgrp,
     -                         umu0,tauerr,pi0err,phferr,surferr,
     -                         taumn,surfgrp,albgrp,
     -                         taugrp,pi0grp,ggrp,pmomgrp,
     -                         alb_b,sur_b,phiw,ws,
     -                         dtau_b,copi0_b,g_b,pmom_b,
     -                         dx_b_i,dalb_b_i)
c
ccccccccccccccccccccccc  l o a d _ o p t i c s   ccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this routine loads optical depths, single scattering albedos,   cc
cc    surface reflectance and surface pressures into arrays from use  cc
cc    in the routine sm_eq_trn.                                       cc
cc    NOTE: this version implements a perturbation that is            cc
cc          related to the actual range in each bin.                  cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        ng0 - index  of current spectral bin                        cc
cc       nlyr - number of computational model layers                  cc
cc     nstate - number of variable elements in the state vecror       cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc              1 - surface pressure                                  cc
cc              2 - surface/atmospheric temperature                   cc
cc              3 - gas absorption coeffient                          cc
cc              4 - cloud/aerosol optical depth                       cc
cc              5 - surface albedo                                    cc
cc      n_rad - radiance calculations flag for each bin               cc
c       n_rad = 1 : radiances for nominal state structure
c             = 2 : radiances for perturbed optical depths
c             = 3 : radiances for perturbed single scattering albedos
c             = 4 : radiances for perturbed phase functions
c             = 5 : radiances for perturbed surface reflectances
cc       nref - number of surface optical properties specified at     cc
cc              each wavelength.                                      cc
cc     tauerr - optical depth relative binning error (0. to ~0.8)     cc
cc     pi0err - co-single scattering albedo absolute binning error    cc
cc    surferr - surface optical property binning error                cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      sur_b - surface albedo for this bin                           cc
cc     dtau_b - differential optical depth in each layer for this bin cc
cc    copi0_b - single scattering co-albedo in each layer in this bin cc
cc        g_b - asymmetry parameter in each layer for this bin        cc
cc                                                                    cc
ccccccccccccccccccccccc  l o a d _ o p t i c s   ccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lamber
      integer ng0,nlyr,nza
      integer nref,iref
      integer nstate,istate(nex)
      integer n_rad(5,ngrp)
      integer k,n,ibn,ist,nr,mom,nz
      integer nmomgrp(ngrp)
c
      real taugrp(kp,ngrp,3),pi0grp(kp,ngrp,3),
     -     ggrp(kp,ngrp,3),pmomgrp(0:mxmom,kp,ngrp,3),
     -     surfgrp(ngrp,4,3),albgrp(ngrp,nsol,3)
c
      real tauerr,pi0err,phferr,surferr
      real taumn,umu0(nsol),phiw,ws
c
c*****   state vector variables.
c
      real sur_b(4,5,ngrp),alb_b(nsol,5,ngrp),
     -     dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp),pmom_b(0:mxmom,kp,5,ngrp),
     -     dx_b_i(kp,5,ngrp),dalb_b_i(nsol,ngrp)
c
c****   internal variables
c
      real tau_range,ss_range,g_range,dref_i(nsol),dref0,dsurf
      real taumin,ssamin,gmin,surfmin,deltau_b,delpi0_b,delg_b
c
      real deltaumin
      real dref,surf_pr(4)
c
c****    set group counter
c
      ibn = ng0
c
c****    define a minimum change for delta-tau
c
c       deltaumin = 1.e-5
       deltaumin = 1.e-3
c
c****   set optical depths for the nomical case and for each variable 
c       component of the optical depth structure and set 
c       the radiance calculation flag, n_rad, where
c       n_rad = 1 : radiances for nominal state structure
c             = 2 : radiances for perturbed optical depths
c             = 3 : radiances for perturbed single scattering albedos
c             = 4 : radiances for perturbed phase functions
c             = 5 : radiances for perturbed surface reflectances
c
c****       define the column-integrated optcal depth, single scatting
c           albedo, and particle phase function ranges in each bin
c
      tau_range = 0.0 
      ss_range = 0.0
      g_range = 0.0
c
      do 1001 k=1,nlyr
c
c****       define the column-integrated optcal depth, single scatting
c           albedo, and particle phase function ranges
c
          tau_range = tau_range + taugrp(k,ibn,3) - taugrp(k,ibn,2)
          ss_range = ss_range + pi0grp(k,ibn,3) - pi0grp(k,ibn,2)
          g_range = g_range + ggrp(k,ibn,3) - ggrp(k,ibn,2)
1001  continue
c
      if(istate(1) .eq. 1) then
c
c****      add the difference between the nominal surface layer 
c          optical depth and the optical depth for the layer with 
c          a perturbed surface pressure.
c
          tau_range = tau_range + abs(taugrp(nlyr+1,ibn,1) - 
     -                                taugrp(nlyr,ibn,1))
          ss_range =  ss_range + abs(pi0grp(nlyr+1,ibn,1)  - 
     -                               pi0grp(nlyr,ibn,1))
          g_range = g_range + abs(ggrp(nlyr+1,ibn,1) - 
     -                            ggrp(nlyr,ibn,1))
c
      endif
c
c****     n o m i n a l    c a s e
c
      n_rad(1,ibn) = 1
c
c****    define the optical properties for the nominal state vector
c     
      do 1021 k=1,nlyr
          dtau_b(k,1,ibn) = taugrp(k,ibn,1)
          copi0_b(k,1,ibn) = pi0grp(k,ibn,1)
          g_b(k,1,ibn) = ggrp(k,ibn,1)
          pmom_b(0,k,1,ibn) = 1.0
          do 1011 mom=1,nmomgrp(ibn)
              pmom_b(mom,k,1,ibn) =  pmomgrp(mom,k,ibn,1)
1011      continue
1021  continue
c
c****   define the nominal surface albedo
c
      do 1041 nr=1,nref
          sur_b(nr,1,ibn) = surfgrp(ibn,nr,1)
1041  continue
      do 1061 nz=1,nza
          alb_b(nz,1,ibn) = albgrp(ibn,nz,1)
1061  continue
c
c****     o p t i c a l   p r o p e r t y    p e r t u r b a t i o n s
c
c***      check to see if optical depth or single scattering albedo
c         are variable parts of the state structure
c
      ist = 0
      do 1301 n=1,nstate
          if(iabs(istate(n)) .ge. 1 .and. iabs(istate(n)) .le. 4) 
     -     ist = 1
1301  continue
c
c****   determine if optical depth varies across the bin 
c
      taumin = 1.0e-3*tauerr*tauerr + taumn
c      if(ibn .eq. 1) write(*,'(/,1a,2(1pe12.4))') 
c     - 'load_optics: taumin,tau_range: ',
c     - taumin,tau_range
      if(tau_range .gt. taumin .or. ist .eq. 1) then
c
        n_rad(2,ibn) = 1
c
        do 1351 k=1,nlyr
c
c****         note: add at least deltaumin to this quantity to 
c                   ensure that radiance difference is no-zero, and
c                   that the jacobian denominator is not too small
c
            dtau_b(k,2,ibn) = deltaumin +
     -                       (1.0 + tauerr*tauerr)*taugrp(k,ibn,1) +      
     -                        0.25*(taugrp(k,ibn,3) - taugrp(k,ibn,2))
            copi0_b(k,2,ibn) = pi0grp(k,ibn,1)
            g_b(k,2,ibn) = ggrp(k,ibn,1)
            pmom_b(0,k,2,ibn) = 1.0
            do 1321 mom=1,nmomgrp(ibn)
                pmom_b(mom,k,2,ibn) =  pmomgrp(mom,k,ibn,1)
1321        continue
c
1351    continue
c      if(ibn .eq. 1) then
c     - write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dtau_b(1)   ',
c     - (dtau_b(k,1,ibn),k=1,nlyr)
c     - write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dtau_b(2)   ',
c     - (dtau_b(k,2,ibn),k=1,nlyr)
c       endif
c
        do 1341 nr=1,nref
            sur_b(nr,2,ibn) = surfgrp(ibn,nr,1)
1341    continue
        do 1361 nz=1,nza
            alb_b(nz,2,ibn) = albgrp(ibn,nz,1)
1361    continue
c
      else
c
        n_rad(2,ibn) = 0
c
      endif
c
c****     s i n g l e   s c a t t e r i n g   a l b e d o 
c
c****   check single scattering co-albedo range in this bin
c
      ssamin = 1.0e-3*pi0err*pi0err + 1.0e-5
c
      if(ss_range .gt. ssamin .or. ist .eq. 1) then
c
        n_rad(3,ibn) = 1
c
        do 1411 k=1,nlyr
            dtau_b(k,3,ibn) = taugrp(k,ibn,1)
c
            copi0_b(k,3,ibn) = 0.01 + 
     -                         (1.0 + pi0err*pi0err)*pi0grp(k,ibn,1) + 
     -                         0.25*(pi0grp(k,ibn,3) - pi0grp(k,ibn,2))
c
            if(copi0_b(k,3,ibn) .gt. 1.00)
     -         copi0_b(k,3,ibn) = -0.01 + 
     -                        pi0grp(k,ibn,1)/(1.0 + pi0err*pi0err) -
     -                        0.25*(pi0grp(k,ibn,3) - pi0grp(k,ibn,2))
c
            g_b(k,3,ibn) = ggrp(k,ibn,1)
            pmom_b(0,k,3,ibn) = 1.0
            do 1401 mom=1,nmomgrp(ibn)
                pmom_b(mom,k,3,ibn) =  pmomgrp(mom,k,ibn,1)
1401        continue
1411    continue
c
        do 1421 nr=1,nref
            sur_b(nr,3,ibn) = surfgrp(ibn,nr,1)
1421    continue
        do 1441 nz=1,nza
            alb_b(nz,3,ibn) = albgrp(ibn,nz,1)
1441    continue
c
      else
c
        n_rad(3,ibn) = 0
c     
      endif 
c
c****     s c a t t e r i n g   p h a s e    f u n c t i o n
c
c*****      NOTE: phase functions cannot be arbitrarily multiplied by 
c           a constant, or they will lose their normaization.  
c           replace the values by the min and max binned values  
c
c****   check phase function range in this bin
c
      gmin = 1.0e-3*phferr*phferr + 1.0e-5
      if(g_range .gt. gmin .or. ist .eq. 1) then
c
        n_rad(4,ibn) = 1
c
        do 1541 k=1,nlyr
            dtau_b(k,4,ibn) = taugrp(k,ibn,1)
            copi0_b(k,4,ibn) = pi0grp(k,ibn,1)
            g_b(k,4,ibn) = ggrp(k,ibn,3) 
            pmom_b(0,k,4,ibn) = 1.0
            do 1521 mom=1,nmomgrp(ibn)
                pmom_b(mom,k,4,ibn) =  pmomgrp(mom,k,ibn,3) 
1521        continue
1541    continue
c
        do 1561 nr=1,nref
            sur_b(nr,4,ibn) = surfgrp(ibn,nr,1)
1561    continue    
        do 1581 nz=1,nza
            alb_b(nz,4,ibn) = albgrp(ibn,nz,1)
1581    continue
c
      else
c
        n_rad(4,ibn) = 0
c     
      endif 
c
c****    s u r f a c e    o p t i c a l   p r o p e r t i e s 
c
c****   in this version, the surface reflectance is always 
c       a variable part of the state structure if any variable changes
c       because it affects the upward radiances at the surface
c
      ist = 0
      do 1601 n=1,nstate
          if(iabs(istate(n)) .eq. 5) ist = 1
1601  continue
c
c*****  set the optical depths and single scattering albedos to nominal
c       values
c
      do 1621 k=1,nlyr
          dtau_b(k,5,ibn) = taugrp(k,ibn,1)
          copi0_b(k,5,ibn) = pi0grp(k,ibn,1)
          g_b(k,5,ibn) = ggrp(k,ibn,1)
          pmom_b(0,k,5,ibn) =  0.0
          do 1611 mom=1,nmomgrp(ibn)
              pmom_b(mom,k,5,ibn) =  pmomgrp(mom,k,ibn,1)
1611      continue
1621  continue
c
c****   determine the surface reflectance variations in this bin
c
      n_rad(5,ibn) = 0
      surfmin = 1.0e-3*surferr*surferr + 1.0e-5
      do 1661 nr=1,nref
c
          dsurf = surfgrp(ibn,nr,3) - surfgrp(ibn,nr,2)
c
          if(dsurf .gt. surfmin .or. ist .eq. 1) then 
c
c****        radiances must be calculated for 2 surface reflectances
c
            n_rad(5,ibn) = 1
c
            sur_b(nr,5,ibn) = 0.25*dsurf + 0.01 + 
     -                        (1.0 + surferr*surferr)*surfgrp(ibn,nr,1) 
            if(sur_b(nr,5,ibn) .gt. 1.0) 
     -         sur_b(nr,5,ibn) = 0.99*surfgrp(ibn,nr,1) - 0.25*dsurf
c
          else
c
            sur_b(nr,5,ibn) = surfgrp(ibn,nr,1)
c     
          endif
c
1661  continue
c
      if(lamber) then
c
c****     set the perturbed albedo equal to the first perturbed sur_b
c         value
c
        do 1671 nz=1,nza
            alb_b(nz,5,ibn) = sur_b(1,5,ibn)
1671    continue
c
      else
c
c****   use the disort routine dref to calculate the value of alb
c       that corresponds to this difference in the parameters
c
        do 1701 nr=1,nref
            surf_pr(nr) = sur_b(nr,5,n)
1701    continue
        if(iref .eq. 4) then
c
c****      set wind speed and direction for Cox/Munk model
c
          surf_pr(3) = ws
          surf_pr(4) = phiw
        endif
        do 1721 nz=1,nza
            alb_b(nz,5,ibn) = dref(umu0(nz),surf_pr,iref)
1721    continue
c
      endif
c
c****    compute denomenators for flux and radiance partial derivatives
c        optical depth, single scattering co-albedo, and surface albedo
c
      do 1741 nz=1,nza      
          dref0 = alb_b(nz,5,ibn) - alb_b(nz,1,ibn)
          if(abs(dref0) .gt. surfmin) then
            dref_i(nz) = 1.0/dref0
          else
            dref_i(nz) = 0.0
          endif
1741  continue
c
c****   define the amplitudes of the perturbations
c
      do 2021 k=1,nlyr
c
c****       nominal calse
c
         dx_b_i(k,1,ibn) = 0.0
c
c****       optical depth: n_rad(2,ibn) = 1
c
          deltau_b = dtau_b(k,2,ibn) - dtau_b(k,1,ibn)
          if(deltau_b .gt.
     -       (1.0e-3*tauerr*tauerr*dtau_b(k,1,ibn)) + taumn) then
            dx_b_i(k,2,ibn) = 1.0/deltau_b
          else
            dx_b_i(k,2,ibn) = 0.0
          endif    
c
c****       single scattering albedo: n_rad(3,ibn) = 1
c
          delpi0_b = copi0_b(k,3,ibn) - copi0_b(k,1,ibn)
          if(abs(delpi0_b) .gt. ssamin) then
            dx_b_i(k,3,ibn) = 1.0/delpi0_b
          else
            dx_b_i(k,3,ibn) = 0.0
          endif
c
c*****        phase function: n_rad(4,ibn) = 1
c
          delg_b = g_b(k,4,ibn) - g_b(k,1,ibn)
          if(abs(delg_b) .gt. gmin) then
            dx_b_i(k,4,ibn) = 1.0/delg_b
          else
            dx_b_i(k,4,ibn) = 0.0
          endif
c
2021  continue
c
      do 2101 nz=1,nza
          dalb_b_i(nz,ibn) = dref_i(nz)
2101  continue
c
      return
      end
      
