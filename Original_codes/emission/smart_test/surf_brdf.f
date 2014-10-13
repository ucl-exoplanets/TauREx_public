      subroutine surf_brdf(usrang,lamber,
     -                     l0,ng0,nz0,nstr,numu,nphi,
     -                     nmom,n_rad,nref,umu,phi,dalb_b_i,
     -                     alb_b,sur_b,ws,phiw,umu0,phi0,
     -                     brdf_b,dbrdfdx)
c
cccccccccccccccccccccccccc  s u r f _ b r d f   cccccccccccccccccccccccc
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
cc        nstr: number of computational radiance streams.             cc
cc         ng0: spectral bin number                                   cc
cc      dtau_b: layer optical depth for bin, ibn0.                    cc
cc     copi0_b: layer single scattering albedo for bin, ibn0.         cc
cc      pmom_b: layer particle phase function for bin, ibn0.          cc
cc         umu: zenith angle of each stream.                          cc
cc         gwt: gaussian weight for each stream.                      cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    brdf_b: surface bi-directional reflection distribution function cc
cc   dbrdfdx: partial derivative of brdf_b with respect to albedo     cc
cc                                                                    cc
cccccccccccccccccccccccccc  s u r f _ b r d f   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang,lamber
      logical prnt(7),onlyfl,plank,usrtau
      save    prnt,onlyfl,plank,usrtau
c
      integer l0,ng0,nz0,nstr,numu,nphi,nref,nmom
      integer n_rad(5,ngrp)
c
      integer l,m,mom,nze,naz,ibn,nz,nr
c
      integer ibcnd,maxcly,maxmom,maxphi,maxulv,maxumu,nlyr
      save    ibcnd,maxcly,maxmom,maxphi,maxulv,maxumu,nlyr
c
      integer ntau,iref,ifirst,iunits
      save    ntau,iref,ifirst,iunits
c
      real umu(mxumu),phi(mxphi),umu0,phi0
c
c****    local variables used in disort
c
      real accur,alb,btemp,fbeam,fisot,temis,ttemp,wvnm
      save accur,alb,btemp,fbeam,fisot,temis,ttemp,wvnm
c
      real surf_pr(4),temper(0:mxrad),utau(mxulv),pi
      save surf_pr,temper,utau,pi
c
      real dtauc(mxrad),ssalb(mxrad),pmom(0:mxmom,mxrad)
c
c*****   state vector variables.
c
      real sur_b(4,5,ngrp),alb_b(nsol,5,ngrp),ws,phiw,
     -     dalb_b_i(nsol,ngrp)
c
c***    disort output variables
c
      real rfldir( mxulv ), rfldn( mxulv ), flup( mxulv ),
     -       dfdt( mxulv ), uavg( mxulv ),uu( mxumu, mxulv, mxphi ), 
     -       albmed( mxumu ),trnmed( mxumu )
c
c****   binned layer radiance transmittances and absorptances
c
      double precision brdf_b(mxumu,mxphi,5,ngrp),
     -                 dbrdfdx(mxumu,mxphi,ngrp)
c
      data ifirst / 0 /
c
      if(ifirst .eq. 0) then
c
c****      initialzie disort varibles
c
        ifirst = 1
        iunits = 1
        ibcnd = 0
c
c****      set logicals
c
        usrang = .false.
        onlyfl = .false.
        plank = .false.
        usrtau = .false.
        iref = 0
        wvnm = 1.0
        pi = acos(-1.0)
c
c****         set dimensions of arrays
c
        maxcly = 2
        maxulv = 2
        maxumu = mxumu
        maxphi = mxphi
        maxmom = mxmom
c
c****     turn off all discr_ord model internal print flags
c
        do 1001 m=1,7
            prnt(m) = .false.
1001    continue
c
        ntau = 0
        accur = 1.e-4
c
c****     set number of layers to 1
c
        nlyr = 1
c
      endif   
c
      ibn = ng0
      l = l0
      nz = nz0
c
      if(l .eq. 1 .or. (l .eq. 5 .and. n_rad(l,ibn) .ne. 0)) then
c
c****     set surface temperature and atmospheric temperature
c
        btemp = 296.0
        temper(0) = 296.0
        temper(1) = 296.0
        wvnm = 1000.
c
c****             set top boundary condition - no emission or reflection
c
        fisot = 0.0
        ttemp = 0.0
        temis = 1.0
        fbeam = 1.0
c
        dtauc(1) = 0.0
        ssalb(1) = 0.0
        pmom(0,1) = 1.0
        do 2001 mom=1,nmom
           pmom(mom,1) =  0.0
2001    continue        
c
c****    define the surface optical properties
c
        alb = alb_b(nz,l,ibn)
c
c****    set the surface optical properties for disort
c
        do 2021 nr=1,nref
            surf_pr(nr) = sur_b(nr,l,ibn)
2021    continue
c  
        if(iref .eq. 4) then
c
c****      set wind speed and direction for Cox/Munk model
c
          surf_pr(3) = ws
          surf_pr(4) = phiw
        endif
c
        call disort( nlyr, dtauc, ssalb, nmom, pmom, temper, 
     $                   wvnm, usrtau, ntau, utau, nstr, iunits,
     &                   usrang, numu, umu, nphi, phi, ibcnd, fbeam, 
     &                   umu0, phi0, fisot, lamber, iref, 
     &                   surf_pr, alb, btemp, ttemp,temis, 
     &                   plank, onlyfl, accur, prnt, 
     &                   maxcly,maxulv, maxumu, maxphi, maxmom, 
     &                   rfldir, rfldn,flup, dfdt, uavg, uu, 
     &                   albmed, trnmed )
c
        if(l .eq. 1) then
c
c****         load brdf's: note-values are independent of atmospheric
c             optical depth, single scattering albedo, and particle 
c             phase function, so values of l=1,4 are the same.
c
          do 2461 m=1,4
              do 2441 naz=1,nphi
                  do 2421 nze=1,numu
                      brdf_b(nze,naz,m,ibn) = pi*uu(nze,2,naz)/
     -                                         rfldir(2)
2421              continue
2441          continue
2461      continue
c
        else
c
c****        find the brdf jacobians.
c
          do 2821 naz=1,nphi
              do 2801 nze=1,numu
                  brdf_b(nze,naz,l,ibn) = pi*uu(nze,2,naz)/rfldir(2)
                  dbrdfdx(nze,naz,ibn) = (brdf_b(nze,naz,l,ibn) - 
     -             brdf_b(nze,naz,1,ibn))*dalb_b_i(nz,ibn)
2801          continue
2821      continue
c
        endif
c
      endif
      
      write(*,'(/,1a,3i5)') 'surf_brdf: nz,l,ibn: ',nz,l,ibn
      write(*,'(1a,16(1pe12.4))') 'rfldir    ',rfldir(nlyr+1)
      write(*,'(1a,16(1pe12.4))') 'rfldn     ',rfldn(nlyr+1)
      write(*,'(1a,16(1pe12.4))') 'uu(1)     ',(uu(nze,1,1),nze=1,numu)
      write(*,'(1a,16(1pe12.4))') 'uu(2)     ',(uu(nze,2,1),nze=1,numu)
      write(*,'(1a,16(1pe12.4))') 'brdf_b: ',
     -                        (brdf_b(nze,1,l,ibn),nze=1,numu)
c
      return
      end
