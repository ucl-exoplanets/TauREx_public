      subroutine layer_flx(usrang,lamber,l0,nlyr0,ng0,
     -                     nstr,numu,nzdn,nzup,nmom,nphi,iref,
     -                     umu,phi,dtau_b,copi0_b,pmom_b,
     -                     umu0,phi0,accur)
c
cccccccccccccccccccccccccc  l a y e r _ f l x   cccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the radiances and fluxes layer by      cc
cc    layer to test the adding method.                                cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       nlyr0: number of layers on the atmosphere.                   cc
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
cc    trn_rad: layer transmittance for each radiance stream.          cc
cc    abs_rad: layer absorptance for each rediance stream             cc
cc                                                                    cc
cccccccccccccccccccccccccc  l a y e r _ f l x   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang,lamber,prnt(7),onlyfl,plank,usrtau
      save prnt,onlyfl,plank,usrtau
c
      integer l0,nlyr0,nstr,numu,ng0,nmom
c
      integer l,k,m,mom,nze,ibn,naz
c
      integer ibcnd,maxcly,maxmom,maxphi,maxulv,maxumu,nlyr
      save ibcnd,maxcly,maxmom,maxphi,maxulv,maxumu,nlyr
c
      integer nphi,ntau,iref,ifirst,iunits
      save ntau,ifirst,iunits
c
      integer nlev,nzdn,nzup
      save nlev
c
      real accur,alb,btemp,fbeam,fisot,temis,ttemp,umu0,phi0,wvnm,
     -     umu(mxumu)
      save btemp,fbeam,fisot,temis,ttemp,wvnm
c
      real phi(mxphi),surf_pr(4),temper(0:mxrad),utau(mxulv)
c
      real dtauc(mxrad),ssalb(mxrad),pmom(0:mxmom,mxrad)
c
c***   define a dummy cosine of zenith angle variable for call
c      to DISORT for the flux transmission calculation.
c
c*****   state vector variables.
c
      real dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     pmom_b(0:mxmom,kp,5,ngrp)
c
c***    disort output variables
c
      real rfldir( mxulv ), rfldn( mxulv ), flup( mxulv ),
     -       dfdt( mxulv ), uavg( mxulv ),uu( mxumu, mxulv, mxphi ), 
     -       albmed( mxumu ),trnmed( mxumu )
c
      real flxu(kp,5),flxd(kp,5),dir(kp),rad(mxumu,mxphi,kp)
      real deltau,deltaui(kp)
      double precision dflxddx(kp),dflxudx(kp)
c
      data ifirst / 0 /
c
c      write(*,*) 'layer_flx',l0,nlyr0,ng0,
c     -                     nstr,numu,nzdn,nzup,nmom,nphi,iref
      if(ifirst .eq. 0) then
c
c****      initialzie disort varibles
c
        ifirst = 1
        nlev = nlyr0+1
c
c****      set logicals
c
        onlyfl = .false.
        plank = .false.
        lamber = .true.
        usrtau = .false.
        wvnm = 1.0
        fbeam = 1.0
        fisot = 0.0
        iunits = 0
c
c****    set variable for isotropic illumniation from top
c
        ibcnd = 0      
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
c
c****     set number of layers to 1
c
        nlyr = 1
c
      endif   
c
      ibn = ng0
      l = l0
c
c****     set surface albedos to zero
c
      alb = 0.0
      surf_pr(1) = alb
      do 1021 m=2,4
          surf_pr(m) = 0.0
1021  continue
c
c****    find the layer transittances and reflectances for each layer
c
      dir(1) = umu0*fbeam
      flxd(1,l) = 0.0
      do 2441 k=1,nlyr0
c
c****    find the layer transittances and reflectances for each sounding
c
          dtauc(1) = dtau_b(k,l,ibn)
          ssalb(1) = 1.0 - copi0_b(k,l,ibn)
          pmom(0,1) = 1.0
          do 2001 mom=1,nmom
             pmom(mom,1) =  pmom_b(mom,k,l,ibn)
2001      continue
c
          fbeam=1.0
c
          call disort( nlyr, dtauc, ssalb, nmom, pmom, temper, 
     -                    wvnm, usrtau, ntau, utau, nstr, iunits,
     -                    usrang, numu, umu, nphi, phi, ibcnd, fbeam, 
     -                    umu0, phi0, fisot, lamber, iref, 
     -                    surf_pr, alb, btemp,ttemp, temis, 
     -                    plank,onlyfl, accur, prnt, 
     -                    maxcly, maxulv, maxumu, maxphi, maxmom, 
     -                    rfldir, rfldn, flup, dfdt, uavg, uu, 
     -                    albmed, trnmed )
c
          flxu(k,l) = flup(1)
          flxd(k+1,l) = rfldn(2)
          dir(k+1) = rfldir(2)
          do 2241 naz=1,nphi
              do 2201 nze=1,nzdn
                  rad(nze,naz,k+1) = uu(nze,2,naz)
2201          continue
              do 2221 nze=nzup,numu
                  rad(nze,naz,k) = uu(nze,1,naz)
2221          continue
2241      continue
c
2441  continue
      k = nlev
          do 2281 naz=1,nphi
              do 2261 nze=nzup,numu
                  rad(nze,naz,k) = uu(nze,1,naz)
2261          continue
2281      continue
      flxu(nlev,l) = flup(2)
      if(ibn .eq. 1) then
      write(*,'(/,1a,i5)') 'from layer_flx: l = ',l
c      write(*,'(1a,16(1pe12.4))') 'dir:        ',(dir(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dtau_b:     ',
     - (dtau_b(k,l,ibn),k=1,nlev-1)
      write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'flxd:       ',
     - (flxd(k,l),k=1,nlev)
      write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'flxu:       ',
     - (flxu(k,l),k=1,nlev)
c      write(*,*)
c      do 1 nze=1,numu
c      naz = 1
c1     write(*,'(1a,i2,16(1pe12.4))') 'layer rad ',
c     -  nze,((rad(nze,naz,k)),k=1,nlev)
      if(l .gt. 1) then
        do 3001 k=1,nlev-1
          deltau = (dtau_b(k,l,ibn) - dtau_b(k,1,ibn))
          if(deltau .ne. 0.0) then
            deltaui(k) = 1./deltau
            dflxddx(k+1) = (flxd(k+1,l) - flxd(k+1,1))*deltaui(k)
            dflxudx(k) = (flxu(k,l) - flxu(k,1))*deltaui(k)
          else
            deltaui(k) = 0.0
            dflxddx(k) = 0.0
            dflxudx(k) = 0.0
          endif
3001    continue
      write(*,'(/,1a,i5)') 'from layer_flx partials: l = ',l
      write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dflxddx:   ',
     - (dflxddx(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'dflxudx:   ',
     - (dflxudx(k),k=1,nlev)
      write(*,'(1a,16(1pe12.4),4(/16(1pe12.4)))') 'deltaui:   ',
     - (deltaui(k),k=1,nlev)
      endif
      endif
c
      return
      end
