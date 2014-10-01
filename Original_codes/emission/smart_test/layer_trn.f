      subroutine layer_trn(usrang0,l0,ng0,nlyr0,
     -                     nstr,numu0,nmom,nphi,n_rad,
     -                     dtau_b,copi0_b,pmom_b,dx_b_i,
     -                     umu,phi,umu_f,gwt_f,
     -                     trnflx_b,dtrnflxdx,refflx_b,drefflxdx,
     -                     absflx_b,dabsflxdx,trnrad_b,dtrnraddx,
     -                     refrad_b,drefraddx,absrad_b,dabsraddx)
c
cccccccccccccccccccccccccc  l a y e r _ t r n   cccccccccccccccccccccccc
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
cc   trnrad_b: layer transmittance for each radiance stream.          cc
cc   refrad_b: layer reflectance for each radiance stream.            cc
cc   absrad_b: layer absorptance for each rediance stream             cc
cc                                                                    cc
cccccccccccccccccccccccccc  l a y e r _ t r n   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical usrang0
      logical usrang,lamber,prnt(7),onlyfl,plank,usrtau
      save usrang,lamber,prnt,onlyfl,plank,usrtau
c
      integer l0,nlyr0,nstr,numu0,nphi,ng0
      integer nmom,n_rad(5,ngrp)
c
      integer l,k,m,mom,nze,naz,ibn
c
      integer ibcnd,maxcly,maxmom,maxphi,maxulv,maxumu,nlyr
      save    ibcnd,maxcly,maxmom,maxphi,maxulv,maxumu,nlyr
c
      integer ntau,iref,ifirst,iunits
      save    ntau,iref,ifirst,iunits
c
      integer numu,numuf
c
      real umu_f(mxumu), gwt_f(mxumu),umu(mxumu),phi(mxphi)
      real umu1(mxumu),phi1(mxphi)
c
c****    local variables used in disort
c
      real accur,alb,btemp,fbeam,fisot,temis,ttemp,umu0,phi0,wvnm
      save accur,alb,btemp,fbeam,fisot,temis,ttemp,umu0,phi0,wvnm
c
      real surf_pr(4),temper(0:mxrad),utau(mxulv)
      save surf_pr,temper,utau
c
      real dtauc(mxrad),ssalb(mxrad),pmom(0:mxmom,mxrad)
c
c*****   state vector variables.
c
      real dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     pmom_b(0:mxmom,kp,5,ngrp),dx_b_i(kp,5,ngrp)
c
c***    disort output variables
c
      real rfldir( mxulv ), rfldn( mxulv ), flup( mxulv ),
     -       dfdt( mxulv ), uavg( mxulv ),uu( mxumu, mxulv, mxphi ), 
     -       albmed( mxumu ),trnmed( mxumu )
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),dtrnflxdx(kp,3,ngrp),
     -                 refflx_b(kp,5,ngrp),drefflxdx(kp,3,ngrp),
     -                 absflx_b(kp,5,ngrp),dabsflxdx(kp,3,ngrp)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision trnrad_b(mxumu,kp,5,ngrp),
     -                 dtrnraddx(mxumu,kp,3,ngrp),
     -                 refrad_b(mxumu,kp,5,ngrp),
     -                 drefraddx(mxumu,kp,3,ngrp),
     -                 absrad_b(mxumu,kp,5,ngrp),
     -                 dabsraddx(mxumu,kp,3,ngrp)
c
      double precision dtau
c
      data ifirst / 0 /
c
      if(ifirst .eq. 0) then
c
c****      initialzie disort varibles
c
        ifirst = 1
        numu = numu0
        iunits = 1
c
c****      set logicals
c
        usrang = .false.
        onlyfl = .false.
        plank = .false.
        lamber = .true.
        usrtau = .false.
        iref = 0
        wvnm = 1.0
c
c****    set variable for isotropic illumniation from top
c
        ibcnd = 1      
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
c****     set surface albedos to zero
c
        alb = 0.0
        do 1021 m=1,4
            surf_pr(m) = 0.0
1021    continue
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
c
      if(l .lt. 5 .and. n_rad(l,ibn) .ne. 0) then
c
c****    find the layer transittances and reflectances for each layer
c
c      write(*,'(/,1a)') 'above disort layer loop: '
c      write(*,'(1a,16(1pe12.4))') 'umu: ',(umu(nze),nze=1,nstr)
c      write(*,'(1a,16(1pe12.4))') 'umu_f:',(umu_f(nze),nze=1,nstr)
        do 2461 k=1,nlyr0
c
c****       find layer transittances and reflectances for each sounding
c
            dtauc(1) = dtau_b(k,l,ibn)
            ssalb(1) = 1.0 - copi0_b(k,l,ibn)
            if(ssalb(1) .lt. 0.0) ssalb(1) = 0.0
            pmom(0,1) = 1.0
            do 2001 mom=1,nmom
               pmom(mom,1) =  pmom_b(mom,k,l,ibn)
2001        continue
c
            numuf = nstr
            usrang = .false.
c
c****         call disort to find layer transmissivities and albedos.
c             note: in this mode, it changes the values of umu, so a 
c             dummy variable should be used.
c
            call disort( nlyr, dtauc, ssalb, nmom, pmom, temper, 
     -                    wvnm, usrtau, ntau, utau, nstr, iunits,
     -                    usrang, numuf, umu1, nphi, phi1, ibcnd, fbeam,
     -                    umu0, phi0, fisot, lamber, iref, 
     -                    surf_pr, alb, btemp,ttemp, temis, 
     -                    plank,onlyfl, accur, prnt, 
     -                    maxcly, maxulv, maxumu, maxphi, maxmom, 
     -                    rfldir, rfldn, flup, dfdt, uavg, uu, 
     -                    albmed, trnmed )
c
c****           repack radiance variables
c
            do 2201 nze=1,nstr/2
                trnmed(nze+nstr/2) = trnmed(nze)
                albmed(nze+nstr/2) = albmed(nze)
2201        continue
            do 2221 nze=1,nstr/2
                trnmed(nze) = trnmed(nstr-nze+1)
                albmed(nze) = albmed(nstr-nze+1)
2221        continue
c
c****         find the flux transmittance and absorptance by 
c             integrating transmittance over zenith angle using 
c             gaussian quadrature
c
            trnflx_b(k,l,ibn) = 0.0
            refflx_b(k,l,ibn) = 0.0
            absflx_b(k,l,ibn) = 0.0
            do 2241 nze=1,nstr/2
                trnflx_b(k,l,ibn) = trnflx_b(k,l,ibn) + 
     -                            2.0*gwt_f(nze)*umu_f(nze)*trnmed(nze)
                refflx_b(k,l,ibn) = refflx_b(k,l,ibn) + 
     -                            2.0*gwt_f(nze)*umu_f(nze)*albmed(nze)
                absflx_b(k,l,ibn) = absflx_b(k,l,ibn) + 
     -                            2.0*gwt_f(nze)*umu_f(nze)*
     -                            (1.0 - (trnmed(nze) + albmed(nze)))
2241        continue
c
            if(usrang0) then
c
c****         find layer transmissivities and reflectivities for
c             the user-specified zenith angles
c
              numu = numu0
              usrang = .true.
              do 2301 nze=1,numu
                  umu1(nze) = umu(nze)
2301          continue
              do 2321 naz=1,nphi
                  phi1(naz) = phi(nze)
2321          continue
c
c****         call disort to find layer transmissivities and albedos.
c             note: in this mode, it changes the values of umu, so a 
c             dummy variable should be used.
c
              call disort( nlyr, dtauc, ssalb, nmom, pmom, temper, 
     $                   wvnm, usrtau, ntau, utau, nstr, iunits,
     &                   usrang, numu, umu1, nphi, phi1, ibcnd, fbeam, 
     &                   umu0, phi0, fisot, lamber, iref, 
     &                   surf_pr, alb, btemp, ttemp,temis, 
     &                   plank, onlyfl, accur, prnt, 
     &                   maxcly,maxulv, maxumu, maxphi, maxmom, 
     &                   rfldir, rfldn,flup, dfdt, uavg, uu, 
     &                   albmed, trnmed )
c
            endif
c
c****         load layer radiance tranmittances
c
            do 2401 nze=1,numu0
                dtau = dtau_b(k,l,ibn)/abs(umu(nze))
                trnrad_b(nze,k,l,ibn) = dexp(-dtau)
                refrad_b(nze,k,l,ibn) = albmed(nze)
                absrad_b(nze,k,l,ibn) = (1.0d0 - (trnmed(nze) + 
     -                                   albmed(nze)))
2401        continue
2461    continue
c      write(*,'(/,1a)') 'below disort layer loop: '
c      write(*,'(1a,16(1pe12.4))') 'umu: ',(umu(nze),nze=1,nstr)
c      write(*,'(1a,16(1pe12.4))') 'umu_f:',(umu_f(nze),nze=1,nstr)
c      do 1 nze=1,numu
c       write(*,'(/,1a,2i5,1pe12.4)') 
c     -      'layer_trn: l,nze,umu(nze)',l,nze,umu(nze)
c1       write(*,'(1a,16(1pe12.4))') 'trnrad_b',
c     - (trnrad_b(nze,k,l,ibn),k=1,nlyr0)
c
      else
c
c****      either l=5 (perturbed surface albedos) or nrad=0.  
c          in either case, the layer transmittance, reflectance, 
c          and absorptance are set to values for nominal case
c
        do 2621 k=1,nlyr0
            trnflx_b(k,5,ibn) = trnflx_b(k,1,ibn)
            refflx_b(k,5,ibn) = refflx_b(k,1,ibn)
            absflx_b(k,5,ibn) = absflx_b(k,1,ibn)
c
            do 2601 nze=1,numu
                trnrad_b(nze,k,5,ibn) = trnrad_b(nze,k,1,ibn)
                refrad_b(nze,k,5,ibn) = refrad_b(nze,k,1,ibn)
                absrad_b(nze,k,5,ibn) = absrad_b(nze,k,1,ibn)
2601        continue
2621    continue
c
      endif
c
      if(l .gt. 1 .and. l .lt. 5) then
c
c****        find the transmittance and absorptance jacobians.
c            these are only needed for optical depth (l=2), 
c            single scattering albedo (l=3), and scattering
c            phase function (l=4)
c
        do 2821 k=1,nlyr0
            dtrnflxdx(k,l-1,ibn) = (trnflx_b(k,l,ibn) - 
     -                              trnflx_b(k,1,ibn))*dx_b_i(k,l,ibn)
            drefflxdx(k,l-1,ibn) = (refflx_b(k,l,ibn) - 
     -                              refflx_b(k,1,ibn))*dx_b_i(k,l,ibn)
            dabsflxdx(k,l-1,ibn) = (absflx_b(k,l,ibn) - 
     -                              absflx_b(k,1,ibn))*dx_b_i(k,l,ibn)
            do 2801 nze=1,numu
                dtrnraddx(nze,k,l-1,ibn) = (trnrad_b(nze,k,l,ibn) - 
     -                                      trnrad_b(nze,k,1,ibn))*
     -                                      dx_b_i(k,l,ibn)
                drefraddx(nze,k,l-1,ibn) = (refrad_b(nze,k,l,ibn) - 
     -                                      refrad_b(nze,k,1,ibn))*
     -                                      dx_b_i(k,l,ibn)
                dabsraddx(nze,k,l-1,ibn) = (absrad_b(nze,k,l,ibn) - 
     -                                      absrad_b(nze,k,1,ibn))*
     -                                      dx_b_i(k,l,ibn)
2801        continue
c
2821    continue
c
      endif
c
      return
      end
