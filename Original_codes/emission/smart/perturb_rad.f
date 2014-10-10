      subroutine perturb_rad(lsolar,lplanck,lamber,usrang,ni0,nz0,
     -         nlyr,iref,nlout,levout,k_out,ibin,nphi,nstr,
     -         numu,nzup,nzdn,npd,npd1,ipd,igs,modepd,
     -         wn_io,dvi,umu_f,gwt_f,umu,phi,umu0nz,phi0nz,pd_pert,dx,
     -         alb,ts,t,p,dp_dp,rmix,dtauaer,tau_ext,tau_sca,
     -         g_sca,dtauex,co_pi0,g,alb_b,dtau_b,copi0_b,g_b,
     -         trnflx_b,dtrnflxdx,refflx_b,drefflxdx,absflx_b,
     -         dabsflxdx,refrad_b,drefraddx,absrad_b,dabsraddx,
     -         brdf_b,dbrdfdx,dnsflxsrc_b,upsflxsrc_b,sradsrc_b,
     -         dntflxsrc_b,uptflxsrc_b,tradsrc_b,
     -         ddntflxdx,duptflxdx,dtraddx,ddnsflxdx,dupsflxdx,dsraddx,
     -         dalb0,deltau0,delpi00,delg0,bb0,bb_sur0,
     -         bb_flx_up0,bb_flx_dn0,bb_rad0,dn_s_src0,up_s_src0,
     -         rad_s_src0,dn_t_src0,up_t_src0,rad_t_src0,
     -         trndir0,trnflx0,refflx0,absflx0,trnrad0,refrad0,absrad0,
     -         trndir1_i,trnflx1_i,refflx1_i,absflx1_i,
     -         dtrndir1_i,dtrnflx1_i,drefflx1_i,dabsflx1_i,
     -         dn_s_src1,up_s_src1,dup_s_src1,ddn_s_src1,
     -         dn_t_src1,up_t_src1,dup_t_src1,ddn_t_src1,
     -         sol_rad1,dsol_rad1,th_rad1,dth_rad1)
c
cccccccccccccccccccccccc  p e r t u r b _ r a d   cccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine finds evaluates the partial derivatives of the  cc
cc    radiances with respect to each variable member of the state     cc
cc    vector.                                                         cc
cc                                                                    cc
cc    NOTE: this version has been modified to produce layer dependent cc
cc    flux source functions and their partial derivatives, rather     cc
cc    than the corresponding flux values.                             cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        nz0 - index of this solar zenith angle (1 to nsol)          cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc       nlyr - number of computational model layers                  cc
cc       nphi - number of output azimuth angles                       cc
cc       numu - number of zenith angles in input file                 cc
cc     uptflx - wn-dependent upward thermal flux                      cc
cc     dntflx - wn-dependent downward thermal flux                    cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc     upsflx - wn-dependent upward solar flux                        cc
cc     dnsflx - wn-dependent downward diffuse + direct solar flux     cc
cc    dirsflx - wn-dependent downward direct solar flux               cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cccccccccccccccccccccccc  p e r t u r b _ r a d   cccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar,lamber,usrang
c
      integer ni0,nz0,nlyr,nphi,nstr,numu,nlout,levout,nzup,nzdn
      integer npd,npd1,ipd(mxpd),modepd(mxpd),igs(mxpd)
      integer k_out(mxlout)
      integer naz,nze,iref
c
c****   counters used in map_back and write_mono
c
      integer nz,ni,nlev,l1,l2,lpd,k,n,l_1,l_2
c
c****   spectral binning parameters
c
      integer ibin,ng0
c
c****   gaussian angles for flux integration
c
      real umu_f(mxumu),gwt_f(mxumu)
c
c****   spectral binning parameters
c
      real alb_b(nsol,5,ngrp),dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp)
c
c****   binned layer radiance transmittances and absorptances
c
      double precision refrad_b(mxumu,kp,5,ngrp),
     -                 drefraddx(mxumu,kp,3,ngrp),
     -                 absrad_b(mxumu,kp,5,ngrp),
     -                 dabsraddx(mxumu,kp,3,ngrp),
     -                 brdf_b(mxumu,mxphi,5,ngrp),
     -                 dbrdfdx(mxumu,mxphi,ngrp)
c
c****   binned layer flux transmittances and absorptances
c
      double precision trnflx_b(kp,5,ngrp),dtrnflxdx(kp,3,ngrp),
     -                 refflx_b(kp,5,ngrp),drefflxdx(kp,3,ngrp),
     -                 absflx_b(kp,5,ngrp),dabsflxdx(kp,3,ngrp)
c
c****    layer radiance transmission values for simplified adding method
c
      double precision trnrad(mxumu,kp),refrad(mxumu,kp),
     -                 absrad(mxumu,kp)
      double precision trn_cmu(mxumu,kp),abs_cmu(mxumu,kp)
c
c****    layer flux transmission values for simplified adding method
c
      double precision trnflx(kp),refflx(kp),absflx(kp),trndir(kp)
c
      real dtauaer(nmode,kp)
c
c****   binned solar source function variables and jacobians
c
      double precision upsflxsrc_b(kp,ngrp),dnsflxsrc_b(kp,ngrp),
     -                 sradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddnsflxdx(kp,4,ngrp),dupsflxdx(kp,4,ngrp),
     -                 dsraddx(mxumu,mxphi,kp,4,ngrp)
c
c****   binned thermal source function variables and jacobians
c
      double precision uptflxsrc_b(kp,ngrp),dntflxsrc_b(kp,ngrp),
     -                 tradsrc_b(mxumu,mxphi,kp,ngrp)
c
      double precision ddntflxdx(kp,4,ngrp),duptflxdx(kp,4,ngrp),
     -                 dtraddx(mxumu,mxphi,kp,4,ngrp)
c
      real umu0nz,phi0nz,umu(mxumu),phi(mxphi)
      real dp_dp(mxlout),rmix(kp,ngas)
c
      real ts,t(kp),p(kp)
      real dtauex(kp,2),co_pi0(kp,2),g(kp)
      real pd_pert(mxpd),dx(mxrad,mxpd)
      real tau_ext(kp,3,mxpd),tau_sca(kp,mxpd),g_sca(kp,mxpd)
      real dtausca0,dtausca
c
c****   vertical structure variables at this wavelength
c
      real alb(nsol),tsurf,tatm(kp),dtau(kp),copi0(kp),g0(kp)
c
      real fbeam,btemp,ttemp,temis,alb0
c
      double precision wn_io,dvi
c
c*****   variables passed to int_rad_lev
c
      real ups(mxlout),dns(mxlout),dirs(mxlout),
     -        upth(mxlout),dnth(mxlout)
      real rad_s(mxumu,mxphi,mxlout),rad_th(mxumu,mxphi,mxlout)
c
c****   spectrally-dependent output flux and radiance
c
      real upsflx(mxrad),dnsflx(mxrad),dirsflx(mxrad),
     -     sol_rad(mxumu,mxphi,mxrad),
     -     uptflx(mxrad),dntflx(mxrad),
     -     th_rad(mxumu,mxphi,mxrad)
c
c****    optical property differences for this wn
c
      real dalb,deltau(kp),delpi0(kp),delg(kp)
c
c****    planck functions used in interp_rad
c
      double precision bb_sur,bb(kp),bb_flx_up(kp),bb_flx_dn(kp),
     -                 bb_rad(mxumu,mxphi,kp)
c
c****   interpolated solar and thermal source functions
c
      double precision up_s_src(kp),dn_s_src(kp),
     -                 rad_s_src(mxumu,mxphi,kp)
      double precision up_t_src(kp),dn_t_src(kp),
     -                 rad_t_src(mxumu,mxphi,kp)
c
c****    baseline optical property differences for this wn
c
      real dalb0,deltau0(kp),delpi00(kp),delg0(kp)
c
c****    baseline planck functions used in interp_rad
c
      double precision bb_sur0,bb0(kp),bb_flx_up0(kp),bb_flx_dn0(kp),
     -                 bb_rad0(mxumu,mxphi,kp)
c
c****    baseline layer transmission values for simplified adding method
c
      double precision trndir0(kp),trnflx0(kp),refflx0(kp),absflx0(kp)
      double precision trnrad0(mxumu,kp),refrad0(mxumu,kp),
     -                 absrad0(mxumu,kp)
c
c****   baseline flux and radiance source functions
c
      double precision up_s_src0(kp),dn_s_src0(kp),
     -                 rad_s_src0(mxumu,mxphi,kp)
      double precision up_t_src0(kp),dn_t_src0(kp),
     -                 rad_t_src0(mxumu,mxphi,kp)
c
c****     perturrbed monochormatic radiance and wn interpolation values
c 
      double precision sol_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,2,nsol),
     -                 dsol_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,nsol)
      double precision th_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,2),
     -                 dth_rad1(mxumu,mxphi,mxlout,mxrad,mxpd)
c
c****   flux source functions at wavenumber wn_io for perturbed state
c
      double precision up_s_src1(mxrad,mxpd,2,nsol),
     -                 dn_s_src1(mxrad,mxpd,2,nsol)
      double precision up_t_src1(mxrad,mxpd,2),dn_t_src1(mxrad,mxpd,2)
c
c****   flux source functions wavenumber derivatives
c
      double precision dup_s_src1(mxrad,mxpd,nsol),
     -                 ddn_s_src1(mxrad,mxpd,nsol)
      double precision dup_t_src1(mxrad,mxpd),ddn_t_src1(mxrad,mxpd)
c
c****    flux transmission and absorption needed for simplified
c        adding method at wavenumber wn_io
c
      double precision trndir1_i(mxrad,mxpd,2,nsol),
     -                 trnflx1_i(mxrad,mxpd,2),
     -                 refflx1_i(mxrad,mxpd,2),
     -                 absflx1_i(mxrad,mxpd,2),
     -                 dtrndir1_i(mxrad,mxpd,nsol),
     -                 dtrnflx1_i(mxrad,mxpd),
     -                 drefflx1_i(mxrad,mxpd),
     -                 dabsflx1_i(mxrad,mxpd)
c
      real dtau_l1,dtau_l2,tautest,dtaumin
c
c****    specify the solar zenith angle and gas index
c
      nz = nz0
      ni = ni0
      nlev = nlyr+1
c
c****       define a small value of delta tau.  If the optical depth
c           perturbation is smaller than this, but greater than taumn,
c           set the value to this value
c
      dtaumin = 0.0005
c
c****       evaluate the radiances and partial derivatives
c           at this wavenumber for each optical property 
c           that can vary
c
      do 5001 n=1,npd
c
c****        determine the number of levels at which the state
c            vector variable chages (T=nlev, tau=nlyr, alb=1)
c
          if(abs(ipd(n)) .eq. 1) then
c
c****         surface pressure
c
            l1 = nlev
            l2 = nlev
c
          else
c
            if(abs(ipd(n)) .eq. 2) then
c
c****            atmospheric/surface temperature
c
              l1 = 1
              l2 = nlev
            else
c
              if(abs(ipd(n)) .eq. 3) then
c
c****             gas optical depth
c
                l1 = 1
                l2 = nlev
c
              else
c
                if(abs(ipd(n)) .eq. 4) then
c
c****               aerosol optical depth
c
                  l1 = 1
                  l2 = nlev
c
                else
c
c****               surface albedo
c
                  l1 = nlev
                  l2 = nlev
c
                endif
              endif
            endif
          endif
c
c****      enter the loop over perturbation level
c
          do 4001 lpd=l1,l2
c
c****           reload baseline arrays for optical properties, 
c               planck functions, and source functions
c
              alb0 = alb(nz)
              dalb = dalb0
              do 1001 k=1,nlev-1
                  dtau(k) = dtauex(k,1)
                  copi0(k) = co_pi0(k,1)
                  g0(k) = g(k)
                  deltau(k) = deltau0(k)
                  delpi0(k) = delpi00(k)
                  delg(k) = delg0(k) 
1001          continue
c
              do 1261 k=1,nlev-1    
                  trndir(k) = trndir0(k)
                  trnflx(k) = trnflx0(k)
                  refflx(k) = refflx0(k)
                  absflx(k) = absflx0(k)
                  do 1201 nze=1,numu
                      trnrad(nze,k) = trnrad0(nze,k)
                      refrad(nze,k) = refrad0(nze,k)
                      absrad(nze,k) = absrad0(nze,k)
1201              continue
 1261          continue
c
              if(lsolar) then
                do 1481 k=1,nlev
                    dn_s_src(k) = dn_s_src0(k)
                    up_s_src(k) = up_s_src0(k)
                    do 1461 naz=1,nphi
                        do 1441 nze=1,numu
                            rad_s_src(nze,naz,k) = rad_s_src0(nze,naz,k)
1441                    continue
1461                continue
1481            continue
              endif
c
              if(lplanck) then
                tsurf = ts
                bb_sur = bb_sur0
                do 1681 k=1,nlev
                    tatm(k) = t(k)
                    bb(k) = bb0(k)
                    bb_flx_up(k) = bb_flx_up0(k)
                    bb_flx_dn(k) = bb_flx_dn0(k)
                    dn_t_src(k) = dn_t_src0(k)
                    up_t_src(k) = up_t_src0(k)
                    do 1661 naz=1,nphi
                        do 1641 nze=1,numu
                            bb_rad(nze,naz,k) = bb_rad0(nze,naz,k)
                            rad_t_src(nze,naz,k) = rad_t_src0(nze,naz,k)
1641                    continue
1661                continue
1681            continue
              endif              
c
c****           define the perturbation increments, dx
c
              if(abs(ipd(n)) .eq. 1) then
c
c                define the surface pressure change
c
                dx(lpd,n) = pd_pert(n)*p(nlev)
c
c*****           calculate radiances for perturbed surface pressure:
c                substitute dtauex(nlev,1) for dtauex(nlyr,1)
c
c****             replace the differential optical depth and
c                 single scattering co-albedo in the lowest layer
c
                dtau(nlyr) = dtauex(nlev,1)
                copi0(nlyr) = co_pi0(nlev,1)
c
c****             set the levels for recalculating source functions
c
                l_1 = nlev-1
                l_2 = nlev
c
              else
c
                if(abs(ipd(n)) .eq. 2) then
c
c****               define atmospheric or surface temperature change
c
                  if(lpd .le. nlev) then
                    dx(lpd,n) = pd_pert(n)*t(lpd)
                  else
                    dx(lpd,n) = pd_pert(n)*ts
                  endif
c
c*****             calculate radiances for perturbed temperatures
c
                  if(lpd .le. nlev) then
                    tatm(lpd) = (1.0 + pd_pert(n))*t(lpd)
                  else
                    tsurf = (1.0 + pd_pert(n))*ts
                  endif
c
c****             set the levels for recalculating source functions
c
                  if(lpd .eq. 1) then
                    l_1 = 1
                  else
                    l_1 = lpd - 1
                  endif
                  l_2 = lpd
c
c****               modify optical depths and single scattering abledos
c                   for temperature dependence
c
                  dtau(l_1) = dtauex(l_1,2)
                  copi0(l_1) = co_pi0(l_1,2)
                  if(l_2 .lt. nlev) then
                    dtau(l_2) = dtauex(l_2,2)
                    copi0(l_2) = co_pi0(l_2,2)
                  endif
c
                else
c
                  if(abs(ipd(n)) .eq. 3) then
c
c****               gas mixing ratio perturbation
c
c****               set the levels for recalculating source functions
c
                    if(lpd .eq. 1) then
                      l_1 = 1
                    else
                      l_1 = lpd - 1
                    endif
c
                    l_2 = lpd
c
c****                 define the gas optical depth change
c
                    dx(lpd,n) = pd_pert(n)*rmix(lpd,igs(n-npd1))
c       
c*****               find optical depths for perturbed gas mixing ratios
c
c****                     note: there should be no perturbation in
c                               g(k) since the asymmetry parameter
c                               for a gas, g_sca, should be 0.0
c
c****                    a change in mixing ratio at level lpd affects
c                        the optical depth in layer l_1 and l_2
c
c****                    modify gas rmix at bottom of layer l_2.  
c                        subtract off the nominal optical depth of 
c                        the gas, tau_ext(l_1,1,n-npd1), and add back
c                        the perturbed value tau_ext(l_1,3,n-npd1)
c
                    dtau_l1 = -tau_ext(l_1,1,n-npd1) +
     -                             tau_ext(l_1,3,n-npd1)
c
c****                    modify gas mixing ratio at top of layer l_2
c                        subtract off the nominal optical depth of 
c                        the gas, tau_ext(l_3,1,n-npd1), and add back
c                        the perturbed value tau_ext(l_3,3,n-npd1)
c
                    if(l_2 .lt. nlev) then
                      dtau_l2 = -tau_ext(l_2,1,n-npd1) +
     -                             tau_ext(l_2,2,n-npd1)
                    else 
                      dtau_l2 = 0.0
                    endif
c
c****                  define a delta-tau test variable that will be 
c                      used to define dx in the case where dtau_l1 and
c                      dtau_l2 are very small
c
                    tautest = 0.5*(dtau_l1 + dtau_l2)
c
                    if(tautest .ge. dtaumin) then
c
c****                   adjust optical depths in layers l_1 and l_2 
c                       by an amount consistent with 
c                       pd_pert(n)*rmix(lpd)
c
                      dtau(l_1) = dtauex(l_1,1) + dtau_l1
                      if(l_2 .lt. nlev) 
     -                  dtau(l_2) = dtauex(l_2,1) + dtau_l2
c
                    else
c
c****                     don't adjust the optical depth
c
                        dx(lpd,n) = 0.0
                        dtau(l_1) = dtauex(l_1,1)    
                        if(l2 .lt. nlev) 
     -                    dtau(l_2) = dtauex(l_2,1)                  
c
                    endif
c
                    if(dtau(l_1) .ne. 0.0) then
c
c****                  find the total scattering optical depth of the
c                      unperturbed layer, l_1
c
                      dtausca0 = (1.0 - co_pi0(l_1,1))*dtauex(l_1,1)
c
c****                  then next statement assumes that perturbing the
c                      mixing ratio does not change the scattering 
c                      optical depth, while it does change the 
c                      single scattering albedo.
c 
                      copi0(l_1) = 1.0 - (dtausca0/dtau(l_1))
c
                    else
c
                      copi0(l_1) = 0.0
c
                    endif
c
                    if(dtau(l_2) .ne. 0.0 .and. l_2 .lt. nlev) then
c
c****                  find the total scattering optical depth of the
c                      unperturbed layer l_2
c
                      dtausca0 = (1.0 - co_pi0(l_2,1))*dtauex(l_2,1)
c
                      copi0(l_2) = 1.0 - (dtausca0/dtau(l_2))
c
                    else
c
                      copi0(l_2) = 0.0
c
                    endif
c      if(wn_io .gt. 5307.45 .and. wn_io .lt. 5307.5 .and. n .eq. 1 
c     - .and. lpd .eq. 12 .and. nz .eq. 1) 
c     - write(*,'(1a,i5,1pe14.6,14(1pe12.4))') 
c     - 'perturb_rad: ',ibin,wn_io,
c     - tau_ext(l_1,1,n-npd1),tau_ext(l_1,3,n-npd1),
c     - dtauex(l_1,1),dtau(l_1),
c     - tau_ext(l_2,1,n-npd1),tau_ext(l_2,2,n-npd1),
c     - dtauex(l_2,1),dtau(l_2),
c     - copi0(l_1),copi0(l_2),tautest,dx(lpd,n)
c
                  else
c
                    if(abs(ipd(n)) .eq. 4) then
c
c****                   define the cloud/aerosol optical depth change
c
                      dx(lpd,n) = pd_pert(n)*dtauaer(modepd(n-npd1),lpd)
c
c****                   set levels for recalculating source functions
c
                      l_1 = lpd
                      l_2 = lpd
c       
c*****                  calculate radiances for perturbed aerosol 
c                       optical depth
c
                      if(lpd .lt. nlev) then
c
                        dtau(lpd) = dtauex(lpd,1) + 
     -                              pd_pert(n)*tau_ext(lpd,1,n-npd1)
c
                        if(dtau(lpd) .ne. 0.0) then
c
c****                       define the nominal scattering optical depth
c
                          dtausca0 = (1.0 - co_pi0(lpd,1))*
     -                               dtauex(lpd,1)
c
c****                       perturb the scattering optical depth
c
                          dtausca = dtausca0 + 
     -                              pd_pert(n)*tau_sca(lpd,n-npd1)
c
c****                      recompute perturbed single scattering albedo
c
                          copi0(lpd) = 1.0 - (dtausca/dtau(lpd))
c
c****                        calculate  perturbed asymmetry parameter
c
                          if(dtausca .ne. 0.0) then
                            g0(lpd) = (g(lpd)*dtausca0 + 
     -                                pd_pert(n)*tau_sca(lpd,n-npd1)*
     -                                g_sca(lpd,n-npd1))/dtausca
                          else
                            g0(lpd) = 0.0
                          endif
c
                        else
                          copi0(lpd) = 0.0
                          g0(lpd) = 0.0
                        endif
c
                      endif
c
                    else
c
c****                  set levels for recalculating source functions
c
                      l_1 = nlev
                      l_2 = nlev
c
c****                   define the surface albedo change
c
                      dx(lpd,n) = pd_pert(n)*alb(nz)
c
c***                   calculate radiances for perturbed albedos
c
                      alb0 = (1.0 + pd_pert(n))*alb(nz)
c
                    endif
c
                  endif
c
                endif       
c
              endif       
c
c****           Layer transmittances are used in mapping fluxes and
c               radiances back to the high resolution spectrum.
c               find the transmission through each layer for each 
c               stream.  
c
              ng0 = ibin
c
              if(lsolar) then
c
c****             find the direct beam transmission for sun
c
                do 3011 k=1,nlyr
                    trndir(k) = exp(-dtau(k)/umu0nz)
3011            continue
c
              endif
c
c      if(wn_io .gt. 6320.5 .and. wn_io .lt. 6320.6 .and. n .eq. 1 
c     - .and. lpd .eq. 10 .and. nz .eq. 1) 
c     - write(*,'(/,1a,i5,1pe14.6,14(1pe12.4))') 
c     - 'perturb_rad: ',ibin,wn_io,
c     - dtau(lpd),copi0(lpd),g0(lpd)
c
c****           calculate radiances for this case
c          
              if(ibin .gt. 0) then
c
                call interp_rad(lsolar,lplanck,
     -           ng0,l_1,l_2,nlev,nz0,nphi,numu,nzup,nzdn,
     -           wn_io,umu0nz,umu,alb_b,dtau_b,copi0_b,g_b,
     -           trnflx_b,dtrnflxdx,refflx_b,drefflxdx,
     -           absflx_b,dabsflxdx,refrad_b,drefraddx,
     -           absrad_b,dabsraddx,brdf_b,dbrdfdx,
     -           tsurf,tatm,alb0,dtau,copi0,g0,dalb,deltau,delpi0,delg,
     -           dnsflxsrc_b,ddnsflxdx,upsflxsrc_b,dupsflxdx,
     -           dntflxsrc_b,ddntflxdx,uptflxsrc_b,duptflxdx,
     -           sradsrc_b,dsraddx,tradsrc_b,dtraddx,
     -           bb,bb_sur,bb_flx_up,bb_flx_dn,bb_rad,
     -           dn_s_src,up_s_src,rad_s_src,
     -           dn_t_src,up_t_src,rad_t_src,
     -           trndir,trnflx,refflx,absflx,trnrad,refrad,absrad,
     -           upsflx,dnsflx,dirsflx,sol_rad,uptflx,dntflx,th_rad)
c
c****            initialize flux transmision values for simplified
c                adding sourse function integration
c
                trnflx1_i(lpd,n,ni) = trnflx(lpd)
                refflx1_i(lpd,n,ni) = refflx(lpd)
                absflx1_i(lpd,n,ni) = absflx(lpd)
                if(lsolar) then
                  if(lpd .eq. 1) then
                    trndir1_i(lpd,n,ni,nz) = trndir(lpd)
                  else
                    trndir1_i(lpd,n,ni,nz) = trndir(lpd)*
     -                                       trndir1_i(lpd-1,n,ni,nz)
                  endif
                endif
c
              else
c
c****             use the pure absorption model to calculate radiances
c
                fbeam = 1.0
                temis = 1.0
                btemp = tsurf
                ttemp = 0.0
c
                call rad_trn(usrang,lplanck,l_1,l_2,ng0,nstr,numu,
     -                       umu,umu_f,gwt_f,dtau,copi0,g0,
     -                       trnrad,absrad,trn_cmu,abs_cmu)
c
                call grey_eq_trn(lsolar,lplanck,usrang,lamber,
     -                       nlev,nz0,nstr,numu,nphi,
     -                       nzup,nzdn,iref,wn_io,dtau,copi0,tatm,
     -                       umu,umu_f,gwt_f,phi,umu0nz,phi0nz,
     -                       fbeam,alb0,btemp,ttemp,temis,
     -                       trndir,trnflx,absflx,trn_cmu,abs_cmu,
     -                       trnrad,absrad,
     -                       upsflx,dnsflx,dirsflx,sol_rad,
     -                       uptflx,dntflx,th_rad,
     -                       dn_s_src,up_s_src,dn_t_src,up_t_src)
c
c****            initialize flux transmision values for simplified
c                adding sourse function integration
c
                trnflx1_i(lpd,n,ni) = trnflx(lpd)
                refflx1_i(lpd,n,ni) = refflx(lpd)
                absflx1_i(lpd,n,ni) = absflx(lpd)
c
                if(lsolar) then
                  if(lpd .eq. 1) then
                    trndir1_i(lpd,n,ni,nz) = trndir(lpd)
                  else
                    trndir1_i(lpd,n,ni,nz) = trndir(lpd)*
     -                                       trndir1_i(lpd-1,n,ni,nz)
                  endif
                endif
c
              endif
c
c****           find the spectral gradients in the flux transmission and
c               absorption values used for the simplified adding method
c
              if(lsolar) dtrndir1_i(lpd,n,nz) = 
     -                            (trndir1_i(lpd,n,2,nz) - 
     -                                   trndir1_i(lpd,n,1,nz))*dvi 
              dtrnflx1_i(lpd,n) = (trnflx1_i(lpd,n,2) - 
     -                                   trnflx1_i(lpd,n,1))*dvi
              drefflx1_i(lpd,n) = (refflx1_i(lpd,n,2) - 
     -                                   refflx1_i(lpd,n,1))*dvi
              dabsflx1_i(lpd,n) = (absflx1_i(lpd,n,2) - 
     -                                   absflx1_i(lpd,n,1))*dvi
c
c****             interpolate flux and radiance values to their 
c                 output levels
c
              call int_rad_lev(lsolar,lplanck,nphi,numu,nz,
     -                       nlyr,nlout,levout,k_out,dp_dp,
     -                       upsflx,dnsflx,dirsflx,uptflx,dntflx,
     -                       sol_rad,th_rad,ups,dns,dirs,
     -                       upth,dnth,rad_s,rad_th)
c
c****             find the spectral gradients in the flux and radiance
c
c
              if(lsolar) then
c
c****             define layer dependent solar source terms at wn_io
c
                dn_s_src1(lpd,n,ni,nz) = dn_s_src(lpd)
                up_s_src1(lpd,n,ni,nz) = up_s_src(lpd)
c      
c****           define solar source term wavelength derivatives
c
                dup_s_src1(lpd,n,nz) = (up_s_src1(lpd,n,2,nz) - 
     -                                  up_s_src1(lpd,n,1,nz))*dvi 
                ddn_s_src1(lpd,n,nz) = (dn_s_src1(lpd,n,2,nz) - 
     -                                  dn_s_src1(lpd,n,1,nz))*dvi
c
c****              define solar radiance wavelength derivatives
c                  at the output levels, k_out
c
                do 3261 k=1,nlout
                    do 3241 naz=1,nphi
                        do 3221 nze=1,numu
c
c****                         define solar radiances and wavelength 
c                             derivatives at output levels, k_out
c
                            sol_rad1(nze,naz,k,lpd,n,ni,nz) = 
     -                                rad_s(nze,naz,k)
                            dsol_rad1(nze,naz,k,lpd,n,nz) = 
     -                           (sol_rad1(nze,naz,k,lpd,n,2,nz) - 
     -                            sol_rad1(nze,naz,k,lpd,n,1,nz))*dvi
3221                    continue
3241                continue
3261            continue
c      if(wn_io .gt. 6320.5 .and. wn_io .lt. 6320.6 .and. n .eq. 1 
c     - .and. lpd .eq. 10 .and. nz .eq. 1) then
c      write(*,'(1a,16(1pe12.4))') 'sol_rad1(1)',
c     - (sol_rad1(numu,1,1,k,n,1,nz),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'sol_rad1(2)',
c     - (sol_rad1(numu,1,1,k,n,2,nz),k=1,nlev)
c      write(*,'(1a,16(1pe12.4))') 'dsol_rad1: ',
c     - (dsol_rad1(numu,1,1,k,n,nz),k=1,nlev)
c      endif
c
              endif
c
              if(lplanck .and. nz .eq. 1) then      
c
c****             define layer dependent solar source terms at wn_io
c
                dn_t_src1(lpd,n,ni) = dn_t_src(lpd)
                up_t_src1(lpd,n,ni) = up_t_src(lpd)
c      
c****           define solar source term wavelength derivatives
c
                dup_t_src1(lpd,n) = (up_t_src1(lpd,n,2) - 
     -                               up_t_src1(lpd,n,1))*dvi 
                ddn_t_src1(lpd,n) = (dn_t_src1(lpd,n,2) - 
     -                               dn_t_src1(lpd,n,1))*dvi
c
c****              define thermal radiance wavelength derivatives
c                  at the output levels, k_out
c
                do 3461 k=1,nlout
                    do 3441 naz=1,nphi
                        do 3421 nze=1,numu
c
c****                         define solar radiances and wavelength 
c                             derivatives at output levels, k_out
c
                            th_rad1(nze,naz,k,lpd,n,ni) = 
     -                                rad_th(nze,naz,k)
                            dth_rad1(nze,naz,k,lpd,n) = 
     -                               (th_rad1(nze,naz,k,lpd,n,2) - 
     -                                th_rad1(nze,naz,k,lpd,n,1))*dvi
3421                    continue
3441                continue
3461            continue
c 
              endif
c
4001      continue
c
c****      exit loop over output levels
c
5001  continue
c
c****   exit loop over constituent
c
      return
      end
