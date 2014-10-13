      subroutine map_back(lsolar,lplanck,lamber,usrang,iugrp,iuthrm,
     -          iusol0,iusol2,iutrn,iu_flux,iuout,npd,ipd,iunits,irad,
     -          iu_pd,iutpd,iuspd,iupdrad,ifrmout,nza,nz0,nstr,numu,
     -          nzup,nzdn,nphi,nlyr,nlay,levout,nlout,k_out,igs,iref,
     -          modepd,nt_pd,nstate,istate,ntau_pd,isptype,islit,
     -          width,dwn,points,wnmin,wnmax,wnsmin,wnsmax,wn_tol,
     -          units,umu0nz,phi0nz,umu_f,gwt_f,umu,phi,pd_pert,
     -          p,t,ts,dp_dp,rmix,dtauaer,alb_b,dtau_b,copi0_b,g_b,
     -          trnflx_b,dtrnflxdx,refflx_b,drefflxdx,
     -          absflx_b,dabsflxdx,refrad_b,drefraddx,
     -          absrad_b,dabsraddx,brdf_b,dbrdfdx,
     -          dnsflxsrc_b,ddnsflxdx,upsflxsrc_b,dupsflxdx,
     -          dntflxsrc_b,ddntflxdx,uptflxsrc_b,duptflxdx,
     -          sradsrc_b,dsraddx,tradsrc_b,dtraddx,
     -          dirsoflx,dnsoflx,upsoflx,dnthflx,upthflx)
c
cccccccccccccccccccccccccc   m a p _ b a c k   ccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine maps binned radiances back to a high-           cc
cc    resolution spectral grid.                                       cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      iugrp - unit number of file with binned radiances             cc
cc      iutrn - unit number for output transmission file              cc
cc    iu_flux - unit number for output level-dependent fluxes         cc
cc     iusol0 - unit number for solar radiance scratch file           cc
cc     iusol2 - unit number with solar radiances                      cc
cc     iuthrm - unit number for thermal radiances                     cc
cc      iuout - unit number of output radiance file                   cc
cc      iutpd - unit numbers for output level-dependent               cc
cc              thermal fluxes and thier partial derivatives          cc
cc      iuspd - unit numbers for output level-dependent               cc
cc              solar fluxes and thier partial derivatives            cc
cc    iupdrad - unit numbers for output level-dependent radiances     cc
cc              and thier partial derivatives                         cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) ergs/s/cm**2/sr/cm-1                               cc
cc              5) photons/s/m**2/sr/micron                           cc
cc        nz0 - index of this solar zenith angle (1 to nsol)          cc
cc        nza - number of solar zenith angles                         cc
cc    ifrmout - output file format (1) ascii or (2) binary            cc
cc     levout - output level index (1) top of atmosphere,             cc
cc              (2) surface, (3) arbitrary level                      cc
cc      dp_dp - vertical pressure gradient (p-pout)/dp                cc
cc      k_out - index of output level (1 - nlev)                      cc
cc       nlyr - number of computational model layers                  cc
cc       nphi - number of output azimuth angles                       cc
cc       numu - number of zenith angles in input file                 cc
cc    isptype - output spectrum type:                                 cc
cc              1) Full-resolution spectrum                           cc
cc              2) Sampled spectrum smoothed with a slit function     cc
cc       irad - index of output file type:                            cc
cc              1) fluxes, radiances, and heating rates,              cc
cc                 at computational azimuths and zenith angles,       cc
cc              2) fluxes, radiances, heating rates, and transmission cc
cc                functions at computational zenith angles,           cc
cc              3) fluxes, radiances, heating rates, and contribution cc
cc                 functions at computational zenith angles,          cc
cc              4) fluxes, radiances, heating rates, transmission     cc
cc                 functions and and contribution functions           cc
cc                 at computational zenith angles,                    cc
cc              5) fluxes, radiances, and heating rates,              cc
cc                 at computational azimuths and zenith angles,       cc
cc              6) fluxes, radiances and transmission functions       cc
cc                 at arbitrary zenith angles,                        cc
cc              7) fluxes, radiances, and contribution functions      cc
cc                 at arbitrary zenith angles,                        cc
cc              8) fluxes, radiances, transmission functions,and      cc
cc                 contribution functions at arbitrary zenith angles. cc
cc     nstate - number of variable elements in the state vecror       cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc              1 - surface pressure                                  cc
cc              2 - surface/atmospheric temperature                   cc
cc              3 - gas absorption coeffient                          cc
cc              4 - cloud/aerosol optical depth                       cc
cc              5 - surface albedo                                    cc
cc      islit - index of slit function type:                          cc
cc              1) boxcar                                             cc
cc              2) triangular (approximate slit spectrometer)         cc
cc      width - half-width at half max of slit function  (cm**-1)     cc
cc        dwn - output sampling resolution (cm**-1)                   cc
cc         wn - current wavenumber (cm**-1)                           cc
cc     solflx - solar flux at this wavenumber                         cc
cc      wnmin - minimum wavenumber in output spectral grid (cm**-1)   cc
cc      wnmax - maximum wavenumber in output spectral grid (cm**-1)   cc
cc     wnsmin - minimum wavenumber in this binned segment (cm**-1)    cc
cc     wnsmax - maximum wavenumber in this binned segment (cm**-1)    cc
cc         ts - surface temperature (K)                               cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc         wl - current wavenumber (cm**-1)                           cc
cc         wn - current wavenumber (cm**-1)                           cc
cc    dir_flx - downward direct flux (irradiance) at wavenumber wn    cc
cc     dn_flx - total downward flux (irradiance) at wavenumber wn     cc
cc     up_flx - upward flux (irradiance) at wavenumber wn             cc
cc        rad - radiance at each output zenith and azimuth angle      cc
cc  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc
cc  trn_gas_0 - normal-incidence gas transmission                     cc
cc  trn_ray_0 - normal-incidence aerosol transmission                 cc
cc                                                                    cc
cccccccccccccccccccccccccc   m a p _ b a c k   ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nfmx
      parameter (nfmx = (mxlout+mxpd+1)*nsol+1)
c
      logical lplanck,lsolar,lamber,usrang
c
c****   spectrally-dependent output flux and radiance
c
      integer iugrp,iutrn,iu_flux,iusol0,iusol2,iuthrm,iuout,ifrmout,
     -        iunits,levout,nlout,irad,nlyr,nlay,nza,nz0,nphi,nstr,numu,
     -        isptype,islit,nt_pd,iref
      integer k_out(mxlout)
      integer k,nout,n1,ne,ne0,ninter,n
      integer nzup,nzdn,nz,iuin
c
c****    number and type of partial derivatives
c
      integer igs(mxpd),ntau_pd,modepd(mxpd)
      integer npd1
      save npd1
c
c****   counters used in map_back
c
      integer iflout(nex,nsol),nsiout(2,nex,nsol)
      save iflout,nsiout
      integer nsivar(2),iend(nex)
      save nsivar,iend
      integer ifl, iend_ne
c
c****    units for flux and radiance jacobians
c 
      integer nstate,istate(nex),npd,ipd(mxpd),
     -        iu_pd,iutpd(mxpd,nsol),iuspd(mxpd,nsol),
     -        iupdrad(mxumu,mxphi,mxpd,mxlout,nsol)
c
c****   spectral binning parameters
c
      integer ibin
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
c****    aerosol optical depth for each particle mode
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
      real umu(mxumu),phi(mxphi),umu0nz,phi0nz
      real dp_dp(mxlout)
      real p(kp),t(kp),ts,rmix(kp,ngas)
c
      real umu_f(mxumu),gwt_f(mxumu)
c
      double precision wn,wnmin,wnmax,width,dwn,
     -                 wnsmin,wnsmax,wn_tol,delnu
c
c****     monochormatic radiance interpolation variables
c 
      double precision wnout(2,nex,nsol),wn_io,wn0,wn1,dvi
      save wnout
      double precision dwn_max,dist(nex),distmin,dist0,wntst,wn_next
c
c****   solar variables 
c
      real solflx,alb(nsol)
      real sol0(2),dsoldv,units
      save sol0,dsoldv
c
      double precision wnvar(2),wns1
      save wnvar,wns1
c
c****    output fluxes and radiances
c
      real dir_s_flx(mxrad,nsol,2),up_s_flx(mxrad,nsol,2),
     -     dn_s_flx(mxrad,nsol,2),up_t_flx(mxrad,2),dn_t_flx(mxrad,2)
c
      real up_flx(mxlout),dn_flx(mxlout),
     -     dir_flx(mxlout),rad(mxumu,mxphi,mxlout)
c
      real p_ray0(mxumu),p_gas0(mxumu),p_aer0(mxumu),
     -        tau_ray_0,tau_gas_0,tau_aer_0,
     -        p_ray_0,p_gas_0,p_aer_0
c
      save p_ray0,p_gas0,p_aer0,
     -     tau_ray_0,tau_gas_0,tau_aer_0,
     -     p_ray_0,p_gas_0,p_aer_0
c
c****    flux transmission and absorption functions at wavenuber, wn
c
      real trn_dir(mxrad),trn_flx(mxrad),ref_flx(mxrad),abs_flx(mxrad)
c
c****    flux transmission and absorption partial 
c        derivatives for simplified adding method at wavenumber wn
c
      real pd_trndir(mxrad,mxpd),
     -     pd_trnflx(mxrad,mxpd),
     -     pd_refflx(mxrad,mxpd),
     -     pd_absflx(mxrad,mxpd)
c
c****   perturbations used for partial derivatives
c
      real pd_pert(mxpd),tau_ext(kp,3,mxpd)
c
c****   layer-dependent source terms for solar and thermal fluxes
c
      real dns_src(mxrad),ups_src(mxrad),dnt_src(mxrad),upt_src(mxrad)
c
c****   output flux and radiance partial derivatives
c
      real pd_ups_src(mxrad,mxpd),pd_dns_src(mxrad,mxpd),
     -     pd_upt_src(mxrad,mxpd),pd_dnt_src(mxrad,mxpd),
     -     pd_rad(mxumu,mxphi,mxlout,mxrad,mxpd)
c
c****   slit function variables
c
      real points(nfmx)
c
c****   double precision variables for flux integration
c
      double precision dirsoflx(mxulv,nsol),dnsoflx(mxulv,nsol),
     -       upsoflx(mxulv,nsol),dnthflx(mxulv),upthflx(mxulv)
c
c****   specify the solar zenith angle index
c
      nz = nz0
c
c**** rewind thermal flux/radiance file if necessary
c
      if(lplanck) rewind(iuthrm)
c
c****    set a maximum spectral interval step size
c
      dwn_max = 0.01d0*(wnmax - wnmin)
      if(dwn_max .lt. 0.01d0*wnmin) dwn_max = 0.01d0*wnmin

c*****   initialize the spectral grid flags for the binned
c        radiances and the solar fluxes
c
      if(lsolar) then
c
c         Two spectral grids (solar fluxes and spectral radiance) 
c         must be combined, so the number of outoput grirds, nout = 2.
c          
        nout = 2
        if(wnsmin .eq. wnmin) then
c
c****      set the value of the variable n1 = 1 so that the
c          index counters, nsiout, for both the solar fluxes 
c          and the binned radiances are reset
c
          n1 = 1
c
c****      since this is the first spectral interval, set the
c          "last solar flux" quantity, wns1 = wnsmin
c
          wns1 = wnsmin
        else
c
c****      Set the the variable n1 = 2 such that the index
c          counters for the binned radiances are reset, but
c          the index counters for the solar flux are not reset
c 
          n1 = 2
c
c****       test to see if the last solar flux value is beyond
c           the current wavenumber, wn.  If it is not, set the
c           mininum wavenumber value for this interval to a 
c           value that is just beyond the last solar flux value.
c
          if(nz .eq. 1) then
            wn_next = (1.0d0 + wn_tol)*wns1
            if(wn_next .lt. wnsmin) wnsmin =wn_next
          endif
        endif
c
c****     wavnumber dependent solar fluxes are written to 
c         unit iusol0 as the program processes the first solar
c         zenith angle.  These values are used for all
c         other solar zenith angles. If this is the first 
c         solar zenith angle, close and reopen the solar 
c         flux unit. if multiple solar zenith angles are 
c         used rewind this unit
c
        if(nza .gt. 1) then
          if(nz .eq. 1) then
c
            close(iusol0)
            open(iusol0,form='unformatted',status='scratch')
c
          else
            rewind(iusol0)
          endif
        endif
c
      else 
c
c****     solar fluxes are not included.  Thermal fluxes are the only
c         output quantity, so nout = 1, and n1 = 1
c
        nout = 1
        n1 = 1
c
      endif
c
c****  initialize wavelength pointers
c
      wn = wnsmin
      wn0 = wnsmin
      wn1 = wnsmin
      wn_io = -999.0d0
c
      if(wnsmin .eq. wnmin) then
c
c****     intialize output flags and spectral grid pointers.
c
        do 1001 ne = n1,nout
c
c****      set output flag
c
            iflout(ne,nz) = 1
c
c****         set end of bin file flag
c
            iend(ne) = 0
c
c****      initialize spectral counters for radiances and solar flux.
c
            nsiout(1,ne,nz) = 0
            nsiout(2,ne,nz) = 0
1001    continue
c
c****     set partial derivative pointers
c
        do 1221 n=1,nstate
            if(istate(n) .ne. 0) then
c
              do 1201 k=1,mxrad
                  pd_ups_src(k,n) = 0.0
                  pd_dns_src(k,n) = 0.0
                  pd_upt_src(k,n) = 0.0
                  pd_dnt_src(k,n) = 0.0
1201          continue        
c
            endif
1221    continue
c
c****     set the index variable for optical depths: if ipd = 1, 
c         the first varible component of the state vector is pressure.
c         if ipd = 2, the second variable component of the state vector 
c         is temperature.  
c
        npd1 = 0
        if(abs(ipd(1)) .gt. 0 .and. abs(ipd(1)) .lt. 3) npd1 = 1
        if(abs(ipd(2)) .gt. 0 .and. abs(ipd(2)) .lt. 3) npd1 = npd1 + 1
c       
      endif
c     
      ninter = 0
c
c****        s p e c t r a l    i n t e r v a l    l o o p
c
c****   initialize spectral interval counters
c
2001  ne = 0
c
          ninter = ninter + 1
c
c****       f i n d    t o a   s o l a r   f l u x e s
c
          if(lsolar) then
c
c****          increment spectral counter and check wn grid flag
c
            ne = ne + 1
c
c****          if this is the first solar zenith angle, read the 
c              toa solar flux values.
c
            if(nz .eq. 1) then
c
              if(iflout(ne,nz) .eq. 1) then
                nsivar(1) = nsiout(1,ne,nz)
                nsivar(2) = nsiout(2,ne,nz)
                wnvar(1) = wnout(1,ne,nz)
                wnvar(2) = wnout(2,ne,nz)
c
c****            read solar fluxes at next wavenumber
c
                call solar(iusol2,ne,nsivar,iend,wnmax,wn,wnvar,
     -                     sol0,dsoldv) 
c
                 nsiout(1,ne,nz) = nsivar(1)
                 nsiout(2,ne,nz) = nsivar(2)
                 wnout(1,ne,nz) = wnvar(1)
                 wnout(2,ne,nz) = wnvar(2)
c
c****             turn off the wavenumber grid flag.
c    
                iflout(ne,nz) = 0
c
              endif
c
c****           interpolate toa solar fluxes to this wn
c
              solflx = sol0(1) + dsoldv*real(wn - wnout(1,ne,nz))
              if(solflx .lt. 0.0) solflx = 0.0
c
c****           if there is more than one solar zenith angle,save
c              toa solar fluxes and spectral derivative test values
c
              if(nza .gt. 1)
     -          write(iusol0) 
     -           nsiout(1,ne,nz),nsiout(2,ne,nz),
     -           wnout(nsiout(1,ne,nz),ne,nz),
     -           wnout(nsiout(2,ne,nz),ne,nz),
     -           wn,solflx
c
              else
c
c****           read previously mapped solar fluxes
c
                read(iusol0,end=2201) 
     -            nsiout(1,ne,nz),nsiout(2,ne,nz),
     -            wnout(nsiout(1,ne,nz),ne,nz),
     -            wnout(nsiout(2,ne,nz),ne,nz),
     -            wntst,solflx
                go to 2221
c
c****             print end of file error if necessary and stop
c
2201            write(*,'(1a,i5,8(1pe14.7))') 
     -             'map_back: read past EOF for solar flux unit',
     -              iusol0,wn,wntst,wn0,wnsmax,wns1
c  
                stop
c
c****            check wavelength grid alignment
c
2221            if(abs(wntst-wn) .gt. 2.0d0*wn_tol*wn) then
                  write(*,'(/,1a,/,2a,2(1pe14.6))') 
     -            'Error in map_back:',
     -            'wavelength grids not aligned for different ',
     -            'solar zenith angles: wn, wntst: ',wn,wntst
                  write(*,'(/,7a)') 'nz,ninter',
     -            'nsiout(1,1,nz),nsiout(2,1,nz)',
     -            'nsiout(1,2,nz),nsiout(2,2,nz)',
     -            'wnout(nsiout(1,1,nz),1,nz)',
     -            'wnout(nsiout(2,1,nz),1,nz)',
     -            'wnout(nsiout(1,2,nz),2,nz)',
     -            'wnout(nsiout(2,2,nz),2,nz)',
     -            'wntst,solflx'
                  write(*,'(/,6i5,6(1pe14.6))') nz,ninter,
     -              nsiout(1,1,nz),nsiout(2,1,nz),
     -              nsiout(1,2,nz),nsiout(2,2,nz),
     -              wnout(nsiout(1,1,nz),1,nz),
     -              wnout(nsiout(2,1,nz),1,nz),
     -              wnout(nsiout(1,2,nz),2,nz),
     -              wnout(nsiout(2,2,nz),2,nz),
     -              wntst,solflx
                  stop
c
                endif
              endif
            endif 
c
c****        m a p    s p e c t r a l    b i n s   t o   t h i s   wn
c
            ne = ne + 1
            ifl = iflout(ne,nz)
            iend_ne = iend(ne)
c
            if(wn_io .lt. wnsmax) then
              ne0 = ne
c
****            interpolate binned properties to wavenumber wn
c
              call map_rad(lsolar,lplanck,lamber,usrang,iugrp,iuthrm,
     -         ifl,iend_ne,ne0,nz0,nphi,nstr,numu,nzup,nzdn,nza,
     -         nlyr,nlay,nlout,levout,k_out,iref,ibin,nsiout,
     -         igs,npd,npd1,ipd,ntau_pd,nt_pd,modepd,irad,
     -         wnmin,wnmax,wn,wn_io,wnout,wn_tol,dvi,umu0nz,phi0nz,
     -         solflx,umu_f,gwt_f,umu,phi,ts,t,p,dp_dp,rmix,alb,
     -         dtauaer,tau_ext,pd_pert,alb_b,dtau_b,copi0_b,g_b,
     -         p_ray_0,p_gas_0,p_aer_0,p_ray0,p_gas0,p_aer0,tau_ray_0,
     -         tau_gas_0,tau_aer_0,trnflx_b,dtrnflxdx,refflx_b,
     -         drefflxdx,absflx_b,dabsflxdx,refrad_b,drefraddx,
     -         absrad_b,dabsraddx,brdf_b,dbrdfdx,
     -         dnsflxsrc_b,ddnsflxdx,upsflxsrc_b,dupsflxdx,
     -         dntflxsrc_b,ddntflxdx,uptflxsrc_b,duptflxdx,
     -         tradsrc_b,dtraddx,sradsrc_b,dsraddx,
     -         up_s_flx,dn_s_flx,dir_s_flx,up_t_flx,dn_t_flx,
     -         dns_src,ups_src,pd_dns_src,pd_ups_src,
     -         dnt_src,upt_src,pd_dnt_src,pd_upt_src,
     -         pd_rad,pd_trndir,pd_trnflx,pd_refflx,pd_absflx,up_flx,
     -         dn_flx,dir_flx,rad,trn_dir,trn_flx,ref_flx,abs_flx)
c
c****           turn off the wavenumber grid flag.
c
              if(iend(ne) .eq. 0) then
                iflout(ne,nz) = 0
              else
c
c*****           end if iugrp file encountered
c
                write(*,'(1a,8(1pe14.7))')
     -              'End of bin file: wn,wn_io',wn,wn_io,wnsmax
                wnsmax = wn_io
                if(nz .eq. nza) wns1 = wn
                return
c
              endif
c
            endif
c
c****           e f f e c t i v e    i n t e r v a l    w i d t h
c
c*****       find the next wavenumber where monochromatic properties
c            are specified.  Check all input data sets and select the
c            wavenumber that is closest to the current wavenumber.
c
c            note: the minimum distance must exceed the round-off
c                tolerance of the computer, wn_tol
c
            distmin = 1.0d30
            dist0 = wn_tol*wn
c
            do 6001 ne0=1,nout
                dist(ne0) =  wnout(nsiout(2,ne0,nz),ne0,nz) - wn
                if(dist(ne0) .lt. distmin) then
                  if(dist(ne0) .gt. dist0) then
                    distmin = dist(ne0)
                  else
                    distmin = dist0
                  endif
                endif
6001        continue
c
            if(distmin .gt. dwn_max) distmin = dwn_max
c
c****        set the next wavenumber and the effective spectral 
c            interval width, dnu.  The spectral interval is assumed
c            to extend between the current wavenumber, wn, and half-way
c            to the previous value, wn0, and the next value, wn1.
c
            wn1 = wn + distmin
            delnu = wn - wn0
c
c****         add this contribution to spectrally-integrated fluxes
c
            call flux_int(lsolar,lplanck,nlyr,nz,ninter,delnu,
     -                    dir_s_flx,dn_s_flx,up_s_flx,
     -                    dn_t_flx,up_t_flx,
     -                    dirsoflx,dnsoflx,upsoflx,
     -                    dnthflx,upthflx)
c
c****          find the factor needed to convert fluxes and radiances
c              from w/m**2/cm-1 to the desired output units
c
            iuin = 1
c
            call find_units(iuin,iunits,wn,units)
c
c****           print these radiances
c          
            call rad_out(lsolar,iuout,iutrn,iu_flux,
     -                   ifrmout,irad,levout,nlout,nlyr,
     -                   nz0,nza,nphi,numu,nzup,isptype,islit,
     -                   wn,wnmin,wnmax,width,dwn,units,p,solflx,
     -                   up_s_flx,dn_s_flx,dir_s_flx,
     -                   up_t_flx,dn_t_flx,
     -                   up_flx,dn_flx,dir_flx,rad,
     -                   p_ray0,p_gas0,p_aer0,
     -                   tau_ray_0,tau_gas_0,tau_aer_0,
     -                   p_ray_0,p_gas_0,p_aer_0,points)
c
            if(nstate .gt. 0) then
c
              call pd_out(nz0,npd,npd1,ipd,igs,ifrmout,
     -                    levout,nlout,nlyr,nphi,numu,nzup,
     -                    iu_pd,iutpd,iuspd,modepd,iupdrad,
     -                    iunits,isptype,islit,lsolar,lplanck,
     -                    wn,wnmin,wnmax,width,dwn,points,units,
     -                    umu0nz,ts,t,p,rmix,alb,dtauaer,
     -                    solflx,up_s_flx,dn_s_flx,dir_s_flx,
     -                    up_t_flx,dn_t_flx,pd_pert,
     -                    dns_src,pd_dns_src,ups_src,pd_ups_src,
     -                    dnt_src,pd_dnt_src,upt_src,pd_upt_src,
     -                    trn_dir,trn_flx,ref_flx,abs_flx,
     -                    pd_trndir,pd_trnflx,pd_refflx,pd_absflx,
     -                    rad,pd_rad)
c
            endif
c
c*****      specify the next wavenumber grid point and reset
c           wavenumber grid flags.
c
            wn0 = wn
            wn = wn1
c   
            do 7001 ne=1,nout
                dist(ne) = wnout(nsiout(2,ne,nz),ne,nz) - wn
                if(dist(ne) .le. 0.0 .and. iend(ne) .eq. 0) then
c
c****              wn is beyond latest input point for this constituent.
c                  set flag to read next point.
c
                  iflout(ne,nz) = 1
                else
                  iflout(ne,nz) = 0
                endif
7001        continue
c
c****         get monochromatic optical properties for next wavenumber.
c
      if(wn0 .lt. wnsmax) go to 2001
c
c****   if this is the last zenith angle, set the end-of interval wn.
c
      if(nz .eq. nza) wns1 = wn1
c
      return
      end
