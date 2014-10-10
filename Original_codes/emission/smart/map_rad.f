      subroutine map_rad(lsolar,lplanck,lamber,usrang,iugrp,iuthrm,
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
cccccccccccccccccccccccccc   m a p _ r a d  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads the wavenumber, interval width, bin       cc
cc    bin index, and albedo of the next spectral point.               cc
cc                                                                    cc
cc    NOTE: this version has been modified to produce layer dependent cc
cc    flux source functions and their partial derivatives, as well    cc
cc    as the corresponding flux values.                               cc
cc                                                                    cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      iugrp - unit number of file with binned radiances             cc
cc         ne - index of extinction source                            cc
cc        nz0 - index of this solar zenith angle (1 to nsol)          cc
cc    iend_ne - end of file flag (0 - okay, 1 - end of file)          cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc       nlyr - number of computational model layers                  cc
cc       nphi - number of output azimuth angles                       cc
cc       numu - number of zenith angles in input file                 cc
cc      wnmin - minimum wavenumber of desired spectral window         cc
cc      wnmax - maximum wavenumber of desired spectral window         cc
cc         wn - current wavenumber                                    cc
cc     uptflx - wn-dependent upward thermal flux                      cc
cc     dntflx - wn-dependent downward thermal flux                    cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc     upsflx - wn-dependent upward solar flux                        cc
cc     dnsflx - wn-dependent downward diffuse + direct solar flux     cc
cc    dirsflx - wn-dependent downward direct solar flux               cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc    tau_ray - wn-dependent rayleigh-scattering optical depth        cc
cc    tau_gas - wn-dependent gas optical depth                        cc
cc    tau_aer - wn-dependent aerosol extinction optical depth         cc
cc     pray_0 - pressure of normial-incident rayleigh optical depth   cc
cc              unity                                                 cc
cc     pgas_0 - pressure of normial-incident gas optical depth unity  cc
cc     paer_0 - pressure of normial-incident aerosol optical depth    cc
cc              unity                                                 cc
cc       pray - pressure of rayleigh optical depth unity along each   cc
cc              output stream                                         cc
cc       pgas - pressure of gas optical depth unity along each stream cc
cc       paer - pressure of aerosol optical depth unity along each    cc
cc              output stream                                         cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    uptflx_i - wn-dependent upward thermal flux                      cc
cc    dntflx_i - wn-dependent downward thermal flux                    cc
cc    th_rad0 - wn-depndent, angle-dependent thermal radiances        cc
cc    upsflx_i - wn-dependent upward solar flux                        cc
cc    dnsflx_i - wn-dependent downward diffuse+direct solar flux       cc
cc   dirsflx_i - wn-dependent downward direct solar flux               cc
cc   sol_rad0 - wn-depndent, angle-dependent solar radiances          cc
cc   tau_ray0 - wn-dependent rayleigh-scattering optical depth        cc
cc   tau_gas0 - wn-dependent gas optical depth                        cc
cc   tau_aer0 - wn-dependent aerosol extinction optical depth         cc
cc    pray_00 - pressure of normial-incident rayleigh optical depth   cc
cc              unity                                                 cc
cc    pgas_00 - pressure of normial-incident gas optical depth unity  cc
cc    paer_00 - pressure of normial-incident aerosol optical depth    cc
cc              unity                                                 cc
cc      pray0 - pressure of rayleigh optical depth unity along each   cc
cc              output stream                                         cc
cc      pgas0 - pressure of gas optical depth unity along each stream cc
cc      paer0 - pressure of aerosol optical depth unity along each    cc
cc              output stream                                         cc
cc    duptflx - wn-dependent upward thermal flux                      cc
cc    ddntflx - wn-dependent downward thermal flux                    cc
cc    dth_rad - wn-depndent, angle-dependent thermal radiances        cc
cc    dupsflx - wn-dependent upward solar flux                        cc
cc    ddnsflx - wn-dependent downward diffuse + direct solar flux     cc
cc   ddirsflx - wn-dependent downward direct solar flux               cc
cc   dsol_rad - wn-depndent, angle-dependent solar radiances          cc
cc                                                                    cc
cccccccccccccccccccccccccc   m a p _ r a d  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar,lamber,usrang
c
      integer iugrp,iuthrm,ifl,iend_ne,ne0,nz0,nza,nphi,nstr,numu,
     -        nzup,nzdn,nlyr,nlout,levout,iref
c     
      integer l_1,l_2,nlev,ni,nz,ne,k,l,n,nze,naz,ni0,nt_pd,ii,m,nn
      integer k_out(mxlout)
      integer npd,npd1,ipd(mxpd),igs(mxpd),ntau_pd,modepd(mxpd)
      integer nlay,irad
c
c****   counters used in map_back
c
      integer nsiout(2,nex,nsol)
c
c****   internal spectral flag for starting wavenumber loop
c
      integer istart
      save istart
c
c****   spectral binning parameters
c
      integer ibin,ng0
c
c****   spectral binning parameters
c
      real alb_b(nsol,5,ngrp),dtau_b(kp,5,ngrp),copi0_b(kp,5,ngrp),
     -     g_b(kp,5,ngrp)
c
      double precision wnmin,wnmax,wn_tol
c
      double precision wn,wn_io,wnout(2,nex,nsol),dv,dvi
c
      real ts,t(kp),p(kp),rmix(kp,ngas),dp_dp(mxlout),solflx
c
c****   gaussian angles for flux integration
c
      real umu_f(mxumu),gwt_f(mxumu)
c
      real umu0nz,phi0nz,umu(mxumu),phi(mxphi)
c
c****   perturbations for finite derivative jacobians
c
      real pd_pert(mxpd),dx(mxrad,mxpd)
c
c****    aerosol optical depth for each particle mode
c
      real dtauaer(nmode,kp)
c
c****    stored monochromatic optical properties
c
      real surf_opt(4),dtauex(kp,2),co_pi0(kp,2),g(kp),
     -     tau_ext(kp,3,mxpd),tau_sca(kp,mxpd),g_sca(kp,mxpd)
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
c****  pressures and optical depths of tau = 1 level
c
      real pray_0,pgas_0,paer_0
      real tau_ray,tau_gas,tau_aer
      real pray(mxumu),pgas(mxumu),paer(mxumu)
c
c****   temporary variables used in interp_rad
c
      real alb(nsol),tsurf,tatm(kp),dtau(kp),copi0(kp),g0(kp)
      real fbeam,btemp,ttemp,temis,alb0
c
c****   spectrally-dependent output flux and radiance from interp_rad
c
      real upsflx(mxrad),dnsflx(mxrad),dirsflx(mxrad),
     -     sol_rad(mxumu,mxphi,mxrad),
     -     uptflx(mxrad),dntflx(mxrad),
     -     th_rad(mxumu,mxphi,mxrad)
c
c****     monochormatic fluxes at wavenumber wn_io
c 
      real upsflx_i(mxrad,2,nsol),dnsflx_i(mxrad,2,nsol),
     -     dirsflx_i(mxrad,2,nsol),uptflx_i(mxrad,2),dntflx_i(mxrad,2)
c
c****   downward and upward solar fluxes gradients at wavenumber wn_io
c
      real dupsflx(mxrad,nsol),ddnsflx(mxrad,nsol),ddirsflx(mxrad,nsol),
     -     duptflx(mxrad),ddntflx(mxrad)
c
      save upsflx_i,dnsflx_i,dirsflx_i,uptflx_i,dntflx_i,
     -     dupsflx,ddnsflx,ddirsflx,duptflx,ddntflx
c
c****   monochormatic radiances and wavelength gradients at wn_io 
c       at output levels, k_out
c
      real sol_rad0(mxumu,mxphi,mxlout,2,nsol),
     -     th_rad0(mxumu,mxphi,mxlout,2)
      real dsol_rad(mxumu,mxphi,mxlout,nsol),dth_rad(mxumu,mxphi,mxrad)
c
      save sol_rad0,th_rad0,dsol_rad,dth_rad
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
      save rad_s_src0,rad_t_src0
c
c****   variables for interpolating source functions to wn
c
      double precision up_s_src_i(mxrad,2,nsol),dn_s_src_i(mxrad,2,nsol)
      double precision up_t_src_i(mxrad,2),dn_t_src_i(mxrad,2)
c
      save up_s_src_i,dn_s_src_i,up_t_src_i,dn_t_src_i
c
c****   wavelength gradients in source terms
c
      double precision dup_s_src(mxrad,nsol),ddn_s_src(mxrad,nsol)
      double precision dup_t_src(mxrad),ddn_t_src(mxrad)
c
      save dup_s_src,ddn_s_src,dup_t_src,ddn_t_src
c
c****   flux source functions at wavenumber wn_io for perturbed state
c
      double precision up_s_src1(mxrad,mxpd,2,nsol),
     -                 dn_s_src1(mxrad,mxpd,2,nsol)
      double precision up_t_src1(mxrad,mxpd,2),dn_t_src1(mxrad,mxpd,2)
c
      save up_s_src1,dn_s_src1,up_t_src1,dn_t_src1
c
c****   flux source functions wavenumber derivatives
c
      double precision dup_s_src1(mxrad,mxpd,nsol),
     -                 ddn_s_src1(mxrad,mxpd,nsol)
      double precision dup_t_src1(mxrad,mxpd),ddn_t_src1(mxrad,mxpd)
c
      save dup_s_src1,ddn_s_src1,dup_t_src1,ddn_t_src1
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
      save trndir1_i,trnflx1_i,refflx1_i,absflx1_i,
     -     dtrndir1_i,dtrnflx1_i,drefflx1_i,dabsflx1_i
c
c****    flux transmission and absorption needed for simplified
c        adding method at wavenumber wn_io
c
      double precision trndir_i(mxrad,2,nsol),trnflx_i(mxrad,2),
     -                 refflx_i(mxrad,2),absflx_i(mxrad,2),
     -                 dtrndir_i(mxrad,nsol),dtrnflx_i(mxrad),
     -                 drefflx_i(mxrad),dabsflx_i(mxrad)
c
      save trndir_i,trnflx_i,refflx_i,absflx_i,
     -     dtrndir_i,dtrnflx_i,drefflx_i,dabsflx_i
c
c****   fluxes and radiances at output levels k_out
c
      real ups0(mxlout,2,nsol),dns0(mxlout,2,nsol),dirs0(mxlout,2,nsol),
     -     dups(mxlout,nsol),ddns(mxlout,nsol),ddirs(mxlout,nsol)
      save ups0,dns0,dirs0,dups,ddns,ddirs
c
      real upt0(mxlout,2),dnt0(mxlout,2),dupt(mxlout),ddnt(mxlout)
      save upt0,dnt0,dupt,ddnt
c
c****    flux transmission and absorption functions at wavenuber, wn
c
      real trn_dir(mxrad),trn_flx(mxrad),ref_flx(mxrad),abs_flx(mxrad)
c
c****     perturrbed monochormatic radiance and wn interpolation values
c 
      double precision sol_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,2,nsol),
     -                 dsol_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,nsol)
      double precision th_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,2),
     -                 dth_rad1(mxumu,mxphi,mxlout,mxrad,mxpd)
c
      save sol_rad1,dsol_rad1,th_rad1,dth_rad1
c
c****    flux transmission and absorption partial 
c        derivatives for simplified adding method at wavenumber wn
c
      real pd_trndir(mxrad,mxpd),
     -     pd_trnflx(mxrad,mxpd),
     -     pd_refflx(mxrad,mxpd),
     -     pd_absflx(mxrad,mxpd)
c
c****    optical depth variables for the output wavenumber grid
c
      real pray0(mxumu,2),pgas0(mxumu,2),paer0(mxumu,2),
     -     dpray(mxumu),dpgas(mxumu),dpaer(mxumu),
     -     tau_ray0(2),tau_gas0(2),tau_aer0(2),
     -     pray_00(2),pgas_00(2),paer_00(2)
c
      save pray0,pgas0,paer0,dpray,dpgas,dpaer,
     -     tau_ray0,tau_gas0,tau_aer0,
     -     pray_00,pgas_00,paer_00
c
      real dpray0,dpgas0,dpaer0
c
      save dpray0,dpgas0,dpaer0
c
      real p_ray0(mxumu),p_gas0(mxumu),p_aer0(mxumu),
     -     tau_ray_0,tau_gas_0,tau_aer_0,
     -     p_ray_0,p_gas_0,p_aer_0
c
      real dtau_ray,dtau_gas,dtau_aer
      save dtau_ray,dtau_gas,dtau_aer
c
c*****   variables passed to int_rad_lev
c
      real ups(mxlout),dns(mxlout),dirs(mxlout),
     -        upth(mxlout),dnth(mxlout)
      real rad_s(mxumu,mxphi,mxlout),rad_th(mxumu,mxphi,mxlout)
c
c****    output fluxes and radiances
c
      real dir_s_flx(mxrad,nsol,2),up_s_flx(mxrad,nsol,2),
     -     dn_s_flx(mxrad,nsol,2),up_t_flx(mxrad,2),dn_t_flx(mxrad,2)
c
      real up_flx(mxlout),dn_flx(mxlout),
     -     dir_flx(mxlout),rad(mxumu,mxphi,mxlout)
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
      double precision wnio_old,d_wn,wnavg
c
      real dvi_s,dwn_s
c
c****    specify the solar zenith angle and gas index
c
      nz = nz0
c
      ne = ne0
      nlev = nlyr+1   
c
      if(nsiout(2,ne,nz) .eq. 0) then
c
c****    initialize spectral interval counters
c
        nsiout(1,ne,nz) = 1
        nsiout(2,ne,nz) = 2
c
c****   initialize spectral radiance interpolation variables
c
        ibin = -1
        iend_ne = 0
        istart = 0
        wn_io = -9999.0d0
c
        call init_interp(nlev,nlout,numu,nphi,nz0,wnmin,wnout,
     -                   upsflx,dnsflx,dirsflx,sol_rad,
     -                   uptflx,dntflx,th_rad,
     -                   up_s_src_i,dn_s_src_i,dup_s_src,ddn_s_src,
     -                   up_t_src_i,dn_t_src_i,dup_t_src,ddn_t_src,
     -                   upsflx_i,dnsflx_i,dirsflx_i,dupsflx,
     -                   ddnsflx,ddirsflx,sol_rad0,dsol_rad,
     -                   uptflx_i,dntflx_i,duptflx,ddntflx,
     -                   th_rad0,dth_rad,ups,dns,dirs,
     -                   upth,dnth,rad_s,rad_th,pray0,pgas0,paer0,
     -                   pray_0,pgas_0,paer_0,tau_ray,tau_gas,tau_aer,
     -                   dpray,dpgas,dpaer,tau_ray0,tau_gas0,tau_aer0,
     -                   pray_00,pgas_00,paer_00,dtau_ray,
     -                   dtau_gas,dtau_aer,dpray0,dpgas0,dpaer0,
     -                   pgas,pray,paer)
c
      endif
c
c****  end of initialization block
c
c****  determine if new spectral radiances and fluxes are needed 
c      (ifl = 1) or if the existing ones span the range of including
c      wavenumber wn, and can be used to interpolate the values
c
      if(ifl .eq. 1) then
c
c****   read the binned optical properties at the next wavelength
c
c****    s p e c t r a l    i n t e r v a l    l o o p
c
2002    ni = nsiout(1,ne,nz)
        nsiout(1,ne,nz) = nsiout(2,ne,nz)
        nsiout(2,ne,nz) = ni
c
c****     read the next wavenumber, group index, and surface albedo. 
c
          wnio_old = wn_io
          read(iugrp,end=2601,err=6001) wn_io,ibin,
     -           (surf_opt(ii),ii=1,4),(alb(nn),nn=1,nza),
     -           ((dtauex(k,l),k=1,nlay),l=1,nt_pd),
     -           ((co_pi0(k,l),k=1,nlay),l=1,nt_pd),
     -           (g(k),k=1,nlyr),
     -           (((tau_ext(k,m,n),k=1,nlev),m=1,3),n=1,ntau_pd),
     -           ((tau_sca(k,n),k=1,nlev),n=1,ntau_pd),
     -           ((g_sca(k,n),k=1,nlev),n=1,ntau_pd)
          if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -           .or. irad .eq. 8) read(iugrp,end=2601,err=6001)
     -           pray_0,(pray(nze),nze=nzup,numu),tau_ray,
     -           pgas_0,(pgas(nze),nze=nzup,numu),tau_gas,
     -           paer_0,(paer(nze),nze=nzup,numu),tau_aer
c
        go to 2621
c
2601    iend_ne = 1
        wn_io = wnmax
        ni = nsiout(2,ne,nz)
        write(*,'(1a,5(1pe15.8),4i5)')
     -   'End of iugrp file encountered:',
     -   wn,wnio_old,wnmax,wnout(1,ne,nz),wnout(2,ne,nz),
     -   nz,nsiout(1,ne,nz),nsiout(2,ne,nz),ni
c
        return
c
c*****    set the radiance wavelength marker to wn_io (which is 
c         presumably beyond the current wavelength, wn)
c
2621    wnout(ni,ne,nz) = wn_io
c
c***      define wavenumber offset of wn
c
        d_wn = wn - wnout(1,ne,nz)
        dwn_s = real(d_wn)
c
c****     find the spectral interval width
c
        dv = wnout(2,ne,nz) - wnout(1,ne,nz)
        wnavg = 0.5d0*(wnout(2,ne,nz) + wnout(1,ne,nz))
c
        if(abs(dv) .ge. wn_tol*wnavg) then
          dvi = 1.0d0/dv
        else
          dvi = 0.0d0
        endif
        dvi_s = real(dvi)
c
c****     load tau=1 arrays for spectral interpolation
c
        if(nz .eq. 1) then
          if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -           .or. irad .eq. 8) then
            tau_gas0(ni) = tau_gas
            tau_ray0(ni) = tau_ray
            tau_aer0(ni) = tau_aer
            pgas_00(ni) = pgas_0
            pray_00(ni) = pray_0
            paer_00(ni) = paer_0
            do 2801 nze=1,numu
                pgas0(nze,ni) = pgas(nze)
                pray0(nze,ni) = pray(nze)
                paer0(nze,ni) = paer(nze)
2801        continue
c
c****        find optical depth and p(tau=1) spectral gradients
c    
            dtau_gas = (tau_gas0(2) - tau_gas0(1))*dvi_s
            dtau_ray = (tau_ray0(2) - tau_ray0(1))*dvi_s
            dtau_aer = (tau_aer0(2) - tau_aer0(1))*dvi_s
            dpgas0 = (pgas_00(2) - pgas_00(1))*dvi_s
            dpray0 = (pray_00(2) - pray_00(1))*dvi_s
            dpaer0 = (paer_00(2) - paer_00(1))*dvi_s
            do 2821 nze=1,numu
                dpgas(nze) = (pgas0(nze,2) - pgas0(nze,1))*dvi_s
                dpray(nze) = (pray0(nze,2) - pray0(nze,1))*dvi_s
                dpaer(nze) = (paer0(nze,2) - paer0(nze,1))*dvi_s
2821        continue
c
          endif
c
        endif
c
c****         evaluate the radiances and partial derivatives
c             at this wavenumber for each optical property 
c             that can vary
c
        do 3001 k=1,nlyr
            tatm(k) = t(k)
            dtau(k) = dtauex(k,1)
            copi0(k) = co_pi0(k,1)
            g0(k) = g(k)
3001    continue
c
        tatm(nlev) = t(nlev)
        tsurf = ts
        alb0 = alb(nz)
c
c****       set layer index range, l_1 and L_2
c
        l_1 = 1
        l_2 = nlev
        ng0 = ibin
        if(lsolar) then
c
c****       find the direct beam transmission for sun
c
          do 3011 k=1,nlyr
              trndir(k) = exp(-dtau(k)/umu0nz)
3011      continue
c
        endif
c
c****  set the surface pressure interpolation index
c
        if(ibin .gt. 0) then
c
c****       calculate radiances for wavenumber wn_io
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
        else
c
c****       use the pure absorption model to calculate the radiances
c           at wn_io
c
          fbeam = 1.0
          temis = 1.0
          btemp = tsurf
          ttemp = 0.0
          l_1 = 1
          l_2 = nlev
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
        endif
c
c****  load arrays for wavelength interpolation
c
        do 3021 k=1,nlev
            trnflx_i(k,ni) = trnflx(k)
            refflx_i(k,ni) = refflx(k)
            absflx_i(k,ni) = absflx(k)
3021    continue
c
        if(lsolar) then
          trndir_i(1,ni,nz) = trndir(1)
          do 3031 k=2,nlev
              trndir_i(k,ni,nz) = trndir(k)*trndir_i(k-1,ni,nz)
3031      continue
        endif
c
c****     find the spectral gradients in the flux transmission and
c         absorption values used for the simplified adding method
c
        do 3101 k=1,nlev
c
c****         define solar flux wavelength derivatives
c
            dtrndir_i(k,nz) = (trndir_i(k,2,nz) - 
     -                           trndir_i(k,1,nz))*dvi 
            dtrnflx_i(k) = (trnflx_i(k,2) - trnflx_i(k,1))*dvi
            drefflx_i(k) = (refflx_i(k,2) - refflx_i(k,1))*dvi
            dabsflx_i(k) = (absflx_i(k,2) - absflx_i(k,1))*dvi
3101    continue
c
c****       interpolate flux and radiance values to their output levels
c
        call int_rad_lev(lsolar,lplanck,nphi,numu,nz,
     -                       nlyr,nlout,levout,k_out,dp_dp,
     -                       upsflx,dnsflx,dirsflx,uptflx,dntflx,
     -                       sol_rad,th_rad,ups,dns,dirs,
     -                       upth,dnth,rad_s,rad_th)
c
c****       find the spectral gradients in the fluxes and
c           flux source functions
c
c
        if(lsolar) then
          do 3121 k=1,nlev
c
c****         define layer dependent solar source terms at wn_io
c
              up_s_src_i(k,ni,nz) = up_s_src(k) 
              dn_s_src_i(k,ni,nz) = dn_s_src(k)
c 
c****         define solar source term wavelength derivatives
c
              dup_s_src(k,nz) = (up_s_src_i(k,2,nz) -
     -                           up_s_src_i(k,1,nz))*dvi_s 
              ddn_s_src(k,nz) = (dn_s_src_i(k,2,nz) - 
     -                           dn_s_src_i(k,1,nz))*dvi_s
3121      continue
c
          do 3201 k=1,nlev
c
c****           define solar fluxes at wn_io
c
              upsflx_i(k,ni,nz) = upsflx(k) 
              dnsflx_i(k,ni,nz) = dnsflx(k)
              dirsflx_i(k,ni,nz) = dirsflx(k)
c
c****           define solar flux wavelength derivatives
c
              dupsflx(k,nz) = (upsflx_i(k,2,nz) - 
     -                         upsflx_i(k,1,nz))*dvi_s
              ddnsflx(k,nz) = (dnsflx_i(k,2,nz) - 
     -                         dnsflx_i(k,1,nz))*dvi_s
              ddirsflx(k,nz) = (dirsflx_i(k,2,nz) - 
     -                          dirsflx_i(k,1,nz))*dvi_s
3201      continue
c
          do 3261 k=1,nlout
c
c****           load fluxes at the output levels, k_out, and 
c               define solar flux wavelength derivatives
c
              ups0(k,ni,nz) = ups(k)
              dups(k,nz) = (ups0(k,2,nz) - ups0(k,1,nz))*dvi_s
              dns0(k,ni,nz) = dns(k)
              ddns(k,nz) = (dns0(k,2,nz) - dns0(k,1,nz))*dvi_s
              dirs0(k,ni,nz) = dirs(k)
              ddirs(k,nz) = (dirs0(k,2,nz) - dirs0(k,1,nz))*dvi_s
        
              do 3241 naz=1,nphi
                  do 3221 nze=1,numu
c
c****                   define solar radiances and wavelength 
c                       derivatives at output levels, k_out
c
                      sol_rad0(nze,naz,k,ni,nz) = rad_s(nze,naz,k)
                      dsol_rad(nze,naz,k,nz) = 
     -                        (sol_rad0(nze,naz,k,2,nz) - 
     -                         sol_rad0(nze,naz,k,1,nz))*dvi_s
3221              continue
3241          continue
3261      continue

        endif
c
        if(lplanck .and. nz .eq. 1) then      
c
          do 3321 k=1,nlev
c
c****         define layer dependent solar source terms at wn_io
c
              up_t_src_i(k,ni) = up_t_src(k) 
              dn_t_src_i(k,ni) = dn_t_src(k)
c 
c****         define solar source term wavelength derivatives
c
              dup_t_src(k) = (up_t_src_i(k,2) -
     -                           up_t_src_i(k,1))*dvi_s 
              ddn_t_src(k) = (dn_t_src_i(k,2) - 
     -                           dn_t_src_i(k,1))*dvi_s
3321      continue
c
          do 3401 k=1,nlev
c
c****           define thermal fluxes at wn_io
c
              uptflx_i(k,ni) = uptflx(k)
              dntflx_i(k,ni) = dntflx(k)
c
c****           define thermal flux wavelength derivatives
c
              duptflx(k) = (uptflx_i(k,2) - uptflx_i(k,1))*dvi_s
              ddntflx(k) = (dntflx_i(k,2) - dntflx_i(k,1))*dvi_s
3401      continue  
c
          do 3461 k=1,nlout
c
c****           load fluxes at the output levels, k_out, and 
c               define solar flux wavelength derivatives
c
              upt0(k,ni) = upth(k)
              dupt(k) = (upt0(k,2) - upt0(k,1))*dvi_s
              dnt0(k,ni) = dnth(k)
              ddnt(k) = (dnt0(k,2) - dnt0(k,1))*dvi_s
        
c
c****           define thermal radiance wavelength derivatives
c
              do 3441 naz=1,nphi
                  do 3421 nze=1,numu
c
c****                   define thermal radiances and wavelength 
c                       derivatives at output levels, k_out
c
                      th_rad0(nze,naz,k,ni) = rad_th(nze,naz,k)
                      dth_rad(nze,naz,k) = (th_rad0(nze,naz,k,2) - 
     -                                      th_rad0(nze,naz,k,1))*dvi_s
3421              continue
3441          continue
3461      continue  
c
        endif
c
        if(npd .gt. 0) then 
c
c****       initialze optical properties, flux and radiance transmission
c           and absorption functions in each layer
c
          dalb0 = dalb
          do 3631 k=1,nlev-1
              deltau0(k) = deltau(k)
              delpi00(k) = delpi0(k)
              delg0(k) = delg(k) 
              trndir0(k) = trndir(k)
              trnflx0(k) = trnflx(k)
              refflx0(k) = refflx(k)
              absflx0(k) = absflx(k)
              do 3601 nze=1,numu
                  trnrad0(nze,k) = trnrad(nze,k)
                  refrad0(nze,k) = refrad(nze,k)
                  absrad0(nze,k) = absrad(nze,k)
3601          continue
3631      continue
c
          if(lsolar) then
            do 3681 k=1,nlev
                dn_s_src0(k) = dn_s_src(k)
                up_s_src0(k) = up_s_src(k)
                do 3661 naz=1,nphi
                    do 3641 nze=1,numu
                        rad_s_src0(nze,naz,k) = rad_s_src(nze,naz,k)
3641                continue
3661            continue
3681        continue
          endif
c
          if(lplanck) then
            bb_sur0 = bb_sur
            do 3781 k=1,nlev
                bb0(k) = bb(k)
                bb_flx_up0(k) = bb_flx_up(k)
                bb_flx_dn0(k) = bb_flx_dn(k)
                dn_t_src0(k) = dn_t_src(k)
                up_t_src0(k) = up_t_src(k)
                do 3721 nze=1,numu
                    do 3701 naz=1,nphi
                        bb_rad0(nze,naz,k) = bb_rad(nze,naz,k)
                        rad_t_src0(nze,naz,k) = rad_t_src(nze,naz,k)
3701                continue
3721            continue
3781        continue
          endif
c
c****        find flux and radiance jacobians at wavenumer wn_io
c
          ni0 = ni
c      
          call perturb_rad(lsolar,lplanck,lamber,usrang,ni0,nz0,
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
        endif
c
c******  if this is is the first point, get the second one for the 
c        spectral interpolation
c
        if(istart .eq. 0) then
          istart = 1
          go to 2002
        endif
c
c****     ensure that at least one valid radiance point beyond wn has  
c         been found.  if not, read the next point.
c
        if(wnout(ni,ne,nz) .lt. wn .and. wnout(ni,ne,nz) .lt. wnmax
     -     .and. iend_ne .ne. 1) go to 2002
c
      endif
c
      if(nz .eq. 1) then
        if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -          .or. irad .eq. 8) then
c
c****       save optical depths and pressure of tau = 1
c
          tau_gas_0 = tau_gas0(1) + dtau_gas*dwn_s
          tau_ray_0 = tau_ray0(1) + dtau_ray*dwn_s
          tau_aer_0 = tau_aer0(1) + dtau_aer*dwn_s
          p_gas_0 = pgas_00(1) + dpgas0*dwn_s
          p_ray_0 = pray_00(1) + dpray0*dwn_s
          p_aer_0 = paer_00(1) + dpaer0*dwn_s
          do 4001 nze=1,numu
              p_gas0(nze) = pgas0(nze,1) + dpgas(nze)*dwn_s
              p_ray0(nze) = pray0(nze,1) + dpray(nze)*dwn_s
              p_aer0(nze) = paer0(nze,1) + dpaer(nze)*dwn_s
4001      continue
c
        endif
c
      endif
c
c****           interpolate flux transmision and absorption functions 
c               to wavenumber, wn
c
      do 4021 k=1,nlev
          if(lsolar) trn_dir(k) = 
     -                 real(trndir_i(k,1,nz) + dtrndir_i(k,nz)*d_wn)
          trn_flx(k) = real(trnflx_i(k,1) + dtrnflx_i(k)*d_wn)
          ref_flx(k) = real(refflx_i(k,1) + drefflx_i(k)*d_wn)
          abs_flx(k) = real(absflx_i(k,1) + dabsflx_i(k)*d_wn)
4021  continue
c
      if(lsolar) then
c
c****     Interpolate solar fluxes, flux source functions, and
c          radiances to wavenumber, wn
c
        do 4101 k=1,nlyr+1
c
c****        interpolate solar source functions to this wn
c
            ups_src(k) = real(up_s_src_i(k,1,nz) + dup_s_src(k,nz)*d_wn)
            dns_src(k) = real(dn_s_src_i(k,1,nz) + ddn_s_src(k,nz)*d_wn)
c
c****        interpolate solar fluxes and radiances to this wn
c            Note: fluxes must be retained at wn adjacent wavenumbers
c            for the trapezoid integration over wavenumber in the 
c            heating rate calculation
c
            up_s_flx(k,nz,2) = solflx*(upsflx_i(k,1,nz) + 
     -                           dupsflx(k,nz)*dwn_s)
            dn_s_flx(k,nz,2) = solflx*(dnsflx_i(k,1,nz) + 
     -                           ddnsflx(k,nz)*dwn_s)
            dir_s_flx(k,nz,2) = solflx*(dirsflx_i(k,1,nz) + 
     -                            ddirsflx(k,nz)*dwn_s)
c
            if(up_s_flx(k,nz,2) .lt. 0.0) up_s_flx(k,nz,2) = 0.0
            if(dn_s_flx(k,nz,2) .lt. 0.0) dn_s_flx(k,nz,2) = 0.0
            if(dir_s_flx(k,nz,2) .lt. 0.0) dir_s_flx(k,nz,2) = 0.0
4101    continue
c
c****        interpolate solar fluxes and radiances at levels nlout
c            to this wn
c
        do 4161 k=1,nlout
c
c****         find the solar fluxes at levels k_out
c
            ups(k) = solflx*(ups0(k,1,nz) + dups(k,nz)*dwn_s)
            dns(k) = solflx*(dns0(k,1,nz) + ddns(k,nz)*dwn_s)
            dirs(k) = solflx*(dirs0(k,1,nz) + ddirs(k,nz)*dwn_s)
c
c****         find radiances at levels k_out
c
            do 4141 naz=1,nphi
                do 4121 nze=1,numu
                    rad_s(nze,naz,k) =
     -                       solflx*(sol_rad0(nze,naz,k,1,nz) +
     -                       dsol_rad(nze,naz,k,nz)*dwn_s)
c
                    if(rad_s(nze,naz,k) .lt. 0.0) 
     -                 rad_s(nze,naz,k) = 0.0
4121            continue
4141        continue
4161    continue
c
      endif
c
c****       t h e r m a l    f l u x e s    a n d    r a d i a n c e s 
c           (do only for nz = 1)
c
      if(lplanck .and. nz .eq. 1) then
c
c****     interpolate thermal fluxes and radiances to this wn
c
        do 4201 k=1,nlyr+1
c
c****        interpolate solar source functions to this wn
c
            upt_src(k) = real(up_t_src_i(k,1) + dup_t_src(k)*d_wn)
            dnt_src(k) = real(dn_t_src_i(k,1) + ddn_t_src(k)*d_wn)
c
c****     interpolate thermal fluxes to this wn
c
            up_t_flx(k,2) = uptflx_i(k,1) + duptflx(k)*dwn_s
            dn_t_flx(k,2) = dntflx_i(k,1) + ddntflx(k)*dwn_s
c
            if(up_t_flx(k,2) .lt. 0.0) up_t_flx(k,2) = 0.0
            if(dn_t_flx(k,2) .lt. 0.0) dn_t_flx(k,2) = 0.0
c
4201    continue
c
c****     interpolate thermal fluxes and radiances at levels
c         k_out to this wn
c
        do 4261 k=1,nlout
c
c****         interpolate thermal fluxes at levels k_out to this wn
c
            upth(k) = upt0(k,1) + dupt(k)*dwn_s
            dnth(k) = dnt0(k,1) + ddnt(k)*dwn_s
        
            do 4241 naz=1,nphi
                do 4221 nze=1,numu
                    rad_th(nze,naz,k) = th_rad0(nze,naz,k,1) +
     -                                     dth_rad(nze,naz,k)*dwn_s
                    if(rad_th(nze,naz,k) .lt. 0.0) 
     -                 rad_th(nze,naz,k) = 0.0
4221            continue
4241        continue
4261    continue
c
        write(iuthrm) nsiout(1,ne,nz),nsiout(2,ne,nz),
     -                    wnout(nsiout(1,ne,nz),ne,nz),
     -                    wnout(nsiout(2,ne,nz),ne,nz),
     -                    wn_io,(upth(k),dnth(k),
     -                    ((rad_th(nze,naz,k),nze=1,numu),
     -                     naz=1,nphi),k=1,nlout)
c
      else
c
c****      nz ne 1.  Read stored values
c
        if(lplanck) then
          read(iuthrm) nsiout(1,ne,nz),nsiout(2,ne,nz),
     -                   wnout(nsiout(1,ne,nz),ne,nz),
     -                   wnout(nsiout(2,ne,nz),ne,nz),
     -                   wn_io,(upth(k),dnth(k),
     -                   ((rad_th(nze,naz,k),nze=1,numu),
     -                   naz=1,nphi),k=1,nlout)
        endif
c
      endif
c
c****    c o m b i n e d   f l u x e s   a n d   r a d i a n c e s
c
c****     define the combined solar+thermal upward and downward 
c         fluxes at the output levels
c
      do 4641 k=1,nlout
c
c****               store fluxes at top or bottom of atmosphere:
c
            up_flx(k) = upth(k) + ups(k)
            dn_flx(k) = dnth(k) + dns(k)
            dir_flx(k) = dirs(k)
c
c****          define the combined solar+thermal radiances at the
c              output level, 
c    
          do 4621 naz=1,nphi
              do 4601 nze=1,numu
                    rad(nze,naz,k) = rad_s(nze,naz,k) +
     -                               rad_th(nze,naz,k)
4601          continue
4621      continue
4641  continue
c
      if(npd .gt. 0) then 
c
c****      find partial derivatives
c
        ni0 = ni
c      
        call find_pd(lsolar,lplanck,nz0,nlyr,
     -                   nlout,nphi,numu,npd,ipd,
     -                   wn,d_wn,dx,solflx,rad,
     -                   trn_dir,trn_flx,ref_flx,abs_flx,
     -                   dns_src,ups_src,dnt_src,upt_src,
     -                   trndir1_i,trnflx1_i,refflx1_i,absflx1_i,
     -                   dtrndir1_i,dtrnflx1_i,drefflx1_i,dabsflx1_i,
     -                   dn_s_src1,up_s_src1,ddn_s_src1,dup_s_src1,
     -                   dn_t_src1,up_t_src1,ddn_t_src1,dup_t_src1,
     -                   sol_rad1,dsol_rad1,th_rad1,dth_rad1,
     -                   pd_trndir,pd_trnflx,pd_refflx,pd_absflx,
     -                   pd_dns_src,pd_ups_src,pd_dnt_src,pd_upt_src,
     -                   pd_rad)
c
      endif
c
      return
6001  write(*,'(/,1a,1pe14.6)') 
     - 'Input Error on unit iugrp in subroutine map_rad after: wn =',
     -  wnio_old
      stop
c
      end
