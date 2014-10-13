      program smart_spectra
c
cccccccccccccccccccccc    s m a r t _ s p e c t r a  ccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This program combines spectral mapping methods and the          cc
cc    discrete ordinate model (Stamnes et al. 1988) to generate       cc
cc    high-resolution synthetic monochromatic radiances in            cc
cc    vertically-inhomogeneous scattering absorbing, emitting,        cc
cc    planetary atmospheres.                                          cc
cc                                                                    cc
cc    note:  This version of the program uses albedo gradients.       cc
cc           It also outputs wavelength-dependent, level-dependent    cc
cc           contribution functions.                                  cc
cc                                                                    cc
cc   9/28/97:   this version prints the direct solar flux as well as  cc
cc              the total downward and upward fluxes, and allows      cc
cc              wavelength-dependent radiances to be printed          cc
cc              at mxlout levels.                                     cc
cc                                                                    cc
cc              as part of this modification, the slit function       cc
cc              routine, slit, was modified such that it only allows  cc
cc              boxcar or triangular slits (other functions required  cc
cc              far too much RAM.                                     cc
cc                                                                    cc
cc    5/3/98:   the aerosol, gas, and rayliegh scattering optical     cc
cc              depth sections have been extracted from the main      cc
cc              program and are called as subroutines.                cc
cc                                                                    cc
cc    6/3/98:   smart.f was modified to include a 0th order           cc
cc              correction for the solar zenith angle to account for  cc
cc              the sphericity of the planet.                         cc
cc                                                                    cc
cc    5/28/01   smart modified to evaluate radiances for the          cc
cc              minimum and maximum optical depth within each         cc
cc              bin, and interpolate linearly in tau                  cc
cc                                                                    cc
cc    5/29/03:  smart modified to use DISORT 2.0                      cc
cc                                                                    cc
cc    12/26/03: smart modified to provide an option to elminate       cc
cc              headers from binary out files (ifrmout>2)             cc
cc    11/14/04: bdrf mode (iref) propagated through program           cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      iuabc - unit numbers for gas line parameters                  cc
cc     iusol0 - unit number as a scratch file for solar fluxes        cc
cc     iusol1 - unit number for solar fluxes                          cc
cc     iusol2 - unit number for solar flux scratch file               cc
cc     iuthrm - unit number for thermal radiances                     cc
cc      iuout - unit number of output radiance file                   cc
cc      iutrn - unit number of output transmission/pressure file      cc
cc      iutpd - unit numbers for output level-dependent               cc
cc              thermal fluxes and thier partial derivatives          cc
cc      iuspd - unit numbers for output level-dependent               cc
cc              solar fluxes and thier partial derivatives            cc
cc    iupdrad - unit numbers for output level-dependent radiances     cc
cc              and thier partial derivatives                         cc
cc    atmfile - name of input atmospheric structure file              cc
cc    aerfile - name of input aerosol vertical structure file         cc
cc    miefile - name of file with aerosol optical properties vs. wn   cc
cc    solfile - name of file with wn-dependent solar fluxes           cc
cc    surfile - name of file with wn-dependent surface albedos        cc
cc    mixfile - name of file with gas mixing ratios                   cc
cc    gasfile - name fo file with gas absorption coeffiecients vs. wn cc
cc     nmodes - number of discrete aerosol partical modes             cc
cc      ncomp - number of rayleigh-scattering constituentes           cc
cc      icomp - index of each rayleigh scattering constituent         cc
cc     volmix - volume mixing ratio of each rayleigh scatterer        cc
cc         ts - surface temperature (K)                               cc
cc         au - distance to the sun (in AU's)                         cc
cc     tauerr - optical depth relative binning error (0. to ~0.8)     cc
cc     pi0err - co-single scattering albedo absolute binning error    cc
cc     phferr - asymmetry factor absolute binning error               cc
cc    surferr - surface albedo absolute binning error                 cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc       pout - output pressure level (bars)                          cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, use a BRDF is used.              cc
cc    isource - index of source function type: (1) solar only         cc
cc              (2) thermal only, (3) both                            cc
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
cc    ifmrout - index of output file format (1) ascii, (2) binary,    cc
cc                                          (3) binary, no header     cc
cc       iref - bidirectional reflectance options                     cc
cc              0 - Lambert                                           cc
cc              1 - Hapke's BDR model                                 cc
cc              2 - Breon's BDR model; combination of Li + Roujean    cc
cc              3 - Roujean's BDR model                               cc
cc              4 - Cox and Munk glint model                          cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) Watts/m**2/sr/Angstrom                             cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       nphi - number of output azimuth angles                       cc
cc       nlev - number of output levels                               cc
cc       nlyr - number of computational model layers                  cc
cc     levout - output level index (1) top of atmosphere,             cc
cc              (2) surface, (3) arbitrary level                      cc
cc        nza - number of solar zenith angles                         cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc       nref - number of surface optical properties specified at     cc
cc              each wavelength.                                      cc
cc  SURF_PR : Wavelength dependent surface properties array           cc
cc            IREF= 0 - Lambert albedo                                cc
cc            IREF= 1 - Hapke : HH, W                                 cc
cc            IREF= 2 - Breon's BDR model: k0, k1, k2                 cc
cc            IREF= 3 - Roujean's BDR model: k0, k1, k2               cc
cc            IREF= 4 - Cox and Munk glint model: n, k, ws, phiw      cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc        phi - emission azimuth angles (degrees)                     cc
cc       umu0 - solar zenith angle cosines                            cc
cc       phi0 - solar azimuth angle cosines                           cc
cc      wnmin - minimum wavenumber of desired spectral window         cc
cc      wnmax - maximum wavenumber of desired spectral window         cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc         wl - current wavenumber (cm**-1)                           cc
cc         wn - current wavenumber (cm**-1)                           cc
cc    dir_flx - direct downward solar flux at wavenumber, wn          cc
cc     dn_flx - total downward flux (irradiance) at wavenumber wn     cc
cc     up_flx - upward flux (irradiance) at wavenumber wn             cc
cc        rad - radiance at each output zenith and azimuth angle      cc
cc  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc
cc  trn_gas_0 - normal-incidence gas transmission                     cc
cc  trn_ray_0 - normal-incidence aerosol transmission                 cc
cc                                                                    cc
cccccccccccccccccccccc    s m a r t _ s p e c t r a  ccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      character*132 atmfile,io_file(nex),mixfile(ngas),
     -              gasfile(ngtmax,ngas),
     -              aerfile(nmode),miefile(nmode),solfile,surfile
      character*40 frmmix(ngas),frmatm,frms0,frmsur
      character*1  name(132)
c
c****   define variables used by discrete ordinate method
c
      logical usrang,lamber,lplanck,lsolar
      logical lcsh(nmode)
c
      integer iuatm, iumix,iumie, iusur
      integer iustat,iuaer,iugas,iuheat
      integer iugrp,iutrn,iuflx,iusol0,iusol1,iusol2,iuthrm,iuout
      integer iuabc(ngtmax,ngas)
c
      integer new_pt,new_mix
      integer iskatm,ifrmatm,icp,ict,nmodes,icwnsur
      integer ifrms0,ioffs0,ixs0,icwns0,icsol,iuins0
      integer ncomp,nlev,nlyr,k,ng,ngt,m,mode,iopen,len,ii,nlb,ntb
      integer next,n,nalbmx,isource
c      
      integer numu,nstr,nphi,k_out(mxlout)
      integer icomp(6),ne
      integer igas(ngas),ngases,ngastyp(ngas)
      integer igtrn(ngtmax,ngas),igpab(ngtmax,ngas),
     -        ifrmmix(ngas),ioffmix(ngas),izmix(ngas),imix(ngas),
     -        icpmix(ngas),icmix(ngas)
c
      integer ioffmie(nmode),icwlq(nmode),icqext(nmode),icqsca(nmode),
     -        icg1(nmode),icg2(nmode),icf(nmode),iofftau(nmode),
     -        iztau(nmode),icptau(nmode),ictau(nmode),
     -        iccsh(nmode),ioffmom(nmode),iang(nmode)
c
      integer io_err(nex),io_end(nex)
      integer ifrmsur,ioffsur,iwnsur,icalb(4)
      integer ifrmout,iunits,levout,nlout,irad,
     -        nza,isptype,islit
      integer iref,nref
c
      integer nze
c
c****    units for flux and radiance jacobians
c 
      integer nstate,istate(nex),ntau_pd,igs(mxpd),
     -        iu_pd,iutpd(mxpd,nsol),iuspd(mxpd,nsol),
     -        iupdrad(mxumu,mxphi,mxpd,mxlout,nsol)
c
      double precision wn_tol,wnmin,wnmax,width,dwn
      double precision wn_eof(2,nex),wneof(2)
      double precision wnaer(nmode)
      real dtauaer(nmode,kp)
      real scptau(nmode),sctau(nmode),sccsh(nmode)
      real pi,scp,scwalb,scalb(4),ws,phiw,scwns0,scsol
      real ratm,tau_tot,dtau_tot,pbar,ang
      real phi(mxphi),umu(mxumu),phi0(nsol),umu0(nsol)
      real scpmix(ngas),scmix(ngas),rmix(kp,ngas),wgtgas(ngas)
      real volmix(6),ts,au,wgtatm
      real tauerr,pi0err,phferr,surferr
      real pout(mxlout),dp_dp(mxlout),accur
      real pd_frac(nex)
      real umu_o(mxumu),gwt_o(mxumu),umu_1(mxumu)
c
c****   atmospheric structure variables
c
      real p(kp),t(kp),alt(kp),grav(kp),z(kp),radius,sgrav
c
c****   double precision variables for heating rates
c
      double precision aid_lr(mxulv),dirsoflx(mxulv,nsol),
     -       dnsoflx(mxulv,nsol),upsoflx(mxulv,nsol),dnthflx(mxulv),
     -       upthflx(mxulv),soheat(mxulv,nsol),thheat(mxulv)
c
c****   define the minimum spectral interval width:
c
      wn_tol = 5.0d-7
c
c****   define the input and output unit numbers for all i/o
c
      call init_iu(iuheat,iustat,iuatm,iumie,iumix,iusur,
     -             iugrp,iuthrm,iutrn,iuaer,iuflx,iuout,
     -             iugas,iusol0,iusol1,iusol2,iu_pd)
c
c****   read in all user-supplied input
c
      call smart_in
     -     (iuatm,atmfile,ifrmatm,frmatm,iskatm,icp,ict,scp,
     -     ngases,ngastyp,igas,iugas,iuabc,igtrn,igpab,wgtgas,
     -     gasfile,iumix,mixfile,ifrmmix,frmmix,ioffmix,
     -     izmix,imix,icpmix,icmix,scpmix,scmix,
     -     nmodes,iuaer,aerfile,iumie,miefile,ioffmie,
     -     ioffmom,iang,icwlq,icqext,icqsca,icg1,icg2,icf,wnaer,
     -     iofftau,iztau,icptau,ictau,scptau,sctau,lcsh,iccsh,sccsh,
     -     iusur,surfile,ifrmsur,frmsur,ioffsur,ws,phiw,
     -     iwnsur,icwnsur,icalb,scwalb,scalb,
     -     iusol1,solfile,ifrms0,frms0,ioffs0,ixs0,
     -     icwns0,icsol,scwns0,scsol,iuins0,
     -     ncomp,icomp,volmix,ts,au,radius,sgrav,wgtatm,
     -     wnmin,wnmax,isptype,islit,width,dwn,
     -     nstr,usrang,numu,umu,nphi,phi,nza,umu0,phi0,
     -     levout,nlout,pout,accur,lamber,lplanck,lsolar,
     -     isource,irad,iref,nref,iuout,iunits,ifrmout,iuheat,
     -     iustat,iutrn,iuflx,tauerr,pi0err,phferr,surferr,
     -     pd_frac,nstate,istate,iu_pd,iutpd,iuspd,iupdrad)
c
c****  read the atmospheric structure file and the gas mixing ratio
c      files and set up the model atmosphere:
c
      nlev = 0
      new_pt = 1
      new_mix = 1
      do 1001 k=1,kp
          p(k) = 0.0
          t(k) = 0.0
          alt(k) = 0.0
          z(k) = 0.0
          grav(k) = 0.0
1001  continue          
c
      call mod_atm(iuatm,iumix,new_pt,new_mix,istate,
     -            atmfile,ifrmatm,frmatm,iskatm,icp,ict,scp,
     -            nlev,nlyr,levout,nlout,k_out,pout,wgtatm,wgtgas,
     -            ngases,igas,mixfile,ifrmmix,frmmix,ioffmix,
     -            izmix,imix,icpmix,icmix,scpmix,scmix,pd_frac,
     -            radius,sgrav,ratm,p,t,alt,z,grav,rmix,dp_dp)
c
c****   define optical depth counters for 
c       rayleigh scattering (ne = 1) gas absorption
c
      ne = 1
      wn_eof(1,ne) = wnmin
      wn_eof(2,ne) = wnmax
c
c****  initialize the counter for perturbations in optical depth 
c
      ntau_pd = 0
c
c****   increment optical depth and optical depth perturbation 
c       counter, ntau_pd for each radiatively active gas that is a 
c       variable part of the state vector
c
      do 1821 ng=1,ngases
          if(iabs(istate(ng+2)) .eq. 3) then
            ntau_pd = ntau_pd + 1
            igs(ntau_pd) = ng
          endif          
          do 1801 ngt = 1,ngastyp(ng)
              ne = ne + 1
              wn_eof(1,ne) = -999.0d0
              wn_eof(2,ne) = -999.0d0
              io_file(ne) = gasfile(ngt,ng)
1801      continue
1821  continue
c
c****   define aerosol optical properties and vertical distribution
c
      do 2221 m=1,nmodes
          ne = ne + 1
          io_file(ne) = miefile(m)
          wn_eof(1,ne) = -999.0d0
          wn_eof(2,ne) = -999.0d0
          mode = m
          iopen = 1
c
c*****       if this particle mode is a variable part of the state
c            vector, update the state vector counter and optical depth
c            counter
c
          if(iabs(istate(m+ngases+2)) .eq. 4) then
            ntau_pd = ntau_pd + 1
          endif          
c
c****      read wavelength-dependent optical properties of aerosols.
c
          call readmie(mode,iumie,iuaer,miefile(m),iang(m),ioffmie(m),
     -                 icwlq(m),icqext(m),icqsca(m),icg1(m),icg2(m),
     -                 icf(m),ioffmom(m),nstr,
     -                 io_end(ne),io_err(ne),wnaer(m),wneof)
c
          wn_eof(1,ne) = wneof(1)
          wn_eof(2,ne) = wneof(2)
c
c****      read the aerosol vertical structure
c
          call cldstr(nlev,iopen,mode,iuaer,aerfile(m),
     -                iofftau(m),iztau(m),icptau(m),ictau(m),
     -                scptau(m),sctau(m),lcsh(m),iccsh(m),sccsh(m),
     -                p,alt,z,dtauaer)
c
2221  continue
c
c****   print the integrated cloud optical depth
c
      write(*,'(/,1a,/,/,1a)') 
     -' Total Cloud and Aerosol Optical Depth at reference wavelength:',
     -'    alt(km)     p (bar)      T(K)       dtau        tau'
      tau_tot = 0.0
      do 2461 k=1,nlyr
          dtau_tot = 0.0
          do 2441 m=1,nmodes          
             dtau_tot = dtau_tot + dtauaer(m,k)
2441      continue
          pbar = 1.e-5*p(k+1)
          tau_tot = tau_tot + dtau_tot
          write(*,'(5(1pe12.4))') alt(k+1),pbar,t(k+1),dtau_tot,tau_tot
2461  continue
c
c****   read surface albedos and bi-directional reflection functions
c
      nalbmx = 32766
      iopen = 1
      ne = ne + 1
      wn_eof(1,ne) = -999.0d0
      wn_eof(2,ne) = -999.0d0
      io_file(ne) = surfile
c
      call readsur(surfile,iusur,iopen,ifrmsur,frmsur,ioffsur,
     -             nalbmx,iwnsur,icwnsur,nref,icalb,scwalb,scalb,
     -             wnmin,wnmax,io_end(ne),io_err(ne),wneof)
c
      wn_eof(1,ne) = wneof(1)
      wn_eof(2,ne) = wneof(2)
c
c****   read solar flux input file.
c
      if(lsolar) then
        ne = ne + 1
        wn_eof(1,ne) = -999.
        wn_eof(2,ne) = -999.
        iopen = 1
        io_file(ne) = solfile
c
        call readsol(solfile,iusol1,iusol2,iopen,ifrms0,frms0,
     -                ioffs0,ixs0,icwns0,icsol,scwns0,scsol,
     -                iuins0,wnmin,wnmax,au,
     -                io_end(ne),io_err(ne),wneof)
c
        wn_eof(1,ne) = wneof(1)
        wn_eof(2,ne) = wneof(2)
c
      else
        nza = 1
      endif
c
c****  save input file names and other data in a header file.
c  
      if(ifrmout .le. 2) then
c
c****      find gaussian points and weights. These quantities are 
c          needed to define the number of upward streams, nzup,
c          and the number of downward streams, nzdn,and
c          for flux integration for grey calculations
c
        call qgausn(nstr/2,umu_o,gwt_o)   
c
c****      restack the values into order used by disort
c
        do 2601 nze=1,nstr/2
            umu_1(nze+nstr/2) = umu_o(nze)
            umu_1(nze) = -umu_o(nstr/2-nze+1)
2601    continue
c
        call smart_hdr
     -     (iuout,iutrn,iuflx,atmfile,aerfile,miefile,solfile,
     -     surfile,mixfile,gasfile,iunits,nmodes,ncomp,icomp,volmix,
     -     ts,au,wgtatm,wnmin,wnmax,isptype,islit,width,dwn,
     -     tauerr,pi0err,phferr,surferr,nstr,numu,umu_1,
     -     nphi,phi,nza,umu0,phi0,levout,nlout,pout,accur,
     -     lamber,isource,irad,ifrmout,radius,sgrav,
     -     igas,ngases,ngastyp)
c
      endif
c
c****   initialze th number of sources of extinction
c
      next = 0
c
c****  Call the spectral mapping Algorithm
c
      call smart(iugrp,iuthrm,iuabc,iuaer,iusur,iusol0,
     -     iusol2,iuout,iuheat,iustat,iutrn,iuflx,nmodes,
     -     ngases,ngastyp,igtrn,igpab,gasfile,iref,nref,
     -     nlev,nlyr,next,k_out,dp_dp,ncomp,icomp,volmix,
     -     radius,wgtatm,ts,grav,p,t,alt,z,rmix,dtauaer,
     -     wnmin,wnmax,wn_tol,isptype,islit,width,dwn,
     -     usrang,nstr,numu,umu,nphi,phi,nza,umu0,phi0,
     -     levout,nlout,ws,phiw,lamber,lplanck,lsolar,irad,
     -     iunits,ifrmout,io_end,io_err,wn_eof,accur,ratm,
     -     io_file,tauerr,pi0err,phferr,surferr,aid_lr,
     -     dirsoflx,dnsoflx,upsoflx,dnthflx,upthflx,
     -     soheat,thheat,pd_frac,nstate,istate,ntau_pd,
     -     igs,iu_pd,iutpd,iuspd,iupdrad)
c
c****      Print final output
c
      write(*,'(/,/,1a,/,/,1a)') 
     - ' O u t p u t    R a d i a n c e    A n g l e s    U s e d: ',
     - '   i   ang(deg)     cos(ang)'
      pi = acos(-1.0)
      do 7001 n=1,numu
          ang = 180.*acos(umu(n))/pi
          write(*,'(1x,i3,2(1pe13.5))') n,ang,umu(n)
7001  continue    
c
c****    spectral ranges used for each wavelength-dependent file:
c
      write(*,'(/,/,1a,/,/,2a)')
     - ' I/O Status and Final Wavenumber for Each Constituent:',
     - '       n     io_end  io_err',
     - ' wavenumber extent (cm**-1)     filename'
      do 7101 ne=1,ngases + 1
c
          call charsp(io_file(ne),name,len,132,nlb,ntb) 
c
          write(*,'(3i8,2x,2(1pe14.6),5x,132a)') 
     -         ne,io_end(ne),io_err(ne),wn_eof(1,ne),
     -         wn_eof(2,ne),(name(ii),ii=1,len)
7101  continue
      do 7121 ne=ngases+2,next
c
          call charsp(io_file(ne),name,len,132,nlb,ntb) 
c
          write(*,'(3i8,2x,2(1pe14.6),5x,132a)') 
     -         ne,io_end(ne),io_err(ne),wn_eof(1,ne),
     -         wn_eof(2,ne),(name(ii),ii=1,len)
7121  continue
c
      if(lsolar) then
        ne = next + 1
c
        call charsp(io_file(ne),name,len,132,nlb,ntb) 
c
        write(*,'(3i8,2x,2(1pe14.6),5x,132a)') 
     -         ne,io_end(ne),io_err(ne),wn_eof(1,ne),
     -         wn_eof(2,ne),(name(ii),ii=1,len)
      endif
c
      stop
      end
     
