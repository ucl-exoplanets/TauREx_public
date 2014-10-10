      subroutine rad_out(lsolar,iuout,iutrn,iu_flux,
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
cccccccccccccccccccccccccc   r a d _ o u t   ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine prints backmapped radiances  and other          cc
cc    spectral quantities back to a high-resolution spectral grid.    cc
cc    It can print monochromatic variables, or smooth these variables cc
cc    with a slit function.                                           cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc      iuout - unit number of output radiance file                   cc
cc      iutrn - unit number of output transmission file               cc
cc    iu_flux - unit number for output level-dependent fluxes         cc
cc         nz - index of this solar zenith angle (1 to nsol)          cc
cc        nza - number of solar zenith angles                         cc
cc    ifrmout - output file format (1) ascii or (2) binary            cc
cc     levout - output level index (1) top of atmosphere,             cc
cc              (2) surface, (3) arbitrary level                      cc
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
cc      islit - index of slit function type:                          cc
cc              1) boxcar                                             cc
cc              2) triangular (approximate slit spectrometer)         cc
cc      width - half-width at half max of slit function  (cm**-1)     cc
cc        dwn - output sampling resolution (cm**-1)                   cc
cc         wn - current wavenumber (cm**-1)                           cc
cc     solflx - solar flux at this wavenumber                         cc
cc      wnmin - minimum wavenumber in output spectral grid (cm**-1)   cc
cc      wnmax - maximum wavenumber in output spectral grid (cm**-1)   cc
cc    dir_flx - downward direct flux at wavenumber, wn                cc
cc     dn_flx - total downward flux (irradiance) at wavenumber wn     cc
cc     up_flx - upward flux (irradiance) at wavenumber wn             cc
cc        rad - radiance at each output zenith and azimuth angle      cc
cc  tau_ray_0 - normal-incidence rayleigh-scattering optical depth    cc
cc  tau_gas_0 - normal-incidence gas absorption optical depth         cc
cc  tau_ray_0 - normal-incidence aerosol extinction optical depth     cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc         wl - current wavenumber (cm**-1)                           cc
cc       wnio - current wavenumber (cm**-1)                           cc
cc    dir_flx - downward direct flux at wavenumber, wn                cc
cc     dn_flx - total downward flux (irradiance) at wavenumber wn     cc
cc     up_flx - upward flux (irradiance) at wavenumber wn             cc
cc        rad - radiance at each output zenith and azimuth angle      cc
cc  trn_ray_0 - normal-incidence rayleigh-scattering transmission     cc
cc  trn_gas_0 - normal-incidence gas transmission                     cc
cc  trn_ray_0 - normal-incidence aerosol transmission                 cc
cc                                                                    cc
cccccccccccccccccccccccccc   r a d _ o u t   ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nfmx,ncd
      parameter (ncd=mxrad*mxrad+mxumu*mxphi+6, 
     -           nfmx = (mxlout+mxpd+1)*nsol+1)
c
      logical lsolar
c
      integer iuout,iutrn,iu_flux,nlyr,nz0,nza,ifrmout,
     -        irad,levout,nlout,nphi,numu,nzup,isptype,islit
      integer nze_1(mxlout),iscwt,nlev,nz
      integer nfl,ntrn,iquit_rad(nsol),iquit_trn,iquit_flx(nsol)
      integer nl,k,iurad,nze,naz,no,nn,iord
      save iquit_rad,iquit_trn,iquit_flx
c
      real p(kp),solflx
c
      double precision wn,wnmin,wnmax,width,dwn
c
      real dn_dir_flx(kp),dn_diff_flx(kp),up_diff_flx(kp)
      real var_1,trn_ray_0,trn_gas_0,trn_aer_0
      real wl,wnio
c
c****    output radiances
c
      real dir_s_flx(mxrad,nsol,2),up_s_flx(mxrad,nsol,2),
     -        dn_s_flx(mxrad,nsol,2),
     -        up_t_flx(mxrad,2),dn_t_flx(mxrad,2),
     -        up_flx(mxlout),dn_flx(mxlout),
     -        dir_flx(mxlout),rad(mxumu,mxphi,mxlout)
c
c****    scaled output radiance (units converted to iunits)
c
      real units,up_flx0(mxlout),dn_flx0(mxlout),
     -        dir_flx0(mxlout),rad0(mxumu,mxphi,mxlout)
c
      real p_ray0(mxumu),p_gas0(mxumu),p_aer0(mxumu),
     -        tau_ray_0,tau_gas_0,tau_aer_0,
     -        p_ray_0,p_gas_0,p_aer_0
c
c****   slit function variables
c
      real spect(ncd),points(nfmx)
c
c****    number of levels
c
      nlev = nlyr + 1
c
c****    index of solar zenith angle
c
       nz = nz0
c
c****   specify the file index for the first flux output file 
c       (required for slit function programs)
c
      nfl = nlout*nza + nz
c
c****    specify the file index for the transmission output file:
c
      ntrn = (nlout + 1)*nza + 1
c
c****   set io variables
c
      if(isptype .eq. 2) then
        if(wn .le. wnmin) then
          iquit_rad(nz) = 0
          iquit_trn = 0
          iquit_flx(nz) = 0
        else
          if(wn .ge. wnmax) then
            if(iquit_rad(nz) .eq. 0) iquit_rad(nz) = 1
            if(iquit_trn .eq. 0) iquit_trn = 1
            if(iquit_flx(nz) .eq. 0) iquit_flx(nz) = 1
          endif
        endif
      endif
c
c*****  define the wavelength
c
      wnio = real(wn)
      if(wn .ne. 0.0d0) then
        wl = real(1.0d4/wn)
      else
        wl = 1.0e9
      endif
c
c****    define the first output zenith angle for radiances 
c        (only upward radiances are saved at the top fo the atmosphere)
c
      do 1021 nl=1,nlout
          if(levout .eq. 1 .or. (levout .eq. 3 .and. nl .eq. 1)) then
             nze_1(nl) = nzup
          else
             nze_1(nl) = 1
          endif
1021  continue
c
      trn_ray_0 = exp(-tau_ray_0)
      trn_gas_0 = exp(-tau_gas_0)
      trn_aer_0 = exp(-tau_aer_0)
c
      if(lsolar) then
        var_1 = units*solflx
      else
        var_1 = tau_ray_0 + tau_gas_0 + tau_aer_0
      endif
c
      if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -         .or. irad .eq. 8) then
        if(wn .le. wnmin) then
            if(ifrmout .eq. 1) then
              write(iu_flux,'(9(1pe14.6))') (p(k),k=1,nlev)
            else
              write(iu_flux) wl,wnio,(p(k),k=1,nlev)
            endif
        endif
      endif
c
c****     find  total downward and upward diffuse fluxes at each level
c         and convert to W/m**2/cm-1 to iunits
c
      do 1041 k=1,nlev
          dn_dir_flx(k) = units*dir_s_flx(k,nz,2)
          dn_diff_flx(k) = units*(dn_s_flx(k,nz,2) + dn_t_flx(k,2))
          up_diff_flx(k) = units*(up_s_flx(k,nz,2) + up_t_flx(k,2))
1041  continue
c
c***     convert radiances and fluxes form W/m**2/cm-1 to iunits
c
      do 1141 k=1,nlout
          dir_flx0(k) = units*dir_flx(k)
          dn_flx0(k) = units*dn_flx(k)
          up_flx0(k) = units*up_flx(k)
          do 1121 naz=1,nphi
              do 1101 nze=nze_1(k),numu
                  rad0(nze,naz,k) = units*rad(nze,naz,k)
1101          continue
1121      continue
1141  continue
c
      if(isptype .eq. 1) then
c
c****     write monochromatic data at full resolution
c
        do 1221 k=1,nlout    
c
c****         compute the unit number
c    
            iurad = iuout + (nz-1)*nlout + k - 1
c
c****         if only top-of-atmosphere radiances are to be saved,
c             print out only upward radiances and fluxes.  
c             For output at other levels, print out both upward 
c             and downward fluxes and radiances
c
            if(levout .eq. 1 .or. (levout .eq. 3 .and. k .eq. 1)) then
c
c****          print out only upward streams at top of atmosphere
c
              if(ifrmout .eq. 1) then
c
c****               print a formatted output file
c
                write(iurad,'(9(1pe14.6))') 
     -                wl,wnio,var_1,up_flx0(k),
     -                ((rad0(nze,naz,k),nze=nze_1(k),numu),naz=1,nphi)
              else
c
c****               print an unformatted output file
c
                write(iurad) wl,wnio,var_1,up_flx0(k),
     -                ((rad0(nze,naz,k),nze=nze_1(k),numu),naz=1,nphi)
              endif
c
            else
c
c****           print out both upward and downward streams.
c
              if(ifrmout .eq. 1) then
c
c****               print a formatted output file
c
                write(iurad,'(6(1pe14.6))') wl,wnio,var_1,
     -               dir_flx0(k),dn_flx0(k),up_flx0(k)
                do 1201 naz=1,nphi
                    write(iurad,'(14x,9(1pe14.6))')
     -                   (rad0(nze,naz,k),nze=nze_1(k),numu)
1201            continue
              else
c
c****               print an unformatted output file
c
                write(iurad) wl,wnio,var_1,
     -               dir_flx0(k),dn_flx0(k),up_flx0(k),
     -               ((rad0(nze,naz,k),nze=nze_1(k),numu),naz=1,nphi)
              endif
            endif
c
1221    continue
c
c****       print other output radiance quantities:
c
        if(ifrmout .eq. 1) then
c
c****       print formatted output files
c
          if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -       .or. irad .eq. 8) then
c
c
c****         print level-dependent fluxes at the wavelength
c
            if(lsolar) then
              write(iu_flux,'(9(1pe14.6))')
     -              wl,wnio,(dn_dir_flx(k),k=1,nlev),
     -                    (dn_diff_flx(k),k=1,nlev),
     -                    (up_diff_flx(k),k=1,nlev)
            else
              write(iu_flux,'(9(1pe14.6))')
     -              wl,wnio,(dn_diff_flx(k),k=1,nlev),
     -                    (up_diff_flx(k),k=1,nlev)
            endif
          endif
c
          if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -       .or. irad .eq. 8) then   
c
c****         print the transmission to space and the p(tau=1)
c         
            if(nz .eq. 1) write(iutrn,'(9(1pe14.6))') wl,wnio,
     -        trn_ray_0,p_ray_0,(p_ray0(nze),nze=nzup,numu),
     -        trn_gas_0,p_gas_0,(p_gas0(nze),nze=nzup,numu),
     -        trn_aer_0,p_aer_0,(p_aer0(nze),nze=nzup,numu)
          endif
c
        else
c
          if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -       .or. irad .eq. 8) then
c
c****         print level-dependent fluxes at the wavelength
c
            if(lsolar) then
              write(iu_flux)
     -              wl,wnio,(dn_dir_flx(k),k=1,nlev),
     -                    (dn_diff_flx(k),k=1,nlev),
     -                    (up_diff_flx(k),k=1,nlev)
            else
              write(iu_flux)
     -              wl,wnio,(dn_diff_flx(k),k=1,nlev),
     -                    (up_diff_flx(k),k=1,nlev)
            endif
          endif
c
          if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -       .or. irad .eq. 8) then
c
c****         print the transmission to space and the p(tau=1)
c 
            if(nz .eq. 1)  write(iutrn) wl,wnio,
     -        trn_ray_0,p_ray_0,(p_ray0(nze),nze=nzup,numu),
     -        trn_gas_0,p_gas_0,(p_gas0(nze),nze=nzup,numu),
     -        trn_aer_0,p_aer_0,(p_aer0(nze),nze=nzup,numu)
          endif
        endif
c
      else
c
c****      call slit function program
c
        do 2021 k=1,nlout
c
c****         compute the unit number
c    
            iurad = iuout + (nz-1)*nlout + k - 1
            no = (nz-1)*nlout + k
            points(no) = points(no) + 1.0
c
c****         load variables into slit arrays
c
            nn = 1
            spect(nn) = var_1
            if(levout .eq. 2 .or. levout .ge. 4 .or. 
     -         (levout .eq. 3 .and. k .gt. 1)) then
              nn = nn + 1
              spect(nn) = dir_flx0(k)
              nn = nn + 1
              spect(nn) = dn_flx0(k)
            endif
            nn = nn + 1                
            spect(nn) = up_flx0(k)
            do 2011 naz=1,nphi
                do 2001 nze=nze_1(k),numu
                    nn = nn + 1
                    spect(nn) = rad0(nze,naz,k)
2001            continue
2011        continue
c
            iord = 2
            iscwt = 0
c
            call rad_slit(iurad,no,ifrmout,islit,iscwt,iord,nn,
     -                    wnmin,wnmax,dwn,width,
     -                    points,wn,spect,iquit_rad(nz))
2021    continue
c
        if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -     .or. irad .eq. 8) then
          if(nz .eq. 1) then
c
c*****         add pressure levels of optical depth unity and column-
c              integrated transmission of gases and aerosols.
c
            points(ntrn) =  points(ntrn) + 1.0
            nn = 1
            spect(nn) = trn_ray_0
            nn = nn + 1
            spect(nn) = p_ray_0
            do 3001 nze=nzup,numu
                nn = nn + 1
                spect(nn) = p_ray0(nze)
3001        continue
            nn = nn + 1
            spect(nn) = trn_gas_0
            nn = nn + 1
            spect(nn) = p_gas_0
            do 3021 nze=nzup,numu
                nn = nn + 1
                spect(nn) = p_gas0(nze)
3021        continue
            nn = nn + 1
            spect(nn) = trn_aer_0
            nn = nn + 1
            spect(nn) = p_aer_0
            do 3041 nze=nzup,numu
                nn = nn + 1
                spect(nn) = p_aer0(nze)
3041        continue
c
            iord = 2
            iscwt = 0
c
            call rad_slit(iutrn,ntrn,ifrmout,islit,iscwt,iord,nn,
     -                    wnmin,wnmax,dwn,width,
     -                    points,wn,spect,iquit_trn)
c
          endif
        endif
c
c****      call slit function program for level-dependent net fluxes
c
        if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -     .or. irad .eq. 8) then
c
c****        increment counting variables
c
          points(nfl) =  points(nfl) + 1.0
c
          if(lsolar) then
            nn = 3*nlev
            do 4021 k=1,nlev
                spect(k) = dn_dir_flx(k)
                spect(k+nlev)= dn_diff_flx(k)
                spect(k+2*nlev) = up_diff_flx(k)
4021        continue
          else
            nn = 2*nlev
            do 4041 k=1,nlev
                spect(k)= dn_diff_flx(k)
                spect(k+nlev) = up_diff_flx(k)
4041        continue
          endif
c
          iord = 2
          iscwt = 0
c
          call rad_slit(iu_flux,nfl,ifrmout,islit,iscwt,iord,nn,
     -                  wnmin,wnmax,dwn,width,
     -                  points,wn,spect,iquit_flx(nz))
c
        endif
c
      endif
c
      return
      end
