      subroutine smart_hdr
     -     (iuout,iutrn,iuflx,atmfile,aerfile,miefile,solfile,
     -     surfile,mixfile,gasfile,iunits,nmodes,ncomp,icomp,volmix,
     -     ts,au,wgtatm,wnmin,wnmax,isptype,islit,width,dwn,
     -     tauerr,pi0err,phferr,surferr,nstr,numu,umu,
     -     nphi,phi,nza,umu0,phi0,levout,nlout,pout,accur,
     -     lamber,isource,irad,ifrmout,radius,sgrav,
     -     igas,ngases,ngastyp)
c
ccccccccccccccccccccccccc  s m a r t _ h d r  cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This subroutine creates a header that documents the input       cc
cc    parameters used in the program smart                            cc
cc                                                                    cc
cc    i n p u t:                                                      cc
cc                                                                    cc
cc      iusol - unit number with solar radiances                      cc
cc     iusol1 - unit number for solar radiance scratch file           cc
cc     iuthrm - unit number for thermal radiances                     cc
cc      iuout - unit number of output radiance file                   cc
cc      iutrn - unit number of output transmission/pressure file      cc
cc    atmfile - name of input atmospheric structure file              cc
cc    aerfile - name of input aerosol vertical structure file         cc
cc    miefile - name of file with aerosol optical properties vs. wn   cc
cc    solfile - name of file with wn-dependent solar fluxes           cc
cc    surfile - name of file with wn-dependent surface optics         cc
cc    mixfile - name of file with gas mixing ratios                   cc
cc    gasfile - name fo file with gas absorption coeffiecients vs. wn cc
cc     nmodes - number of discrete aerosol partical modes             cc
cc      ncomp - number of rayleigh-scattering constituentes           cc
cc      icomp - index of each rayleigh scattering constituent         cc
cc     volmix - volume mixing ratio of each rayleigh scatterer        cc
cc         ts - sufrace temperature (K)                               cc
cc         au - distance to the sun (in AU's)                         cc
cc     tauerr - optical depth relative binning error (0. to ~0.8)     cc
cc     pi0err - co-single scattering albedo absolute binning error    cc
cc     phferr - asymmetry factor absolute binning error               cc
cc    surferr - surface optical property binning error                cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc       pout - output pressure level (bars)                          cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, a bi-directional reflection      cc
cc              function is required.                                 cc
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
cc    ifmrout - index of output file format (1) ascii, (2) binary     cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       nphi - number of output azimuth angles                       cc
cc       nlyr - number of computational model layers                  cc
cc        nza - number of solar zenith angles                         cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) Watts/m**2/sr/Angstrom                             cc
cc        phi - emission azimuth angles (degrees)                     cc
cc       umu0 - solar zenith angle cosines                            cc
cc       phi0 - solar azimuth angle cosines                           cc
cc      wnmin - minimum wavenumber of desired spectral window         cc
cc      wnmax - maximum wavenumber of desired spectral window         cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    header file with smart/dart header format                       cc
cc                                                                    cc
ccccccccccccccccccccccccc  s m a r t _ h d r  cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
      integer lnrec
      parameter (lnrec = 80)
c
      character*132 atmfile,mixfile(ngas),gasfile(ngtmax,ngas),
     -             aerfile(nmode),miefile(nmode),solfile,surfile
      character*20 val
      character*80 record
c
      logical lamber
c
      integer icomp(6),igas(ngas),ngases,ngastyp(ngas)
      integer iuout,iutrn,iuflx,iunits,nmodes,ncomp,isptype,islit
      integer nstr,numu,nphi,nza,levout,nlout,isource,irad,ifrmout
      integer nrec,n,ngt,i,nl,nz,iuf,lenrec,m,iuo
c
      double precision wnmin,wnmax,width,dwn
      real ts,au,wgtatm,accur
      real tauerr,pi0err,phferr,surferr
      real phi(mxphi),umu(mxumu),phi0(nsol),umu0(nsol),pout(mxlout)
      real sza0
c
      real volmix(6),radius,sgrav
c
      lenrec = lnrec
      nrec = 0
c
c****  input atmospheric structure file.
c
      call header(iuout,ifrmout,'atmfile',atmfile,
     -            ' Atmospheric structure file ',8,30,lenrec)
      nrec = nrec + 1
c
      do 1021 m=1,nmodes
          call header(iuout,ifrmout, 'aerfile',aerfile(m),
     -    ' aerosol vert. str. file',8,30,lenrec)
          nrec = nrec + 1
          call header(iuout,ifrmout, 'miefile',miefile(m),
     -    ' aerosol mie scat. file',8,30,lenrec)
          nrec = nrec + 1
1021  continue
c
c****    a t m o s p h e r i c    g a s e s 
c
      write(val,'(1x,i5)') ngases
      call header(iuout,ifrmout,'ngases',val,
     - ' number of absorbing gases ',8,30,lenrec)
      nrec = nrec + 1
c
      do 1081 n=1,ngases
          val = ' '
          write(val,'(1x,i5)') igas(n)
          call header(iuout,ifrmout,'igas',val,
     -    ' afgl gas code for gas',8,30,lenrec)
          do 1041 ngt =1,ngastyp(n)
              nrec = nrec + 1
              call header(iuout,ifrmout, 'gasfile',gasfile(ngt,n),
     -        ' Name of abs. coeff.',8,30,lenrec)
1041      continue
          nrec = nrec + 1
          call header(iuout,ifrmout, 'mixfile',mixfile(n),
     -    ' gas mixing ratio file',8,30,lenrec)
          nrec = nrec + 1
1081  continue
c
c****    r a y l e i g h    s c a t t e r i n g
c
      write(val,'(1x,i5)') ncomp        
      call header(iuout,ifrmout,'ncomp',val,
     - ' Number of Rayleigh scatterers',8,30,lenrec)
      nrec = nrec + 1
      do 1661 i=1,ncomp
          write(val,'(1x,i5)') icomp(i)      
          call header(iuout,ifrmout,'icomp',val,
     -     ' Index: (1) air (2) co2 (3) n2 (4) o2',8,30,lenrec)
          nrec = nrec + 1
          write(val,'(1x,1pe15.6)') volmix(i)
          call header(iuout,ifrmout,'volmix',val,
     -                ' Volume mixing ratio ',8,30,lenrec)
          nrec = nrec + 1
1661  continue
c
c****    s u r f a c e    c h a r a c t e r i s t i c s
c
      write(val,'(1x,1pe12.4)') ts
      call header(iuout,ifrmout,'ts',val, ' Surface temperature ',
     -            8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,l7)') lamber
      call header(iuout,ifrmout,'lamber',val,
     -' Use Lambert surface?',8,30,lenrec)
      nrec = nrec + 1
c
        call header(iuout,ifrmout,'surfile',surfile,
     -   ' surface optical property file ',8,30,lenrec)
      nrec = nrec + 1
c
      if(isource .eq. 1 .or. isource .eq. 3) then
c
c****     s o l a r    f l u x e s
c
        call header(iuout,ifrmout,'solfile',solfile,
     -   ' Name of solar flux file ',8,30,lenrec)
c
      endif
c
c****  read physical properties of planet and atmosphere.
c
      write(val,'(1x,1pe15.6)') au
      call header(iuout,ifrmout,'au',val,
     - ' Planets distance from sun (au)',8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,1pe15.6)') sgrav
      call header(iuout,ifrmout,'sgrav',val,
     -' surface gravity (m/s**2) ',8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,1pe15.6)') radius
      call header(iuout,ifrmout,'radius',val,
     -            ' Radius of the planet (km) ',8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,1pe15.6)') wgtatm
      call header(iuout,ifrmout,'wgtatm',val,
     - ' Atmospheric molecular weight (kg/kmole) ',8,30,lenrec)
      nrec = nrec + 1
c
c****    o u t p u t   s p e c t r a l   g r i d.
c
      write(val,'(1x,1pe15.6)') wnmin
      call header(iuout,ifrmout,'wnmin',val,
     - 'minimum output wavenumber',8,30,lenrec)
      nrec = nrec + 1
      write(val,'(1x,1pe15.6)') wnmax
      call header(iuout,ifrmout,
     - 'wnmax',val,'maximum output wavenumber',8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,i5)') isptype
      call header(iuout,ifrmout,'isptype',val,
     - ' Index of output spectrum type ',8,30,lenrec)
      nrec = nrec + 1
c
      if(isptype .eq. 2) then
c
        write(val,'(1x,i5)') islit        
        call header(iuout,ifrmout,'islit',val,
     - ' Type of spectral response function ',8,30,lenrec)
        nrec = nrec + 1
c
        write(val,'(1x,1pe15.6)') width  
        call header(iuout,ifrmout,'width',val,
     -       ' Half-width-at-half-max ',8,30,lenrec)
        nrec = nrec + 1
c
        write(val,'(1x,1pe15.6)') dwn 
        call header(iuout,ifrmout,'dwn',val,
     -   ' output sampling resolution (cm**-1)',8,30,lenrec)
        nrec = nrec + 1
      endif
c
      write(val,'(1x,i5)') iunits        
      call header(iuout,ifrmout,'iunits',val,
     - ' Index of output radiance units ',8,30,lenrec)
      nrec = nrec + 1
c
c****   set error limits for smt method
c
      write(val,'(1x,1pe15.6)') tauerr
      call header(iuout,ifrmout,'tauerr',val,
     - ' Fractional error for smt tau binning ',8,30,lenrec)
      nrec = nrec + 1
      write(val,'(1x,1pe15.6)') pi0err
      call header(iuout,ifrmout,'pi0err',val,
     - ' error for smt pi0 binning ',8,30,lenrec)
      nrec = nrec + 1
      write(val,'(1x,1pe15.6)') phferr
      call header(iuout,ifrmout,'phferr',val,
     - ' error for smt <cos> binning ',8,30,lenrec)
      nrec = nrec + 1
      write(val,'(1x,1pe15.6)') surferr
      call header(iuout,ifrmout,'surferr',val,
     - ' error for smt surface optical property binning ',8,30,lenrec)
      nrec = nrec + 1
c
c****    d i s c r e t e   o r d i n a n t   m e t h o d
c
      write(val,'(1x,i5)') nstr      
      call header(iuout,ifrmout,'nstr',val,
     - ' Number of streams for dom (> 2) ',8,30,lenrec)
      nrec = nrec + 1
c
c****   o u t p u t    r a d i a n c e    a n g l e s
c
      write(val,'(1x,i5)') irad 
      call header(iuout,ifrmout,'irad',val,
     - ' Index of output radiances/fluxes dist.',8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,i5)') numu  
      call header(iuout,ifrmout,'numu',val,
     -  ' # of arbitrary emission zenith angles ',8,30,lenrec)
      nrec = nrec + 1
c
      do 3281 n=1,numu
          write(val,'(1x,1pe15.6)') umu(n)
          call header(iuout,ifrmout,'umu',val,
     -         ' Emission zenith angle',8,30,lenrec)
          nrec = nrec + 1
3281  continue
c
      write(val,'(1x,i5)') nphi
      call header(iuout,ifrmout,'nphi',val,
     -' # of emission azimuth angles ',8,30,lenrec)
      nrec = nrec + 1
c
      do 3641 n=1,nphi
          write(val,'(1x,1pe15.6)') phi(n)
          call header(iuout,ifrmout,'phi',val,
     -                ' Emission azimuth angle',8,30,lenrec)
          nrec = nrec + 1
3641  continue
c
c****    o u t p u t    l e v e l s
c
      write(val,'(1x,i5)') levout
      call header(iuout,ifrmout,'levout',val,
     - ' output level index: ',8,30,lenrec)
      nrec = nrec + 1
c
      write(val,'(1x,i5)') nlout
      call header(iuout,ifrmout,'nlout',val,
     - ' number of output levels: ',8,30,lenrec)
      nrec = nrec + 1
c
      do 3801 nl=1,nlout
          write(val,'(1x,1pe15.6)') pout(nl)
          call header(iuout,ifrmout,'pout',val,
     -     ' output pressure level (bars): ',8,30,lenrec)
          nrec = nrec + 1
3801  continue
c
c****   s o u r c e    f u n c t i o n s
c
      write(val,'(1x,i5)') isource   
      call header(iuout,ifrmout,'isource',val,
     - ' Index of source function type ',8,30,lenrec)
      nrec = nrec + 1
c
      if(isource .eq. 1 .or. isource .eq. 3) then
c
        write(val,'(1x,1pe15.6)') accur
        call header(iuout,ifrmout,'accur',val,
     -   ' Convergence criteria for azimuthal series',8,30,lenrec)
      endif
c
c****   rewind output file and write header for optical depth file 
c       and each solar zenith angle.
c
      rewind(iuout)
4001  if(ifrmout .eq. 1) then
         read(iuout,'(1a)',end=4201) record
      else
         read(iuout,end=4201) record
      endif
c
c****   print out the headers for the files for the remainder of the
c       levels at the first zenith angle.
c
      do 4101 nl=2,nlout
          iuo = iuout + nl - 1
          if(ifrmout .eq. 1) then
            write(iuo,'(1a80)') record
          else
            write(iuo) record
          endif
4101  continue
c
c*****   if the number of zenith angles is greater than 1, write 
c        out headers for remaining output zenith angles and levels
c
      do 4141 nz=2,nza
          do 4121 nl=1,nlout
              iuo = iuout + (nz-1)*nlout + nl - 1
              if(ifrmout .eq. 1) then
                write(iuo,'(1a80)') record
              else
                write(iuo) record
              endif
4121      continue
4141  continue
c
      if(irad .eq. 3 .or. irad .eq. 4 .or. irad .eq. 7 
     -         .or. irad .eq. 8) then
        do 4161 nz=1,nza
            iuf = iuflx + nz - 1
            if(ifrmout .eq. 1) then
              write(iuf,'(1a80)') record
            else
              write(iuf) record
            endif
4161    continue
      endif
c
c*****   write out to optical depth file
c
      if(irad .eq. 2 .or. irad .eq. 4 .or. irad .eq. 6 
     -   .or. irad .eq. 8) then
        if(ifrmout .eq. 1) then
          write(iutrn,'(1a80)') record
        else
          write(iutrn) record
        endif
      endif
c
      go to 4001
c
4201  continue 
c
      do 5021 nz=1,nza
          do 5001 nl=1,nlout
              iuo = iuout + (nz-1)*nlout + nl - 1
c
              if(isource .eq. 1 .or. isource .eq. 3) then
c
                sza0 = 180.*acos(umu0(nz))/3.14159
c                write(*,*) 'nz,umu0(nz),sza0',nz,umu0(nz),sza0
                write(val,'(1x,1pe15.6)') sza0
                call header(iuo,ifrmout,'sza0',val,
     -            ' solar zenith angle (deg)',8,30,lenrec)
                write(val,'(1x,1pe15.6)') phi0(nz)
                call header(iuo,ifrmout,'phi0',val,
     -                  ' Solar azimuth angle (deg)',8,30,lenrec)
              endif
c
c****            e n d   o f    h e a d e r
c
              call header(iuo,ifrmout,'end',' 0 ',
     -           ' End of file',8,30,lenrec)
c
5001      continue
5021  continue



c
      return
      end
