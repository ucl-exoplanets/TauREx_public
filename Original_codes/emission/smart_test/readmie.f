      subroutine readmie(mode,iumie,iuaer0,miefile,iang,ioffmie,
     -                   icwlq,icqext,icqsca,icg1,icg2,
     -                   icf,ioffmom,nstr,
     -                   io_end,io_err,wnaer,wneof)
c
cccccccccccccccccccccccccc  r e a d m i e  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    *******************   p u r p o s e   **********************    cc
cc                                                                    cc
cc    this subroutine reads  the aerosol extinction, absorption       cc
cc    and scattering efficiencies and the scattering asymmetry        cc
cc    parameters and phase function moments at all wavelengths.       cc
cc    these quantities are normalized to their values at a reference  cc
cc    wavelength and then written out into a binary file that is      cc
cc    monochromatically increasing in wavenumber.                     cc
cc                                                                    cc
cc    *******************     i n p u t     **********************    cc
cc                                                                    cc
cc          ns0 - number of wavelengths where scattering properties   cc
cc                are specified                                       cc
cc       wl0(i) - ith wavelength at which optical properties are      cc
cc                specified                                           cc
cc      qext(i) - extinction efficiency for the ith wavelength        cc
cc      qsca(i) - scattering efficiency for the ith wavelength        cc
cc        g1(i) - asymmetry parameter for forward scattering for      cc
cc                the  ith wavelength                                 cc
cc        g2(i) - asymmetry parameter for back scattering for         cc
cc                the  ith wavelength                                 cc
cc         f(i) - forward scattering fraction for two-term            cc
cc                Henyey-Greenstein phase function.                   cc
cc                                                                    cc
cc                      - calling arguments -                         cc
cc                                                                    cc
cc                                                                    cc
cc    *******************    o u t p u t    **********************    cc
cc                                                                    cc
cc       wl0(i) - ith wavelength at which optical properties are      cc
cc                specified                                           cc
cc      qext(i) - extinction efficiency for the ith wavelength        cc
cc      qsca(i) - scattering efficiency for the ith wavelength        cc
cc        g1(i) - asymmetry parameter for forward scattering for      cc
cc                the  ith wavelength                                 cc
cc        g2(i) - asymmetry parameter for back scattering for         cc
cc                the  ith wavelength                                 cc
cc         f(i) - forward scattering fraction for two-term            cc
cc                Henyey-Greenstein phase function.                   cc
cc          ns0 - number of wavelengths where scattering properties   cc
cc                are specified                                       cc
cc                                                                    cc
cccccccccccccccccccccccccc  r e a d m i e  ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer npmx
      parameter (npmx=8192)
      character*(132) miefile
      character*(132) miefil,momfil
      character*1 fileout(132)
c
      integer ioffmie,icwlq,icqext,icqsca,icg1,icg2,icf
      integer mode,iumie,iuaer0,iang,nstr,io_end,io_err
      integer nmomstd(npmx),nmom0
      integer i,k,l,m,length,nlb,iuaer,icmax,ia,ns0,mom,nmom,ntb,ioffmom
c
      double precision wn0d(npmx),dwnmn,dwn0,wn0p,wl0p
      double precision wnaer,wneof(2)
      real qext(npmx),qsca(npmx),g1(npmx),g2(npmx),f(npmx),
     -     pmomstd(0:mxmom,npmx),vector(100),pmom0(0:1024)
      real qexts
c
c*****    a e r o s o l    p a r t i c l e    m o d e    l o o p
c
      m = mode
c
c****    open the input mie parameter file
c
      call charsp(miefile,fileout,length,132,nlb,ntb)
c
      miefil = ' '
      write(miefil,'(132a)') (fileout(i),i=1,length),'.mie'
      write(*,'(/,132a)') ' mie file: ',(fileout(i),i=1,length),'.mie'
c
      open(iumie,file=miefil,form='formatted',status='old')
c
c****   create a scatch unit for aerosol optical properties
c
      iuaer = iuaer0 + m
      close(iuaer)
      open(iuaer,form='unformatted',status='scratch')
c
c****    skip unnecessary descriptive information at top of file
c
      do 1201 l=1,ioffmie
          read(iumie,*,end=1241)
1201  continue
      go to 1301
1241  io_end = 1
c
c****    find the last column to read:
c
1301  icmax = 1
      if(icwlq .gt. icmax) icmax = icwlq
      if(icqext .gt. icmax) icmax = icqext
      if(icqsca .gt. icmax) icmax = icqsca
      if(icg1 .gt. icmax) icmax = icg1
      if(icg2 .gt. icmax) icmax = icg2
      if(icf .gt. icmax) icmax = icf
c
c****   read the wavelength, and the extinction, and scattering
c       cross-sections or efficiencies and the scattering
c       assymetry parameter
c
      dwnmn = 1.0d20
      k = 0
c
1401  read(iumie,*,err=1601,end=1621) (vector(ia),ia=1,icmax)
c
            k = k + 1
            wn0d(k) = 1.0d4/vector(icwlq)
            qext(k) = vector(icqext)
            qsca(k) = vector(icqsca)
            g1(k) = vector(icg1)
            if(icg2 .gt. 0) then
              g2(k) = vector(icg2)
              f(k) = vector(icf) 
            else
              g2(k) = 0.0
              f(k) = 1.0
            endif
c              
            dwn0 = abs(wnaer - wn0d(k))
            if(dwn0 .lt. dwnmn) then
              dwnmn = dwn0
              qexts = vector(icqext)
            endif
c
      go to 1401
c
c****   the following statement lets the program stumble over 
c       unnecessary descriptive lines at the top of the file
c
1601  if(k .eq. 0) go to 1401
      io_err = 1
c
c****   set number of spectral intervals, and eof markers
c
1621  ns0 = k
      if(io_err .eq. 0) io_end = 1
      if(wn0d(1) .lt. wn0d(k)) then
        wneof(1) = wn0d(1)
        wneof(2) = wn0d(k)
      else
        wneof(1) = wn0d(k)
        wneof(2) = wn0d(1)
      endif
c
c****   normalize extinction and scattering efficiencies by qexts.
c
      do 2001 k=1,ns0
          qext(k) =  qext(k)/qexts
          qsca(k) =  qsca(k)/qexts
2001  continue
c
c****    close efficiency factor file
c
      close(iumie)
c
c****   initialize phase function legendre coefficients
c
      do 2221 k=1,ns0
c
          pmomstd(0,k) = 1.0
          do 2201 mom=1,mxmom
              pmomstd(mom,k) = 0.0
2201      continue
2221  continue
c
      if(iang .eq. 1) then
c
c****     create name of file with legendre-polynomial coefficients
c         for the scattering phase function
c
        momfil = ' '
        write(momfil,'(132a)') (fileout(i),i=1,length),'.mom'
        write(*,'(132a)') ' pmom file: ',(fileout(i),i=1,length),'.mom'
c
c****     open file with legendre coefficients of phase function
c
        open(iumie,file=momfil,status='old',form='formatted')
c
c****     skip over unnecessary descriptive records at top of file
c
        do 3001 i=1,ioffmom
            read(iumie,*)
3001    continue
c
        do 3041 k=1,ns0
c
c****         read wavelength and number of phase function coefficients
c
            read(iumie,*,end=3061) wl0p,nmom0
c
            if(wl0p .ne. 0.0) wn0p = 1.0d4/wl0p
            if(abs(wn0p - wn0d(k)) .gt. 1.0d-4*wn0d(k)) then
               write(*,'(/,1a,/,1a,1pe14.6,1a,1pe14.6)') 
     -         ' Error reading phase function moments: ',
     -         ' wn0p =',wn0p,'  wn0d =',wn0d(k)
              stop
            endif
            nmomstd(k) = nmom0
            if(nmomstd(k) .gt. mxmom) then
              write(*,'(/,1a,i10,/,1a,i5,1a,1pe12.4)') 
     -        ' readmie warning: Number of phase function moments, ',
     -         nmomstd(k),' exceeds dimension bounds: mxmom=',mxmom,
     -         ' at wavenumber, ',wn0p 
               nmomstd(k)= mxmom
            endif
c
c****         read the legendre polynomial coeficients:
c
            read(iumie,*) (pmom0(mom),mom=0,nmom0)
            do 3021 mom=0,nmomstd(k)
                pmomstd(mom,k) = pmom0(mom)
3021        continue
c
3041    continue
c
c****         close phase function momemnt file
c
3061    close(iumie)
c
      else
c
c****     use a henyey-greenstein phase function
c
        do 3241 k=1,ns0
            if(iang .ne. 3) then
              do 3201 mom=1,nstr
                  pmomstd(mom,k) = g1(k)**mom
3201          continue
            else
c
c****           use a double Henyey-Greenstein phase function
c
              do 3221 mom=0,nstr
                  pmomstd(mom,k) = f(k)*g1(k)**mom + 
     -                             (1.0 - f(k))*g2(k)**mom
3221          continue
            endif
            nmomstd(k) = nstr
3241    continue
c         
      endif
c
c****   write these values to iuaer unit - in order of increasing 
c       wavenumber
c
      if(wn0d(1) .lt. wn0d(ns0)) then
        do 3641 k=1,ns0
            if(nmomstd(k) .gt. mxmom) then
              write(*,'(1a,i5,1a,1pe14.6,1a,i5)') 
     -         'Warning, dimension mxmom=',mxmom,
     -         ' is not adequate to resolve phase function at wn=',
     -         wn0d(k),' where nmomstd=',nmomstd(k)
              nmom = mxmom
            else
              nmom = nmomstd(k)
            endif 
            write(iuaer) nmom,wn0d(k),qext(k),qsca(k),g1(k),
     -                   (pmomstd(mom,k),mom=0,nmom) 
3641    continue
      else
        do 3661 k=ns0,1,-1
            if(nmomstd(k) .gt. mxmom) then
              write(*,'(1a,i5,1a,1pe14.6,1a,i5)') 
     -         'Warning, dimension mxmom=',mxmom,
     -         ' is not adequate to resolve phase function at wn=',
     -         wn0d(k),' where nmomstd=',nmomstd(k)
              nmom = mxmom
            else
              nmom = nmomstd(k)
            endif 
            write(iuaer) nmom,wn0d(k),qext(k),qsca(k),g1(k),
     -                   (pmomstd(mom,k),mom=0,nmom) 
3661    continue
      endif
c
      write(*,'(/,1a,i5)') ' Number of input aerosol wavelengths = ',ns0
c
      rewind(iuaer)
c
      return
      end
