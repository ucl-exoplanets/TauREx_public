      subroutine smart_in
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
ccccccccccccccccccccccccccc  s m a r t _ i n  cccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This subroutine reads list-directed input for the program smart cc
cc                                                                    cc
cc    i n p u t:                                                      cc
cc                                                                    cc
cc      iusol1 - unit number with solar radiances                     cc
cc     iuthrm - unit number for thermal radiances                     cc
cc      iuout - unit number for output radiance file                  cc
cc      iutrn - unit number for output transmission/pressure file     cc
cc      iuflx - unit number for output level-dependent fluxes         cc
cc     iuheat - unit number for output solar heating rates            cc
cc     iustat - unit number for output binning statistics             cc
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
cc    surfile - name of file with wn-dependent surface optical prop.  cc
cc    mixfile - name of file with gas mixing ratios                   cc
cc    gasfile - name of file with gas absorption coeffiecients vs. wn cc
cc    radfile - name output file for flux/radiance spectra            cc
cc   heatfile - name of output file with heating/cooling rates        cc
cc   statfile - name of output file with binning statistics           cc
cc    trnfile - name of output file with tau=1 trans/pressure         cc
cc    flxfile - name of output file with level-dependent flux spectra cc
cc     nstate - number of variable elements in the state vecror       cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc              0 - not a variable component of the state vector      cc
cc              1 - surface pressure                                  cc
cc              2 - surface/atmospheric temperature                   cc
cc              3 - gas absorption coeffient                          cc
cc              4 - cloud/aerosol optical depth                       cc
cc              5 - surface albedo                                    cc
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
cc       pout - output pressure level (bars)                          cc
cc    isource - index of source function type: (1) solar only         cc
cc              (2) thermal only, (3) both                            cc
cc       irad - index of output file type:                            cc
cc              1) wavelength-dependent fluxes and radiances at the   cc
cc                 computational azimuths and zenith angles at the    cc
cc                 specified output levels, and spectrally-integrated cc
cc                 fluxes and heating rates at each computational     cc
cc                 level.                                             cc
cc              2) wavelength-dependent fluxes, radiances, and        cc
cc                 transmission values at computational zenith angles cc
cc                 and specified output levels, and spectrally-       cc
cc                 integrated fluxes and heating rates at each        cc
cc                 computational level.                               cc
cc              3) wavelength-dependent fluxes and radiances at the   cc
cc                 computational azimuths and zenith angles at the    cc
cc                 specified output levels, wavelength-dependent,     cc
cc                 level-dependent pressure-weighted, flux            cc
cc                 divergences, and spectrally integrated fluxes      cc
cc                 and heating rates at each computational level.     cc
cc              4) wavelength-dependent fluxes, radiances,            cc
cc                 transmission values, wavelength-dependent,         cc
cc                 level-dependent pressure-weighted, flux            cc
cc                 divergences, and spectrally-integrated fluxes      cc
cc                 and heating rates at each computational level.     cc
cc              5) fluxes, radiances, and heating rates,              cc
cc                 at arbitrary azimuths and zenith angles,           cc
cc              6) fluxes, radiances and transmission functions       cc
cc                 at arbitrary zenith angles,                        cc
cc              7) fluxes, radiances, and contribution functions      cc
cc                 at arbitrary zenith angles,                        cc
cc              8) fluxes, radiances, transmission functions,and      cc
cc                 contribution functions at arbitrary zenith angles. cc
cc     ifmout - index of output file format (1) ascii, (2) binary,    cc
cc                                          (3) binary, no header     cc
cc       iref - bidirectional reflectance options                     cc
cc              0 - Lambert                                           cc
cc              1 - Hapke's BDR model                                 cc
cc              2 - Breon's BDR model; combination of Li + Roujean    cc
cc              3 - Roujean's BDR model                               cc
cc              4 - Cox and Munk glint model                          cc
cc       sza0 - solar zenith angles (degrees)                         cc
cc       umu0 - cosine of solar zenith angles                         cc
cc       phi0 - solar azimuth angles (degrees)                        cc
cc       nstr - number of gaussian zenith angles used in D/O code     cc
cc       nphi - number of output azimuth angles                       cc
cc       nlev - number of output levels                               cc
cc       nlyr - number of computational model layers                  cc
cc      nlout - number of levels where radiance spectra are output    cc
cc     levout - output level index (1) top of atmosphere,             cc
cc              (2) surface, (3) arbitrary level                      cc
cc        nza - number of solar zenith angles                         cc
cc        phi - emission azimuth angles (degrees)                     cc
cc        umu - emission zenith angle cosines                         cc
cc      accur - azimuth convergence accuracy for D/O routine          cc
cc     lamber - Include a lambertian surface? (Logical: T/F)          cc
cc              note: if lamber = F, a bi-directional reflection      cc
cc              function is required.                                 cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc     usrang - output radiances at user angles? (logical: T/F)       cc
cc     iunits - index of output radiance units:                       cc
cc              1) Watts/m**2/sr/cm**-1                               cc
cc              2) Watts/m**2/sr/micron                               cc
cc              3) Watts/m**2/sr/nanometer                            cc
cc              4) ergs/s/cm**2/sr/cm-1                               cc
cc              5) photons/s/m**2/sr/micron                           cc
cc      wnmin - minimum wavenumber of desired spectral window         cc
cc      wnmax - maximum wavenumber of desired spectral window         cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    same as input                                                   cc
cc                                                                    cc
cc    m o d i f i c a t i o n    h i s t o r y                        cc
cc                                                                    cc
cc    8/19/00: an aerosol optical depth scaling factor, sctau(nmode), cc
cc             was added to the input prompt for the optical depth    cc
cc             pressure scaling factor. In its current form, the code cc
cc             should work whether sctau is specified or not, as long cc
cc             as the input record includes a non-numeric character.  cc
cc    4/15/06  modified to support radiance and flux jacobian files   cc  
cc                                                                    cc
ccccccccccccccccccccccccccc  s m a r t _ i n  cccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      character*132 filein,atmfile,mixfile(ngas),gasfile(ngtmax,ngas),
     -             aerfile(nmode),miefile(nmode),solfile,surfile
      character*132 radfile(mxlout,nsol),heatfile,statfile,trnfile,
     -             flxfile(nsol)
      character*40 frmmix(ngas),frmatm,frms0,frmsur
      character*1 name(132)
      character*4 clev(mxlout)
      character*6 cgas(ngas)
      character*9 j_ext(nex)
c
      logical usrang,lamber,lplanck,lsolar
      logical lcsh(nmode)
c
      integer ifrmsur,ifrmatm
      integer icomp(6)
      integer iref,nref,icalb(4)
      integer igas(ngas),ngases,ngastyp(ngas)
      integer iumix,iuabc(ngtmax,ngas),
     -        igtrn(ngtmax,ngas),igpab(ngtmax,ngas),
     -        ifrmmix(ngas),ioffmix(ngas),izmix(ngas),imix(ngas),
     -        icpmix(ngas),icmix(ngas)
      integer ioffmie(nmode),icwlq(nmode),icqext(nmode),icqsca(nmode),
     -        icg1(nmode),icg2(nmode),icf(nmode),iofftau(nmode),
     -        iztau(nmode),icptau(nmode),ictau(nmode),
     -        iccsh(nmode),ioffmom(nmode),iang(nmode)
      integer iuatm,iskatm,icp,ict,iugas,nmodes,iuaer,iumie,
     -        iusur,ioffsur,iwnsur,icwnsur,iusol1,
     -        ifrms0,ioffs0,icwns0,ixs0,iuins0,icsol,ncomp,
     -        isptype,islit,nstr,numu,nphi,nza,levout,nlout,
     -        isource,irad,iuout,iunits,ifrmout,iuheat,iustat,
     -        iutrn,iuflx
      integer icmx,ic,ist,len,nlb,ntb,j,npabs,n,ngt,m,isurtyp,
     -        ir,i,nl,nze,k,nzup,nzdn,nza_1(mxlout)
c
c****    units for flux and radiance jacobians
c 
      integer nstate,istate(nex),iu_pd,
     -        iutpd(mxpd,nsol),iuspd(mxpd,nsol),
     -        iupdrad(mxumu,mxphi,mxpd,mxlout,nsol)
c
      double precision wnmin,wnmax,width,dwn
      double precision wnaer(nmode)
      real phi(mxphi),umu(mxumu),phi0(nsol),umu0(nsol)
      real sza0(nsol),volmix(6),scpmix(ngas),scmix(ngas),
     -     scptau(nmode),sctau(nmode),sccsh(nmode)
      real pout(mxlout),ts,au,radius,wgtatm
      real pd_frac(nex)
      real ws,phiw
c
      real wgtgas(ngas),scp,sgrav,scwalb,scsol,accur,scwns0,scalb(4)
      real tauerr,pi0err,phferr,surferr
      real pi,sctau0,ang
c
      data cgas /
     1 '   H2O','   CO2','    O3','   N2O','    CO','   CH4','    O2',
     2 '    NO','   SO2','   NO2','   NH3','  HNO3','    OH','    HF',
     3 '   HCl','   HBr','    HI','   ClO','   OCS','  H2CO','  HOCl',
     4 '    N2','   HCN',' CH3Cl','  H2O2','  C2H2','  C2H6','   PH3',
     5 '  COF2','   SF6','   H2S',' HCOOH','   HO2','     O','ClONO2',
     6 '   NO+','  HOBr','  C2H4', 'CH3OH','H2','He',9*'Other'/
c
      pi = acos(-1.0)
c
c****  set maximimum number of retries for input errors
c
      icmx = 10
c
      do 1001 n=1,nex
          istate(n) = 0
1001  continue
c
c****   initialize the number of values in the state structure:
c       The first 2 values are pressure and temperature
c
      ic = 0 
1101  ic = ic + 1
      if(ic .gt. icmx) stop
      nstate = 1
      istate(1) = 0
      write(*,'(/,1a)') ' Compute Jacobians for Pressure?'
      write(*,'(1a,3(/,1a))') ' Enter index of Jacobian Type: ',
     -      ' 0) No Jacobians',
     -      ' 1) Radiance Jacobians',
     -      ' 2) Flux Jacobians'
      read(*,*,err=1101) ist
c
      if(ist .ne. 0) then
        j_ext(1) = '.j_pres'
        if(ist .eq. 1) then
          istate(nstate) = 1
        else
          istate(nstate) = -1
        endif
c
        write(*,'(3(1x,1a,i5),2a)') 'ist =',ist,' nstate=',nstate,
     -           ' istate =',istate(nstate),' j_ext =',j_ext(nstate)
c
        write(*,'(/,2a)')
     -   ' Enter the fractional change in pressure (0.0-1.0): '
        read(*,*) pd_frac(nstate)
        write(*,'(1x,1a,1pe13.5)') 'pd_frac =',pd_frac(nstate)
        else
          write(*,'(1x,1a,i5)') 'ist =',ist        
      endif
c
      nstate = 2
      istate(nstate) = 0
      pd_frac(nstate) = 0.0
1111  write(*,'(/,2a,/,1a)') 
     -' Compute Jacobians for Temperature?'
      write(*,'(1a,3(/,1a))') ' Enter index of Jacobian Type: ',
     -      ' 0) No Jacobians',
     -      ' 1) Radiance Jacobians',
     -      ' 2) Flux Jacobians'
      read(*,*,err=1111) ist
c       
      if(ist .ne. 0) then
        j_ext(nstate) = '.j_temp'
        if(ist .eq. 1) then
          istate(nstate) = 2
        else 
          istate(nstate) = -2
        endif
        write(*,'(3(1x,1a,i5),2a)') 'ist =',ist,' nstate=',nstate,
     -           ' istate =',istate(nstate),' j_ext =',j_ext(nstate)
        write(*,'(/,2a)')
     -   ' Enter the fractional change in Temperature (0.0-1.0): '
        read(*,*) pd_frac(nstate)
        write(*,'(1x,1a,1pe13.5)') 'pd_frac =',pd_frac(nstate)        
        else
          write(*,'(1x,1a,i5)') 'ist =',ist        
      endif         
c
c****  enter input file names and other data for this run.
c
      icmx = 5
      ic = 0
1121  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(2(/,1a))') 
     - ' Enter atmospheric structure file format index (1 - 3):  ',
     - '  1) formatted 2) unformatted 3) list directed '
      read(*,*,err=1121) ifrmatm
      write(*,'(1x,1a,i5)') 'ifrmatm =',ifrmatm
c
      if(ifrmatm .eq. 1) then
        write(*,'(2(/,1a))')
     -    ' Enter format for reading pressure and temperature',
     -    ' [enclosed in parenthesis, ie; (2f3.5) ]: '
        read(*,'(1a)') frmatm
        write(*,'(1x,2a)') 'frmatm = ',frmatm
      else
        frmatm = ' '
      endif
c
      ic = 0
1132  ic = ic + 1
      if(ic .gt. icmx) stop
      atmfile = ' '
      write(*,'(/,1a)') 
     -' Enter name of the Atmospheric Structure file: '
      read(*,'(1a)') atmfile
c
      call charsp(atmfile, name, len,132,nlb,ntb) 
c
      filein = ' '
      write(filein,'(132a)') (name(j),j=1,len)
      write(*,'(/,1a,i5,2a)') ' unit: ',iuatm,' atmfile:',filein(1:len)
c
      if(ifrmatm .ne. 2) then
        open(iuatm,file=filein(1:len),form='formatted',status='old',
     -        err=1132)
      else
        open(iuatm,file=filein(1:len),form='unformatted',status='old',
     -        err=1132)
      endif
c
      close(iuatm)
c
      ic = 0
1141  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') 
     - ' Enter number of lines to skip at top of this file:'
      read(*,*,err=1141) iskatm
      write(*,'(1x,1a,i5)') 'iskatm =',iskatm
      ic = 0
1161  ic = ic + 1
      if(ic .gt. icmx) stop
        write(*,'(/,1a)') 
     - ' Enter columns containing p and T:'
      read(*,*,err=1161) icp,ict
      write(*,'(1x,2(1a,i5))') ' icp =',icp,'  ict =',ict 
      ic = 0 
1181  ic = ic + 1
      if(ic .gt. icmx) stop
        write(*,'(/,1a)') 
     - ' Enter scaling factor to convert pressure to Pascals:'
      read(*,*,err=1181) scp
      write(*,'(1x,1a,1pe13.5)') 'scp =',scp        
c
      ic = 0
1201  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter the surface temperature, ts (K): '
      read(*,*,err=1201) ts
      write(*,'(1x,1a,1pe13.5)') 'ts = ',ts
c
c****                  a b s o r b i n g    g a s e s
c
      ic = 0
1301  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(4(/,1a))')
     - ' Enter the number of different absorbing gases: '
      read(*,*,err=1301) ngases
      write(*,'(1x,1a,i5)') 'ngases =',ngases 
      if(ngases .gt. ngas) then
        write(*,*) 'Number of gases exceeds dimension bound: ngas=',ngas
        stop
      endif
c
      npabs = 0
c
      do 1841 n=1,ngases
c
c****       increment state structure counter
c
          nstate = nstate + 1
c
          ic = 0
1401      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,i5)') 
     -    ' Enter the HITRAN gas code for gas # ',n
          read(*,*,err=1401) igas(n)
c
          call charsp(cgas(igas(n)),name,len,6,nlb,ntb) 
c
          write(*,'(1x,1a,i5,10a)') 
     -         'igas =',igas(n),' name = ',(name(j),j=1,len)
          ic = 0
1421      ic = ic + 1
          if(ic .gt. icmx) stop
          istate(nstate) = 0
          write(*,'(/,10a)') 
     -    ' Compute Jacobians for ',(name(j),j=1,len)
          write(*,'(1a,3(/,1a))') ' Enter index of Jacobian Type: ',
     -      ' 0) No Jacobians',
     -      ' 1) Radiance Jacobians',
     -      ' 2) Flux Jacobians'
          read(*,*,err=1421) ist
c
          if(ist .ne. 0) then
            write(j_ext(nstate),'(1a3,6a)') '.j_',(name(j),j=1,len)
            if(ist .eq. 1) then
              istate(nstate) = 3
            else 
              istate(nstate) = -3
            endif
            write(*,'(3(1x,1a,i5),2a)') 'ist =',ist,' nstate=',nstate,
     -           ' istate =',istate(nstate),' j_ext =',j_ext(nstate)
c
            write(*,'(/,2a)')
     -       ' Enter fractional change in the mixing ratio (0.0-1.0): '
            read(*,*) pd_frac(nstate)
            write(*,'(1x,1a,1pe13.5)') 'pd_frac =',pd_frac(nstate)
          else
            write(*,'(1x,1a,i5)') 'ist =',ist
          endif         
c
          ic = 0
1432      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,10a)') 
     -    ' Enter  number of absorption coefficient types for ',
     -     (name(j),j=1,len)
          read(*,*,err=1432) ngastyp(n)
          if(ngastyp(n) .gt. ngtmax) then
            write(*,*) ' Error: number of absorption coeffienct types ',
     -                'exceeds maximum dimension bound: ngtmax=',ngtmax
            go to 1432
          endif 
          write(*,'(1x,1a,i5)') 'igas =',ngastyp(n)
c
          do 1471 ngt =1,ngastyp(n)
c
              ic = 0
1441          ic = ic + 1
              if(ic .gt. icmx) stop
              write(*,'(/,2a,/,3(/,1a))') 
     -        ' Enter index of transition type for ',cgas(igas(n)),
     -        '  (1) vibration-rotation transitions (line absorbers), ',
     -        '  (2) pressure-induced transitions (density-squared), ',
     -        '  (3) cross-sections per molecule (UV and VIS continuum)'
              read(*,*,err=1441) igtrn(ngt,n)
              write(*,'(1x,1a,i5)') 'igtrn =',igtrn(ngt,n)
c
c****           create a gas absorption coefficient input unit number.
c
              iuabc(ngt,n) = iugas + (n - 1)*ngtmax + ngt - 1
c
              ic = 0
1462          ic = ic + 1
              if(ic .gt. icmx) stop
              write(*,'(/,10a)')
     -         ' Enter name of absorption coefficient file for ',
     -           (name(j),j=1,len)
              read(*,'(1a)') gasfile(ngt,n)
c
              write(*,'(/,1a,i5)') ' opening gas unit:iuabc(n)=',
     -                              iuabc(ngt,n)
              if(igtrn(ngt,n) .eq. 1) then
                open(iuabc(ngt,n),file=gasfile(ngt,n),
     -               form='unformatted',status='old',err=1462)
              else
                open(iuabc(ngt,n),file=gasfile(ngt,n),form='formatted',
     -               status='old',err=1462)
c
                if(igtrn(ngt,n) .eq. 2) then
                  npabs = npabs + 1
                  igpab(ngt,n) = npabs
                endif
              endif
c
              write(*,'(/,2a,/,2x,1a)') 
     -       ' Absorption coefficient file for ',cgas(igas(n)),
     -         gasfile(ngt,n)
c
1471      continue
c
          ic = 0
1481      ic = ic + 1
          if(ic .gt. icmx) stop
            write(*,'(2(/,1a))') 
     -     ' Enter gas mixing ratio file format index (1 - 3):  ',
     -     '  1) formatted 2) unformatted 3) list directed '
          read(*,*,err=1481) ifrmmix(n)
          write(*,'(1x,1a,i5)') 'ifrmmix =',ifrmmix(n)
c
          if(ifrmmix(n) .eq. 1) then
            write(*,'(2(/,1a))')
     -      ' Enter format for reading pressure and rmix',
     -      ' [enclosed in parenthesis, ie; (2f3.5) ]: '
            read(*,'(1a)') frmmix(n)
            write(*,'(1x,2a)') 'frmmix = ',frmmix(n)
          else
            frmmix(n) = ' '
          endif
c
          ic = 0
1602      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,2a)') 
     -     ' Enter name of gas mixing ratio file for ',cgas(igas(n))
          read(*,'(1a)') mixfile(n)
c
          write(*,'(/,1a,i5)') 
     -      ' Opening gas mixing ratio unit: iumix=',iumix
          if(ifrmmix(n) .ne. 2) then
            open(iumix,file=mixfile(n),form='formatted',status='old',
     -           err=1602)
          else
            open(iumix,file=mixfile(n),form='unformatted',status='old',
     -           err=1602)
          endif
c
          write(*,'(1x,10a)') (name(j),j=1,len),
     -        ' mixfile(n) = ',mixfile(n)
c
          close(iumix)
c
          ic = 0
1621      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') 
     -     ' Enter number of lines to skip at top of file:'
          read(*,*,err=1621) ioffmix(n)
          write(*,'(1x,1a,i5)') 'ioffmix =',ioffmix(n)
c
          ic = 0
1641      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,/,1a)') 
     -    ' Enter type of vertical coordinate: ',
     -    '  1) pressure,   2) altitude (km)'
          read(*,*,err=1641) izmix(n)
          write(*,'(1x,1a,i5)') 'izmix =',izmix(n)
          ic = 0
1661      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,10a)') 
     -    ' Enter number of column containing z and rmix for ',
     -    (name(j),j=1,len)
          read(*,*,err=1661) icpmix(n),icmix(n)
          write(*,'(1x,1a,i5)') 'icpmix =',icpmix(n)
          write(*,'(1x,1a,i5)') 'icmix =',icmix(n)
          ic = 0
1681      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(2a,/,2(/,1a))') 
     -     ' Enter type of mixing ratios for ',cgas(igas(n)),
     -     '  1) volume mixing ratio,',
     -     '  2) mass mixing ratio'
          read(*,*,err=1681) imix(n)
          write(*,'(1x,1a,i5)') 'imix =',imix(n)
          ic = 0
1801      ic = ic + 1
          if(ic .gt. icmx) stop
          if(izmix(n) .eq. 1) then
            write(*,'(/,1a,/,1a)')
     -     ' Enter factors to convert pressure to Pascals ',
     -     ' and mixing ratios to range 0-1: '
          else 
            write(*,'(/,1a,/,1a)')
     -     ' Enter multiplicative factors to convert altitudes to km ',
     -     ' and mixing ratios to range 0-1: '
          endif
          read(*,*,err=1801) scpmix(n),scmix(n)
          write(*,'(1x,1a,1pe13.6)') 'scpmix =',scpmix(n)
          write(*,'(1x,1a,1pe13.6)') 'scmix  =',scmix(n)
          if(imix(n) .eq. 1) then
            ic = 0
1821        ic = ic + 1
            if(ic .gt. icmx) stop
            write(*,'(/,1a)') 
     -       ' Enter molecular weight of gas (kg/Kmole): '
             read(*,*,err=1821) wgtgas(n)
             write(*,'(1x,1a,1pe13.6)') 'wgtgas =',wgtgas(n)
          endif
c
1841  continue
c
c****       a e r o s o l    o p t i c a l    p r o p e r t i e s
c
c****   enter the number of particle modes.
c
      ic = 0
2001  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,/,1a)') ' Enter number of aerosol particle modes: '
      read(*,*,err=2001) nmodes
      write(*,'(1x,1a,i5)') 'nmodes =',nmodes
      if(nmodes .gt. nmode) then
        write(*,'(1a,i5)') 
     -  ' Error: Number of modes exceeds dimension bound: nmode=',nmode
        stop
      endif
c
c****   initialize i/o variables for each mode
c
      do 2101 m=1,nmode
          miefile(m) = ' '
          ioffmie(m) = 0
          ioffmom(m) = 0
          icwlq(m) = 0
          icqext(m) = 0
          icqsca(m) = 0
          icg1(m) = 0
          icg2(m) = 0
          icf(m) = 0
2101  continue
c
c****   read in file names for aerosol optical properties
c
      do 2861 m=1,nmodes
          nstate = nstate + 1
          ic = 0
2201      ic = ic + 1
          if(ic .gt. icmx) stop
            istate(nstate) = 0
            write(*,'(/,1a,i5)') 
     -      ' Compute Jacobians for aerosol mode: ',m
            write(*,'(1a,3(/,1a))') ' Enter index of Jacobian Type: ',
     -      ' 0) No Jacobians',
     -      ' 1) Radiance Jacobians',
     -      ' 2) Flux Jacobians'
            read(*,*,err=2201) ist
c
            if(ist .ne. 0) then
c
              write(j_ext(nstate),'(1a7,i2.2)') '.j_aer',m
              if(ist .eq. 1) then
                istate(nstate) = 4
              else 
                istate(nstate) = -4
              endif
              write(*,'(3(1x,1a,i5),2a)') 'ist =',ist,' nstate=',nstate,
     -           ' istate =',istate(nstate),' j_ext =',j_ext(nstate)
              write(*,'(/,2a)')
     -         ' Enter fractional change in optical depth (0.0-1.0): '
              read(*,*) pd_frac(nstate)
              write(*,'(1x,1a,1pe13.5)') 'pd_frac =',pd_frac(nstate)
            else
              write(*,'(1x,1a,i5)') 'ist =',ist 
            endif
c
          ic = 0
2221      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,i5)') 
     -    ' Enter prefix of the mie scattering files for mode # ',m
          read(*,'(1a)') miefile(m)
c
          call charsp(miefile(m), name, len,132,nlb,ntb)
c
          filein = ' '
          write(filein,'(132a)') (name(j),j=1,len),'.mie'
          write(*,'(/,2a)') ' miefile = ',filein
c
          write(*,'(1a,i5)') ' opening unit iumie=',iumie
          open(iumie,file=filein,form='formatted',
     -         status='old',err=2221)
          close(iumie)
c
          ic = 0
2241      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') 
     -    ' Enter number of records to skip above mie file: '
          read(*,*,err=2241) ioffmie(m)
          write(*,'(1x,1a,i5)') 'ioffmie =',ioffmie(m)
          ic = 0
2261      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') 
     -    ' Enter column numbers for wl, qext, qsca, and g1 : '
          read(*,*,err=2261) icwlq(m),icqext(m),icqsca(m),icg1(m)
          write(*,'(1x,1a,4i5)') 
     -         'icwlq, icqext, icqsca, icg1 =',icwlq(m),icqext(m),
     -                           icqsca(m),icg1(m)
c
          ic = 0
2281      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') ' Choose phase function type '
          write(*,'(3(/1a))') 
     -    ' 1) full phase function',
     -    ' 2) Henyey Greenstein,  ',
     -    ' 3) Double Henyey-Greenstein,  '
          read(*,*,err=2281) iang(m)
          write(*,'(1x,1a,i5)') 'iang =',iang(m)
          ioffmom(m) = 0
          if(iang(m) .eq. 1) then
            ic = 0
2401        ic = ic + 1
            if(ic .gt. icmx) stop
            write(*,'(/,1a)') 
     -      ' Enter number of records to skip above phase fcn file: '
              read(*,*,err=2401) ioffmom(m)
            write(*,'(1x,1a,i5)') 'ioffmom =',ioffmom(m)
          else
            if(iang(m) .eq. 3) then
              ic = 0
2421          ic = ic + 1
              if(ic .gt. icmx) stop
              write(*,'(/,1a)') 
     -        ' Enter column numbers for g2 and f (0-1): '
              read(*,*,err=2421) icg2(m),icf(m)
              write(*,'(1x,1a,2i5)') 
     -         'icg2, icf =',icg2(m),icf(m)
            else
              icg2(m) = 0
              icf(m) = 0
            endif
          endif
c  
c****       a e r o s o l    v e r t i c a l    s t r u c t u r e
c
          ic = 0
2461      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,i5)')
     -     ' Enter name of file with optical depths for mode: ',m
          read(*,'(1a)') aerfile(m)
          write(*,'(1x,2a)') 'aerfile =',aerfile(m)
c
          write(*,'(/,1a,i5)') ' opening unit iuaer=',iuaer
          open(iuaer,file=aerfile(m),form='formatted',
     -          status='old',err=2461)
          close(iuaer)
c
c****       a e r o s o l   o p t i c a l    d e p t h s
c 
          ic = 0
2601      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,2a)') 
     -      ' Enter standard wavenumber where optical',
     -      ' depths are given: '
          read(*,*,err=2601) wnaer(m)
          write(*,'(1x,1a,1pe13.5)') 'wnaer =',wnaer(m)
          write(*,'(/,1a)') 
     -    ' Enter number of records to skip above optical depths: '
          read(*,*) iofftau(m)
          write(*,'(1x,1a,i5)') 'iofftau =',iofftau(m)
          ic = 0
2621      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,/,1a)') 
     -     ' Enter type of vertical coordinate: ',
     -     '  1) pressure,   2) altitude (km)'
          read(*,*,err=2621) iztau(m)
          write(*,'(1x,1a,i5)') 'iztau =',iztau(m)
          ic = 0
2641      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,/,1a)')
     -    ' Enter column numbers with vertical coordinate ',
     -    ' and differential optical depth:'
          read(*,*,err=2641) icptau(m),ictau(m)
          write(*,'(1x,1a,2i5)') 'icptau, ictau =',icptau(m),ictau(m)
          ic = 0
          ic = ic + 1
          if(ic .gt. icmx) stop
          if(iztau(m) .eq. 1) then
            write(*,'(/,1a,/,1a)')
     -     ' Enter multiplicative factors to convert',
     -     ' pressure to Pascals and to scale optical depth: '
          else 
            write(*,'(/,1a)')
     -     ' Enter multiplicative factors to convert',
     -     ' altitudes to km and to scale optical depth: '
          endif
          scptau(m) = 1.0
          sctau(m) = 1.0
          read(*,*,err=2661) scptau(m), sctau0
          sctau(m) = sctau0
2661      write(*,'(1x,1a,1pe13.5)') 'scptau =',scptau(m),
     -                               'sctau =',sctau(m)
c
c****       a e r o s o l    s c a l e    h e i g h t s 
c
          ic = 0
2801      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,2a)') 
     -     ' Are aerosol scale heights specified for',
     -     ' each layer (true/false)? '
          read(*,*,err=2801) lcsh(m)
          write(*,'(1x,1a,l6)') ' lcsh =',lcsh(m)
          if(lcsh(m)) then
            ic = 0
2821        ic = ic + 1
            if(ic .gt. icmx) stop
            write(*,'(/,1a)')
     -      ' Enter column numbers with aerosol scale height:'
            read(*,*,err=2821) iccsh(m)
            write(*,'(1x,1a,i5)') 'iccsh =',iccsh(m)
            ic = 0
2841        ic = ic + 1
            if(ic .gt. icmx) stop
            write(*,'(/,1a)')
     -     ' Enter factor to convert scale heights to km: '
            read(*,*,err=2841) sccsh(m)
            write(*,'(1x,1a,1pe13.5)') 'sccsh =',sccsh(m)
          endif
c
2861  continue
c
c****       s u r f a c e    o p t i c a l    p r o p e r t i e s 
c
      ic = 0
3001  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a,/,1a)') 
     - ' Enter the index of the surface BRDF type: ',
     - ' (0) Lambertian, (1) Non-Lambertian: ' 
      read(*,*,err=3001) isurtyp
      if(isurtyp .eq. 0) then
        lamber = .true.
      else
        lamber = .false.
        write(*,'(/,1a,4(/,1a))') 
     -   ' Enter the index of the surface BRDF mode: ',
     -   '  1) Hapke BDR model',
     -   '  2) Breon BDR model; combination of Li + Roujean',
     -   '  3) Roujean BDR model',
     -   '  4) Cox and Munk glint model'
        read(*,*,err=3001) iref
        write(*,'(1x,1a,i5)') 'iref =',iref
      endif
      write(*,'(1x,1a,l6)') 'lamber =',lamber
c
      nstate = nstate + 1
      ic = 0
3011  ic = ic + 1
      if(ic .gt. icmx) stop
      istate(nstate) = 0
      write(*,'(/,2a,/,1a)') 
     -  ' Compute Jacobians for Surface Reflectance?'
      write(*,'(1a,3(/,1a))') ' Enter index of Jacobian Type: ',
     -      ' 0) No Jacobians',
     -      ' 1) Radiance Jacobians',
     -      ' 2) Flux Jacobians'
      read(*,*,err=3011) ist
c
      if(ist .ne. 0) then
        write(j_ext(nstate),'(1a7)') '.j_surf'
        if(ist .eq. 1) then
          istate(nstate) = 5
        else 
          istate(nstate) = -5
        endif
        write(*,'(3(1x,1a,i5),2a)') 'ist =',ist,' nstate=',nstate,
     -           ' istate =',istate(nstate),' j_ext =',j_ext(nstate)
c
3021    write(*,'(/,2a)')
     -  ' Enter the fractional change in reflectance (0.0-1.0): '
        read(*,*,err=3021) pd_frac(nstate)
        write(*,'(1x,1a,1pe13.5)') 'pd_frac =',pd_frac(nstate)        
        else
          write(*,'(1x,1a,i5)') 'ist =',ist        
      endif         
c
      ic = 0
3031  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a,/,1a)') 
     -   ' Enter the surface reflectance file format index (1 - 3): ',
     -   '  1) formatted 2) unformatted 3) list directed '
      read(*,*,err=3031) ifrmsur
      write(*,'(1x,1a,i5)') 'ifrmsur =',ifrmsur
c
      if(ifrmsur .eq. 1) then
        write(*,'(2(/,1a))')
     -    ' Enter format for reading surface properties',
     -    ' [enclosed in parenthesis, ie; (2f3.5) ]: '
        read(*,'(1a)') frmsur
        write(*,'(1x,2a)') 'frmsur =',frmsur
      else
        frmsur = ' '
      endif
c
      ic = 0
3042  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter name surface optical property file: '
      read(*,'(1a)') surfile
c
      write(*,'(/,1a,i5)') ' opening unit iusur=',iusur
      if(ifrmsur .ne. 2) then
        open(iusur,file = surfile,form='formatted',status='old',
     -      err=3042)
      else
        open(iusur,file = surfile,form='unformatted',status='old',
     -      err=3042)
      endif
c
      close(iusur)
      write(*,'(2a)') ' surfile =',surfile
c
      ic = 0
3061  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') 
     - ' Enter number of lines to skip at top of this file:'
      read(*,*,err=3061) ioffsur
      write(*,'(1x,1a,i5)') 'ioffsur =',ioffsur
c
c****    initialize the number of surface reflectance properties 
c        at each wavelength used in the subroutine bdrf
c             iref= 0, nref = 1: Lambert albedo
c             iref= 1, nref = 1: Hapke : single scatter albedo, W, and 
c                                angular width factor HH 
c             iref= 2, nref = 3: - Breon's BDR model: k0, k1, k2
c             iref= 3, nref = 3: - Roujean's BDR model: k0, k1, k2
c             iref= 4, nref = 2: - Cox and Munk glint model: n, k
c
      if(iref .eq. 0) then
        nref = 1
      else
        if(iref .eq. 1) then
          nref = 2
        else
          if(iref .eq. 2 .or. iref .eq. 3) then
            nref = 3
          else
            nref = 2
          endif
        endif
      endif
      ic = 0
3081  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,2a)') 
     - ' Enter columns with spectral quantity and each',
     - ' of the surface optical properties:'
      read(*,*,err=3081) icwnsur,(icalb(ir),ir=1,nref)
      write(*,'(1x,1a,5i5)') 'icwnsur, icalb =',
     -                        icwnsur,(icalb(ir),ir=1,nref)
c
      ic = 0
3221  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(2(/,1a))') 
     - ' Enter type of spectral quantity:',
     - '  1) wavelength,  2) wavenumber'
      read(*,*,err=3221) iwnsur
      write(*,'(1x,1a,i5)') 'iwnsur =',iwnsur
c
      if(iwnsur .eq. 1) then
        ic = 0
3241    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a)') 
     -  ' Enter factor needed to convert wavelengths to microns: '
        read(*,*,err=3241) scwalb
        write(*,'(1x,1a,1pe13.5)') 'scwalb =',scwalb
      else
        scwalb = 1.0
      endif
c
      ic = 0
3261  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') 
     - ' Enter factor needed to scale surface optical properties: '
      read(*,*,err=3261) (scalb(ir),ir=1,nref)
      write(*,'(1x,1a,4(1pe13.5))') 'scalb =',(scalb(ir),ir=1,nref)
c
      if(iref .eq. 4) then
c
c****    set wind speed and direction for Cox/Munk model
c
3301    write(*,'(/,1x,1a)') 
     -           'Enter the wind speed (m/s) and azimuth (deg): '
        read(*,*,err=3301) ws,phiw
        write(*,'(1x,1a,1pe12.4,1a,1pe12.4,1a)') 
     -        'Wind speed =',ws,' m/s   Wind Azimuth =',phiw,' deg'
      else
        ws = 0.0
        phiw = 0.0
      endif
c
c****  read physical properties of planet and atmosphere.
c
      ic = 0
3401  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter the planets distance from sun (au):'
      read(*,*,err=3401) au
      write(*,'(1x,1a,1pe13.5)') 'au =',au
c
      ic = 0
3441  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)')
     - ' Enter gravitational acceleration at surface (m/s**2): '
      read(*,*,err=3441) sgrav
      write(*,'(1x,1a,1pe13.5)') 'sgrav =',sgrav
c
      ic = 0
3461  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter the radius of the planet (km): '
      read(*,*,err=3461) radius
      write(*,'(1x,1a,1pe13.5)') 'radius =',radius
c
      ic = 0
3481  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)')
     - ' Enter the atmospheric mean molecular weight (kg/kmole): '
      read(*,*,err=3481) wgtatm
      write(*,'(1x,1a,1pe13.5)') 'wgtatm =',wgtatm
c
      ic = 0
3601  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)')
     - ' Enter the number of major atmospheric components:'
      read(*,*,err=3601) ncomp
      write(*,'(1x,1a,i5)') 'ncomp =',ncomp
      if(ncomp .gt. ngas) then
        write(*,*) 'ncomp exceeds dimension bound, ngas=',ngas
        stop
      endif
c
      do 3661 i=1,ncomp
          write(*,'(/,1a)')
     -     ' component types: (1) air  (2) co2  (3) n2  (4) o2'
          ic = 0
3621      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,i5)') ' Enter the index of component # ',i
          read(*,*,err=3621) icomp(i)
          write(*,'(1x,1a,i5)') 'icomp =',icomp(i)
          ic = 0
3641      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a,i5)')
     -     ' Enter the volume mixing ratio of component ',i
          read(*,*,err=3641) volmix(i)
          write(*,'(1x,1a,1pe13.5)') 'volmix =',volmix(i)
3661  continue
c
c****    d i s c r e t e    o r d i n a t e    m e t h o d
c
      ic = 0
3802  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a,/,1a)')
     - ' Enter number of streams for discrete ordinate method (> 2): ',
     - ' [half negative (downward) and half positive (upward)]'
      read(*,*,err= 3802) nstr
      write(*,'(1x,1a,i5)') 'nstr =',nstr
c
      if(nstr .gt. mxumu) then
        write(*,'(/,1a)') 
     -  ' number of streams exceeds dimension bound: mxumu=',mxumu
        stop
      endif
c
c****   s o u r c e    f u n c t i o n s
c
      ic = 0
4001  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(5(/,1a))') 
     - ' Choose types of source functions: ',
     - '   1) direct solar beam',
     - '   2) internal thermal sources',
     - '   3) both solar and thermal sources',
     - ' Select index of source type(s): '
      read(*,*,err=4001) isource
      write(*,'(1x,1a,i5)') 'isource =',isource
c
      if(isource .eq. 1 .or. isource .eq. 3) then
        if(isource .eq. 3) lplanck = .true.
c
c****     s o l a r    f l u x e s
c
        lsolar = .true.
        ic = 0
4041    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a,/,1a)') 
     -   ' Enter the solar flux file format index (1 - 3):  ',
     -   '  1) formatted 2) unformatted 3) list directed '
        read(*,*,err=4041) ifrms0
        write(*,'(1x,1a,i5)') 'ifrms0 =',ifrms0
c
        if(ifrms0 .eq. 1) then
          write(*,'(2(/,1a))')
     -      ' Enter format for reading solar fluxes',
     -      ' [enclosed in parenthesis, ie; (2f3.5) ]: '
          read(*,'(1a)') frms0
          write(*,'(1x,2a)') 'frms0 =',frms0
        else
          frms0 = ' '
        endif
c
        ic = 0
4062    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a)') ' Enter name of solar flux file: '
        read(*,'(1a)') solfile
        write(*,'(1x,2a)') 'solfile =',solfile
c
        write(*,'(/,1a,i5)') ' opening unit, iusol1=',iusol1
        if(ifrms0 .ne. 2) then
          open(iusol1,file=solfile,form='formatted',status='old',
     -          err=4062)
        else
          open(iusol1,file=solfile,form='unformatted',status='old',
     -         err=4062)
        endif
c
        close(iusol1)
c
        ic = 0
4081    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a)') 
     -    ' Enter number of lines to skip at top of file:'
        read(*,*,err=4081) ioffs0
        write(*,'(1x,1a,i5)') 'ioffs0 =',ioffs0
c
        ic = 0
4201    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(5(/,1a))') 
     -   ' Enter input radiance units of input solar fluxes: ',
     -   '  1) Watts/m**2/cm**-1        ',
     -   '  2) Watts/m**2/micron        ',
     -   '  3) Watts/m**2/nanometer     ',
     -   '  4) Watts/m**2/Angstrom      '
        read(*,*,err=4201) iuins0  
        write(*,'(1x,1a,i5)') 'iuins0 =',iuins0
c
        ic = 0
4221    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(4(/,1a))') 
     -   ' Specify type of abscissa for solar spectrum: ',
     -   '  1) wavelength (microns)',
     -   '  2) wavenumber',
     -   ' Enter index of choice: '
        read(*,*,err=4221) ixs0
        write(*,'(1x,1a,i5)') 'ixs0 =',ixs0
c
        ic = 0
        if(ixs0 .eq. 1) then
4241      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') 
     -    ' Enter factor to convert wavelengths to microns:'
          read(*,*,err=4241) scwns0
          write(*,'(1x,1a,1pe13.5)') 'scwns0 =',scwns0 
          ic = 0
4261      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') 
     -    ' Enter column numbers of wavelength and solar flux: '
          read(*,*,err=4261) icwns0,icsol
          write(*,'(1x,1a,2i5)') 'icwns0, icsol =',icwns0,icsol
        else
4281      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(/,1a)') 
     -    ' Enter column numbers of wavenumber and solar flux: '
          read(*,*,err=4281) icwns0,icsol
          write(*,'(1x,1a,2i5)') 'icwns0, icsol =',icwns0,icsol
          scwns0 = 1.0
        endif 
        scsol = 1.0  
c
        ic = 0
4402    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a)') 
     -  ' Enter number of solar incidence angles (> 0): '
        read(*,*,err=4402) nza
        write(*,'(1x,1a,i5)') 'nza =',nza
c
        if(nza .gt. nsol) then
          write(*,'(/,2a,i5)') 
     -    ' number solar zenith angles exceeds dimension',
     -    ' bound: nsol=',nsol
            stop
        endif
c
        write(*,'(/,2a,/,1a,/,1a,/1a)') 
     -  ' Enter solar zenith and azimuth angles (degrees)',
     -  ' for each in incidence angle: ',
     -  ' NOTE:  Solar zenith angles less than ~2 degrees and ',
     -  '        angles between ~58 and 62 degrees often cause ',
     -  '        near-singular matrices in the discrete ordinate model'
        do 4441 n=1,nza
            ic = 0
4421        ic = ic + 1
            if(ic .gt. icmx) stop
            write(*,'(/,1a,i5)') 
     -      ' Enter solar zenith and azimuth angles for path # ',n 
            read(*,*,err=4421) sza0(n),phi0(n)
            write(*,'(1x,1a,2(1pe13.5))') 'sza0, phi0', sza0(n),phi0(n)
            umu0(n) = cos(pi*sza0(n)/180.)
4441    continue
c
        write(*,'(/,1a)') 
     -   ' Enter convergence criteria for azimuthal series (0.0-0.01):'
        read(*,*) accur
        write(*,'(1x,1a,1pe13.5)') 'accur =',accur
      else
        lplanck = .true.
        lsolar = .false.
        nza = 1
      endif
c
      ic = 0
4621  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(16(/,1a))')
     - ' Specify quantities to be included in output files: ',
     - '   1) wavelength-dependent fluxes and radiances at the',
     - '      computational azimuths and zenith angles at the ',
     - '      specified output levels, and spectrally-integrated ',
     - '      fluxes and heating rates at each computational level.  ',
     - '   2) wavelength-dependent fluxes, radiances, and ',
     - '      transmission values at computational zenith angles ',
     - '      and specified output levels, and spectrally- ',
     - '      integrated fluxes and heating rates at each ',
     - '      computational level.  ',
     - '   3) wavelength-dependent fluxes and radiances at the ',
     - '      computational azimuths and zenith angles at the ',
     - '      specified output levels, wavelength-dependent,',
     - '      level-dependent, pressure-weighted flux divergences, ',
     - '      and spectrally-integrated fluxes and heating rates ',
     - '      at each computational level. '
      write(*,'(1a,19(/,1a))')
     - '   4) wavelength-dependent fluxes, radiances, ',
     - '      transmission values, wavelength-dependent, ',
     - '      level-dependent, pressure-weighted flux divergences, ',
     - '      and spectrally-integrated fluxes and heating rates ',
     - '      at each computational level. ',
     - '   5) same as (1) with radiances saved at arbitrary ',
     - '      azimuths and zenith angles. ',
     - '   6) same as (2) with radiances saved at arbitrary ',
     - '      azimuths and zenith angles. ',
     - '   7) same as (3) with radiances saved at arbitrary ',
     - '      azimuths and zenith angles. ',
     - '   8) same as (4) with radiances saved at arbitrary ',
     - '      azimuths and zenith angles.',
     - ' Enter index of choice: '
      read(*,*,err=4621) irad
      write(*,'(1x,1a,i5)') 'irad =',irad
c
      if(irad .le. 4) then 
        usrang = .false.
        numu = nstr
      else
c
c****      set disort flag that creates output at user-specified angles
c
        usrang = .true.
c
        ic = 0
4641    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a,/,1a)')
     -  ' Enter number of arbitrary emission zenith angles (> 0): ',
     -  ' [half negative (downward) and half positive (upward)]'
        read(*,*,err=4641) numu
        write(*,'(1x,1a,i5)') 'numu =',numu
c
        if(numu .gt. mxumu) then
          write(*,'(/,2a,i5)') 
     -    ' number of streams exceeds dimension',
     -    ' bound: mxumu=',mxumu
          stop
        endif
c
        write(*,'(/,1a,/,1a)') 
     -  ' Enter output emission zenith angles (180 - 0 degrees): ',
     -  ' (Note: Angles must be monotonically DECREASING) '
        do 4681 n=1,numu
            ic = 0
4661        ic = ic + 1
            if(ic .gt. icmx) stop
            write(*,'(1a,i5)') ' Enter emission zenith angle #',n
            read(*,*,err= 4661) ang
            umu(n) = cos(pi*ang/180.)
            write(*,'(1x,1a,1pe13.5)') 'ang =',ang
            write(*,'(1x,1a,1pe13.5)') 'umu = ',umu(n)
4681    continue
      endif
c
c****    specify emission azimuth angles
c
      ic = 0
4802  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)')
     -' Enter number of output emission azimuth angles (> 0): '
      read(*,*,err=4802) nphi
      write(*,'(1x,1a,i5)') 'nphi =',nphi
c
      if(nphi .gt. mxphi) then
        write(*,'(/,2a,i5)') 
     -  ' number of output azimuths exceeds dimension',
     -  ' bounds: mxphi=',mxphi
        stop
      endif
c
      write(*,'(/,1a)') 
     -' Enter output emission azimuth angles (degrees): '
      do 4841 n=1,nphi
          ic = 0
4821      ic = ic + 1
          if(ic .gt. icmx) stop
          write(*,'(1a,i5)') ' Enter emission azimuth angle #',n
          read(*,*,err= 4821) phi(n)
          write(*,'(1x,1a,1pe13.5)') 'phi =',phi(n)
4841  continue
c
c****    o u t p u t    l e v e l s
c
      ic = 0
5002  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(6(/1a))')
     - ' Specify level at which output radiances are to be printed: ',
     - '   1) top of the atmosphere only, ',
     - '   2) surface only, ',
     - '   3) top of atmosphere and surface,',
     - '   4) one or more (monotonically-increasing) pressure levels',
     - ' Choose index of output level: '
      read(*,*,err=5002) levout
      write(*,'(1x,1a,i5)') 'levout =',levout
c
      if(levout .le. 2) then
        nlout = 1
        if(levout .eq. 1) then
          clev(1) = '_toa'
        else
          clev(1) = '_sur'
        endif
      else
        if(levout .eq. 3) then
          nlout = 2
          clev(1) = '_toa'
          clev(2) = '_sur'
        else
          write(*,'(/,1a,/,1a,i5)') 
     -   ' Enter the number or arbitrary output levels: ',
     -   ' NOTE: the current maximum is mxlout = ',mxlout
          read(*,*,err=5002) nlout
          write(*,*) 'nlout = ',nlout
          if(nlout .gt. mxlout) then
            write(*,*) 'nlout exceeds dimension bound, mxlout=',mxlout
            stop
          endif
c
          do 5011 k=1,nlout
              write(*,'(/,1a,i5)') 
     -        ' Enter the output pressure (bars) for level #: ',k
              read(*,*,err=5002) pout(k)
              write(*,'(1x,1a,1pe13.5)') 'pout =',pout(k)
              write(clev(k),'(1a2,i2.2)') '_l',k
5011      continue
        endif
      endif
c
      ic = 0
5021  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(6(/,1a))') 
     - ' Enter index of desired output radiance units: ',
     - '  1) Watts/m**2/sr/cm**-1        ',
     - '  2) Watts/m**2/sr/micron        ',
     - '  3) Watts/m**2/sr/nanometer     ',
     - '  4) ergs/s/cm**2/sr/cm-1        ',
     - '  5) photons/s/m**2/sr/cm-1      '
      read(*,*,err=5021) iunits
      write(*,'(1x,1a,i5)') 'iunits =',iunits
c
c****    enter characteristics of the output spectral grid.
c
      ic = 0
6001  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,2a)')
     - ' Enter min and max wavenumbers of the ',
     - 'output spectral grid (cm**-1): '
      read(*,*,err=6001) wnmin,wnmax
      write(*,'(1x,1a,2(1pe14.6))') 'wnmin, wnmax =',wnmin,wnmax
c
      ic = 0
6201  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(4(/,1a))') 
     - ' Choose index of output spectrum type: ',
     - '  1) Full-resolution spectrum (this creates enourmous files)',
     - '  2) Sampled spectrum convolved with a slit function '
c     - '  3) SMT binned data with spectral map: '
      read(*,*,err=6201) isptype
      write(*,'(1x,1a,i5)') 'isptype =',isptype
c
      if(isptype .eq. 2) then
c
        ic = 0
6221    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(9(/,1a))') 
     - ' Enter type of spectral response function: ',
     - '  1) boxcar',
     - '  2) triangular (approximate slit spectrometer)'
c
        read(*,*,err=6221) islit
        write(*,'(1x,1a,i5)') 'islit =',islit
c
        ic = 0
6241    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a)') ' Enter the half-width-at-half-max: '
        read(*,*,err=6241) width
        write(*,'(1x,1a,1pe13.5)') 'width =',width
c
        ic = 0
6281    ic = ic + 1
        if(ic .gt. icmx) stop
        write(*,'(/,1a)') 
     -   ' Enter sampling resolution of output file (cm**-1):'
        read(*,*,err=6281) dwn
         write(*,'(1x,1a,1pe13.5)') 'dwn =',dwn
      endif
c
c****   set error limits for spectral mapping method
c
      ic = 0
6801  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter fractional error for tau binning: '
      read(*,*,err=6801) tauerr
      write(*,'(1x,1a,1pe13.5)') 'tauerr =',tauerr
      ic = 0
6821  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter error for pi0 binning: '
      read(*,*,err=6821) pi0err
      write(*,'(1x,1a,1pe13.5)') 'pi0err =',pi0err
      ic = 0
6841  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter error for <cos> binning: '
      read(*,*,err=6841) phferr
      write(*,'(1x,1a,1pe13.5)') 'phferr =',phferr
      ic = 0
6861  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a)') ' Enter error for surface optics binning: '
      read(*,*,err=6861) surferr
      write(*,'(1x,1a,1pe13.5)') 'surferr =',surferr
c
c****   n a m e    o f    o u t p u t    f i l e
c
      ic = 0
7001  ic = ic + 1
      if(ic .gt. icmx) stop
      write(*,'(/,1a,/,1a)') 
     - ' Enter index of output file type: ',
     - '  1) ASCII,   2) Binary with header,   3) binary no header'
      read(*,*,err=7001) ifrmout
      write(*,'(1a,i5)') ' ifrmout = ',ifrmout
c
      call init_spect_io(nza,nlout,irad,ifrmout,sza0,clev,
     -                         iuout,iuflx,iuheat,iustat,iutrn,
     -                         name,len,radfile,heatfile,statfile,
     -                         trnfile,flxfile) 
c
c****    define the first output zenith angle for radiances
c        (only upward radiances are saved at the top fo the atmosphere)
c
      if(usrang) then
        nzdn = 0
        do 8021 nze=1,numu
            if(umu(nze) .lt. 0.0) nzdn = nze
8021    continue
        nzup = nzdn + 1
      else
        nzdn = nstr/2
        nzup = nzdn + 1
      endif
c
      do 8011 nl=1,nlout
          if(levout .eq. 1 .or. (levout .eq. 3 .and. nl .eq. 1)) then
             nza_1(nl) = nzup
          else
             nza_1(nl) = 1
          endif
8011  continue
c
c****       create an output files for each partial derivative
c
      call init_pd_iu(nza,numu,nphi,nlout,ifrmout,
     -                      nza_1,radfile,lplanck,lsolar,
     -                      clev,nstate,istate,iu_pd,
     -                      j_ext,iutpd,iuspd,iupdrad)
c
      return
      end
