      subroutine cldstr(nlev,iopen,mode,iuaer,aerfile,
     -                  iofftau,iztau,icptau,ictau,
     -                  scptau,sctau,lcsh,iccsh,sccsh,
     -                  p,alt,z,dtauaer)
c
cccccccccccccccccccccccccccc  c l d s t r  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    ********************   p u r p o s e   *********************    cc
cc                                                                    cc
cc    this subroutine initializes the aerosol properties of the       cc
cc    atmosphere for radiative transfer models.                       cc
cc                                                                    cc
cc    ********************     i n p u t     *********************    cc
cc                                                                    cc
cc        iztau - vertical coordinate type 1) pressure, 2) altitude   cc
cc      iofftau - number of records to skip above optical depths      cc
cc       icptau - column with vertical coordinate                     cc
cc        ictau - column with differential optical depth              cc
cc       scptau - scaling factor for optical depth vertical           cc
cc                coordinate (pressures -> pascals, altitudes -> km)  cc
cc        sctau - scaling factor for optical depth                    cc
cc         lsch - scale height flag (true -> use scale heights)       cc
cc        icsch - column with scale heights.                          cc
cc        sccsh - scaling factor to convert scale heights to km       cc
cc        nclev - number of discrete cloud layers.                    cc
cc       cpl(l) - pressure (bars) at homogeneous aerosol layer boun-  cc
cc                daries                                              cc
cc     csh(m,l) - cloud particle scale height (km) for particle mode  cc
cc                m in homogeneous layer l                            cc
cc   dtaus(m,l) - aerosol differential optical depth for mode m in    cc
cc                layer l                                             cc
cc       ncloud - number of homogeneous aerosol layers                cc
cc       nmodes - number of aerosol particle modes                    cc
cc         tcld - total fractional cloudiness by all layers           cc
cc                                                                    cc
cc    ********************    o u t p u t    *********************    cc
cc                                                                    cc
cc                            - unit 6 -                              cc
cc                                                                    cc
cc   dtaus(m,l) - aerosol differential optical depth for mode m in    cc
cc                layer l                                             cc
cc                                                                    cc
cc                       - calling arguments -                        cc
cc                                                                    cc
cc                                                                    cc
cc                    - other computed quantities -                   cc
cc                                                                    cc
cc cpnd(m,l,na) - number density between model half-levels l and l+1  cc
cc                for mode m and model atmosphere na (particles/cm**3)cc
cc         nstd - spectral interval index where standard atmosphere   cc
cc                aerosol optical depths are specified                cc
cc dtauaer(m,l) - fractional optical depth for particle mode m, at    cc
cc                lth level in the nath model atmosphere              cc
cc                                                                    cc
cccccccccccccccccccccccccccc  c l d s t r cccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ncd
      parameter(ncd = 80)
c
      logical lcsh
c
c*****   file structure variables
c
      integer iofftau,iztau,icptau,ictau,iccsh
      integer nlev,iopen,mode,iuaer
c
      integer i,k,l,m,n,ierr,nzdim,ncdim,iform,invz,iskip,ncol,
     -        ncmax,nptj,nxmax,nymax,ncloud,nclev,ndiv,icl,ic
c
      real scptau,sctau,sccsh,p(kp),alt(kp),z(kp)
      real altc0(100*kp),tauc0(100*kp),cldabs(100*kp),
     -       altaer(kp),tauaer(kp), dtauaer(nmode,kp)
      real z1,dtaus1,fndiv,dz,dz10,cshm,cldab0,tautot,pbar
c
c****   readfil variables
c
      character*132 aerfile
      character*40 formxy
      integer icol(ncd)
      real scz(ncd),add(ncd),array(1,ncd)
c
      real zaer(kp),altc(kp),dtaus(kp),csh(kp)
c
      m = mode
c
c*****    initialize i/o variables
c
      ierr = 0
      nzdim = 1
      ncdim = ncd
      formxy = ' '
      iform = 3
      invz = 2
      iskip = 0
      if(lcsh) then
        ncol = 3
        icol(3) = iccsh
        scz(3) = sccsh
      else
        ncol = 2
        do 1001 k=1,kp
            csh(k) = 0.0
1001    continue
      endif
      icol(1) = icptau
      icol(2) = ictau
      scz(1) = scptau
      scz(2) = sctau
      add(1) = 0.
      add(2) = 0.
      ncmax = icol(1)
      do 1201 i=2,ncol
          if(icol(i) .gt. ncmax) ncmax = icol(i)
1201  continue
c
      if(iopen .eq. 1) then
c
c****    o p e n   a e r o s o l   s t r u c t u r e    f i l e 
c
        call readfil(0,iuaer,aerfile,formxy,iform,nzdim,ncdim,
     -               iofftau,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
        if(ierr .ne. 0) then
          write(*,'(/,2a)') ' cldstr error opening file: ',aerfile
          stop
        endif
      endif
c
      if(iofftau .ge. 1) then
c
c****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s .
c
        call readfil(-1,iuaer,aerfile,formxy,iform,nzdim,ncdim,
     -               iofftau,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) then
          write(*,*) 
     -    'End-of-file Encountered while skipping records in file ',
     -             aerfile
          stop
c
        endif
      endif
c
c*****   read each p-tau record 
c
      k = 0
      dtaus(1) = 0.0
c
2402  nptj = 1
      call readfil(1,iuaer,aerfile,formxy,iform,nzdim,ncdim,
     -               iofftau,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) go to 2801
        k = k + 1
c
c***      load vectors for vertical coordinate and aerosol amounts.
c
        if(iztau .eq. 1) then
c
c****      the vertical coordinate is pressure
c          define the log of the cloud-level pressure
c
          if(array(1,1) .gt. 0.) then
            zaer(k) = alog(array(1,1))
          else
            zaer(k) = 1.e10
          endif
c
        else
c
c****     define altitude of input
c
          zaer(k) = array(1,1)
c
        endif
c
c****     define the differential optical depth above level k
c
        if(k .gt. 1) dtaus(k-1) = array(1,2)
c
        if(lcsh) csh(k) = array(1,3)
c
      if(k .lt. nclev .or. nclev .eq. 0) go to 2402
c
2801  nclev = k
      ncloud = nclev - 1
c
      write(*,'(1a,i5,1a,i5)') 
     - ' Number of input levels for aerosol mode ',mode,' = ',nclev
c
c****   find altitudes corresponding to the aerosol layer boundaries
c
      if(iztau .eq. 1) then
c
c****     the input vertical coordinate, zaer, is log pressure.
c         interpolate to the altitude grid
c
        nxmax = kp
        nymax = 1
c
        call xyinterp(z,alt,zaer,altc,nxmax,nymax,nlev,nclev,1)
c
      else
c
        if(zaer(nclev) .lt. zaer(1)) then
          do 3001 k=1,nclev
              altc(k) = zaer(k)
3001      continue
        else
c
c****       the array must be turned around
c
          do 3021 k=1,nclev/2
              z1 = zaer(k)
              dtaus1 = dtaus(k)
              zaer(k) = zaer(nclev-k+1)
              dtaus(k) = dtaus(nclev-k+1)
              zaer(nclev-k+1) = z1
              dtaus(nclev-k+1) = dtaus1
3021      continue
c
          do 3041 k=1,nclev
              altc(k) = zaer(k)
3041      continue
        endif
      endif
c
c****   interpolate grid to 10 times the vertical resolution.  use the
c       cloud particle scale height data to find the optical depth
c       profile in each sub-layer.
c
      ndiv = 100*kp/ncloud - 1
      fndiv = float(ndiv)
      icl = 0
      do 3221 l=1,ncloud
          do 3201 ic=1,ndiv
              icl = icl + 1
              tauc0(icl) = 0.0
              cldabs(icl) = 0.0
3201      continue
3221  continue
c
      icl = 1
      if(alt(1) .gt. altc(1)) then
        altc0(icl) = alt(1)
        tauc0(icl) = 0.0
        cldabs(icl) = 0.0
        icl = icl + 1
      endif
c
      altc0(icl) = altc(1)
      tauc0(icl) = 0.0
      cldabs(icl) = 0.0
c
      do 4021 n=1,ncloud
c
c****        determine the thickness of the cloud layer
c
          dz = altc(n) - altc(n+1)
c
c****       determine if the layer has a constant or
c           variable number density, and find the absorption
c           coefficient per unit length at the base of the layer.
c
          dz10 = dz/fndiv
          cshm = csh(n+1)
          if(abs(cshm) .gt. 1.e-5) then
c
c****          use the cloud particle scale heights to distribute 
c              cloud particles
c
            cldab0 = dtaus(n)/(cshm*(1. - exp(-dz/cshm)))
          endif
c
          do 4001 i=1,ndiv
              icl = icl + 1
              altc0(icl) = altc(n) - float(i)*dz10
              if(abs(cshm) .gt. 1.e-5) then
c
c****             use the particle scale height to distribute
c                 the particles in the layer
c
                cldabs(icl) = cldab0*exp(-(altc0(icl)-altc(n+1))/cshm)
              else
c
c****             the scale height is specified as zero -assume
c                 particles have a constant number density.
c
                cldabs(icl) = dtaus(n)/dz
c
              endif
              tauc0(icl) = tauc0(icl-1) + 0.5*(cldabs(icl) +
     -                     cldabs(icl-1))*(altc0(icl-1) - altc0(icl))
4001      continue
4021  continue
c
c****   interpolate the cloud absorption coefficients
c       to the pressure grid.
c
      do 4041 k=1,nlev
          altaer(k) = alt(k)
4041  continue
c
      nxmax = 100*kp
      nymax = 1
c
      call xyinterp(altc0,tauc0,altaer,tauaer,
     -                  nxmax,nymax,icl,nlev,nymax)
c
c****   find the differential optical depth in each layer.
c
      do 4061 k=1,nlev-1
          dtauaer(m,k) = tauaer(k+1) - tauaer(k)
4061  continue
c
c****   print the cumulative optical depth for each particle mode
c
      write(*,'(/,1a,i5,/,/,1a)')
     - ' Optical depth for particle mode:',mode,
     -    '    alt(km)     p (bar)      dtau        tau'
      tautot = 0.0
      do 4281 k=1,nlev-1
           pbar = 1.e-5*p(k+1)
           tautot = tautot + dtauaer(m,k)
           write(*,'(4(1pe12.4))') alt(k+1),pbar,dtauaer(m,k),tautot
4281  continue
c
      return
      end
