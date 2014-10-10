      subroutine readmix(mixfile,iumix,ng,iopen,ifrmmix,
     -                   frmmix,ioffmix,nlmax,izmix,
     -                   imix,icpmix,icmix,
     -                   scpmix,scmix,wgtatm,wgtgas,
     -                   nlev,alt,z,rmix)
c
cccccccccccccccccccccccccccc r e a d m i x ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    ********************   p u r p o s e   **********************   cc
cc                                                                    cc
cc    this subroutine reads reads gas mixing ratios as a function of  cc
cc    pressure, converts pressures to pascals, and mixing ratios to   cc
cc    mass mixing ratios, if necessary, and interpolates them to a    cc
cc    standard pressure grid.                                         cc
cc                                                                    cc
cc    ********************     i n p u t     *********************    cc
cc                                                                    cc
cc     iumix - input unit number for atmospheric variables.           cc
cc        ng - hitran gas index for this gas.                         cc
cc     iopen - open unit flag - 0) don't open unit, 1) open unit.     cc
cc   mixfile - name of input file                                     cc
cc     izmix - idex of vertical coordinate: 1) pressure, 2) altitude. cc
cc      imix - type of mixing ratios: 1) volume, 2) mass (kg/kg)      cc
cc    icpmix - column index of vertical coordinate in each record.    cc
cc     icmix - column index of mixing ratio in each record.           cc
cc     ifrmmix - format of input file:                                cc
cc             1) formatted, 2) unformatted, 3) list directed)        cc
cc    frmmix - format for formatted data files (ifrmmix=1 only)       cc
cc    scpmix - multiplicative factor to convert pressure to pascals   cc
cc             or altitudes to km.                                    cc
cc     scmix - multiplicative factor to convert mixing ratio to       cc
cc             range (0.0 to 1.0)                                     cc
cc      ioffmix - number of records to skip at top of file               cc
cc      nlev - number of levels in the standard atmosphere p grid.    cc
cc     nlmax - maximum number of levels to read from input file.      cc
cc     ps(l) - pressure at level l of input rmix atmosphere.          cc
cc  rmixs(l) - mixing ratio at level l of input atmosphere.           cc
cc    wgtgas - molecular weight of gas (kg/kmole).                    cc
cc    wgtatm - molecular weight of background atmosphere (kg/kmole).  cc
cc                                                                    cc
cc    ********************   o u t p u t     *********************    cc
cc                                                                    cc
cc      rmix(l) - mass mixing ratio at model level l (kg/kg).         cc
cc                                                                    cc
cccccccccccccccccccccccccccc r e a d m i x ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ncd
      parameter(ncd = 20)
c
      character*132 mixfile
      character*40 frmmix
c
      integer icol(ncd)
      integer iumix,ng,iopen,ifrmmix,ioffmix,nlmax,
     -        izmix,imix,icpmix,icmix
      integer ierr,ncol,ncdim,nzdim,invz,iskip,ncmax,nptj,k,
     -        nlstd,nxmax,nymax
c
      real scpmix,scmix,wgtatm
      real scz(ncd),add(ncd),array(1,ncd)
c
      real zs(8096),zmix(kp),rmix0(8096),rmix1(kp)
c
c****   atmospheric structure variables
c
      real alt(kp),z(kp),rmix(kp,ngas),wgtgas(ngas)
      integer nlev
c
c*****    initialize i/o variables
c
      ierr = 0
      ncol = 2
      ncdim = ncd
      nzdim= 1
      invz = 2
      iskip = 0
      icol(1) = icpmix
      icol(2) = icmix
      ncmax = icol(2)
      if(icol(1) .gt. icol(2)) ncmax = icol(1)
      scz(1) = scpmix
      scz(2) = scmix
      add(1) = 0.
      add(2) = 0.
c
      if(iopen .eq. 1) then
c
c****    o p e n   a t m o s p h e r i c   s t r u c t u r e    f i l e 
c
        call readfil(0,iumix,mixfile,frmmix,ifrmmix,nzdim,ncdim,
     -               ioffmix,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
        if(ierr .ne. 0) then
          write(*,*) 'readmix error opening file: ',mixfile
          stop
        endif
      endif
c
c****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s .
c
      if(ioffmix .ge. 1) then
        call readfil(-1,iumix,mixfile,frmmix,ifrmmix,nzdim,ncdim,
     -               ioffmix,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) then
          write(*,*) 
     -        'End-of-file Encountered while skipping records in file ',
     -             mixfile
          stop
c
        endif
      endif
c
c*****   read each p-T record 
c
      k = 0
c
2401  nptj = 1
      call readfil(1,iumix,mixfile,frmmix,ifrmmix,nzdim,ncdim,
     -               ioffmix,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) go to 2801
        k = k + 1
c
c***      load vectors for vertical coordinate and mixing ratios.
c
        if(izmix .eq. 1) then
c
c****      define the log of pressure
c
          if(array(1,1) .gt. 0.) then
            zs(k) = alog(array(1,1))
          else
            zs(k) = -1.e10
          endif
        else
c
c****     define altitude of input
c
          zs(k) = array(1,1)
        endif
c
c****     define the log of the mass mixing ratio
c
        if(array(1,2) .gt. 0.0) then
          if(imix .eq. 1) then
            rmix0(k) = alog(wgtgas(ng)*array(1,2)/wgtatm)
          else
            rmix0(k) = alog(array(1,2))
          endif
        else
          rmix0(k) = -100.
        endif
c
      if(k .lt. nlmax .or. nlmax .eq. 0) go to 2401
c
2801  nlstd = k
c
      write(*,'(1a,i5,1a,i5)') 
     - ' Number of input levels for gas ',ng,' = ',nlstd
c
c****   interpolate the input mixing ratios to the input
c       pressure grid.  Assume log of mixing ratio varies linearly
c       with the log of pressure.
c
      if(izmix .eq. 1) then
        do 3201 k=1,nlev
            zmix(k) = z(k)
3201    continue
      else
        do 3221 k=1,nlev
            zmix(k) = alt(k)
3221    continue
      endif
c
      nxmax = kp
      nymax = 1
c
      call xyinterp(zs,rmix0,zmix,rmix1,nxmax,nymax,nlstd,nlev,1)
c
c****     convert output back to mixing ratio
c
      do 4001 k=1,nlev
          rmix(k,ng) = exp(rmix1(k))
4001  continue
c
      return
      end
