      subroutine atmstr(atmfile,iuatm,iopen,iform,formxy,
     -                   ioff,nlmax,icp,ict,scp,ratm,
     -                   radius,sgrav,p,t,z,alt,grav,np_mx,nlev)
c
cccccccccccccccccccccccccccc  a t m s t r  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    ********************   p u r p o s e   *********************    cc
cc                                                                    cc
cc    this subroutine reads in standard atmosphere pressures and      cc
cc    temperatures, and uses the hydrostatic equation to find the     cc
cc    altitudes and gravitational acceleration at each pressure level.cc
cc                                                                    cc
cc    ********************     i n p u t     **********************   cc
cc                                                                    cc
cc     iuatm - input unit number for atmospheric variables.           cc
cc     iopen - open unit flag - 0) don't open unit, 1) open unit.     cc
cc   atmfile - name of input file                                     cc
cc       icp - column index of p value in each record.                cc
cc       ict - column index of t value in each record.                cc
cc     iform - format of input file:                                  cc
cc             1) formatted, 2) unformatted, 3) list directed)        cc
cc    formxy - format for formatted data files (iform=1 only)         cc
cc       scp - multiplicative factor to convert pressure to pascals.  cc
cc      ioff - number of records to skip at top of file               cc
cc      nlev - number of levels in the input atmosphere:              cc
cc             (if nlev=0, the ps in adopted as the input p field.    cc
cc             otherwise, values are interpolated to input p-grid.    cc
cc     nlmax - maximum number of input levels to read from file.      cc
cc             if this variable is set to zero, the program reads     cc
cc             all levels in this file.                               cc
cc     ps(l) - pressure at level l of standard atmosphere (pascals)   cc
cc     ts(l) - temperature at level l of standard atmosphere (k)      cc
cc      ratm - gas constant of atmosphere (J/kg/K)                    cc
cc     sgrav - surface gravity (m/s**2)                               cc
cc    radius - radius of planet (km)                                  cc
cc                                                                    cc
cc    ********************   o u t p u t     ********************     cc
cc                                                                    cc
cc       alt(l) - altitude above planet's surface of kth level (km)   cc
cc         t(l) - temperature at model level l (k)                    cc
cc         p(l) - pressure at model level l (pascals)                 cc
cc         z(l) - log of pressure at each level, l                    cc
cc      grav(l) - gravitational acceleraton at each level l (m/s**2)  cc
cc                                                                    cc
cccccccccccccccccccccccccccc  a t m s t r cccccccccccccccccccccccccccccc
c
      implicit none
c
      integer ncd, npmax 
      parameter(ncd = 20, npmax = 200)
c
      character*132 atmfile
      character*40 formxy
c
      integer np_mx
c
      integer icol(ncd),nl0
c
      integer iuatm,ioff,iopen,iform,nlmax,icp,ict,nlev
c
      integer ierr,ncol,ncdim,nzdim,iskip,invz,nptj,
     -        ncmax,nxmax,nymax,k,l
c
      real scz(ncd),add(ncd),array(1,ncd)
c
      real zs(npmax),ps(npmax),ts(npmax)
c
c****   atmospheric structure variables
c
      real scp,p(np_mx),t(np_mx),alt(np_mx),grav(np_mx),z(np_mx),
     -     ratm,radius,sgrav
c
      real p5ratm
c
c*****    initialize i/o variables
c
      ierr = 0
      ncol = 2
      ncdim = ncd
      nzdim = 1
      iskip = 0
      invz = 2
      nptj = 1
      icol(1) = icp
      icol(2) = ict
      ncmax = icol(2)
      if(icol(1) .gt. icol(2)) ncmax = icol(1)
      scz(1) = scp
      scz(2) = 1.0
      add(1) = 0.
      add(2) = 0.
c
      if(iopen .eq. 1) then
c
c****    o p e n   a t m o s p h e r i c   s t r u c t u r e    f i l e 
c
        call readfil(0,iuatm,atmfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
        if(ierr .ne. 0) then
          write(*,*) 'atmstr error opening file: ',atmfile
          stop
        endif
      endif
c
c****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s .
c
      if(ioff .ge. 1) then
        call readfil(-1,iuatm,atmfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) then
          write(*,*) 
     -        'End-of-file Encountered while skipping records in file ',
     -             atmfile
          stop
c
        endif
      endif
c
c*****   read each p-T record 
c
      nlmax = np_mx
      nl0 = 0
c
2401  nptj = 1
      call readfil(1,iuatm,atmfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) go to 2801
        nl0 = nl0 + 1
c
c****       determine if the number of vertical levels exceeds
c           the dimension bound:  Note, some of the partial 
c           deriviative arrays use nlev + 1 levels, so the
c           maximum number of levels, nlev must be less than kp
c
        if(nl0 .ge. np_mx) then
          write(*,'(/,/,5(1a,/),1a,i5,/,1a,2e13.5)') 
     -    ' ****  E R R O R  ****    E R R O R  ****  E R R O R  ****', 
     -    ' The number of levels in the model atmosphere exceeds the',
     -    ' vertical level dimension, np_mx, in the subroutine ATMSTR.',
     -    ' Reduce the number of levels in the model atmosphere, or', 
     -    ' change np_mx in the calling program and recompile.',
     -    ' current dimension bound, np_mx =',np_mx,
     -    ' pressure and temperature at last level read =',
     -      array(1,1),array(1,2)
c
          stop
        endif
c
c***      load p-t vectors
c
        ps(nl0) = array(1,1)
        ts(nl0) = array(1,2)
c
c****      define the log of pressure
c
        if(ps(nl0) .gt. 0.) then
          zs(nl0) = alog(ps(nl0))
        else
          zs(nl0) = -1.e10
        endif
c
      if(nl0 .le. nlmax .or. nlmax .eq. 0) go to 2401
c
2801  continue
      write(*,'(1a,i5)') 
     - ' Number of input levels in p-T file = ',nl0
c
      if(nlev .eq. 0) then
c
c****     use the input grid as the default pressure grid
c
        nlev = nl0
        do 3001 k=1,nlev
            p(k) = ps(k)
            t(k) = ts(k)
            z(k) = zs(k)
3001    continue
      else
c
c****     interpolate the input temperature field to the input
c         pressure grid.  Assume temperatures vary linearly
c         with the log of pressure.
c
        do 3201 k=1,nlev
          if(p(k) .gt. 0.) then
            z(k) = alog(p(k))
          else
            z(k) = -1.e10
          endif
3201    continue
c
        nxmax = np_mx
        nymax = 1
c
        call xyinterp(zs,ts,z,t,nxmax,nymax,nl0,nlev,1)
c
      endif
c
c****   solve the hydrostatic equation to find the altitude and
c       gravitational acceleration at each pressure level.
c
      p5ratm = 0.0005*ratm
      alt(nlev) = 0.0d0
      grav(nlev) = sgrav
      do 4241 l=nlev-1,1,-1
          alt(l) = alt(l+1) + p5ratm*(t(l+1) + t(l))*
     -            (z(l+1) - z(l))/grav(l+1)
          grav(l) = sgrav*(radius/(radius+alt(l)))**2
4241  continue
c
c****  if top pressure is zero, set altitude there to a nominal value.
c
      if(p(1) .eq. 0.0) then
        alt(1) = 2.*alt(2)
        grav(1) = sgrav*(radius/(radius+alt(1)))**2
      endif
c
      return
      end
