      subroutine readsur(surfile,iusur,iopen,iform,formxy,ioff,
     -                   nalbmx,iwn,icwn,nref,icalb,scwn,scalb,
     -                   wnmin,wnmax,io_end,io_err,wneof)
c
ccccccccccccccccccccccccccc  r e a d s u r  ccccccccccccccccccccccccccc
cc                                                                   cc
cc    ********************   p u r p o s e   *********************   cc
cc                                                                   cc
cc    this subroutine reads surface optical properties as a function cc
cc    ofwavelength or wavenumber.  These quantities are then         cc
cc    rewritten into a binary file that lists albedos monotonically  cc 
cc    with wavenumber.                                               cc
cc                                                                   cc
cc    ********************     i n p u t     *********************   cc
cc                                                                   cc
cc    iusur - input unit number for surface albedos and              cc
cc             moments of surface phase functions.                   cc
cc     iusur - output unit number for surface albedos.               cc
cc     iopen - open unit flag - 0) don't open unit, 1) open unit.    cc
cc   surfile - name of input file                                    cc
cc       iwn - type of frequency coordinate for input albedos:       cc
cc              1) wavelength, 2) wavenumber.                        cc
cc      icwn - column index of spectral quantity in each record.     cc
cc     iclab - column index of surface optics in each record.        cc
cc     iform - format of input file:                                 cc
cc             1) formatted, 2) unformatted, 3) list directed)       cc
cc    formxy - format for formatted data files (iform=1 only)        cc
cc      scwn - multiplicative factor to convert wavelength - microns cc
cc     scalb - multiplicative factor to convert albedo to range(0-1).cc
cc      ioff - number of records to skip at top of file              cc
cc    nalbmx - maximum number of albedos in input file               cc
cc   wn_0(l) - wavelength or wavenumber of each input albedo.        cc
cc    surf_0 - input surface optical properties                      cc
cc     wnmin - minimum wavenumber where albedos are needed (cm**-1)  cc
cc     wnmax - maximum wavenumber where albedos are needed (cm**-1)  cc
cc                                                                   cc
cc    ********************   o u t p u t     *********************   cc
cc                                                                   cc
cc  wnalb(l) - wavenumber of each input albedo.                      cc
cc    surf_0 - input surface optical properties at each wavelength   cc
cc                                                                   cc
ccccccccccccccccccccccccccc  r e a d s u r cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer ncd
      parameter(ncd = 80)
c
      character*(*) surfile
      character*(*) formxy
      integer iusur,iopen,iform,ioff,nalbmx,iwn,icwn,
     -        io_end,io_err,nref,icalb(nref)
      integer ierr,nr,nptj,k
      integer ncdim,nzdim,iskip,invz,icol(ncd),ncmax,ncol
c
      double precision wnmin,wnmax,wneof(2)
c
      real scwn,scalb(nref)
      real scz(ncd),add(ncd),array(1,ncd)
c
      real surf_0(32766,4)
c
      double precision wn_0(32766)
c
c*****    initialize i/o variables
c
      ierr = 0
c
      ncdim = ncd
      nzdim = 1
      iskip = 0
      invz = 2
      icol(1) = icwn
      ncmax = icol(1)
      scz(1) = scwn
      add(1) = 0.
      ncol = nref + 1
      do 1001 nr=1,nref
          icol(nr+1) = icalb(nr)
          if(icol(nr+1) .gt. ncmax) ncmax = icol(nr+1)
          scz(nr+1) = scalb(nr)
          add(nr+1) = 0.0
1001  continue
c
      if(iopen .eq. 1) then
c
c****    o p e n   s u r f a c e    a l b e d o    f i l e 
c
        call readfil(0,iusur,surfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) then
          write(*,'(/,2a)') ' readsur: error opening file: ',surfile
          stop
        endif
      endif
c
c****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s .
c
      if(ioff .ge. 1) then
c
        call readfil(-1,iusur,surfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) then
          write(*,'(/,2a)') 
     -    'End-of-file Encountered while skipping records in file ',
     -             surfile
          stop
c
        endif
      endif
c
c*****   read each surface albedo record 
c
      k = 0
c
2001  nptj = 1
c
      call readfil(1,iusur,surfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
      if(ierr .eq. 0) then
c
c***      load wavenumber vector
c
        if(iwn .eq. 1) then
          if(array(1,1) .ne. 0.) then
            wn_0(k+1) = 1.0d4/array(1,1)
          else
            wn_0(k+1) = 1.0d10
          endif
        else
          wn_0(k+1) = array(1,1)
        endif
c
c****     determine if this value is needed.
c
        if(k .eq. 0) then
c
          k = k + 1
c
c****      initialize the very first point.
c
          do 2281 nr=1,nref
              surf_0(k,nr) = array(1,nr+1)
2281      continue
c
        else  !  if(k .eq. 0) k is not equal to zero
c
          if(wn_0(k+1) .gt. wn_0(k)) then
c
c****        input data is ordered in order of increasing wavenumber.
c
            if(wn_0(k+1) .lt. wnmin) then
c
c****              wavenumber is too small - reset intial value 
c
              k = 1
              wn_0(k) = wn_0(k+1)
              do 2201 nr=1,nref
                  surf_0(k,nr) = array(1,nr+1)
2201          continue
c
            else
c
c****              wavenumber is too large - add this point and quit
c
              k = k + 1
              do 2221 nr=1,nref
                  surf_0(k,nr) = array(1,nr+1)
2221          continue
c
              if(wn_0(k) .gt. wnmax) go to 2402
c
            endif ! if(wn_0(k+1) .lt. wnmin)
c
          else ! wn_0(k+1) < wn_0(k)
c
c****        input data is ordered in order of increasing wavelength
c
            if(wn_0(k+1) .gt. wnmax) then
c
c****          wavenumber is too large - reset intial value
c
              k = 1
              wn_0(k) = wn_0(k+1)
              do 2241 nr=1,nref
                  surf_0(k,nr) = array(1,nr+1)
2241          continue
c
            else ! if(wn_0(k+1) .gt. wnmax)
c
c****              add this point 
c
                k = k + 1
                do 2261 nr=1,nref
                    surf_0(k,nr) = array(1,nr+1)
2261            continue
c
              if(wn_0(k) .lt. wnmin) go to 2402
c
            endif ! if(wn_0(k+1) .gt. wnmax)
c
          endif!  if(wn_0(k+1) .gt. wn_0(k))
c
        endif !  if(k .eq. 0)
c
        if(k .lt. nalbmx .or. nalbmx .eq. 0) go to 2001
c
      endif  ! if(ierr .eq. 0)
c
c****    You are done.  reaorder arrays if necessary and close
c
2402  nalbmx = k
      if(ierr .eq. 1) then 
        io_end = 1
      else
        io_err = 1
      endif
c
      write(*,'(/,1a,i5,/)') 
     - ' Number of input surface optical properties = ',nalbmx
c
      close(iusur)
c
c****   open a scratch unit for albedos, monotonically ordered in 
c       increasing wavenumber.
c
      open(iusur,form='unformatted',status='scratch')
c
      if(wn_0(1) .lt. wn_0(nalbmx)) then
        do 3001 k=1,nalbmx
           write(iusur) wn_0(k),(surf_0(k,nr),nr=1,nref)
3001    continue
        wneof(1) = wn_0(1)
        wneof(2) = wn_0(nalbmx)
      else
        do 3021 k=nalbmx,1,-1
            write(iusur) wn_0(k),(surf_0(k,nr),nr=1,nref)
3021    continue
        wneof(1) = wn_0(nalbmx)
        wneof(2) = wn_0(1)
      endif
c
      rewind(iusur)
c
      return
      end
