      subroutine readsol(solfile,iusol1,iusol2,iopen,iform,formxy,
     -                   ioff,iwn,icwn,icsol,scwn,scsol,
     -                   iuin,wnmin,wnmax,au,
     -                   io_end,io_err,wneof)
c
cccccccccccccccccccccccccccc r e a d s o l cccccccccccccccccccccccccccc
cc                                                                   cc
cc    ********************   p u r p o s e   *********************   cc
cc                                                                   cc
cc    this subroutine reads solar fluxes as a function of wavelength cc
cc    or wavenumber.  These quantities are then rewritten into       cc
cc    a binary file that lists fluxes monotonically with wavenumber. cc
cc                                                                   cc
cc    ********************     i n p u t     *********************   cc
cc                                                                   cc
cc    iusol1 - input unit number for solar fluxes.                   cc
cc     iopen - open unit flag - 0) don't open unit, 1) open unit.    cc
cc   solfile - name of input file                                    cc
cc       iwn - type of frequency coordinate for input solar fluxes   cc
cc              1) wavelength, 2) wavenumber.                        cc
cc      iuin - units index for input solar fluxes:                   cc
cc              1) Watts/m**2/sr/cm**-1                              cc
cc              2) Watts/m**2/sr/micron                              cc
cc              3) Watts/m**2/sr/nanometer                           cc
cc              4) ergs/s/cm**2/sr/cm-1                              cc
cc              5) photons/s/m**2/sr/micron                          cc
cc      icwn - column index of spectral quantity in each record.     cc
cc     iclab - column index of solar flux in each record.            cc
cc     iform - format of input file:                                 cc
cc             1) formatted, 2) unformatted, 3) list directed)       cc
cc    formxy - format for formatted data files (iform=1 only)        cc
cc      scwn - multiplicative factor to convert wavelength - microns cc
cc      ioff - number of records to skip at top of file              cc
cc     nsol0 - maximum number of solar fluxe in input file           cc
cc    wn0(l) - wavelength or wavenumber of each input solar flux.    cc
cc   sol0(l) - input solar flux value                                cc
cc     wnmin - minimum wavenumber where fluxes are needed (cm**-1)   cc
cc     wnmax - maximum wavenumber where fluxes are needed (cm**-1)   cc
cc        au - relative distance from the sun in astronomical units  cc
cc                                                                   cc
cc    ********************   o u t p u t     *********************   cc
cc                                                                   cc
cc    wn0(l) - wavenumber of each input solar flux.                  cc
cc   sol0(l) - solar flux in units w/m/m/sr/cm-1 at distance, au     cc
cc                                                                   cc
cccccccccccccccccccccccccccc r e a d s o l cccccccccccccccccccccccccccc
c
      implicit none
c
      integer ncd,nsp
      parameter(ncd = 20, nsp = 32767)
c
      character*132 solfile
      character*40 formxy
c
      integer iusol1,iusol2,iopen,iform,icol(ncd),iwn,icwn,icsol,ioff
      integer iuin,io_end,io_err
      integer iunits,ius0,nsol0,nstot,ncol,ierr,ncdim,nzdim,iskip,invz,
     -        ncmax,nptj,k,iu
c
      real scwn,scsol,scz(ncd),add(ncd),array(1,ncd)
c
      double precision wnmin,wnmax,wneof(2)
c
      double precision wn0(nsp)
c
      real au,au2,sol0(nsp),units
c
      ius0 = 0
      nsol0 = 0
      nstot = 0
c
c****  set the default units index iunits = 1 to yield internal solar
c      flux units of W/m**2/cm-1
c
      iunits = 1
c
c****   open a scratch unit for solar flux, monotonically ordered in 
c       increasing wavenumber.
c
      open(iusol2,form='unformatted',status='scratch')
c
c****   find the square of the relative distance to the sun:
c
      au2 = au*au
      if(au2 .eq. 0.0) au2 = 1.0
c
c*****    initialize i/o variables
c
      ierr = 0
      ncol = 2
      ncdim = ncd
      nzdim = 1
      iskip = 0
      invz = 2
      icol(1) = icwn
      icol(2) = icsol
      ncmax = icol(2)
      if(icol(1) .gt. ncmax) ncmax = icol(1)
      if(scwn .ne. 0.0) then
        scz(1) = scwn
      else
        scz(1) = 1.0
      endif
      scz(2) = scsol/au2
      add(1) = 0.
      add(2) = 0.
c
      if(iopen .eq. 1) then
c
c****    o p e n   a t m o s p h e r i c   s t r u c t u r e    f i l e 
c
        call readfil(0,iusol1,solfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
        if(ierr .ne. 0) then
          write(*,*) 'readsol error opening file: ',solfile
          stop
        endif
      endif
c
c****   s k i p   o v e r   u n n e c e s s a r y   r e c o r d s .
c
      if(ioff .ge. 1) then
        call readfil(-1,iusol1,solfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) then
          write(*,*) 
     -        'End-of-file Encountered while skipping records in file ',
     -             solfile
          stop
c
        endif
      endif
c
c*****   read each solar flux record 
c
      k = 0
c
2403  nptj = 1
      call readfil(1,iusol1,solfile,formxy,iform,nzdim,ncdim,
     -               ioff,iskip,invz,ncol,icol,nptj,scz,add,
     -               array,ncmax,ierr)
c
        if(ierr .ne. 0) go to 2803
c
c***      load wavenumber vector
c
        if(iwn .eq. 1) then
          if(array(1,1) .ne. 0.) then
            wn0(k+1) = 1.0d4/array(1,1)
          else
            wn0(k+1) = 1.0d10
          endif
        else
          wn0(k+1) = array(1,1)
          if(k .gt. 0) then
            if(wn0(k+1) .eq. wn0(k)) go to 2403
          endif
        endif
c
c****     determine if this value is needed.
c
        if(k .ge. 1) then
          if(wn0(k+1) .gt. wn0(k)) then
c
c****        input data is ordered in order of increasing wavenumber.
c
            if(wn0(k+1) .lt. wnmin) then
c
c****            wavenumber is too small - discard
c
               k = 1
               wn0(k) = wn0(k+1)
               sol0(k) = array(1,2)
               go to 2403
            else
              if(wn0(k+1) .gt. wnmax) then
c
c****              wavenumber is too large - quit
c
                 k = k + 1
                 sol0(k) = array(1,2)
                 go to 2803
              endif
            endif
          else
c
c****      input data is ordered in order of increasing wavelength
c
            if(wn0(k+1) .gt. wnmax) then
c
c****            wavenumber is too large - discard
c
              k = 1
              wn0(k) = wn0(k+1)
              sol0(k) = array(1,2)
              go to 2403
            else
              if(wn0(k+1) .lt. wnmin) then
c
c****              wavenumber is too small - quit
c
                k = k + 1
                sol0(k) = array(1,2)
                go to 2803
              endif
            endif
          endif                
        endif
c
        k = k + 1
c
        sol0(k) = array(1,2)
c
      if(k .lt. nsp - 1) go to 2403
c
2803  nsol0 = k
c
      if(ierr .eq. 1) then 
        io_end = 1
      else
        io_err = 1
      endif
c
      nstot = nstot + nsol0
      write(*,'(1x,1a,i10)') 'Number of input solar fluxes = ',nstot
      if(wn0(1) .lt. wn0(k)) then
        if(wneof(1) .lt. 0.0) wneof(1) = wn0(1)
        wneof(2) = wn0(k)
      else
        wneof(1) = wn0(k)
        if(wneof(2) .le. 0.0) wneof(2) = wn0(1)
      endif
c
c****    convert the input solar fluxes to w/m/m/sr/cm-1:
c
      if(iuin .ne. 1) then
c
c****    find the conversoin factor needed to convert solar
c        fluxes from the input units to W/m/m/sr/cm-1
c
        do 3001 k=1,nsol0
c
            call find_units(iuin,iunits,wn0(k),units)
c
            sol0(k) = units*sol0(k)
3001    continue
      endif
c
      if(nsol0 .eq. 0) then
        rewind (iusol2)
        return
      endif
c
      if(wn0(1) .lt. wn0(nsol0)) then
        wneof(1) = wn0(1)
        wneof(2) = wn0(nsol0)
        do 4001 k=1,nsol0
           write(iusol2) wn0(k),sol0(k)
4001    continue
c
c****     read additional values if necessary
c
        if(wn0(nsol0) .lt. wnmax .and. ierr .eq. 0) then
          k = 0
          nsol0 = 0
          go to 2403
        endif
      else
        ius0 = ius0 + 1
        iu = iusol2 + ius0
c
          open(iu,form='unformatted',status='scratch')
c
        do 4021 k=nsol0,1,-1
            write(iu) wn0(k),sol0(k)
4021    continue
c
        if(wn0(nsol0) .gt. wnmin .and. ierr .eq. 0) then
c
c****       get the next spectral segment
c
          k = 0
          nsol0 = 0
          go to 2403
        else
c
c****     put all solar data segments in one file
c
          do 5041 k=ius0,1,-1
              iu = iusol2 + k
              write(*,*) 'readsol: rewinding and reading iu=',iu
              rewind (iu)
5001          read(iu,end=5021) wn0(1),sol0(1)
              write(iusol2) wn0(1),sol0(1)
              if(wneof(1) .le. 0.0d0) wneof(1) = wn0(1)
              go to 5001
5021          close (iu, status='delete')
              write(*,*) 'readsol: closing iu=',iu
5041      continue
          wneof(2) = wn0(1)
        endif
      endif
c
      rewind(iusol2)
c
c****   set file i/o termination status
c
      if(ierr .gt. 0) io_end = 1
      if(ierr .lt. 0) io_err = 1
c
      return
      end
