      subroutine readfil(i,iunit,filez,formz,iform,nzd,ncd,ioff,
     -                   iskip,invz,nval,icol,npt,scz,addz,z,
     -                   ncmax,ierr)
c
cccccccccccccccccccccccccc  r e a d f i l  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads data records from a formatted or          cc
cc    unformatted file.                                               cc
cc                                                                    cc
cc         i - i/o mode -1) skip records, 0) open file, +1) read file cc
cc     iunit - unit number for input                                  cc
cc     formz - format for formatted data files                        cc
cc     filez - name of input file                                     cc
cc       nzd - maximum number of input records                        cc
cc       ncd - maximum number of input columns                        cc
cc         i - data point index                                       cc
cc     iunit - unit number of coordinate grid input                   cc
cc       npt - number of points in z array                            cc
cc     ncmax - number of variables to read in each record             cc
cc      nval - number of variables to keep in each record             cc
cc      icol - column index of each value to keep in each record      cc
cc      invz - data order inversion flag (1) yes (2) no               cc
cc     iform - format mode (formatted, unfor., list directed)         cc
cc       scz - multiplicative scaling factor for each input quantity  cc
cc      addz - additive offset factor for each input quantity         cc
cc      ioff - number of records to skip at top of file               cc
cc     iskip - number of records to skip between records              cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       z(npt,icol) - value of quantity at each point                cc
cc                                                                    cc
cccccccccccccccccccccccccc  r e a d f i l  ccccccccccccccccccccccccccccc
c
      implicit none
c
      integer ncmx
      parameter (ncmx = 80)
c
      character*(40) formz
      character*(132) filez,file
      character*1 fl,name(132)
c
      integer nzd,ncd
      integer icol(ncd),iarray(ncmx),ncmax
      integer nlb,ntb,len,ii,l,j,k,n,iflag,nptj
      integer i,iunit,iform,ioff,iskip,invz,nval,npt,ierr
c
      real array(ncmx),z(nzd,ncd),scz(ncd),addz(ncd)
      double precision darray(ncmx)
c
      ierr = 0
      ncmax = icol(1)
      do 1001 n=1,nval
          if(icol(n) .gt. ncmax) ncmax = icol(n)
1001  continue
      if(ncmax .gt. ncmx) then
        write(*,*) ' Number of input values on each record exceeds ',
     -             'dimension bound in readfil: ',
     -             'ncmax =',ncmax
        ierr = 1
        return
      endif
c
c****   truncate file name if necessary
c
      call charsp(filez,name,len,132,nlb,ntb)
      write(file,'(132a)') (name(ii),ii=1,len)
c
      if(iform .eq. 1) then
c
c****     f o r m a t t e d    d a t a    f i l e
c
        if(i .eq. 0) then
            nptj = 0
            write(*,'(/,2a)') 
     -          ' readfil reading formatted file: ',file
            close(iunit)
            open(iunit,file=file,form='formatted',
     -           status='old',err= 5203)
            iflag = 0
c
c*****        determine if the quantites are real or integer
c
            do 1201 l=1,40
                fl = formz(l:l)
                if(fl .eq. 'i' .or. fl .eq. 'i') then
                  iflag = 1
                endif
1201        continue
        else
          if(i .lt. 0) then
c
c****         skip over records at beginning of file
c
            do 2001 l=1,ioff
              read(iunit,*,end=5010)
2001        continue
          else
c
c****         read in data values and scale.
c
            do 2441 l=1,npt
               j=l
               if(invz.ne.2) j=npt-l+1
               if(iflag .eq. 1) then
                 read(iunit,formz,err=5110,end=5010)
     -                   (iarray(k),k=1,ncmax)
                 do 2201 k=1,nval
                     z(j,k)=  scz(k)*float(iarray(icol(k))) + addz(k)
2201             continue
               else
                 read(iunit,formz,err=5110,end=5010)
     -               (array(k),k=1,ncmax)
                 do 2221 k=1,nval
                     z(j,k)=  scz(k)*array(icol(k)) + addz(k)
2221             continue
               endif
               do 2421 k=1,iskip
                  read(iunit,formz,err=5110,end=5010)
2421           continue
               nptj=nptj+1
2441        continue
          endif
        endif
c
      else
        if(iform .eq. 2) then
c
c****      read an unformatted file
c
          if(i .eq. 0) then
              nptj = 0
              iflag = 2
              write(*,'(2a)') 'formz = ',formz(1:6)
              if(formz(1:7) .eq. 'integer') iflag = 1
              if(formz(1:6) .eq. 'real*8' .or. 
     -             formz(1:6) .eq. 'double') iflag = 3
c
              write(*,'(/,132a)') 
     -             ' readfil reading from unformatted file: ',
     -             (name(ii),ii=1,len)
              close(iunit)
              open(iunit,file=file,form='unformatted',
     -             status='old',err=5203)
c
          else
            if(i .lt. 0) then
c
c****         skip records at top of file.
c
              do 3001 l=1,ioff
                 read(iunit,end=5010)
3001          continue
            else
              do 3441 l=1,npt
                  j=l
                  if(invz.ne.2) j=npt-l+1
c
c****               read the next data point
c
                  if(iflag .eq. 1) then
                    read(iunit,err=5110,end=5010) 
     -                   (iarray(k),k=1,ncmax)
                    do 3101 k=1,nval
                        z(j,k)=  scz(k)*float(iarray(icol(k))) + 
     -                           addz(k)
3101                continue
                  else
                    if(iflag .eq. 2) then
                      read(iunit,err=5110,end=5010) 
     -                    (array(k),k=1,ncmax)
                      do 3201 k=1,nval
                          z(j,k)=  scz(k)*array(icol(k)) + addz(k)
3201                  continue
                    else
                      read(iunit,err=5110,end=5010) 
     -                    (darray(k),k=1,ncmax)
                      do 3221 k=1,nval
                          z(j,k)=  scz(k)*darray(icol(k)) + addz(k)
3221                  continue
                    endif
                  endif
                  do 3421 k=1,iskip
                     read(iunit,end=5010,err=5110)
3421              continue
                  nptj=nptj+1
3441          continue
            endif
          endif
c
        else
c
c****     read a list-directed file
c
          if(i .eq. 0) then
c
c****         open the file
c
            nptj = 0
            write(*,'(/,132a)') 
     -            ' readfil reading from list directed file: ',
     -            (name(ii),ii=1,len)
            close(iunit)
            open(iunit,file=file,form='formatted',
     -           status='old',err=5203)
c
          else
            if(i .lt. 0) then
c
c****         skip records at top of file.
c
              do 4001 l=1,ioff
                read(iunit,*,end=5010)
4001          continue
            else
              do 4441 l=1,npt
                 j=l
                 if(invz.ne.2) j=npt-l+1
                   read(iunit,*,err=5110,end=5010)
     -                  (array(k),k=1,ncmax)
                   do 4201 k=1,nval
                       z(j,k)=  scz(k)*array(icol(k)) + addz(k)
4201               continue
                   do 4421 k=1,iskip
                      read(iunit,*,end=5010,err=5110)
4421               continue
                   nptj = nptj + 1
4441             continue
            endif
          endif
        endif
      endif
c
      return
c
c****   end of file encountered
c
5010  ierr = 1
      npt = nptj
      close(iunit)
c
      return
c
c****    error encountered in file
c
5110  ierr = -1
      close(iunit)
c
      return
c
c****  could not open file
c
5203   ierr = -2
      write(*,'(/,1a,i5,132a)') 
     - 'I/0 Error - could not open unit ',iunit,' for file: ',
     -            (name(ii),ii=1,len)
      close(iunit)
c
      return
      end
