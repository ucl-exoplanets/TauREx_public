      subroutine readpab(iu,wnmin,wnmax,npab0,ntpab,
     -                   tpab0,wnpab0,pab0,
     -                   ioend,ioerr,wneof)
c
ccccccccccccccccccccccccccc  r e a d p a b  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine reads wavelength-dependent pressure-induced     cc
cc    absorption coefficients of gases at one or more temperatures    cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        iu - unit number for absorption coefficeints.               cc
cc     wnpab - wavenumbers where absorption coefficients are given    cc
cc       pab - absorption coefficient at each wavenumber              cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      npab0 - number of wavenumbers where absorption coefficients   cc
cc             are specified                                          cc
cc     wnpab - wavenumbers where absorption coefficients are given    cc
cc       pab - absorption coefficient at each wavenumber              cc
cc                                                                    cc
ccccccccccccccccccccccccccc  r e a d p a b  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nsp
      parameter (nsp=10000)
c
      character*255 str
      character*1 str1(255)
c
      integer iu,npab0,ntpab,ioerr,ioend
      integer ic,len,nlb,ntb,itst,i,n
c      
      double precision wnmin,wnmax,wneof(2)
      real pab0(kp,nsp),tpab0(kp)
      double precision wn0,wnpab0(nsp)
c
c****    read the number of temperatures
c
      ic=0
1021  if(ic .lt. 100) then
        ic = ic + 1
c
c****     skip over any ascii header at the top of the file, and
c         find the first record with numbers in it
c
        read(iu,'(1a)') str
c
c*****   find the first non-blank character
c
        call charsp(str,str1,len,255,nlb,ntb)
c
        itst = ichar(str1(1))
        if(itst .ge. 48 .and. itst .le. 57 .and. len .gt. 0) then 
c
c         this record includes numerical values
c             
          read(str,*,err=1021) ntpab,(tpab0(i),i=1,ntpab)
          read(iu,'(/,/)')
        else
          go to 1021
        endif
      else
        write(*,'(/,1a,i5,1a)') 
     -   ' Error reading unit',iu,' in readpab'
        stop
      endif
c
      n = 0
2001  read(iu,*,end=2021,err=2041) wn0,(pab0(i,n+1),i=1,ntpab)
c
      wnpab0(n+1) = wn0
      if(wnpab0(n+1) .ge. wnmin) n = n + 1
      if(n .gt. nsp) go to 2061
c
      if(wn0 .le. wnmax) go to 2001
c
2021  close(iu)
      npab0 = n
c
      ioend = 1
      ioerr = 0
      wneof(1) = wn0
      write(*,'(/,1a,1pe14.6)') 'EOF in readpab after wn = ',wnpab0(n)
      return
c
2041  npab0 = n
      ioend = 0
      ioerr = 1
      wneof(1) = wn0
      write(*,'(/,1a,1pe14.6)') 'Error in readpab after wn = ',wnpab0(n)
      return
c
2061  npab0 = n
      ioend = 0
      ioerr = 1
      wneof(1) = wn0
      write(*,'(/,1a,1pe14.6,/2a,i5)') 
     -    'Error in readpab after wn = ',wnpab0(n),
     -    'Number of elements in array exceeds dimension bound:',
     -    ' nsp = ',nsp
c
      return
      end
