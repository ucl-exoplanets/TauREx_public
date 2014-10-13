      program gaus_int
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                     cc
cc    p u r p o s e :                                                  cc
cc                                                                     cc
cc    this program uses gaussian quadrature schemes to evaulate        cc
cc    integrals of sample functions.                                   cc
cc                                                                     cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      parameter (maxmu = 16)
c
      real gmu(maxmu),gwt(maxmu)
c
      pi = 3.14159
c
      write(*,*) 'enter number of gaussian points: '
      read(*,*) nstr
c
c****   call standard gaussian quadrature routine:
c
      call qgausn(nstr,gmu,gwt)
c
      write(*,'(/,/,1a,/,1a)') ' qgausn output ','  n     gmu      gwt'
      gsum1 = 0.0
      do 1001 n=1,nstr
c          gmu(nstr-n+1) = -gmu(n)
c          gwt(nstr-n+1) = gwt(n)
          write(*,'(i5,2e12.4)') n,gmu(n),gwt(n)
          gsum1 = gsum1 + gwt(n)*cos(pi*gmu(n))
1001  continue
c
      write(*,'(/,/,1a,1pe12.4)') 'gsum1 =',gsum1
c
      stop
      end


