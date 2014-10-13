      subroutine norm_pi0_ph(nmommx,nlay,nt_pd,small_tau,
     -                       dtauex,dtausc,co_pi0,g,phmom)
c
ccccccccccccccccccccccccc  n o r m _ pi0_ph  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine normalizes the single scattering albedo and     cc
cc    scattering phase function in each atmospheric layer.            cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc       nlay - number of layers in the model atmosphere.             cc
cc     nmommx - maximum number of legendre polynomial moments         cc
cc      nt_pd - number of temperature profiles needed                 cc
cc              1 - radiances only, 2 - partial derivaties            cc
cc  small_tau - very small optical depth (small_tau*small_tau)        cc
cc     dtauex - differential extintion optical depth in each layer    cc
cc     dtausc - differential scattering optical depth in each layer   cc
cc     co_pi0 - single scattering co-albedo in each layer             cc
cc      phmom - scattering phase function moments in each layer.      cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      phmom - scattering phase function moments in each layer.      cc
cc     co_pi0 - single scattering co-albedo in each layer             cc
cc                                                                    cc
ccccccccccccccccccccccccc  n o r m _ pi0_ph  ccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlay,nt_pd,nmommx
      integer mom,k,l
c
      real dtausc(kp),small_tau
c
c*****   monochromatic optical properties
c
      real dtauex(kp,2),co_pi0(kp,2),g(kp),phmom(0:mxmom,kp)
c
c****   normalize the effective asymmetry parameter
c       and moments of phase function.  
c
c       NOTE:  The values of these parameters must be normalized by
c              dtausc rather than dtauex because the rayleigh scattering
c              contribution is explicitly included in these values. 
c
      do 4441 k=1,nlay
          if(dtausc(k) .gt. small_tau) then
            g(k) = g(k)/dtausc(k)
            do 4401 mom=0,nmommx
                phmom(mom,k) = phmom(mom,k)/dtausc(k)
4401        continue
          else
            g(k) = 0.0
            phmom(0,k) = 1.0
            do 4421 mom=1,nmommx
                phmom(mom,k) = 0.0
4421        continue
          endif
4441  continue
c
c****   define the single scattering CO-albedo (1. - pi0) and 
c       normalize the asymmetry parameter.  Note: the asymmetry 
c       parameter is normalized by 
c
      do 4621 l=1,nt_pd
          do 4601 k=1,nlay
              if(dtauex(k,l) .gt. small_tau) then
                co_pi0(k,l) = 1.0 - dtausc(k)/dtauex(k,l)
              else
                co_pi0(k,l) = 1.0
              endif
4601      continue
4621  continue
c
      return
      end
