      function cp_atm(ncomp,icomp,volmix,t)
c
cccccccccccccccccccccccccccccc  c p _ a t m  ccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the specific heat at constant pressure cc
cc    for a specified temperature and atmospheric temperature.        cc
cc    specific heats for specific gases are determined by evaluating  cc
cc    a cubic polynomial of the form:                                 cc
cc                                                                    cc
cc                cp(t) = cp1*t**3 + cp2*t**2 + cp3*t + cp4           cc
cc                                                                    cc
cc    These coefficients are derived from cp/R vs T data (Hilsenrath) cc
cc    with the aid of the program: /home/dc/util/fitcp.f              cc
cc                                                                    cc
cc    The net atmospheric specific heat of a mixture of gases is      cc
cc    then determined by taking the volume-mixing-ratio-weighted      cc
cc    average of the specific heats of all components.                cc
cc                                                                    cc
cc    r e f e r e n c e s :                                           cc
cc                                                                    cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    ncomp = number of major atmosphere constituents                 cc
cc    icomp = atmosphere constituent index                            cc
cc            (1) air  (2) co2  (3) n2  (4) o2  (5) h2  (6) he        cc
cc   volmix = volume mixing ratio of each constituent                 cc
cc        t = temperature at which cp is needed (K)                   cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc       cp = specific heat at constant pressure at temperature t.    cc
cc                                                                    cc
cccccccccccccccccccccccccccccc  c p _ a t m  ccccccccccccccccccccccccccc
c
      implicit none
c
      integer ncomp,icomp(ncomp)
c
      integer i
c
      real volmix(ncomp),t,cp_atm
c
      real cp1(2),cp2(2),cp3(2),cp4(2),vmix,cp
c
      data cp1 /-2.629e-7,3.47e-7/
      data cp2 /6.059e-4,-1.269e-3/
      data cp3 /-2.232e-1,1.688e+0/
      data cp4 /1.024e+3,4.431e+2/
c
      vmix = 0.
      cp = 0.
c
      do 1001 i=1,ncomp
          vmix = vmix + volmix(i)
          if(icomp(i) .le. 2) then
            cp = cp + volmix(i)*(cp1(icomp(i))*t*t*t + 
     -                cp2(icomp(i))*t*t + 
     -                       cp3(icomp(i))*t + cp4(icomp(i)))
          endif
c
c****          use temperature-independent values
c
          if(icomp(i) .eq. 3) then
c
c****          the following N2 value comes from the CRC 
c
              cp = cp + 1059.
c
          endif
c
          if(icomp(i) .eq. 4) then
c
c****            the following O2 value comes from the CRC pg D165
c
            cp = cp + 916.29
c
          endif
c
          if(icomp(i) .eq. 5) then
c
c****        the following value comes from the hilsenrath tables
c            for H2 at 80K
c 
            cp = cp + 10933.
c
          endif
c
          if(icomp(i) .eq. 6) then
c
c****               the following He value comes from the CRC pg D165
c
                  cp = cp + 5196.0
c
          endif
1001  continue
c
      cp_atm = cp/vmix
c
      return
      end
