      subroutine map_spect(nlyr,nlay,nstr,nmom,
     -                     levout,levtau1,nref,nza,
     -                     wn,delnu,dtauex,co_pi0,g,
     -                     phmom,surf_opt,alb,taumn,
     -                     tauerr,pi0err,phferr,surferr,
     -                     dnugrp,wngrp,surfgrp,albgrp,
     -                     taugrp,pi0grp,ggrp,pmomgrp,
     -                     nmomgrp,iwngrp,ngroup,ismterr,
     -                     itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt)
c
cccccccccccccccccccccccccc  m a p _ s p e c t  ccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine takes monochromatic optical properties and      cc
cc    uses a spectral binning technique to identify monochromatic     cc
cc    spectral regions with similar absorption properties at all      cc
cc    levels of the atmosphere.  These intervals are combined into a  cc
cc    smaller number of quasi-monochrmatic bins.                      cc
cc                                                                    cc
cc    ***  this version of the program uses an exponetial-sum adding  cc
cc         method to combine simlilar optical depths.                 cc
cc         it also uses an error criterion that depends on the        cc
cc         layer absorption optical depth - layers with optical       cc
cc         depths less than unity have the most stringent constraint. cc
cc         Other values are weighted as:                              cc
cc            tauerr =  err*(1 + 0.1dtau + 0.001*dtau**2)             cc
cc								      cc
cc    note:  this version of the program uses surface optical         cc
cc           property gradients.                                      cc
cc    note:  this version of the program uses optcial depth           cc
cc           gradients and finds corresponding values of pi0grp.      cc
cc                                                                    cc
cc    i n p u t    p a r a m e t e r s :                              cc
cc        wn : wavenumber of input monochromatic interval.            cc
cc      nlyr : number of layers in which coefficients are found       cc
cc      nlay : number of layers were optical depths are needed:       cc
cc              nlyr if pressure is not a varible part of the state   cc
cc              vector, nlev if it is (istate(1) = 1)                 cc
cc      nmom : number of legendre polynomial coefficients for the     cc
cc              scattering phase function expansion.                  cc
cc       nref - number of surface optical properties specified at     cc
cc              each wavelength.                                      cc
cc    tauerr : fractional dtau error allowed in interval bins         cc
cc    pi0err : fractional pi0 error allowed in interval bins          cc
cc    phferr : maximum fractional asymmetry parameter error allowed   cc
cc   surferr : fractional surface optical property error allowed      cc
cc             in interval bins                                       cc
cc     delnu : width of monochromatic spectral interval               cc
cc      dtau : extinction optical depth in spectral interval          cc
cc       pi0 : single scattering albedo in spectral interval          cc
cc                                                                    cc
cc    o u t p u t    v a r i a b l e s :                              cc
cc                                                                    cc
cc     taugrp - array of total optical depth bins at model            cc
cc            levels for each temperature profile                     cc
cc    ismterr - error flag: ngroup > ngrp                             cc
cc     ngroup - number of spectral bins (0 to ngrp)                   cc
cc     iwngrp - index of spectral bin                                 cc
cc   itau1cnt - number of bins rejected by poor fit at tau=1          cc
cc    itaucnt - number of bins rejected by poor fit away from tau=1   cc
cc    ipi0cnt - number of bins rejected by poor pi0 fit               cc
cc      igcnt - number of bins rejected by poor g fit                 cc
cc    isurcnt - number of bins rejected by poor albedo fit            cc
cc    levtau1 - index of tau=1 level                                  cc
cc    nmomgrp - maximum number of phase function moments in bin       cc
cc     taugrp - mean, min, and max optical depth in each bin          cc
cc     piogrp - mean, min, and max single scattering albedo in  bin   cc
cc    pmomgrp - mean, min, and max phase function in each bin         cc
cc       ggrp - mean, min, and max asymmetry parmeter in each bin     cc
cc                                                                    cc
cccccccccccccccccccccccccc  m a p _ s p e c t  ccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlyr,nstr,nmom,nref,nlay,levout,ismterr,nza
      integer itau1cnt,itaucnt,ipi0cnt,igcnt,isurcnt
c
c****   index of optical depth unity for each bin
c
      integer levtau1(ngrp,2)
c
c****   spectral binning parameters
c
      integer nmomgrp(ngrp),iwngrp,ngroup
c
c****    internal variables
c
      integer i,j,nbin,m,nr,k,mom,nz
c
      real taugrp(kp,ngrp,3),pi0grp(kp,ngrp,3),
     -     ggrp(kp,ngrp,3),pmomgrp(0:mxmom,kp,ngrp,3),
     -     surfgrp(ngrp,4,3),albgrp(ngrp,nsol,3)
c
      real taumn,tauerr,pi0err,phferr,surferr
c
c       monochromatic optical properties
c
      real dtauex(kp,2),co_pi0(kp,2),g(kp),phmom(0:mxmom,kp),
     -     surf_opt(4),alb(nsol)
c
c****   internal variables
c
      real albmin0,albmax0,smallpi0,smallg,smallalb,errtst,abserr,
     -     dtmn,dtmx,pi0mn,pi0mx,ph_mn,ph_mx,tauabs,tauerror,
     -     taumin0,taumax0,pi0min0,pi0max0,gmax0,gmin0,gmn,gmx,
     -     ph_min0,ph_max0
c
c****   variables defining bin limits
c
      real tautot(kp,2),dtaumn(kp,ngrp),tau_rng(kp,ngrp),
     -        pi0_rng(kp,ngrp),g_rng(kp,ngrp),ph_rng(mxmom,kp,ngrp),
     -        alb_rng(ngrp,4)
      save tautot,dtaumn,tau_rng,pi0_rng,g_rng,ph_rng,alb_rng
c
      double precision wn,delnu,wngrp(ngrp,3),dnugrp(ngrp)
      double precision dng0,dng1
c
c****    define a small number as a minimum tau and pi0 error
c
      smallpi0 = 5.e-3*pi0err*pi0err + 1.0e-4
      smallg = 5.e-3*phferr*phferr + 1.0e-4
      smallalb = 5.e-3*surferr*surferr + 1.e-4
c
c****   define the reference level for integrated absorption optical
c       depths.  This reference level will be the first level tested
c       by the binning scheme.  If radiances are needed only at the top 
c       of the atmosphere, the optical depth reference is the top 
c       of the atmosphere. If radiances are needed at the surface, 
c       the reference level is the surface.
c
      j = 1
      if(levout .eq. 2) j = 2
c
c*****  compare this monochromatic segment to each smt bin
c
      nbin = 0
      errtst = 1.e9
      do 1712 m=ngroup,1,-1
          abserr = 0.0
c
c****         check errors at tau = 1
c
          dtmn = taugrp(levtau1(m,j),m,3) - tau_rng(levtau1(m,j),m)
          if(dtmn .lt. dtaumn(levtau1(m,j),m)) 
     -       dtmn = dtaumn(levtau1(m,j),m)
          if(dtauex(levtau1(m,j),1) .lt. dtmn) then
            itau1cnt = itau1cnt + 1
            go to 1712
          else
            dtmx = taugrp(levtau1(m,j),m,2) + tau_rng(levtau1(m,j),m)
            if(dtauex(levtau1(m,j),1) .gt. dtmx) then
              itau1cnt = itau1cnt + 1
              go to 1712
            endif
          endif
          abserr = abserr + abs(dtauex(levtau1(m,j),1) - 
     -             taugrp(levtau1(m,j),m,1))/tau_rng(levtau1(m,j),m)
c
c****       check surface albedo match
c              
          do 1001 nz=1,nza
              if(alb(nz) .lt. 
     -          (albgrp(m,nz,3) - alb_rng(m,nz))) then
                isurcnt = isurcnt + 1
                go to 1712
              else
                if(alb(nz) .gt. 
     -            (albgrp(m,nz,2) + alb_rng(m,nz))) then
                  isurcnt = isurcnt + 1
                  go to 1712
                endif
              endif
1001      continue
c
c****       check optical depth values in all layers below tau = 1
c
          do 1201 k=levtau1(m,j)+1,nlyr
c
c****           define the optical depth limits of this bin
c
              dtmn = taugrp(k,m,3) - tau_rng(k,m)
              if(dtmn .lt. dtaumn(k,m)) dtmn = dtaumn(k,m)
              dtmx = taugrp(k,m,2) + tau_rng(k,m)
              if(dtauex(k,1) .lt. dtmn .or. dtauex(k,1) .gt. dtmx) then
c
c****               this value does not belong in this bin.
c
                itaucnt = itaucnt + 1
                go to 1712
              endif
              abserr = abserr + abs(dtauex(k,1) - taugrp(k,m,1))/
     -                 tau_rng(k,m)
1201      continue
c
c****       check optical depth values in all layers above tau = 1
c
          do 1221 k=levtau1(m,j)-1,1,-1
c
c****           define the optical depth limits of this bin
c
              dtmn = taugrp(k,m,3) - tau_rng(k,m)
              if(dtmn .lt. dtaumn(k,m)) dtmn = dtaumn(k,m)
              dtmx = taugrp(k,m,2) + tau_rng(k,m)
              if(dtauex(k,1) .lt. dtmn .or. dtauex(k,1) .gt. dtmx) then
c
c****               this value does not belong in this bin.
c
                itaucnt = itaucnt + 1
                go to 1712
              endif
              abserr = abserr + abs(dtauex(k,1) - taugrp(k,m,1))/
     -                 tau_rng(k,m)
1221      continue
c
c****        check single scattering co-albedos in layers below tau = 1
c
          do 1401 k=levtau1(m,j),nlyr
c
c****            if the single scattering co-albedo is larger than 
c                a minimum value -scattering is important
c
              if((1.0 - co_pi0(k,1))*dtauex(k,1) .gt. taumn) then
c
c****             the scattering optical depth is finite.  
c                 define the single-scattering co-albedo limits
c
                pi0mn = pi0grp(k,m,3) - pi0_rng(k,m)
                pi0mx = pi0grp(k,m,2) + pi0_rng(k,m)
                if(co_pi0(k,1) .lt. pi0mn .or. 
     -             co_pi0(k,1) .gt. pi0mx) then
c
c****             this value does not belong in this bin.
c
                  ipi0cnt = ipi0cnt + 1      
                  go to 1712
                endif
                abserr = abserr + abs(co_pi0(k,1) -
     -                  pi0grp(k,m,1))/pi0_rng(k,m)
             endif
1401      continue
c
c****        check single scattering co-albedos in layers above tau = 1
c
          do 1421 k=levtau1(m,j)-1,1,-1
c
c****            if the single scattering co-albedo is larger than 
c                a minimum value -scattering is important
c
              if((1.0 - co_pi0(k,1))*dtauex(k,1) .gt. taumn) then
c
c****             the scattering optical depth is finite.  
c                 define the single-scattering co-albedo limits
c
                pi0mn = pi0grp(k,m,3) - pi0_rng(k,m)
                pi0mx = pi0grp(k,m,2) + pi0_rng(k,m)
                if(co_pi0(k,1) .lt. pi0mn .or. 
     -             co_pi0(k,1) .gt. pi0mx) then
c
c****             this value does not belong in this bin.
c
                  ipi0cnt = ipi0cnt + 1      
                  go to 1712
                endif
                abserr = abserr + abs(co_pi0(k,1) -
     -                  pi0grp(k,m,1))/pi0_rng(k,m)
              endif
1421      continue
c
c****        check scattering phase function in layers below tau = 1
c
          do 1521 k=levtau1(m,j),nlyr
c
c*****          determine if scattering is important
c
              if((1.0 - co_pi0(k,1))*dtauex(k,1) .gt. taumn) then
c
                gmn = ggrp(k,m,3) - g_rng(k,m)
                gmx = ggrp(k,m,2) + g_rng(k,m)
                if(g(k) .lt. gmn .or. g(k) .gt. gmx) then
c
c****             this value does not belong in this bin.
c
                  igcnt = igcnt + 1      
                  go to 1712
                endif
c
c****                check each term in the phase function expansion
c
                do 1501 mom=1,nstr
                    ph_mn = pmomgrp(mom,k,m,3) - ph_rng(mom,k,m)
                    ph_mx = pmomgrp(mom,k,m,2) + ph_rng(mom,k,m)
                    if(phmom(mom,k) .lt. ph_mn .or. 
     -                  phmom(mom,k) .gt. ph_mx) then
c
c****                       this value does not belong in this bin.
c
                      igcnt = igcnt + 1
                      go to 1712
                    endif
1501            continue
c
              endif
1521      continue
c
c****        check scattering phase function in layers above tau = 1
c
          do 1621 k=levtau1(m,j)-1,1,-1
c
c*****          determine if scattering is important
c
              if((1.0 - co_pi0(k,1))*dtauex(k,1) .gt. taumn) then
c
                gmn = ggrp(k,m,3) - g_rng(k,m)
                gmx = ggrp(k,m,2) + g_rng(k,m)
                if(g(k) .lt. gmn .or. g(k) .gt. gmx) then
c
c****             this value does not belong in this bin.
c
                  igcnt = igcnt + 1      
                  go to 1712
                endif
c
c****                check each term in the phase function expansion
c
                do 1601 mom=1,nstr-1
                    ph_mn = pmomgrp(mom,k,m,3) - ph_rng(mom,k,m)
                    ph_mx = pmomgrp(mom,k,m,2) + ph_rng(mom,k,m)
                    if(phmom(mom,k) .lt. ph_mn .or. 
     -                  phmom(mom,k) .gt. ph_mx) then
c
c****                       this value does not belong in this bin.
c
                      igcnt = igcnt + 1
                      go to 1712
                    endif
1601            continue
c
              endif
1621      continue
c
c****       congratulations! you passed all of the tests and you
c           qualify to enter this bin
c
          if(abserr .lt. errtst) then
            errtst = abserr
            nbin = m
          endif
1712  continue
c
      if(nbin .ne. 0) then
        m = nbin
c
c****       congratulations! you passed all of the tests and you
c           qualify to enter this bin
c
          iwngrp = m
          dng0 = dnugrp(m) + 1.0d-7
          dnugrp(m) = dnugrp(m) + delnu
          dng1 = dnugrp(m) + 1.0d-7
          wngrp(m,1) = (dng0*wngrp(m,1) + delnu*wn)/dng1
          wngrp(m,3) = wn
c
c****       modify min and max albedo values
c
          do 2001 nz=1,nza
              albgrp(m,nz,1) = real((dng0*albgrp(m,nz,1) + 
     -                          delnu*alb(nz))/dng1)
              if(alb(nz) .lt. albgrp(m,nz,2)) 
     -            albgrp(m,nz,2) = alb(nz)
              if(alb(nz) .gt. albgrp(m,nz,3)) 
     -            albgrp(m,nz,3) = alb(nz)
2001      continue
          do 2021 nr=1,nref
              surfgrp(m,nr,1) = real((dng0*surfgrp(m,nr,1) + 
     -                          delnu*surf_opt(nr))/dng1)
              if(surf_opt(nr) .lt. surfgrp(m,nr,2))
     -            surfgrp(m,nr,2) = surf_opt(nr)
              if(surf_opt(nr) .gt. surfgrp(m,nr,3))
     -            surfgrp(m,nr,3) = surf_opt(nr)
2021      continue
c
          do 2201 k=1,nlay
c
              taugrp(k,m,1) = real((dng0*taugrp(k,m,1) + delnu*
     -                           dtauex(k,1))/dng1)
c
c****           update min and max optical depth values
c 
              if(dtauex(k,1) .lt. taugrp(k,m,2)) 
     -           taugrp(k,m,2) = dtauex(k,1)
c
c****           load maximum optical depth and appropriate 
c               single scattering abledo
c
              if(dtauex(k,1) .gt. taugrp(k,m,3))  
     -           taugrp(k,m,3) = dtauex(k,1)
c
c               use a running average to add co-pi0 contribution.
c
              pi0grp(k,m,1) = real((dng0*pi0grp(k,m,1) + delnu*
     -                       co_pi0(k,1))/dng1)
c
c****           update min and max co-albedo values
c 
              if(co_pi0(k,1) .lt. pi0grp(k,m,2)) 
     -           pi0grp(k,m,2) = co_pi0(k,1)
              if(co_pi0(k,1) .gt. pi0grp(k,m,3)) 
     -           pi0grp(k,m,3) = co_pi0(k,1)
c
c****            update the asymmetry factor
c
              ggrp(k,m,1) = real((dng0*ggrp(k,m,1) + delnu*
     -                       g(k))/dng1)
c
c****            update moments of phase function
c
              do 2101 mom=0,nmom
                  pmomgrp(mom,k,m,1) = real((dng0*pmomgrp(mom,k,m,1) + 
     -                                  delnu*phmom(mom,k))/dng1)
2101          continue 
c
c****           update min and max asymmetry factor values
c               NOTE: the phase function moments are slaved to the
c               variations in the asymmetry parameter
c 
              if(g(k) .lt. ggrp(k,m,2)) then
                ggrp(k,m,2) = g(k)
                do 2121 mom=1,nmom
                    pmomgrp(mom,k,m,2) = phmom(mom,k)
2121            continue
              endif
              if(g(k) .gt. ggrp(k,m,3)) then
                ggrp(k,m,3) = g(k)
                do 2141 mom=1,nmom
                    pmomgrp(mom,k,m,3) = phmom(mom,k)
2141            continue
              endif
c
2201      continue
c
c****       check to see if the number of phase function moments
c           must be updated.
c
          if(nmom .gt. nmomgrp(m)) nmomgrp(m) = nmom
c
c****       you are done. go to the next monochromatic spectral interval
c
          return
c
      endif
c
c****   sorry, you're a real misfit. start a bin of your own.
c
c****   determine if number of bins exceeds dimension bound.
c
      if(ngroup .eq. ngrp) then
c
c****     no more groups will fit - return
c
        ismterr = 1
        return
c
      endif
c
c****  increment the group counter
c
      ngroup = ngroup + 1
c
c****   initialize bin properties
c
      m = ngroup
      iwngrp = ngroup
      nmomgrp(m) = nmom  
      dnugrp(m) = delnu
      tautot(1,1) = 1.e-7
      tautot(nlyr+1,2) = 1.e-7
      do 3041 i=1,3
          wngrp(m,i) = wn
          do 3001 nr=1,nref
              surfgrp(m,nr,i) = surf_opt(nr)
3001      continue
          do 3021 nz=1,nza
              albgrp(m,nz,i) = alb(nz)
3021      continue
3041  continue
c
      do 3201 k=1,nlay
c
c****       define the column-integrated absorption optical depth
c
          tautot(k+1,1) = tautot(k,1) + co_pi0(k,1)*dtauex(k,1)
          tautot(nlay-k+1,2) = tautot(nlay-k+2,2) + co_pi0(nlay-k+1,1)*
     -                       dtauex(nlay-k+1,1)
3201  continue
c
      do 3421 k=1,nlay
c
c****        initialize optical depth, single scattering co-albedo,
c            assymetry parameter and phase function moments for bin
c
          taugrp(k,m,1) = dtauex(k,1)
c
          pi0grp(k,m,1) = co_pi0(k,1)
c
          ggrp(k,m,1) = g(k)
c
          do 3401 mom=0,nmom
              pmomgrp(mom,k,m,1) = phmom(mom,k)
3401      continue
3421  continue
c
c****   find the index of the tau(abs) =1 level for this bin
c
      levtau1(m,1) = nlyr
      levtau1(m,2) = 1
      do 3601 k=1,nlyr
          if(levtau1(m,1) .eq. nlyr .and. tautot(k+1,1) .ge. 1.) 
     -           levtau1(m,1) = k
          if(levtau1(m,2) .eq. 1 .and. tautot(nlyr-k+1,2) .ge. 1.) 
     -           levtau1(m,2) = nlyr-k+1
3601  continue
c
c*****  define the optical depth range for this bin.
c       the criteria is much weaker for large absorption 
c       optical depths.
c
      do 4021 k=1,nlyr
c
c****        define the absorption optical depth
c
          tauabs = pi0grp(k,m,1)*taugrp(k,m,1) 
c
c****       define the layer-dependent optical depth error for bin.
c           the following function varies from ~tauerr at tau < 1,
c           to a much larger value at large layer absorption optical 
c           depths (ie for completely opaque layers).
c 
          tauerror = taugrp(k,m,1)*tauerr*(1.0 + 0.1*tauabs + 
     -                    0.01*tauabs*tauabs) + 10.*taumn
c
c****       define the nominal optical depth limits
c
          taumax0 = taugrp(k,m,1) + tauerror
          taumin0 = taugrp(k,m,1)*(1.0 - tauerr) - taumn
          tau_rng(k,m) = taumax0 - taumin0
          taugrp(k,m,3) = taugrp(k,m,1)
          taugrp(k,m,2) = taugrp(k,m,1)
          dtaumn(k,m) = taugrp(k,m,1)*(1.0 - tauerr) - taumn
c
c****       define the single scattering co-albedo limits
c
          pi0max0 = pi0grp(k,m,1)*(1.0 + pi0err) + smallpi0
          pi0min0 = pi0grp(k,m,1)*(1.0 - pi0err) - smallpi0
          pi0_rng(k,m) = pi0max0 - pi0min0
          pi0grp(k,m,3) = pi0grp(k,m,1)
          pi0grp(k,m,2) = pi0grp(k,m,1)
c
c****       define the asymmetry parameter limits
c
          gmax0 = ggrp(k,m,1)*(1.0 + phferr) + smallg
          gmin0 = ggrp(k,m,1)*(1.0 - phferr) - smallg
          g_rng(k,m) = gmax0 - gmin0
          ggrp(k,m,3) = ggrp(k,m,1)
          ggrp(k,m,2) = ggrp(k,m,1)
c
c****       define the phase function limits.
c
          pmomgrp(0,k,m,2) = 1.0
          pmomgrp(0,k,m,3) = 1.0 
          do 4001 mom=1,nmom
              ph_max0 = pmomgrp(mom,k,m,1)*(1.0 + phferr) + smallg
              ph_min0 = pmomgrp(mom,k,m,1)*(1.0 - phferr) - smallg
              ph_rng(mom,k,m) = ph_max0 - ph_min0
              pmomgrp(mom,k,m,3) = pmomgrp(mom,k,m,1)
              pmomgrp(mom,k,m,2) = pmomgrp(mom,k,m,1)
4001      continue
c
4021  continue
c
c****     define surface optical property limits
c
      do 4041 nr=1,nref
          albmax0 = albgrp(m,nr,1)*(1.0 + surferr) + smallalb
          albmin0 = albgrp(m,nr,1)*(1.0 - surferr) - smallalb
          alb_rng(m,nr) = albmax0 - albmin0
4041  continue
c
c****  get the next spectral interval
c
      return
      end
