      subroutine delted(nlay,nza,alb,u0,solflx,dtau,om,g,
     -                  flxu,flxd,dir)
c
ccccccccccccccccccccccccccc  d e l t e d  cccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the upward, downward and net solar     cc
cc    fluxes in an inhomogeneous absorbing, scattering atmosphere     cc
cc    the delta-eddington approximation is used to find the diffuse   cc
cc    reflectivity and transmissivity and the total upward and        cc
cc    downward fluxes for each of the nlay homogeneous layers.        cc
cc    the adding method is then used to combine these layers.  if     cc
cc    any layer is thicker than dtau = *emax1*, it is assumed to be   cc
cc    semi-infinite.  layers thicker than dtau = *emax2*, are treated cc
cc    as infinite layers.                                             cc
cc                                                                    cc
cc    note: to account for a diffuse flux at the top of the atmos-    cc
cc          phere, the user must initialize flxd(1,m) to the appro-   cc
cc          priate value at all m zenith angles.  if there is no      cc
cc          downward diffuse flux at the top of the atmosphere, the   cc
cc          user must initialize flxd(1,m) to zero in the calling     cc
cc          program.                                                  cc
cc                                                                    cc
cc      nlay - number of homogeneous model layers.                    cc
cc       nza - number of zenith angles                                cc
cc       alb - surface albedo                                         cc
cc        u0 - array of cosines solar zenith angles                   cc
cc    solflx - normal incidence solar flux at top of atmosphere.      cc
cc      dtau - array of normal incidence optical depths in each       cc
cc             homogeneous model layer.                               cc
cc        om - array of single scattering albedos for each homo-      cc
cc             geneous model layer.                                   cc
cc         g - array of assymetry parameters for each homogeneous     cc
cc             model layer.                                           cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc      flxu - upward flux at nlay+1 layer boundries.                 cc
cc             (flxu(l) refers to the upward flux at the top          cc
cc              of layer l)                                           cc
cc      flxd - downward flux at nlay+1 layer boundries.               cc
cc             (flxd(l) refers to the downward flux at the bottom     cc
cc              of layer l-1)                                         cc
cc       dir - direct solar flux at the top of each layer and at      cc
cc             surface.                                               cc
cc                                                                    cc
ccccccccccccccccccccccccccc  d e l t e d  cccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      parameter (kp = 75, nsol =10)
c
c********** common blocks used in delta-eddington routines.
c
      dimension dtau(kp),tau(kp),om(kp),g(kp),u0(nsol)
c
      dimension dir(kp,nsol),flxu(kp,nsol),flxd(kp,nsol)
c
      dimension gp(kp),omp(kp),dtaup(kp),ru(kp),dd(kp),dflux(kp)
c
      dimension g1(kp),g2(kp),g3(kp),g4(kp),
     *          ufl(kp),dfl(kp),uflux(kp),gsq(kp),ekt(kp),
     *          etu(kp),da(kp),db(kp),dc(kp),du(kp),
     *          sk(kp),rl(kp),tl(kp),rs(kp),a1(kp),
     *          a2(kp),tdmu(kp)
c
      data prec/1.e-10/
c
      data emax1,emax2/8.0d0,70.0d0/
c
      nlev = nlay + 1
      np2 = nlay+ 2
      nlm = nlay- 1
      vlarge =  dexp(emax2)
      tau(1) = 0.
c
c****   scale the optical depths, single scattering albedos and the
c       scattering assymetry factors for use in the delta-eddington
c       approximation.  initialize other quantities.  use wiscombe's
c       trick to subtract a small value from the single scattering
c       for the case of a conservative atmosphere
c
      do 1081 l=1,nlay
          gsq(l) = g(l)*g(l)
          gp(l) = g(l)/(1. + g(l))
          if(1. - om(l) .lt. prec) om(l) = 1. - prec
          dtaup(l) = (1. - om(l)*gsq(l))*dtau(l)
          tau(l+1) = tau(l) + dtaup(l)
          omp(l) = (1. - gsq(l))*om(l)/(1. - om(l)*gsq(l))
          g1(l) = .25*(7. - omp(l)*(4. + 3.*gp(l)))
          g2(l) = -.25*(1. - omp(l)*(4. - 3.*gp(l)))
          sk(l) = dsqrt(g1(l)*g1(l) - g2(l)*g2(l))
1081  continue
c
c**** diffuse transmittance and reflectance of each homogeneous layer
c
c     the equations for rl and tl were derived by solving the homogeneous
c     part of the equation of transfer (ie. no source term)
c
      do 1163 l=1,nlay
          skt = sk(l)*dtaup(l)
          if(skt .gt.emax2) go to 1141
          ekt(l) = dexp(skt)
          if(skt .gt. emax1) go to 1121
          e2ktm = ekt(l)*ekt(l) - 1.
          denom = g1(l)*e2ktm + sk(l)*(e2ktm + 2.)
          rl(l) = g2(l)*e2ktm/denom
          tl(l) = 2.*sk(l)*ekt(l)/denom
          go to 1163
c
c     semi-infinite layers
c
1121      rl(l) = g2(l)/(g1(l) + sk(l))
          tl(l) = 2.*sk(l)/(ekt(l)*(g1(l) + sk(l)))
          go to 1163
c
c     infintie layers
c
1141      rl(l) = g2(l)/(g1(l) + sk(l))
          tl(l) = 0.
          ekt(l) = vlarge
1163  continue
c
c****   set the "reflectivity" and "transmissivity" of the surface.
c
      rl(nlev) = alb
      tl(nlev) = 0.
c
c****   use adding method to find the reflectance and transmittance
c       of combined layers.  add downward from the top and upward
c       from the bottom at the same time.
c
      rs(1) = rl(1)
      ru(nlev) =alb
      do 1221 l=1,nlay
          dd(l) = 1./(1. - rs(l)*rl(l+1))
          rs(l+1) = rl(l+1) + tl(l+1)*tl(l+1)*rs(l)*dd(l)
          du(nlev-l) = 1./(1. - rl(nlev-l)*ru(np2-l))
          ru(nlev-l) = rl(nlev-l) + tl(nlev-l)*tl(nlev-l)*ru(np2-l)*
     *                 du(nlev-l)
1221  continue
c
c****    z e n i t h    a n g l e    l o o p
c
      do 1782 m=1,nza
c
c****       initialize the direct solar flux at the top of each layer
c           and other quantities needed to find the diffuse fluxes.
c
          dir(1,m) = u0(m)*solflx
          do 1502 l=1,nlay
              g3(l) = .25*(2. - 3.*gp(l)*u0(m))
              g4(l) = 1. - g3(l)
              db(l) = 1. + sk(l)*u0(m)
              a1(l) = g1(l)*g4(l) + g2(l)*g3(l)
              a2(l) = g1(l)*g3(l) + g2(l)*g4(l)
              dc(l) = 2. - db(l)
              ufl(l) = 0.
              dfl(l) = 0.
              dir(l+1,m) = 0.
              tdmu(l) = dtaup(l)/u0(m)
              etu(l) = 1.0/vlarge
              if(tdmu(l) .le. emax2) etu(l) = dexp(-tdmu(l))
              da(l) = db(l)*dc(l)*((sk(l) + g1(l))*ekt(l) +
     *                (sk(l) - g1(l))/ekt(l))
1502      continue
          do 1521 l=1,nlev
              tdmu(l) = tau(l)/u0(m)
1521      continue
1551      do 1561 l=2,nlev
              if(tdmu(l) .le. emax2) dir(l,m) = solflx*u0(m)*
     -                                          dexp(-tdmu(l))
1561      continue
c
c****       find the upward flux at the top and the downward
c           flux at the bottom of each homogeneous model layer.
c
          ufl(nlev) = alb*dir(nlev,m)
          dfl(nlev) = 0.
          do 1621 l=1,nlay
              ufl(l) = omp(l)*dir(l,m)/da(l)*
     *                 (dc(l)*(a2(l) + sk(l)*g3(l))*ekt(l) -
     *                 db(l)*(a2(l) - sk(l)*g3(l))/ekt(l) -
     *                 2.*sk(l)*(g3(l) - a2(l)*u0(m))*etu(l))
              dfl(l) = -dir(l+1,m)*omp(l)/da(l)*
     *                 (db(l)*(a1(l) + sk(l)*g4(l))*ekt(l) -
     *                 dc(l)*(a1(l) - sk(l)*g4(l))/ekt(l) -
     *                 2.*sk(l)*(g4(l) + a1(l)*u0(m))/etu(l))
1621      continue
c
c****       use adding method to find upward and downward fluxes
c           for combined layers.  start at top
c
          dflux(1) = dfl(1) + tl(1)*flxd(1,m)
          do 1701 l=1,nlm
c              dflux(l+1) = tl(l+1)*(rs(l)*(rl(l+1)*dflux(l) +
c     *                     ufl(l+1))*dd(l) + dflux(l)) + dfl(l+1)
              dflux(l+1) = dfl(l+1) + tl(l+1)*
     *                     (dflux(l + rs(l)*ufl(l+1))*dd(l))
1701      continue
          flxd(nlev,m) = (dflux(nlay) + rs(nlay)*dir(nlev,m)*alb)/
     *                   (1. - rs(nlay)*alb)
c
c****       use adding method to find upward and downward fluxes
c           for combined layers.  start at bottom.
c
          uflux(nlev) = ufl(nlev)
          do 1721 l=1,nlay
              uflux(nlev-l) = tl(nlev-l)*(uflux(np2-l) +
     *                        ru(np2-l)*dfl(nlev-l))*du(nlev-l) +
     *                        ufl(nlev-l)
1721  continue
c
c****       find the total upward and downward fluxes at interfaces
c           between inhomogeneous layers.
c
          flxu(1,m) = uflux(1) + ru(1)*flxd(1,m)
          flxu(nlev,m) = alb*(flxd(nlev,m) + dir(nlev,m))
          do 1741 l=1,nlm
              flxu(nlev-l,m) = (uflux(nlev-l) + ru(nlev-l)*
     *                         dflux(nlay-l))/(1. - ru(nlev-l)*
     *                         rs(nlay-l))
1741      continue
          do 1782 l=1,nlm
              flxd(l+1,m) = dflux(l) + rs(l)*flxu(l+1,m)
1782  continue
c
      return
      end
