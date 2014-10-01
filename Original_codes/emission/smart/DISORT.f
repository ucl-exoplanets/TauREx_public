c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c rcs version control information:
c  disort.f,v 2.1 2000/04/04 18:21:55 laszlo exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine disort( nlyr, dtauc, ssalb, nmom, pmom, temper, 
     $                   wvnm, usrtau, ntau, utau, nstr, iunits,
     &                   usrang, numu, umu, nphi, phi, ibcnd, fbeam, 
     &                   umu0, phi0, fisot, lamber, iref, 
     &                   surf_pr, albedo, btemp, ttemp,temis, 
     &                   plank, onlyfl, accur, prnt, 
     &                   maxcly,maxulv, maxumu, maxphi, maxmom, 
     &                   rfldir, rfldn,flup, dfdt, uavg, uu, 
     &                   albmed, trnmed )

c *******************************************************************
c       plane-parallel discrete ordinates radiative transfer program
c             ( see disort.doc for complete documentation )
c *******************************************************************
c
c***** this version passes the surface reflection properties 
c      through to bdref
c
c +------------------------------------------------------------------+
c  calling tree (omitting calls to errmsg):
c  (routines in parentheses are not in this file)
c
c  disort-+-(r1mach)
c         +-slftst-+-(tstbad)
c         +-zeroit
c         +-chekin-+-(wrtbad)
c         |        +-(wrtdim)
c         |        +-dref
c         +-zeroal
c         +-setdis-+-qgausn-+-(d1mach)
c         +-prtinp
c         +-albtrn-+-lepoly
c         |        +-zeroit
c         |        +-soleig-+-asymtx-+-(d1mach)
c         |        +-terpev
c         |        +-setmtx-+-zeroit
c         |        +-(sgbco)
c         |        +-solve1-+-zeroit
c         |        |        +-(sgbsl)
c         |        +-altrin
c         |        +-spaltr
c         |        +-praltr
c         +-planck-+-(r1mach)
c         +-lepoly
c         +-surfac-+-qgausn-+-(d1mach)
c         |        +-bdref
c         |        +-zeroit
c         +-soleig-+-asymtx-+-(d1mach)
c         +-upbeam-+-(sgeco)
c         |        +-(sgesl)
c         +-upisot-+-(sgeco)
c         |        +-(sgesl)
c         +-terpev
c         +-terpso
c         +-setmtx-+-zeroit
c         +-solve0-+-zeroit
c         |        +-(sgbco)
c         |        +-(sgbsl)
c         +-fluxes--zeroit
c         +-zeroit
c         +-usrint
c         +-cmpint
c         +-pravin
c         +-zeroit
c         +-ratio--(r1mach)
c         +-intcor-+-sinsca
c         |        +-secsca-+-xifunc
c         +-prtint
c
c *** intrinsic functions used in disort package which take
c     non-negligible amount of time:
c
c    exp :  called by- albtrn, altrin, cmpint, fluxes, setdis,
c                      setmtx, spaltr, usrint, planck
c
c    sqrt : called by- asymtx, soleig
c
c +-------------------------------------------------------------------+
c
c  index conventions (for all do-loops and all variable descriptions):
c
c     iu     :  for user polar angles
c
c  iq,jq,kq  :  for computational polar angles ('quadrature angles')
c
c   iq/2     :  for half the computational polar angles (just the ones
c               in either 0-90 degrees, or 90-180 degrees)
c
c     j      :  for user azimuthal angles
c
c     k,l    :  for legendre expansion coefficients or, alternatively,
c               subscripts of associated legendre polynomials
c
c     lu     :  for user levels
c
c     lc     :  for computational layers (each having a different
c               single-scatter albedo and/or phase function)
c
c    lev     :  for computational levels
c
c    mazim   :  for azimuthal components in fourier cosine expansion
c               of intensity and phase function
c +-------------------------------------------------------------------+
c
c   variables added to the subroutine call list
c
c     iref   :  surface reflectance mode
c             0 - Lambert
c             1 - Hapke's BDR model
c             2 - Breon's BDR model; combination of Li + Roujean
c             3 - Roujean's BDR model
c             4 - Cox and Munk glint model
c   SURF_PR : Wavelength dependent surface properties array 
c             IREF= 0 - Lambert albedo
c             IREF= 1 - Hapke : HH, W
c             IREF= 2 - Breon's BDR model: k0, k1, k2
c             IREF= 3 - Roujean's BDR model: k0, k1, k2
c             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw
c
c +------------------------------------------------------------------+
c
c               i n t e r n a l    v a r i a b l e s
c
c   amb(iq/2,iq/2)    first matrix factor in reduced eigenvalue problem
c                     of eqs. ss(12), stwj(8e), stwl(23f)
c                     (used only in soleig)
c
c   apb(iq/2,iq/2)    second matrix factor in reduced eigenvalue problem
c                     of eqs. ss(12), stwj(8e), stwl(23f)
c                     (used only in soleig)
c
c   array(iq,iq)      scratch matrix for soleig, upbeam and upisot
c                     (see each subroutine for definition)
c
c   b()               right-hand side vector of eq. sc(5) going into
c                     solve0,1;  returns as solution vector
c                     vector  l, the constants of integration
c
c   bdr(iq/2,0:iq/2)  bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  first index always
c                     refers to a computational angle.  second index:
c                     if zero, refers to incident beam angle umu0;
c                     if non-zero, refers to a computational angle.
c
c   bem(iq/2)         bottom-boundary directional emissivity at compu-
c                     tational angles.
c
c   bplank            intensity emitted from bottom boundary
c
c   cband()           matrix of left-hand side of the linear system
c                     eq. sc(5), scaled by eq. sc(12);  in banded
c                     form required by linpack solution routines
c
c   cc(iq,iq)         c-sub-ij in eq. ss(5)
c
c   cmu(iq)           computational polar angles (gaussian)
c
c   cwt(iq)           quadrature weights corresponding to cmu
c
c   corint            when set true, correct intensities for
c                     delta-scaling effects (see nakajima and tanaka,
c                     1988). when false, intensities are not corrected.
c                     in general, corint should be set true when beam
c                     source is present (fbeam is not zero) asurf_prnd deltam
c                     is true in a problem including scattering.
c                     however, execution is faster when corint is false,
c                     and intensities outside the aureole may still be
c                     accurate enough.  when corint is true, it is
c                     important to have a sufficiently high order of
c                     legendre approximation of the phase function. this
c                     is because the intensities are corrected by
c                     calculating the single-scattered radiation, for
c                     which an adequate representation of the phase
c                     function is crucial.  in case of a low order
c                     legendre approximation of an otherwise highly
c                     anisotropic phase function, the intensities might
c                     actually be more accurate when corint is false.
c                     when only fluxes are calculated (onlyfl is true),
c                     or there is no beam source (fbeam=0.0),surf_pr or there
c                     is no scattering (ssalb=0.0 for all layers) corint
c                     is set false by the code.
c
c   delm0             kronecker delta, delta-sub-m0, where m = mazim
c                     is the number of the fourier component in the
c                     azimuth cosine expansion
c
c   deltam            true,  use delta-m method ( see wiscombe, 1977 );
c                     false, do not use delta-m method. in general, for
c                     a given number of streams, intensities and
c                     fluxes will be more accurate for phase functions
c                     with a large forward peak if deltam is set true.
c                     intensities close to the forward scattering
c                     direction are often less accurate, however, when
c                     the delta-m method is applied. the intensity
c                     correction of nakajima and tanaka is used to
c                     improve the accuracy of the intensities.
c
c   dither            small quantity subtracted from single-scattering
c                     albedos of unity, in order to avoid usurf_prsing special
c                     case formulas;  prevents an eigenvalue of exactly
c                     zero from occurring, which would cause an
c                     immediate overflow
c
c   dtaucp(lc)        computational-layer optical depths (delta-m-scaled
c                     if deltam = true, otherwise equal to dtauc)
c
c   emu(iu)           bottom-boundary directional emissivity at user
c                     angles.
c
c   eval(iq)          temporary storage for eigenvalues of eq. ss(12)
c
c   evecc(iq,iq)      complete eigenvectors of ss(7) on return from
c                     soleig; stored permanently in  gc
c
c   expbea(lc)        transmission of direct beam in delta-m optical
c                     depth coordinates
c
c   flyr(lc)          separated fraction in delta-m method
csurf_pr
c   gl(k,lc)          phase function legendre polynomial expansion
c                     coefficients, calculated from pmom by
c                     including single-scattering albedo, factor
c                     2k+1, and (if deltam=true) the delta-m
c                     scaling
c
c   gc(iq,iq,lc)      eigenvectors at polar quadrature angles,
c                     g  in eq. sc(1)
c
c   gu(iu,iq,lc)      eigenvectors interpolated to user polar angles
c                     ( g  in eqs. sc(3) and s1(8-9), i.e.
c                       g without the l factor )
c
c   ipvt(lc*iq)       integer vector of pivot indices for linpack
c                     routines
csurf_pr
c   kk(iq,lc)         eigenvalues of coeff. matrix in eq. ss(7)
c
c   kconv             counter in azimuth convergence test
c
c   layru(lu)         computational layer in which user output level
c                     utau(lu) is located
c
c   ll(iq,lc)         constants of integration l in eq. sc(1),
c                     obtained by solving scaled version of eq. sc(5)
c
c   lyrcut            true, radiation is assumed zero below layer
c                     ncut because of almost complete absorption
c
c   naz               number of azimuthal components considered
c
c   ncut              computational layer number in which absorption
c                     optical depth first exceeds abscut
c
c   oprim(lc)         single scattering albedo after delta-m scaling
c
c   pass1             true on first entry, false thereafter
c
c   pkag(0:lc)        integrated planck function for internal emission
c
c   prntu0(l)         logical flag to trigger printing of azimuthally-
c                     averaged intensities:
c                       l    quantities printed
c                      --    ------------------
c                       1    azimuthally-averaged intensities at user
c                               levels and computational polar angles
c                       2    azimuthally-averaged intensities at user
c                               levels and user polar angles
c
c   psi0(iq)          sum just after square bracket in  eq. sd(9)
c
c   psi1(iq)          sum in  eq. stwl(31d)
c
c   rmu(iu,0:iq)      bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  first index always
c                     refers to a user angle.  second index:
c                     if zero, refers to incident beam angle umu0;
c                     if non-zero, refers to a computational angle.
c
c   sqt(k)            square root of k (used only in lepoly for
c                     computing associated legendre polynomials)
c
c   tauc(0:lc)        cumulative optical depth (un-delta-m-scaled)
c
c   taucpr(0:lc)      cumulative optical depth (delta-m-scaled if
c                     deltam = true, otherwise equal to tauc)
c
c   tplank            intensity emitted from top boundary
c
c   uum(iu,lu)        expansion coefficients when the intensity
c                     (u-super-m) is expanded in fourier cosine series
c                     in azimuth angle
c
c   u0c(iq,lu)        azimuthally-averaged intensity at quadrature
c                     angle
c
c   u0u(iu,lu)        if onlyfl = false, azimuthally-averaged intensity
c                     at user angles and user levels
c
c                     if onlyfl = true and maxumu.ge.nstr,
c                     azimuthally-averaged intensity at computational
c                     (gaussian quadrature) angles and user levels;
c                     the corresponding quadrature angle cosines are
c                     returned in umu.  if maxumu.lt.nstr, u0u will be
c                     zeroed, and umu, numu will not be set.
c
c   utaupr(lu)        optical depths of user output levels in delta-m
c                     coordinates;  equal to  utau(lu) if no delta-m
c
c   wk()              scratch array
c
c   xr0(lc)           x-sub-zero in expansion of thermal source func-
c                     tion preceding eq. ss(14)(has no mu-dependence);
c                     b-sub-zero in eq. stwl(24d)
c
c   xr1(lc)           x-sub-one in expansion of thermal source func-
c                     tion; see  eqs. ss(14-16); b-sub-one in stwl(24d)
c
c   ylm0(l)           normalized associated legendre polynomial
c                     of subscript l at the beam angle (not saved
c                     as function of superscipt m)
c
c   ylmc(l,iq)        normalized associated legendre polynomial
c                     of subscript l at the computational angles
c                     (not saved as function of superscipt m)
c
c   ylmu(l,iu)        normalized associated legendre polynomial
c                     of subscript l at the user angles
c                     (not saved as function of superscipt m)
c
c   z()               scratch array used in solve0, albtrn to solve
c                     a linear system for the constants of integration
c
c   z0(iq)            solution vectors z-sub-zero of eq. ss(16)
c
c   z0u(iu,lc)        z-sub-zero in eq. ss(16) interpolated to user
c                     angles from an equation derived from ss(16)
c
c   z1(iq)            solution vectors z-sub-one  of eq. ss(16)
c
c   z1u(iu,lc)        z-sub-one in eq. ss(16) interpolated to user
c                     angles from an equation derived from ss(16)
c
c   zbeam(iu,lc)      particular solution for beam source
c
c   zj(iq)            right-hand side vector  x-sub-zero in
c                     eq. ss(19), also the solution vector
c                     z-sub-zero after solving that system
c
c   zz(iq,lc)         permanent storage for the beam source vectors zj
c
c   zplk0(iq,lc)      permanent storage for the thermal source
c                     vectors  z0  obtained by solving  eq. ss(16)
c
c   zplk1(iq,lc)      permanent storage for the thermal source
c                     vectors  z1  obtained by solving  eq. ss(16)
c
c +-------------------------------------------------------------------+
c
c  local symbolic dimensions (have big effect on storage requirements):
c
c       mxcly  = max no. of computational layers
c       mxulv  = max no. of output levels
c       mxcmu  = max no. of computation polar angles
c       mxumu  = max no. of output polar angles
c       mxphi  = max no. of output azimuthal angles
c       mxsqt  = max no. of square roots of integers (for lepoly)
c +-------------------------------------------------------------------+

c     .. parameters ..

      integer   mxcly, mxulv, mxcmu, mxumu, mxphi, mi, mi9m2, nnlyri,
     &          mxsqt
      parameter ( mxcly = 70, mxulv = 70, mxcmu = 16, mxumu = 16,
     &          mxphi = 16, mi = mxcmu / 2, mi9m2 = 9*mi - 2,
     &          nnlyri = mxcmu*mxcly, mxsqt = 1000 )
c     ..
c     .. scalar arguments ..

      logical   lamber, onlyfl, plank, usrang, usrtau
      integer   ibcnd, maxcly, maxmom, maxphi, maxulv, maxumu, nlyr,
     &          nmom, nphi, nstr, ntau, numu, iunits
      real      accur, albedo, btemp, fbeam, fisot, phi0, temis, ttemp,
     &          umu0, wvnm
c     ..
c     .. array arguments ..

      logical   prnt( 5 )
      real      albmed( maxumu ), dfdt( maxulv ), dtauc( maxcly ),
     &          flup( maxulv ), phi( maxphi ), pmom( 0:maxmom, maxcly ),
     &          rfldir( maxulv ), rfldn( maxulv ), ssalb( maxcly ),
     &          temper( 0:maxcly ), trnmed( maxumu ), uavg( maxulv ),
     &          umu( maxumu ), utau( maxulv ),
     &          uu( maxumu, maxulv, maxphi )
c     ..
c     .. local scalars ..

      logical   compar, corint, deltam, lyrcut, pass1
      integer   iq, iu, j, kconv, l, lc, lev, lu, mazim, naz, ncol,
     &          ncos, ncut, nn, ns
      real      azerr, azterm, bplank, cosphi, delm0, dither,
     &          dum, pi, rpd, sgn, tplank
      real      angcos(1)
c     ..
c     .. local arrays ..

      logical   prntu0( 2 )
      integer   ipvt( nnlyri ), layru( mxulv )
      integer   iref

      real      amb( mi, mi ), apb( mi, mi ), array( mxcmu, mxcmu ),
     &          b( nnlyri ), bdr( mi, 0:mi ), bem( mi ),
     &          cband( mi9m2, nnlyri ), cc( mxcmu, mxcmu ),
     &          cmu( mxcmu ), cwt( mxcmu ), dtaucp( mxcly ),
     &          emu( mxumu ), eval( mi ), evecc( mxcmu, mxcmu ),
     &          expbea( 0:mxcly ), fldir( mxulv ), fldn( mxulv ),
     &          flyr( mxcly ), gc( mxcmu, mxcmu, mxcly ),
     &          gl( 0:mxcmu, mxcly ), gu( mxumu, mxcmu, mxcly ),
     &          kk( mxcmu, mxcly ), ll( mxcmu, mxcly ), oprim( mxcly ),
     &          phasa( mxcly ), phast( mxcly ), phasm( mxcly ),
     &          phirad( mxphi ), pkag( 0:mxcly ), psi0( mxcmu ),
     &          psi1( mxcmu ), rmu( mxumu, 0:mi ), sqt( mxsqt ),
     &          tauc( 0:mxcly ), taucpr( 0:mxcly ), u0c( mxcmu, mxulv ),
     &          u0u( mxumu, mxulv ), utaupr( mxulv ),
     &          uum( mxumu, mxulv ), wk( mxcmu ), xr0( mxcly ),
     &          xr1( mxcly ), ylm0( 0:mxcmu ), ylmc( 0:mxcmu, mxcmu ),
     &          ylmu( 0:mxcmu, mxumu ), z( nnlyri ), z0( mxcmu ),
     &          z0u( mxumu, mxcly ), z1( mxcmu ), z1u( mxumu, mxcly ),
     &          zbeam( mxumu, mxcly )
      real      zj( mxcmu ), zplk0( mxcmu, mxcly ),
     &          zplk1( mxcmu, mxcly ), zz( mxcmu, mxcly )
      real      surf_pr( * )
c
      real tempbb( mxcly )
      double precision bb( mxcly )

      double precision aad( mi, mi ), evald( mi ), eveccd( mi, mi ),
     &                 wkd( mxcmu )
c     ..
c     .. external functions ..

      real      r1mach, ratio
      external  r1mach, ratio
c     ..
c     .. external subroutines ..

      external  albtrn, chekin, cmpint, fluxes, intcor, lepoly, pravin,
     &          prtinp, prtint, setdis, setmtx, slftst, soleig, solve0,
     &          surfac, terpev, terpso, upbeam, upisot, usrint, zeroal,
     &          zeroit
c     ..
c     .. intrinsic functions ..

      intrinsic abs, asin, cos, float, max, sqrt
c     ..
      save      dither, pass1, pi, rpd, sqt
      data      pass1 / .true. /, prntu0 / 2*.false. /

      deltam = .true.
      corint = .true.

      if( pass1 ) then

         pi     = 2.*asin( 1.0 )
         dither = 10.*r1mach( 4 )

c                            ** must dither more on high (14-digit)
c                            ** precision machine

         if( dither.lt.1.e-10 ) dither = 10.*dither

         rpd  = pi / 180.0

         do 10 ns = 1, mxsqt
            sqt( ns ) = sqrt( float( ns ) )
   10    continue
c                            ** set input values for self-test.
c                            ** be sure slftst sets all print flags off.
         compar = .false.

         call slftst( corint, accur, albedo, btemp, deltam, dtauc( 1 ),
     &                fbeam, fisot, ibcnd, lamber, nlyr, plank, nphi,
     &                numu, nstr, ntau, onlyfl, phi( 1 ), phi0, nmom,
     &                pmom( 0,1 ), prnt, prntu0, ssalb( 1 ), temis,
     &                temper( 0 ), ttemp, umu( 1 ), usrang, usrtau,
     &                utau( 1 ), umu0, wvnm, compar, dum,
     &                dum, dum, dum )

      end if


   20 continue
c
c                                  ** calculate cumulative optical depth
c                                  ** and dither single-scatter albedo
c                                  ** to improve numerical behavior of
c                                  ** eigenvalue/vector computation
      call zeroit( tauc, mxcly + 1 )

      do 30 lc = 1, nlyr

         if( ssalb( lc ).eq.1.0 ) ssalb( lc ) = 1.0 - dither
         tauc( lc ) = tauc( lc - 1 ) + dtauc( lc )

   30 continue
c                                ** check input dimensions and variables

      call chekin( nlyr, dtauc, ssalb, nmom, pmom, temper, wvnm,
     &             usrtau, ntau, iref, utau, nstr, usrang,
     &             numu, umu, nphi, phi, ibcnd, fbeam, umu0,
     &             phi0, fisot, lamber, albedo, btemp, ttemp,temis, 
     &             surf_pr, plank, onlyfl, deltam, corint, accur,
     &             tauc, maxcly, maxulv, maxumu, maxphi, maxmom,
     &             mxcly, mxulv, mxumu, mxcmu, mxphi, mxsqt )

c                                 ** zero internal and output arrays

      call  zeroal( mxcly, expbea(1), flyr, oprim, phasa, phast, phasm,
     &                     taucpr(1), xr0, xr1,
     &              mxcmu, cmu, cwt, psi0, psi1, wk, z0, z1, zj,
     &              mxcmu+1, ylm0,
     &              mxcmu**2, array, cc, evecc,
     &              (mxcmu+1)*mxcly, gl,
     &              (mxcmu+1)*mxcmu, ylmc,
     &              (mxcmu+1)*mxumu, ylmu,
     &              mxcmu*mxcly, kk, ll, zz, zplk0, zplk1,
     &              mxcmu**2*mxcly, gc,
     &              mxulv, layru, utaupr,
     &              mxumu*mxcmu*mxcly, gu,
     &              mxumu*mxcly, z0u, z1u, zbeam,
     &              mi, eval,
     &              mi**2, amb, apb,
     &              nnlyri, ipvt, z,
     &              maxulv, rfldir, rfldn, flup, uavg, dfdt,
     &              maxumu, albmed, trnmed,
     &              mxumu*mxulv, u0u,
     &              maxumu*maxulv*maxphi, uu )

c                                 ** perform various setup operations

      call setdis( cmu, cwt, deltam, dtauc, dtaucp, expbea, fbeam, flyr,
     &             gl, ibcnd, layru, lyrcut, maxmom, maxumu, mxcmu,
     &             ncut, nlyr, ntau, nn, nstr, plank, numu, onlyfl,
     &             corint, oprim, pmom, ssalb, tauc, taucpr, utau,
     &             utaupr, umu, umu0, usrtau, usrang )

c                                 ** print input information
      if( prnt( 1 ) )
     &    call prtinp( nlyr, dtauc, dtaucp, ssalb, nmom, pmom, temper,
     &                 wvnm, ntau, utau, nstr, numu, umu,
     &                 nphi, phi, ibcnd, fbeam, umu0, phi0, fisot,
     &                 lamber, albedo, btemp, ttemp, temis, deltam,
     &                 plank, onlyfl, corint, accur, flyr, lyrcut,
     &                 oprim, tauc, taucpr, maxmom, prnt( 5 ) )

c                              ** handle special case for getting albedo
c                              ** and transmissivity of medium for many
c                              ** beam angles at once
      if( ibcnd.eq.1 ) then

         call albtrn( albedo, amb, apb, array, b, bdr, cband, cc, cmu,
     &                cwt, dtaucp, eval, evecc, gl, gc, gu, ipvt, kk,
     &                ll, nlyr, nn, nstr, numu, prnt, taucpr, umu, u0u,
     &                wk, ylmc, ylmu, z, aad, evald, eveccd, wkd, mi,
     &                mi9m2, maxumu, mxcmu, mxumu, nnlyri, sqt, albmed,
     &                trnmed )
         return

      end if
c                                   ** calculate planck functions
      if( .not.plank ) then

         bplank = 0.0
         tplank = 0.0
         call zeroit( pkag,  mxcly + 1 )

      else

         tempbb(1) = ttemp
         call planck(1,1,1,iunits,wvnm,tempbb,bb)
         tplank = real(temis*bb(1))
         tempbb(1) = btemp
         call planck(1,1,1,1,wvnm,tempbb,bb)
         bplank = real(bb(1))

         do 39 lev = 0, nlyr
             tempbb(lev+1) = temper(lev)
39       continue
c
         call planck(1,nlyr+1,1,iunits,wvnm,tempbb,bb)
c
         do 40 lev = 0, nlyr
            pkag( lev ) = real(bb(lev+1))
   40    continue

      end if

c
c ========  begin loop to sum azimuthal components of intensity  =======
c           (eq stwj 5, stwl 6)

      kconv  = 0
      naz    = nstr - 1
c                                    ** azimuth-independent case

      if( fbeam.eq.0.0 .or. abs(1.-umu0).lt.1.e-5 .or. onlyfl .or.
     &   ( numu.eq.1 .and. abs(1.-umu(1)).lt.1.e-5 ) .or.
     &   ( numu.eq.1 .and. abs(1.+umu(1)).lt.1.e-5 ) .or.
     &   ( numu.eq.2 .and. abs(1.+umu(1)).lt.1.e-5 .and.
     &     abs(1.-umu(2)).lt.1.e-5 ) )
     &   naz = 0


      do 180 mazim = 0, naz

         if( mazim.eq.0 ) delm0  = 1.0
         if( mazim.gt.0 ) delm0  = 0.0

c                             ** get normalized associated legendre
c                             ** polynomials for
c                             ** (a) incident beam angle cosine
c                             ** (b) computational and user polar angle
c                             **     cosines
         if( fbeam.gt.0.0 ) then

            ncos   = 1
            angcos(1) = -umu0

            call lepoly( ncos, mazim, mxcmu, nstr-1, angcos, sqt, ylm0 )

         end if


         if( .not.onlyfl .and. usrang )
     &       call lepoly( numu, mazim, mxcmu, nstr-1, umu, sqt, ylmu )

         call lepoly( nn, mazim, mxcmu, nstr-1, cmu, sqt, ylmc )

c                       ** get normalized associated legendre polys.
c                       ** with negative arguments from those with
c                       ** positive arguments; dave/armstrong eq. (15),
c                       ** stwl(59)
         sgn  = -1.0

         do 70 l = mazim, nstr - 1

            sgn  = -sgn

            do 60 iq = nn + 1, nstr
               ylmc( l, iq ) = sgn*ylmc( l, iq - nn )
   60       continue

   70    continue
c                                 ** specify users bottom reflectivity
c                                 ** and emissivity properties
         if( .not.lyrcut )
     &       call surfac( albedo, delm0, cmu, fbeam, lamber, iref, mi, 
     &                    mazim, mxumu, nn, numu, onlyfl, pi, umu, umu0,
     &                    usrang, surf_pr, bdr, emu,bem, rmu )


c ===================  begin loop on computational layers  =============

         do 80 lc = 1, ncut

c                      ** solve eigenfunction problem in eq. stwj(8b),
c                      ** stwl(23f); return eigenvalues and eigenvectors

            call soleig( amb, apb, array, cmu, cwt, gl( 0,lc ), mi,
     &                   mazim, mxcmu, nn, nstr, ylmc, cc, evecc, eval,
     &                   kk( 1,lc ), gc( 1,1,lc ), aad, eveccd, evald,
     &                   wkd )

c                                  ** calculate particular solutions of
c                                  ** eq. ss(18), stwl(24a) for incident
c                                  ** beam source
            if( fbeam.gt.0.0 )
     &          call upbeam( array, cc, cmu, delm0, fbeam, gl( 0,lc ),
     &                       ipvt, mazim, mxcmu, nn, nstr, pi, umu0, wk,
     &                       ylm0, ylmc, zj, zz( 1,lc ) )

c                              ** calculate particular solutions of eq.
c                              ** ss(15), stwl(25) for thermal emission
c                              ** source
c
            if( plank .and. mazim.eq.0 ) then

               xr1( lc ) = 0.0

               if( dtaucp( lc ).gt.0.0 ) xr1( lc ) =
     &             ( pkag( lc ) - pkag( lc-1 ) ) / dtaucp( lc )

               xr0( lc ) = pkag( lc-1 ) - xr1( lc )*taucpr( lc-1 )

               call upisot( array, cc, cmu, ipvt, mxcmu, nn, nstr,
     &                      oprim( lc ), wk, xr0( lc ), xr1( lc ),
     &                      z0, z1, zplk0( 1,lc ), zplk1( 1,lc ) )
            end if


            if( .not.onlyfl .and. usrang ) then

c                                            ** interpolate eigenvectors
c                                            ** to user angles

               call terpev( cwt, evecc, gl( 0,lc ), gu( 1,1,lc ), mazim,
     &                      mxcmu, mxumu, nn, nstr, numu, wk, ylmc,
     &                      ylmu )
c                                            ** interpolate source terms
c                                            ** to user angles

               call terpso( cwt, delm0, fbeam, gl( 0,lc ), mazim, mxcmu,
     &                      plank, numu, nstr, oprim( lc ), pi, ylm0,
     &                      ylmc, ylmu, psi0, psi1, xr0( lc ),
     &                      xr1( lc ), z0, z1, zj, zbeam( 1,lc ),
     &                      z0u( 1,lc ), z1u( 1,lc ) )

            end if

   80    continue

c ===================  end loop on computational layers  ===============


c                      ** set coefficient matrix of equations combining
c                      ** boundary and layer interface conditions

         call setmtx( bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk,
     &                lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut,
     &                nnlyri, nn, nstr, taucpr, wk )

c                      ** solve for constants of integration in homo-
c                      ** geneous solution (general boundary conditions)

         call solve0( b, bdr, bem, bplank, cband, cmu, cwt, expbea,
     &                fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim, mi,
     &                mi9m2, mxcmu, ncol, ncut, nn, nstr, nnlyri, pi,
     &                tplank, taucpr, umu0, z, zz, zplk0, zplk1 )

c                                  ** compute upward and downward fluxes

         if( mazim.eq.0 )
     &       call fluxes( cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut,
     &                    maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau,
     &                    pi, prnt, prntu0( 1 ), ssalb, taucpr, umu0,
     &                    utau, utaupr, xr0, xr1, zz, zplk0, zplk1,
     &                    dfdt, flup, fldn, fldir, rfldir, rfldn, uavg,
     &                    u0c )

         if( onlyfl ) then

            if( maxumu.ge.nstr ) then
c                                     ** save azimuthal-avg intensities
c                                     ** at quadrature angles
               do 100 lu = 1, ntau

                  do 90 iq = 1, nstr
                     u0u( iq, lu ) = u0c( iq, lu )
   90             continue

  100          continue

            end if

            go to  190

         end if


         call zeroit( uum, mxumu*mxulv )

         if( usrang ) then
c                                     ** compute azimuthal intensity
c                                     ** components at user angles

            call usrint( bplank, cmu, cwt, delm0, dtaucp, emu, expbea,
     &                   fbeam, fisot, gc, gu, kk, lamber, layru, ll,
     &                   lyrcut, mazim, mxcmu, mxulv, mxumu, ncut, nlyr,
     &                   nn, nstr, plank, numu, ntau, pi, rmu, taucpr,
     &                   tplank, umu, umu0, utaupr, wk, zbeam, z0u, z1u,
     &                   zz, zplk0, zplk1, uum )

         else
c                                     ** compute azimuthal intensity
c                                     ** components at quadrature angles

            call cmpint( fbeam, gc, kk, layru, ll, lyrcut, mazim, mxcmu,
     &                   mxulv, mxumu, ncut, nn, nstr, plank, ntau,
     &                   taucpr, umu0, utaupr, zz, zplk0, zplk1, uum )
         end if


         if( mazim.eq.0 ) then
c                               ** save azimuthally averaged intensities

            do 130 lu = 1, ntau

               do 120 iu = 1, numu
                  u0u( iu, lu ) = uum( iu, lu )

                  do 110 j = 1, nphi
                     uu( iu, lu, j ) = uum( iu, lu )
  110             continue

  120          continue

  130       continue
c                              ** print azimuthally averaged intensities
c                              ** at user angles

            if( prntu0( 2 ) )
     &          call pravin( umu, numu, mxumu, utau, ntau, u0u )

            if( naz.gt.0 ) then

               call zeroit( phirad, mxphi )
               do 140 j = 1, nphi
                  phirad( j ) = rpd*( phi( j ) - phi0 )
  140          continue

            end if


         else
c                                ** increment intensity by current
c                                ** azimuthal component (fourier
c                                ** cosine series);  eq sd(2), stwl(6)
            azerr  = 0.0

            do 170 j = 1, nphi

               cosphi = cos( mazim*phirad( j ) )

               do 160 lu = 1, ntau

                  do 150 iu = 1, numu
                     azterm = uum( iu, lu )*cosphi
                     uu( iu, lu, j ) = uu( iu, lu, j ) + azterm
                     azerr  = max( azerr,
     &                        ratio( abs(azterm), abs(uu(iu,lu,j)) ) )
  150             continue

  160          continue

  170       continue

            if( azerr.le.accur ) kconv  = kconv + 1

            if( kconv.ge.2 ) go to  190

         end if

  180 continue

c ===================  end loop on azimuthal components  ===============


  190 continue

c                                    ** apply nakajima/tanaka intensity
c                                    ** corrections

      if( corint )
     &    call intcor( dither, fbeam, flyr, layru, lyrcut, maxmom,
     &                 maxulv, maxumu, nmom, ncut, nphi, nstr, ntau,
     &                 numu, oprim, phasa, phast, phasm, phirad, pi,
     &                 rpd, pmom, ssalb, dtauc, tauc, taucpr, umu,
     &                 umu0, utau, utaupr, uu )

c                                          ** print intensities

      if( prnt( 3 ) .and. .not.onlyfl )
     &    call prtint( uu, utau, ntau, umu, numu, phi, nphi, maxulv,
     &                 maxumu )


      if( pass1 ) then
c                                    ** compare test case results with
c                                    ** correct answers and abort if bad
         compar = .true.

         call slftst( corint, accur, albedo, btemp, deltam, dtauc( 1 ),
     &                fbeam, fisot, ibcnd, lamber, nlyr, plank, nphi,
     &                numu, nstr, ntau, onlyfl, phi( 1 ), phi0, nmom,
     &                pmom( 0,1 ), prnt, prntu0, ssalb( 1 ), temis,
     &                temper( 0 ), ttemp, umu( 1 ), usrang, usrtau,
     &                utau( 1 ), umu0, wvnm, compar,
     &                flup( 1 ), rfldir( 1 ), rfldn( 1 ), uu( 1,1,1 ) )

         pass1  = .false.
         go to  20

      end if


      return
      end

      subroutine asymtx( aa, evec, eval, m, ia, ievec, ier, wkd, aad,
     &                   evecd, evald )

c    =======  d o u b l e    p r e c i s i o n    v e r s i o n  ======

c       solves eigenfunction problem for real asymmetric matrix
c       for which it is known a priori that the eigenvalues are real.

c       this is an adaptation of a subroutine eigrf in the imsl
c       library to use real instead of complex arithmetic, accounting
c       for the known fact that the eigenvalues and eigenvectors in
c       the discrete ordinate solution are real.  other changes include
c       putting all the called subroutines in-line, deleting the
c       performance index calculation, updating many do-loops
c       to fortran77, and in calculating the machine precision
c       tol instead of specifying it in a data statement.

c       eigrf is based primarily on eispack routines.  the matrix is
c       first balanced using the parlett-reinsch algorithm.  then
c       the martin-wilkinson algorithm is applied.

c       there is a statement 'j  = wkd( i )' that converts a double
c       precision variable to an integer variable, that seems dangerous
c       to us in principle, but seems to work fine in practice.

c       references:
c          dongarra, j. and c. moler, eispack -- a package for solving
c             matrix eigenvalue problems, in cowell, ed., 1984:
c             sources and development of mathematical software,
c             prentice-hall, englewood cliffs, nj
c         parlett and reinsch, 1969: balancing a matrix for calculation
c             of eigenvalues and eigenvectors, num. math. 13, 293-304
c         wilkinson, j., 1965: the algebraic eigenvalue problem,
c             clarendon press, oxford


c   i n p u t    v a r i a b l e s:
c
c       aa    :  input asymmetric matrix, destroyed after solved
c
c        m    :  order of  aa
c
c       ia    :  first dimension of  aa
c
c    ievec    :  first dimension of  evec
c
c
c   o u t p u t    v a r i a b l e s:
c
c       evec  :  (unnormalized) eigenvectors of  aa
c                   ( column j corresponds to eval(j) )
c
c       eval  :  (unordered) eigenvalues of aa ( dimension at least m )
c
c       ier   :  if .ne. 0, signals that eval(ier) failed to converge;
c                   in that case eigenvalues ier+1,ier+2,...,m  are
c                   correct but eigenvalues 1,...,ier are set to zero.
c
c
c   s c r a t c h   v a r i a b l e s:
c
c       wkd   :  work area ( dimension at least 2*m )
c       aad   :  double precision stand-in for aa
c       evecd :  double precision stand-in for evec
c       evald :  double precision stand-in for eval
c
c   called by- soleig
c   calls- d1mach, errmsg
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical tf
      integer   ia, ier, ievec, m
c     ..
c     .. array arguments ..

      real      aa( ia, m ), eval( m ), evec( ievec, m )
      double precision aad( ia, m ), evald( m ), evecd( ia, m ),
     &                 wkd( * )
c     ..
c     .. local scalars ..

      logical   noconv, notlas
      integer   i, ii, in, j, k, ka, kkk, l, lb, lll, n, n1, n2
      double precision c1, c2, c3, c4, c5, c6, col, discri, f, g, h,
     &                 one, p, q, r, repl, rnorm, row, s, scale, sgn, t,
     &                 tol, uu, vv, w, x, y, z, zero
c     ..
c     .. external functions ..

      double precision d1mach
      external  d1mach
c     ..
c     .. external subroutines ..

      external  errmsg
c     ..
c     .. intrinsic functions ..

      intrinsic abs, min, sign, sqrt
c     ..
      data      c1 / 0.4375d0 / , c2 / 0.5d0 / , c3 / 0.75d0 / ,
     &          c4 / 0.95d0 / , c5 / 16.d0 / , c6 / 256.d0 / ,
     &          zero / 0.d0 / , one / 1.d0 /


      ier  = 0
      tol  = d1mach( 4 )

      if( m.lt.1 .or. ia.lt.m .or. ievec.lt.m ) then
        tf = .true.
        call errmsg( 'asymtx--bad input variable(s)', tf )
      endif

c                           ** handle 1x1 and 2x2 special cases
      if( m.eq.1 ) then

         eval( 1 )   = aa( 1,1 )
         evec( 1,1 ) = 1.0
         return

      else if( m.eq.2 ) then

         discri = ( aa( 1,1 ) - aa( 2,2 ) )**2 + 4.*aa( 1,2 )*aa( 2,1 )

         tf = .true.
         if( discri .lt. 0.0 )
     &      call errmsg( 'asymtx--complex evals in 2x2 case',tf )

         sgn  = one

         if( aa( 1,1 ) .lt. aa( 2,2 ) ) sgn  = - one

         eval( 1 ) = 0.5*( aa( 1,1 ) + aa( 2,2 ) + sgn*sqrt( discri ) )
         eval( 2 ) = 0.5*( aa( 1,1 ) + aa( 2,2 ) - sgn*sqrt( discri ) )
         evec( 1,1 ) = 1.0
         evec( 2,2 ) = 1.0

         if( aa( 1,1 ) .eq. aa( 2,2 ) .and.
     &       ( aa( 2,1 ).eq.0.0 .or. aa( 1,2 ).eq.0.0 ) ) then

            rnorm = abs( aa( 1,1 ) ) + abs( aa( 1,2 ) ) +
     &              abs( aa( 2,1 ) ) + abs( aa( 2,2 ) )
            w     = tol * rnorm
            evec( 2,1 ) =   aa( 2,1 ) / w
            evec( 1,2 ) = - aa( 1,2 ) / w

         else

            evec( 2,1 ) = aa( 2,1 ) / ( eval( 1 ) - aa( 2,2 ) )
            evec( 1,2 ) = aa( 1,2 ) / ( eval( 2 ) - aa( 1,1 ) )

         end if

         return

      end if

c                               ** convert single-prec. matrix to double
      do 20 j = 1, m

         do 10 k = 1, m
            aad( j,k ) = aa( j,k )
   10    continue

   20 continue

c                                ** initialize output variables
      ier  = 0

      do 40 i = 1, m

         evald( i ) = zero

         do 30 j = 1, m
            evecd( i, j ) = zero
   30    continue

         evecd( i, i ) = one

   40 continue

c                  ** balance the input matrix and reduce its norm by
c                  ** diagonal similarity transformation stored in wk;
c                  ** then search for rows isolating an eigenvalue
c                  ** and push them down
      rnorm  = zero
      l  = 1
      k  = m

   50 continue
      kkk  = k

      do 90 j = kkk, 1, -1

         row  = zero

         do 60 i = 1, k
            if( i.ne.j ) row  = row + abs( aad( j,i ) )
   60    continue

         if( row.eq.zero ) then

            wkd( k ) = j

            if( j.ne.k ) then

               do 70 i = 1, k
                  repl        = aad( i, j )
                  aad( i, j ) = aad( i, k )
                  aad( i, k ) = repl
   70          continue

               do 80 i = l, m
                  repl        = aad( j, i )
                  aad( j, i ) = aad( k, i )
                  aad( k, i ) = repl
   80          continue

            end if

            k  = k - 1
            go to  50

         end if

   90 continue
c                                ** search for columns isolating an
c                                ** eigenvalue and push them left
  100 continue
      lll  = l

      do 140 j = lll, k

         col  = zero

         do 110 i = l, k
            if( i.ne.j ) col  = col + abs( aad( i,j ) )
  110    continue

         if( col.eq.zero ) then

            wkd( l ) = j

            if( j.ne.l ) then

               do 120 i = 1, k
                  repl        = aad( i, j )
                  aad( i, j ) = aad( i, l )
                  aad( i, l ) = repl
  120          continue

               do 130 i = l, m
                  repl        = aad( j, i )
                  aad( j, i ) = aad( l, i )
                  aad( l, i ) = repl
  130          continue

            end if

            l  = l + 1
            go to  100

         end if

  140 continue

c                           ** balance the submatrix in rows l through k
      do 150 i = l, k
         wkd( i ) = one
  150 continue

  160 continue
      noconv = .false.

      do 220 i = l, k

         col  = zero
         row  = zero

         do 170 j = l, k

            if( j.ne.i ) then
               col  = col + abs( aad( j,i ) )
               row  = row + abs( aad( i,j ) )
            end if

  170    continue

         f  = one
         g  = row / c5
         h  = col + row

  180    continue
         if( col.lt.g ) then

            f    = f*c5
            col  = col*c6
            go to  180

         end if

         g  = row*c5

  190    continue
         if( col.ge.g ) then

            f    = f / c5
            col  = col / c6
            go to  190

         end if
c                                                ** now balance
         if( ( col + row ) / f.lt.c4*h ) then

            wkd( i ) = wkd( i )*f
            noconv = .true.

            do 200 j = l, m
               aad( i, j ) = aad( i, j ) / f
  200       continue

            do 210 j = 1, k
               aad( j, i ) = aad( j, i )*f
  210       continue

         end if

  220 continue


      if( noconv ) go to  160
c                                   ** is a already in hessenberg form?
      if( k-1 .lt. l+1 ) go to  370

c                                   ** transfer a to a hessenberg form
      do 310 n = l + 1, k - 1

         h  = zero
         wkd( n + m ) = zero
         scale  = zero
c                                                 ** scale column
         do 230 i = n, k
            scale  = scale + abs( aad( i,n - 1 ) )
  230    continue

         if( scale.ne.zero ) then

            do 240 i = k, n, -1
               wkd( i + m ) = aad( i, n - 1 ) / scale
               h  = h + wkd( i + m )**2
  240       continue

            g    = - sign( sqrt( h ), wkd( n + m ) )
            h    = h - wkd( n + m )*g
            wkd( n + m ) = wkd( n + m ) - g
c                                            ** form (i-(u*ut)/h)*a
            do 270 j = n, m

               f  = zero

               do 250 i = k, n, -1
                  f  = f + wkd( i + m )*aad( i, j )
  250          continue

               do 260 i = n, k
                  aad( i, j ) = aad( i, j ) - wkd( i + m )*f / h
  260          continue

  270       continue
c                                    ** form (i-(u*ut)/h)*a*(i-(u*ut)/h)
            do 300 i = 1, k

               f  = zero

               do 280 j = k, n, -1
                  f  = f + wkd( j + m )*aad( i, j )
  280          continue

               do 290 j = n, k
                  aad( i, j ) = aad( i, j ) - wkd( j + m )*f / h
  290          continue

  300       continue

            wkd( n + m ) = scale*wkd( n + m )
            aad( n, n - 1 ) = scale*g

         end if

  310 continue


      do 360 n = k - 2, l, -1

         n1   = n + 1
         n2   = n + 2
         f  = aad( n + 1, n )

         if( f.ne.zero ) then

            f  = f*wkd( n + 1 + m )

            do 320 i = n + 2, k
               wkd( i + m ) = aad( i, n )
  320       continue

            if( n + 1.le.k ) then

               do 350 j = 1, m

                  g  = zero

                  do 330 i = n + 1, k
                     g  = g + wkd( i + m )*evecd( i, j )
  330             continue

                  g  = g / f

                  do 340 i = n + 1, k
                     evecd( i, j ) = evecd( i, j ) + g*wkd( i + m )
  340             continue

  350          continue

            end if

         end if

  360 continue


  370 continue

      n  = 1

      do 390 i = 1, m

         do 380 j = n, m
            rnorm  = rnorm + abs( aad( i,j ) )
  380    continue

         n  = i

         if( i.lt.l .or. i.gt.k ) evald( i ) = aad( i, i )

  390 continue

      n  = k
      t  = zero

c                                      ** search for next eigenvalues
  400 continue
      if( n.lt.l ) go to  550

      in  = 0
      n1  = n - 1
      n2  = n - 2
c                          ** look for single small sub-diagonal element
  410 continue

      do 420 i = l, n

         lb  = n + l - i

         if( lb.eq.l ) go to  430

         s  = abs( aad( lb - 1,lb - 1 ) ) + abs( aad( lb,lb ) )

         if( s.eq.zero ) s  = rnorm

         if( abs( aad( lb, lb-1 ) ).le. tol*s ) go to  430

  420 continue


  430 continue
      x  = aad( n, n )

      if( lb.eq.n ) then
c                                        ** one eigenvalue found
         aad( n, n ) = x + t
         evald( n ) = aad( n, n )
         n  = n1
         go to  400

      end if

      y  = aad( n1, n1 )
      w  = aad( n, n1 )*aad( n1, n )

      if( lb.eq.n1 ) then
c                                        ** two eigenvalues found
         p  = ( y - x )*c2
         q  = p**2 + w
         z  = sqrt( abs( q ) )
         aad( n, n ) = x + t
         x  = aad( n, n )
         aad( n1, n1 ) = y + t
c                                        ** real pair
         z  = p + sign( z, p )
         evald( n1 ) = x + z
         evald( n ) = evald( n1 )

         if( z.ne.zero ) evald( n ) = x - w / z

         x  = aad( n, n1 )
c                                  ** employ scale factor in case
c                                  ** x and z are very small
         r  = sqrt( x*x + z*z )
         p  = x / r
         q  = z / r
c                                             ** row modification
         do 440 j = n1, m
            z  = aad( n1, j )
            aad( n1, j ) = q*z + p*aad( n, j )
            aad( n, j ) = q*aad( n, j ) - p*z
  440    continue
c                                             ** column modification
         do 450 i = 1, n
            z  = aad( i, n1 )
            aad( i, n1 ) = q*z + p*aad( i, n )
            aad( i, n ) = q*aad( i, n ) - p*z
  450    continue
c                                          ** accumulate transformations
         do 460 i = l, k
            z  = evecd( i, n1 )
            evecd( i, n1 ) = q*z + p*evecd( i, n )
            evecd( i, n ) = q*evecd( i, n ) - p*z
  460    continue

         n  = n2
         go to  400

      end if


      if( in.eq.30 ) then

c                    ** no convergence after 30 iterations; set error
c                    ** indicator to the index of the current eigenvalue
         ier  = n
         go to  700

      end if
c                                                  ** form shift
      if( in.eq.10 .or. in.eq.20 ) then

         t  = t + x

         do 470 i = l, n
            aad( i, i ) = aad( i, i ) - x
  470    continue

         s  = abs( aad( n,n1 ) ) + abs( aad( n1,n2 ) )
         x  = c3*s
         y  = x
         w  = -c1*s**2

      end if


      in  = in + 1

c                ** look for two consecutive small sub-diagonal elements

      do 480 j = lb, n2
         i  = n2 + lb - j
         z  = aad( i, i )
         r  = x - z
         s  = y - z
         p  = ( r*s - w ) / aad( i + 1, i ) + aad( i, i + 1 )
         q  = aad( i + 1, i + 1 ) - z - r - s
         r  = aad( i + 2, i + 1 )
         s  = abs( p ) + abs( q ) + abs( r )
         p  = p / s
         q  = q / s
         r  = r / s

         if( i.eq.lb ) go to  490

         uu   = abs( aad( i, i-1 ) )*( abs( q ) + abs( r ) )
         vv   = abs( p ) * ( abs( aad( i-1, i-1 ) ) + abs( z ) +
     &                       abs( aad( i+1, i+1 ) ) )

         if( uu .le. tol*vv ) go to  490

  480 continue

  490 continue
      aad( i+2, i ) = zero

      do 500 j = i + 3, n
         aad( j, j - 2 ) = zero
         aad( j, j - 3 ) = zero
  500 continue

c             ** double qr step involving rows k to n and columns m to n

      do 540 ka = i, n1

         notlas = ka.ne.n1

         if( ka.eq.i ) then

            s  = sign( sqrt( p*p + q*q + r*r ), p )

            if( lb.ne.i ) aad( ka, ka - 1 ) = -aad( ka, ka - 1 )

         else

            p  = aad( ka, ka - 1 )
            q  = aad( ka + 1, ka - 1 )
            r  = zero

            if( notlas ) r  = aad( ka + 2, ka - 1 )

            x  = abs( p ) + abs( q ) + abs( r )

            if( x.eq.zero ) go to  540

            p  = p / x
            q  = q / x
            r  = r / x
            s  = sign( sqrt( p*p + q*q + r*r ), p )
            aad( ka, ka - 1 ) = -s*x

         end if

         p  = p + s
         x  = p / s
         y  = q / s
         z  = r / s
         q  = q / p
         r  = r / p
c                                              ** row modification
         do 510 j = ka, m

            p  = aad( ka, j ) + q*aad( ka + 1, j )

            if( notlas ) then

               p  = p + r*aad( ka + 2, j )
               aad( ka + 2, j ) = aad( ka + 2, j ) - p*z

            end if

            aad( ka + 1, j ) = aad( ka + 1, j ) - p*y
            aad( ka, j ) = aad( ka, j ) - p*x

  510    continue
c                                                 ** column modification
         do 520 ii = 1, min( n, ka + 3 )

            p  = x*aad( ii, ka ) + y*aad( ii, ka + 1 )

            if( notlas ) then

               p  = p + z*aad( ii, ka + 2 )
               aad( ii, ka + 2 ) = aad( ii, ka + 2 ) - p*r

            end if

            aad( ii, ka + 1 ) = aad( ii, ka + 1 ) - p*q
            aad( ii, ka ) = aad( ii, ka ) - p

  520    continue
c                                          ** accumulate transformations
         do 530 ii = l, k

            p  = x*evecd( ii, ka ) + y*evecd( ii, ka + 1 )

            if( notlas ) then

               p  = p + z*evecd( ii, ka + 2 )
               evecd( ii, ka + 2 ) = evecd( ii, ka + 2 ) - p*r

            end if

            evecd( ii, ka + 1 ) = evecd( ii, ka + 1 ) - p*q
            evecd( ii, ka ) = evecd( ii, ka ) - p

  530    continue

  540 continue

      go to  410
c                     ** all evals found, now backsubstitute real vector
  550 continue

      if( rnorm.ne.zero ) then

         do 580 n = m, 1, -1

            n2   = n
            aad( n, n ) = one

            do 570 i = n - 1, 1, -1

               w  = aad( i, i ) - evald( n )

               if( w.eq.zero ) w  = tol*rnorm

               r  = aad( i, n )

               do 560 j = n2, n - 1
                  r  = r + aad( i, j )*aad( j, n )
  560          continue

               aad( i, n ) = -r / w
               n2   = i

  570       continue

  580    continue
c                      ** end backsubstitution vectors of isolated evals
         do 600 i = 1, m

            if( i.lt.l .or. i.gt.k ) then

               do 590 j = i, m
                  evecd( i, j ) = aad( i, j )
  590          continue

            end if

  600    continue
c                                   ** multiply by transformation matrix
         if( k.ne.0 ) then

            do 630 j = m, l, -1

               do 620 i = l, k

                  z  = zero

                  do 610 n = l, min( j, k )
                     z  = z + evecd( i, n )*aad( n, j )
  610             continue

                  evecd( i, j ) = z

  620          continue

  630       continue

         end if

      end if


      do 650 i = l, k

         do 640 j = 1, m
            evecd( i, j ) = evecd( i, j ) * wkd( i )
  640    continue

  650 continue

c                           ** interchange rows if permutations occurred
      do 670 i = l-1, 1, -1

         j  = wkd( i )

         if( i.ne.j ) then

            do 660 n = 1, m
               repl   = evecd( i, n )
               evecd( i, n ) = evecd( j, n )
               evecd( j, n ) = repl
  660       continue

         end if

  670 continue


      do 690 i = k + 1, m

         j  = wkd( i )

         if( i.ne.j ) then

            do 680 n = 1, m
               repl   = evecd( i, n )
               evecd( i, n ) = evecd( j, n )
               evecd( j, n ) = repl
  680       continue

         end if

  690 continue

c                         ** put results into output arrays
  700 continue

      do 720 j = 1, m

         eval( j ) = evald( j )

         do 710 k = 1, m
            evec( j, k ) = evecd( j, k )
  710    continue

  720 continue


      return
      end

      subroutine cmpint( fbeam, gc, kk, layru, ll, lyrcut, mazim, mxcmu,
     &                   mxulv, mxumu, ncut, nn, nstr, plank, ntau,
     &                   taucpr, umu0, utaupr, zz, zplk0, zplk1, uum )

c          calculates the fourier intensity components at the quadrature
c          angles for azimuthal expansion terms (mazim) in eq. sd(2),
c          stwl(6)
c
c
c    i n p u t    v a r i a b l e s:
c
c       kk      :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b)
c
c       gc      :  eigenvectors at polar quadrature angles in eq. sc(1)
c
c       ll      :  constants of integration in eq. sc(1), obtained
c                  by solving scaled version of eq. sc(5);
c                  exponential term of eq. sc(12) not included
c
c       lyrcut  :  logical flag for truncation of computational layer
c
c       mazim   :  order of azimuthal component
c
c       ncut    :  number of computational layer where absorption
c                  optical depth exceeds abscut
c
c       nn      :  order of double-gauss quadrature (nstr/2)
c
c       taucpr  :  cumulative optical depth (delta-m-scaled)
c
c       utaupr  :  optical depths of user output levels in delta-m
c                  coordinates;  equal to utau if no delta-m
c
c       zz      :  beam source vectors in eq. ss(19), stwl(24b)
c
c       zplk0   :  thermal source vectors z0, by solving eq. ss(16),
c                  y-sub-zero in stwl(26ab)
c
c       zplk1   :  thermal source vectors z1, by solving eq. ss(16),
c                  y-sub-one in stwl(26ab)
c
c       (remainder are 'disort' input variables)
c
c
c    o u t p u t   v a r i a b l e s:
c
c       uum     :  fourier components of the intensity in eq. sd(12)
c                    (at polar quadrature angles)
c
c
c    i n t e r n a l   v a r i a b l e s:
c
c       fact    :  exp( - utaupr / umu0 )
c       zint    :  intensity of m=0 case, in eq. sc(1)
c
c   called by- disort
c +--------------------------------------------------------------------

c     .. scalar arguments ..

      logical   lyrcut, plank
      integer   mazim, mxcmu, mxulv, mxumu, ncut, nn, nstr, ntau
      real      fbeam, umu0
c     ..
c     .. array arguments ..

      integer   layru( * )
      real      gc( mxcmu, mxcmu, * ), kk( mxcmu, * ), ll( mxcmu, * ),
     &          taucpr( 0:* ), utaupr( mxulv ), uum( mxumu, mxulv ),
     &          zplk0( mxcmu, * ), zplk1( mxcmu, * ), zz( mxcmu, * )
c     ..
c     .. local scalars ..

      integer   iq, jq, lu, lyu
      real      zint
c     ..
c     .. intrinsic functions ..

      intrinsic exp
c     ..

c                                       ** loop over user levels
      do 40 lu = 1, ntau

         lyu  = layru( lu )

         if( lyrcut .and. lyu.gt.ncut ) go to  40

         do 30 iq = 1, nstr

            zint = 0.0

            do 10 jq = 1, nn
               zint = zint + gc( iq, jq, lyu ) * ll( jq, lyu ) *
     &                       exp( -kk( jq,lyu )*
     &                     ( utaupr( lu ) - taucpr( lyu ) ) )
   10       continue

            do 20 jq = nn + 1, nstr
               zint = zint + gc( iq, jq, lyu ) * ll( jq, lyu ) *
     &                       exp( -kk( jq,lyu )*
     &                     ( utaupr( lu ) - taucpr( lyu-1 ) ) )
   20       continue

            uum( iq, lu ) = zint

            if( fbeam.gt.0.0 ) uum( iq, lu ) = zint +
     &                         zz( iq, lyu )*exp( -utaupr( lu )/umu0 )

            if( plank .and. mazim.eq.0 )
     &          uum( iq, lu ) = uum( iq, lu ) + zplk0( iq,lyu ) +
     &                          zplk1( iq,lyu ) * utaupr( lu )
   30    continue

   40 continue


      return
      end

      subroutine fluxes( cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut,
     &                   maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau,
     &                   pi, prnt, prntu0, ssalb, taucpr, umu0, utau,
     &                   utaupr, xr0, xr1, zz, zplk0, zplk1, dfdt,
     &                   flup, fldn, fldir, rfldir, rfldn, uavg, u0c )

c       calculates the radiative fluxes, mean intensity, and flux
c       derivative with respect to optical depth from the m=0 intensity
c       components (the azimuthally-averaged intensity)
c
c
c    i n p u t     v a r i a b l e s:
c
c       cmu      :  abscissae for gauss quadrature over angle cosine
c
c       cwt      :  weights for gauss quadrature over angle cosine
c
c       gc       :  eigenvectors at polar quadrature angles, sc(1)
c
c       kk       :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b)
c
c       layru    :  layer number of user level utau
c
c       ll       :  constants of integration in eq. sc(1), obtained
c                   by solving scaled version of eq. sc(5);
c                   exponential term of eq. sc(12) not included
c
c       lyrcut   :  logical flag for truncation of comput. layer
c
c       nn       :  order of double-gauss quadrature (nstr/2)
c
c       ncut     :  number of computational layer where absorption
c                   optical depth exceeds abscut
c
c       prntu0   :  true, print azimuthally-averaged intensity at
c                   quadrature angles
c
c       taucpr   :  cumulative optical depth (delta-m-scaled)
c
c       utaupr   :  optical depths of user output levels in delta-m
c                   coordinates;  equal to utau if no delta-m
c
c       xr0      :  expansion of thermal source function in eq. ss(14),
c                   stwl(24c)
c
c       xr1      :  expansion of thermal source function eq. ss(16),
c                   stwl(24c)
c
c       zz       :  beam source vectors in eq. ss(19), stwl(24b)
c
c       zplk0    :  thermal source vectors z0, by solving eq. ss(16),
c                   y0 in stwl(26b)
c
c       zplk1    :  thermal source vectors z1, by solving eq. ss(16),
c                   y1 in stwl(26a)
c
c       (remainder are disort input variables)
c
c
c    o u t p u t     v a r i a b l e s:
c
c       u0c      :  azimuthally averaged intensities
c                   ( at polar quadrature angles )
c
c       (rfldir, rfldn, flup, dfdt, uavg are disort output variables)
c
c
c    i n t e r n a l       v a r i a b l e s:
c
c       dirint   :  direct intensity attenuated
c       fdntot   :  total downward flux (direct + diffuse)
c       fldir    :  direct-beam flux (delta-m scaled)
c       fldn     :  diffuse down-flux (delta-m scaled)
c       fnet     :  net flux (total-down - diffuse-up)
c       fact     :  exp( - utaupr / umu0 )
c       plsorc   :  planck source function (thermal)
c       zint     :  intensity of m = 0 case, in eq. sc(1)
c
c   called by- disort
c   calls- zeroit
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   lyrcut, prntu0
      integer   maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau
      real      fbeam, pi, umu0
c     ..
c     .. array arguments ..

      logical   prnt( * )
      integer   layru( mxulv )
      real      cmu( mxcmu ), cwt( mxcmu ), dfdt( maxulv ),
     &          fldir( mxulv ), fldn( mxulv ), flup( maxulv ),
     &          gc( mxcmu, mxcmu, * ), kk( mxcmu, * ), ll( mxcmu, * ),
     &          rfldir( maxulv ), rfldn( maxulv ), ssalb( * ),
     &          taucpr( 0:* ), u0c( mxcmu, mxulv ), uavg( maxulv ),
     &          utau( maxulv ), utaupr( mxulv ), xr0( * ), xr1( * ),
     &          zplk0( mxcmu, * ), zplk1( mxcmu, * ), zz( mxcmu, * )
c     ..
c     .. local scalars ..

      integer   iq, jq, lu, lyu
      real      ang1, ang2, dirint, fact, fdntot, fnet, plsorc, zint
c     ..
c     .. external subroutines ..

      external  zeroit
c     ..
c     .. intrinsic functions ..

      intrinsic exp
c     ..


      if( prnt( 2 ) ) write ( *, '(//,21x,a,/,2a,/,2a,/)' )
     &    '<----------------------- fluxes ----------------------->',
     &    '   optical  compu    downward    downward    downward     ',
     &    ' upward                    mean      planck   d(net flux)',
     &    '     depth  layer      direct     diffuse       total     ',
     &    'diffuse         net   intensity      source   / d(op dep)'

c                                        ** zero disort output arrays
      call zeroit( u0c, mxulv*mxcmu )
      call zeroit( fldir, mxulv )
      call zeroit( fldn, mxulv )

c                                        ** loop over user levels
      do 80 lu = 1, ntau

         lyu  = layru( lu )

         if( lyrcut .and. lyu.gt.ncut ) then
c                                                ** no radiation reaches
c                                                ** this level
            fdntot = 0.0
            fnet   = 0.0
            plsorc = 0.0
            go to  70

         end if


         if( fbeam.gt.0.0 ) then

            fact         = exp( -utaupr( lu ) / umu0 )
            dirint       = fbeam*fact
            fldir( lu )  = umu0*( fbeam*fact )
            rfldir( lu ) = umu0*fbeam * exp( -utau( lu ) / umu0 )

         else

            dirint       = 0.0
            fldir( lu )  = 0.0
            rfldir( lu ) = 0.0

         end if


         do 30 iq = 1, nn

            zint = 0.0

            do 10 jq = 1, nn
               zint = zint + gc( iq, jq, lyu )*ll( jq, lyu )*
     &                exp( -kk( jq,lyu )*( utaupr( lu ) -
     &                taucpr( lyu ) ) )
   10       continue

            do 20 jq = nn + 1, nstr
               zint = zint + gc( iq, jq, lyu )*ll( jq, lyu )*
     &                exp( -kk( jq,lyu )*( utaupr( lu ) -
     &                taucpr( lyu-1 ) ) )
   20       continue

            u0c( iq, lu ) = zint

            if( fbeam.gt.0.0 ) u0c( iq, lu ) = zint + zz( iq, lyu )*fact

            u0c( iq, lu ) = u0c( iq, lu ) + zplk0( iq,lyu ) +
     &                      zplk1( iq,lyu )*utaupr( lu )
            uavg( lu ) = uavg( lu ) + cwt( nn + 1 - iq )*u0c( iq, lu )
            fldn( lu ) = fldn( lu ) + cwt( nn + 1 - iq )*
     &                   cmu( nn + 1 - iq )*u0c( iq, lu )
   30    continue


         do 60 iq = nn + 1, nstr

            zint = 0.0

            do 40 jq = 1, nn
               zint = zint + gc( iq, jq, lyu )*ll( jq, lyu )*
     &                exp( -kk( jq,lyu )*( utaupr( lu ) -
     &                taucpr( lyu ) ) )
   40       continue

            do 50 jq = nn + 1, nstr
               zint = zint + gc( iq, jq, lyu )*ll( jq, lyu )*
     &                exp( -kk( jq,lyu )*( utaupr( lu ) -
     &                taucpr( lyu-1 ) ) )
   50       continue

            u0c( iq, lu ) = zint

            if( fbeam.gt.0.0 ) u0c( iq, lu ) = zint + zz( iq, lyu )*fact

            u0c( iq, lu ) = u0c( iq, lu ) + zplk0( iq,lyu ) +
     &                      zplk1( iq,lyu )*utaupr( lu )
            uavg( lu ) = uavg( lu ) + cwt( iq - nn )*u0c( iq, lu )
            flup( lu ) = flup( lu ) + cwt( iq - nn )*cmu( iq - nn )*
     &                   u0c( iq, lu )
   60    continue


         flup( lu )  = 2.*pi*flup( lu )
         fldn( lu )  = 2.*pi*fldn( lu )
         fdntot      = fldn( lu ) + fldir( lu )
         fnet        = fdntot - flup( lu )
         rfldn( lu ) = fdntot - rfldir( lu )
         uavg( lu )  = ( 2.*pi*uavg( lu ) + dirint ) / ( 4.*pi )
         plsorc      = xr0( lyu ) + xr1( lyu )*utaupr( lu )
         dfdt( lu )  = ( 1. - ssalb( lyu ) ) * 4.*pi *
     &                 ( uavg( lu ) - plsorc )

   70    continue
         if( prnt( 2 ) ) write ( *, '(f10.4,i7,1p,7e12.3,e14.3)' )
     &       utau( lu ), lyu, rfldir( lu ), rfldn( lu ), fdntot,
     &       flup( lu ), fnet, uavg( lu ), plsorc, dfdt( lu )

   80 continue


      if( prntu0 ) then

         write ( *, '(//,2a)' ) ' ******** azimuthally averaged ',
     &     'intensities ( at polar quadrature angles ) *******'

         do 100 lu = 1, ntau

            write ( *, '(/,a,f10.4,//,2a)' )
     &        ' optical depth =', utau( lu ),
     &        '     angle (deg)   cos(angle)     intensity',
     &        '     angle (deg)   cos(angle)     intensity'

            do 90 iq = 1, nn
               ang1 = ( 180./pi )*acos( cmu( 2 *nn-iq+1 ) )
               ang2 = ( 180./pi )*acos( cmu( iq ) )
               write ( *, '(2(0p,f16.4,f13.5,1p,e14.3))' )
     &           ang1, cmu(2*nn-iq+1), u0c(iq,lu),
     &           ang2, cmu(iq),        u0c(iq+nn,lu)
   90       continue

  100    continue

      end if


      return
      end

      subroutine intcor( dither, fbeam, flyr, layru, lyrcut, maxmom,
     &                   maxulv, maxumu, nmom, ncut, nphi, nstr, ntau,
     &                   numu, oprim, phasa, phast, phasm, phirad, pi,
     &                   rpd, pmom, ssalb, dtauc, tauc, taucpr, umu,
     &                   umu0, utau, utaupr, uu )

c       corrects intensity field by using nakajima-tanaka algorithm
c       (1988). for more details, see section 3.6 of stwl nasa report.

c                i n p u t   v a r i a b l e s
c
c       dither  10 times machine precision
c
c       dtauc   computational-layer optical depths
c
c       fbeam   incident beam radiation at top
c
c       flyr    separated fraction in delta-m method
c
c       layru   index of utau in multi-layered system
c
c       lyrcut  logical flag for truncation of computational layer
c
c       nmom    number of phase function legendre coefficients supplied
c
c       ncut    total number of computational layers considered
c
c       nphi    number of user azimuthal angles
c
c       nstr    number of polar quadrature angles
c
c       ntau    number of user-defined optical depths
c
c       numu    number of user polar angles
c
c       oprim   delta-m-scaled single-scatter albedo
c
c       phirad  azimuthal angles in radians
c
c       pmom    phase function legendre coefficients (k, lc)
c                   k = 0 to nmom, lc = 1 to nlyr with pmom(0,lc)=1
c
c       rpd     pi/180
c
c       ssalb   single scattering albedo at computational layers
c
c       tauc    optical thickness at computational levels
c
c       taucpr  delta-m-scaled optical thickness
c
c       umu     cosine of emergent angle
c
c       umu0    cosine of incident zenith angle
c
c       utau    user defined optical depths
c
c       utaupr  delta-m-scaled version of utau
c
c                o u t p u t   v a r i a b l e s
c
c       uu      corrected intensity field; uu(iu,lu,j)
c                         iu=1,numu; lu=1,ntau; j=1,nphi
c
c                i n t e r n a l   v a r i a b l e s
c
c       ctheta  cosine of scattering angle
c       dtheta  angle (degrees) to define aureole region as
c                    direction of beam source +/- dtheta
c       phasa   actual (exact) phase function
c       phasm   delta-m-scaled phase function
c       phast   phase function used in tms correction; actual phase
c                    function divided by (1-flyr*ssalb)
c       pl      ordinary legendre polynomial of degree l, p-sub-l
c       plm1    ordinary legendre polynomial of degree l-1, p-sub-(l-1)
c       plm2    ordinary legendre polynomial of degree l-2, p-sub-(l-2)
c       theta0  incident zenith angle (degrees)
c       thetap  emergent angle (degrees)
c       ussndm  single-scattered intensity computed by using exact
c                   phase function and scaled optical depth
c                   (first term in stwl(68a))
c       ussp    single-scattered intensity from delta-m method
c                   (second term in stwl(68a))
c       duims   intensity correction term from ims method
c                   (delta-i-sub-ims in stwl(a.19))
c
c   called by- disort
c   calls- sinsca, secsca

c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   lyrcut
      integer   maxmom, maxulv, maxumu, ncut, nmom, nphi, nstr, ntau,
     &          numu
      real      dither, fbeam, pi, rpd, umu0
c     ..
c     .. array arguments ..

      integer   layru( * )
      real      dtauc( * ), flyr( * ), oprim( * ), phasa( * ),
     &          phast( * ), phasm( * ), phirad( * ),
     &          pmom( 0:maxmom, * ), ssalb( * ), tauc( 0:* ),
     &          taucpr( 0:* ), umu( * ), utau( * ), utaupr( * ),
     &          uu( maxumu, maxulv, * )
c     ..
c     .. local scalars ..

      integer   iu, jp, k, lc, ltau, lu
      real      ctheta, dtheta, duims, pl, plm1, plm2, theta0, thetap,
     &          ussndm, ussp
c     ..
c     .. external functions ..

      real      secsca, sinsca
      external  secsca, sinsca
c     ..
c     .. intrinsic functions ..

      intrinsic abs, acos, cos, sqrt
c     ..


      dtheta = 10.

c                                ** start loop over zenith angles

      do 110 iu = 1, numu

         if( umu( iu ).lt.0. ) then

c                                ** calculate zenith angles of icident
c                                ** and emerging directions

            theta0 = acos( -umu0 ) / rpd
            thetap = acos( umu( iu ) ) / rpd

         end if

c                                ** start loop over azimuth angles

         do 100 jp = 1, nphi

c                                ** calculate cosine of scattering
c                                ** angle, eq. stwl(4)

            ctheta = -umu0*umu( iu ) + sqrt( ( 1.-umu0**2 )*
     &               ( 1.-umu( iu )**2 ) )*cos( phirad( jp ) )

c                                ** initialize phase function
            do 10 lc = 1, ncut

               phasa( lc ) = 1.
               phasm( lc ) = 1.

   10       continue
c                                ** initialize legendre poly. recurrence
            plm1 = 1.
            plm2 = 0.

            do 40 k = 1, nmom
c                                ** calculate legendre polynomial of
c                                ** p-sub-l by upward recurrence

               pl   = ( ( 2 *k-1 )*ctheta*plm1 - ( k-1 )*plm2 ) / k
               plm2 = plm1
               plm1 = pl
c                                ** calculate actual phase function
               do 20 lc = 1, ncut

                  phasa( lc ) = phasa( lc ) +
     &                          ( 2*k + 1 )*pl*pmom( k, lc )

   20          continue

c                                ** calculate delta-m transformed
c                                ** phase function
               if( k.le.nstr - 1 ) then

                  do 30 lc = 1, ncut

                     phasm( lc ) = phasm( lc ) + ( 2*k + 1 ) * pl *
     &                             ( pmom( k,lc ) - flyr( lc ) ) /
     &                             ( 1. - flyr( lc ) )
   30             continue

               end if

   40       continue


c                                ** apply tms method, eq. stwl(68)
            do 70 lc = 1, ncut

               phast( lc ) = phasa(lc) / ( 1. - flyr(lc) * ssalb(lc) )

   70       continue

            do 80 lu = 1, ntau

               if( .not.lyrcut .or. layru( lu ).lt.ncut ) then

                   ussndm  = sinsca( dither, layru( lu ), ncut, phast,
     &                               ssalb, taucpr, umu( iu ), umu0,
     &                               utaupr( lu ), fbeam, pi )

                   ussp    = sinsca( dither, layru( lu ), ncut, phasm,
     &                               oprim, taucpr, umu( iu ), umu0,
     &                               utaupr( lu ), fbeam, pi )

                   uu( iu, lu, jp ) = uu( iu, lu, jp ) + ussndm - ussp

               end if

   80       continue

            if( umu(iu).lt.0. .and. abs( theta0-thetap ).le.dtheta) then

c                                ** emerging direction is in the aureole
c                                ** (theta0 +/- dtheta). apply ims
c                                ** method for correction of secondary
c                                ** scattering below top level.

               ltau = 1

               if( utau( 1 ).le.dither ) ltau = 2

               do 90 lu = ltau, ntau

                  if( .not.lyrcut .or. layru( lu ).lt.ncut ) then

                      duims = secsca( ctheta, flyr, layru( lu ), maxmom,
     &                                nmom, nstr, pmom, ssalb, dtauc,
     &                                tauc, umu( iu ), umu0, utau( lu ),
     &                                fbeam, pi )

                      uu( iu, lu, jp ) = uu( iu, lu, jp ) - duims

                  end if

   90          continue

            end if
c                                ** end loop over azimuth angles
  100    continue

c                                ** end loop over zenith angles
  110 continue


      return
      end

      real function  secsca( ctheta, flyr, layru, maxmom, nmom, nstr,
     &                       pmom, ssalb, dtauc, tauc, umu, umu0, utau,
     &                       fbeam, pi )

c          calculates secondary scattered intensity of eq. stwl (a7)

c                i n p u t   v a r i a b l e s

c        ctheta  cosine of scattering angle
c
c        dtauc   computational-layer optical depths
c
c        flyr    separated fraction f in delta-m method
c
c        layru   index of utau in multi-layered system
c
c        maxmom  maximum number of phase function moment coefficients
c
c        nmom    number of phase function legendre coefficients supplied
c
c        nstr    number of polar quadrature angles
c
c        pmom    phase function legendre coefficients (k, lc)
c                k = 0 to nmom, lc = 1 to nlyr, with pmom(0,lc)=1
c
c        ssalb   single scattering albedo of computational layers
c
c        tauc    cumulative optical depth at computational layers
c
c        umu     cosine of emergent angle
c
c        umu0    cosine of incident zenith angle
c
c        utau    user defined optical depth for output intensity
c
c        fbeam   incident beam radiation at top
c
c        pi       3.1415...
c
c   local variables
c
c        pspike  2*p"-p"**2, where p" is the residual phase function
c        wbar    mean value of single scattering albedo
c        fbar    mean value of separated fraction f
c        dtau    layer optical depth
c        stau    sum of layer optical depths between top of atmopshere
c                and layer layru
c
c   called by- intcor
c   calls- xifunc
c +-------------------------------------------------------------------+

c     .. scalar arguments ..
      integer   layru, maxmom, nmom, nstr
      real      ctheta, fbeam, pi, umu, umu0, utau
c     ..
c     .. array arguments ..
      real      dtauc( * ), flyr( * ), pmom( 0:maxmom, * ), ssalb( * ),
     &          tauc( 0:* )
c     ..
c     .. local scalars ..
      integer   k, lyr
      real      dtau, fbar, gbar, pl, plm1, plm2, pspike, stau, umu0p,
     &          wbar, zero
c     ..
c     .. external functions ..
      real      xifunc
      external  xifunc
c     ..

      zero = 1e-4

c                          ** calculate vertically averaged value of
c                          ** single scattering albedo and separated
c                          ** fraction f, eq. stwl (a.15)

      dtau = utau - tauc( layru - 1 )
      wbar = ssalb( layru ) * dtau
      fbar = flyr( layru ) * wbar
      stau = dtau

      do 10 lyr = 1, layru - 1

         wbar = wbar + ssalb( lyr ) * dtauc( lyr )
         fbar = fbar + ssalb( lyr ) * dtauc( lyr ) * flyr( lyr )
         stau = stau + dtauc( lyr )

   10 continue

      if( wbar.le.zero .or.
     &    fbar.le.zero .or. stau.le.zero .or.fbeam.le.zero ) then

          secsca = 0.0
          return

      end if

      fbar  = fbar / wbar
      wbar  = wbar / stau


c                          ** calculate pspike=(2p"-p"**2)
      pspike = 1.
      gbar   = 1.
      plm1    = 1.
      plm2    = 0.
c                                   ** pspike for l<=2n-1
      do 20 k = 1, nstr - 1

         pl   = ( ( 2 *k-1 )*ctheta*plm1 - ( k-1 )*plm2 ) / k
         plm2  = plm1
         plm1  = pl

         pspike = pspike + ( 2.*gbar - gbar**2 )*( 2*k + 1 )*pl

   20 continue
c                                   ** pspike for l>2n-1
      do 40 k = nstr, nmom

         pl   = ( ( 2 *k-1 )*ctheta*plm1 - ( k-1 )*plm2 ) / k
         plm2  = plm1
         plm1  = pl

         dtau = utau - tauc( layru - 1 )

         gbar = pmom( k, layru ) * ssalb( layru ) * dtau

         do 30 lyr = 1, layru - 1
            gbar = gbar + pmom( k, lyr ) * ssalb( lyr ) * dtauc( lyr )
   30    continue

         if( fbar*wbar*stau .le. zero ) then
            gbar   = 0.0
         else
            gbar   = gbar / ( fbar*wbar*stau )
         end if

         pspike = pspike + ( 2.*gbar - gbar**2 )*( 2*k + 1 )*pl

   40 continue

      umu0p = umu0 / ( 1. - fbar*wbar )

c                              ** calculate ims correction term,
c                              ** eq. stwl (a.13)

      secsca = fbeam / ( 4.*pi ) * ( fbar*wbar )**2 / ( 1.-fbar*wbar ) *
     &         pspike * xifunc( -umu, umu0p, umu0p, utau )


      return
      end

      subroutine setdis( cmu, cwt, deltam, dtauc, dtaucp, expbea, fbeam,
     &                   flyr, gl, ibcnd, layru, lyrcut, maxmom, maxumu,
     &                   mxcmu, ncut, nlyr, ntau, nn, nstr, plank, numu,
     &                   onlyfl, corint, oprim, pmom, ssalb, tauc,
     &                   taucpr, utau, utaupr, umu, umu0, usrtau,
     &                   usrang )

c          perform miscellaneous setting-up operations
c
c    input :  all are disort input variables (see doc file)
c
c
c    o u t p u t     v a r i a b l e s:
c
c       ntau,utau   if usrtau = false (defined in disort.doc)
c       numu,umu    if usrang = false (defined in disort.doc)
c
c       cmu,cwt     computational polar angles and
c                   corresponding quadrature weights
c
c       expbea      transmission of direct beam
c
c       flyr        separated fraction in delta-m method
c
c       gl          phase function legendre coefficients multiplied
c                   by (2l+1) and single-scatter albedo
c
c       layru       computational layer in which utau falls
c
c       lyrcut      flag as to whether radiation will be zeroed
c                   below layer ncut
c
c       ncut        computational layer where absorption
c                   optical depth first exceeds  abscut
c
c       nn          nstr / 2
c
c       oprim       delta-m-scaled single-scatter albedo
c
c       taucpr      delta-m-scaled optical depth
c
c       utaupr      delta-m-scaled version of  utau
c
c   called by- disort
c   calls- qgausn, errmsg
c ---------------------------------------------------------------------

c     .. scalar arguments ..

      logical   corint, deltam, lyrcut, onlyfl, plank, usrang, usrtau
      logical   tf
      integer   ibcnd, maxmom, maxumu, mxcmu, ncut, nlyr, nn, nstr,
     &          ntau, numu
      real      fbeam, umu0
c     ..
c     .. array arguments ..

      integer   layru( * )
      real      cmu( mxcmu ), cwt( mxcmu ), dtauc( * ), dtaucp( * ),
     &          expbea( 0:* ), flyr( * ), gl( 0:mxcmu, * ), oprim( * ),
     &          pmom( 0:maxmom, * ), ssalb( * ), tauc( 0:* ),
     &          taucpr( 0:* ), umu( maxumu ), utau( * ), utaupr( * )
c     ..
c     .. local scalars ..

      integer   iq, iu, k, lc, lu
      real      abscut, abstau, f, yessct
c     ..
c     .. external subroutines ..

      external  errmsg, qgausn
c     ..
c     .. intrinsic functions ..

      intrinsic abs, exp
c     ..
      data      abscut / 10. /


      if( .not.usrtau ) then
c                              ** set output levels at computational
c                              ** layer boundaries
         ntau  = nlyr + 1

         do 10 lc = 0, ntau - 1
            utau( lc + 1 ) = tauc( lc )
   10    continue

      end if
c                        ** apply delta-m scaling and move description
c                        ** of computational layers to local variables
      expbea( 0 ) = 1.0
      taucpr( 0 ) = 0.0
      abstau      = 0.0
      yessct      = 0.0

      do 40 lc = 1, nlyr

         yessct = yessct + ssalb( lc )

         pmom( 0, lc ) = 1.0

         if( abstau.lt.abscut ) ncut  = lc

         abstau = abstau + ( 1. - ssalb( lc ) )*dtauc( lc )

         if( .not.deltam ) then

            oprim( lc )  = ssalb( lc )
            dtaucp( lc ) = dtauc( lc )
            taucpr( lc ) = tauc( lc )

            do 20 k = 0, nstr - 1
               gl( k, lc ) = ( 2*k + 1 )*oprim( lc )*pmom( k, lc )
   20       continue

            f  = 0.0


         else
c                                    ** do delta-m transformation

            f  = pmom( nstr, lc )
            oprim( lc )  = ssalb( lc )*( 1. - f ) / ( 1. - f*ssalb(lc) )
            dtaucp( lc ) = ( 1. - f*ssalb( lc ) )*dtauc( lc )
            taucpr( lc ) = taucpr( lc - 1 ) + dtaucp( lc )

            do 30 k = 0, nstr - 1
               gl( k, lc ) = ( 2*k + 1 )*oprim( lc )*
     &                       ( pmom( k,lc ) - f ) / ( 1. - f )
   30       continue

         end if

         flyr( lc ) = f
         expbea( lc ) = 0.0

         if( fbeam.gt.0.0 ) expbea( lc ) = exp( -taucpr( lc )/umu0 )

   40 continue
c                      ** if no thermal emission, cut off medium below
c                      ** absorption optical depth = abscut ( note that
c                      ** delta-m transformation leaves absorption
c                      ** optical depth invariant ).  not worth the
c                      ** trouble for one-layer problems, though.
      lyrcut = .false.

      if( abstau.ge.abscut .and. .not.plank .and. ibcnd.ne.1 .and.
     &    nlyr.gt.1 ) lyrcut = .true.

      if( .not.lyrcut ) ncut = nlyr

c                             ** set arrays defining location of user
c                             ** output levels within delta-m-scaled
c                             ** computational mesh
      do 70 lu = 1, ntau

         do 50 lc = 1, nlyr

            if( utau( lu ).ge.tauc( lc-1 ) .and.
     &          utau( lu ).le.tauc( lc ) ) go to  60

   50    continue
         lc   = nlyr

   60    continue
         utaupr( lu ) = utau( lu )
         if( deltam ) utaupr( lu ) = taucpr( lc - 1 ) +
     &                               ( 1. - ssalb( lc )*flyr( lc ) )*
     &                               ( utau( lu ) - tauc( lc-1 ) )
         layru( lu ) = lc

   70 continue
c                      ** calculate computational polar angle cosines
c                      ** and associated quadrature weights for gaussian
c                      ** quadrature on the interval (0,1) (upward)
      nn   = nstr / 2

      call qgausn( nn, cmu, cwt )
c                                  ** downward (neg) angles and weights
      do 80 iq = 1, nn
         cmu( iq + nn ) = -cmu( iq )
         cwt( iq + nn ) = cwt( iq )
   80 continue


      if( fbeam.gt.0.0 ) then
c                     ** compare beam angle to comput. angles
         tf = .true.
         do 90 iq = 1, nn

            if( abs( umu0-cmu( iq ) )/umu0.lt.1.e-4 ) call errmsg(
     &          'setdis--beam angle=computational angle; change nstr',
     &          tf )

   90    continue

      end if


      if( .not.usrang .or. ( onlyfl.and.maxumu.ge.nstr ) ) then

c                                   ** set output polar angles to
c                                   ** computational polar angles
         numu = nstr

         do 100 iu = 1, nn
            umu( iu ) = -cmu( nn + 1 - iu )
  100    continue

         do 110 iu = nn + 1, nstr
            umu( iu ) = cmu( iu - nn )
  110    continue

      end if


      if( usrang .and. ibcnd.eq.1 ) then

c                               ** shift positive user angle cosines to
c                               ** upper locations and put negatives
c                               ** in lower locations
         do 120 iu = 1, numu
            umu( iu + numu ) = umu( iu )
  120    continue

         do 130 iu = 1, numu
            umu( iu ) = -umu( 2*numu + 1 - iu )
  130    continue

         numu = 2*numu

      end if

c                               ** turn off intensity correction when
c                               ** only fluxes are calculated, there
c                               ** is no beam source, no scattering,
c                               ** or delta-m transformation is not
c                               ** applied
c
      if( onlyfl .or. fbeam.eq.0.0 .or. yessct.eq.0.0 .or.
     &   .not.deltam )  corint = .false.


      return
      end

      subroutine setmtx( bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk,
     &                   lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut,
     &                   nnlyri, nn, nstr, taucpr, wk )

c        calculate coefficient matrix for the set of equations
c        obtained from the boundary conditions and the continuity-
c        of-intensity-at-layer-interface equations;  store in the
c        special banded-matrix format required by linpack routines
c
c
c    i n p u t      v a r i a b l e s:
c
c       bdr      :  surface bidirectional reflectivity
c
c       cmu,cwt     abscissae, weights for gauss quadrature
c                   over angle cosine
c
c       delm0    :  kronecker delta, delta-sub-m0
c
c       gc       :  eigenvectors at polar quadrature angles, sc(1)
c
c       kk       :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b)
c
c       lyrcut   :  logical flag for truncation of computational layers
c
c       nn       :  number of streams in a hemisphere (nstr/2)
c
c       ncut     :  total number of computational layers considered
c
c       taucpr   :  cumulative optical depth (delta-m-scaled)
c
c       (remainder are disort input variables)
c
c
c   o u t p u t     v a r i a b l e s:
c
c       cband    :  left-hand side matrix of linear system eq. sc(5),
c                   scaled by eq. sc(12); in banded form required
c                   by linpack solution routines
c
c       ncol     :  number of columns in cband
c
c
c   i n t e r n a l    v a r i a b l e s:
c
c       irow     :  points to row in cband
c       jcol     :  points to position in layer block
c       lda      :  row dimension of cband
c       ncd      :  number of diagonals below or above main diagonal
c       nshift   :  for positioning number of rows in band storage
c       wk       :  temporary storage for exp evaluations
c
c
c   band storage
c
c      linpack requires band matrices to be input in a special
c      form where the elements of each diagonal are moved up or
c      down (in their column) so that each diagonal becomes a row.
c      (the column locations of diagonal elements are unchanged.)
c
c      example:  if the original matrix is
c
c          11 12 13  0  0  0
c          21 22 23 24  0  0
c           0 32 33 34 35  0
c           0  0 43 44 45 46
c           0  0  0 54 55 56
c           0  0  0  0 65 66
c
c      then its linpack input form would be:
c
c           *  *  *  +  +  +  , * = not used
c           *  * 13 24 35 46  , + = used for pivoting
c           * 12 23 34 45 56
c          11 22 33 44 55 66
c          21 32 43 54 65  *
c
c      if a is a band matrix, the following program segment
c      will convert it to the form (abd) required by linpack
c      band-matrix routines:
c
c               n  = (column dimension of a, abd)
c               ml = (band width below the diagonal)
c               mu = (band width above the diagonal)
c               m = ml + mu + 1
c               do j = 1, n
c                  i1 = max(1, j-mu)
c                  i2 = min(n, j+ml)
c                  do i = i1, i2
c                     k = i - j + m
c                     abd(k,j) = a(i,j)
c                  end do
c               end do
c
c      this uses rows  ml+1  through  2*ml+mu+1  of abd.
c      the total number of rows needed in abd is  2*ml+mu+1 .
c      in the example above, n = 6, ml = 1, mu = 2, and the
c      row dimension of abd must be >= 5.
c
c
c   called by- disort, albtrn
c   calls- zeroit
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   lamber, lyrcut
      integer   mi, mi9m2, mxcmu, ncol, ncut, nn, nnlyri, nstr
      real      delm0
c     ..
c     .. array arguments ..

      real      bdr( mi, 0:mi ), cband( mi9m2, nnlyri ), cmu( mxcmu ),
     &          cwt( mxcmu ), dtaucp( * ), gc( mxcmu, mxcmu, * ),
     &          kk( mxcmu, * ), taucpr( 0:* ), wk( mxcmu )
c     ..
c     .. local scalars ..

      integer   iq, irow, jcol, jq, k, lc, lda, ncd, nncol, nshift
      real      expa, sum
c     ..
c     .. external subroutines ..

      external  zeroit
c     ..
c     .. intrinsic functions ..

      intrinsic exp
c     ..


      call zeroit( cband, mi9m2*nnlyri )

      ncd    = 3*nn - 1
      lda    = 3*ncd + 1
      nshift = lda - 2*nstr + 1
      ncol   = 0
c                         ** use continuity conditions of eq. stwj(17)
c                         ** to form coefficient matrix in stwj(20);
c                         ** employ scaling transformation stwj(22)
      do 60 lc = 1, ncut

         do 10 iq = 1, nn
            wk( iq ) = exp( kk( iq,lc )*dtaucp( lc ) )
   10    continue

         jcol  = 0

         do 30 iq = 1, nn

            ncol  = ncol + 1
            irow  = nshift - jcol

            do 20 jq = 1, nstr
               cband( irow + nstr, ncol ) =   gc( jq, iq, lc )
               cband( irow, ncol )        = - gc( jq, iq, lc )*wk( iq )
               irow  = irow + 1
   20       continue

            jcol  = jcol + 1

   30    continue


         do 50 iq = nn + 1, nstr

            ncol  = ncol + 1
            irow  = nshift - jcol

            do 40 jq = 1, nstr
               cband( irow + nstr, ncol ) =   gc( jq, iq, lc )*
     &                                          wk( nstr + 1 - iq )
               cband( irow, ncol )        = - gc( jq, iq, lc )
               irow  = irow + 1
   40       continue

            jcol  = jcol + 1

   50    continue

   60 continue
c                  ** use top boundary condition of stwj(20a) for
c                  ** first layer
      jcol  = 0

      do 80 iq = 1, nn

         expa  = exp( kk( iq,1 )*taucpr( 1 ) )
         irow  = nshift - jcol + nn

         do 70 jq = nn, 1, -1
            cband( irow, jcol + 1 ) = gc( jq, iq, 1 )*expa
            irow  = irow + 1
   70    continue

         jcol  = jcol + 1

   80 continue


      do 100 iq = nn + 1, nstr

         irow  = nshift - jcol + nn

         do 90 jq = nn, 1, -1
            cband( irow, jcol + 1 ) = gc( jq, iq, 1 )
            irow  = irow + 1
   90    continue

         jcol  = jcol + 1

  100 continue
c                           ** use bottom boundary condition of
c                           ** stwj(20c) for last layer

      nncol = ncol - nstr
      jcol  = 0

      do 130 iq = 1, nn

         nncol  = nncol + 1
         irow   = nshift - jcol + nstr

         do 120 jq = nn + 1, nstr

            if( lyrcut .or. ( lamber .and. delm0.eq.0 ) ) then

c                          ** no azimuthal-dependent intensity if lam-
c                          ** bert surface; no intensity component if
c                          ** truncated bottom layer

               cband( irow, nncol ) = gc( jq, iq, ncut )

            else

               sum  = 0.0

               do 110 k = 1, nn
                  sum  = sum + cwt( k )*cmu( k )*bdr( jq - nn, k )*
     &                     gc( nn + 1 - k, iq, ncut )
  110          continue

               cband( irow, nncol ) = gc( jq, iq, ncut ) -
     &                                ( 1.+ delm0 )*sum
            end if

            irow  = irow + 1

  120    continue

         jcol  = jcol + 1

  130 continue


      do 160 iq = nn + 1, nstr

         nncol  = nncol + 1
         irow   = nshift - jcol + nstr
         expa   = wk( nstr + 1 - iq )

         do 150 jq = nn + 1, nstr

            if( lyrcut .or. ( lamber .and. delm0.eq.0 ) ) then

               cband( irow, nncol ) = gc( jq, iq, ncut )*expa

            else

               sum  = 0.0

               do 140 k = 1, nn
                  sum  = sum + cwt( k )*cmu( k )*bdr( jq - nn, k )*
     &                         gc( nn + 1 - k, iq, ncut )
  140          continue

               cband( irow, nncol ) = ( gc( jq,iq,ncut ) -
     &                                ( 1.+ delm0 )*sum )*expa
            end if

            irow  = irow + 1

  150    continue

         jcol  = jcol + 1

  160 continue


      return
      end

      real function  sinsca( dither, layru, nlyr, phase, omega, tau,
     &                       umu, umu0, utau, fbeam, pi )

c        calculates single-scattered intensity from eqs. stwl (65b,d,e)

c                i n p u t   v a r i a b l e s

c        dither   10 times machine precision
c
c        layru    index of utau in multi-layered system
c
c        nlyr     number of sublayers
c
c        phase    phase functions of sublayers
c
c        omega    single scattering albedos of sublayers
c
c        tau      optical thicknesses of sublayers
c
c        umu      cosine of emergent angle
c
c        umu0     cosine of incident zenith angle
c
c        utau     user defined optical depth for output intensity
c
c        fbeam   incident beam radiation at top
c
c        pi       3.1415...
c
c   called by- intcor
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      integer   layru, nlyr
      real      dither, fbeam, pi, umu, umu0, utau
c     ..
c     .. array arguments ..

      real      omega( * ), phase( * ), tau( 0:* )
c     ..
c     .. local scalars ..

      integer   lyr
      real      exp0, exp1
c     ..
c     .. intrinsic functions ..

      intrinsic abs, exp
c     ..


      sinsca = 0.
      exp0 = exp( -utau/umu0 )

      if( abs( umu+umu0 ).le.dither ) then

c                                 ** calculate downward intensity when
c                                 ** umu=umu0, eq. stwl (65e)

         do 10 lyr = 1, layru - 1
            sinsca = sinsca + omega( lyr ) * phase( lyr ) *
     &               ( tau( lyr ) - tau( lyr-1 ) )
   10    continue

         sinsca = fbeam / ( 4.*pi * umu0 ) * exp0 * ( sinsca +
     &            omega( layru )*phase( layru )*( utau-tau(layru-1) ) )

         return

      end if


      if( umu.gt.0. ) then
c                                 ** upward intensity, eq. stwl (65b)
         do 20 lyr = layru, nlyr

            exp1 = exp( -( ( tau( lyr )-utau )/umu + tau( lyr )/umu0 ) )
            sinsca = sinsca + omega( lyr )*phase( lyr )*( exp0 - exp1 )
            exp0 = exp1

   20    continue

      else
c                                 ** downward intensity, eq. stwl (65d)
         do 30 lyr = layru, 1, -1

            exp1 = exp( -( ( tau(lyr-1)-utau )/umu + tau(lyr-1)/umu0 ) )
            sinsca = sinsca + omega( lyr )*phase( lyr )*( exp0 - exp1 )
            exp0 = exp1

   30    continue

      end if

      sinsca = fbeam / ( 4.*pi * ( 1. + umu/umu0 ) ) * sinsca


      return
      end

      subroutine soleig( amb, apb, array, cmu, cwt, gl, mi, mazim,
     &                   mxcmu, nn, nstr, ylmc, cc, evecc, eval, kk, gc,
     &                   aad, eveccd, evald, wkd )

c         solves eigenvalue/vector problem necessary to construct
c         homogeneous part of discrete ordinate solution; stwj(8b),
c         stwl(23f)
c         ** note ** eigenvalue problem is degenerate when single
c                    scattering albedo = 1;  present way of doing it
c                    seems numerically more stable than alternative
c                    methods that we tried
c
c
c   i n p u t     v a r i a b l e s:
c
c       gl     :  delta-m scaled legendre coefficients of phase function
c                 (including factors 2l+1 and single-scatter albedo)
c
c       cmu    :  computational polar angle cosines
c
c       cwt    :  weights for quadrature over polar angle cosine
c
c       mazim  :  order of azimuthal component
c
c       nn     :  half the total number of streams
c
c       ylmc   :  normalized associated legendre polynomial
c                 at the quadrature angles cmu
c
c       (remainder are disort input variables)
c
c
c   o u t p u t    v a r i a b l e s:
c
c       cc     :  c-sub-ij in eq. ss(5); needed in ss(15&18)
c
c       eval   :  nn eigenvalues of eq. ss(12), stwl(23f) on return
c                 from asymtx but then square roots taken
c
c       evecc  :  nn eigenvectors  (g+) - (g-)  on return
c                 from asymtx ( column j corresponds to eval(j) )
c                 but then  (g+) + (g-)  is calculated from ss(10),
c                 g+  and  g-  are separated, and  g+  is stacked on
c                 top of  g-  to form nstr eigenvectors of ss(7)
c
c       gc     :  permanent storage for all nstr eigenvectors, but
c                 in an order corresponding to kk
c
c       kk     :  permanent storage for all nstr eigenvalues of ss(7),
c                 but re-ordered with negative values first ( square
c                 roots of eval taken and negatives added )
c
c
c   i n t e r n a l   v a r i a b l e s:
c
c       amb,apb :  matrices (alpha-beta), (alpha+beta) in reduced
c                    eigenvalue problem
c       array   :  complete coefficient matrix of reduced eigenvalue
c                    problem: (alfa+beta)*(alfa-beta)
c       gpplgm  :  (g+) + (g-) (cf. eqs. ss(10-11))
c       gpmigm  :  (g+) - (g-) (cf. eqs. ss(10-11))
c       wkd     :  scratch array required by asymtx
c
c   called by- disort, albtrn
c   calls- asymtx, errmsg
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   tf
      integer   mazim, mi, mxcmu, nn, nstr
c     ..
c     .. array arguments ..

      real      amb( mi, mi ), apb( mi, mi ), array( mi, * ),
     &          cc( mxcmu, mxcmu ), cmu( mxcmu ), cwt( mxcmu ),
     &          eval( mi ), evecc( mxcmu, mxcmu ), gc( mxcmu, mxcmu ),
     &          gl( 0:mxcmu ), kk( mxcmu ), ylmc( 0:mxcmu, mxcmu )
      double precision aad( mi, mi ), evald( mi ), eveccd( mi, mi ),
     &                 wkd( mxcmu )
c     ..
c     .. local scalars ..

      integer   ier, iq, jq, kq, l
      real      alpha, beta, gpmigm, gpplgm, sum
c     ..
c     .. external subroutines ..

      external  asymtx, errmsg
c     ..
c     .. intrinsic functions ..

      intrinsic abs, sqrt
c     ..

c                             ** calculate quantities in eqs. ss(5-6),
c                             ** stwl(8b,15,23f)
      do 40 iq = 1, nn

         do 20 jq = 1, nstr

            sum  = 0.0
            do 10 l = mazim, nstr - 1
               sum  = sum + gl( l )*ylmc( l, iq )*ylmc( l, jq )
   10       continue

            cc( iq, jq ) = 0.5*sum*cwt( jq )

   20    continue

         do 30 jq = 1, nn
c                             ** fill remainder of array using symmetry
c                             ** relations  c(-mui,muj) = c(mui,-muj)
c                             ** and        c(-mui,-muj) = c(mui,muj)

            cc( iq + nn, jq ) = cc( iq, jq + nn )
            cc( iq + nn, jq + nn ) = cc( iq, jq )

c                                       ** get factors of coeff. matrix
c                                       ** of reduced eigenvalue problem

            alpha  = cc( iq, jq ) / cmu( iq )
            beta   = cc( iq, jq + nn ) / cmu( iq )
            amb( iq, jq ) = alpha - beta
            apb( iq, jq ) = alpha + beta

   30    continue

         amb( iq, iq ) = amb( iq, iq ) - 1.0 / cmu( iq )
         apb( iq, iq ) = apb( iq, iq ) - 1.0 / cmu( iq )

   40 continue
c                      ** finish calculation of coefficient matrix of
c                      ** reduced eigenvalue problem:  get matrix
c                      ** product (alfa+beta)*(alfa-beta); ss(12),
c                      ** stwl(23f)
      do 70 iq = 1, nn

         do 60 jq = 1, nn

            sum  = 0.
            do 50 kq = 1, nn
               sum  = sum + apb( iq, kq )*amb( kq, jq )
   50       continue

            array( iq, jq ) = sum

   60    continue

   70 continue
c                      ** find (real) eigenvalues and eigenvectors

      call asymtx( array, evecc, eval, nn, mi, mxcmu, ier, wkd, aad,
     &             eveccd, evald )

      tf = .true.
      if( ier.gt.0 ) then

         write( *, '(//,a,i4,a)' ) ' asymtx--eigenvalue no. ',
     &      ier, '  didnt converge.  lower-numbered eigenvalues wrong.'

         call errmsg( 'asymtx--convergence problems',tf)

      end if


      do 80 iq = 1, nn
         eval( iq )    = sqrt( abs( eval( iq ) ) )
         kk( iq + nn ) = eval( iq )
c                                      ** add negative eigenvalue
         kk( nn + 1 - iq ) = -eval( iq )
   80 continue

c                          ** find eigenvectors (g+) + (g-) from ss(10)
c                          ** and store temporarily in apb array
      do 110 jq = 1, nn

         do 100 iq = 1, nn

            sum  = 0.
            do 90 kq = 1, nn
               sum  = sum + amb( iq, kq )*evecc( kq, jq )
   90       continue

            apb( iq, jq ) = sum / eval( jq )

  100    continue

  110 continue


      do 130 jq = 1, nn

         do 120 iq = 1, nn

            gpplgm = apb( iq, jq )
            gpmigm = evecc( iq, jq )
c                                ** recover eigenvectors g+,g- from
c                                ** their sum and difference; stack them
c                                ** to get eigenvectors of full system
c                                ** ss(7) (jq = eigenvector number)

            evecc( iq,      jq ) = 0.5*( gpplgm + gpmigm )
            evecc( iq + nn, jq ) = 0.5*( gpplgm - gpmigm )

c                                ** eigenvectors corresponding to
c                                ** negative eigenvalues (corresp. to
c                                ** reversing sign of 'k' in ss(10) )
            gpplgm = - gpplgm
            evecc(iq,   jq+nn) = 0.5 * ( gpplgm + gpmigm )
            evecc(iq+nn,jq+nn) = 0.5 * ( gpplgm - gpmigm )
            gc( iq+nn,   jq+nn )   = evecc( iq,    jq )
            gc( nn+1-iq, jq+nn )   = evecc( iq+nn, jq )
            gc( iq+nn,   nn+1-jq ) = evecc( iq,    jq+nn )
            gc( nn+1-iq, nn+1-jq ) = evecc( iq+nn, jq+nn )

  120    continue

  130 continue


      return
      end

      subroutine solve0( b, bdr, bem, bplank, cband, cmu, cwt, expbea,
     &                   fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim,
     &                   mi, mi9m2, mxcmu, ncol, ncut, nn, nstr, nnlyri,
     &                   pi, tplank, taucpr, umu0, z, zz, zplk0, zplk1 )

c        construct right-hand side vector b for general boundary
c        conditions stwj(17) and solve system of equations obtained
c        from the boundary conditions and the continuity-of-
c        intensity-at-layer-interface equations.
c        thermal emission contributes only in azimuthal independence.
c
c
c    i n p u t      v a r i a b l e s:
c
c       bdr      :  surface bidirectional reflectivity
c
c       bem      :  surface bidirectional emissivity
c
c       bplank   :  bottom boundary thermal emission
c
c       cband    :  left-hand side matrix of linear system eq. sc(5),
c                   scaled by eq. sc(12); in banded form required
c                   by linpack solution routines
c
c       cmu,cwt  :  abscissae, weights for gauss quadrature
c                   over angle cosine
c
c       expbea   :  transmission of incident beam, exp(-taucpr/umu0)
c
c       lyrcut   :  logical flag for truncation of computational layers
c
c       mazim    :  order of azimuthal component
c
c       ncol     :  number of columns in cband
c
c       nn       :  order of double-gauss quadrature (nstr/2)
c
c       ncut     :  total number of computational layers considered
c
c       tplank   :  top boundary thermal emission
c
c       taucpr   :  cumulative optical depth (delta-m-scaled)
c
c       zz       :  beam source vectors in eq. ss(19), stwl(24b)
c
c       zplk0    :  thermal source vectors z0, by solving eq. ss(16),
c                   y0 in stwl(26b)
c
c       zplk1    :  thermal source vectors z1, by solving eq. ss(16),
c                   y1 in stwl(26a)
c
c       (remainder are disort input variables)
c
c
c    o u t p u t     v a r i a b l e s:
c
c       b        :  right-hand side vector of eq. sc(5) going into
c                   sgbsl; returns as solution vector of eq. sc(12),
c                   constants of integration without exponential term
c
c      ll        :  permanent storage for b, but re-ordered
c
c
c   i n t e r n a l    v a r i a b l e s:
c
c       ipvt     :  integer vector of pivot indices
c       it       :  pointer for position in  b
c       ncd      :  number of diagonals below or above main diagonal
c       rcond    :  indicator of singularity for cband
c       z        :  scratch array required by sgbco
c
c   called by- disort
c   calls- zeroit, sgbco, errmsg, sgbsl
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   lamber, lyrcut
      logical   tf
      integer   mazim, mi, mi9m2, mxcmu, ncol, ncut, nn, nnlyri, nstr
      real      bplank, fbeam, fisot, pi, tplank, umu0
c     ..
c     .. array arguments ..

      integer   ipvt( * )
      real      b( nnlyri ), bdr( mi, 0:mi ), bem( mi ),
     &          cband( mi9m2, nnlyri ), cmu( mxcmu ), cwt( mxcmu ),
     &          expbea( 0:* ), ll( mxcmu, * ), taucpr( 0:* ),
     &          z( nnlyri ), zplk0( mxcmu, * ), zplk1( mxcmu, * ),
     &          zz( mxcmu, * )
c     ..
c     .. local scalars ..

      integer   ipnt, iq, it, jq, lc, ncd
      real      rcond, sum
c     ..
c     .. external subroutines ..

      external  errmsg, sgbco, sgbsl, zeroit
c     ..

      call zeroit( b, nnlyri )
c                              ** construct b,  stwj(20a,c) for
c                              ** parallel beam + bottom reflection +
c                              ** thermal emission at top and/or bottom

      if( mazim.gt.0 .and. fbeam.gt.0.0 ) then

c                                         ** azimuth-dependent case
c                                         ** (never called if fbeam = 0)
         if( lyrcut .or. lamber ) then

c               ** no azimuthal-dependent intensity for lambert surface;
c               ** no intensity component for truncated bottom layer

            do 10 iq = 1, nn
c                                                  ** top boundary
               b( iq ) = -zz( nn + 1 - iq, 1 )
c                                                  ** bottom boundary

               b( ncol - nn + iq ) = -zz( iq + nn, ncut )*expbea( ncut )

   10       continue


         else

            do 30 iq = 1, nn

               b( iq ) = -zz( nn + 1 - iq, 1 )

               sum  = 0.
               do 20 jq = 1, nn
                  sum  = sum + cwt( jq )*cmu( jq )*bdr( iq, jq )*
     &                         zz( nn + 1 - jq, ncut )*expbea( ncut )
   20          continue

               b( ncol - nn + iq ) = sum
               if( fbeam.gt.0.0 ) b( ncol - nn + iq ) = sum +
     &             ( bdr( iq,0 )*umu0*fbeam/pi - zz( iq+nn,ncut ) )*
     &             expbea( ncut )

   30       continue

         end if
c                             ** continuity condition for layer
c                             ** interfaces of eq. stwj(20b)
         it  = nn

         do 50 lc = 1, ncut - 1

            do 40 iq = 1, nstr
               it  = it + 1
               b( it ) = ( zz( iq,lc+1 ) - zz( iq,lc ) )*expbea( lc )
   40       continue

   50    continue


      else
c                                   ** azimuth-independent case

         if( fbeam.eq.0.0 ) then

            do 60 iq = 1, nn
c                                      ** top boundary

               b( iq ) = -zplk0( nn + 1 - iq, 1 ) + fisot + tplank

   60       continue


            if( lyrcut ) then
c                               ** no intensity component for truncated
c                               ** bottom layer
               do 70 iq = 1, nn
c                                      ** bottom boundary

                  b( ncol - nn + iq ) = - zplk0( iq + nn, ncut ) -
     &                                    zplk1( iq + nn, ncut ) *
     &                                    taucpr( ncut )
   70          continue


            else

               do 90 iq = 1, nn

                  sum  = 0.
                  do 80 jq = 1, nn
                     sum  = sum + cwt( jq )*cmu( jq )*bdr( iq, jq )*
     &                        ( zplk0( nn+1-jq, ncut ) +
     &                          zplk1( nn+1-jq, ncut ) *taucpr( ncut ) )
   80             continue

                  b( ncol - nn + iq ) = 2.*sum + bem( iq )*bplank -
     &                                  zplk0( iq + nn, ncut ) -
     &                                  zplk1( iq + nn, ncut ) *
     &                                  taucpr( ncut )
   90          continue

            end if
c                             ** continuity condition for layer
c                             ** interfaces, stwj(20b)
            it  = nn
            do 110 lc = 1, ncut - 1

               do 100 iq = 1, nstr
                  it  = it + 1
                  b( it ) =   zplk0( iq, lc + 1 ) - zplk0( iq, lc ) +
     &                      ( zplk1( iq, lc + 1 ) - zplk1( iq, lc ) )*
     &                      taucpr( lc )
  100          continue

  110       continue


         else

            do 120 iq = 1, nn
               b( iq ) = -zz( nn + 1 - iq, 1 ) -
     &                   zplk0( nn + 1 - iq, 1 ) + fisot + tplank
  120       continue

            if( lyrcut ) then

               do 130 iq = 1, nn
                  b( ncol-nn+iq ) = - zz(iq+nn, ncut) * expbea(ncut)
     &                              - zplk0(iq+nn, ncut)
     &                              - zplk1(iq+nn, ncut) * taucpr(ncut)
  130          continue


            else

               do 150 iq = 1, nn

                  sum  = 0.
                  do 140 jq = 1, nn
                     sum = sum + cwt(jq) * cmu(jq) * bdr(iq,jq)
     &                          * ( zz(nn+1-jq, ncut) * expbea(ncut)
     &                            + zplk0(nn+1-jq, ncut)
     &                            + zplk1(nn+1-jq, ncut) * taucpr(ncut))
  140             continue

                  b(ncol-nn+iq) = 2.*sum + ( bdr(iq,0) * umu0*fbeam/pi
     &                                - zz(iq+nn, ncut) ) * expbea(ncut)
     &                            + bem(iq) * bplank
     &                            - zplk0(iq+nn, ncut)
     &                            - zplk1(iq+nn, ncut) * taucpr(ncut)
  150          continue

            end if


            it  = nn

            do 170 lc = 1, ncut - 1

               do 160 iq = 1, nstr

                  it  = it + 1
                  b(it) = ( zz(iq,lc+1) - zz(iq,lc) ) * expbea(lc)
     &                    + zplk0(iq,lc+1) - zplk0(iq,lc) +
     &                    ( zplk1(iq,lc+1) - zplk1(iq,lc) ) * taucpr(lc)
  160          continue

  170       continue

         end if

      end if
c                     ** find l-u (lower/upper triangular) decomposition
c                     ** of band matrix cband and test if it is nearly
c                     ** singular (note: cband is destroyed)
c                     ** (cband is in linpack packed format)
      rcond  = 0.0
      ncd    = 3*nn - 1

      call sgbco( cband, mi9m2, ncol, ncd, ncd, ipvt, rcond, z )

      tf = .false.
      if( 1.0 + rcond.eq.1.0 )
     &    call errmsg('solve0--sgbco says matrix near singular',tf)

c                   ** solve linear system with coeff matrix cband
c                   ** and r.h. side(s) b after cband has been l-u
c                   ** decomposed.  solution is returned in b.

      call sgbsl( cband, mi9m2, ncol, ncd, ncd, ipvt, b, 0 )

c                   ** zero cband (it may contain 'foreign'
c                   ** elements upon returning from linpack);
c                   ** necessary to prevent errors

      call zeroit( cband, mi9m2*nnlyri )

      do 190 lc = 1, ncut

         ipnt  = lc*nstr - nn

         do 180 iq = 1, nn
            ll( nn + 1 - iq, lc ) = b( ipnt + 1 - iq )
            ll( iq + nn,     lc ) = b( iq + ipnt )
  180    continue

  190 continue


      return
      end

      subroutine surfac( albedo, delm0, cmu, fbeam, lamber, iref, mi,
     &                   mazim, mxumu, nn, numu, onlyfl, pi, umu, umu0,
     &                   usrang, surf_pr, bdr, emu, bem, rmu )

c       computes user's surface bidirectional properties, stwl(41)
c
c   i n p u t     v a r i a b l e s:
c
c       cmu    :  computational polar angle cosines (gaussian)
c
c       delm0  :  kronecker delta, delta-sub-m0
c
c       mazim  :  order of azimuthal component
c
c       nn     :  order of double-gauss quadrature (nstr/2)
c
c       iref   : bidirectional reflectance options
c             0 - Lambert
c             1 - Hapke's BDR model
c             2 - Breon's BDR model; combination of Li + Roujean
c             3 - Roujean's BDR model
c             4 - Cox and Munk glint model
c
c   SURF_PR : Wavelength dependent surface properties array 
c             IREF= 0 - Lambert albedo
c             IREF= 1 - Hapke : HH, W
c             IREF= 2 - Breon's BDR model: k0, k1, k2
c             IREF= 3 - Roujean's BDR model: k0, k1, k2
c             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw
c
c       (remainder are 'disort' input variables)
c
c    o u t p u t     v a r i a b l e s:
c
c       bdr :  fourier expansion coefficient of surface bidirectional
c                 reflectivity (computational angles)
c
c       rmu :  surface bidirectional reflectivity (user angles)
c
c       bem :  surface directional emissivity (computational angles)
c
c       emu :  surface directional emissivity (user angles)
c
c    i n t e r n a l     v a r i a b l e s:

c       dref   :  directional reflectivity
c
c       nmug   :  number of angle cosine quadrature points on (-1,1)
c                 for integrating bidirectional reflectivity to get
c                 directional emissivity (it is necessary to use a
c                 quadrature set distinct from the computational angles,
c                 because the computational angles may not be dense
c                 enough -- i.e. 'nstr' may be too small-- to give an
c                 accurate approximation for the integration).
c
c       gmu    :  the 'nmug' angle cosine quadrature points on (0,1)
c
c       gwt    :  the 'nmug' angle cosine quadrature weights on (0,1)
c
c   called by- disort
c   calls- qgausn, bdref, zeroit
c+---------------------------------------------------------------------+

c     .. parameters ..

      integer   nmug
      parameter ( nmug = 50 )
c     ..
c     .. scalar arguments ..

      logical   lamber, onlyfl, usrang
      integer   mazim, mi, mxumu, nn, numu
      integer   iref
      real      albedo, delm0, fbeam, pi, umu0
c     ..
c     .. array arguments ..
c
      real      bdr( mi, 0:mi ), bem( mi ), cmu( * ), emu( mxumu ),
     &          rmu( mxumu, 0:mi ), umu( * )
      real      surf_pr( * )
c     ..
c     .. local scalars ..

      logical   pass1
      integer   iq, iu, jg, jq, k
      real      dref, sum, pigmu
c     ..
c     .. local arrays ..

      real      gmu( nmug ), gwt( nmug )
c     ..
c     .. external functions ..

      real      bdref
      external  bdref
c     ..
c     .. external subroutines ..

      external  qgausn, zeroit
c     ..
c     .. intrinsic functions ..

      intrinsic cos
c     ..
      save      pass1, gmu, gwt
      data      pass1 / .true. /


      if( pass1 ) then

         pass1  = .false.

         call qgausn( nmug/2, gmu, gwt )

         do 10 k = 1, nmug / 2
            gmu( k + nmug/2 ) = -gmu( k )
            gwt( k + nmug/2 ) = gwt( k )
   10    continue

      end if


      call zeroit( bdr, mi*( mi+1 ) )
      call zeroit( bem, mi )

c                             ** compute fourier expansion coefficient
c                             ** of surface bidirectional reflectance
c                             ** at computational angles eq. stwl (41)

      if( lamber .and. mazim.eq.0 ) then

         do 30 iq = 1, nn

            bem( iq ) = 1.0 - albedo

            do 20 jq = 0, nn
               bdr( iq, jq ) = albedo
   20       continue

   30    continue

      else if( .not.lamber ) then

         do 70 iq = 1, nn

            do 50 jq = 1, nn

               sum  = 0.0
               do 40 k = 1, nmug
                  pigmu = pi*gmu( k )
                  sum  = sum + gwt( k ) *
     &                   bdref( cmu(iq), cmu(jq),
     &                          pigmu, surf_pr, iref ) * 
     &                          cos( mazim*pigmu )
   40          continue

               bdr( iq, jq ) = 0.5 * ( 2. - delm0 ) * sum

   50       continue


            if( fbeam.gt.0.0 ) then

               sum  = 0.0
               do 60 k = 1, nmug
                  pigmu = pi*gmu( k )
                  sum  = sum + gwt( k ) *
     &                   bdref( cmu(iq), umu0,
     &                          pigmu, surf_pr, iref ) * 
     &                          cos( mazim*pigmu )
   60          continue

               bdr( iq, 0 ) = 0.5 * ( 2. - delm0 ) * sum

            end if

   70    continue


         if( mazim.eq.0 ) then

c                             ** integrate bidirectional reflectivity
c                             ** at reflection polar angle cosines -cmu-
c                             ** and incident angle cosines -gmu- to get
c                             ** directional emissivity at computational
c                             ** angle cosines -cmu-.
            do 100 iq = 1, nn

               dref  = 0.0

               do 90 jg = 1, nmug
                  pigmu = pi*gmu( jg )

                  sum  = 0.0
                  do 80 k = 1, nmug / 2
                     sum  = sum + gwt( k ) * gmu( k ) *
     &                      bdref( cmu(iq), gmu(k),
     &                             pigmu, surf_pr, iref )
   80             continue

                  dref  = dref + gwt( jg )*sum

   90          continue

               bem( iq ) = 1.0 - dref

  100       continue

         end if

      end if
c                             ** compute fourier expansion coefficient
c                             ** of surface bidirectional reflectance
c                             ** at user angles eq. stwl (41)

      if( .not.onlyfl .and. usrang ) then

         call zeroit( emu, mxumu )
         call zeroit( rmu, mxumu*( mi+1 ) )

         do 170 iu = 1, numu

            if( umu( iu ).gt.0.0 ) then

               if( lamber .and. mazim.eq.0 ) then

                  do 110 iq = 0, nn
                     rmu( iu, iq ) = albedo
  110             continue

                  emu( iu ) = 1.0 - albedo

               else if( .not.lamber ) then

                  do 130 iq = 1, nn

                     sum  = 0.0
                     do 120 k = 1, nmug
                        pigmu = pi*gmu( k )
                        sum  = sum + gwt( k ) *
     &                         bdref( umu(iu), cmu(iq),
     &                                pigmu, surf_pr, iref ) *
     &                           cos( mazim*pigmu )
  120                continue

                     rmu( iu, iq ) = 0.5 * ( 2. - delm0 ) * sum

  130             continue

                  if( fbeam.gt.0.0 ) then

                     sum  = 0.0
                     do 140 k = 1, nmug
                        pigmu = pi*gmu(k)
                        sum  = sum + gwt( k ) *
     &                         bdref( umu(iu), umu0,
     &                                pigmu,surf_pr, iref ) *
     &                           cos( mazim*pigmu )
  140                continue

                     rmu( iu, 0 ) = 0.5 * ( 2. - delm0 ) * sum

                  end if


                  if( mazim.eq.0 ) then

c                               ** integrate bidirectional reflectivity
c                               ** at reflection angle cosines -umu- and
c                               ** incident angle cosines -gmu- to get
c                               ** directional emissivity at
c                               ** user angle cosines -umu-.
                     dref  = 0.0

                     do 160 jg = 1, nmug
                        pigmu = pi*gmu( jg )

                        sum  = 0.0
                        do 150 k = 1, nmug / 2
                           sum  = sum + gwt( k )*gmu( k )*
     &                            bdref( umu(iu),
     &                                   gmu(k), pigmu, surf_pr, iref )
  150                   continue

                        dref  = dref + gwt( jg ) * sum

  160                continue

                     emu( iu ) = 1.0 - dref

                  end if

               end if

            end if

  170    continue

      end if


      return
      end

      subroutine terpev( cwt, evecc, gl, gu, mazim, mxcmu, mxumu, nn,
     &                   nstr, numu, wk, ylmc, ylmu )

c         interpolate eigenvectors to user angles; eq sd(8)

c   called by- disort, albtrn
c --------------------------------------------------------------------+

c     .. scalar arguments ..

      integer   mazim, mxcmu, mxumu, nn, nstr, numu
c     ..
c     .. array arguments ..

      real      cwt( mxcmu ), evecc( mxcmu, mxcmu ), gl( 0:mxcmu ),
     &          gu( mxumu, mxcmu ), wk( mxcmu ), ylmc( 0:mxcmu, mxcmu ),
     &          ylmu( 0:mxcmu, mxumu )
c     ..
c     .. local scalars ..

      integer   iq, iu, jq, l
      real      sum
c     ..


      do 50 iq = 1, nstr

         do 20 l = mazim, nstr - 1
c                                   ** inner sum in sd(8) times all
c                                   ** factors in outer sum but plm(mu)
            sum  = 0.0
            do 10 jq = 1, nstr
               sum  = sum + cwt( jq )*ylmc( l, jq )*evecc( jq, iq )
   10       continue

            wk( l + 1 ) = 0.5*gl( l )*sum

   20    continue
c                                    ** finish outer sum in sd(8)
c                                    ** and store eigenvectors
         do 40 iu = 1, numu

            sum  = 0.
            do 30 l = mazim, nstr - 1
               sum  = sum + wk( l + 1 )*ylmu( l, iu )
   30       continue

            if( iq.le.nn ) gu( iu, iq + nn )       = sum
            if( iq.gt.nn ) gu( iu, nstr + 1 - iq ) = sum

   40    continue

   50 continue


      return
      end

      subroutine terpso( cwt, delm0, fbeam, gl, mazim, mxcmu, plank,
     &                   numu, nstr, oprim, pi, ylm0, ylmc, ylmu, psi0,
     &                   psi1, xr0, xr1, z0, z1, zj, zbeam, z0u, z1u )

c         interpolates source functions to user angles, eq. stwl(30)
c
c
c    i n p u t      v a r i a b l e s:
c
c       cwt    :  weights for gauss quadrature over angle cosine
c
c       delm0  :  kronecker delta, delta-sub-m0
c
c       gl     :  delta-m scaled legendre coefficients of phase function
c                 (including factors 2l+1 and single-scatter albedo)
c
c       mazim  :  order of azimuthal component
c
c       oprim  :  single scattering albedo
c
c       xr0    :  expansion of thermal source function, eq. stwl(24d)
c
c       xr1    :  expansion of thermal source function eq. stwl(24d)
c
c       ylm0   :  normalized associated legendre polynomial
c                 at the beam angle
c
c       ylmc   :  normalized associated legendre polynomial
c                 at the quadrature angles
c
c       ylmu   :  normalized associated legendre polynomial
c                 at the user angles
c
c       z0     :  solution vectors z-sub-zero of eq. ss(16), stwl(26a)
c
c       z1     :  solution vectors z-sub-one  of eq. ss(16), stwl(26b)
c
c       zj     :  solution vector z-sub-zero after solving eq. ss(19),
c                 stwl(24b)
c
c       (remainder are disort input variables)
c
c
c    o u t p u t     v a r i a b l e s:
c
c       zbeam  :  incident-beam source function at user angles
c
c       z0u,z1u:  components of a linear-in-optical-depth-dependent
c                 source (approximating the planck emission source)
c
c
c   i n t e r n a l    v a r i a b l e s:
c
c       psi0  :  sum just after square bracket in  eq. sd(9)
c       psi1  :  sum in eq. stwl(31d)
c
c   called by- disort
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   plank
      integer   mazim, mxcmu, nstr, numu
      real      delm0, fbeam, oprim, pi, xr0, xr1
c     ..
c     .. array arguments ..

      real      cwt( mxcmu ), gl( 0:mxcmu ), psi0( mxcmu ),
     &          psi1( mxcmu ), ylm0( 0:mxcmu ), ylmc( 0:mxcmu, mxcmu ),
     &          ylmu( 0:mxcmu, * ), z0( mxcmu ), z0u( * ), z1( mxcmu ),
     &          z1u( * ), zbeam( * ), zj( mxcmu )
c     ..
c     .. local scalars ..

      integer   iq, iu, jq
      real      fact, psum, psum0, psum1, sum, sum0, sum1
c     ..


      if( fbeam.gt.0.0 ) then
c                                  ** beam source terms; eq. sd(9)

         do 20 iq = mazim, nstr - 1

            psum   = 0.
            do 10 jq = 1, nstr
               psum  = psum + cwt( jq )*ylmc( iq, jq )*zj( jq )
   10       continue

            psi0( iq + 1 ) = 0.5*gl( iq )*psum

   20    continue

         fact   = ( 2. - delm0 )*fbeam / ( 4.0*pi )

         do 40 iu = 1, numu

            sum    = 0.
            do 30 iq = mazim, nstr - 1
               sum  = sum + ylmu( iq, iu )*
     &                    ( psi0( iq+1 ) + fact*gl( iq )*ylm0( iq ) )
   30       continue

            zbeam( iu ) = sum

   40    continue

      end if


      if( plank .and. mazim.eq.0 ) then

c                          ** thermal source terms, stwj(27c), stwl(31c)
c
         do 60 iq = mazim, nstr - 1

            psum0  = 0.0
            psum1  = 0.0
            do 50 jq = 1, nstr
               psum0  = psum0 + cwt( jq )*ylmc( iq, jq )*z0( jq )
               psum1  = psum1 + cwt( jq )*ylmc( iq, jq )*z1( jq )
   50       continue

            psi0( iq + 1 ) = 0.5*gl( iq ) * psum0
            psi1( iq + 1 ) = 0.5*gl( iq ) * psum1

   60    continue

         do 80 iu = 1, numu

            sum0   = 0.0
            sum1   = 0.0
            do 70 iq = mazim, nstr - 1
               sum0  = sum0 + ylmu( iq, iu ) * psi0( iq + 1 )
               sum1  = sum1 + ylmu( iq, iu ) * psi1( iq + 1 )
   70       continue

            z0u( iu ) = sum0 + ( 1. - oprim ) * xr0
            z1u( iu ) = sum1 + ( 1. - oprim ) * xr1

   80    continue

      end if


      return
      end

      subroutine upbeam( array, cc, cmu, delm0, fbeam, gl, ipvt, mazim,
     &                   mxcmu, nn, nstr, pi, umu0, wk, ylm0, ylmc, zj,
     &                   zz )

c         finds the incident-beam particular solution of ss(18),
c         stwl(24a)
c
c   i n p u t    v a r i a b l e s:
c
c       cc     :  c-sub-ij in eq. ss(5)
c
c       cmu    :  abscissae for gauss quadrature over angle cosine
c
c       delm0  :  kronecker delta, delta-sub-m0
c
c       gl     :  delta-m scaled legendre coefficients of phase function
c                 (including factors 2l+1 and single-scatter albedo)
c
c       mazim  :  order of azimuthal component
c
c       ylm0   :  normalized associated legendre polynomial
c                 at the beam angle
c
c       ylmc   :  normalized associated legendre polynomial
c                 at the quadrature angles
c
c       (remainder are disort input variables)
c
c
c   o u t p u t    v a r i a b l e s:
c
c       zj     :  right-hand side vector x-sub-zero in ss(19),stwl(24b);
c                 also the solution vector z-sub-zero after solving
c                 that system
c
c       zz     :  permanent storage for zj, but re-ordered
c
c
c   i n t e r n a l    v a r i a b l e s:
c
c       array  :  coefficient matrix in left-hand side of eq. ss(19),
c                   stwl(24b)
c       ipvt   :  integer vector of pivot indices required by linpack
c       wk     :  scratch array required by linpack
c
c   called by- disort
c   calls- sgeco, errmsg, sgesl
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   tf
      integer   mazim, mxcmu, nn, nstr
      real      delm0, fbeam, pi, umu0
c     ..
c     .. array arguments ..

      integer   ipvt( * )
      real      array( mxcmu, mxcmu ), cc( mxcmu, mxcmu ), cmu( mxcmu ),
     &          gl( 0:mxcmu ), wk( mxcmu ), ylm0( 0:mxcmu ),
     &          ylmc( 0:mxcmu, * ), zj( mxcmu ), zz( mxcmu )
c     ..
c     .. local scalars ..

      integer   iq, job, jq, k
      real      rcond, sum
c     ..
c     .. external subroutines ..

      external  errmsg, sgeco, sgesl
c     ..


      do 30 iq = 1, nstr

         do 10 jq = 1, nstr
            array( iq, jq ) = -cc( iq, jq )
   10    continue

         array( iq, iq ) = 1.+ cmu( iq ) / umu0 + array( iq, iq )

         sum  = 0.
         do 20 k = mazim, nstr - 1
            sum  = sum + gl( k )*ylmc( k, iq )*ylm0( k )
   20    continue

         zj( iq ) = ( 2.- delm0 )*fbeam*sum / ( 4.*pi )
   30 continue

c                  ** find l-u (lower/upper triangular) decomposition
c                  ** of array and see if it is nearly singular
c                  ** (note:  array is altered)
      rcond  = 0.0

      call sgeco( array, mxcmu, nstr, ipvt, rcond, wk )

      tf = .false.
      if( 1.0 + rcond.eq.1.0 )
     &    call errmsg('upbeam--sgeco says matrix near singular',tf)

c                ** solve linear system with coeff matrix array
c                ** (assumed already l-u decomposed) and r.h. side(s)
c                ** zj;  return solution(s) in zj
      job  = 0

      call sgesl( array, mxcmu, nstr, ipvt, zj, job )


      do 40 iq = 1, nn
         zz( iq + nn )     = zj( iq )
         zz( nn + 1 - iq ) = zj( iq + nn )
   40 continue


      return
      end

      subroutine upisot( array, cc, cmu, ipvt, mxcmu, nn, nstr, oprim,
     &                   wk, xr0, xr1, z0, z1, zplk0, zplk1 )

c       finds the particular solution of thermal radiation of stwl(25)
c
c
c
c    i n p u t     v a r i a b l e s:
c
c       cc     :  c-sub-ij in eq. ss(5), stwl(8b)
c
c       cmu    :  abscissae for gauss quadrature over angle cosine
c
c       oprim  :  delta-m scaled single scattering albedo
c
c       xr0    :  expansion coefficient b-sub-zero of thermal source
c                   function, eq. stwl(24c)
c
c       xr1    :  expansion coefficient b-sub-one of thermal source
c                   function eq. stwl(24c)
c
c       (remainder are disort input variables)
c
c
c    o u t p u t    v a r i a b l e s:
c
c       z0     :  solution vectors z-sub-zero of eq. ss(16), stwl(26a)
c
c       z1     :  solution vectors z-sub-one  of eq. ss(16), stwl(26b)
c
c       zplk0, :  permanent storage for z0,z1, but re-ordered
c        zplk1
c
c
c   i n t e r n a l    v a r i a b l e s:
c
c       array  :  coefficient matrix in left-hand side of eq. ss(16)
c       ipvt   :  integer vector of pivot indices required by linpack
c       wk     :  scratch array required by linpack
c
c   called by- disort
c   calls- sgeco, errmsg, sgesl
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   tf
      integer   mxcmu, nn, nstr
      real      oprim, xr0, xr1
c     ..
c     .. array arguments ..

      integer   ipvt( * )
      real      array( mxcmu, mxcmu ), cc( mxcmu, mxcmu ), cmu( mxcmu ),
     &          wk( mxcmu ), z0( mxcmu ), z1( mxcmu ), zplk0( mxcmu ),
     &          zplk1( mxcmu )
c     ..
c     .. local scalars ..

      integer   iq, jq
      real      rcond
c     ..
c     .. external subroutines ..

      external  errmsg, sgeco, sgesl
c     ..


      do 20 iq = 1, nstr

         do 10 jq = 1, nstr
            array( iq, jq ) = -cc( iq, jq )
   10    continue

         array( iq, iq ) = 1.0 + array( iq, iq )

         z1( iq ) = ( 1. - oprim ) * xr1

   20 continue
c                       ** solve linear equations: same as in upbeam,
c                       ** except zj replaced by z1 and z0
      rcond  = 0.0

      call sgeco( array, mxcmu, nstr, ipvt, rcond, wk )

      tf = .false.
      if( 1.0 + rcond.eq.1.0 )
     &    call errmsg('upisot--sgeco says matrix near singular',tf)

      call sgesl( array, mxcmu, nstr, ipvt, z1, 0 )

      do 30 iq = 1, nstr
         z0( iq ) = ( 1. - oprim ) * xr0 + cmu( iq ) * z1( iq )
   30 continue

      call sgesl( array, mxcmu, nstr, ipvt, z0, 0 )

      do 40 iq = 1, nn
         zplk0( iq + nn ) = z0( iq )
         zplk1( iq + nn ) = z1( iq )
         zplk0( nn + 1 - iq ) = z0( iq + nn )
         zplk1( nn + 1 - iq ) = z1( iq + nn )
   40 continue


      return
      end

      subroutine usrint( bplank, cmu, cwt, delm0, dtaucp, emu, expbea,
     &                   fbeam, fisot, gc, gu, kk, lamber, layru, ll,
     &                   lyrcut, mazim, mxcmu, mxulv, mxumu, ncut, nlyr,
     &                   nn, nstr, plank, numu, ntau, pi, rmu, taucpr,
     &                   tplank, umu, umu0, utaupr, wk, zbeam, z0u, z1u,
     &                   zz, zplk0, zplk1, uum )

c       computes intensity components at user output angles
c       for azimuthal expansion terms in eq. sd(2), stwl(6)
c
c
c   i n p u t    v a r i a b l e s:
c
c       bplank :  integrated planck function for emission from
c                 bottom boundary
c
c       cmu    :  abscissae for gauss quadrature over angle cosine
c
c       cwt    :  weights for gauss quadrature over angle cosine
c
c       delm0  :  kronecker delta, delta-sub-m0
c
c       emu    :  surface directional emissivity (user angles)
c
c       expbea :  transmission of incident beam, exp(-taucpr/umu0)
c
c       gc     :  eigenvectors at polar quadrature angles, sc(1)
c
c       gu     :  eigenvectors interpolated to user polar angles
c                    (i.e., g in eq. sc(1) )
c
c       kk     :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b)
c
c       layru  :  layer number of user level utau
c
c       ll     :  constants of integration in eq. sc(1), obtained
c                 by solving scaled version of eq. sc(5);
c                 exponential term of eq. sc(12) not included
c
c       lyrcut :  logical flag for truncation of computational layer
c
c       mazim  :  order of azimuthal component
c
c       ncut   :  total number of computational layers considered
c
c       nn     :  order of double-gauss quadrature (nstr/2)
c
c       rmu    :  surface bidirectional reflectivity (user angles)
c
c       taucpr :  cumulative optical depth (delta-m-scaled)
c
c       tplank :  integrated planck function for emission from
c                 top boundary
c
c       utaupr :  optical depths of user output levels in delta-m
c                 coordinates;  equal to utau if no delta-m
c
c       z0u    :  z-sub-zero in eq. ss(16) interpolated to user
c                 angles from an equation derived from ss(16),
c                 y-sub-zero on stwl(26b)
c
c       z1u    :  z-sub-one in eq. ss(16) interpolated to user
c                 angles from an equation derived from ss(16),
c                 y-sub-one in stwl(26a)
c
c       zz     :  beam source vectors in eq. ss(19), stwl(24b)
c
c       zplk0  :  thermal source vectors z0, by solving eq. ss(16),
c                 y-sub-zero in stwl(26)
c
c       zplk1  :  thermal source vectors z1, by solving eq. ss(16),
c                 y-sub-one in stwl(26)
c
c       zbeam  :  incident-beam source vectors
c
c       (remainder are disort input variables)
c
c
c    o u t p u t    v a r i a b l e s:
c
c       uum    :  azimuthal components of the intensity in eq. stwj(5),
c                 stwl(6)
c
c
c    i n t e r n a l    v a r i a b l e s:
c
c       bnddir :  direct intensity down at the bottom boundary
c       bnddfu :  diffuse intensity down at the bottom boundary
c       bndint :  intensity attenuated at both boundaries, stwj(25-6)
c       dtau   :  optical depth of a computational layer
c       lyrend :  end layer of integration
c       lyrstr :  start layer of integration
c       palint :  intensity component from parallel beam
c       plkint :  intensity component from planck source
c       wk     :  scratch vector for saving exp evaluations
c
c       all the exponential factors ( exp1, expn,... etc.)
c       come from the substitution of constants of integration in
c       eq. sc(12) into eqs. s1(8-9).  they all have negative
c       arguments so there should never be overflow problems.
c
c   called by- disort
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   lamber, lyrcut, plank
      integer   mazim, mxcmu, mxulv, mxumu, ncut, nlyr, nn, nstr, ntau,
     &          numu
      real      bplank, delm0, fbeam, fisot, pi, tplank, umu0
c     ..
c     .. array arguments ..

      integer   layru( * )
      real      cmu( mxcmu ), cwt( mxcmu ), dtaucp( * ), emu( mxumu ),
     &          expbea( 0:* ), gc( mxcmu, mxcmu, * ),
     &          gu( mxumu, mxcmu, * ), kk( mxcmu, * ), ll( mxcmu, * ),
     &          rmu( mxumu, 0:* ), taucpr( 0:* ), umu( * ),
     &          utaupr( mxulv ), uum( mxumu, mxulv ), wk( mxcmu ),
     &          z0u( mxumu, * ), z1u( mxumu, * ), zbeam( mxumu, * ),
     &          zplk0( mxcmu, * ), zplk1( mxcmu, * ), zz( mxcmu, * )
c     ..
c     .. local scalars ..

      logical   negumu
      integer   iq, iu, jq, lc, lu, lyrend, lyrstr, lyu
      real      bnddfu, bnddir, bndint, denom, dfuint, dtau, dtau1,
     &          dtau2, exp0, exp1, exp2, expn, f0n, f1n, fact, palint,
     &          plkint, sgn
c     ..
c     .. intrinsic functions ..

      intrinsic abs, exp
c     ..

c                          ** incorporate constants of integration into
c                          ** interpolated eigenvectors
      do 30 lc = 1, ncut

         do 20 iq = 1, nstr

            do 10 iu = 1, numu
               gu( iu, iq, lc ) = gu( iu, iq, lc ) * ll( iq, lc )
   10       continue

   20    continue

   30 continue
c                           ** loop over levels at which intensities
c                           ** are desired ('user output levels')
      do 160 lu = 1, ntau

         if( fbeam.gt.0.0 ) exp0  = exp( -utaupr( lu ) / umu0 )
         lyu  = layru( lu )
c                              ** loop over polar angles at which
c                              ** intensities are desired
         do 150 iu = 1, numu

            if( lyrcut .and. lyu.gt.ncut ) go to  150

            negumu = umu( iu ) .lt. 0.0

            if( negumu ) then

               lyrstr = 1
               lyrend = lyu - 1
               sgn    = -1.0

            else

               lyrstr = lyu + 1
               lyrend = ncut
               sgn    = 1.0

            end if
c                          ** for downward intensity, integrate from top
c                          ** to lyu-1 in eq. s1(8); for upward,
c                          ** integrate from bottom to lyu+1 in s1(9)
            palint = 0.0
            plkint = 0.0

            do 60 lc = lyrstr, lyrend

               dtau = dtaucp( lc )
               exp1 = exp( ( utaupr(lu) - taucpr(lc-1) ) / umu( iu ) )
               exp2 = exp( ( utaupr(lu) - taucpr(lc)   ) / umu( iu ) )

               if( plank .and. mazim.eq.0 ) then

c                          ** eqs. stwl(36b,c, 37b,c)
c
                  f0n = sgn * ( exp1 - exp2 )

                  f1n = sgn * ( ( taucpr( lc-1 ) + umu( iu ) ) * exp1 -
     &                          ( taucpr( lc )   + umu( iu ) ) * exp2 )

                  plkint = plkint + z0u( iu,lc )*f0n + z1u( iu,lc )*f1n

               end if


               if( fbeam.gt.0.0 ) then

                  denom  = 1. + umu( iu ) / umu0

                  if( abs( denom ).lt.0.0001 ) then
c                                                   ** l'hospital limit
                     expn   = ( dtau / umu0 )*exp0

                  else

                     expn   = ( exp1*expbea( lc-1 ) -
     &                          exp2*expbea( lc ) ) * sgn / denom

                  end if

                  palint = palint + zbeam( iu, lc )*expn

               end if

c                                                   ** kk is negative
               do 40 iq = 1, nn

                  wk( iq ) = exp( kk( iq,lc )*dtau )
                  denom  = 1.0 + umu( iu )*kk( iq, lc )

                  if( abs( denom ).lt.0.0001 ) then
c                                                   ** l'hospital limit
                     expn   = dtau / umu( iu )*exp2

                  else

                     expn   = sgn*( exp1*wk( iq ) - exp2 ) / denom

                  end if

                  palint = palint + gu( iu, iq, lc )*expn

   40          continue

c                                                   ** kk is positive
               do 50 iq = nn + 1, nstr

                  denom  = 1.0 + umu( iu )*kk( iq, lc )

                  if( abs( denom ).lt.0.0001 ) then
c                                                   ** l'hospital limit
                     expn  = -dtau / umu( iu )*exp1

                  else

                     expn  = sgn*( exp1 - exp2*wk( nstr+1-iq ) ) / denom

                  end if

                  palint = palint + gu( iu, iq, lc )*expn

   50          continue


   60       continue
c                           ** calculate contribution from user
c                           ** output level to next computational level

            dtau1  = utaupr( lu ) - taucpr( lyu - 1 )
            dtau2  = utaupr( lu ) - taucpr( lyu )

            if( abs( dtau1 ).lt.1.e-6 .and. negumu ) go to  90
            if( abs( dtau2 ).lt.1.e-6 .and. (.not.negumu ) ) go to  90

            if( negumu )      exp1  = exp( dtau1/umu( iu ) )
            if( .not.negumu ) exp2  = exp( dtau2/umu( iu ) )

            if( fbeam.gt.0.0 ) then

               denom  = 1. + umu( iu ) / umu0

               if( abs( denom ).lt.0.0001 ) then

                  expn   = ( dtau1 / umu0 )*exp0

               else if( negumu ) then

                  expn  = ( exp0 - expbea( lyu-1 )*exp1 ) / denom

               else

                  expn  = ( exp0 - expbea( lyu )*exp2 ) / denom

               end if

               palint = palint + zbeam( iu, lyu )*expn

            end if

c                                                   ** kk is negative
            dtau  = dtaucp( lyu )

            do 70 iq = 1, nn

               denom  = 1. + umu( iu )*kk( iq, lyu )

               if( abs( denom ).lt.0.0001 ) then

                  expn = -dtau2 / umu( iu )*exp2

               else if( negumu ) then

                  expn = ( exp( -kk( iq,lyu ) * dtau2 ) -
     &                     exp(  kk( iq,lyu ) * dtau  ) * exp1 ) / denom

               else

                  expn = ( exp( -kk( iq,lyu ) * dtau2 ) - exp2 ) / denom

               end if

               palint = palint + gu( iu, iq, lyu )*expn

   70       continue

c                                                   ** kk is positive
            do 80 iq = nn + 1, nstr

               denom  = 1. + umu( iu )*kk( iq, lyu )

               if( abs( denom ).lt.0.0001 ) then

                  expn   = -dtau1 / umu( iu )*exp1

               else if( negumu ) then

                  expn = ( exp( -kk( iq,lyu ) * dtau1 ) - exp1 ) / denom

               else

                  expn = ( exp( -kk( iq,lyu ) * dtau1 ) -
     &                     exp( -kk( iq,lyu ) * dtau  ) * exp2 ) / denom

               end if

               palint = palint + gu( iu, iq, lyu )*expn

   80       continue


            if( plank .and. mazim.eq.0 ) then

c                            ** eqs. stwl (35-37) with tau-sub-n-1
c                            ** replaced by tau for upward, and
c                            ** tau-sub-n replaced by tau for downward
c                            ** directions

               if( negumu ) then

                  expn  = exp1
                  fact  = taucpr( lyu - 1 ) + umu( iu )

               else

                  expn  = exp2
                  fact  = taucpr( lyu ) + umu( iu )

               end if

               f0n  = 1. - expn
               f1n  = utaupr( lu ) + umu( iu ) - fact * expn

               plkint = plkint + z0u( iu, lyu )*f0n + z1u( iu, lyu )*f1n

            end if

c                            ** calculate intensity components
c                            ** attenuated at both boundaries.
c                            ** note: no azimuthal intensity
c                            ** component for isotropic surface
   90       continue
            bndint = 0.0

            if( negumu .and. mazim.eq.0 ) then

               bndint = (fisot + tplank) * exp( utaupr(lu ) / umu(iu) )


            else if( .not.negumu ) then

               if( lyrcut .or. ( lamber.and.mazim.gt.0 ) ) go to  140

               do 100 jq = nn + 1, nstr
                  wk( jq ) = exp( -kk( jq,nlyr )*dtaucp( nlyr ) )
  100          continue

               bnddfu = 0.0

               do 130 iq = nn, 1, -1

                  dfuint = 0.0
                  do 110 jq = 1, nn
                     dfuint = dfuint + gc( iq, jq, nlyr )*ll( jq, nlyr )
  110             continue

                  do 120 jq = nn + 1, nstr
                     dfuint = dfuint + gc( iq, jq, nlyr )*
     &                                 ll( jq, nlyr )*wk( jq )
  120             continue

                  if( fbeam.gt.0.0 ) dfuint = dfuint +
     &                                     zz( iq, nlyr )*expbea( nlyr )

                  dfuint = dfuint + delm0 * ( zplk0( iq, nlyr ) +
     &                              zplk1( iq,nlyr ) *taucpr( nlyr ) )
                  bnddfu = bnddfu + ( 1.+delm0 ) * rmu(iu,nn+1-iq)
     &                            * cmu(nn+1-iq) * cwt(nn+1-iq)* dfuint
  130          continue

               bnddir = 0.0
               if( fbeam.gt.0.0 ) bnddir = umu0*fbeam / pi*rmu( iu, 0 )*
     &                                     expbea( nlyr )

               bndint = ( bnddfu + bnddir + delm0 * emu(iu) * bplank )
     &                  * exp( (utaupr(lu)-taucpr(nlyr)) / umu(iu) )

            end if

  140       continue

            uum( iu, lu ) = palint + plkint + bndint

  150    continue

  160 continue


      return
      end

      real function  xifunc( umu1, umu2, umu3, tau )

c          calculates xi function of eq. stwl (72)

c                    i n p u t   v a r i a b l e s

c        tau         optical thickness of the layer
c
c        umu1,2,3    cosine of zenith angle_1, _2, _3
c
c   called by- secsca
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      real      tau, umu1, umu2, umu3
c     ..
c     .. local scalars ..

      real      exp1, x1, x2
c     ..
c     .. intrinsic functions ..

      intrinsic exp
c     ..


      x1     = 1. / umu1 - 1. / umu2
      x2     = 1. / umu1 - 1. / umu3

      exp1 = exp( -tau/umu1 )

      if( umu2.eq.umu3 .and. umu1.eq.umu2 ) then

         xifunc = tau*tau * exp1 / ( 2.*umu1*umu2 )

      else if( umu2.eq.umu3 .and. umu1.ne.umu2 ) then

         xifunc = ( ( tau - 1./x1 ) * exp( -tau/umu2 ) + exp1 / x1 )
     &            / ( x1*umu1*umu2 )

      else if( umu2.ne.umu3 .and. umu1.eq.umu2 ) then

         xifunc = ( ( exp( -tau/umu3 ) - exp1 ) / x2 - tau * exp1 )
     &            / ( x2*umu1*umu2 )

      else if( umu2.ne.umu3 .and. umu1.eq.umu3 ) then

         xifunc = ( ( exp( -tau/umu2 ) - exp1 ) / x1 - tau * exp1 )
     &            / ( x1*umu1*umu2 )

      else

         xifunc = ( ( exp( -tau/umu3 ) - exp1 ) / x2 -
     &            (   exp( -tau/umu2 ) - exp1 ) / x1 ) /
     &            ( x2*umu1*umu2 )

      end if


      return
      end

c ******************************************************************
c ********** disort service routines ************************
c ******************************************************************

      subroutine chekin( nlyr, dtauc, ssalb, nmom, pmom, temper, wvnm,
     &                   usrtau, ntau, iref, utau, nstr, usrang,
     &                   numu, umu, nphi, phi, ibcnd, fbeam, umu0, phi0,
     &                   fisot, lamber, albedo, btemp, ttemp,temis, 
     &                   surf_pr, plank, onlyfl, deltam, corint, accur,
     &                   tauc, maxcly, maxulv, maxumu, maxphi, maxmom,
     &                   mxcly, mxulv, mxumu, mxcmu, mxphi, mxsqt )

c           checks the input dimensions and variables

c   calls- wrtbad, wrtdim, dref, errmsg
c   called by- disort
c --------------------------------------------------------------------

c     .. scalar arguments ..

      logical   corint, deltam, lamber, onlyfl, plank, usrang, usrtau
      integer   ibcnd, maxcly, maxmom, maxphi, maxulv, maxumu, mxcly,
     &          mxcmu, mxphi, mxsqt, mxulv, mxumu, nlyr, nmom, nphi,
     &          nstr, ntau, numu
      real      accur, albedo, btemp, fbeam, fisot, phi0, temis, ttemp,
     &          umu0, wvnm
c     ..
c     .. array arguments ..

      real      dtauc( maxcly ), phi( maxphi ),
     &          pmom( 0:maxmom, maxcly ), ssalb( maxcly ),
     &          tauc( 0:mxcly ), temper( 0:maxcly ), umu( maxumu ),
     &          utau( maxulv )
c     ..
c     .. local scalars ..

      logical   inperr, tf
c      integer   irmu
      integer iu, j, k, lc, lu, numsqt
c      real      flxalb, rmu
      real yessct
      real      surf_pr( * )
c     ..
c     .. external functions ..

      logical   wrtbad, wrtdim
      real      dref
      external  wrtbad, wrtdim, dref
c     ..
c     .. external subroutines ..

      external  errmsg
c     ..
c     .. intrinsic functions ..

      intrinsic abs, max, mod
c     ..


      inperr = .false.

      if( nstr.lt.2 .or. mod( nstr,2 ).ne.0 ) inperr = wrtbad( 'nstr' )

      tf = .true.
      if( nstr.eq.2 )
     &    call errmsg( 'chekin--2 streams not recommended; '//
     &                 'use specialized 2-stream code twostr instead',
     &                 tf)

      if( nlyr.lt.1 ) inperr = wrtbad( 'nlyr' )

      if( nlyr.gt.maxcly ) inperr = wrtbad( 'maxcly' )

      yessct = 0.0

      do 10 lc = 1, nlyr

         if( dtauc( lc ).lt.0.0 ) inperr = wrtbad( 'dtauc' )

         if( ssalb( lc ).lt.0.0 .or. ssalb( lc ).gt.1.0 )
     &       inperr = wrtbad( 'ssalb' )

         yessct = yessct + ssalb( lc )

         if( plank .and. ibcnd.ne.1 ) then

            if( lc.eq.1 .and. temper( 0 ).lt.0.0 )
     &          inperr = wrtbad( 'temper' )

            if( temper( lc ).lt.0.0 ) inperr = wrtbad( 'temper' )

         end if

   10 continue

      if( nmom.lt.0 .or. ( yessct.gt.0.0 .and. nmom.lt.nstr ) )
     &    inperr = wrtbad( 'nmom' )

      if( maxmom.lt.nmom ) inperr = wrtbad( 'maxmom' )

      do 30 lc = 1, nlyr

         do 20 k = 0, nmom

            if( pmom( k,lc ).lt.-1.0 .or. pmom( k,lc ).gt.1.0 )
     &          inperr = wrtbad( 'pmom' )

   20    continue

   30 continue

      if( ibcnd.eq.1 ) then

         if( maxulv.lt.2 ) inperr = wrtbad( 'maxulv' )

      else if( usrtau ) then

         if( ntau.lt.1 ) inperr = wrtbad( 'ntau' )

         if( maxulv.lt.ntau ) inperr = wrtbad( 'maxulv' )

         do 40 lu = 1, ntau

            if( abs( utau( lu )-tauc( nlyr ) ).le.1.e-4 )
     &          utau( lu ) = tauc( nlyr )

            if( utau( lu ).lt.0.0 .or. utau( lu ).gt.tauc( nlyr ) )
     &          inperr = wrtbad( 'utau' )

   40    continue

      else

         if( maxulv.lt.nlyr + 1 ) inperr = wrtbad( 'maxulv' )

      end if


      if( usrang ) then

         if( numu.lt.0 ) inperr = wrtbad( 'numu' )

         if( .not.onlyfl .and. numu.eq.0 ) inperr = wrtbad( 'numu' )

         if( numu.gt.maxumu ) inperr = wrtbad( 'maxumu' )

         if( ibcnd.eq.1 .and. 2*numu.gt.maxumu )
     &       inperr = wrtbad( 'maxumu' )

         do 50 iu = 1, numu

            if( umu( iu ).lt.-1.0 .or. umu( iu ).gt.1.0 .or.
     &          umu( iu ).eq.0.0 ) inperr = wrtbad( 'umu' )

            if( ibcnd.eq.1 .and. umu( iu ).lt.0.0 )
     &          inperr = wrtbad( 'umu' )

            if( iu.gt.1 ) then

               if( umu( iu ).lt.umu( iu-1 ) ) inperr = wrtbad( 'umu' )

            end if

   50    continue

      else

         if( maxumu.lt.nstr ) inperr = wrtbad( 'maxumu' )

      end if


      if( .not.onlyfl .and. ibcnd.ne.1 ) then

         if( nphi.le.0 ) inperr = wrtbad( 'nphi' )

         if( nphi.gt.maxphi ) inperr = wrtbad( 'maxphi' )

         do 60 j = 1, nphi

            if( phi( j ).lt.0.0 .or. phi( j ).gt.360.0 )
     &          inperr = wrtbad( 'phi' )

   60    continue

      end if


      if( ibcnd.lt.0 .or. ibcnd.gt.1 ) inperr = wrtbad( 'ibcnd' )

      if( ibcnd.eq.0 ) then

         if( fbeam.lt.0.0 ) inperr = wrtbad( 'fbeam' )

         if( fbeam.gt.0.0 .and. ( umu0.le.0.0 .or. umu0.gt.1.0 ) )
     &       inperr = wrtbad( 'umu0' )

         if( fbeam.gt.0.0 .and. ( phi0.lt.0.0 .or. phi0.gt.360.0 ) )
     &       inperr = wrtbad( 'phi0' )

         if( fisot.lt.0.0 ) inperr = wrtbad( 'fisot' )

         if( lamber ) then

            if( albedo.lt.0.0 .or. albedo.gt.1.0 )
     &          inperr = wrtbad( 'albedo' )

         else
c                    ** make sure flux albedo at dense mesh of incident
c                    ** angles does not assume unphysical values
c
c  ********note, the following check was commented out because it 
c                dramatically increases the amount of time needed
c                in cases with a non-lambertian surface

c            do 70 irmu = 3, 100
            do 70 irmu = 1,10
c
               rmu  = irmu*0.1
               flxalb = dref(rmu, surf_pr, iref )
c
               if( flxalb.lt.0.0 .or. flxalb.gt.1.0 )
     &             inperr = wrtbad( 'function bdref' )
               if( inperr ) write(*,*) 'at rmu=',rmu,
     &             'function bdref gives flxalb=',flxalb
   70       continue

         end if


      else if( ibcnd.eq.1 ) then

         if( albedo.lt.0.0 .or. albedo.gt.1.0 )
     &       inperr = wrtbad( 'albedo' )
         if( inperr ) write(*,*) 'albedo=',albedo
      end if


      if( plank .and. ibcnd.ne.1 ) then

         if( wvnm.lt.0.0 )
     &       inperr = wrtbad( 'wvnm' )
         if( inperr ) write(*,*) 'wvnm',wvnm
         if( temis.lt.0.0 .or. temis.gt.1.0 ) inperr = wrtbad( 'temis' )
         if( inperr ) write(*,*) 'temis =',temis
         if( btemp.lt.0.0 ) inperr = wrtbad( 'btemp' )
         if( inperr ) write(*,*) 'btemp =',btemp
         if( ttemp.lt.0.0 ) inperr = wrtbad( 'ttemp' )
         if( inperr ) write(*,*) 'ttemp =',ttemp
      end if


      if( accur.lt.0.0 .or. accur.gt.1.e-2 ) inperr = wrtbad( 'accur' )
      if( inperr ) write(*,*) 'accur =',accur

      if( mxcly.lt.nlyr ) inperr = wrtdim( 'mxcly', nlyr )

      if( ibcnd.ne.1 ) then

         if( usrtau .and. mxulv.lt.ntau )
     &       inperr = wrtdim( 'mxulv', ntau )

         if( .not.usrtau .and. mxulv.lt.nlyr + 1 )
     &       inperr = wrtdim( 'mxulv', nlyr + 1 )

      else

         if( mxulv.lt.2 ) inperr = wrtdim( 'mxulv', 2 )

      end if

      if( mxcmu.lt.nstr ) inperr = wrtdim( 'mxcmu', nstr )
      if( inperr ) write(*,*) 'mxcmu,nstr',mxcmu,nstr
      if( usrang .and. mxumu.lt.numu ) inperr = wrtdim( 'mxumu', numu )
      if( inperr ) write(*,*) 'mxumu, numu',mxumu,numu
      if( usrang .and. ibcnd.eq.1 .and. mxumu.lt.2*numu )
     &    inperr = wrtdim( 'mxumu', 2*numu )
      if( .not.usrang .and. mxumu.lt.nstr )
     &    inperr = wrtdim( 'mxumu', nstr )
      if( inperr ) write(*,*) 'mxumu, nstr',mxumu,nstr

      if( .not.onlyfl .and. ibcnd.ne.1 .and. mxphi.lt.nphi )
     &    inperr = wrtdim( 'mxphi', nphi )
      if( inperr ) write(*,*) 'mxphi, nphi',mxphi,nphi
      numsqt = 2*max( 100, nstr )
      if( mxsqt.lt.numsqt ) inperr = wrtdim( 'mxsqt', numsqt )
      if( inperr ) write(*,*) 'mxsqt,numsqt',mxsqt,numsqt
      tf = .true.
      if( inperr )
     &    call errmsg( 'disort--input and/or dimension errors', tf )

      if( plank ) then
         tf = .false.
         do 80 lc = 1, nlyr

c            if( abs( temper( lc )-temper( lc-1 ) ).gt. 10.0 )
            if( abs( temper( lc )-temper( lc-1 ) ).gt. 50.0 )
     &          call errmsg('chekin--vertical temperature step may'
     &                      //' be too large for good accuracy',
     &                      tf )
   80    continue

      end if

      tf = .false.
      if( .not.corint .and. .not.onlyfl .and. fbeam.gt.0.0 .and.
     &    yessct.gt.0.0 .and. deltam )
     &     call errmsg( 'chekin--intensity correction is off; '//
     &                  'intensities may be less accurate', tf )


      return
      end

      real function  dref( mu, surf_pr, iref )

c        flux albedo for given angle of incidence, given
c        a bidirectional reflectivity.
c
c  input :   mu      cosine of incidence angle
c
c            wvnm  wavenumber (inv-cm) of spectral interval
c
c            iref  : bidirectional reflectance options
c             0 - Lambert
c             1 - Hapke's BDR model
c             2 - Breon's BDR model; combination of Li + Roujean
c             3 - Roujean's BDR model
c             4 - Cox and Munk glint model
c
c   SURF_PR : Wavelength dependent surface properties array 
c             IREF= 0 - Lambert albedo
c             IREF= 1 - Hapke : HH, W
c             IREF= 2 - Breon's BDR model: k0, k1, k2
c             IREF= 3 - Roujean's BDR model: k0, k1, k2
c             IREF= 4 - Cox and Munk glint model: n, k, ws, phiw
c
c  internal variables :
c
c       nmug   :  number of angle cosine quadrature points on (-1,1)
c                 for integrating bidirectional reflectivity to get
c                 directional emissivity (it is necessary to use a
c                 quadrature set distinct from the computational angles,
c                 because the computational angles may not be dense
c                 enough -- i.e. 'nstr' may be too small -- to give an
c                 accurate approximation for the integration).
c
c       gmu    :  the 'nmug' angle cosine quadrature points on (0,1)
c
c       gwt    :  the 'nmug' angle cosine quadrature weights on (0,1)
c
c   called by- chekin
c   calls- qgausn, errmsg, bdref
c +--------------------------------------------------------------------+

c     .. parameters ..

      integer   nmug
      parameter ( nmug = 50 )
c     ..
c     .. scalar arguments ..

      real      mu
c
c     .. array arguments
c
      real surf_pr( 4 )
c     ..
c     .. local scalars ..

      logical   pass1
      logical   tf
      integer   iref
      integer   jg, k
      real      pi, sum, pigmu
c     ..
c     .. local arrays ..

      real      gmu( nmug ), gwt( nmug )
c     ..
c     .. external functions ..

      real      bdref
      external  bdref
c     ..
c     .. external subroutines ..

      external  errmsg, qgausn
c     ..
c     .. intrinsic functions ..

      intrinsic abs, asin
c     ..
      save      pass1, gmu, gwt, pi
      data      pass1 / .true. /


      if( pass1 ) then

         pass1 = .false.
         pi   = 2.*asin( 1.0 )

         call qgausn( nmug/2, gmu, gwt )

         do 10 k = 1, nmug / 2
            gmu( k + nmug/2 ) = -gmu( k )
            gwt( k + nmug/2 ) = gwt( k )
   10    continue

      end if

      tf = .true.
      if( abs( mu ).gt.1.0 )
     &    call errmsg( 'dref--input argument error(s)',tf )

      dref = 0.0

c                       ** loop over azimuth angle difference
      do 30 jg = 1, nmug

         sum  = 0.0
c                       ** loop over angle of reflection
         do 20 k = 1, nmug / 2
            pigmu = pi*gmu( jg )
            sum  = sum + gwt( k )*gmu( k )*
     &             bdref( gmu( k ), mu, pigmu, surf_pr, iref )
   20    continue

         dref = dref + gwt( jg )*sum

   30 continue

      if( dref.lt.0.0 .or. dref.gt.1.0 ) then
        tf = .false.
        write(*,'(/,1a,1pe12.4)') ' dref = ',dref
        call errmsg( 'dref--albedo value not in (0,1)',tf )
      endif

      return
      end

      subroutine lepoly( nmu, m, maxmu, twonm1, mu, sqt, ylm )

c       computes the normalized associated legendre polynomial,
c       defined in terms of the associated legendre polynomial
c       plm = p-sub-l-super-m as
c
c             ylm(mu) = sqrt( (l-m)!/(l+m)! ) * plm(mu)
c
c       for fixed order m and all degrees from l = m to twonm1.
c       when m.gt.0, assumes that y-sub(m-1)-super(m-1) is available
c       from a prior call to the routine.
c
c       reference: dave, j.v. and b.h. armstrong, computations of
c                  high-order associated legendre polynomials,
c                  j. quant. spectrosc. radiat. transfer 10,
c                  557-562, 1970.  (hereafter d/a)
c
c       method: varying degree recurrence relationship.
c
c       notes:
c       (1) the d/a formulas are transformed by setting m=n-1; l=k-1.
c       (2) assumes that routine is called first with  m = 0, then with
c           m = 1, etc. up to  m = twonm1.
c
c
c  i n p u t     v a r i a b l e s:
c
c       nmu    :  number of arguments of ylm
c
c       m      :  order of ylm
c
c       maxmu  :  first dimension of ylm
c
c       twonm1 :  max degree of ylm
c
c       mu(i)  :  arguments of ylm (i = 1 to nmu)
c
c       sqt(k) :  square root of k
c
c       if m.gt.0, ylm(m-1,i) for i = 1 to nmu is assumed to exist
c       from a prior call.
c
c
c  o u t p u t     v a r i a b l e:
c
c       ylm(l,i) :  l = m to twonm1, normalized associated legendre
c                   polynomials evaluated at argument mu(i)
c
c   called by- disort, albtrn
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      integer   m, maxmu, nmu, twonm1
c     ..
c     .. array arguments ..

      real      mu( * ), ylm( 0:maxmu, * ), sqt( * )
c     ..
c     .. local scalars ..

      integer   i, l
      real      tmp1, tmp2
c     ..


      if( m.eq.0 ) then
c                             ** upward recurrence for ordinary
c                             ** legendre polynomials
         do 20 i = 1, nmu
            ylm( 0, i ) = 1.0
            ylm( 1, i ) = mu( i )
   20    continue


         do 40 l = 2, twonm1

            do 30 i = 1, nmu
               ylm( l, i ) = ( ( 2*l - 1 )*mu( i )*ylm( l-1, i ) -
     &                         ( l - 1 )*ylm( l-2, i ) ) / l
   30       continue

   40    continue


      else

         do 50 i = 1, nmu
c                               ** y-sub-m-super-m; derived from
c                               ** d/a eqs. (11,12), stwl(58c)

            ylm( m, i ) = - sqt( 2*m - 1 ) / sqt( 2*m )*
     &                      sqrt( 1.- mu(i)**2 )*ylm( m-1, i )

c                              ** y-sub-(m+1)-super-m; derived from
c                              ** d/a eqs.(13,14) using eqs.(11,12),
c                              ** stwl(58f)

            ylm( m+1, i ) = sqt( 2*m + 1 )*mu( i )*ylm( m, i )

   50    continue

c                                   ** upward recurrence; d/a eq.(10),
c                                   ** stwl(58a)
         do 70 l = m + 2, twonm1

            tmp1  = sqt( l - m )*sqt( l + m )
            tmp2  = sqt( l - m - 1 )*sqt( l + m - 1 )

            do 60 i = 1, nmu
               ylm( l, i ) = ( ( 2*l - 1 )*mu( i )*ylm( l-1, i ) -
     &                         tmp2*ylm( l-2, i ) ) / tmp1
   60       continue

   70    continue

      end if


      return
      end


      subroutine pravin( umu, numu, mxumu, utau, ntau, u0u )

c        print azimuthally averaged intensities at user angles

c   called by- disort

c     lenfmt   max number of polar angle cosines umu that can be
c              printed on one line, as set in format statement
c --------------------------------------------------------------------

c     .. scalar arguments ..

      integer   mxumu, ntau, numu
c     ..
c     .. array arguments ..

      real      u0u( mxumu, * ), umu( numu ), utau( ntau )
c     ..
c     .. local scalars ..

      integer   iu, iumax, iumin, lenfmt, lu, np, npass
c     ..
c     .. intrinsic functions ..

      intrinsic min
c     ..


      if( numu.lt.1 )  return

      write( *, '(//,a)' )
     &   ' *******  azimuthally averaged intensities ' //
     &   '(at user polar angles)  ********'

      lenfmt = 8
      npass  = 1 + (numu-1) / lenfmt

      write( *,'(/,a,/,a)') '   optical   polar angle cosines',
     &                      '     depth'

      do 20 np = 1, npass

         iumin  = 1 + lenfmt * ( np - 1 )
         iumax  = min( lenfmt*np, numu )
         write( *,'(/,10x,8f14.5)') ( umu(iu), iu = iumin, iumax )

         do 10 lu = 1, ntau
            write( *, '(0p,f10.4,1p,8e14.4)' ) utau( lu ),
     &           ( u0u( iu,lu ), iu = iumin, iumax )
   10    continue

   20 continue


      return
      end

      subroutine prtinp( nlyr, dtauc, dtaucp, ssalb, nmom, pmom, temper,
     &                   wvnm, ntau, utau, nstr, numu, umu,
     &                   nphi, phi, ibcnd, fbeam, umu0, phi0, fisot,
     &                   lamber, albedo, btemp, ttemp, temis, deltam,
     &                   plank, onlyfl, corint, accur, flyr, lyrcut,
     &                   oprim, tauc, taucpr, maxmom, prtmom )

c        print values of input variables

c   called by- disort
c --------------------------------------------------------------------

c     .. scalar arguments ..

      logical   corint, deltam, lamber, lyrcut, onlyfl, plank, prtmom
      integer   ibcnd, maxmom, nlyr, nmom, nphi, nstr, ntau, numu
      real      accur, albedo, btemp, fbeam, fisot, phi0, temis, ttemp,
     &          umu0, wvnm
c     ..
c     .. array arguments ..

      real      dtauc( * ), dtaucp( * ), flyr( * ), oprim( * ),
     &          phi( * ), pmom( 0:maxmom, * ), ssalb( * ), tauc( 0:* ),
     &          taucpr( 0:* ), temper( 0:* ), umu( * ), utau( * )
c     ..
c     .. local scalars ..

      integer   iu, j, k, lc, lu
      real      yessct
c     ..


      write( *, '(/,a,i4,a,i4)' ) ' no. streams =', nstr,
     &       '     no. computational layers =', nlyr

      if( ibcnd.ne.1 ) write( *, '(i4,a,10f10.4,/,(26x,10f10.4))' )
     &    ntau, ' user optical depths :', ( utau(lu), lu = 1, ntau )

      if( .not.onlyfl ) write( *, '(i4,a,10f9.5,/,(31x,10f9.5))' )
     &    numu, ' user polar angle cosines :', ( umu(iu), iu = 1, numu )

      if( .not.onlyfl .and. ibcnd.ne.1 )
     &    write( *, '(i4,a,10f9.2,/,(28x,10f9.2))' )
     &           nphi,' user azimuthal angles :',( phi(j), j = 1, nphi )

      if( .not.plank .or. ibcnd.eq.1 )
     &    write( *, '(a)' ) ' no thermal emission'


      write( *, '(a,i2)' ) ' boundary condition flag: ibcnd =', ibcnd

      if( ibcnd.eq.0 ) then

         write( *, '(a,1p,e11.3,a,0p,f8.5,a,f7.2,/,a,1p,e11.3)' )
     &          '    incident beam with intensity =', fbeam,
     &          ' and polar angle cosine = ', umu0,
     &          '  and azimuth angle =', phi0,
     &          '    plus isotropic incident intensity =', fisot

         if( lamber ) write( *, '(a,0p,f8.4)' )
     &                '    bottom albedo (lambertian) =', albedo

         if( .not.lamber ) write( *, '(a)' )
     &       '    bidirectional reflectivity at bottom'

         if( plank ) write( *, '(a,f14.4,/,a,f10.2,a,f10.2,a,f8.4)' )
     &       '    thermal emission at wavenumber: ', wvnm,
     &       '    bottom temperature =', btemp,
     &       '    top temperature =', ttemp,
     &       '    top emissivity =', temis

      else if( ibcnd.eq.1 ) then

         write( *, '(a)' )
     &          '    isotropic illumination from top and bottom'
         write( *, '(a,0p,f8.4)' )
     &          '    bottom albedo (lambertian) =', albedo

      end if


      if( deltam ) write( *, '(a)' ) ' uses delta-m method'
      if( .not.deltam ) write( *, '(a)' ) ' does not use delta-m method'

      if( corint ) write( *, '(a)' ) ' uses tms/ims method'
      if( .not.corint ) write( *,'(a)' ) ' does not use tms/ims method'


      if( ibcnd.eq.1 ) then

         write( *, '(a)' ) ' calculate albedo and transmissivity of'//
     &                     ' medium vs. incident beam angle'

      else if( onlyfl ) then

         write( *, '(a)' )
     &          ' calculate fluxes only'

      else

         write( *, '(a)' ) ' calculate fluxes and intensities'

      end if

      write( *, '(a,1p,e11.2)' )
     &       ' relative convergence criterion for azimuth series =',
     &       accur

      if( lyrcut ) write( *, '(a)' )
     &    ' sets radiation = 0 below absorption optical depth 10'


c                                    ** print layer variables
c                                    ** (to read, skip every other line)

      if( plank ) write( *, '(/,37x,a,3(/,2a))' )
     &  '<------------- delta-m --------------->',
     &  '                   total    single                           ',
     &  'total    single',
     &  '       optical   optical   scatter   separated   ',
     &  'optical   optical   scatter    asymm',
     &  '         depth     depth    albedo    fraction     ',
     &  'depth     depth    albedo   factor   temperature'

      if( .not.plank ) write( *, '(/,37x,a,3(/,2a))' )
     &  '<------------- delta-m --------------->',
     &  '                   total    single                           ',
     &  'total    single',
     &  '       optical   optical   scatter   separated   ',
     &  'optical   optical   scatter    asymm',
     &  '         depth     depth    albedo    fraction     ',
     &  'depth     depth    albedo   factor'


      yessct = 0.0

      do 10 lc = 1, nlyr

         yessct = yessct + ssalb( lc )
c                                       ** f90 nonadvancing i/o would
c                                       ** simplify this a lot (also the
c                                       ** two writes above)
         if( plank )
     &       write( *,'(i4,2f10.4,f10.5,f12.5,2f10.4,f10.5,f9.4,f14.3)')
     &             lc, dtauc( lc ), tauc( lc ), ssalb( lc ), flyr( lc ),
     &             dtaucp( lc ), taucpr( lc ), oprim( lc ),
     &             pmom( 1,lc ), temper( lc-1 )

         if( .not.plank )
     &       write( *,'(i4,2f10.4,f10.5,f12.5,2f10.4,f10.5,f9.4)' )
     &             lc, dtauc( lc ), tauc( lc ), ssalb( lc ), flyr( lc ),
     &             dtaucp( lc ), taucpr( lc ), oprim( lc ), pmom( 1,lc )
   10 continue

      if( plank ) write( *, '(85x,f14.3)' ) temper( nlyr )


      if( prtmom .and. yessct.gt.0.0 ) then

         write( *, '(/,a,i5)' ) ' number of phase function moments = ',
     &        nmom + 1
         write( *, '(a)' ) ' layer   phase function moments'

         do 20 lc = 1, nlyr

            if( ssalb( lc ).gt.0.0 )
     &          write( *, '(i6,10f11.6,/,(6x,10f11.6))' )
     &                 lc, ( pmom( k, lc ), k = 0, nmom )
   20    continue

      end if


      return
      end

      subroutine prtint( uu, utau, ntau, umu, numu, phi, nphi, maxulv,
     &                   maxumu )

c         prints the intensity at user polar and azimuthal angles

c     all arguments are disort input or output variables

c   called by- disort

c     lenfmt   max number of azimuth angles phi that can be printed
c                on one line, as set in format statement
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      integer   maxulv, maxumu, nphi, ntau, numu
c     ..
c     .. array arguments ..

      real      phi( * ), umu( * ), utau( * ), uu( maxumu, maxulv, * )
c     ..
c     .. local scalars ..

      integer   iu, j, jmax, jmin, lenfmt, lu, np, npass
c     ..
c     .. intrinsic functions ..

      intrinsic min
c     ..


      if( nphi.lt.1 )  return

      write( *, '(//,a)' )
     &   ' *********  i n t e n s i t i e s  *********'

      lenfmt = 10
      npass  = 1 + (nphi-1) / lenfmt

      write( *, '(/,a,/,a,/,a)' )
     &   '             polar   azimuth angles (degrees)',
     &   '   optical   angle',
     &   '    depth   cosine'

      do 30 lu = 1, ntau

         do 20 np = 1, npass

            jmin   = 1 + lenfmt * ( np - 1 )
            jmax   = min( lenfmt*np, nphi )

            write( *, '(/,18x,10f11.2)' ) ( phi(j), j = jmin, jmax )

            if( np.eq.1 ) write( *, '(f10.4,f8.4,1p,10e11.3)' )
     &             utau(lu), umu(1), (uu(1, lu, j), j = jmin, jmax)
            if( np.gt.1 ) write( *, '(10x,f8.4,1p,10e11.3)' )
     &                       umu(1), (uu(1, lu, j), j = jmin, jmax)

            do 10 iu = 2, numu
               write( *, '(10x,f8.4,1p,10e11.3)' )
     &                 umu( iu ), ( uu( iu, lu, j ), j = jmin, jmax )
   10       continue

   20    continue

   30 continue


      return
      end

      subroutine qgausn( m, gmu, gwt )

c       compute weights and abscissae for ordinary gaussian quadrature
c       on the interval (0,1);  that is, such that

c           sum(i=1 to m) ( gwt(i) f(gmu(i)) )

c       is a good approximation to

c           integral(0 to 1) ( f(x) dx )

c   input :    m       order of quadrature rule

c   output :  gmu(i)   array of abscissae (i = 1 to m)
c             gwt(i)   array of weights (i = 1 to m)

c   reference:  davis, p.j. and p. rabinowitz, methods of numerical
c                   integration, academic press, new york, pp. 87, 1975

c   method:  compute the abscissae as roots of the legendre
c            polynomial p-sub-m using a cubically convergent
c            refinement of newton's method.  compute the
c            weights from eq. 2.7.3.8 of davis/rabinowitz.  note
c            that newton's method can very easily diverge; only a
c            very good initial guess can guarantee convergence.
c            the initial guess used here has never led to divergence
c            even for m up to 1000.

c   accuracy:  relative error no better than tol or computer
c              precision (machine epsilon), whichever is larger

c   internal variables:

c    iter      : number of newton method iterations
c    maxit     : maximum allowed iterations of newton method
c    pm2,pm1,p : 3 successive legendre polynomials
c    ppr       : derivative of legendre polynomial
c    p2pri     : 2nd derivative of legendre polynomial
c    tol       : convergence criterion for legendre poly root iteration
c    x,xi      : successive iterates in cubically-convergent version
c                of newtons method (seeking roots of legendre poly.)

c   called by- dref, setdis, surfac
c   calls- d1mach, errmsg
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical tf
      integer   m
c     ..
c     .. array arguments ..

      real      gmu( m ), gwt( m )
c     ..
c     .. local scalars ..

      integer   iter, k, lim, maxit, nn, np1
      real      cona, pi, t
      double precision en, nnp1, one, p, p2pri, pm1, pm2, ppr, prod,
     &                 tmp, tol, two, x, xi
c     ..
c     .. external functions ..

      double precision d1mach
      external  d1mach
c     ..
c     .. external subroutines ..

      external  errmsg
c     ..
c     .. intrinsic functions ..

      intrinsic abs, asin, cos, float, mod, tan
c     ..
      save      pi, tol

      data      pi / 0.0 / , maxit / 1000 / , one / 1.d0 / ,
     &          two / 2.d0 /


      if( pi.eq.0.0 ) then

         pi   = 2.*asin( 1.0 )
         tol  = 10.*d1mach( 4 )

      end if

      tf = .true.
      if( m.lt.1 ) call errmsg( 'qgausn--bad value of m',tf)

      if( m.eq.1 ) then

         gmu( 1 ) = 0.5
         gwt( 1 ) = 1.0
         return

      end if

      en   = m
      np1  = m + 1
      nnp1 = m*np1
      cona = float( m - 1 ) / ( 8*m**3 )

      lim  = m / 2

      do 30 k = 1, lim
c                                        ** initial guess for k-th root
c                                        ** of legendre polynomial, from
c                                        ** davis/rabinowitz (2.7.3.3a)
         t  = ( 4*k - 1 )*pi / ( 4*m + 2 )
         x  = cos( t + cona / tan( t ) )
         iter = 0
c                                        ** upward recurrence for
c                                        ** legendre polynomials
   10    continue
         iter   = iter + 1
         pm2    = one
         pm1    = x

         do 20 nn = 2, m
            p    = ( ( 2*nn - 1 )*x*pm1 - ( nn - 1 )*pm2 ) / nn
            pm2  = pm1
            pm1  = p
   20    continue
c                                              ** newton method
         tmp    = one / ( one - x**2 )
         ppr    = en*( pm2 - x*p )*tmp
         p2pri  = ( two*x*ppr - nnp1*p )*tmp
         xi     = x - ( p / ppr )*( one +
     &            ( p / ppr )*p2pri / ( two*ppr ) )

c                                              ** check for convergence
         if( abs( xi - x ).gt.tol ) then
            tf = .true.
            if( iter.gt.maxit )
     &          call errmsg( 'qgausn--max iteration count',tf)

            x  = xi
            go to  10

         end if
c                             ** iteration finished--calculate weights,
c                             ** abscissae for (-1,1)
         gmu( k ) = -x
         gwt( k ) = two / ( tmp*( en*pm2 )**2 )
         gmu( np1 - k ) = -gmu( k )
         gwt( np1 - k ) = gwt( k )
   30 continue
c                                    ** set middle abscissa and weight
c                                    ** for rules of odd order
      if( mod( m,2 ).ne.0 ) then

         gmu( lim + 1 ) = 0.0
         prod   = one

         do 40 k = 3, m, 2
            prod   = prod * k / ( k - 1 )
   40    continue

         gwt( lim + 1 ) = two / prod**2
      end if

c                                        ** convert from (-1,1) to (0,1)
      do 50 k = 1, m
         gmu( k ) = 0.5*gmu( k ) + 0.5
         gwt( k ) = 0.5*gwt( k )
   50 continue


      return
      end

      real function ratio( a, b )

c        calculate ratio  a/b  with over- and under-flow protection
c        (thanks to prof. jeff dozier for some suggestions here).
c        since this routine takes two logs, it is no speed demon,
c        but it is invaluable for comparing results from two runs
c        of a program under development.
c
c        note:  in fortran90, built-in functions tiny and huge
c               can replace the r1mach calls.
c
c   called by- disort
c   calls- r1mach
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      real      a, b
c     ..
c     .. local scalars ..

      logical   pass1
      real      absa, absb, huge, powa, powb, powmax, powmin, tiny
c     ..
c     .. external functions ..

      real      r1mach
      external  r1mach
c     ..
c     .. intrinsic functions ..

      intrinsic abs, log10, sign
c     ..
c     .. save statement ..

      save      pass1, tiny, huge, powmax, powmin
c     ..
c     .. data statements ..

      data      pass1 / .true. /
c     ..


      if( pass1 ) then

         tiny   = r1mach( 1 )
         huge   = r1mach( 2 )
         powmax = log10( huge )
         powmin = log10( tiny )
         pass1  = .false.

      end if


      if( a.eq.0.0 ) then

         if( b.eq.0.0 ) then

            ratio  = 1.0

         else

            ratio  = 0.0

         end if


      else if( b.eq.0.0 ) then

         ratio  = sign( huge, a )

      else

         absa   = abs( a )
         absb   = abs( b )
         powa   = log10( absa )
         powb   = log10( absb )

         if( absa.lt.tiny .and. absb.lt.tiny ) then

            ratio  = 1.0

         else if( powa - powb.ge.powmax ) then

            ratio  = huge

         else if( powa - powb.le.powmin ) then

            ratio  = tiny

         else

            ratio  = absa / absb

         end if
c                      ** dont use old trick of determining sign
c                      ** from a*b because a*b may (over/under)flow

         if( ( a.gt.0.0 .and. b.lt.0.0 ) .or.
     &       ( a.lt.0.0 .and. b.gt.0.0 ) ) ratio = -ratio

      end if


      return
      end

      subroutine slftst( corint, accur, albedo, btemp, deltam, dtauc,
     &                   fbeam, fisot, ibcnd, lamber, nlyr, plank, nphi,
     &                   numu, nstr, ntau, onlyfl, phi, phi0, nmom,
     &                   pmom, prnt, prntu0, ssalb, temis, temper,
     &                   ttemp, umu, usrang, usrtau, utau, umu0, wvnm, 
     &                   compar, flup, rfldir, rfldn, uu )

c       if  compar = false, save user input values that would otherwise
c       be destroyed and replace them with input values for self-test.
c       if  compar = true, compare self-test case results with correct
c       answers and restore user input values if test is passed.
c
c       (see file 'disort.doc' for variable definitions.)
c
c
c     i n t e r n a l    v a r i a b l e s:
c
c         acc     relative accuracy required for passing self-test
c
c         errorn  relative errors in disort output variables
c
c         ok      logical variable for determining failure of self-test
c
c         all variables ending in 's' are temporary 's'torage for input
c
c   called by- disort
c   calls- tstbad, errmsg
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical   compar, corint, deltam, lamber, onlyfl, plank, usrang,
     &          usrtau
      integer   ibcnd, nlyr, nmom, nphi, nstr, ntau, numu
      real      accur, albedo, btemp, dtauc, fbeam, fisot, flup, phi,
     &          phi0, rfldir, rfldn, ssalb, temis, ttemp, umu, umu0,
     &          utau, uu, wvnm
c     ..
c     .. array arguments ..

      logical   prnt( * ), prntu0( * )
      logical   tf
      real      pmom( 0:* ), temper( 0:* )
c     ..
c     .. local scalars ..

      logical   corins, deltas, lambes, ok, onlyfs, planks, usrans,
     &          usrtas
      integer   i, ibcnds, n, nlyrs, nmoms, nphis, nstrs, ntaus, numus
      real      acc, accurs, albeds, btemps, dtaucs, error1, error2,
     &          error3, error4, fbeams, fisots, phi0s, phis, ssalbs,
     &          temiss, ttemps, umu0s, umus, utaus, wvnms
c     ..
c     .. local arrays ..

      logical   prnts( 5 ), prnu0s( 2 )
      real      pmoms( 0:4 ), tempes( 0:1 )
c     ..
c     .. external functions ..

      logical   tstbad
      external  tstbad
c     ..
c     .. external subroutines ..

      external  errmsg
c     ..
c     .. intrinsic functions ..

      intrinsic abs
c     ..
      save
      data      acc / 1.e-4 /


      if( .not.compar ) then
c                                     ** save user input values
         nlyrs  = nlyr
         dtaucs = dtauc
         ssalbs = ssalb

         do 10 n = 0, 4
            pmoms( n ) = pmom( n )
   10    continue

         nstrs  = nstr
         nmoms  = nmom
         usrans = usrang
         numus  = numu
         umus   = umu
         usrtas = usrtau
         ntaus  = ntau
         utaus  = utau
         nphis  = nphi
         phis   = phi
         ibcnds = ibcnd
         fbeams = fbeam
         umu0s  = umu0
         phi0s  = phi0
         fisots = fisot
         lambes = lamber
         albeds = albedo
         deltas = deltam
         onlyfs = onlyfl
         corins = corint
         accurs = accur
         planks = plank
         wvnms = wvnm
         btemps = btemp
         ttemps = ttemp
         temiss = temis
         tempes( 0 ) = temper( 0 )
         tempes( 1 ) = temper( 1 )

         do 20 i = 1, 5
            prnts( i ) = prnt( i )
   20    continue

         do 30 i = 1, 2
            prnu0s( i ) = prntu0( i )
   30    continue

c                                     ** set input values for self-test
         nstr   = 4
         nlyr   = 1
         dtauc  = 1.0
         ssalb  = 0.9
         nmom   = 4
c                          ** haze l moments
         pmom( 0 ) = 1.0
         pmom( 1 ) = 0.8042
         pmom( 2 ) = 0.646094
         pmom( 3 ) = 0.481851
         pmom( 4 ) = 0.359056
         usrang = .true.
         numu   = 1
         umu    = 0.5
         usrtau = .true.
         ntau   = 1
         utau   = 0.5
         nphi   = 1
         phi    = 90.0
         ibcnd  = 0
         fbeam  = 3.14159265
         umu0   = 0.866
         phi0   = 0.0
         fisot  = 1.0
         lamber = .true.
         albedo = 0.7
         deltam = .true.
         onlyfl = .false.
         corint = .true.
         accur  = 1.e-4
         plank  = .true.
         wvnm = 0.0
         btemp  = 300.0
         ttemp  = 100.0
         temis  = 0.8
         temper( 0 ) = 210.0
         temper( 1 ) = 200.0

         do 40 i = 1, 5
            prnt( i ) = .false.
   40    continue

         do 50 i = 1, 2
            prntu0( i ) = .false.
   50    continue


      else
c                                    ** compare test case results with
c                                    ** correct answers and abort if bad
         ok     = .true.

         error1 = ( uu - 47.865571 ) / 47.865571
         error2 = ( rfldir - 1.527286 ) / 1.527286
         error3 = ( rfldn - 28.372225 ) / 28.372225
         error4 = ( flup - 152.585284 ) / 152.585284

         if( abs( error1 ).gt.acc ) ok = tstbad( 'uu', error1 )

         if( abs( error2 ).gt.acc ) ok = tstbad( 'rfldir', error2 )

         if( abs( error3 ).gt.acc ) ok = tstbad( 'rfldn', error3 )

         if( abs( error4 ).gt.acc ) ok = tstbad( 'flup', error4 )

         tf = .true.
         if( .not.ok ) call errmsg( 'disort--self-test failed', tf )

c                                      ** restore user input values
         nlyr   = nlyrs
         dtauc  = dtaucs
         ssalb  = ssalbs

         do 60 n = 0, 4
            pmom( n ) = pmoms( n )
   60    continue

         nstr   = nstrs
         nmom   = nmoms
         usrang = usrans
         numu   = numus
         umu    = umus
         usrtau = usrtas
         ntau   = ntaus
         utau   = utaus
         nphi   = nphis
         phi    = phis
         ibcnd  = ibcnds
         fbeam  = fbeams
         umu0   = umu0s
         phi0   = phi0s
         fisot  = fisots
         lamber = lambes
         albedo = albeds
         deltam = deltas
         onlyfl = onlyfs
         corint = corins
         accur  = accurs
         plank  = planks
         wvnm = wvnms
         btemp  = btemps
         ttemp  = ttemps
         temis  = temiss
         temper( 0 ) = tempes( 0 )
         temper( 1 ) = tempes( 1 )

         do 70 i = 1, 5
            prnt( i ) = prnts( i )
   70    continue

         do 80 i = 1, 2
            prntu0( i ) = prnu0s( i )
   80    continue

      end if


      return
      end

      subroutine zeroal( nd1, expbea, flyr, oprim, phasa, phast, phasm,
     &                        taucpr, xr0, xr1,
     &                   nd2, cmu, cwt, psi0, psi1, wk, z0, z1, zj,
     &                   nd3, ylm0,
     &                   nd4, array, cc, evecc,
     &                   nd5, gl,
     &                   nd6, ylmc,
     &                   nd7, ylmu,
     &                   nd8, kk, ll, zz, zplk0, zplk1,
     &                   nd9, gc,
     &                   nd10, layru, utaupr,
     &                   nd11, gu,
     &                   nd12, z0u, z1u, zbeam,
     &                   nd13, eval,
     &                   nd14, amb, apb,
     &                   nd15, ipvt, z,
     &                   nd16, rfldir, rfldn, flup, uavg, dfdt,
     &                   nd17, albmed, trnmed,
     &                   nd18, u0u,
     &                   nd19, uu )

c         zero arrays; ndn is dimension of all arrays following
c         it in the argument list

c   called by- disort
c --------------------------------------------------------------------

c     .. scalar arguments ..

      integer   nd1, nd10, nd11, nd12, nd13, nd14, nd15, nd16, nd17,
     &          nd18, nd19, nd2, nd3, nd4, nd5, nd6, nd7, nd8, nd9
c     ..
c     .. array arguments ..

      integer   ipvt( * ), layru( * )
      real      albmed( * ), amb( * ), apb( * ), array( * ), cc( * ),
     &          cmu( * ), cwt( * ), dfdt( * ), eval( * ), evecc( * ),
     &          expbea( * ), flup( * ), flyr( * ), gc( * ), gl( * ),
     &          gu( * ), kk( * ), ll( * ), oprim( * ), phasa( * ),
     &          phast( * ), phasm( * ), psi0( * ), psi1( * ),
     &          rfldir( * ), rfldn( * ), taucpr( * ), trnmed( * ),
     &          u0u( * ), uavg( * ), utaupr( * ), uu( * ), wk( * ),
     &          xr0( * ), xr1( * ), ylm0( * ), ylmc( * ), z( * ),
     &          z0( * ), z0u( * ), z1( * ), z1u( * ), ylmu( * ),
     &          zbeam( * ), zj( * ), zplk0( * ), zplk1( * ), zz( * )
c     ..
c     .. local scalars ..

      integer   n
c     ..


      do 10 n = 1, nd1
         expbea( n ) = 0.0
         flyr( n )   = 0.0
         oprim( n )  = 0.0
         phasa( n )  = 0.0
         phast( n )  = 0.0
         phasm( n )  = 0.0
         taucpr( n ) = 0.0
         xr0( n )    = 0.0
         xr1( n )    = 0.0
   10 continue

      do 20 n = 1, nd2
         cmu( n )  = 0.0
         cwt( n )  = 0.0
         psi0( n ) = 0.0
         psi1( n ) = 0.0
         wk( n )   = 0.0
         z0( n )   = 0.0
         z1( n )   = 0.0
         zj( n )   = 0.0
   20 continue

      do 30 n = 1, nd3
         ylm0( n ) = 0.0
   30 continue

      do 40 n = 1, nd4
         array( n ) = 0.0
         cc( n )    = 0.0
         evecc( n ) = 0.0
   40 continue

      do 50 n = 1, nd5
         gl( n ) = 0.0
   50 continue

      do 60 n = 1, nd6
         ylmc( n ) = 0.0
   60 continue

      do 70 n = 1, nd7
         ylmu( n ) = 0.0
   70 continue

      do 80 n = 1, nd8
         kk( n )    = 0.0
         ll( n )    = 0.0
         zz( n )    = 0.0
         zplk0( n ) = 0.0
         zplk1( n ) = 0.0
   80 continue

      do 90 n = 1, nd9
         gc( n ) = 0.0
   90 continue

      do 100 n = 1, nd10
         layru( n )  = 0
         utaupr( n ) = 0.0
  100 continue

      do 110 n = 1, nd11
         gu( n ) = 0.0
  110 continue

      do 120 n = 1, nd12
         z0u( n )   = 0.0
         z1u( n )   = 0.0
         zbeam( n ) = 0.0
  120 continue

      do 130 n = 1, nd13
         eval( n ) = 0.0
  130 continue

      do 140 n = 1, nd14
         amb( n ) = 0.0
         apb( n ) = 0.0
  140 continue

      do 150 n = 1, nd15
         ipvt( n ) = 0
         z( n )    = 0.0
  150 continue

      do 160 n = 1, nd16
         rfldir( n ) = 0.
         rfldn( n )  = 0.
         flup( n )   = 0.
         uavg( n )   = 0.
         dfdt( n )   = 0.
  160 continue

      do 170 n = 1, nd17
         albmed( n ) = 0.
         trnmed( n ) = 0.
  170 continue

      do 180 n = 1, nd18
         u0u( n ) = 0.
  180 continue

      do 190 n = 1, nd19
         uu( n ) = 0.
  190 continue


      return
      end

      subroutine zeroit( a, length )

c         zeros a real array a having length elements
c
c   called by- disort, albtrn, solve1, surfac, setmtx, solve0, fluxes
c --------------------------------------------------------------------

c     .. scalar arguments ..

      integer   length
c     ..
c     .. array arguments ..

      real      a( length )
c     ..
c     .. local scalars ..

      integer   l
c     ..


      do 10 l = 1, length
         a( l ) = 0.0
   10 continue


      return
      end

c ******************************************************************
c ********** end of disort service routines ************************
c ******************************************************************

c ******************************************************************
c ********** ibcnd=1 special case routines *************************
c ******************************************************************

      subroutine albtrn( albedo, amb, apb, array, b, bdr, cband, cc,
     &                   cmu, cwt, dtaucp, eval, evecc, gl, gc, gu,
     &                   ipvt, kk, ll, nlyr, nn, nstr, numu, prnt,
     &                   taucpr, umu, u0u, wk, ylmc, ylmu, z, aad,
     &                   evald, eveccd, wkd, mi, mi9m2, maxumu,
     &                   mxcmu, mxumu, nnlyri, sqt, albmed, trnmed )

c    disort special case to get only albedo and transmissivity
c    of entire medium as a function of incident beam angle
c    (many simplifications because boundary condition is just
c    isotropic illumination, there are no thermal sources, and
c    particular solutions do not need to be computed).  see
c    ref. s2 and references therein for details.

c    the basic idea is as follows.  the reciprocity principle leads to
c    the following relationships for a plane-parallel, vertically
c    inhomogeneous medium lacking thermal (or other internal) sources:
c
c       albedo(theta) = u_0(theta) for unit-intensity isotropic
c                       illumination at *top* boundary
c
c       trans(theta) =  u_0(theta) for unit-intensity isotropic
c                       illumination at *bottom* boundary
c
c    where
c
c       albedo(theta) = albedo for beam incidence at angle theta
c       trans(theta) = transmissivity for beam incidence at angle theta
c       u_0(theta) = upward azim-avg intensity at top boundary
c                    at angle theta


c   o u t p u t    v a r i a b l e s:
c
c       albmed(iu)   albedo of the medium as a function of incident
c                    beam angle cosine umu(iu)
c
c       trnmed(iu)   transmissivity of the medium as a function of
c                    incident beam angle cosine umu(iu)


c    i n t e r n a l   v a r i a b l e s:

c       ncd         number of diagonals below/above main diagonal

c       rcond       estimate of the reciprocal condition of matrix
c                   cband; for system  cband*x = b, relative
c                   perturbations in cband and b of size epsilon may
c                   cause relative perturbations in x of size
c                   epsilon/rcond.  if rcond is so small that
c                          1.0 + rcond .eq. 1.0
c                   is true, then cband may be singular to working
c                   precision.

c       cband       left-hand side matrix of linear system eq. sc(5),
c                   scaled by eq. sc(12); in banded form required
c                   by linpack solution routines

c       ncol        number of columns in cband matrix

c       ipvt        integer vector of pivot indices

c       (most others documented in disort)

c   called by- disort
c   calls- lepoly, zeroit, sgbco, soleig, terpev, setmtx, solve1,
c          altrin, spaltr, praltr
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      logical tf
      integer   maxumu, mi, mi9m2, mxcmu, mxumu, nlyr, nn, nnlyri,
     &          nstr, numu
      real      albedo
c     ..
c     .. array arguments ..

      logical   prnt( * )
      integer   ipvt( * )
      real      albmed( maxumu ), amb( mi, mi ), apb( mi, mi ),
     &          array( mxcmu, mxcmu ), b( nnlyri ), bdr( mi, 0:mi ),
     &          cband( mi9m2, nnlyri ), cc( mxcmu, mxcmu ),
     &          cmu( mxcmu ), cwt( mxcmu ), dtaucp( * ), eval( mi ),
     &          evecc( mxcmu, mxcmu ), gc( mxcmu, mxcmu, * ),
     &          gl( 0:mxcmu, * ), gu( mxumu, mxcmu, * ), kk( mxcmu, * ),
     &          ll( mxcmu, * ), sqt( * ), taucpr( 0:* ),
     &          trnmed( maxumu ), u0u( mxumu, * ), umu( maxumu ),
     &          wk( mxcmu ), ylmc( 0:mxcmu, mxcmu ), ylmu( 0:mxcmu, * ),
     &          z( nnlyri )

      double precision aad( mi, mi ), evald( mi ), eveccd( mi, mi ),
     &                 wkd( mxcmu )
c     ..
c     .. local scalars ..

      logical   lamber, lyrcut
      integer   iq, iu, l, lc, mazim, ncd, ncol, ncut
      real      delm0, fisot, rcond, sgn, sphalb, sphtrn
c     ..
c     .. external subroutines ..

      external  altrin, errmsg, lepoly, praltr, setmtx, sgbco, soleig,
     &          solve1, spaltr, terpev, zeroit
c     ..
c     .. intrinsic functions ..

      intrinsic exp
c     ..

      mazim  = 0
      delm0  = 1.0
c                    ** set disort variables that are ignored in this
c                    ** special case but are needed below in argument
c                    ** lists of subroutines shared with general case
      ncut   = nlyr
      lyrcut = .false.
      fisot  = 1.0
      lamber = .true.
c                          ** get legendre polynomials for computational
c                          ** and user polar angle cosines

      call lepoly( numu, mazim, mxcmu, nstr-1, umu, sqt, ylmu )

      call lepoly( nn, mazim, mxcmu, nstr-1, cmu, sqt, ylmc )

c                       ** evaluate legendre polynomials with negative
c                       ** arguments from those with positive arguments;
c                       ** dave/armstrong eq. (15), stwl(59)
      sgn  = -1.0

      do 20 l = mazim, nstr - 1

         sgn  = -sgn

         do 10 iq = nn + 1, nstr
            ylmc( l, iq ) = sgn*ylmc( l, iq - nn )
   10    continue

   20 continue
c                                  ** zero out bottom reflectivity
c                                  ** (albedo is used only in analytic
c                                  ** formulae involving albedo = 0
c                                  ** solutions; eqs 16-17 of ref s2)

      call zeroit( bdr, mi*( mi+1 ) )


c ===================  begin loop on computational layers  =============

      do 30 lc = 1, nlyr

c                                       ** solve eigenfunction problem
c                                       ** in eq. stwj(8b), stwl(23f)

         call soleig( amb, apb, array, cmu, cwt, gl( 0,lc ), mi, mazim,
     &                mxcmu, nn, nstr, ylmc, cc, evecc, eval,
     &                kk( 1,lc ), gc( 1,1,lc ), aad, eveccd, evald,
     &                wkd )

c                          ** interpolate eigenvectors to user angles

         call terpev( cwt, evecc, gl( 0,lc ), gu( 1,1,lc ), mazim,
     &                mxcmu, mxumu, nn, nstr, numu, wk, ylmc, ylmu )

   30 continue

c ===================  end loop on computational layers  ===============


c                      ** set coefficient matrix (cband) of equations
c                      ** combining boundary and layer interface
c                      ** conditions (in band-storage mode required by
c                      ** linpack routines)

      call setmtx( bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk,
     &             lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut,
     &             nnlyri, nn, nstr, taucpr, wk )

c                      ** lu-decompose the coeff. matrix (linpack)

      ncd  = 3*nn - 1
      call sgbco( cband, mi9m2, ncol, ncd, ncd, ipvt, rcond, z )
      tf = .false.
      if( 1.0+rcond .eq. 1.0 )
     &    call errmsg('albtrn--sgbco says matrix near singular',tf)

c                             ** first, illuminate from top; if only
c                             ** one layer, this will give us everything

c                             ** solve for constants of integration in
c                             ** homogeneous solution

      call solve1( b, cband, fisot, 1, ipvt, ll, mi9m2, mxcmu,
     &             ncol, nlyr, nn, nnlyri, nstr )

c                             ** compute azimuthally-averaged intensity
c                             ** at user angles; gives albedo if multi-
c                             ** layer (eq. 9 of ref s2); gives both
c                             ** albedo and transmissivity if single
c                             ** layer (eqs. 3-4 of ref s2)

      call altrin( gu, kk, ll, mxcmu, mxumu, maxumu, nlyr, nn, nstr,
     &             numu, taucpr, umu, u0u, wk )

c                               ** get beam-incidence albedos from
c                               ** reciprocity principle
      do 40 iu = 1, numu / 2
         albmed( iu ) = u0u( iu + numu/2, 1 )
   40 continue


      if( nlyr.eq.1 ) then

         do 50 iu = 1, numu / 2
c                               ** get beam-incidence transmissivities
c                               ** from reciprocity principle (1 layer);
c                               ** flip them end over end to correspond
c                               ** to positive umu instead of negative

            trnmed( iu ) = u0u( numu/2 + 1 - iu, 2 )
     &                     + exp( -taucpr( nlyr ) / umu( iu + numu/2 ) )

   50    continue

      else
c                             ** second, illuminate from bottom
c                             ** (if multiple layers)

         call solve1( b, cband, fisot, 2, ipvt, ll, mi9m2, mxcmu,
     &                ncol, nlyr, nn, nnlyri, nstr )

         call altrin( gu, kk, ll, mxcmu, mxumu, maxumu, nlyr, nn, nstr,
     &                numu, taucpr, umu, u0u, wk )

c                               ** get beam-incidence transmissivities
c                               ** from reciprocity principle
         do 60 iu = 1, numu / 2
            trnmed( iu ) = u0u( iu + numu/2, 1 )
     &                     + exp( -taucpr( nlyr ) / umu( iu + numu/2 ) )
   60    continue

      end if


      if( albedo.gt.0.0 ) then

c                             ** get spherical albedo and transmissivity
         if( nlyr.eq.1 ) then

            call spaltr( cmu, cwt, gc, kk, ll, mxcmu, nlyr,
     &                    nn, nstr, taucpr, sphalb, sphtrn )
         else

            call spaltr( cmu, cwt, gc, kk, ll, mxcmu, nlyr,
     &                    nn, nstr, taucpr, sphtrn, sphalb )
         end if

c                                ** ref. s2, eqs. 16-17 (these eqs. have
c                                ** a simple physical interpretation
c                                ** like that of adding-doubling eqs.)
         do 70 iu = 1, numu

            albmed(iu) = albmed(iu) + ( albedo / (1.-albedo*sphalb) )
     &                                * sphtrn * trnmed(iu)

            trnmed(iu) = trnmed(iu) + ( albedo / (1.-albedo*sphalb) )
     &                                * sphalb * trnmed(iu)
   70    continue

      end if
c                          ** return umu to all positive values, to
c                          ** agree with ordering in albmed, trnmed
      numu  = numu / 2
      do 80 iu = 1, numu
         umu( iu ) = umu( iu + numu )
   80 continue

      if( prnt(4) ) call praltr( umu, numu, albmed, trnmed )


      return
      end

      subroutine altrin( gu, kk, ll, mxcmu, mxumu, maxumu, nlyr, nn,
     &                   nstr, numu, taucpr, umu, u0u, wk )

c       computes azimuthally-averaged intensity at top and bottom
c       of medium (related to albedo and transmission of medium by
c       reciprocity principles; see ref s2).  user polar angles are
c       used as incident beam angles. (this is a very specialized
c       version of usrint)
c
c       ** note **  user input values of umu (assumed positive) are
c                   temporarily in upper locations of  umu  and
c                   corresponding negatives are in lower locations
c                   (this makes gu come out right).  i.e. the contents
c                   of the temporary umu array are:
c
c                     -umu(numu),..., -umu(1), umu(1),..., umu(numu)
c
c
c   i n p u t    v a r i a b l e s:
c
c       gu     :  eigenvectors interpolated to user polar angles
c                   (i.e., g in eq. sc(1), stwl(31ab))
c
c       kk     :  eigenvalues of coeff. matrix in eq. ss(7), stwl(23b)
c
c       ll     :  constants of integration in eq. sc(1), obtained
c                   by solving scaled version of eq. sc(5);
c                   exponential term of eq. sc(12) not included
c
c       nn     :  order of double-gauss quadrature (nstr/2)
c
c       taucpr :  cumulative optical depth (delta-m-scaled)
c
c       (remainder are disort input variables)
c
c
c   o u t p u t    v a r i a b l e:
c
c       u0u  :    diffuse azimuthally-averaged intensity at top and
c                 bottom of medium (directly transmitted component,
c                 corresponding to bndint in usrint, is omitted).
c
c
c   i n t e r n a l    v a r i a b l e s:
c
c       dtau   :  optical depth of a computational layer
c       palint :  non-boundary-forced intensity component
c       utaupr :  optical depths of user output levels (delta-m scaled)
c       wk     :  scratch vector for saving 'exp' evaluations
c       all the exponential factors (i.e., exp1, expn,... etc.)
c       come from the substitution of constants of integration in
c       eq. sc(12) into eqs. s1(8-9).  all have negative arguments.
c
c   called by- albtrn
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      integer   maxumu, mxcmu, mxumu, nlyr, nn, nstr, numu
c     ..
c     .. array arguments ..

      real      gu( mxumu, mxcmu, * ), kk( mxcmu, * ), ll( mxcmu, * ),
     &          taucpr( 0:* ), u0u( mxumu, * ), umu( maxumu ),
     &          wk( mxcmu )
c     ..
c     .. local scalars ..

      integer   iq, iu, iumax, iumin, lc, lu
      real      denom, dtau, exp1, exp2, expn, mu, palint, sgn
c     ..
c     .. local arrays ..

      real      utaupr( 2 )
c     ..
c     .. intrinsic functions ..

      intrinsic abs, exp
c     ..


      utaupr( 1 ) = 0.0
      utaupr( 2 ) = taucpr( nlyr )

      do 50 lu = 1, 2

         if( lu.eq.1 ) then

            iumin  = numu / 2 + 1
            iumax  = numu
            sgn    = 1.0

         else

            iumin  = 1
            iumax  = numu / 2
            sgn    = - 1.0

         end if
c                                   ** loop over polar angles at which
c                                   ** albedos/transmissivities desired
c                                   ** ( upward angles at top boundary,
c                                   ** downward angles at bottom )
         do 40 iu = iumin, iumax

            mu   = umu( iu )
c                                     ** integrate from top to bottom
c                                     ** computational layer
            palint = 0.0

            do 30 lc = 1, nlyr

               dtau   = taucpr( lc ) - taucpr( lc - 1 )
               exp1   = exp( ( utaupr( lu ) - taucpr( lc - 1 ) ) / mu )
               exp2   = exp( ( utaupr( lu ) - taucpr( lc ) ) / mu )

c                                      ** kk is negative
               do 10 iq = 1, nn

                  wk( iq ) = exp( kk( iq,lc )*dtau )
                  denom  = 1.0 + mu*kk( iq, lc )

                  if( abs( denom ).lt.0.0001 ) then
c                                                   ** l'hospital limit
                     expn   = dtau / mu*exp2

                  else

                     expn   = ( exp1*wk( iq ) - exp2 )*sgn / denom

                  end if

                  palint = palint + gu( iu, iq, lc )*ll( iq, lc )*expn

   10          continue

c                                        ** kk is positive
               do 20 iq = nn + 1, nstr

                  denom  = 1.0 + mu*kk( iq, lc )

                  if( abs( denom ).lt.0.0001 ) then

                     expn   = - dtau / mu * exp1

                  else

                     expn = ( exp1 - exp2 * wk(nstr+1-iq) ) *sgn / denom

                  end if

                  palint = palint + gu( iu, iq, lc )*ll( iq, lc )*expn

   20          continue

   30       continue

            u0u( iu, lu ) = palint

   40    continue

   50 continue


      return
      end

      subroutine praltr( umu, numu, albmed, trnmed )

c        print planar albedo and transmissivity of medium
c        as a function of incident beam angle

c   called by- albtrn
c --------------------------------------------------------------------

c     .. parameters ..

      real      dpr
      parameter ( dpr = 180.0 / 3.14159265 )
c     ..
c     .. scalar arguments ..

      integer   numu
c     ..
c     .. array arguments ..

      real      albmed( numu ), trnmed( numu ), umu( numu )
c     ..
c     .. local scalars ..

      integer   iu
c     ..
c     .. intrinsic functions ..

      intrinsic acos
c     ..


      write( *, '(///,a,//,a)' )
     &   ' *******  flux albedo and/or transmissivity of ' //
     &   'entire medium  ********',
     &  ' beam zen ang   cos(beam zen ang)      albedo   transmissivity'

      do 10 iu = 1, numu
         write( *, '(0p,f13.4,f20.6,f12.5,1p,e17.4)' )
     &      dpr*acos( umu( iu ) ), umu( iu ), albmed( iu ), trnmed( iu )
   10 continue


      return
      end

      subroutine solve1( b, cband, fisot, ihom, ipvt, ll, mi9m2, mxcmu,
     &                   ncol, ncut, nn, nnlyri, nstr )

c        construct right-hand side vector b for isotropic incidence
c        (only) on either top or bottom boundary and solve system
c        of equations obtained from the boundary conditions and the
c        continuity-of-intensity-at-layer-interface equations
c
c
c     i n p u t      v a r i a b l e s:
c
c       cband    :  left-hand side matrix of banded linear system
c                   eq. sc(5), scaled by eq. sc(12); assumed already
c                   in lu-decomposed form, ready for linpack solver
c
c       ihom     :  direction of illumination flag (1, top; 2, bottom)
c
c       ncol     :  number of columns in cband
c
c       nn       :  order of double-gauss quadrature (nstr/2)
c
c       (remainder are disort input variables)
c
c
c    o u t p u t     v a r i a b l e s:
c
c       b        :  right-hand side vector of eq. sc(5) going into
c                   sgbsl; returns as solution vector of eq.
c                   sc(12), constants of integration without
c                   exponential term
c
c       ll      :   permanent storage for b, but re-ordered
c
c
c    i n t e r n a l    v a r i a b l e s:
c
c       ipvt     :  integer vector of pivot indices
c       ncd      :  number of diagonals below or above main diagonal
c
c   called by- albtrn
c   calls- zeroit, sgbsl
c +-------------------------------------------------------------------+

c     .. scalar arguments ..

      integer   ihom, mi9m2, mxcmu, ncol, ncut, nn, nnlyri, nstr
      real      fisot
c     ..
c     .. array arguments ..

      integer   ipvt( nnlyri )
      real      b( nnlyri ), cband( mi9m2, nnlyri ), ll( mxcmu, * )
c     ..
c     .. local scalars ..

      integer   i, ipnt, iq, lc, ncd
c     ..
c     .. external subroutines ..

      external  sgbsl, zeroit
c     ..


      call zeroit( b, nnlyri )

      if( ihom.eq.1 ) then
c                             ** because there are no beam or emission
c                             ** sources, remainder of b array is zero
         do 10 i = 1, nn
            b( i )             = fisot
            b( ncol - nn + i ) = 0.0
   10    continue

      else if( ihom.eq.2 ) then

         do 20 i = 1, nn
            b( i )             = 0.0
            b( ncol - nn + i ) = fisot
   20    continue

      end if


      ncd  = 3*nn - 1
      call sgbsl( cband, mi9m2, ncol, ncd, ncd, ipvt, b, 0 )

      do 40 lc = 1, ncut

         ipnt  = lc*nstr - nn

         do 30 iq = 1, nn
            ll( nn + 1 - iq, lc ) = b( ipnt + 1 - iq )
            ll( iq + nn,     lc ) = b( iq + ipnt )
   30    continue

   40 continue


      return
      end

      subroutine spaltr( cmu, cwt, gc, kk, ll, mxcmu, nlyr, nn, nstr,
     &                   taucpr, sflup, sfldn )

c       calculates spherical albedo and transmissivity for the entire
c       medium from the m=0 intensity components
c       (this is a very specialized version of fluxes)
c
c
c    i n p u t    v a r i a b l e s:
c
c       cmu,cwt    abscissae, weights for gauss quadrature
c                  over angle cosine
c
c       kk      :  eigenvalues of coeff. matrix in eq. ss(7)
c
c       gc      :  eigenvectors at polar quadrature angles, sc(1)
c
c       ll      :  constants of integration in eq. sc(1), obtained
c                  by solving scaled version of eq. sc(5);
c                  exponential term of eq. sc(12) not included
c
c       nn      :  order of double-gauss quadrature (nstr/2)
c
c       (remainder are disort input variables)
c
c
c    o u t p u t   v a r i a b l e s:
c
c       sflup   :  up-flux at top (equivalent to spherical albedo due to
c                  reciprocity).  for illumination from below it gives
c                  spherical transmissivity
c
c       sfldn   :  down-flux at bottom (for single layer, equivalent to
c                  spherical transmissivity due to reciprocity)
c
c
c    i n t e r n a l   v a r i a b l e s:
c
c       zint    :  intensity of m=0 case, in eq. sc(1)
c
c   called by- albtrn
c +--------------------------------------------------------------------

c     .. scalar arguments ..

      integer   mxcmu, nlyr, nn, nstr
      real      sfldn, sflup
c     ..
c     .. array arguments ..

      real      cmu( mxcmu ), cwt( mxcmu ), gc( mxcmu, mxcmu, * ),
     &          kk( mxcmu, * ), ll( mxcmu, * ), taucpr( 0:* )
c     ..
c     .. local scalars ..

      integer   iq, jq
      real      zint
c     ..
c     .. intrinsic functions ..

      intrinsic exp
c     ..


      sflup  = 0.0

      do 30 iq = nn + 1, nstr

         zint   = 0.0
         do 10 jq = 1, nn
            zint  = zint + gc( iq, jq, 1 )*ll( jq, 1 )*
     &                     exp( kk( jq,1 )*taucpr( 1 ) )
   10    continue

         do 20 jq = nn + 1, nstr
            zint  = zint + gc( iq, jq, 1 )*ll( jq, 1 )
   20    continue

         sflup  = sflup + cwt( iq - nn )*cmu( iq - nn )*zint

   30 continue


      sfldn  = 0.0

      do 60 iq = 1, nn

         zint   = 0.0
         do 40 jq = 1, nn
            zint  = zint + gc( iq, jq, nlyr )*ll( jq, nlyr )
   40    continue

         do 50 jq = nn + 1, nstr
            zint  = zint + gc( iq, jq, nlyr )*ll( jq, nlyr )*
     &                     exp( - kk( jq,nlyr ) *
     &                     ( taucpr( nlyr ) - taucpr( nlyr-1 ) ) )
   50    continue

         sfldn  = sfldn + cwt( nn + 1 - iq )*cmu( nn + 1 - iq )*zint

   60 continue

      sflup  = 2.0*sflup
      sfldn  = 2.0*sfldn


      return
      end

c ******************************************************************
c ********** end of ibcnd=1 special case routines ******************
c ******************************************************************
