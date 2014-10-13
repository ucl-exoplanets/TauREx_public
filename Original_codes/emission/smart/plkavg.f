      real function plkavg( wnumlo, wnumhi, t )

c        computes planck function integrated between two wavenumbers
c
c  input :  wnumlo : lower wavenumber (inv cm) of spectral interval
c
c           wnumhi : upper wavenumber
c
c           t      : temperature (k)
c
c  output : plkavg : integrated planck function ( watts/sq m )
c                      = integral (wnumlo to wnumhi) of
c                        2h c**2  nu**3 / ( exp(hc nu/kt) - 1)
c                        (where h=plancks constant, c=speed of
c                         light, nu=wavenumber, t=temperature,
c                         and k = boltzmann constant)
c
c  reference : specifications of the physical world: new value
c                 of the fundamental constants, dimensions/n.b.s.,
c                 jan. 1974
c
c  method :  for wnumlo close to wnumhi, a simpson-rule quadrature
c            is done to avoid ill-conditioning; otherwise
c
c            (1)  for wnumlo or wnumhi small,
c                 integral(0 to wnumlo/hi) is calculated by expanding
c                 the integrand in a power series and integrating
c                 term by term;
c
c            (2)  otherwise, integral(wnumlo/hi to infinity) is
c                 calculated by expanding the denominator of the
c                 integrand in powers of the exponential and
c                 integrating term by term.
c
c  accuracy :  at least 6 significant digits, assuming the
c              physical constants are infinitely accurate
c
c  errors which are not trapped:
c
c      * power or exponential series may underflow, giving no
c        significant digits.  this may or may not be of concern,
c        depending on the application.
c
c      * simpson-rule special case is skipped when denominator of
c        integrand will cause overflow.  in that case the normal
c        procedure is used, which may be inaccurate if the
c        wavenumber limits (wnumlo, wnumhi) are close together.
c
c  local variables
c
c        a1,2,... :  power series coefficients
c        c2       :  h * c / k, in units cm*k (h = plancks constant,
c                      c = speed of light, k = boltzmann constant)
c        d(i)     :  exponential series expansion of integral of
c                       planck function from wnumlo (i=1) or wnumhi
c                       (i=2) to infinity
c        epsil    :  smallest number such that 1+epsil .gt. 1 on
c                       computer
c        ex       :  exp( - v(i) )
c        exm      :  ex**m
c        mmax     :  no. of terms to take in exponential series
c        mv       :  multiples of v(i)
c        p(i)     :  power series expansion of integral of
c                       planck function from zero to wnumlo (i=1) or
c                       wnumhi (i=2)
c        pi       :  3.14159...
c        sigma    :  stefan-boltzmann constant (w/m**2/k**4)
c        sigdpi   :  sigma / pi
c        smallv   :  number of times the power series is used (0,1,2)
c        v(i)     :  c2 * (wnumlo(i=1) or wnumhi(i=2)) / temperature
c        vcut     :  power-series cutoff point
c        vcp      :  exponential series cutoff points
c        vmax     :  largest allowable argument of exp function
c
c   called by- disort
c   calls- r1mach, errmsg
c ----------------------------------------------------------------------

c     .. parameters ..

      real      a1, a2, a3, a4, a5, a6
      parameter ( a1 = 1. / 3., a2 = -1. / 8., a3 = 1. / 60.,
     &          a4 = -1. / 5040., a5 = 1. / 272160.,
     &          a6 = -1. / 13305600. )
c     ..
c     .. scalar arguments ..

      real      t, wnumhi, wnumlo
c     ..
c     .. local scalars ..

      integer   i, k, m, mmax, n, smallv
      real      c2, conc, del, epsil, ex, exm, hh, mv, oldval, pi,
     &          sigdpi, sigma, val, val0, vcut, vmax, vsq
c     ..
c     .. local arrays ..

      real      d( 2 ), p( 2 ), v( 2 ), vcp( 7 )
c     ..
c     .. external functions ..

      real      r1mach
      external  r1mach
c     ..
c     .. external subroutines ..

      external  errmsg
c     ..
c     .. intrinsic functions ..

      intrinsic abs, asin, exp, log, mod
c     ..
c     .. statement functions ..

      real      plkf
c     ..
      save      pi, conc, vmax, epsil, sigdpi

      data      c2 / 1.438786 / , sigma / 5.67032e-8 / , vcut / 1.5 / ,
     &          vcp / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      data      pi / 0.0 /

c     .. statement function definitions ..

      plkf( x ) = x**3 / ( exp( x ) - 1.0 )
c     ..


      if( pi .eq. 0.0 ) then

         pi     = 2.*asin( 1.0 )
         vmax   = log( r1mach( 2 ) )
         epsil  = r1mach( 4 )
         sigdpi = sigma / pi
         conc   = 15. / pi**4

      end if


      if( t.lt.0.0 .or. wnumhi.le.wnumlo .or. wnumlo.lt.0. )
     &    call errmsg('plkavg--temperature or wavenums. wrong',.true.)


      if( t .lt. 1.e-4 ) then

         plkavg = 0.0
         return

      end if


      v( 1 ) = c2*wnumlo / t
      v( 2 ) = c2*wnumhi / t

      if( v( 1 ).gt.epsil .and. v( 2 ).lt.vmax .and.
     &    ( wnumhi - wnumlo ) / wnumhi .lt. 1.e-2 ) then

c                          ** wavenumbers are very close.  get integral
c                          ** by iterating simpson rule to convergence.

         hh     = v( 2 ) - v( 1 )
         oldval = 0.0
         val0   = plkf( v( 1 ) ) + plkf( v( 2 ) )

         do 20 n = 1, 10

            del  = hh / ( 2*n )
            val  = val0

            do 10 k = 1, 2*n - 1
               val  = val + 2*( 1 + mod( k,2 ) )*
     &                      plkf( v( 1 ) + k*del )
   10       continue

            val  = del / 3.*val
            if( abs( ( val - oldval ) / val ).le.1.e-6 ) go to  30
            oldval = val

   20    continue

         call errmsg( 'plkavg--simpson rule didnt converge',.false.)

   30    continue

         plkavg = sigdpi * t**4 * conc * val

         return

      end if

c                          *** general case ***
      smallv = 0

      do 60 i = 1, 2

         if( v( i ).lt.vcut ) then
c                                   ** use power series
            smallv = smallv + 1
            vsq    = v( i )**2
            p( i ) = conc*vsq*v( i )*( a1 +
     &               v( i )*( a2 + v( i )*( a3 + vsq*( a4 + vsq*( a5 +
     &               vsq*a6 ) ) ) ) )

         else
c                      ** use exponential series
            mmax  = 0
c                                ** find upper limit of series
   40       continue
            mmax  = mmax + 1

            if( v(i) .lt. vcp( mmax ) ) go to  40

            ex     = exp( - v(i) )
            exm    = 1.0
            d( i ) = 0.0

            do 50 m = 1, mmax
               mv     = m*v( i )
               exm    = ex*exm
               d( i ) = d( i ) + exm*( 6.+ mv*( 6.+ mv*( 3.+ mv ) ) )
     &                  / m**4
   50       continue

            d( i ) = conc*d( i )

         end if

   60 continue

c                              ** handle ill-conditioning
      if( smallv.eq.2 ) then
c                                    ** wnumlo and wnumhi both small
         plkavg = p( 2 ) - p( 1 )

      else if( smallv.eq.1 ) then
c                                    ** wnumlo small, wnumhi large
         plkavg = 1.- p( 1 ) - d( 2 )

      else
c                                    ** wnumlo and wnumhi both large
         plkavg = d( 1 ) - d( 2 )

      end if

      plkavg = sigdpi * t**4 * plkavg

c      if( plkavg.eq.0.0 )
c     &    call errmsg('plkavg--returns zero; possible underflow',
c     &    .false.)


      return
      end
