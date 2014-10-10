      subroutine lpmom(func,gwt,gmu,numu,nmom,pmom)
c
cccccccccccccccccccccccccccccc  l p m o m  cccccccccccccccccccccccccccccc
cc                                                                     cc
cc    p u r p o s e :                                                  cc
cc                                                                     cc
cc      compute the legendre polynomial coefficient moments for an     cc
cc      arbitrary function, func, that is defined only at the          cc
cc      gaussian angles, umu.                                          cc
cc                                                                     cc
cc     i n p u t :                                                     cc
cc                                                                     cc
cc   func( n ) :  input angle-dependent quantity at n-th angle         cc
cc    gmu( n ) :  cosine of nth gaussian angle.                        cc
cc    wgt( n ) :  gaussian weights at each angle                       cc
cc        numu :  number of gaussian angles.                           cc
cc        nmom :  number of legendre-polynomial momemts (numu - 1)     cc
cc                                                                     cc
cc                            *** specification of local variables     cc
cc                                                                     cc
cc   pl( n )   : legendre polynomial p-sub-mom  of argument  gmu(n)    cc
cc   plm1( n ) : legendre polynomial p-sub-(mom-1) of argument gmu(n)  cc
cc                                                                     cc
cccccccccccccccccccccccccccccc  l p m o m  cccccccccccccccccccccccccccccc
c
      parameter (maxang = 1001)
c
      integer  nmom, numu
      real     pmom( 0:nmom ), gmu( numu ), gwt( numu ), func( numu )
      real  pl( maxang ), plm1( maxang )
c
      pi = acos(-1.0)
c
      do  40  mom = 0, nmom
c
c****        calculate legendre polys.
c
           if( mom.eq.0 )  then
             do  22  n = 1, numu
                  pl( n ) = 1.0
 22          continue
           else if( mom.eq.1 )  then
                  do  24  n = 1, numu
                       plm1( n ) = 1.0
                       pl( n ) = gmu( n )
 24               continue
           else
             do  26  n = 1, numu
                  pol = ( (2*mom-1) * gmu(n) * pl(n)
     $                     - (mom-1) * plm1(n) ) / mom
                  plm1( n ) = pl( n )
                  pl( n ) = pol
 26          continue
           end if
c
c****        get legendre moments by angular quadrature
c
           pmom(mom) = 0.0
           do  30  n = 1, numu
                pmom(mom) = pmom(mom) + 0.5*gwt(n)*pl(n)*func(n)
 30        continue
c
 40   continue
c
      return
      end
