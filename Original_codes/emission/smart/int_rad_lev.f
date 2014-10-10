      subroutine int_rad_lev(lsolar,lplanck,nphi,numu,nz,
     -                       nlyr,nlout,levout,k_out,dp_dp,
     -                       upsflx,dnsflx,dirsflx,uptflx,dntflx,
     -                       sol_rad,th_rad,ups,dns,dirs,
     -                       upth,dnth,rad_s,rad_th)
c
cccccccccccccccccccccccc  i n t _ r a d _ l e v  ccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine interpolates radiances and fluxes to their      cc
cc    output levels, k_out.                                           cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc       nlyr - number of computational model layers                  cc
cc       nphi - number of output azimuth angles                       cc
cc       numu - number of zenith angles in input file                 cc
cc     uptflx - wn-dependent upward thermal flux                      cc
cc     dntflx - wn-dependent downward thermal flux                    cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc     upsflx - wn-dependent upward solar flux                        cc
cc     dnsflx - wn-dependent downward diffuse + direct solar flux     cc
cc    dirsflx - wn-dependent downward direct solar flux               cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    uptflx0 - wn-dependent upward thermal flux                      cc
cc    dntflx0 - wn-dependent downward thermal flux                    cc
cc    th_rad0 - wn-depndent, angle-dependent thermal radiances        cc
cc    upsflx0 - wn-dependent upward solar flux                        cc
cc    dnsflx0 - wn-dependent downward diffuse+direct solar flux       cc
cc   dirsflx0 - wn-dependent downward direct solar flux               cc
cc   sol_rad0 - wn-depndent, angle-dependent solar radiances          cc
cc              output stream                                         cc
cc                                                                    cc
cccccccccccccccccccccccc  i n t _ r a d _ l e v  ccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar
c
      integer nphi,numu,nlyr,nlout,levout
     
      integer nz,k,nze,naz
      integer k_out(mxlout)
c
      real dp_dp(mxlout)
c
c*****   fluxes interpolated to the output levels
c
      real ups(mxlout),dns(mxlout),dirs(mxlout),
     -        upth(mxlout),dnth(mxlout)
      real rad_s(mxumu,mxphi,mxlout),rad_th(mxumu,mxphi,mxlout)
c 
c****    output fluxes and radiances
c
      real dirsflx(mxrad),upsflx(mxrad),
     -     dnsflx(mxrad),uptflx(mxrad),dntflx(mxrad)
c
      real sol_rad(mxumu,mxphi,mxrad),th_rad(mxumu,mxphi,mxrad)
c
      if(lsolar) then
c
c****        interpolate solar radiances and fluxes to output levels
c
        do 4281 k=1,nlout
            if(levout .le. 3) then
c
c****             load fluxes at the appropriate level
c
              ups(k) = upsflx(k_out(k))
              dns(k) = dnsflx(k_out(k))
              dirs(k) = dirsflx(k_out(k))               
c
c****             load solar radiances at the appropriate level
c
              do 4211 naz=1,nphi
                  do 4201 nze=1,numu
                      rad_s(nze,naz,k) = sol_rad(nze,naz,k_out(k))
4201              continue
4211          continue
c
            else
c
c****             interpolate fluxes to the output level
c
              if(k_out(k) .ne. 1 .and. k_out(k) .ne. nlyr+1) then
                ups(k) = upsflx(k_out(k)) - dp_dp(k)*
     -                     (upsflx(k_out(k)) - 
     -                      upsflx(k_out(k)-1))
                dns(k) = dnsflx(k_out(k)) - dp_dp(k)*
     -                     (dnsflx(k_out(k)) - 
     -                      dnsflx(k_out(k)-1))
                dirs(k) = dirsflx(k_out(k)) - dp_dp(k)*
     -                     (dirsflx(k_out(k)) - 
     -                      dirsflx(k_out(k)-1))
c
c****               interpolate solar radiances to the output level
c
                do 4241 naz=1,nphi
                    do 4231 nze=1,numu
                        rad_s(nze,naz,k) = 
     -                        sol_rad(nze,naz,k_out(k)) - 
     -                        dp_dp(k)*(sol_rad(nze,naz,k_out(k)) -
     -                                  sol_rad(nze,naz,k_out(k)-1))
4231                continue
4241            continue
c
              else
c
c****              load solar fluxes and radiances into the
c                  top and/or bottlm level
c
                ups(k) = upsflx(k_out(k))
                dns(k) = dnsflx(k_out(k))
                dirs(k) = dirsflx(k_out(k))
c
c****               load solar radiances at the appropriate level
c
                do 4261 naz=1,nphi
                    do 4251 nze=1,numu
                        rad_s(nze,naz,k) = sol_rad(nze,naz,k_out(k))
4251                continue
4261            continue
              endif
c
            endif
4281    continue
c
      endif
c
c****       t h e r m a l    f l u x e s    a n d    r a d i a n c e s 
c           (do only for nz = 1)
c
      if(lplanck .and. nz .eq. 1) then
c
c****        interpolate thermal radiances and fluxes to output levels
c
        do 4461 k=1,nlout
c
            if(levout .le. 3) then
c
c****             load thermal fluxes to the appropriate level
c
              upth(k) = uptflx(k_out(k))
              dnth(k) = dntflx(k_out(k))
c
c****           interpolate thermal radiances
c
              do 4411 naz=1,nphi
                  do 4401 nze=1,numu
                      rad_th(nze,naz,k) =  th_rad(nze,naz,k_out(k))
4401              continue
4411          continue
c
            else
c
c****         interpolate fluxes and radiances to the output level
c
              if(k_out(k) .ne. 1 .and. k_out(k) .ne. nlyr+1) then
c
                upth(k) = uptflx(k_out(k)) - dp_dp(k)*
     -                   (uptflx(k_out(k)) - uptflx(k_out(k)-1))
                dnth(k) = dntflx(k_out(k)) - dp_dp(k)*
     -                   (dntflx(k_out(k)) - dntflx(k_out(k)-1))
c
c****               interpolate thermal radiances
c
                do 4431 naz=1,nphi
                    do 4421 nze=1,numu
                        rad_th(nze,naz,k) = 
     -                        th_rad(nze,naz,k_out(k)) -
     -                        dp_dp(k)*(th_rad(nze,naz,k_out(k)) - 
     -                                  th_rad(nze,naz,k_out(k)-1))
4421                continue
4431            continue
c
              else
c
c****              load thermal fluxes and radiances into the
c                  top and/or bottom level
c
                upth(k) = uptflx(k_out(k))
                dnth(k) = dntflx(k_out(k))
c
                do 4451 naz=1,nphi
                    do 4441 nze=1,numu
                        rad_th(nze,naz,k) = th_rad(nze,naz,k_out(k))
4441                continue
4451            continue
c
              endif
            endif
4461    continue
c
      endif
c
      return
      end
