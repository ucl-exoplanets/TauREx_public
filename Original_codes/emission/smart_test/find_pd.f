      subroutine find_pd(lsolar,lplanck,nz0,nlyr,
     -                   nlout,nphi,numu,npd,ipd,
     -                   wn,d_wn,dx,solflx,rad,
     -                   trn_dir,trn_flx,ref_flx,abs_flx,
     -                   dns_src,ups_src,dnt_src,upt_src,
     -                   trndir1_i,trnflx1_i,refflx1_i,absflx1_i,
     -                   dtrndir1_i,dtrnflx1_i,drefflx1_i,dabsflx1_i,
     -                   dn_s_src1,up_s_src1,ddn_s_src1,dup_s_src1,
     -                   dn_t_src1,up_t_src1,ddn_t_src1,dup_t_src1,
     -                   sol_rad1,dsol_rad1,th_rad1,dth_rad1,
     -                   pd_trndir,pd_trnflx,pd_refflx,pd_absflx,
     -                   pd_dns_src,pd_ups_src,pd_dnt_src,pd_upt_src,
     -                   pd_rad)
c
ccccccccccccccccccccccccccc  f i n d _ p d  cccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine finds evaluates the partial derivatives of the  cc
cc    radiances with respect to each variable member of the state     cc
cc    vector.                                                         cc
cc                                                                    cc
cc    NOTE: this version has been modified to produce layer dependent cc
cc    flux source functions and their partial derivatives, rather     cc
cc    than the corresponding flux values.                             cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc        nz0 - index of this solar zenith angle (1 to nsol)          cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc       nlyr - number of computational model layers                  cc
cc       nphi - number of output azimuth angles                       cc
cc       numu - number of zenith angles in input file                 cc
cc     uptflx - wn-dependent upward thermal flux                      cc
cc     dntflx - wn-dependent downward thermal flux                    cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc   up_s_src - layer upward flux source function                     cc
cc   dn_s_src - layer downward flux source function                   cc
cc     th_rad - wn-depndent, angle-dependent thermal radiances        cc
cc    sol_rad - wn-depndent, angle-dependent solar radiances          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc                                                                    cc
ccccccccccccccccccccccccccc  f i n d _ p d  cccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      logical lplanck,lsolar
c
      integer nz0,nlyr,nphi,numu,nlout
      integer npd,ipd(mxpd)
      integer naz,nze
c
c****   counters used in map_back and write_mono
c
      integer nz,nlev,l_1,l_2,lpd,k,n
c
      real dx(mxrad,mxpd),solflx
c
      double precision d_wn
      double precision wn
c
c****    flux transmission and absorption functions at wavenuber, wn
c
      real trn_dir(mxrad),trn_flx(mxrad),ref_flx(mxrad),abs_flx(mxrad)
c
      real dns_src(mxrad),ups_src(mxrad),dnt_src(mxrad),upt_src(mxrad)
c
c****   flux source functions at wavenumber wn_io for perturbed state
c
      double precision up_s_src1(mxrad,mxpd,2,nsol),
     -                 dn_s_src1(mxrad,mxpd,2,nsol)
      double precision up_t_src1(mxrad,mxpd,2),dn_t_src1(mxrad,mxpd,2)
c
c****   flux source functions wavenumber derivatives
c
      double precision dup_s_src1(mxrad,mxpd,nsol),
     -                 ddn_s_src1(mxrad,mxpd,nsol)
      double precision dup_t_src1(mxrad,mxpd),ddn_t_src1(mxrad,mxpd)
c
c****     monochormatic radiance for basic state
c 
      real rad(mxumu,mxphi,mxlout)
c
c****    flux transmission and absorption needed for simplified
c        adding method at wavenumber wnio
c
      double precision trndir1_i(mxrad,mxpd,2,nsol),
     -                 trnflx1_i(mxrad,mxpd,2),
     -                 refflx1_i(mxrad,mxpd,2),
     -                 absflx1_i(mxrad,mxpd,2),
     -                 dtrndir1_i(mxrad,mxpd,nsol),
     -                 dtrnflx1_i(mxrad,mxpd),
     -                 drefflx1_i(mxrad,mxpd),
     -                 dabsflx1_i(mxrad,mxpd)
c
c****     perturrbed monochormatic radiance and wn interpolation values
c 
      double precision sol_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,2,nsol),
     -                 dsol_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,nsol)
      double precision th_rad1(mxumu,mxphi,mxlout,mxrad,mxpd,2),
     -                 dth_rad1(mxumu,mxphi,mxlout,mxrad,mxpd)
c
c****    flux transmission and absorption partial 
c        derivatives for simplified adding method at wavenumber wn
c
      real pd_trndir(mxrad,mxpd),
     -     pd_trnflx(mxrad,mxpd),
     -     pd_refflx(mxrad,mxpd),
     -     pd_absflx(mxrad,mxpd)
c
c****   output flux and radiance partial derivatives
c
      real pd_ups_src(mxrad,mxpd),pd_dns_src(mxrad,mxpd),
     -     pd_upt_src(mxrad,mxpd),pd_dnt_src(mxrad,mxpd),
     -     pd_rad(mxumu,mxphi,mxlout,mxrad,mxpd)
c
      double precision trndir,trnflx,refflx,absflx
      double precision dns_src1,ups_src1
c
      real rad_s1,rad_th1,rad1(mxumu,mxphi)
c
      real dnt_src1,upt_src1
c
      real dxi(kp,mxpd),dwn
c
c****    specify the solar zenith angle and gas index
c
      nz = nz0
      nlev = nlyr+1
      dwn = real(d_wn)
c
c****       evaluate the radiances and partial derivatives
c           at this wavenumber for each optical property 
c           that can vary
c
      do 5001 n=1,npd
c
c****        determine the number of levels at which the state
c            vector variable chages (T=nlev, tau=nlyr, alb=1)
c
          if(abs(ipd(n)) .eq. 1) then
c
c****         surface pressure
c
            l_1 = nlev
            l_2 = nlev
c
          else
c
            if(abs(ipd(n)) .eq. 2) then
c
c****            atmospheric/surface temperature
c
              l_1 = 1
              l_2 = nlev
            else
c
              if(abs(ipd(n)) .eq. 3 .or. abs(ipd(n)) .eq. 4) then
c
c****             gas or aerosol optical depth
c
                l_1 = 1
                l_2 = nlev
c
              else
c
c****             this is surface albedo
c
                l_1 = nlev
                l_2 = nlev
              endif
            endif
          endif
c
c****      enter the loop over perturbation level
c
          do 4981 lpd=l_1,l_2
c
c*****         define the inverse of the state vector change
c
              if(dx(lpd,n) .ne. 0.0) then
                dxi(lpd,n) = 1.0/dx(lpd,n)
              else
                dxi(lpd,n) = 0.0
              endif
c
c****           interpolate flux transmision and absorption functions 
c               and jacobians to this wavenumber
c
              if(lpd .le. nlev) then
                if(ipd(n) .lt. 0) then
                  trnflx = trnflx1_i(lpd,n,1) + dtrnflx1_i(lpd,n)*d_wn
                  refflx = refflx1_i(lpd,n,1) + drefflx1_i(lpd,n)*d_wn
                  absflx = absflx1_i(lpd,n,1) + dabsflx1_i(lpd,n)*d_wn
c
                  pd_trnflx(lpd,n) = 
     -               real((trnflx - trn_flx(lpd))*dxi(lpd,n))
                  pd_refflx(lpd,n) = 
     -               real((refflx - ref_flx(lpd))*dxi(lpd,n))
                  pd_absflx(lpd,n) = 
     -               real((absflx - abs_flx(lpd))*dxi(lpd,n))
                else
                  pd_trnflx(lpd,n) = 0.0
                  pd_refflx(lpd,n) = 0.0
                  pd_absflx(lpd,n) = 0.0
                endif
              else
                pd_trnflx(lpd,n) = 0.0
                pd_refflx(lpd,n) = 0.0
                pd_absflx(lpd,n) = 0.0
              endif
c
 
              if(l_2 .le. nlev) then
                if(lsolar) then
c
                  trndir = trndir1_i(lpd,n,1,nz) + 
     -                     dtrndir1_i(lpd,n,nz)*d_wn
c
                  pd_trndir(lpd,n) =
     -               real((trndir - trn_dir(lpd))*dxi(lpd,n))
c
                  if(ipd(n) .lt. 0) then
c
c*****                  interpolate perturbed layer dependent solar  
c                       source terms to wn
c     
                    dns_src1 = dn_s_src1(lpd,n,1,nz) + 
     -                          ddn_s_src1(lpd,n,nz)*d_wn
                    ups_src1 = up_s_src1(lpd,n,1,nz) + 
     -                          dup_s_src1(lpd,n,nz)*d_wn
c
c*****              find jacobians for solar source functions
c
                    pd_dns_src(lpd,n) = real(dns_src1 - 
     -                                       dns_src(lpd))*dxi(lpd,n)
                    pd_ups_src(lpd,n) = real(ups_src1 - 
     -                                       ups_src(lpd))*dxi(lpd,n)
c
                  else
c                  
                    pd_dns_src(lpd,n) = 0.0
                    pd_ups_src(lpd,n) = 0.0
c
                  endif
c
                endif
c
                if(lplanck) then
c
                  if(ipd(n) .lt. 0) then
c
c*****                  interpolate perturbed layer dependent solar  
c                       source terms to wn
c     
                    dnt_src1 = real(dn_t_src1(lpd,n,1) + 
     -                         ddn_t_src1(lpd,n)*d_wn)
                    upt_src1 = real(up_t_src1(lpd,n,1) + 
     -                         dup_t_src1(lpd,n)*d_wn)
c
c*****              find jacobians for thermal source functions
c
                    pd_dnt_src(lpd,n) = real(dnt_src1 - 
     -                                       dnt_src(lpd))*dxi(lpd,n)
                    pd_upt_src(lpd,n) = real(upt_src1 - 
     -                                       upt_src(lpd))*dxi(lpd,n)
c
                  else
c                  
                    pd_dnt_src(lpd,n) = 0.0
                    pd_upt_src(lpd,n) = 0.0
c
                  endif
c
                endif
c
              else
c
                if(lsolar) then
c
                  pd_dns_src(lpd,n) = 0.0
                  pd_ups_src(lpd,n) = 0.0
                  pd_trndir(lpd,n) = 0.0
c
                endif
                if(lplanck) then
c
                  pd_dns_src(lpd,n) = 0.0
                  pd_ups_src(lpd,n) = 0.0
c
                endif
c
              endif
c
c****          interpolate perturbed radiances to wavenumber, wn,
c              and find partial derivatives
c
              if(ipd(n) .gt. 0) then
                do 4061 k=1,nlout
                    do 4041 naz=1,nphi
                        do 4021 nze=1,numu
c
c*****                         interpolate solar radiance to wn
c
                            rad_s1 = 0.0
                            if(lsolar) rad_s1 = real(solflx*
     -                        (sol_rad1(nze,naz,k,lpd,n,1,nz) +
     -                         dsol_rad1(nze,naz,k,lpd,n,nz)*dwn))
c
c*****                       interpolate thermal radiance to wn
c
                            rad_th1 = 0.0
                            if(lplanck) rad_th1 = 
     -                             real(th_rad1(nze,naz,k,lpd,n,1) +
     -                                 dth_rad1(nze,naz,k,lpd,n)*dwn)
c
c****                         define the total (solar+thermal) radiance
c
                            rad1(nze,naz) = rad_s1 + rad_th1
c
c****                         find radiance jacobian
c
                            pd_rad(nze,naz,k,lpd,n) = (rad1(nze,naz) - 
     -                                        rad(nze,naz,k))*dxi(lpd,n)
        
4021                  continue
4041                continue
4061            continue
c
              endif
c
4981      continue
c
c****      exit loop over output levels
c
5001  continue
c
c****   exit loop over constituent
c
      return
      end
