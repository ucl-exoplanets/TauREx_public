      subroutine init_interp(nlev,nlout,numu,nphi,nz0,wnmin,wnout,
     -                   upsflx,dnsflx,dirsflx,sol_rad,
     -                   uptflx,dntflx,th_rad,
     -                   up_s_src_i,dn_s_src_i,dup_s_src,ddn_s_src,
     -                   up_t_src_i,dn_t_src_i,dup_t_src,ddn_t_src,
     -                   upsflx_i,dnsflx_i,dirsflx_i,dupsflx,
     -                   ddnsflx,ddirsflx,sol_rad0,dsol_rad,
     -                   uptflx_i,dntflx_i,duptflx,ddntflx,
     -                   th_rad0,dth_rad,ups,dns,dirs,
     -                   upth,dnth,rad_s,rad_th,pray0,pgas0,paer0,
     -                   pray_0,pgas_0,paer_0,tau_ray,tau_gas,tau_aer,
     -                   dpray,dpgas,dpaer,tau_ray0,tau_gas0,tau_aer0,
     -                   pray_00,pgas_00,paer_00,dtau_ray,
     -                   dtau_gas,dtau_aer,dpray0,dpgas0,dpaer0,
     -                   pgas,pray,paer)
c
cccccccccccccccccccccccccccc  init_interp  ccccccccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine initializes extinction flags and counters and   cc
cc    variables that record the wavelength range for each wavelength  cc
cc    dependent input file.                                           cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    ngases : number of absorbing gases included in calculation.     cc
cc    nmodes : number of aerosol particle modes used in calculation.  cc
cc                                                                    cc
cc                                                                    cc
cccccccccccccccccccccccccccc  init_interp  ccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'
c
      integer nlev,nlout,numu,nphi,nz0
      integer nz,ne,naz,nze,ni,k
c
      double precision wnmin,wnout(2,nex,nsol)
c
      real upsflx(mxrad),dnsflx(mxrad),dirsflx(mxrad),
     -     sol_rad(mxumu,mxphi,mxrad),
     -     uptflx(mxrad),dntflx(mxrad),
     -     th_rad(mxumu,mxphi,mxrad)
c
c****     monochormatic fluxes at wavenumber wnout
c 
      real upsflx_i(mxrad,2,nsol),dnsflx_i(mxrad,2,nsol),
     -     dirsflx_i(mxrad,2,nsol),uptflx_i(mxrad,2),dntflx_i(mxrad,2)
c
c****   downward and upward solar fluxes gradients at wavenumber wnout
c
      real dupsflx(mxrad,nsol),ddnsflx(mxrad,nsol),ddirsflx(mxrad,nsol),
     -     duptflx(mxrad),ddntflx(mxrad)
c
c****   monochormatic radiances and wavelength gradients at wnout 
c       at output levels, k_out
c
      real sol_rad0(mxumu,mxphi,mxlout,2,nsol),
     -     th_rad0(mxumu,mxphi,mxlout,2)
      real dsol_rad(mxumu,mxphi,mxlout,nsol),dth_rad(mxumu,mxphi,mxrad)
c
c****     monochormatic flux source interpolation variables
c 
      double precision up_s_src_i(mxrad,2,nsol),dn_s_src_i(mxrad,2,nsol)
      double precision up_t_src_i(mxrad,2),dn_t_src_i(mxrad,2)
c
c****   wavelength gradients in source terms
c
      double precision dup_s_src(mxrad,nsol),ddn_s_src(mxrad,nsol)
      double precision dup_t_src(mxrad),ddn_t_src(mxrad)
c
      real ups(mxlout),dns(mxlout),dirs(mxlout),
     -        upth(mxlout),dnth(mxlout)
      real rad_s(mxumu,mxphi,mxlout),rad_th(mxumu,mxphi,mxlout)
c
c****    variables for the output wavenumber grid
c
      real pray_0,pgas_0,paer_0,tau_ray,tau_gas,tau_aer,
     -     pray(mxumu),pgas(mxumu),paer(mxumu)
c
      real dpray0,dpgas0,dpaer0
c
      real pray0(mxumu,2),pgas0(mxumu,2),paer0(mxumu,2),
     -     dpray(mxumu),dpgas(mxumu),dpaer(mxumu),
     -     tau_ray0(2),tau_gas0(2),tau_aer0(2),
     -     pray_00(2),pgas_00(2),paer_00(2)
c
      real dtau_ray,dtau_gas,dtau_aer
c
      nz = nz0
c
c****   initialize spectral variables
c
      do 1001 ne=1,2
           wnout(1,ne,nz) = -999.0d0
           wnout(2,ne,nz) = wnmin
1001  continue
c
c****     initialize the thermal fluxes and radiances
c
      do 1061 k=1,nlev
          upsflx(k) = 0.0
          dnsflx(k) = 0.0
          dirsflx(k) = 0.0
          uptflx(k) = 0.0
          dntflx(k) = 0.0
          do 1041 naz=1,nphi
              do 1021 nze=1,numu
                  sol_rad(nze,naz,k) = 0.0
                  th_rad(nze,naz,k) = 0.0
1021          continue
1041      continue
1061  continue
      do 1161 k=1,nlout
          do 1141 naz=1,nphi
              do 1121 nze=1,numu
                  rad_s(nze,naz,k) = 0.0
                  rad_th(nze,naz,k) = 0.0
1121          continue
1141      continue
1161  continue
c
c****    initialize the optical depths and pressure of tau=1
c
      if(nz .eq. 1) then
        tau_gas = 0.0
        tau_ray = 0.0
        tau_aer = 0.0
        pgas_0 = 0.0
        pray_0 = 0.0
        paer_0 = 0.0
        dtau_gas = 0.0
        dtau_ray = 0.0
        dtau_aer = 0.0
        dpgas0 = 0.0
        dpray0 = 0.0
        dpaer0 = 0.0
c
        do 1201 nze=1,numu
            pgas(nze) = 0.0
            pray(nze) = 0.0
            paer(nze) = 0.0
            dpgas(nze) = 0.0
            dpray(nze) = 0.0
            dpaer(nze) = 0.0
1201    continue
      endif
c
      do 1221 k=1,nlev
          dup_s_src(k,nz) = 0.0d0
          ddn_s_src(k,nz) = 0.0d0
          dupsflx(k,nz) = 0.0
          ddnsflx(k,nz) = 0.0
          ddirsflx(k,nz) = 0.0
1221  continue
c
      do 1361 k=1,nlout
          do 1341 naz=1,nphi
              do 1321 nze=1,numu
                  dsol_rad(nze,naz,k,nz) = 0.0
                  dth_rad(nze,naz,k) = 0.0
1321          continue
1341      continue
1361  continue
c
      do 1801 ni=1,2
          tau_gas0(ni) = 0.0
          tau_ray0(ni) = 0.0
          tau_aer0(ni) = 0.0
          pray_00(ni) = 0.0
          pgas_00(ni) = 0.0
          paer_00(ni) = 0.0
          do 1441 nze=1,numu
              pgas0(nze,ni) = 0.0
              pray0(nze,ni) = 0.0
              paer0(nze,ni) = 0.0
              dpgas(nze) = 0.0
              dpray(nze) = 0.0
              dpaer(nze) = 0.0
1441      continue
c
c****       initialize the solar fluxes and radiances for this solar
c           zenith angle
c
          do 1541 k=1,nlev
              up_s_src_i(k,ni,nz) = 0.0d0
              dn_s_src_i(k,ni,nz) = 0.0d0
              upsflx_i(k,ni,nz) = 0.0
              dnsflx_i(k,ni,nz) = 0.0
              dirsflx_i(k,ni,nz) = 0.0
1541      continue
c
c****     initialize the spectral interpolation variables
c
          do 1641 k=1,nlev
              uptflx_i(k,ni) = 0.0
              dntflx_i(k,ni) = 0.0
1641      continue
c 
          do 1741 k=1,nlout
              ups(k) = 0.0
              dns(k) = 0.0
              dirs(k) = 0.0
              upth(k) = 0.0
              dnth(k) = 0.0
              do 1721 naz=1,nphi
                  do 1701 nze=1,numu
                      sol_rad0(nze,naz,k,ni,nz) = 0.0
                      th_rad0(nze,naz,k,ni) = 0.0
1701              continue
1721          continue
1741      continue
c
1801  continue
c
      return
      end
