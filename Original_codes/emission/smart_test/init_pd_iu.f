      subroutine init_pd_iu(nza,numu,nphi,nlout,ifrmout,
     -                      nza_1,radfile,lplanck,lsolar,
     -                      clev,nstate,istate,iu_pd,
     -                      j_ext,iutpd,iuspd,iupdrad)
c
cccccccccccccccccccccccc  i n i t _ p d _ i u   cccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    This subroutine initializes the output units for radiance and   cc
cc    and flux jacobians.                                             cc
cc                                                                    cc
cc    i n p u t:                                                      cc
cc                                                                    cc
cc        nza - number of solar zenith angles                         cc
cc       numu - number of output zenith angles                        cc
cc       nphi - number of output azimuth angles                       cc
cc        npd - number of variable elements of the state vector       cc 
cc      nlout - number of levels where radiance spectra are output    cc
cc     ifmout - index of output file format (1) ascii, (2) binary,    cc
cc                                          (3) binary, no header     cc
cc      nza_1 - index of first stream to be printed out               cc
cc       name - string vector with parsed name of output file         cc
cc        len - length of name (no leading or trailing blanks)        cc
cc     lsolar - include solar fluxes? (logical: T/F)                  cc
cc    lplanck - include thermal fluxes? (logical: T/F)                cc
cc     nstate - number of elements in the state vector                cc
cc     istate - state vector flag indicating which state variables    cc
cc              are variable components of the state vector.          cc
cc              0 - not a variable component of the state vector      cc
cc              1 - surface pressure                                  cc
cc              2 - surface/atmospheric temperature                   cc
cc              3 - gas absorption coeffient                          cc
cc              4 - cloud/aerosol optical depth                       cc
cc              5 - surface albedo                                    cc
cc      iu_pd - starting index of partial derivative output units     cc
cc       clev - character vector with out level index                 cc
cc      j_ext - character varible with state vector type (istate)     cc 
cc                                                                    cc
cc    o u t p u t:                                                    cc
cc                                                                    cc
cc       iutpd - unit number for each thermal flux 
cc                                                                    cc
cccccccccccccccccccccccc  i n i t _ p d _ i u   cccccccccccccccccccccccc
c
      implicit none
c
      include 'param.inc'     
c
      character*20 pd_ext
      character*132 pdfile,radfile(mxlout,nsol)
      character*1 name(132),pd_name(9),nm1(132)
      character*9 j_ext(nex)
      character*11 sflx
      character*4 clev(mxlout)
c
      logical lplanck,lsolar
      
      integer nza,numu,nphi,nlout,nza_1(mxlout),ifrmout,len
c
c****    units for flux and radiance jacobians
c 
      integer nstate,istate(nex),npd,
     -        iu_pd,iutpd(mxpd,nsol),iuspd(mxpd,nsol),
     -        iupdrad(mxumu,mxphi,mxpd,mxlout,nsol)
c
c****   local integer variables
c
      integer n,nz,nl,naz,nze,ipdct,ir,j,lenpd,nlb,ntb,lenf
c
c****       create an output files for each partial derivative
c
      ipdct = -1
c
c****   enter loop over solar zenith angle
c
      do 1281 nz=1,nza
c
          npd = 0
c
c****       enter loop over partial derivative
c
          do 1261 n=1,nstate
c
c****            flux partial derivtives are only calculated 
c                if istate is less than zero
c
              if(istate(n) .lt. 0) then
c
                call charsp(radfile(1,nz),name,len,132,nlb,ntb)
c
                npd = npd + 1
c
                call charsp(j_ext(n),pd_name,lenpd,9,nlb,ntb) 
c
c****             increment the unit counter.  This counter has to
c                 be incremented even if thermal fluxes aren't found
c                 in this run to support the separate solar and thermal
c                 flux calculations in vpl_climate
c
                if( nz .eq. 1) then 
                  ipdct = ipdct + 1
                  if(lplanck) then
c
c****                 open units for partial derivatives for the  
c                     diffuse thermal flux
c
                    write(pdfile,'(132a)') (name(j),j=1,len),'_tflx',
     -                 (pd_name(j),j=1,lenpd)
c
                    iutpd(npd,nz) = iu_pd + ipdct 
                    if(ifrmout .eq. 1) then
                      open(iutpd(npd,nz),file=pdfile,form='formatted',
     -                          status='unknown')
                    else
                      open(iutpd(npd,nz),file=pdfile,form='unformatted',
     -                        status='unknown')
                    endif
c
                    call charsp(pdfile,nm1,lenf,132,nlb,ntb) 
c
                    write(*,'(i5,a16,132a)') iutpd(npd,nz),
     -                 ' Jacobian  ',(nm1(j),j=1,lenf)
c
                  endif
                endif
c
c****             increment the unit counter.  This counter has to
c                 be incremented even if solar fluxes aren't found
c                 in this run to support the separate solar and thermal
c                 flux calculations in vpl_climate
c
                ipdct = ipdct + 1
                if(lsolar) then
c
c****               open unit for partial derivatives for downward 
c                   solar flux
c
                  write(sflx,'(1a9,i2.2)') '_sflx_sza',nz
                  write(pdfile,'(132a)') (name(j),j=1,len),sflx,
     -                 (pd_name(j),j=1,lenpd)
c
                  iuspd(npd,nz) = iu_pd + ipdct 
                  if(ifrmout .eq. 1) then
                    open(iuspd(npd,nz),file=pdfile,form='formatted',
     -                            status='unknown')
                  else
                    open(iuspd(npd,nz),file=pdfile,form='unformatted',
     -                          status='unknown')
                  endif
c
c
                    call charsp(pdfile,nm1,lenf,132,nlb,ntb) 
c
                    write(*,'(i5,a16,132a)') iuspd(npd,nz),
     -                 ' Jacobian  ',(nm1(j),j=1,lenf)
c
                endif
c
              endif
c
              if(istate(n) .gt. 0) then
c
                npd = npd + 1
c
                call charsp(j_ext(n),pd_name,lenpd,9,nlb,ntb) 
c
c****            r a d i a n c e   p a r t i a l   d e r i v a t i v e s
c
c***             enter loop over output level
c
                do 1241 nl=1,nlout
c
                    call charsp(radfile(nl,nz),name,len,132,nlb,ntb)
c
c****                 enter loops over emission azimuth and zenith angle
c
                    ir = 0
                    do 1221 naz=1,nphi
c
                        do 1201 nze=nza_1(nl),numu
c
                            ir = ir + 1
                            ipdct = ipdct + 1
                            write(pd_ext,'(1a4,i3.3,1a4,9a)') 
     -                             '_rad',ir,clev(nl),
     -                             (pd_name(j),j=1,lenpd)
c  
c****                          compute the unit number
c    
                            iupdrad(nze,naz,npd,nl,nz) = iu_pd + ipdct
c
                            write(pdfile,'(132a)') 
     -                           (name(j),j=1,len),pd_ext
                            if(ifrmout .eq. 1) then
                              open(iupdrad(nze,naz,npd,nl,nz),
     -                             file=pdfile,form='formatted',
     -                             status='unknown')
                            else
                              open(iupdrad(nze,naz,npd,nl,nz),
     -                             file=pdfile,form='unformatted',
     -                             status='unknown')
                            endif
c
                            call charsp(pdfile,nm1,lenf,132,nlb,ntb) 
c
                            write(*,'(i5,a16,132a)') 
     -                            iupdrad(nze,naz,npd,nl,nz),
     -                            ' Jacobian  ',(nm1(j),j=1,lenf)
1201                    continue
1221                continue
1241            continue
c
              endif
c
1261      continue          
1281  continue
c
      return
      end
