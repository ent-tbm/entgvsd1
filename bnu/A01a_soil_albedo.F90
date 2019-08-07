! Carrer_soilalbedo_to_GISS.f
! @auth 10/14/2018 N.Y.Kiang
      
!     Takes Carrer soil albedo annual (min,mean,max, and std) already interpolated to desired spatial resolution and generates GISS GCM boundary conditions input versions.
!     Versions:
!     1) Spectral albedo TOA, GISS 6 bands.  This is the Carrer albedos from MODIS VIS (300-700 nm) and NIR (700-5000 nm) translated to GISS 6 bands, utilizing Judith Lean's top-of-the-atmosphere 1850 preindustrial solar irradiance spectrum(S0=1360.6718 W/m2).
!     2) Spectral albedo Surf, GISS 6 bands.  Like TOA, but using a surface irradiance spectrum from
!     a) NREL's Direct+circumsolar irradiance
!     b) Brian Cairns (GISS) annual averages for: i) 60S-60N, ii) arctic zonal, iii) mid-latitude zonal, and iv) tropical zonal.

!###  To do:  Need to pull out the spectral code from hntr4_Carrer_nc.f
!###  To do:  Doublecheck MODIS NIR 700-4000 nm vs. Carrer 700-5000 nm??      
!###  To do:  May want to correct spectral waiting for 300-400 nm to be lower.
!###  To do:  Get Brian Cairns zonal spectral irradiances.
      
! To compile With netcdf on Mac:
!     sudo port select --set gcc mp-gcc49  
!     gfortran -cpp -fconvert=big-endian -O0 -fno-range-check -I$NETCDFHOME/include -c arrayutil.f convertnc_util.f EntGVSD_util.f Carrer_soilalbedo_to_GISS.f
!     gfortran arrayutil.o convertnc_util.o EntGVSD_util.o Carrer_soilalbedo_to_GISS.o -L$NETCDFHOME/lib -lnetcdf -lnetcdff -o a.out
      
!* To avoid memory problems:
!     > ulimit -s 64000


!------------------------------------------------------------------------------
      
program Carrer_soilalbedo_to_GISS
    use convertnc
    use netcdf
    use EntGVSD_netcdf_util
    
    implicit none
    !!include 'netcdf.inc'

    
!    integer, parameter :: N_BANDS = 6
      
    integer, parameter :: IMK = 7200 !long at 0.05 degrees ~ 1 km
    integer, parameter :: JMK = 3600 !lat at 0.05 degrees ~ 1 im
    integer, parameter :: IMH = 720 !long at 0.5 degrees
    integer, parameter :: JMH = 360 !lat at 0.5 degrees
    integer, parameter :: IM1 = 360 !long at 1 degrees
    integer, parameter :: JM1 = 180 !lat at 1 degrees
    integer, parameter :: IM2 = 144 !long at 2.5 degrees
    integer, parameter :: JM2 = 90 !lat at 2 degrees
    integer, parameter :: IM4X5 = 72 !long at 5 degrees
    integer, parameter :: JM4X5 = 46 !lat at 4 degrees

    !integer, parameter :: NSPEC = 3785 !Number of rows in Lean spectrum
    !character*118 :: GISSBANDStxt = !len = 6*19 + 4
    !Below parameters are moved to EntGVSD_util.f
!     integer, parameter :: N_BANDS = 6 !GISS 6 spectral bands
!      character*101, parameter :: GISSBANDStxt =  &
!          'VIS (330-770), NIR1 (770-860), NIR2 (860-1250), '// &
!          'NIR3 (1250-1500), NIR4 (1500-2200), NIR5 (2200-4000)'

    integer, parameter :: IM = IMH 
    integer, parameter :: JM = JMH 
    character*3,parameter :: res = 'HXH'
    
    real*4 :: lon(IM),lat(JM)

    !MODIS/Carrer bands
    integer, parameter :: MODISVIS = 1
    integer, parameter :: MODISNIR = 2

    !GISS/MODIS band sections
    integer, parameter :: nm300_400 = 1
    integer, parameter :: nm400_700 = 2
    integer, parameter :: nm700_770 = 3
    integer, parameter :: nm770_860 = 4
    integer, parameter :: nm860_1250 = 5
    integer, parameter :: nm1250_1500 = 6
    integer, parameter :: nm1500_2200 = 7
    integer, parameter :: nm2200_4000 = 8

    !GISS ModelE grey albedo fractions
    integer, parameter :: BRIGHT = 1
    integer, parameter :: DARK = 2
    
    ! netcdf names in Carrer soil albedo files
    character*4 :: varnames(4) = (/ &
        'max ','mean','min ','std ' /)
    character*16 :: vark
    integer, parameter :: MMAX=1
    integer, parameter :: MMEAN=2
    integer, parameter :: MMIN=3
    integer, parameter :: MSTD=4
    
    !real*4 :: solargiss(N_BANDS) !W m-2 per GISS 6 bands
    !real*4 :: solarmodis(2)   !W m-2 per MODIS VIS and NIR band
    !real*4 :: solarsplit(5)   !W m-2 per split MODIS/GISS bands
    !real*8 :: solarr8(5)   !W m-2 per split MODIS/GISS bands

    !From Judith Lean 0 km solar surface irradiance (from 2006 version)
    !Fraction of shortwave 300-4000 nm in MODIS & GISS band portions.
    real*4 :: fracSW_MG(8) =  (/ &
        0.105992883, &        !1 300-400 nm
        0.992030205,&	!2 400-700 nm
        0.174912449,&	!3 700-770 nm
        0.203020351,&	!4 770-860 nm
        0.477761295,&	!5 860-1250 nm
        0.089450965,&	!6 1250-1500 nm
        0.17866227,&  !7 1500-2200 nm
        0.051105118 & !8 2200-4000 nm
      /)

    real*4 :: lc_ice(IM,JM)   !Permanent ice cover to be masked out
    real*4 :: lc_water(IM,JM) !Water cover to be masked out
    
    real*4 :: fracSW_MODIS(2) !fraction of 1-VIS and 2-NIR in SW 300-4000 nm
    real*4 :: fracSW_GISS(N_BANDS) !fraction of SW W m-2 nm-1 300-4000 nm

    real*4 :: albmodis(IM,JM,2,4) !1-VIS, 2-NIR; 1-4 MMAX,MMEAN,MMIN,MSTD
    real*4 :: albSW(IM,JM)
    real*4 :: albgiss(IM, JM, N_BANDS)
    real*4 :: fracbd(IM,JM,N_BANDS,2) !1-bright,2-dark
    real*4 :: fracgrey(IM,JM,2) !1-bright,2-dark
    integer :: band(N_BANDS)

    real*4 :: varin(IM,JM)
    
    character*80 :: TITLE
    character*256 :: filein, fileout !, fileoutij - goes to fort.999
    character*80 :: DIR
    character*256 ::PATHin, PATHout
    character*10 :: version = 'v1.0b'  !* UPDATE ME *!
    
    integer :: i, j, k, m, b, dimlon, dimlat, dim(2), dimband
    integer :: ii,jj,dg
    real*4 :: s
    integer :: ncidin, ncidout, status, varid
    integer :: I0i,I1i, J0i,J1i
    integer :: I0o,I1o, J0o,J1o

    
    !Set up if need to do parallel or tiles for large files
    I0i = 1
    I1i = IM
    J0i = 1
    J1i = JM
    I0o = 1
    I1o = IM
    J0o = 1
    J1o = JM

    !* Select spectral band irradiance fractions
    fracSW_GISS(:) = fracSW_Lean_0km_2006(:)
    
    !* Paths.   !* UPDATE ME *!
    DIR = '/Users/nkiang/NancyResearch/GISS/Models/Ent/Datasets/'
    PATHin = trim(DIR)//'Soil/Carrer/'// &
         'Carrer2014_means_2004_Qingsong/'// &
        'Montes_Carrer_temporal_avg_NKfix/'  !Fixed flipped latitude
    PATHout ='/Users/nkiang/NancyResearch/GISS/Models/Ent/Datasets/'// &
        'Soil/Carrer/CarrerGISS_soil_albedo/'

    !* Make Output file - netcdf - ##UPDATE ME##
    fileout = trim(PATHout)// &
        'CarrerGISS_soil_albedo_multiband_annual_2004_'// &
        trim(version)//"_"//trim(res)//'.nc'
    call my_nf_create_ij(trim(fileout),IM,JM, &
        ncidout,dimlon,dimlat)
    write(0,*) "dimlon, dimlat", dimlon, dimlat
    call my_nf_defglobal(trim(fileout), &
        'Soil albedo derived from D. Carrer et al. (2014) '// &
        'MODIS soil albedos.  Created at same as v1.0b Ent GVSD.', &
        'Includes copy of source VIS and NIR annual averages. '// &
        'Partial coastline and permanent ice cells filled by '// &
        'nearest-5 neighbor. fringeice-uses cell albedo if no '// &
        'neighbor; shows islands.  noice-undef if no neighbor. '// &
        ' Fortran program:  Carrer_soilalbedo_to_GISS.f')
    status = nf_close(ncidout)
    status = nf_open(trim(fileout),NF_WRITE,ncidout)
    write(*,*) status, 'nf_open out ',trim(fileout)
    !* Add some more def attributes
    status = nf_redef(ncidout)
    status=nf_def_dim(ncidout, 'band', N_BANDS, dimband)
    write(0,*) 'GISSBANDStxt:', GISSBANDStxt
    status = my_nf_inq_put_att_any(ncidout, &
        'GISS_bands', 'description',GISSBANDstxt)
    write(0,*) 'def band status',status, dimband
    status = nf_enddef(ncidout)

    !* Input file - Permanent ice and water -----------------------------
    filein = 'lc_lai_ent/EntMM_lc_max_'//trim(res)//'.'// &
        trim(version)//'.nc'
    print *, trim(filein)
    status = nf_open(trim(filein),0,ncidin)
    status = my_nf_inq_get_var_real32_2(ncidin,'permanent_ice' &
           ,varid,varin)
    lc_ice(:,:) = varin(:,:)
    status = my_nf_inq_get_var_real32_2(ncidin,'water' &
           ,varid,varin)
    lc_water(:,:) = varin(:,:)
    status = nf_close(ncidin)
    
    !* Input file - VIS -------------------------------------------------
    !filein = trim(PATHin)//'VIS_Alb_soil_yearly.006.2004.NK.nc' !## UPDATE ME
    filein = trim(PATHin)// &
        'Carrer_soil_albedo_VIS_NIR_annual_2004_upscale_HXH.nc' !## UPDATE ME
    write(*,*) 'filein: ', filein
    status = nf_open(trim(filein),0,ncidin)
    write(*,*) status, 'nf_open in ',trim(filein)

     !* Get VIS
    do k=1,4
       vark = 'soilalb_VIS_'//trim(varnames(k))
       write(*,*) 'VIS ', vark
       status = my_nf_inq_get_var_real32_2(ncidin,trim(vark) &
           ,varid,varin)
       write(*,*) 'Got here 1', status
       albmodis(:,:,MODISVIS, k) = varin(:,:)
       write(*,*) 'Got here'
       status = my_nf_inq_def_put_var_real32_2(ncidout, &
             IM,JM,1,IM,1,JM,dimlon,dimlat, &
!             "soilalb_VIS_"//trim(vark), &
             trim(vark), &
             "Carrer soil albedo VIS (300-700 nm) annual " &
             //trim(vark), &
             "fraction", varin)
       !status = nf_close(ncidin)
       write(*,*) status, 'nf_close in ',trim(filein)
    enddo

    !* Input file - NIR --------------------------------------------
    !filein = trim(PATHin)//'NIR_Alb_soil_yearly.006.2004.NK.nc' !## UPDATE ME
!   filein = trim(PATHin)//!'NIR_Alb_soil_yearly.006.2004.NK.nc' !## UPDATE ME &
!        'Carrer_soil_albedo_VIS_NIR_annual_2004_upscale_HXH.nc' !## UPDATE ME
!     write(*,*) 'filein: ', filein
!     status = nf_open(trim(filein),0,ncidin)
!     write(*,*) status, 'nf_open in ',trim(filein)
     !write(0,*) trim(nf90_strerror(status)), trim(filein)

    !* Get NIR
    do k=1,4
       vark = 'soilalb_NIR_'//trim(varnames(k))
       status = my_nf_inq_get_var_real32_2(ncidin,trim(vark) &
           ,varid,varin)
       albmodis(:,:,MODISNIR,k) = varin
       
       !* Before writing, check against land mask for any undef on land.
       status = my_nf_inq_def_put_var_real32_2(ncidout, &
             IM,JM,1,IM,1,JM,dimlon,dimlat, &
!             "soilalb_NIR_"//trim(vark), &
             trim(vark), &
             "Carrer soilalbedo NIR (700-5000 nm) annual " &
             //trim(vark), &
             "fraction", varin)
      end do
      status = nf_close(ncidin)
      write(*,*) status, 'nf_close in ',trim(filein)
      
      !* Calculate total SW albedo --------------------------------------
      albSW(:,:) = FillValue
      do i=I0o,I1o
         do j=J0o,J1o
            if ((albmodis(i,j,MODISVIS,MMEAN).eq.FillValue).or. &
                (albmodis(i,j,MODISNIR,MMEAN).eq.FillValue)) then
               albSW(i,j) = FillValue
            else
               albSW(i,j) =  &
                   ( albmodis(i,j,MODISVIS,MMEAN)*sum(fracSW_MG(1:2)) !( fracSW_MG(nm300_400) + fracSW_MG(nm400_700)) + &
                   + albmodis(i,j, MODISNIR,MMEAN)* &
                   sum(fracSW_MG(3:8)) )/sum(fracSW_MG(1:8)) 
            endif
         enddo
      enddo
      write(*,*) 'Got here SW'
      status = my_nf_inq_def_put_var_real32_2(ncidout, &
          IM,JM,1,IM,1,JM,dimlon,dimlat, &
          'albedo_soil_SW', &
          'soil albedo shortwave', &
          'fraction', albSW)
       status = my_nf_inq_put_att_any(ncidout, &
          'albedo_soil_SW', 'description', &
          'Carrer annual mean soil albedo shortwave 400-4000 nm')

      
      !* Put spectral breakdown into GISS bands ------------------------
      albgiss(:,:,:) = FillValue
      do i=I0o,I1o
         do j=J0o,J1o
            if ((albmodis(i,j,MODISVIS,MMEAN).eq.FillValue).or. &
                (albmodis(i,j,MODISNIR,MMEAN).eq.FillValue)) then
               albgiss(i,j,1:N_BANDS) = FillValue
            else
               albgiss(i,j,GISS300_770NM) = &
                   ( albmodis(i,j,MODISVIS,MMEAN) * &
                   (fracSW_MG(nm300_400) + fracSW_MG(nm400_700)) + &
                   albmodis(i,j,MODISNIR,MMEAN)*fracSW_MG(nm700_770)) &
                   /( fracSW_MG(nm300_400) + fracSW_MG(nm400_700) + &
                   fracSW_MG(nm700_770) )
               do k=2,N_BANDS
                  albgiss(i,j,k) = albmodis(i,j,MODISNIR,MMEAN) !All NIR bands
               end do
            endif
         enddo
      enddo
      status = my_nf_inq_def_put_var_real32_3t(ncidout, &
          IM,JM,N_BANDS, &
          1,IM,1,JM,dimlon,dimlat,dimband, &
          'albedo_soil_giss', &
          'soil albedo GISS bands', 'fraction', albgiss)
      status = my_nf_inq_put_att_any(ncidout, &
          'albedo_soil_giss', 'description', &
          'Carrer annual mean soil albedo, GISS GCM spectral bands')

!      !------ Old VEG bare_bright vs. bare_dark soil fractions
!      !* Calculate GISS ModelE grey albedo bright vs. dark fractions
!      !Same netcdf names as in Ent_pfts.
!      !1) Non-grey albedos: each GISS band has its own bright|dark fractions.
!      !  For partial-land cells (that have ice and/or water) use
       ! soil albedo of nearest all-land neighbor.  If none, default to
       ! to cell's mixed albedo.
       fracbd(:,:,:,:) = FillValue
       do i=I0o,I1o
         do j=J0o,J1o
            do k=1,N_BANDS
               if (albgiss(i,j,k).ne.FillValue) then
                  if ((lc_ice(i,j)+lc_water(i,j)).eq.1.0) then
                     !All permanent ice, no soil.
                     fracbd(i,j,k,BRIGHT) = FillValue
                     fracbd(i,j,k,DARK) = FillValue
                  elseif (((lc_ice(i,j)+lc_water(i,j)).gt.0.).and. &
                        ((lc_ice(i,j)+lc_water(i,j)).lt.1.0)) then  !Partial
                     !Search for nearest grid cell up to 5 that is all ground
                     s = FillValue
                     do dg=1,5
                      do jj=j-dg,j+dg
                       do ii=i-dg,i+dg
                        if ((lc_ice(ii,jj)+lc_water(ii,jj)).eq.0.) then
                         !All ground found
                           s = albgiss(ii,jj,k)
                           exit
                        endif
                     enddo
                     if (s.ne.FillValue) exit
                    enddo
                    if (s.ne.FillValue) exit
                   enddo
                   if (s.eq.FillValue) then !No nearby ground cells found 
                      !fracbd(i,j,k,BRIGHT) = FillValue
                      !fracbd(i,j,k,DARK) = FillValue
!                     !*noice - If cannot assign land albedo, cut the fractional cells.
!                     !  Loses Alaska peninsula
!                     !fracgrey(i,j,:) = FillValue
!                     !*fringeice - Or, since there is some ground, assign i,j value there.
!                     !  Leaves some coastal ice fringe.
                      fracgrey(i,j,BRIGHT) = &
                          min(0.5, albgiss(i,j,k))/0.5 !
                      fracgrey(i,j,DARK) = 1.0 - fracgrey(i,j,BRIGHT)
                   else
                      fracbd(i,j,k,BRIGHT) = min(s,0.5)/0.5
                      fracbd(i,j,k,DARK) = 1.0 - fracbd(i,j,k,BRIGHT)
                   endif
                else
                  !All ground
                   fracbd(i,j,k,BRIGHT) = min(albgiss(i,j,k),0.5)/0.5
                   fracbd(i,j,k,DARK) = 1.0 - fracbd(i,j,k,BRIGHT)
                endif
             else
                fracbd(i,j,k,BRIGHT) = FillValue
                fracbd(i,j,k,DARK) = FillValue
             endif
          enddo
        enddo
       enddo
      
      status = my_nf_inq_def_put_var_real32_3t(ncidout, &
          IM,JM,N_BANDS, &
          1,IM,1,JM,dimlon,dimlat,dimband, &
          'bare_bright_band', &
          'soil albedo bright fraction in GISS spectral bands', &
          'fraction', fracbd(:,:,:,BRIGHT))
      status = my_nf_inq_put_att_any(ncidout, &
          'bare_bright_band', 'description', &
          'Annual mean soil albedo band bright fraction')
      
      status = my_nf_inq_def_put_var_real32_3t(ncidout, &
          IM,JM,N_BANDS, &
          1,IM,1,JM,dimlon,dimlat,dimband, &
          'bare_dark_band', &
          'soil albedo dark fraction in GISS spectral bands', &
          'fraction', fracbd(:,:,:,DARK))
      status = my_nf_inq_put_att_any(ncidout, &
          'bare_dark_band', 'description', &
          'Annual mean soil albedo band dark fraction')

!     !2) Grey shortwave albedos - exclude 100% ice+water, nearest-neighbor fill
!      !  For partial-land cells (that have water) use
       ! soil albedo of nearest all-land/non-ice neighbor.  If none, default to
       ! to cell's mixed albedo.
      fracgrey(:,:,:) = 0.0
      do i=I0o,I1o
         do j=J0o,J1o
            if (albSW(i,j).ne.FillValue) then 
               if ((lc_ice(i,j)+lc_water(i,j)).eq.1.0) then !No ground
                  fracgrey(i,j,:) = FillValue
               elseif (((lc_ice(i,j)+lc_water(i,j)).gt.0.).and. &
                      ((lc_ice(i,j)+lc_water(i,j)).lt.1.0)) then !Partial
                  !Search for nearest grid cell up to 5 that is all ground
                  s = FillValue
                  do dg=1,5
                  do ii=i-dg,i+dg
                     do jj=j-dg,j+dg
                        if ((lc_ice(ii,jj)+lc_water(ii,jj)).eq.0.) then
                          !All ground found
                          s = albSW(ii,jj)
                          exit
                       endif
                    enddo
                    if (s.ne.FillValue) exit
                  enddo
                  if (s.ne.FillValue) exit
                  enddo
                  if (s.eq.FillValue) then !No nearby ground cells found
!                    !*If cannot assign land albedo, cut the fractional cells.
!                    !Loses Alaska peninsula
                     !fracgrey(i,j,:) = FillValue
!                    !*Or, since there is some ground, assign i,j value there.
!                    !Leaves some coastal ice fringe.
                     fracgrey(i,j,BRIGHT) = min(0.5, albSW(i,j))/0.5 !
                     fracgrey(i,j,DARK) = 1.0 - fracgrey(i,j,BRIGHT)
                  else !Found nearby all-ground cell
                     fracgrey(i,j,BRIGHT) = min(0.5,s)/0.5
                     fracgrey(i,j,DARK) = 1.0 - fracgrey(i,j,BRIGHT)
                  endif
               else    !Cell is all ground with SW
                  fracgrey(i,j,BRIGHT) = min(0.5, albSW(i,j))/0.5
                  fracgrey(i,j,DARK) = 1.0 - fracgrey(i,j,BRIGHT)
               endif
            else
               fracgrey(i,j,:) = FillValue
            endif
         enddo
      enddo
      status = my_nf_inq_def_put_var_real32_2(ncidout, &
             IM,JM,I0o,I1o,J0o,J1o, &
          dimlon,dimlat, &
          'bare_bright_grey','soil albedo grey SW bright fraction', &
          'fraction', fracgrey(:,:,BRIGHT))
      status = my_nf_inq_put_att_any(ncidout, &
          'bare_bright_grey', 'description', &
          'Carrer annual mean SW soil albedo grey bright fraction')
      status = my_nf_inq_def_put_var_real32_2(ncidout, &
          IM,JM,I0o,I1o,J0o,J1o, &
          dimlon,dimlat, &
          'bare_dark_grey', 'soil albedo grey SW dark fraction', &
          'fraction', fracgrey(:,:,DARK))
      status = my_nf_inq_put_att_any(ncidout, &
          'bare_dark_grey', 'description',&
          'Carrer annual mean SW soil albedo grey dark fraction')
      !status = nf_close(ncidin)
      !write(*,*) status, 'nf_close in ',trim(filein)

      
      status = nf_close(ncidout)
      write(*,*) status, 'nf_close out ',trim(fileout)
         
      end program Carrer_soilalbedo_to_GISS
            
