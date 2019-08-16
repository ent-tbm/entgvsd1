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


    use carrer_mod
    use netcdf
    use chunker_mod
    use chunkparams_mod
    use paths_mod
    use ent_labels_mod
    use geom_mod
    use assign_laimax_mod


!    use convertnc
!    use netcdf
!    use EntGVSD_netcdf_util
    
    implicit none

    
!    integer, parameter :: N_BANDS = 6

    !integer, parameter :: NSPEC = 3785 !Number of rows in Lean spectrum
    !character*118 :: GISSBANDStxt = !len = 6*19 + 4
    !Below parameters are moved to EntGVSD_util.f
!     integer, parameter :: N_BANDS = 6 !GISS 6 spectral bands
!      character*101, parameter :: GISSBANDStxt =  &
!          'VIS (330-770), NIR1 (770-860), NIR2 (860-1250), '// &
!          'NIR3 (1250-1500), NIR4 (1500-2200), NIR5 (2200-4000)'

    integer, parameter :: NBANDS_GISS = 6
    character*4, parameter :: sbands_giss(NBANDS_GISS) = &
        (/"VIS ", "NIR1", "NIR2", "NIR3", "NIR4", "NIR5"/)
    real*8, parameter :: bandpoints_giss(NBANDS_GISS+1) = &
        (/300,770,860,1250,1500,2200,4000/)





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


Convert spectrum in one set of bands to another set of bands while conserving energy
    !From Judith Lean 0 km solar surface irradiance (from 2006 version)
    !Fraction of shortwave 300-4000 nm in MODIS & GISS band portions.
    real*4, parameter :: fracSW_MG(8) =  (/ &
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

    ! Input files
    type(ChunkIO_t), target :: ioall_albcarr(NBANDS_CARRER)
    type(ChunkIO_t), target :: io_albcarr(NSTATS,NBANDS_CARRER)
    integer :: iband,istat
    real, dimension(:,:), allocatable, target :: wta    ! Surmised landmask

    !* Select spectral band irradiance fractions
    fracSW_GISS(:) = fracSW_Lean_0km_2006(:)

    call chunker%init(IMK, JMK, IMH*2,JMH*2, 'qxq', 10, 10, 10, (/6,5/))
    allocate(wta(chunker%chunk_size(1), chunker%chunk_size(2)))

    ! ------------- Open Input Files
    do iband=1,NBANDS_CARRER
        call chunker%nc_open(ioall_albmodis(iband), LC_LAI_ENT_DIR, &
            'carrer/', &
            'soilalb_'//trim(sbands(iband))//'.nc', &
            'soilalb_'//trim(sbands(iband)), 0)
        do istat=1,nstats
            call chunker%nc_reuse_var( &
                ioall_albmodis(iband), io_albmodis(istat,iband), &
                (/1,1,istat/))
        end do
    end do

    ! ---------------- Open Output Files
    call chunker%nc_create(ioall_albgiss, weighting(wta,1d0,0d0), &
        'carrer/', &
        'albgiss', 'albgiss', &
        'Soil albedo in GISS bands', &
        '1', sbands_giss)
    do iband=1,NBANDS_GISS
        call chunker%nc_reuse_var(ioall_albgiss, io_albgiss(iband), (/1,1,iband/), &
            weighting(wta, 1d0, 0d0))
    end do


#ifdef ENTGVSD_DEBUG
    !do jchunk = chunker%nchunk(2)*3/4,chunker%nchunk(2)*3/4+1
    !do ichunk = chunker%nchunk(1)*3/4,chunker%nchunk(1)*3/4+1
    do jchunk = 1,3
    do ichunk = 2,4
#else
    do jchunk = 1,chunker%nchunk(2)
    do ichunk = 1,chunker%nchunk(1)
#endif

        call chunker%move_to(ichunk,jchunk)
        wta = 0

        do jc = 1,chunker%chunk_size(2)
        do ic = 1,chunker%chunk_size(1)

            ! Compute overall NetCDF index of current cell
            ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
            jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

            ! Infer soil mask
            if (io_albmodis(SMEAN,1)%buf(ic,jc) /= FillValue)
                wta(ic,jc) = 1
            else
                sta(ic,jc) = 0
            end if


            !* Calculate total SW albedo --------------------------------------
            if ((albmodis(SMEAN, VIS_MODIS).eq.FillValue).or. &
                (albmodis(SMEAN, NIR_MODIS).eq.FillValue)) then
               io_albsw%buf(ic,jc) = FillValue
            else
               io_albsw%buf(ic,jc) = ( &
                     io_albmodis(SMEAN,VIS_MODIS)%buf(ic,jc) * sum(fracSW_MG(1:2)) &
                   + io_albmodis(SMEAN,NIR_MODIS)%buf(ic,jc) * sum(fracSW_MG(3:8)) &
                ) / sum(fracSW_MG(1:8)) 
            endif

            !* Put spectral breakdown into GISS bands ------------------------
            do iband=1,NBANDS_GISS
                io_albgiss(iband)%buf = FillValue
            end do

            if ((io_albmodis(SMEAN,VIS_MODIS)%buf(ic,jc).eq.FillValue).or. &
                (io_albmodis(SMEAN,NIR_MODIS)%buf(ic,jc).eq.FillValue)) then
                do iband=1,NBANDS_GISS
                   io_albgiss(iband)%buf(ic,jc) = FillValue
                end do
            else
                io_albgiss(VIS_GISS)%buf(ic,jc) = &
                   ( io_albmodis(SMEAN,VIS_MODIS)%buf(ic,jc) * &
                   (fracSW_MG(nm300_400) + fracSW_MG(nm400_700)) + &
                   io_albmodis(SMEAN,NIR_MODS)%buf(ic,jc)*fracSW_MG(nm700_770)) &
                   /( fracSW_MG(nm300_400) + fracSW_MG(nm400_700) + &
                   fracSW_MG(nm700_770) )
               do k=2,NBANDS_GISS
                    io_albgiss(k)%buf(ic,jc) = io_albmodis(SMEAN, NIR_MODIS)%buf(ic,jc)
               end do
            endif

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
ice+water==1 ==> no soil in that gridcell.  So we must interpolate what soil under that ice/water might be.  SO search for adjacent gridcells that have soil in them.  It can go up to 5, but not as big as 5.  dg=5.  Played around with dg, 5 was biggest had to go.

At 1km.. skip cover that is ocean, permanetn lake or permanent ice
Don't want to skip anything that's just water.

First look at output files and input files from Nancy's lo-res version


First do a map of all gricells where lc_ice and lc_water add up to 1, see where they are.  If they're ocean or Antarctica, it's fine if we exclude them.  Don't want to use MODIS water cover for lots of lakes in Arctic.

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
This is where had to cap bright/dark fraction to .5
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
            
