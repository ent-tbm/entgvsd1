! Author 10/14/2018 N.Y.Kiang, Elizabeth Fischer
!      

!     Takes Carrer soil albedo annual (min,mean,max, and std) already
!     interpolated to desired spatial resolution and generates GISS
!     GCM boundary conditions input versions.  Creates downscaled 1 km
!     bare soil bright fraction ('bs_brightratio.nc') for ModelE grey
!     bare soil albedo scheme.
!
!     Versions:
!
!  1) Spectral albedo TOA, GISS 6 bands.  This is the Carrer albedos
!     from MODIS VIS (300-700 nm) and NIR (700-5000 nm) translated to
!     GISS 6 bands, utilizing Judith Lean's top-of-the-atmosphere 1850
!     preindustrial solar irradiance spectrum(S0=1360.6718 W/m2).
!
!  2) Spectral albedo Surf, GISS 6 bands.  Like TOA, but using a
!  surface irradiance spectrum from
!     a) NREL's Direct+circumsolar irradiance
!     b) Brian Cairns (GISS) annual averages for: i) 60S-60N, ii)
!        arctic zonal, iii) mid-latitude zonal, and iv) tropical zonal.
!
!  3) Bare soil bright fraction 'bs_brightratio.nc' for GISS ModelE
!     grey soil albedo scheme.

#ifdef JUST_DEPENDENCIES
#    define THIS_OUTPUTS_DIR MKFILES_DIR
#else
#    define THIS_OUTPUTS_DIR DEFAULT_OUTPUTS_DIR
#endif

! Use the filled-in version of Carrer MODIS albedo
#define USE_FILLED

module B07_mod

    use carrer_mod
    use netcdf
    use chunker_mod
!    use ent_labels_mod
    use ent_params_mod
!    use assign_laimax_mod
    use hntr_mod

    implicit none

    ! GISS Band Defintions
    integer, parameter :: VIS_GISS = 1
    integer, parameter :: NBANDS_GISS = 6
    character*4, parameter :: sbands_giss(NBANDS_GISS) = &
        (/"VIS ", "NIR1", "NIR2", "NIR3", "NIR4", "NIR5"/)
    character*(*), parameter :: sbands_giss_long(NBANDS_GISS) = &
        (/"VIS (300-770nm)   ", &
          "NIR1 (770-860nm)  ", &
          "NIR2 (860-1250nm) ", &
          "NIR3 (1250-1500nm)", &
          "NIR4 (1500-2200nm)", &
          "NIR5 (2200-4000nm)"/)
    real*8, parameter :: bandpoints_giss(NBANDS_GISS+1) = &
        (/300,770,860,1250,1500,2200,4000/)

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
    integer, parameter :: NBRIGHT_DARK = 2
    character*6, parameter :: sbright_dark(2) = (/"bright", "dark  "/)
    character*11, parameter :: sbright_dark_long(2) = (/"bright soil", "dark soil  "/)

    ! Convert spectrum in one set of bands to another set of bands while conserving energy
    !From Judith Lean 0 km solar surface irradiance (from 2006 version)
    !Fraction of shortwave 300-4000 nm in MODIS & GISS band portions.
    !Calculated in solar.lean_TOA.0km_toNancyin2006.xlsx
    real*4, parameter :: fracSW_MG(8) =  (/ &
         0.046632595,&   !1 300-400 nm
         0.436453295,&   !2 400-700 nm
         0.076954426,&   !3 700-770 nm
         0.089320769,&   !4 770-860 nm
         0.210195708,&   !5 860-1250 nm
         0.039354818,&   !6 1250-1500 nm
         0.078604196,&   !7 1500-2200 nm
         0.022484192 &   !8 2200-4000 nm
      /)

      ! Fraction of shortwave (SW) surface irradiance (300-4000 nmm)
      !  in GISS bands.
      ! From solar.lean_TOA.0km_toNancyin2006.xlsx.  Old 1366 W/m2 S0, but at surface
      real*4, parameter :: fracSW_Lean_0km_2006(NBANDS_GISS) = (/ &
          0.56, 0.0893, 0.2102, 0.0394, 0.0786, 0.0225  /)


    !integer, parameter :: NSPEC = 3785 !Number of rows in Lean spectrum
    !character*118 :: GISSBANDStxt = !len = 6*19 + 4
    !Below parameters are moved to EntGVSD_util.f
!     integer, parameter :: NBANDS_GISS = 6 !GISS 6 spectral bands
!      character*101, parameter :: GISSBANDStxt =  &
!          'VIS (330-770), NIR1 (770-860), NIR2 (860-1250), '// &
!          'NIR3 (1250-1500), NIR4 (1500-2200), NIR5 (2200-4000)'
    
    real*4 :: fracSW_GISS(NBANDS_GISS) !fraction of SW W m-2 nm-1 300-4000 nm


type Output_t
    type(Chunker_t) :: chunker
    type(HntrCalc_t) :: hntr    ! Preparation to regrid
    character(8) :: sres    ! Resolution indicator

    logical :: create_lr    ! Create a _forplot version?  Only do this for original 6km...
    real, dimension(:,:), allocatable :: wta1   ! Mask that's 1 everywhere
    real, dimension(:,:), allocatable :: wta   ! Surmised landmask based on albmodis
    real, dimension(:,:), allocatable :: wta_fracbd   ! Surmised landmask based on fracbd

    type(ChunkIO_t) :: io_albsw
    type(ChunkIO_t) :: io_albgiss(NBANDS_GISS)
    type(ChunkIO_t) :: ioall_fracbd(NBANDS_GISS), io_fracbd(NBRIGHT_DARK,NBANDS_GISS)
    type(ChunkIO_t) :: ioall_fracgrey, io_fracgrey(NBRIGHT_DARK)
    type(ChunkIO_t) :: io_bs_brightratio
end type Output_t


CONTAINS

subroutine B07_init
    !* Select spectral band irradiance fractions
    fracSW_GISS(:) = fracSW_Lean_0km_2006(:)

end subroutine B07_init

subroutine init_output(out, ichunker, sres, IM, JM, create_lr)
    type(output_t), target, intent(INOUT) :: out
    type(Chunker_t), intent(IN) :: ichunker
    character*(*), intent(IN) :: sres
    integer, intent(IN) :: IM, JM
    logical :: create_lr
    ! -------------- Local Vars
    type(HntrSpec_t) :: ospec, ispec

    out%sres = sres
    call out%chunker%init(IM, JM, IMH*2,JMH*2, 'forplot', 100, 100, 10,(/18,15/), &
        outputs_dir=THIS_OUTPUTS_DIR)
    ospec = hntr_spec(out%chunker%chunk_size(1), out%chunker%ngrid(2), 0d0, 180d0*60d0 / out%chunker%ngrid(2))
    ispec = hntr_spec(ichunker%chunk_size(1), ichunker%ngrid(2), 0d0, 180d0*60d0 / ichunker%ngrid(2))
    out%hntr = hntr_calc(ospec, ispec, 1d0*FillValue)    ! datmis=FillValue

    ! Only need weightings for original (6km) resolution
    out%create_lr = create_lr
    if (create_lr) then
        allocate(out%wta1(out%chunker%chunk_size(1), out%chunker%chunk_size(2)))
        allocate(out%wta(out%chunker%chunk_size(1), out%chunker%chunk_size(2)))
        allocate(out%wta_fracbd(out%chunker%chunk_size(1), out%chunker%chunk_size(2)))
        out%wta1 = 1
    end if

end subroutine init_output

subroutine open_output_files(out, rw)
    type(output_t), target, intent(INOUT) :: out
    type(ReadWrites_t) :: rw
    ! ------------- Locals
    type(FileInfo_t) :: info
    integer iband,k

    ! ------------ albsw
    call clear_file_info(info)
    info%vname = 'albsw'
    info%long_name = 'Carrer soil albedo shortwave (300-4000 nm) annual mean '//sLAI_YEAR
    info%units = '1'
    info%file_metadata_type = 'carrer'
    call out%chunker%nc_create1(out%io_albsw, weighting(out%wta,1d0,0d0), &
        'soilalbedo/', 'soilalbedo_'//trim(out%sres)//'_EntGVSD_v1.1_CarrerGISS_SW_annual_'//sLAI_YEAR, info, &
        create_lr=out%create_lr)

    ! ----------- albgiss
    do iband=1,NBANDS_GISS
        call clear_file_info(info)
        info%vname = 'albgiss_'//trim(sbands_giss(iband))
        info%long_name = 'Carrer soil albedo '//trim(sbands_giss_long(iband))//' annual mean '//sLAI_YEAR
        info%units = '1'
        info%file_metadata_type = 'carrer'
        call out%chunker%nc_create1(out%io_albgiss(iband), weighting(out%wta,1d0,0d0), &
            'soilalbedo/', &
            'soilalbedo_'//trim(out%sres)//'_EntGVSD_v1.1_CarrerGISS_'//trim(sbands_giss(iband))//'_annual_'//sLAI_YEAR, &
            info, create_lr=out%create_lr)
    end do

    ! ----------- fracbd
    do iband=1,NBANDS_GISS
        call clear_file_info(info)
        info%vname = 'fracbd_'//trim(sbands_giss(iband))
        info%long_name = 'Bright/Dark Soil (GISS '//trim(sbands_giss(iband))//' band)'
        info%units = '1'
        info%file_metadata_type = 'carrer'

        call out%chunker%nc_create1(out%ioall_fracbd(iband), weighting(out%wta_fracbd,1d0,0d0), &
            'soilalbedo/', &
            'soilalbedo_'//trim(out%sres)//'_EntGVSD_v1.1_CarrerGISS_fracbd_'//trim(sbands_giss(iband))//'_annual_'//sLAI_YEAR, &
            info, sbright_dark, sbright_dark_long, create_lr=out%create_lr)
        do k=1,NBRIGHT_DARK
            call out%chunker%nc_reuse_var( &
                out%ioall_fracbd(iband), out%io_fracbd(k,iband), &
                (/1,1,k/), weighting(out%wta_fracbd,1d0,0d0))
        end do

    end do

    ! -------------- fracgrey
    call clear_file_info(info)
    info%vname = 'fracgrey'
    info%long_name = 'Bright/Dark Soil in (grey avg of GISS bands)'
    info%units = '1'
    info%file_metadata_type = 'carrer'
    call out%chunker%nc_create1(out%ioall_fracgrey, weighting(out%wta_fracbd,1d0,0d0), &
        'soilalbedo/', &
        'soilalbedo_'//trim(out%sres)//'_EntGVSD_v1.1_CarrerGISS_fracgrey_annual_'//sLAI_YEAR, info, &
        sbright_dark, sbright_dark_long, create_lr=out%create_lr)
    do k=1,2
        call out%chunker%nc_reuse_var( &
            out%ioall_fracgrey, out%io_fracgrey(k), &
            (/1,1,k/), weighting(out%wta_fracbd,1d0,0d0))
    end do

    ! ---------------- bs_brightratio_hr (1km resolution) 
    call clear_file_info(info)
    info%vname = 'bs_brightratio'
    info%long_name = 'Bright Fraction of Soil'
    info%units = '1'
    info%file_metadata_type = 'carrer'
    call out%chunker%nc_create1( &
        out%io_bs_brightratio, weighting(out%wta_fracbd,1d0,0d0), &
        'soilalbedo/', &
        'soilalbedo_'//trim(out%sres)//'_bs_brightratio', info, create_lr=out%create_lr)

    call out%chunker%nc_check(rw=rw)

end subroutine open_output_files


end module B07_mod

!------------------------------------------------------------------------------
      
program Carrer_soilalbedo_to_GISS

    use B07_mod
    use carrer_mod
    use netcdf
    use chunker_mod
    use ent_labels_mod
    use ent_params_mod
    use assign_laimax_mod
    use hntr_mod

!    use convertnc
!    use netcdf
!    use EntGVSD_netcdf_util
    
    implicit none


    integer :: band(NBANDS_GISS)

    integer :: i,k, m, b, dimlon, dimlat, dim(2), dimband
    integer :: dg
    real*4 :: s

    type(Chunker_t) :: ichunker     ! INPUT: Native 6km resolution
    integer, parameter :: NOUTS = 4
    type(Output_t), target :: outs(NOUTS)    ! OUTPUTS: 6km, 1km, HXH, 2HX2 resolutions
    type(Output_t), pointer :: o6    ! => outs(1), which is the (primary) 6km version
    ! Input files
    type(ChunkIO_t) :: ioall_lc, io_lcice, io_lcwater
    real*4 :: lcice, lcwater
    type(ChunkIO_t) :: io_albmodis(NBANDS_MODIS)
    real*4 :: albmodis(NBANDS_MODIS)
    ! Output Files
    real*4 :: albsw
    real*4 :: albgiss(NBANDS_GISS)
    real*4 :: bsbr

    INTEGER :: ichunk,jchunk,ic,jc
    integer :: i1km,j1km   ! Indexing into hi-res 1km variable
    integer :: iband,istat
    type(EntSet_t) :: ent2
    type(FileInfo_t) :: info
    integer, parameter :: nhr(2) = (/ IM1km/IMK, JM1km/JMK /)   ! # of hi-res gridcells per lo-res
    type(ReadWrites_t) :: rw
    call rw%init(THIS_OUTPUTS_DIR, 'B07_soil_albedo', 30,30)

    call init_ent_labels
    call B07_init

    call ichunker%init(IMK, JMK, IMH*2,JMH*2, 'forplot', 100, 100, 10,(/18,15/), &
        outputs_dir=THIS_OUTPUTS_DIR)    ! READ-ONLY
    call init_output(outs(1), ichunker, '6km',  IMK, JMK,.true.)   ! Special: No regrid for this one
    o6 => outs(1)
#ifdef ENTGVSD_DEBUG
    call init_output(outs(2), ichunker, 'xkm',  IM2,JM2,.false.)
#else
    call init_output(outs(2), ichunker, '1km',  IM1km, JM1km,.false.)
#endif
    call init_output(outs(3), ichunker, 'HXH',  IMH, JMH,.false.)
    call init_output(outs(4), ichunker, '2HX2', IM2, JM2,.false.)

    ! ------------- Open Input Files

    ! --------- LC
    ent2 = make_ent2()
    call ichunker%file_info(info, ent2, LAI_SOURCE, 'M', 'lc', LAI_YEAR, 'ent17', '1.1')
    call ichunker%nc_open(ioall_lc, &
        ichunker%outputs_dir, trim(info%dir), trim(info%leaf)//'.nc', trim(info%vname), 0)
    call ichunker%nc_reuse_var(ioall_lc, io_lcice, (/1,1,ent2%svm(SNOW_ICE)/))
    call ichunker%nc_reuse_var(ioall_lc, io_lcwater, (/1,1,ent2%svm(CV_WATER)/))

    ! ------------ albmodis
#ifdef USE_FILLED
    print *, 'Using grid filled soil albedo mean.'
    do iband=1,NBANDS_MODIS
        ! Read from 2D NetCDF var
        call ichunker%nc_open(io_albmodis(iband), ichunker%outputs_dir, &
            'tmp/carrer/', &
            'albfill_'//trim(sbands_modis(iband))//'.nc', &
            'albfill_'//trim(sbands_modis(iband)), SMEAN)
    end do
#else
    print *, 'Using non-filled soil albedo mean.'
    do iband=1,NBANDS_MODIS
        ! Read from 3D NetCDF var
        call ichunker%nc_open(io_albmodis(iband), ichunker%outputs_dir, &
            'tmp/carrer/', &
            'albmodis_'//trim(sbands_modis(iband))//'.nc', &
            'albmodis_'//trim(sbands_modis(iband)), SMEAN)
    end do
#endif

    ! ===================== Open Output Files

    do i=1,NOUTS
print *,'***************** open_output_files',i
        call open_output_files(outs(i), rw)
    end do

    call ichunker%nc_check(rw=rw)
    call rw%write_mk

#ifdef JUST_DEPENDENCIES
    STOP 0
#endif

    ! ================== Main Loop

#ifdef ENTGVSD_DEBUG
    do jchunk = dbj0,dbj1+1
    do ichunk = 1,ichunker%nchunk(1)
#else
    do jchunk = 1,ichunker%nchunk(2)
    do ichunk = 1,ichunker%nchunk(1)
#endif

        call ichunker%move_to(ichunk,jchunk)
        do i=1,NOUTS
            call outs(i)%chunker%move_to(ichunk,jchunk)
        end do
        o6%wta = 0

        do jc = 1,o6%chunker%chunk_size(2)
        do ic = 1,o6%chunker%chunk_size(1)

            ! ---------- Read inputs
            lcice = io_lcice%buf(ic,jc)
            lcwater = io_lcwater%buf(ic,jc)
            do iband=1,NBANDS_MODIS
                albmodis(iband) = io_albmodis(iband)%buf(ic,jc)
            end do

            ! Compute overall NetCDF index of current cell
            !ii = (ichunk-1)*chunker%chunk_size(1)+(ic-1)+1
            !jj = (jchunk-1)*chunker%chunk_size(2)+(jc-1)+1

            ! Infer soil mask
            ! TODO: Get this from LC instead
            if (albmodis(1) /= FillValue) then
                o6%wta(ic,jc) = 1
            else
                o6%wta(ic,jc) = 0
            end if


            !* Calculate total SW albedo --------------------------------------
            if ((albmodis(VIS_MODIS).eq.FillValue).or. &
                (albmodis(NIR_MODIS).eq.FillValue)) then
                albsw = FillValue
            else
                albsw = ( &
                     albmodis(VIS_MODIS) * sum(fracSW_MG(1:2)) &
                   + albmodis(NIR_MODIS) * sum(fracSW_MG(3:8)) &
                ) / sum(fracSW_MG(1:8)) 
            endif

            !* Put spectral breakdown into GISS bands ------------------------
            do iband=1,NBANDS_GISS
                albgiss(iband) = FillValue
            end do

            if ((albmodis(VIS_MODIS).eq.FillValue).or. &
                (albmodis(NIR_MODIS).eq.FillValue)) then
                do iband=1,NBANDS_GISS
                    albgiss(iband) = FillValue
                end do
            else
                albgiss(VIS_GISS) = &
                   ( albmodis(VIS_MODIS) * &
                   (fracSW_MG(nm300_400) + fracSW_MG(nm400_700)) + &
                   albmodis(NIR_MODIS)*fracSW_MG(nm700_770)) &
                   /( fracSW_MG(nm300_400) + fracSW_MG(nm400_700) + &
                   fracSW_MG(nm700_770) )
                do k=2,NBANDS_GISS
                    albgiss(k) = albmodis(NIR_MODIS)
                end do
            endif

            !------ Old VEG bare_bright vs. bare_dark soil fractions
            !* Calculate GISS ModelE grey albedo bright vs. dark fractions
            !Same netcdf names as in Ent_pfts.
            !1) Non-grey albedos: each GISS band has its own bright|dark fractions.
            !  For partial-land cells (that have ice and/or water) use
            ! soil albedo of nearest all-land neighbor.  If none, default to
            ! to cell's mixed albedo.
            do iband=1,NBANDS_GISS
            do k=1,2
                o6%io_fracbd(k,iband)%buf(ic,jc) = FillValue
            end do
            end do
            do k=1,NBANDS_GISS
                if (albgiss(k) /= FillValue) then
                    if (abs(1d0 - (lcice+lcwater)) < 1d-5) then
                        ! ------ All permanent ice, no soil.
                        o6%io_fracbd(BRIGHT,k)%buf(ic,jc) = FillValue
                        o6%io_fracbd(DARK,k)%buf(ic,jc) = FillValue
                    else if (((lcice+lcwater) > 0.).and. &
                        ((lcice+lcwater) < 1.0)) then
                        ! --------- Partial ground
                        !Search for nearest grid cell up to 5 that is all ground -- NOT DONE
                        s = FillValue

                        if (s.eq.FillValue) then !No nearby ground cells found 
                            !io_fracbd(BRIGHT,k)%buf(ic,jc) = FillValue
                            !io_fracbd(DARK,k)%buf(ic,jc) = FillValue
                            !*noice - If cannot assign land albedo, cut the fractional cells.
                            !  Loses Alaska peninsula
                            !fracgrey(i,j,:) = FillValue
                            !*fringeice - Or, since there is some ground, assign i,j value there.
                            !  Leaves some coastal ice fringe.
                            o6%io_fracgrey(BRIGHT)%buf(ic,jc) = &
                                min(0.5, albgiss(k))/0.5 !
                            o6%io_fracgrey(DARK)%buf(ic,jc) = 1.0 - o6%io_fracgrey(BRIGHT)%buf(ic,jc)
                        else
                            o6%io_fracbd(BRIGHT,k)%buf(ic,jc) = min(s,0.5)/0.5
                            o6%io_fracbd(DARK,k)%buf(ic,jc) = 1.0 - o6%io_fracbd(BRIGHT,k)%buf(ic,jc)
                        end if
                    else
                        ! -------- All ground
                        o6%io_fracbd(BRIGHT,k)%buf(ic,jc) = min(albgiss(k),0.5)/0.5
                        o6%io_fracbd(DARK,k)%buf(ic,jc) = 1.0 - o6%io_fracbd(BRIGHT,k)%buf(ic,jc)
                    end if
                else
                    o6%io_fracbd(BRIGHT,k)%buf(ic,jc) = FillValue
                    o6%io_fracbd(DARK,k)%buf(ic,jc) = FillValue
                end if
            end do   ! k=1,NBANDS

            if (o6%io_fracbd(1,1)%buf(ic,jc) == FillValue) then
                o6%wta_fracbd(ic,jc) = 0
            else
                o6%wta_fracbd(ic,jc) = 1
            end if

            !2) Grey shortwave albedos - exclude 100% ice+water, nearest-neighbor fill
            !  For partial-land cells (that have water) use
            ! soil albedo of nearest all-land/non-ice neighbor.  If none, default to
            ! to cell's mixed albedo.
            o6%io_fracgrey(BRIGHT)%buf(ic,jc) = 0
            o6%io_fracgrey(DARK)%buf(ic,jc) = 0

            if (albsw.ne.FillValue) then 
               if ((lcice+lcwater).eq.1.0) then !No ground
                    o6%io_fracgrey(BRIGHT)%buf(ic,jc) = FillValue
                    o6%io_fracgrey(DARK)%buf(ic,jc) = FillValue
               elseif (((lcice+lcwater).gt.0.).and. &
                      ((lcice+lcwater).lt.1.0)) then !Partial
                  !Search for nearest grid cell up to 5 that is all ground
                  s = FillValue
                  if (s.eq.FillValue) then !No nearby ground cells found
!                    !*If cannot assign land albedo, cut the fractional cells.
!                    !Loses Alaska peninsula
                     !fracgrey(i,j,:) = FillValue
!                    !*Or, since there is some ground, assign i,j value there.
!                    !Leaves some coastal ice fringe.
                     ! *** This is where had to cap bright/dark fraction to .5
                     o6%io_fracgrey(BRIGHT)%buf(ic,jc) = min(0.5, albsw)/0.5 !
                     o6%io_fracgrey(DARK)%buf(ic,jc) = 1d0 - o6%io_fracgrey(BRIGHT)%buf(ic,jc)
                  else !Found nearby all-ground cell
                     o6%io_fracgrey(BRIGHT)%buf(ic,jc) = min(0.5,s)/0.5
                     o6%io_fracgrey(DARK)%buf(ic,jc) = 1d0 - o6%io_fracgrey(BRIGHT)%buf(ic,jc)
                  endif
               else    !Cell is all ground with SW
                  o6%io_fracgrey(BRIGHT)%buf(ic,jc) = min(0.5, albsw)/0.5
                  o6%io_fracgrey(DARK)%buf(ic,jc) = 1d0 - o6%io_fracgrey(BRIGHT)%buf(ic,jc)
               endif
            else
                o6%io_fracgrey(BRIGHT)%buf(ic,jc) = FillValue
                o6%io_fracgrey(DARK)%buf(ic,jc) = FillValue
            endif

            ! ------------ Regrid fracgrey
            o6%io_bs_brightratio%buf(ic,jc) = o6%io_fracgrey(BRIGHT)%buf(ic,jc)

            ! ---------- Store outputs
            o6%io_albsw%buf(ic,jc) = albsw
            do k=1,NBANDS_GISS
                o6%io_albgiss(k)%buf(ic,jc) = albgiss(k)
            end do


        end do
        end do

        ! ------------ Regrid 6km to other resolutions
        do i=2,NOUTS
            ! ------------ Regrid to GISS resolution

            ! albgiss
            do k=1,NBANDS_GISS
                call outs(i)%hntr%regrid4( &
                    outs(i)%io_albgiss(k)%buf, o6%io_albgiss(k)%buf, &
                    o6%wta, 1d0, 0d0, &    ! weighting,
                    outs(i)%io_albgiss(k)%startB(2), outs(i)%io_albgiss(k)%chunker%chunk_size(2))
            end do

            ! fracgrey
            do k=1,2
                call outs(i)%hntr%regrid4( &
                    outs(i)%io_fracgrey(k)%buf, o6%io_fracgrey(k)%buf, &
                    o6%wta_fracbd, 1d0, 0d0, &    ! weighting,
                    outs(i)%io_fracgrey(k)%startB(2), outs(i)%io_fracgrey(k)%chunker%chunk_size(2))
            end do

            ! fracbd
            do iband=1,NBANDS_GISS
            do k=1,NBRIGHT_DARK
                call outs(i)%hntr%regrid4( &
                    outs(i)%io_fracbd(k,iband)%buf, o6%io_fracbd(k,iband)%buf, &
                    o6%wta_fracbd, 1d0, 0d0, &    ! weighting,
                    outs(i)%io_fracbd(k,iband)%startB(2), outs(i)%io_fracbd(k,iband)%chunker%chunk_size(2))
            end do
            end do

            call outs(i)%hntr%regrid4( &
                outs(i)%io_albsw%buf, o6%io_albsw%buf, &
                o6%wta, 1d0, 0d0, &    ! weighting,
                outs(i)%io_albsw%startB(2), outs(i)%io_albsw%chunker%chunk_size(2))

            call outs(i)%hntr%regrid4( &
                outs(i)%io_bs_brightratio%buf, o6%io_bs_brightratio%buf, &
                o6%wta_fracbd, 1d0, 0d0, &    ! weighting,
                outs(i)%io_bs_brightratio%startB(2), outs(i)%io_bs_brightratio%chunker%chunk_size(2))


        end do    ! i=2,NOUTS

        do i=1,NOUTS
            call outs(i)%chunker%write_chunks
        end do
    end do
    end do


    call ichunker%close_chunks
    do i=1,NOUTS
        call outs(i)%chunker%write_chunks
    end do

end program Carrer_soilalbedo_to_GISS
            
