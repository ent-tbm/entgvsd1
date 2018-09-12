module ioutil_mod

use entgvsd_config_mod
use, intrinsic :: iso_fortran_env
use netcdf

implicit none


! Hardcoded directories for input / output

! ---------- Original values (Carlos on discover)
!      --- Input (read-only)
!      CHARACTER (LEN=*), PARAMETER :: DATA_DIR='../../data/'
!      !CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR='../lc_lai_giss/'
!      CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR=
!     &     '../../lc_lai_for_1km1km/'
!      --- Output
!      CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR='../lc_lai_ent/'

! Values for Elizabeth's work on gibbs
!      --- Input (read-only)
CHARACTER (LEN=*), PARAMETER :: DATA_DIR= &
   '/home2/rpfische/entgvsd0/discover/Vegcover_1km/data/'

CHARACTER (LEN=*), PARAMETER :: DATA_INPUTS= &
   '/home2/rpfische/git/entgvsd/inputs/data/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_DIR= &
   '/home2/rpfische/entgvsd0/discover/Vegcover_1km/BNU/'// &
   'lc_lai_giss'

CHARACTER (LEN=*), PARAMETER :: LC_LAI_GISS_INPUTS= &
   '/home2/rpfische/git/entgvsd/inputs/lc_lai_giss/'



CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_DIR= &
    '/home2/rpfische/entgvsd0/discover/Vegcover_1km/'// &
    'lc_lai_for_1kmx1km/'
CHARACTER (LEN=*), PARAMETER :: LC_LAI_FOR_1KM1KM_INPUTS= &
   '/home2/rpfische/git/entgvsd/inputs/lc_lai_giss/lc_lai_for_1kmx1km/'




!      --- Output
CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_DIR= &
   '/home2/rpfische/git/entgvsd/bnu/lc_lai_ent/'

CHARACTER (LEN=*), PARAMETER :: TEMPLATE_DIR= &
   '/home2/rpfische/git/entgvsd/templates/'


! Place of output files originally written by carlos
! Used to create templates.

CHARACTER (LEN=*), PARAMETER :: LC_LAI_ENT_ORIG= &
   '/home2/rpfische/entgvsd0/discover/Vegcover_1km/BNU/lc_lai_ent'

! ===========================================================


! Keep track if we've globally seen an error
integer :: entgvsd_errs = 0

contains
! ---------------------------------------

function create_nc(dir,leaf)
    character(*), intent(in) :: dir
    character(*), intent(in) :: leaf
    integer :: create_nc
    ! -------- Local vars
    character(*), parameter :: ncgen = &
        NETCDF4_C_BINARY_DIR//'/ncgen'
    integer :: ncid,err
    character(4192) :: cmd


    create_nc = -1

    ! ------- Create the file from a template

    ! cdl = ENTGVSD_PROJECT_SOURCE_DIR//"/templates/"//cdl_leaf
    ! cmd = NCGEN//' -k nc4 -o '//ofname//' '//trim(cdl)
    cmd = 'entgvsd_create_nc '//LC_LAI_ENT_DIR//' '//LC_LAI_ENT_ORIG//' '//TEMPLATE_DIR// &
          ' '//trim(dir)//' '//trim(leaf)

    call execute_command_line(cmd, .true., err)

    if (err /= 0) then
        write(ERROR_UNIT,*) 'Error running command to create', &
            LC_LAI_ENT_DIR//'/'//LC_LAI_ENT_ORIG//'/'//TEMPLATE_DIR// &
            '/'//trim(dir)//'/'//trim(leaf)

        entgvsd_errs = entgvsd_errs + 1
        return
    end if

    ! ------- Open the file we just created
    err = NF90_OPEN(LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'.nc', NF90_WRITE, create_nc)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening after create',LC_LAI_ENT_DIR//trim(dir)//trim(leaf)//'.nc',err
        entgvsd_errs = entgvsd_errs + 1
        return
    end if
end function create_nc


subroutine open_input_nc(iroot, oroot, dir, leaf, ncid)
    character*(*), intent(in) :: iroot   ! Read (maybe compressed) file from here
    character*(*), intent(in) :: oroot   ! Write or link to here, then open
    character*(*), intent(in) :: dir
    character*(*), intent(in) :: leaf
    integer, intent(out) :: ncid
    ! --------- Local vars
    character(4192) :: cmd
    integer :: err
    logical :: exist

    ncid = -1

    ! -------- Decompress / link the file if it doesn't exist
    inquire(FILE=oroot//dir//leaf, EXIST=exist)
    if (.not.exist) then
        ! cdl = ENTGVSD_PROJECT_SOURCE_DIR//"/templates/"//cdl_leaf
        ! cmd = NCGEN//' -k nc4 -o '//ofname//' '//trim(cdl)
        cmd = 'entgvsd_link_input '//trim(iroot)//' '//trim(oroot)//' '//trim(dir)//' '//trim(leaf)

        call execute_command_line(cmd, .true., err)
        if (err /= 0) then
            write(ERROR_UNIT,*) 'Error running entgvsd_link_input',leaf,err
            entgvsd_errs = entgvsd_errs + 1
            return
        end if
    end if

    ! ------- Now open the file
    err = nf90_open(oroot//dir//leaf,NF90_NOWRITE,ncid)
    if (err /= NF90_NOERR) then
        write(ERROR_UNIT,*) 'Error opening',trim(leaf),err
        entgvsd_errs = entgvsd_errs + 1
        return
    end if

end subroutine open_input_nc


subroutine check_nf_open_errsors
    if (entgvsd_errs > 0) then
        write(ERROR_UNIT,*) 'Previous errors encountered creating/opening files'
        stop -1
    end if
end subroutine check_nf_open_errsors

end module ioutil_mod