!trim_EntMM_monthly_05x05_nc.f - from trim_EntMM_monthly_05x05_nc.f.  Outputs netcdf.
!     This program converts EntMM 17 PFTs to Ent 16 PFTs + bright + dark.
! Converts lc, laimax, monthly lc and lai, and Simard heights.
! Original from trim_EntMM_monthly_noht.f, which did not do heights.
! It does NOT interpolate from fine to coarse grid -- this should be done
!  prior to using this program, so input and output are the same resolution.
! To change input/output resolution edit lines below 
! "define input file resolution" and "new (interpolated) values"
! To combine C3 and C4 crops: #define COMBINE_CROPS_C3_C4
! 9/12/13 Fixed bright/dark soil:  need to do partitioning after
!         each step of trim/scale/no crops to account for new
!         bare soil cover added, especially in crop grid cells.
! 1/16/13 Added subroutine replace_crops for nocrops to replace cover and LAI 
!         of crops with dominant natural veg in grid cell. If no natural veg,
!         then searches adjacent grid cells.
!         Added subroutine fill_crops for _ext version of crop LAI, to fill
!         in crop LAI in some grid cells that have no crops in MODIS cover
!         but may have crop cover in Pongratz historical cover.  Called once
!         for single grid extension.  Can be called again for further filling.
! 3/17/14 Added ext1 files at end for lai max, lai monthly, and height, by
!         replacing 15-crops in laic, laim, hm, and hsd with the ext values for
!         crops.
! 12/12/17 Fixed merge of C3 and C4 crops for revised Ent17 input files that have WATER at the 
!         end rather than beginning of list of cover types.  Replaced hard-coded array numbers with
!         variable indices c3 and c4 so that these can be updated more easily.
! 12/15/17 Added netcdf outputs, including new subroutines in create netcdf files.
! 12/20/17 Fixed convert_vf so that it keeps original vf2 if zero LAI passed in.
!     This got rid of the mismatch in vfc and vfm after trimming, and solved a lot
!     of other errors down the line.  Still some tiny round-off differences, but ok.
!     Got rid of "NO WATER" layer in heights arrays h and hsd, so that cover indices
!     match up with vfc, vfm, laic, and laim
!12/21/17 Finished adding netcdf outputs. Tried to output monthly LAI to one file
!     with time dimension for months, but cannot get nf_put_var to work, despite
!     status=NO_ERROR.  . ALSO, should have fill_crops extend by more
!     than 1 grid cell for the "ext" files. 10/28/2018 - Abandoning arrays with
!     time dimension because better for GCM to read in all layers from one file
!     at a time point, rather than read one time point from many files.
!12/21/17  Changed fill_crops ext grids from 1 to 5 and generated version ext5.
!10/24/2018  MAJOR OVERHAUL from previous version.  Assumes input files are in order
!     of PFTs first and non-veg including WATER last (no skip of first layer that was water).
!     Uses newly-generated Simard heights and
!     Carrer soil albedo input files.  Adding netcdf functionality.
!      
! Compile the program with:
!
! ifort -cpp convert_VEG5.f -convert big_endian
!
! or
! gfortran -cpp -fconvert=big-endian -O0 -fno-range-check arrayutil.f trim_EntMM_monthly_05x05.f
! BEFORE RUNNING:  mkdir lc_lai_ent16
! 
! With netcdf utilities after 12/15/17, compile as:
! gfortran -cpp -fconvert=big-endian -O0 -fno-range-check -I$NETCDFHOME/include -c convertnc_util.f arrayutil.f trim_EntMM_monthly_05x05_nc.f
! gfortran convertnc_util.o arrayutil.o trim_EntMM_monthly_05x05_nc.o -L$NETCDFHOME/lib -lnetcdf -lnetcdff!
! With netcdf and EntGVSD utilities after 10/24/2018, compile as:
!     gfortran  -cpp -fconvert=big-endian -O0 -fno-range-check -I$NETCDFHOME/include -c convertnc_util.f arrayutil.f EntGVSD_util.f trim_EntMM_monthly_05x05_nc.f
!     gfortran convertnc_util.o arrayutil.o EntGVSD_util.o trim_EntMM_monthly_05x05_nc.o -L$NETCDFHOME/lib -lnetcdf -lnetcdff      

#define COMBINE_CROPS_C3_C4
#define SPLIT_BARE_SOIL

      module conversions
      use convertnc
      use netcdf
      use EntGVSD_netcdf_util
      implicit none

!      !Moved to EntGVSD_util.f
!      character*3, parameter :: MONTH(12) =
!     &     (/
!     &     "Jan","Feb","Mar","Apr","May","Jun",
!     &     "Jul","Aug","Sep","Oct","Nov","Dec"
!     &     /)

      !Need to move to EntGVSD_util.f
!      integer, parameter :: N_COVERTYPES = 18
      !character(len=13), parameter :: ent_cover_names(N_COVERTYPES) = (/
!      character(len=13), parameter :: ent_names(N_COVERTYPES) = (/
!     &     "ever_br_early",
!     &     "ever_br_late ",
!     &     "ever_nd_early",
!     &     "ever_nd_late ",
!     &     "cold_br_early",
!     &     "cold_br_late ",
!     &     "drought_br   ",
!     &     "decid_nd     ",
!     &     "cold_shrub   ",
!     &     "arid_shrub   ",
!     &     "c3_grass_per ",
!     &     "c4_grass     ",
!     &     "c3_grass_ann ",
!     &     "c3_grass_arct",
!     &     "crops_herb   ",
!     &     "crops_woody  ",
!     &     "bare_bright  ",
!     &     "bare_dark    "
!     &     /)

      !Need to move to EntGVSD_util.f
!            character*50, parameter :: Ent_title(N_COVERTYPES) =
!     &     (/
!     &     '1 - evergreen broadleaf early succ               ',
!     &     '2 - evergreen broadleaf late succ                ',
!     &     '3 - evergreen needleleaf early succ              ',
!     &     '4 - evergreen needleleaf late succ               ',
!     &     '5 - cold deciduous broadleaf early succ          ',
!     &     '6 - cold deciduous broadleaf late succ           ',
!     &     '7 - drought deciduous broadleaf                  ',
!     &     '8 - deciduous needleleaf                         ',
!     &     '9 - cold adapted shrub                           ',
!     &     '10 - arid adapted shrub                          ',
!     &     '11 - C3 grass perennial                          ',
!     &     '12 - C4 grass                                    ',
!     &     '13 - C3 grass - annual                           ',
!     &     '14 - arctic C3 grass                             ',
!     &     '15 - crops herb                                  ',
!     &     '16 - crops woody                                 ',
!     &     '17 - bright bare soil                            ',
!     &     '18 - dark bare soil                              '
!     &     /)

      contains

      subroutine convert_vf(vf1, lai1, vf2, lai2, laimin)
      !Converts LAI that is on BARE to a vegetated fraction.
      !Reduces vf1, increases vf2 by increment that maintains lai2,
      ! or only increases lai2 if lai1 is so big that the cover-weighted
      ! lai is bigger than lai2.
      !vf1, lai1 - BARE cover fraction and LAI
      !vf2, lai2 - vegetated cover fraction and max LAI
      real*4 vf1, lai1, vf2, lai2, laimin
      !---
      real*4 tot_la, tot_vf, new_lai, new_vf

      tot_la = vf1*lai1 + vf2*lai2
      tot_vf = vf1 + vf2
      new_lai = tot_la/tot_vf
      !Fix for if 100% conversion to vf2.
      if (new_lai.eq.0.0) then
         new_vf = vf2 !0.0 !Keep original cover
      else
         new_lai = max( new_lai, laimin)
         new_vf = tot_la/new_lai
      endif
      ! get rid of round-off errors
      !!new_vf = min( new_vf, tot_vf )

      vf2 = new_vf
      lai2 = new_lai
      vf1 = max(0.,tot_vf - vf2) !Argh, round-off errors
      lai1 = 0.

      end subroutine convert_vf

!               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
!     &              vfm(m,i,j,9),vfm(m,i,j,9), vfc(i,j,9))

      subroutine convert_vfm(vf1, lai1, vf2, lai2, vfc)
      !convert monthly vf using vfc for new vf2 derived from laimax trim. 
      !vf1, lai1 - BARE cover fraction and LAI
      !vf2, lai2 - vegetated cover fraction and monthly LAI
      !vfc = vf2 + dvf1 from laimax trim (dvf1 is positive)
      !vf1 + vf2 = vfc + (vf1-dvf1)
      !vf1*lai1 = 0*(vf1-dvf1) + (dvf1*dlai1)
      !new_lai2 = ((dvf1*dlai1) + (vf2*lai2))/(dvf1+vf2)
      !         = ((vf1*lai1) + (vf2*lai2))/vfc
      real*4, intent(inout) :: vf1, lai1, vf2, lai2
      real*4, intent(in) :: vfc
      !---
      real*4 :: tot_lai, new_lai, new_vf

      new_vf = vfc
      if (new_vf.le.0.) then  
         print *,'new_vf is .leq. zero',vf1,vf2,lai1,lai2,vfc
         STOP
      endif
      new_lai = (vf1*lai1 + vf2*lai2)/new_vf
!      vf2 = new_vf
      lai2 = new_lai
      vf1 = max(0., vf1 - (vfc - vf2))   !bare, max for neg. round-off error
      vf2 = new_vf
      lai1 = 0.

      end subroutine convert_vfm



      subroutine convert_vfh(vf1, h1, hsd1, vf2, h2, hsd2, vfc)
      !convert height data, h and hsd, based on the converted vfc.
      !Since conversion is to shrub and crops, or, if none, to the
      ! next dominant PFT cover type, there is no cover-weighted height
      ! averaging, but the h and hsd values are simply preserved (this
      ! routine effectively does nothing, since there is no height
      ! over BARE of ICE in EntGVSDmosaic_05x05.ij.
      !However, for the last conversion of any remaining BARE laimax to
      ! arid-adapted shrub cover that does not pre-exist, a new height
      ! has to be assigned.  ##DOUBLEBCHECK simard.f for height assigned
      ! to arid shrubs!
      !vf1, h1, hsd1 - BARE cover fraction, height, and height stdev
      !vf2, h2, hsd2  - vegetated cover fraction, heigth, and height stdev
      !vfc - new vegetated cover fraction that was derived from laimax.
      !vfc = vf2 + dvf1 from laimax trim (dvf1 is positive)
      !vf1 + vf2 = vfc + (vf1-dvf1)
      !vf1*h1 = 0*(vf1-dvf1) + (dvf1*dlai1)
      !new_lai2 = ((dvf1*dlai1) + (vf2*lai2))/(dvf1+vf2)
      !         = ((vf1*lai1) + (vf2*lai2))/vfc
      real*4, intent(inout) :: vf1, h1, hsd1, vf2, h2, hsd2
      real*4, intent(in) :: vfc
      !---
      real*4 tot_lai, new_h, new_hsd, new_vf2

      new_vf2 = vfc
      if (new_vf2>0.) then !vf2 is present, keep h2
         new_h = h2
         new_hsd = hsd2
      else
         new_h = h1
         new_hsd = hsd1
      endif
      vf2 = new_vf2
      h2 = new_h
      vf1 = max(0., vf1 - (vfc - vf2)) !bare, max for neg. round-off error
      h1 = 0.
      hsd1 = 0.

      end subroutine convert_vfh


      !Need to move to convertnc_util.f
!      subroutine handle_nf_error(status, message)
!      use netcdf
!      integer, intent(in) :: status
!      character*(*) :: message
!      
!       write(*,*) status,message, ' '
!     &     ,trim(nf90_strerror(status))!      if(status /= nf90_NoErr)  then
!        write(*,*) 'DONT FORGET TO DELETE THE OLD NETCDF FILES FIRST'
!	STOP
!      endif
!      end subroutine handle_nf_error

      !Keep local to each program to define own values.
      subroutine my_nf_defglobal(file)
      use netcdf
      use convertnc
      character*(*) :: file
      !---Local---
      integer :: status, ncid
      character*1024 :: text
      character(8) :: date
      
           !Enter define mode
      status=nf_open(trim(file), NF_WRITE, ncid)
      write(0,*) 'nf_open ', status, file, ncid
      call handle_nf_error(status,  'nf_open '//file)       
      status=nf_redef(ncid)
      write(0,*) 'nf_redef ', status, ncid
      !call handle_nf_error(status,  'nf_redef')
      !Put metadata global attributes
      text = 'Ent Global Vegetation Structure Dataset '//
     &     '(Ent GVSD) v1.0b  MODIS-Monfreda '
      !write(0,*) 'len, text: ', len(trim(text)), trim(text)
      status=nf_put_att_text(ncid, NF_GLOBAL, 'Description'
     &     ,len(trim(text)), trim(text))
      call handle_nf_error(status, 'nf_global'//trim(text))
      text = 'TEMPORARY working version downscaled '//
     &     'bare/bright soil from 2.5x2 degrees. '//
     &     'With ext1 files having crops LAI extended by 5 grid cells.'
      !write(0,*) 'len, text: ',len(trim(text)), trim(text)
      status=nf_put_att_text(ncid, NF_GLOBAL, 'Comments'
     &     ,len(trim(text)), trim(text))
      call handle_nf_error(status, '')
      text = 'Institution:  NASA Goddard Institute for Space Studies'
      status=nf_put_att_text(ncid, NF_GLOBAL, 'Institution'
     &     ,len(trim(text)), trim(text))
      call handle_nf_error(status, '')
      text = 'Nancy.Y.Kiang@nasa.gov'
      status=nf_put_att_text(ncid, NF_GLOBAL, 'Contact'
     &     ,len(trim(text)), trim(text))
      call handle_nf_error(status, '')
      call DATE_AND_TIME(date)
      text = date//' Created'
      status=nf_put_att_text(ncid, NF_GLOBAL, 'History'
     &     ,len(trim(text)), trim(text))
      call handle_nf_error(status,  'nf_put_att NF_GLOBAL')

!      status=nf_enddef(ncid)      
!      call handle_nf_error(status,  'nf_enddef my_nf_defglobal')

      end subroutine my_nf_defglobal


      subroutine write_output_lai(IM,JM,titlec, laic, n, fileprefix
     &     , MISC, MON, resoutt)
      !Same as write_output, but only lai, not lc.
      integer, intent(in) :: IM, JM
      character*80 :: titlec(:)
      real*4 :: laic(:,:,:)
      integer :: n
      character*(*) :: fileprefix
      character*(*) :: MISC, MON
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      character*4  :: MONstr
      integer :: k
      character*256 :: fileout
      integer :: dimlon,dimlat,status, ncid, varid

      if (trim(MON)=="") then
         MONstr=""
      else
         MONstr="_"//MON
      endif

      open(90,file=fileprefix//"_lai_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!        read( titlec(k), '(a4, a46)' ) dum, pft_name
!        write( title, '(i2, " - ", a46, " (LAI) 4x5")' )
!     &       k, adjustl(pft_name)
!        write(90) title, laic(:,:,k)
         title = titlec(k)(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
         write(90) title, laic(:,:,k)
      enddo

      close(90)

      !* Netcdf output
      fileout = fileprefix//"_lai_"//trim(MISC)//trim(MONstr)//".nc"
      call my_nf_create_Ent(IM,JM
     &     ,trim(fileout)
     &     ,n,ent_names16, Ent_title16,
     &     "leaf area index m^2 leaf / m^2 ground", ncid)
      call my_nf_defglobal(trim(fileout))
!     &     ,'Ent 16 PFTs MODIS-Monfreda 2004 leaf area index (LAI)'
!     &     , '  ')
      status = my_nf_open(trim(fileout), NF_WRITE, ncid)
!      call my_nf_create_ij(
!     &     fileprefix//"_lai_"//trim(MISC)//trim(MONstr)//".nc", IM,JM
!     &     ,ncid,dimlon,dimlat)
      do k=1,n
!         status = my_nf_inq_def_put_var_real32_2(ncid,IM,JM,1,IM,1,JM,
!     &        dimlon,dimlat, trim(ent_names16(k))//'_lai',
!     &        trim(Ent_title16(k))//
!     &        ' leaf area index', 'm^2 leaf / m^ ground', 
!     &        laic(:,:,k))
         status = my_nf_inq_put_var_real32_2(ncid,
     &        trim(ent_names16(k)), varid,laic(:,:,k))
      end do
      status=nf90_close(ncid)
      call handle_nf_error(status, 'write_output_lai '//
     &     trim(fileprefix//"_lai_"//trim(MISC)//trim(MONstr)//".nc"))
      
      end subroutine write_output_lai

      
      subroutine write_output_lai_monthly_nc(IM,JM,ncov
     &     , laim, fileprefix
     &     , MISC, units, resoutt)
      !Create and write to 3D arrays of Ent cover variables, ixjxmonth.
      use netcdf
      integer,intent(in) :: IM,JM,ncov
      real*4 :: laim(:,:,:,:)
      character*(*) :: fileprefix
      character*(*) :: MISC
      character*(*) :: units
      character*(*) :: resoutt !resout for titles
      !---
      character*1024 :: file
      integer :: k, m
      integer :: status, ncid, varid
      real*4 :: laimrev(IM,JM,12) !Need this order for output to netcdf
      
      file=fileprefix//"_lai_"//trim(MISC)//".nc"

      call my_nf_create_Ent_vartime(IM,JM,trim(file),ncov,
     &     ent_names16, Ent_title16, units,ncid)
      status = my_nf_open(trim(file), NF_WRITE, ncid)

      do k=1,ncov
         do m=1,12
            laimrev(:,:,m) = laim(m,:,:,k)
         enddo
         status=my_nf_inq_put_var_real32_3(ncid,
     &     trim(ent_names16(k)),varid,laimrev(:,:,:))
         call handle_nf_error(status, 'my_nf_put_var_rea32_3 '//
     &        trim(ent_names16(k)))
      end do

      status=nf_close(ncid)
      call handle_nf_error(status, 'nf_close '//trim(file))

      end subroutine write_output_lai_monthly_nc


      subroutine write_output_single(IM,JM,titlefoo, varname, vf, lai,
     &     fileprefix, MISC, MON, resoutt)
      !Same as write-output, but arrays are for only one veg type.
      integer,intent(in) :: IM,JM
      character*80 :: titlefoo
      character*(*) :: varname
      real*4 :: vf(:,:), lai(:,:)
      character*(*) :: fileprefix
      character*(*) :: MISC, MON
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title
      character*4  :: MONstr
      character*256 :: fileout
      integer :: k, status,varid, ncid

      if (trim(MON)=="") then
         MONstr=""
      else
         MONstr="_"//MON
      endif

      open(80,file=fileprefix//"_lc_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      title = titlefoo(1:48)//"  "//trim(MON)//" (cover)  "//resoutt
      write(80) title, vf(:,:)
      close(80)

      open(90,file=fileprefix//"_lai_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      title = titlefoo(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
      write(90) title, lai(:,:)

      close(90)

!     !* Netcdf output
      fileout = fileprefix//"_lc_"//trim(MISC)//trim(MONstr)//".nc"
      call my_nf_create_Ent_single(IM,JM
     &     ,trim(fileout)
     &     ,varname,trim(MISC),'cover fraction', ncid)
      call my_nf_defglobal(trim(fileout))
      status = my_nf_open(trim(fileout), NF_WRITE,ncid)
      status = my_nf_inq_put_var_real32_2(ncid,
     &        varname, varid,vf(:,:))
      status=nf90_close(ncid)

      fileout = fileprefix//"_lai_"//trim(MISC)//trim(MONstr)//".nc"
      call my_nf_create_Ent_single(IM,JM
     &     ,trim(fileout)
     &     ,varname,trim(MISC), 'm^2 leaf / m^ground', ncid)
      call my_nf_defglobal(trim(fileout))
      status = my_nf_open(trim(fileout), NF_WRITE,ncid)
      status = my_nf_inq_put_var_real32_2(ncid,
     &        varname, varid,lai(:,:))
      status=nf90_close(ncid)

      end subroutine write_output_single


      subroutine write_output(IM,JM,titlec, vfc, laic, n
     &     , fileprefix, MISC, MON, resoutt)
      !LC and LAI, separate files each with all cover types
      !Output monthly files, single file x all cover for each month.
      use netcdf
      integer, intent(in) :: IM, JM
      character*80 :: titlec(:)
      real*4 :: vfc(:,:,:), laic(:,:,:)
      integer :: n  !Last layer number
      character*(*) :: fileprefix
      character*(*) :: MISC, MON
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      character*4  :: MONstr
      character*256 :: fileout
      integer :: k
      !-- Netcdf vars
      integer :: ncid, status, dimid, varid
      integer, parameter :: start(2) = (/ 1, 1 /)
      integer  :: count(2)
      count  = (/ IM, JM /)
      
      if (trim(MON)=="") then
         MONstr=""
      else
         MONstr="_"//MON
      endif

      open(80,file=fileprefix//"_lc_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!        read( titlec(k), '(a4, a46)' ) dum, pft_name
!        write( title, "(i2, ' - ', a46, ' (cover) 4x5')" )
!     &       k, adjustl(pft_name)
!        titlevar = trim('cover) '//resout)
!        write( title, "(i2, ' - ', a46, titlevar)" )
!     &       k, adjustl(pft_name)
!        write(80) title, vfc(:,:,k)
         title = titlec(k)(1:48)//"  "//trim(MON)//" (cover)  "//resoutt
         write(80) title, vfc(:,:,k)
      enddo
      close(80)

      open(90,file=fileprefix//"_lai_"//
     &     trim(MISC)//trim(MONstr)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!        read( titlec(k), '(a4, a46)' ) dum, pft_name
!        write( title, '(i2, " - ", a46, " (LAI) 4x5")' )
!     &       k, adjustl(pft_name)
!        write(90) title, laic(:,:,k)
         title = titlec(k)(1:48)//"  "//trim(MON)//" (LAI)  "//resoutt
         write(90) trim(title), laic(:,:,k)
      enddo

      close(90)


      !Write same netcdf output -----------------
      !LC
      fileout = fileprefix//"_lc_"//trim(MISC)//trim(MONstr)//".nc"
      call my_nf_create_Ent(IM,JM
     &     ,trim(fileout)
     &     ,n,ent_names16, Ent_title16,'cover fraction', ncid)
      call my_nf_defglobal(trim(fileout))
      status = my_nf_open(trim(fileout), NF_WRITE,ncid)
      do k=1,n
         !status=nf90_inq_varid(ncid, ent_names(k), varid)
         !status=nf90_inq_dimid(ncid, ent_names(k), dimid)
         !status=nf90_put_var(ncid, varid, vfc(:,:,k), start, count)
         !call handle_nf_error(status, 'nf90_put_var'//ent_names(k))
         status = my_nf_inq_put_var_real32_2(ncid,
     &        trim(ent_names16(k)), varid,vfc(:,:,k))

      end do
      status=nf90_close(ncid)
      
      !LAI
      fileout = fileprefix//"_lai_"//trim(MISC)//trim(MONstr)//".nc"
      call my_nf_create_Ent(IM,JM
     &     ,trim(fileout)
     &     ,n,ent_names16, Ent_title16,"m^2 leaf / m^2 ground", ncid)
      call my_nf_defglobal(trim(fileout))
      status = my_nf_open(trim(fileout), NF_WRITE,ncid)
      do k=1,n
         !status=nf90_inq_varid(ncid, ent_names(k), varid)
         !status=nf90_inq_dimid(ncid, ent_names(k), dimid)
         !status=nf90_put_var(ncid, varid, laic(:,:,k), start, count)
         !call handle_nf_error(status, 'nf90_put_var'//ent_names(k))
         status = my_nf_inq_put_var_real32_2(ncid,
     &        trim(ent_names16(k)), varid,laic(:,:,k))
      end do
      status=nf90_close(ncid)
      
      end subroutine write_output

      
      
      subroutine write_output_h_single(IM,JM
     &     ,titleh, entnamek, h,hsd, filename
     &     , MISC, resoutt)
      use netcdf
      use convertnc
      implicit none
      !Same as write_output_h, but single arrays
      integer,intent(in) :: IM,JM
      character*80 :: titleh(:)
      integer :: entnamek
      real*4 :: h(:,:), hsd(:,:)
      character*(*) :: filename
      character*(*) :: MISC
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      character*256 :: fileout
      integer :: status, ncid, varid
      integer :: lonid, latid, dim(2)
      
      open(90,file=filename//trim(MISC)//".ij",
     &     form="unformatted",status="unknown")

      title = titleh(1)(1:63)//"  "//trim(resoutt)
      write(90) title, h(:,:)

      title = titleh(2)(1:63)//"  "//trim(resoutt)
      write(90) title, hsd(:,:)

      close(90)

      !Netcdf output
      fileout = filename//trim(MISC)//".nc"
      call my_nf_create_Ent_single(IM,JM
     &     ,trim(fileout),"h_"//trim(ent_names16(entnamek))
     &     ,Ent_Title16(entnamek)//" height", "meters", ncid)
      call my_nf_defglobal(trim(fileout))
      status = my_nf_open(trim(fileout), NF_WRITE,ncid)
      status = my_nf_inq_put_var_real32_2(ncid,
     &     "h_"//ent_names16(entnamek), varid,h(:,:))
      
      !Need to make a my_nf_create_Entx2 to output two sets of vars.
 !     status = my_nf_inq_put_var_real32_2(ncid,
 !    &        trim(ent_names(k)), varid,hsd(:,:))
      write(*,*) 'Got here'
      status = nf_inq_dimid(ncid,'lon',lonid)
      write(*,*) status, 'lon', lonid
      status = nf_inq_dimid(ncid,'lat',latid)
      write(*,*) status, 'lat', latid
      status = my_nf_inq_def_put_var_real32_2(ncid,
     &     IM,JM,1,IM,1,JM,lonid,latid,
     &     'stdv_h'//ent_names16(entnamek),
     &     Ent_title16(entnamek)//' standard deviation of height'
     &     , 'meters', hsd(:,:))
      status=nf90_close(ncid)
      
      end subroutine write_output_h_single


      subroutine write_output_h(IM,JM,n,titleh, h,hsd, filename
     &     , MISC, resoutt)
      integer,intent(in) :: IM,JM
      character*80 :: titleh(:,:)
      real*4 :: h(:,:,:), hsd(:,:,:)
      integer :: n
      character*(*) :: filename
      character*(*) :: MISC
      character*(*) :: resoutt !resout for titles
      !---
      character*80 :: title, pft_name, dum
      character*20 :: titlevar
      character*256 :: fileout
      integer :: k
      integer :: status, ncid, varid

      open(90,file=filename//trim(MISC)//".ij",
     &     form="unformatted",status="unknown")

      do k=1,n
!         title = titleh(1,k)(1:48)//" "//trim(MISC)//
!     &        "   height (m) "//resout
         title = titleh(1,k)(1:63)//"  "//trim(resoutt)
         write(90) title, h(:,:,k)
      enddo
      do k=1,n
!         title = titleh(2,k)(1:48)//" "//trim(MISC)//
!     &        "   stdev (m) "//resout
         title = titleh(2,k)(1:63)//"  "//trim(resoutt)
         write(90) title, hsd(:,:,k)
      enddo
      close(90)

      !* Netcdf output
      fileout = filename//trim(MISC)//".nc"
      call my_nf_create_Ent(IM,JM
     &     ,trim(fileout),n
     &     ,ent_names16,Ent_title16,'height (m)', ncid)
      call my_nf_defglobal(trim(fileout))
      status = my_nf_open(trim(fileout), NF_WRITE,ncid)
      do k=1,n
         status = my_nf_inq_put_var_real32_2(ncid,
     &        trim(ent_names16(k)), varid,h(:,:,k))
      end do
      !Need to make a my_nf_create_Entx2 to output two sets of vars.
!      do k=1,n
!         status = my_nf_inq_put_var_real32_2(ncid,
!     &        trim(ent_names(k)), varid,hsd(:,:,k))
!      end do
      status=nf90_close(ncid)
 
      end subroutine write_output_h



      subroutine calculate_bs_brightratio(IMn,JMn, bsf_1,bsf_0,
     o     bs_brightratio)
      ! compute ratio of bright/total bare soil.
      ! this will be used to compute bright and dark cover fractions
      !   so that their sum preserves albedo of the bare soil.
      ! extend it to all cells with no data
      implicit none
      integer, intent(in) :: IMn, JMn
      real*4, intent(in) :: bsf_1(IMn,JMn), bsf_0(IMn,JMn) !1-bright cov, 0-dark cov
      real*4, intent(out) :: bs_brightratio(IMn,JMn) !Fraction of bare that is bright.
      !--------------------
      real*4, parameter :: undef = -1.e30
      character*80 :: title_bs
      real*4 bsf(IMn,JMn)
      integer count, i,j, ii,jj, k
      real*4 :: a, s, s1
      
      bsf(:,:) = undef
      bs_brightratio(:,:) = undef
      do j=1,JMn
         do i=1,IMn
            if ((bsf_1(i,j).eq.undef).or.(bsf_0(i,j).eq.undef)) then
               bsf(i,j) = 0.d0
               bs_brightratio(i,j) = undef
            elseif ((bsf_1(i,j)+bsf_0(i,j)).gt.0.0) then
!            if (bsf_1(i,j)<0.0) then
!               write(*,*) 'bsf_1<0.0', bsf_1(i,j)
!            endif
!            if (bsf_0(i,j)<0.0) then
!               write(*,*) 'bsf_0<0.0', bsf_0(i,j)
!            endif
               bsf(i,j) = bsf_1(i,j) + bsf_0(i,j)
               bs_brightratio(i,j) = bsf_1(i,j) /  bsf(i,j) 
               if (bs_brightratio(i,j).gt.1.0) then
                  write(*,*) 'bs_brightratio>1.0,i,j, bratio'
     &                 , i,j,bs_brightratio(i,j)
               endif
            else
               bsf(i,j) = 0.d0
               bs_brightratio(i,j) = undef
            endif
         enddo
      enddo 

      !do
      count = 0
      do j=1,JMn
         do i=1,IMn
!     if ( (bsf(i,j).gt.0.).and.(bsf.lt.1.))  then !Fractional bare soil, only with VEG
            if  (( bsf(i,j).gt.0.).and.(bsf(i,j).lt.1.) ) then !Smear from averaging adjacent cells, 1 grid cell apart.
               a = 0.
               s = 0.
               s1 = 0.
               do jj=max(1,j-1),min(JMn,j+1)
                  do ii=max(1,i-1),min(IMn,i+1)
                     if (bsf(ii,jj).gt.0.0) then
                        a = a + bs_brightratio(ii,jj)*bsf(ii,jj)
                        s = s + bsf(ii,jj)
                        s1 = s1 + 1.
                     endif
                  enddo
               enddo
               if ( s.gt.0.0 ) then
                  bs_brightratio(i,j) = a/s
                  bsf(i,j) = s/s1
               else             !No adjacent bare soil cells
                 !bs_brightratio(i,j) = undef !Do nothing, should already be undef
                  count = count + 1
               endif
            endif
         enddo
      enddo
       
      print *, "count=", count
        !if ( count== 0 ) exit
      !enddo

      title_bs = "bsf_1"
      write(91) title_bs, bsf_1
      title_bs = "bsf_0"
      write(91) title_bs, bsf_0
      title_bs = "bare soil brightratio"
      write(91) title_bs, bs_brightratio
      title_bs = "ground+ice cover"
      write(91) title_bs, bsf

      end subroutine calculate_bs_brightratio

      subroutine get_bare_soil_brightratio(IMn,JMn, filename
     &     ,bs_brightratio)
      ! read bare soil from the "old"  dataset, in GISS layer format
      ! compute ratio of bright/total bare soil.
      ! this will be used to compute bright and dark cover fractions
      !   so that their sum preserves albedo of the bare soil.
      ! extend it to all cells with no data
      implicit none
      integer, intent(in) :: IMn, JMn
      character*(*) :: filename
      real*4, intent(out) :: bs_brightratio(IMn,JMn) !Fraction of bare that is bright.
      !---
      ! bare soil data from "old" dataset
      real*4 :: bsf_1(IMn,JMn), bsf_0(IMn,JMn) !1-bright cov, 0-dark cov
      character*80 :: title_bs
      integer :: k

      open(1,file=filename,
     &     form="unformatted",status="old")

      read(1) title_bs, bsf_1(:,:) !BRIGHT cover
      write(*,*) 'read BRIGHT soil'
      do k=2,9
        read(1)
        write(*,*) k
      enddo
      read(1) title_bs, bsf_0(:,:) !DARK cover

      close(1)

      call calculate_bs_brightratio(IMn, JMn, bsf_1, bsf_0,
     o     bs_brightratio)

      end subroutine get_bare_soil_brightratio

      
      subroutine get_bare_soil_brightratio_Carrer(IMn,JMn, filename
     &     ,bs_brightratio)
      ! Read soil albedo from netcdf file. This can be with bright+dark fractions = 1
      ! for soil albedo only, or from a VEG data set with actual bright and dark cover fractions.
      ! Compute bs_brightratio =  bright/total bare soil.
      ! This will be used to compute bright and dark cover fractions
      !   so that their sum preserves albedo of the bare soil.
      ! Extend it to all cells with no data
      use netcdf
      use convertnc
      implicit none
      integer, intent(in) :: IMn, JMn
      character*(*) :: filename
      real*4, intent(out) :: bs_brightratio(IMn,JMn) !Fraction of bare that is bright.
      !---
      ! bare soil data from "old" dataset
      real*4 :: bsf_1(IMn,JMn), bsf_0(IMn,JMn) !1-bright cov, 0-dark cov
      character*80 :: title_bs
      integer :: k, status, varid, ncid

      status = my_nf_open(filename,0,ncid)

      status = my_nf_inq_get_var_real32_2(ncid,'bare_bright_grey',
     &     varid,bsf_1)
      status = my_nf_inq_get_var_real32_2(ncid,'bare_dark_grey',
     &     varid,bsf_0)
      status = nf_close(ncid)
      
      call calculate_bs_brightratio(IMn, JMn, bsf_1, bsf_0,
     o     bs_brightratio)

!      title_bs = 'Carrer soil bright ratio'
!      write(990) title_bs, bs_brightratio
      
      end subroutine get_bare_soil_brightratio_Carrer

      

      subroutine split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bsbr,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
      !Split BARE soil into BRIGHT and DARK cover to preserve albedo from
      !  "old" ModelE cover.  Should be called after each trim, scale, nocrops.
      !Any LAI on BARE soil should already have been moved to vegetated cover,
      !  so laic(:,:,N_BARE) should be zero.
      !This checks for cases if BARE is original total or was previously split.
      use EntGVSD_netcdf_util, only : MONTH
      implicit none
      integer, intent(in) :: N_VEG, IMn,JMn,KM
      integer, intent(inout) ::N_BARE
      real*4, intent(in) :: bsbr(:,:) !(IMn,JMn) !Fraction of bare that is bright.
      real*4 vfc(:,:,:), laic(:,:,:) !(IMn,JMn,KM)
      real*4 vfm(:,:,:,:), laim(:,:,:,:) !(12,IMn,JMn,KM)
      real*4 vfh(:,:,:) !(IMn,JMn,KM),
      real*4 hm(:,:,:), hsd(:,:,:) !(IMn,JMn,KM)
      character*80 :: titlec(:) !(18)
      character*80 :: titlem(:,:) !(12,18)
      character*80 :: titleh(:,:) !(2,18) !1-h, 2-hsd
      character*(*) :: res_out
      !-----Local----
      real*4 :: vfc_bare(IMn,JMn), vft !vf_tot(IMn,JMn)
      integer :: i,j,m
      integer :: BR, DK !bright, dark indices


      if ((N_BARE-N_VEG).eq.1) then !First time splitting
      ! increase number of fractions by 1 so that
      ! vfc(:,:,N_BARE-1) = fraction of bright soil
      ! vfc(:,:,N_BARE) = fraction of dark soil
         write(*,*) 'N_BARE from bare_sparse:', N_BARE
         BR=N_BARE
         DK=N_BARE+1
         do i=1,IMn
            do j=1,JMn
               vft = vfc(i,j,N_BARE) !All bare soil orig in N_BARE layer
              !if ((vf_tot(i,j).lt.0.0).and.(vf_tot(i,j).ne.-1.e30)) then
              !   vf_tot(i,j) = 0.0
              !endif
              if ((vft.lt.0.0).or.
     &             ((vft.gt.0.0).and.(bsbr(i,j).lt.0.0))) then !Bad data
                 write(*,*)'ERR: i,j,vft,bsbr',i,j,vft,bsbr(i,j)
                 if (vft.gt.0.) then
                    !Need to split to bright and dark.
                    write(*,*) 'No bsbr, assigning gray 0.2,',i,j 
                    vfc(i,j,BR) = max(0.0,vft)*0.4 !Bright fraction GREY
                    vfc(i,j,DK) = max(0.0,vft)*0.6 !Dark fraction   GREY
                    laic(i,j,BR:DK) = 0.
                 else           !else check why negative vft
                    write(*,*) 'Negative vft', i,j, vft
                    STOP
                 endif
              else              !Good data
                 vfc(i,j,BR) = max(0.0,vft)*bsbr(i,j) !Bright fraction
                 vfc(i,j,DK) = max(0.0,vft)*(1.-bsbr(i,j)) !Dark fraction
                 laic(i,j,BR:DK) = 0.
              endif

              do m=1,12         !Monthly
                 vft = vfm(m,i,j,N_BARE)
                 if ((vft.lt.0.0).or.
     &                ((vft.gt.0.).and.(bsbr(i,j).lt.0.))) then !Bad data
                    write(*,*)'ERRm',N_BARE,m,i,j,vft,bsbr(i,j)
                    if (vft.gt.0.) then
                       !Need to split to bright and dark.
                       write(*,*) 'm:No bsbr, assigning gray 0.2,',i,j,m
                       vfm(m,i,j,BR) = max(0.0,vft)*0.4 !Bright fraction GREY
                       vfm(m,i,j,DK) = max(0.0,vft)*0.6 !Dark fraction   GREY
                       laim(m,i,j,BR:DK) = 0.
                    else        !check why negative vft
                       write(*,*) 'Negative vftm:', m,i,j,vft
                       STOP
                    endif
                 else           !Good data
                    vfm(m,i,j,BR) = max(vft,0.0) * bsbr(i,j)
                    !vfm(m,i,j,N_BARE+1) = max(vf_tot(i,j),0.0) - vfm(m,i,j,N_BARE)
                    vfm(m,i,j,DK) = vft*(1.0-bsbr(i,j))
                    laim(m,i,j,BR:DK) = 0.
                 endif
              enddo              
           end do
        end do

        N_BARE = N_BARE+1       !Index of last BARE layer

      else  !N_BARE-NVEG.eq.2, previously split and incremented N_BARE
         write(*,*) 'N_BARE from bright and dark:',N_BARE
         BR=N_BARE-1
         DK=N_BARE
         do i=1,IMn
            do j=1,JMn
               vft = vfc(i,j,BR) + vfc(i,j,DK)
c               vfc(:,:,N_BARE-1) = vf_tot(:,:)*bsbr(:,:)
c              !vfc(:,:,N_BARE) = vf_tot(:,:) - vfc(:,:,N_BARE-1)
c               vfc(:,:,N_BARE) = vf_tot(:,:)*(1.-bsbr(:,:))
c               laic(:,:,N_BARE-1:N_BARE) = 0.
c               do m=1,12
c                  vf_tot(:,:) = vfm(m,:,:,N_BARE-1) + vfm(m,:,:,N_BARE)
c                  vfm(m,:,:,N_BARE-1) = vf_tot(:,:)*bsbr(:,:)
c                  vfm(m,:,:,N_BARE) = vf_tot(:,:) - vfm(m,:,:,N_BARE-1)
c                  laim(m,:,:,N_BARE-1:N_BARE) = 0.
c                enddo

               if ((vft.lt.0.).or.
     &              ((vft.gt.0).and.(bsbr(i,j).lt.0.))) then !Bad data
                  write(*,*) 'ERR2: bsbr',i,j,vft,bsbr(i,j)
                  !Keep existing bright and dark
               else             !Good data
                  vfc(i,j,BR) = vft*bsbr(i,j)
                  !vfc(i,j,N_BARE) = vf_tot(i,j) - vfc(i,j,N_BARE-1)
                  vfc(i,j,DK) = vft*(1.-bsbr(i,j))
                  laic(i,j,BR:DK) = 0.
               endif
               do m=1,12
                  if ((vft.lt.0.).or.
     &                 ((vft.gt.0.).and.(bsbr(i,j).lt.0.))) then !Bad data
                     write(*,*) 'ERR2m: bsbr',i,j,vft
     &                    ,vfc(i,j,BR),vfc(i,j,DK),bsbr(i,j)
                     !Keep existing bright and dark, or check why vft<0.
                  else          !Good data
                     vfm(m,i,j,BR) = vft*bsbr(i,j)
                     !vfm(m,i,j,N_BARE) = vft - vfm(m,i,j,N_BARE-1)
                     vfm(m,i,j,DK) = vft*(1-bsbr(i,j))
                  endif
                  laim(m,i,j,BR:DK) = 0.
               enddo
            enddo
         enddo
      endif

      vfh(:,:,BR) = vfc(:,:,BR)
      vfh(:,:,DK) = vfc(:,:,DK)
      hm(:,:,N_BARE) = 0.
      hsd(:,:,N_BARE) = 0.

      titlec(BR) =   "17 - bright bare soil  "//
     &     "                           (COVER) "//res_out
      titlec(DK) =     "18 - dark bare soil    "//
     &     "                           (COVER) "//res_out
      do m=1,12
            titlem(m,BR) = "17 - bright bare soil  "//
     &           "                          "//MONTH(m)//" (cover) "
     &           //res_out
            titlem(m,DK) =   "18 - dark bare soil    "//
     &           "                          "//MONTH(m)//" (cover) "
     &           //res_out
         enddo
         
      titleh(1,BR) =   "17 - bright bare soil  "//
     &     "height (m)                                "//res_out
      titleh(1,DK) =     "18 - dark bare soil    "//
     &     "height (m)                                "//res_out
      titleh(2,BR) =   "17 - bright bare soil  "//
     &     "height stdev(m)                           "//res_out
      titleh(2,DK) =     "18 - dark bare soil    "//
     &     "height stdev(m)                           "//res_out

      write(*,*) 'Finished split_bare_soil'
      end subroutine split_bare_soil


      subroutine sum_lai_single(iu, IMn,JMn,KM,Ln, tag, vfs, lais)
      !Sum LAI over all cover for a single point in time. Output to iu.
      use arrayutil, only : IJAdd4, IJMult4, IJDiv4
      implicit none
      integer, intent(in) :: iu
      integer, intent(in) :: IMn,JMn,KM
      integer, intent(in) :: Ln !# layers to loop through
      character*(*) :: tag
      real*4 :: vfs(IMn,JMn,KM), lais(IMn,JMn,KM)
      !---Local---

      integer :: i,j,k,m
      character*80 :: TITLE
      real*4 :: L1(IMn,JMn)

      L1(:,:) = 0.0
      do k=1,Ln
            L1(:,:) = IJADD4(1,IMn,1,JMn,L1(:,:),
     &           IJMULT4(1,IMn,1,JMn,vfs(:,:,k),lais(:,:,k)))
      enddo

      TITLE = trim(tag)//'   LAI'
      write(iu) TITLE, L1(:,:)

      end subroutine sum_lai_single

      subroutine sum_lai(iu, IMn,JMn,KM,Ln, tag, vfmon, laimon)
      !Sum LAI over all cover by month.  Output to iu.
      use arrayutil, only : IJAdd4, IJMult4, IJDiv4
      implicit none
      integer, intent(in) :: iu
      integer, intent(in) :: IMn,JMn,KM,Ln
      character*(*) :: tag
      real*4 :: vfmon(12,IMn,JMn,KM), laimon(12,IMn,JMn,KM)
      !---Local---
      character*3, parameter :: MON(12) =
     &     (/
     &     "Jan","Feb","Mar","Apr","May","Jun",
     &     "Jul","Aug","Sep","Oct","Nov","Dec"
     &     /)

      integer :: i,j,k,m
      character*80 :: TITLE
      real*4 :: L1(12,IMn,JMn)

      L1(:,:,:) = 0.0
      do m=1,12
         do k=1,Ln
            L1(m,:,:) = IJADD4(1,IMn,1,JMn,L1(m,:,:),
     &           IJMULT4(1,IMn,1,JMn,vfmon(m,:,:,k),laimon(m,:,:,k)))
         enddo
      enddo

      do m=1,12
         TITLE = trim(tag)//' '//MON(m)//'   LAI'
         write(iu) TITLE, L1(m,:,:)
      enddo

      
      end subroutine sum_lai

      subroutine check_lc_lai_mismatch(KX,IX,JX,vf,lai,label,titlek)
      implicit none
      integer, intent(in) :: KX,IX,JX
      real*4, intent(in) :: vf(KX,IX,JX),lai(KX,IX,JX)
      character*3,intent(in) :: label
      character*80 :: titlek(:)
      !------
      integer :: i,j,k
      real*4 :: LAYER(IX,JX)
      character*80 :: TITLE

      do k=1,KX
         LAYER(:,:) = 0.
         do i=1,IX
            do j=1,JX
               if (((vf(i,j,k).eq.0.0).and.(lai(i,j,k).gt.0.0))
     &              .or.(vf(i,j,k).gt.0.0).and.(lai(i,j,k).eq.0.0))
     &              then
                  print *,label,' mismatch,k,i,j,cov,lai'
     &                 ,k,i,j,vf(i,j,k),lai(i,j,k)
                  LAYER(i,j) = vf(i,j,k)
               endif
            enddo
         enddo
         TITLE = label//' lc_lai mis '//titlek(k)(1:20)
         write(100) TITLE, LAYER
      enddo
      end subroutine check_lc_lai_mismatch


      subroutine replace_crops(IMn,JMn,KM,i,j !,s
     &     ,bs_brightratio
     &     ,vfc,vfm,laic,laim
     &     ,naturalvegfound)  !,flag)
      !Replace crops in (i,j) with main natural cover in cell or adjacent.
      !The check for existence of crops in (i,j) is done before subroutine call
      !Fix 6 cells that are all crops due to MODIS error or two islands.
      !Assign LAI (cover-avg) of natural cover if from adjacent cells.
      !Do not zero crop LAI but keep for historical cover change.
      !NOTES:  Checks of adjacent cells will use replaced cover from previously
      ! processed cells (previous i,j).  I don't see this as a problem,
      ! since there won't be a consistent dominance either side of the cell.
      ! it helps reduce iteration for adjacent cells.
      implicit none
      integer, intent(in) :: IMn,JMn,KM,i, j
      !real*4 :: s
      real*4, intent(in) :: bs_brightratio(:,:) !(IMn,JMn) !Fraction of bare that is bright.
      !real*4 :: vfct(:,:,:) !vfc trim before nocrops
      !real*4 :: vfmt(:,:,:,:)
      real*4 :: vfc(:,:,:) , laic(:,:,:)
      real*4 :: vfm(:,:,:,:), laim(12,IMn,JMn,KM)
      !integer :: flag
      logical,intent(out) :: naturalvegfound
      !----Local----
      integer :: m, k, ii, jj !cell to search for natural (i,j or adjacent)
      real*4 :: covmax
      integer :: covmaxk, covmaxii,covmaxjj
      real*4 :: covsum(14), covavglai(14)
      real*4 :: covsumm(12,14),covavglaim(12,14)
      integer :: dg
      real*4 :: br, brcov   !bright soil ratio, br, to total bare soil, brcov
      
      if (sum(vfc(i,j,15:16)).le.0.0) then
         return !No crops, no need to do
      endif
      
      !--First check just (i,j) cell for natural veg---
      naturalvegfound = .false.
      covmax = 0.d0
      covmaxk = 0
      covmaxii = i
      covmaxjj = j
      do k=1,14         !Find max non-crop, non-bare natural cover type
         if (vfc(i,j,k) > covmax)  then
            covmax = vfc(i,j,k)
            covmaxk = k
         endif
      enddo
      if (covmax.gt.0.d0) then  !Assign dominant natural veg to crop
         naturalvegfound = .true.
         vfc(i,j,covmaxk) = vfc(i,j,covmaxk) 
     &        + vfc(i,j,15)  + vfc(i,j,16)
         vfc(i,j,15:16) = 0.0 !zero out crop cover - done below
         vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk)
     &        + vfm(:,i,j,15)  + vfm(:,i,j,16)
         vfm(:,i,j,15:16) = 0.0 !zero out crop cover - done below
         !If in same cell, keep original LAI in (i,j) of natural veg
         !!DON'T DO ASSIGNMENT BELOW
         !!laic(i,j,covmaxk) = laic(covmaxii,covmaxjj,covmaxk)
         !!laim(:,i,j,covmaxk) = laim(:,covmaxii,covmaxjj,covmaxk)
         !write(*,*) 'Cell has crops + natural',i,j,s
      else
         write(*,*) "No natural veg in (i,j). Checking adjacent",i,j !,s
         dg = 1
         do while ((.not.naturalvegfound).and.(dg.le.9)) !Check up to 9 grid 
         
            covsum(:) = 0.d0
            covavglai(:) = 0.d0
            covsumm(:,:) = 0.
            covavglaim(:,:) = 0.
            do ii=max(1,i-dg),min(i+dg,IMn)
               do jj=max(1,j-dg),min(j+dg,JMn)
                  if ( (ii.ge.1).and.(ii.le.IMn)
     &                 .and.(jj.ge.1).and.(jj.le.JMn) !in grid range
     &                 .and.((ii.ne.i).or.(jj.ne.j)) ) !not the i,j center cell
     &                 then
                     do k=1,14     
                     !Sum adjacent natural veg cover by type.
                        covsum(k) = covsum(k) + vfc(ii,jj,k)
                        covavglai(k) = covavglai(k) + 
     &                       laic(ii,jj,k)*vfc(ii,jj,k)
                        covsumm(:,k) = covsumm(:,k) + vfm(:,ii,jj,k)
                        covavglaim(:,k) = covavglaim(:,k) +
     &                       laim(:,ii,jj,k)*vfm(:,ii,jj,k)
                     enddo
                  endif
               enddo
            enddo
            covmax = 0.
            covmaxk = 0
            covmaxii = 0
            covmaxjj = 0
            do k=1,14           !Find largest adjacent natural cover
               if (covsum(k)>covmax) then
                  covmax = covsum(k)
                  covmaxk = k
               endif
            enddo
            if (covmax>0.) then
               do m=1,12  !Error check
                  if (covsumm(m,covmaxk)<0.) then
                     write(*,*) 'ERRc ',i,j,ii,jj,m
     &                    ,covmax,covsum(covmaxk), covsumm(m,covmaxk)
     &                    ,covavglai(covmaxk),covavglaim(m,covmaxk)
                     STOP
                  endif
               enddo
               ! Assign adjacent natural cover type and LAI to crop
               naturalvegfound = .true.
               !covavglai(k) = covavglai(k)/covsum(k)  !BUG FOUND
               covavglai(covmaxk) = covavglai(covmaxk)/covsum(covmaxk)
               covavglaim(:,covmaxk) = covavglaim(:,covmaxk)/
     &              covsumm(:,covmaxk)
               vfc(i,j,covmaxk) = vfc(i,j,covmaxk) 
     &              + vfc(i,j,15) + vfc(i,j,16)
               vfc(i,j,15:16) = 0.0 !zero out crop cover - done below
               vfm(:,i,j,covmaxk) = vfm(:,i,j,covmaxk)
     &              + vfm(:,i,j,15) + vfm(:,i,j,16)
               vfm(:,i,j,15:16) = 0.0 !zero out crop cover  - done below
               !Assign LAI in (i,j) from adjacent cell
               laic(i,j,covmaxk) = covavglai(covmaxk)
               laim(:,i,j,covmaxk) = covavglaim(:,covmaxk)
               !write(*,*) 'Cell or adjacent has crops + natural',i,j,s
               write(*,*) 'Found natural veg in adjacent cells'
     &              ,i,j,covmaxk,dg
               do m=1,12 !Error check
                  if (covavglaim(m,covmaxk)>10.) then
                     write(*,*) 'ERRc2 bad avg lai',i,j,ii,jj,m,covmaxk
     &                    ,covsumm(m,covmaxk),covavglaim(m,covmaxk)
                     STOP
                  endif
               enddo
            else
               write(*,*) 'No natural veg in adjacent cells'
               dg=dg+1
               write(*,*) 'Increasing dg,',i,j, dg
            endif
         end do                 !while
      endif !Checking for adjacent natural veg
      
      !!!** PRESCRIBED FIXED BY GRID CELL FOR 144X90 ONLY ***!!!
      !For IMn,JMn = 144x90
      !Fix remaining all-crop cells.  These are:  
      !- MODIS error in Antarctic, 
      !     (i,j) = (34,9),(41,9),(48,9),(34,10) -> C3 arctic grass
      !- Islands:  
      !  Mauritius (98,36) - sugar cane,tea,pasture,forest,savanna -> C4 grass 
      !  Nauru (139,45) - grassland bordered by tropicalforest -> C4 grass
      !Antarctic MODIS error
c      if ((IMn.eq.144).and.(JMn.eq.90)) then
c         if ( (((i.eq.34).or.(i.eq.41).or.(i.eq.48)).and.(j.eq.9))
c     &        .or.( (i.eq.34).and.(j.eq.10) ) ) then !Antartic
c            vfc(i,j,14) = vfc(i,j,14) + vfc(i,j,15) + vfc(i,j,16)
c            vfm(:,i,j,14) = vfm(:,i,j,14) 
c     &           + vfm(:,i,j,15) + vfm(:,i,j,16)
c            laic(i,j,14) = (laic(i,j,15)*vfc(i,j,15) 
c     &           + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
c            do m=1,12
c               laim(m,i,j,14) = ( laim(m,i,j,15)*vfm(m,i,j,15) +
c     &              laim(m,i,j,16)*vfm(m,i,j,16) ) / 
c     &              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
c            enddo
c            write(*,*) 'Replaced Antarctic crops.'
c            naturalvegfound=.true.
c         elseif ( ((i.eq.98).and.(j.eq.36))
c     &           .or.((i.eq.139).and.(j.eq.45)) ) then !tropical islands
c            vfc(i,j,12) = vfc(i,j,12) + vfc(i,j,15) + vfc(i,j,16)
c            laic(i,j,12) = (laic(i,j,15)*vfc(i,j,15) 
c     &        + laic(i,j,16)*vfc(i,j,16)) / (vfc(i,j,15)+vfc(i,j,16))
c            do m=1,12
c               laim(m,i,j,12) = ( laim(m,i,j,15)*vfm(m,i,j,15) +
c     &              laim(m,i,j,16)*vfm(m,i,j,16)) / 
c     &              ( vfm(m,i,j,15)+vfm(m,i,j,16) )
c            enddo
c            write(*,*) 'Replaced island crops.'
c            naturalvegfound=.true.
c         endif
c      endif

      !*** PRESCRIBED GRID CELL FIXES FOR 720X360 *********
      if (.not.naturalvegfound) then !Assign bare
         if (bs_brightratio(i,j).eq.-1.e30) then
            !Find bs_brightratio in adjacent cells
            br = 0.
            brcov = 0.
            dg=0
            do while ((brcov.eq.0.).and.(dg.le.5)) !arbitrary neighbor range
               dg=dg+1
               do ii=max(1,i-dg),min(i+dg,IMn)
                  do jj=max(1,j-dg),min(j-dg,JMn)
                     if (bs_brightratio(ii,jj).gt.0.) then
                        br = br + bs_brightratio(ii,jj)
                        brcov = brcov + sum(vfc(ii,jj,17:18))
                     endif
                  enddo
               enddo
            enddo
            if (brcov.gt.0.) then
               br = br/brcov    !cover-weighted average br
               write(*,*) "Adjacent bare soil found",i,j,dg, br
            else
               write(*,*) "Adjacent bare soil NOT found",i,j,dg, br
               if (.FALSE.) then !Done only to diagnose cells w/o bratio
                  write(*,*) "assigning gray 0.2"
                  br=0.4        !This gives 0.5*.4 + 0.*0. = 0.2 albedo
               else
                  write(*,*) "Assign specific grid cells"
                  br=-1.e30     !Assign undef to make sure no wrong values.
               endif
            endif
         else !bs_brightratio was defined
            br = bs_brightratio(i,j)
         endif

         !Assign bare soil fractions
         if (br.ne.-1.e30) then  !Adjacent bs_brightratio found
            vfc(i,j,17) = vfc(i,j,17) +
     &           sum(vfc(i,j,15:16))*br
            vfc(i,j,18) = vfc(i,j,18) +
     &           sum(vfc(i,j,15:16))*(1. - br)
            vfc(i,j,15:16) = 0. !zero out crops
            vfm(:,i,j,17) = vfm(:,i,j,17) +
     &           sum(vfm(:,i,j,15:16))*br
            vfm(:,i,j,18) = vfm(:,i,j,18) +
     &           sum(vfm(:,i,j,15:16))*(1. - br)
            vfm(:,i,j,15:16) = 0. !zero out crops
            write(*,*) 'Assigned nearby bare'
     &           ,i,j, dg, bs_brightratio(i,j),br
         elseif ((IMn.eq.720).and.(JMn.eq.360)) then
      !Adjacent bs_brightratio not found for only these cells
!     Replace specific grid cells ONLY VALID FOR IM=720,JM=360!!!
!     NOTE:  "Mauritius" and "Fiji/Amer.Samoa" might be off-grid and land mask.
!     assigning bare Found for IM=720,JM=360 = (lon, lat)
!     Below:  i,j,bs_brightratio(i,j), br = (lon, lat)
!     Antarctica 1: (pure17 has C3 crop, partial coastal cell/ice)
!     169  36  -1.00000002E+30  0.400000006 = (-95.75,-72.25)
!     168  37  -1.00000002E+30  0.400000006 = (-96.25,-71.75)
!     Antarctica 2: (pure17 has C3 crop, partial coastal cell/ice)
!     239  36  -1.00000002E+30  0.400000006 = (-60.75,-72.75)
!     204  36  -1.00000002E+30  0.400000006 = (-78.25,-72.25)
!     Mauritius (too east of):   (pure17 was everbroad, C4 crop, sparse, coast)
!     487 141 0.461078912 0.461078912  = (63.25, -19.75)
!     Fiji?  American Samoa? (pure17 was everbroad, C4 crop, sparse, coast)
!     694 179 0.353166729 0.353166729  = (163.25,-0.75)
            if (((i.eq.169).and.(j.eq.36)) !Antarctica1
     &           .or.((i.eq.168).and.(j.eq.37)) !Antarctica1
     &           .or.((i.eq.239).and.(j.eq.36)) !Antarctica2
     &           .or.((i.eq.204).and.(j.eq.36))) !Antarctica2
     &           then
                  !* Assign Antarctica cells to arctic grass
               vfc(i,j,14) = vfc(i,j,14) + sum(vfc(i,j,15:16))
               vfc(i,j,15:16) = 0. !zero out crops
               laic(i,j,14) = max(laic(i,j,14),
     &              max(laic(i,j,15),laic(i,j,16)))
               do m=1,12
                  vfm(m,i,j,14)=vfm(m,i,j,14) + sum(vfm(m,i,j,15:16))
                  vfm(m,i,j,15:16) = 0. !zero out crops
                  laim(m,i,j,14)=max(laim(m,i,j,14),
     &                 max(laim(m,i,j,15),laim(m,i,j,16)))
               enddo
               write(*,*) 'Assigned c3 crops cells in Antartica',i,j
            endif
            if (((i.eq.487).and.(j.eq.141)) !Mauritius
     &           .or.((i.eq.694).and.(j.eq.179))) !Fiji or American Samoa
     &           then
               !* Assign tropical rainforest
               vfc(i,j,2) = vfc(i,j,2) + sum(vfc(i,j,15:16))
               vfc(i,j,15:16) = 0. !zero out crops
               laic(i,j,2) = max(laic(i,j,2),
     &              max(laic(i,j,15), laic(i,j,16)))
               do m=1,12
                  vfm(m,i,j,2)=vfm(m,i,j,2) + sum(vfm(m,i,j,15:16))
                  vfm(m,i,j,15:16) = 0. !zero out crops
                  laim(m,i,j,2)=max(laim(m,i,j,2),
     &                 max(laim(m,i,j,15),laim(m,i,j,16)))
               enddo
               write(*,*) 'Assigned c4 crops cells in islands',i,j
            endif
         else
            !Don't do any grid cell fixes
            write(*,*) 'WRONG GRID RES:  PRESCRIBED FIX NOT DONE: '
            write(*,*) '**Check array indices for grid res for nocrops'            endif !Correct grid resolution
         endif                  !Adjacent bs_brightratio not found
      endif                     !.not.(naturalvegfound)
      
      !Zero out crop cover - done above for each vfc assignment.
!      vfc(i,j,15:16) = 0.
!      vfm(:,i,j,15:16) = 0.
         
      end subroutine replace_crops


      subroutine fill_crops(IMn,JMn,vfc15,laic15
     i     ,vfm15,laim15,hm15,hsd15
     o     ,laiccrop, laimcrop, hmcrop, hsdcrop)
      !This performs in-fill once for herb crop (PFT15) LAI to create
      !an extended crop LAI data set for use with historical crop cover.
      integer :: IMn, JMn
      real*4 :: vfc15(:,:), laic15(:,:)  !(i,j)
      real*4 :: vfm15(:,:,:), laim15(:,:,:) !(m,i,j)
      real*4 :: hm15(:,:), hsd15(:,:) !(i,j)
      real*4 :: laiccrop(:,:), laimcrop(:,:,:) !OUTPUT, max monthly
      real*4 :: hmcrop(:,:), hsdcrop(:,:)
      !--Local----
      integer :: i, j, m, ii, jj, nocropcells, dg
      real*4 :: covsum15, laiavg15, covmsum15(12), laimavg15(12)
      real*4 :: hmavg15, hsdavg15
      character*80 :: titlefoo

      titlefoo = 'vfc15 in - crop fill'
      write(93) titlefoo, vfc15
      titlefoo = 'laic15 in - crop fill'
      write(93) titlefoo, laic15

      nocropcells=0
      do i=1,IMn
         do j=1,JMn
            if (vfc15(i,j).gt.0.d0) then !crops in cell, replicate
               laiccrop(i,j) = laic15(i,j)
               laimcrop(:,i,j) = laim15(:,i,j)
               hmcrop(i,j) = hm15(i,j)
               hsdcrop(i,j) = hsd15(i,j)
            else !(vfc15(i,j).eq.0.d0) then !no crops in cell, fill in LAI15
               covsum15 = 0.d0
               laiavg15 = 0.d0
               hmavg15 = 0.d0
               hsdavg15 = 0.d0
               covmsum15(:) = 0.d0
               laimavg15(:) = 0.d0
!               do ii=i-1,i+1  !ext 1 grid cell
!                  do jj=j-1,j+1
               dg=5
               do ii=i-dg,i+dg    !ext 5 grid cells = 2.5 degrees at HXH
                  do jj=j-dg,j+dg
                     if ( (ii.ge.1).and.(ii.le.IMn)
     &                    .and.(jj.ge.1).and.(jj.le.JMn) !in grid range
     &                    .and.((ii.ne.i).or.(jj.ne.j)) ) !not in i,j 
     &                    then
                        if (vfc15(ii,jj).gt.0.d0) then
                           covsum15 = covsum15 + vfc15(ii,jj)
                           laiavg15 = laiavg15 
     &                          + laic15(ii,jj)*vfc15(ii,jj)
                           hmavg15 = hmavg15 + hm15(ii,jj)*vfc15(ii,jj)
                           hsdavg15 = hsdavg15 
     &                          + hsd15(ii,jj)*vfc15(ii,jj)
                           covmsum15(:) = covmsum15(:) + vfm15(:,ii,jj)
                           laimavg15(:) = laimavg15(:) +
     &                          laim15(:,ii,jj)*vfm15(:,ii,jj)
                        endif
                     endif
                  enddo
               enddo
               if (covsum15>0.d0) then !adjacent crops found, assign laiavg
                  laiavg15 = laiavg15/covsum15
                  hmavg15 = hmavg15/covsum15
                  hsdavg15 = hsdavg15/covsum15
                  laimavg15(:) = laimavg15(:)/covmsum15(:)
                  laiccrop(i,j) = laiavg15
                  hmcrop(i,j) = hmavg15
                  hsdcrop(i,j) = hsdavg15
                  laimcrop(:,i,j) = laimavg15(:)
               else
                  nocropcells = nocropcells + 1
               endif
            endif
         enddo
      enddo

      write(*,*) 'Crop fill-in: no-crop cells:',nocropcells

      end subroutine fill_crops

      end module conversions
!------------------------------------------------------------------------

      program convert
      use conversions
      use netcdf
      use convertnc
      use EntGVSD_netcdf_util

      implicit none

      real*4, parameter :: undef = -1.e30

      ! define input file resolution
      !* .5x.5
      integer, parameter :: IM=720, JM=360, KM=N_COVERTYPES17
!      character*(*), parameter :: res_in="05x05"  !degrees
      character*(*), parameter :: res_in="HXH"  !degrees
      character*(*), parameter :: res_in_int="720x360"  !number of grid cells
      character*(*), parameter :: old_veg_file=
     &     "/Users/nkiang/NancyResearch/GISS/Models/Ent/Code/cmrun/"//
     &     "V720x360_no_crops_downsc.ext"
      character*(*), parameter :: soilalb_file = 
     &     "/Users/nkiang/NancyResearch/GISS/Models/Ent/Datasets/"//
!     &     "MODIS/Schaaf/"//
     &     "Soil/Carrer/CarrerGISS_soil_albedo/"//
     &     "CarrerGISS_soil_albedo_multiband_annual_2004_v1.0b_HXH"//
     &     "_fringeice.nc"
      !character*(*), parameter :: fversion_out='_downscg5' !TEMP FOR MAX KELLEY
      character*(*), parameter :: version_in = 'v1.0b'
      character*(*), parameter :: fversion_out='v1.0b'
      
      !* 1x1
      !integer, parameter :: IM=360, JM=180, KM=20
      !character*(*), parameter :: res_in="1x1"
      !character*(*), parameter :: res_in_int="360x180"
      !* 2.5x2
!      integer, parameter :: IM=144, JM=90, KM=20
!      character*(*), parameter :: res_in="2.5x2"
!      character*(*), parameter :: res_in="144x90"
!      character*(*), parameter :: res_in_int="144x90"
!      character*(*), parameter :: old_veg_file=
!     &     "/Users/nkiang/NancyResearch/GISS/Models/Ent/Code/cmrun/"//
!     &     "V144x90_no_crops.ext"
      !* 4x5
!      integer, parameter :: IM=72, JM=46, KM=20
!      character*(*), parameter :: res_in="4x5"
!      character*(*), parameter :: res_in_int="72x46"
!      character*(*), parameter :: old_veg_file=
!     &     "/Users/nkiang/NancyResearch/GISS/Models/Ent/Code/cmrun/"//
!     &     "V72x46.1.cor2_no_crops.ext"

      !Igor's run:   
      !character*(*), parameter :: old_veg_file="V72x46_no_crops.ext"
      !character*(*), parameter :: old_veg_file="V144X90_no_crops.ext"


      integer, parameter :: N_HT = KM !19 ! number of height layers input

!      real*4 vf(IM,JM,KM), lai(IM,JM,KM)
      character*80 :: title(KM),title_tmp,title12(12,KM),titlehn(2,KM)

      ! Input values
      ! new (interpolated) values  ## NO INTERPOLATION IN THIS PROGRAM ##
      integer, parameter :: IMn=IM, JMn=JM
      character*(*), parameter :: res_out=res_in
      character*(*), parameter :: res_out_int=res_in_int   
      character*(*), parameter :: fversion=fversion_out  !Misc file version

      !character*6, parameter :: fileprein = 'EntMM_'
      character*12, parameter :: fileprein = 'EntGVSD17_MM'
      character*12, parameter :: filepreout = 'EntGVSD16_MM'
      character*(*), parameter :: file_checksum = 
     &     'lc_lai_ent16/'//filepreout//'_checksum_'//res_out_int//'.ij'
      character*256 :: filein_h
      
      ! Input values, max, monthly
      real*4 vfn(IMn,JMn,KM), lain(IMn,JMn,KM), area(IMn,JMn), a
      real*4 vfnm(12,IMn,JMn,KM),lainm(12,IMn,JMn,KM),aream(IMn,JMn),am
      real*4 hmn(IMn,JMn,KM),hsdn(IMn,JMn,KM)
      ! Converted values
      real*4 vfc(IMn,JMn,KM), laic(IMn,JMn,KM), dvf, s
      character*80 :: titlec(18)
      ! Converted values - monthly
      real*4 vfm(12,IMn,JMn,KM), laim(12,IMn,JMn,KM)
      character*80 :: titlem(12,18)
      ! Converted values - heights
      real*4 vfh(IMn,JMn,KM),hm(IMn,JMn,KM), hsd(IMn,JMn,KM)
      character*80 :: titleh(2,18) !1-h, 2-hsd
      ! Converted values - crop ext 
      real*4 laicropext(IMn,JMn), laimcropext(12,IMn,JMn)
      real*4 hmcropext(IMn,JMn),hsdcropext(IMn,JMn)
      ! Vars for calculating nocrops
      logical :: naturalfound
      integer :: flag, nonaturalcount      !if no natural veg in vicinity

      real*4 vf_xx(IMn,JMn), lai_xx(IMn,JMn)
      real*4 vf_yy(IMn,JMn), lai_yy(IMn,JMn)
      real*4 LAYER(IMn,JMn)
      character*80 :: title_xx="xx"
      character*80 :: title_yy="yy"
      character*80 :: titlefoo

      ! bs_brightratio = bare soil brightratio
      real*4 :: bs_brightratio(IMn,JMn), vfc_tmp(IMn,JMn)

      integer i, j, k, io, in, jn, maxpft, kx, m
      real*8 lat
      real*4 foo(IMn,JMn)
      integer, parameter :: delta_j = JM/JMn
      integer, parameter :: delta_i2 = (IM*2)/IMn
      real*8,  parameter :: delta_lat = 180.d0/JM
      integer N_VEG             ! number of PFTs in output
      integer N_BARE            ! index of bare soil in output
      integer count
      integer :: c3, c4, cropsw, bare   ! index to 17PFT + 3 nonveg file

      ! Netcdf vars
      integer :: ncidout, status
      integer :: ncih
      character*120 :: filenc
      
      ! lc & laimax file of Ent 17 PFTs
      !open(1,file="lc_lai_ent/EntMM_lc_laimax_"//res_in//
      open(1,file='lc_lai_ent/'//fileprein//'_lc_laimax_'//
     &     trim(res_in)//"_"//version_in//'.bin',
     &     form="unformatted",status="old")

      do k=1,KM 
        read(1) title(k), vfn(:,:,k)
        write(*,*) 'title vfn, ', title(k)
      enddo
      do k=1,KM
        read(1) title_tmp, lain(:,:,k)
        write(*,*) 'title lain, ', title_tmp !title(k)
      enddo
      close(1)
      !Check if mismatch lc or lai values (one is zero and the other not)
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfn,lain,'vfn',title)


      ! monthly lc and lai files
      do m=1,12
         open(2,file='lc_lai_ent/'//fileprein//'_lc_lai_'//MONTH(m)
     &        //"_"//res_in//"_"//version_in//".bin",
     &        form="unformatted",status="old")
         do k=1,KM
            read(2) title12(m,k), vfnm(m,:,:,k)
            write(*,*) 'title12 vfnm, ', title12(m,k)
         enddo
         do k=1,KM
            read(2) title_tmp, lainm(m,:,:,k)
            write(*,*) 'title12 lainm, ', title_tmp !title12(m,k)
         enddo
         close(2)
      enddo

      ! height file - insert dummy WATER layer at beginning toi avoid
      !               confusion in numbering with vfc and vfm.
      filein_h = "/Users/nkiang/NancyResearch/GISS/Models/Ent/"//
!!     &     "Vegcover/Simard/Conversions/EntGVSDmosaic_height_"
!!     &     "Datasets/MODIS/Schaaf/EntGVSD/VEG_EntGVSDmosaic_heights_"
!     &     "Datasets/MODIS/Schaaf/EntGVSD_05x05/EntGVSDmosaic17_height_"
!     &     //res_in//".ij",
     &     'Vegcover/Simard/Conversions/EntGVSD17_Simard_height_'//
     &     trim(version_in)//'/'//
     &     'EntGVSD17_Simard'//'_height_'//
     &     res_in//'_'//trim(version_in)//".ij"
      write(*,*) filein_h
      open(3,file=trim(filein_h),
     &     form="unformatted",status="old")
      !titlehn(1,1) = 'NO WATER LAYER' !Old - now water layer is last

      !* Netcdf version of Simard heights
!      filein_h = trim(PATHEnt)//trim(PATHfile)//
!     &     'EntGVSD17_height_'//trim(RESOUT)//'_'//
!     &     trim(version)//'.nc'
!      status = my_nf_open(trim(filein_h), ncih)
      
      hmn(:,:,:) = 0.
      do k=1,N_HT !height
         !read(3) titlehn(1,k+1), hmn(:,:,k+1)
         read(3) titlehn(1,k), hmn(:,:,k)
         write(*,*) titlehn(1,k)
         write(92) titlehn(1,k), hmn(:,:,k)
      enddo

      hsdn(:,:,:) = 0.
      do k=1,N_HT !stdev
         !read(3) titlehn(2,k+1), hsdn(:,:,k+1)
         read(3) titlehn(2,k), hsdn(:,:,k)
         write(*,*) titlehn(2,k)
      enddo
      close(3)
      
!      call get_bare_soil_brightratio(IMn,JMn, old_veg_file
!     &     ,bs_brightratio)

      call get_bare_soil_brightratio_Carrer(IMn,JMn, soilalb_file
     &     ,bs_brightratio)

      open(100,file=file_checksum,form="unformatted",status="unknown")


      ! Use Gary Russell's HNTR4.




      !* Convert to GISS 16 pfts format


      ! first 14 pfts though grass are the same, ignore WATER
      !  lc laimax
      vfc(:,:,:) = 0.           !zero
      laic(:,:,:) = 0.
      vfc(:,:,1:14) = vfn(:,:,1:14)
      laic(:,:,1:14) = lain(:,:,1:14)
      titlec(1:14) = title(1:14)
      !  lc lai monthly
      vfm(:,:,:,:) = 0.
      laim(:,:,:,:) = 0.
      vfm(:,:,:,1:14) = vfnm(:,:,:,1:14)
      laim(:,:,:,1:14) = lainm(:,:,:,1:14)
      do m=1,12
         titlem(m,1:14) = title12(m,1:14)
         write(*,*) 'titlem,',m,titlem(m,1:14)
      enddo
      !  heights
      vfh(:,:,1:14) = vfn(:,:,1:14) !Should be the same cover from MODIS.
      hm(:,:,1:14) = hmn(:,:,1:14) 
      hsd(:,:,1:14) = hsdn(:,:,1:14)
      titleh(1,1:14) = titlehn(1,1:14)
      titleh(2,1:14) = titlehn(2,1:14)

      ! crops
#ifdef COMBINE_CROPS_C3_C4
      c3=15 !vfn layers
      c4=16 !vfn
      cropsw=17 !vfn
      do j=1,JMn
         do i=1,IMn
            !lc laimax
            a = vfn(i,j,c3) + vfn(i,j,c4)
            if ( a > 0. ) then
               laic(i,j,15) = (vfn(i,j,c3)*lain(i,j,c3)
     &              + vfn(i,j,c4)*lain(i,j,c4)) / a
               vfc(i,j,15) = a
            else
               laic(i,j,15) = 0.
               vfc(i,j,15) = 0.
            endif
            !lc lai monthly
            do m=1,12
               am = vfnm(m,i,j,c3) + vfnm(m,i,j,c4)
               if (am > 0. ) then
                  laim(m,i,j,15) = (vfnm(m,i,j,c3)*lainm(m,i,j,c3)
     &                 + vfnm(m,i,j,c4)*lainm(m,i,j,c4)) / am
                  vfm(m,i,j,15) = am
               else
                  laim(m,i,j,15) = 0.
                  vfm(m,i,j,15) = 0.
               endif
            enddo
            !heights - DO NOT AVERAGE. PRESERVE HEIGHTS. LAI will scale density
            a = vfn(i,j,c3) + vfn(i,j,c4)  !input cover
            if ( a > 0. ) then
!               hm(i,j,15) = (vfn(i,j,c3)*hmn(i,j,c3)
!     &              + vfn(i,j,c4)*hmn(i,j,c4)) / a
               if ((hmn(i,j,c3)>0.).and.(hmn(i,j,c4)>0.)) then
                  !average if both exist
                  hm(i,j,15) = (vfn(i,j,c3)*hmn(i,j,c3)
     &                 + vfn(i,j,c4)*hmn(i,j,c4)) / a
               else
                  !don't average if only one or none exists
                  hm(i,j,15) = max(hmn(i,j,c3),hmn(i,j,c4))
               endif
               !Sum of squares for sd.  Don't weight if only one or less exists
               if ((hmn(i,j,c3)>0.).and.(hmn(i,j,c4)>0.)) then
                  hsd(i,j,15) = sqrt((vfn(i,j,c3)*hsdn(i,j,c3)**2
     &                 + vfn(i,j,c4)*hsdn(i,j,c4)**2) / a)
               else
                  hsd(i,j,15) = max(hsd(i,j,c3),hsd(i,j,c4))
               endif
               vfh(i,j,15) = a
            else
               hm(i,j,15) = 0.
               hsd(i,j,15) = 0.
               vfh(i,j,15) = 0.
            endif
         enddo                  !j=
      enddo                     !i=
      write(*,*) "Re-doing crops.."
      titlec(15) = title(c3)
      titlec(15)(1:18) = "15 - crops herb   "
      do m=1,12
         titlem(m,15) = title12(m,15)
         titlem(m,15)(1:18) = "15 - crops herb   "
      enddo
      titleh(1,15) = titlehn(1,15)
      titleh(2,15) = titlehn(2,15)
      titleh(1,15)(1:18) = "15 - crops herb   "
      titleh(2,15)(1:18) = "15 - crops herb   "

      write(*,*) titlec(15)
      ! crops woody
      vfc(:,:,16) = vfn(:,:,cropsw)
      laic(:,:,16) = lain(:,:,cropsw)
      titlec(16) = '16 - '//title(cropsw)(6:80)
      write(*,*) "titlec: "
      write(*,*) titlec(16)
      do m=1,12
         vfm(m,:,:,16) = vfnm(m,:,:,cropsw)
         laim(m,:,:,16) = lainm(m,:,:,cropsw)
         titlem(m,16) = '16 - '//title12(m,cropsw)(6:80)
      enddo
      vfh(:,:,16) = vfn(:,:,cropsw)
      titleh(1,16) = '16 - '//titlehn(1,cropsw)(6:80)
      titleh(2,16) = '16 - '//titlehn(2,cropsw)(6:80)
      ! bare soil
      bare=19
      vfc(:,:,17) = vfn(:,:,bare)
      laic(:,:,17) = lain(:,:,bare)
      titlec(17) = title(bare)
      vfm(:,:,:,17) = vfnm(:,:,:,bare)
      laim(:,:,:,17) = lainm(:,:,:,bare)
      do m=1,12
         titlem(m,17) = title12(m,bare)
      enddo
      vfh(:,:,17) = vfn(:,:,bare)
      hm(:,:,17) = hmn(:,:,bare)
      hsd(:,:,17) = hsdn(:,:,bare)
      titleh(1,17) = titlehn(1,bare)
      titleh(2,17) = titlehn(2,bare)
      N_VEG = 16
      N_BARE = 17

      titlefoo = 'vfnc3'
      write(92) titlefoo, vfn(:,:,c3)
      titlefoo = 'vfnc4'
      write(92) titlefoo, vfn(:,:,c4)
      titlefoo = 'Crops 15 after combining C3 and C4'
      write(92) titlefoo, vfc(:,:,15)

#else
      !crops
      vfc(:,:,15:17) = vfn(:,:,15:17)
      laic(:,:,15:17) = lain(:,:,15:17)
      titlec(15:17) = title(15:17)
      do m=1,12
         vfm(m,:,:,15:17) = vfnm(m,:,:,15:17)
         laim(m,:,:,15:17) = lainm(m,:,:,15:17)
         titlem(m,15:17) = title12(m,15:17)
      enddo
      vfh(:,:,15:17) = vfn(:,:,15:17)
      hm(:,:,15:17) = hmn(:,:,15:17)
      titleh(1,15:17) = titlehn(1,15:17)
      titleh(2,15:17) = titlehn(2,15:17)
      ! bare soil
      vfc(:,:,18) = vfn(:,:,bare)
      laic(:,:,18) = lain(:,:,bare)
      titlec(18) = title(bare)
      vfm(:,:,:,18) = vfnm(:,:,:,bare)
      laim(:,:,:,18) = lainm(:,:,:,bare)
      titlem(:,18) = title12(:,bare)
      vfh(:,:,18) = vfh(:,:,bare)
      hm(:,:,18) = hmn(:,:,bare)
      hsd(:,:,18) = hsd(:,:,bare)
      titleh(1,18) = titlehn(1,bare)
      titleh(2,18) = titlehn(2,bare)
      N_VEG = 17
      N_BARE = 18
#endif
      
      do k=1,N_BARE
        write(3) titlec(k), vfc(:,:,k)
      enddo
      do k=1,N_BARE
        write(3) titlec(k), laic(:,:,k)
      enddo
     
      ! check if "bare" soil is not bare
      vf_xx(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .01 .and. laic(i,j,N_BARE) > .5 ) then
            vf_xx(i,j) = vfc(i,j,N_BARE)
            lai_xx(i,j) = laic(i,j,N_BARE)
          endif
        enddo
      enddo

      write(4) title_xx, vf_xx(:,:)
      write(4) title_xx, lai_xx(:,:)

      vf_yy(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,10) > .1 .and. laic(i,j,10) < .5 ) then
            vf_yy(i,j) = vfc(i,j,10)
            lai_yy(i,j) = laic(i,j,10)
          endif
        enddo
      enddo

      write(4) title_yy, vf_yy(:,:)
      write(4) title_yy, lai_yy(:,:)

      vf_yy(:,:) = 0.
      lai_yy(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .1 .and. laic(i,j,10) < .01 
     &      .and. laic(i,j,9) < .01 .and. laic(i,j,11) < .01 
     &      .and. laic(i,j,12) < .01 .and. laic(i,j,13) < .01 ) then
            vf_yy(i,j) = vfc(i,j,N_BARE)
            lai_yy(i,j) = laic(i,j,N_BARE)
          endif
        enddo
      enddo

      write(4) title_yy, vf_yy(:,:)
      write(4) title_yy, lai_yy(:,:)


      !!!! do conversions !!!!

      ! convert sparse veg to cold adapted shrub 9 if present
      do j=1,JMn
        do i=1,IMn

          s = sum(vfc(i,j,1:N_BARE))
          if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
             write(*,*) 'ERROR orig:  max and monthly lc different'
     &            ,i,j ,s,sum(vfm(1,i,j,1:N_BARE))
             write(*,*) vfc(i,j,1:N_BARE)
             write(*,*) vfm(1,i,j,1:N_BARE)
          endif

          if( vfc(i,j,N_BARE) > .0 .and. vfc(i,j,N_BARE) < .15
     &         .and. laic(i,j,N_BARE) > .0
     &         .and. vfc(i,j,9) > .0 ) then
            !convert_vf(vf1, lai1, vf2, lai2, laimin)
            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,9), laic(i,j,9), laic(i,j,9) )
                                ! lai >= lai(9)
            do m=1,12
               !subroutine convert_vfm(vf1, lai1, vf2, lai2, vfc)
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
!!     &              vfm(m,i,j,9),vfm(m,i,j,9), vfc(i,j,9))!BUG??
     &              vfm(m,i,j,9),laim(m,i,j,9), vfc(i,j,9))
             
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,9),hm(i,j,9),hsd(i,j,9), vfc(i,j,9))
          endif

          s = sum(vfc(i,j,1:N_BARE))
          if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
             write(*,*) 'ERROR sparse:  max and monthly lc different'
     &            ,i,j, s,sum(vfm(1,i,j,1:N_BARE))
             write(*,*) vfc(i,j,1:N_BARE)
             write(*,*) vfm(1,i,j,1:N_BARE)
          endif
        enddo
      enddo

      ! convert sparse veg to arid adapted shrub 10 if present
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .0 .and. vfc(i,j,N_BARE) < .15
     &         .and. laic(i,j,N_BARE) > .0
     &         .and. vfc(i,j,10) > .0 ) then

            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,10), laic(i,j,10), laic(i,j,10) )
                                ! lai >= lai(10)
            do m=1,12
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &              vfm(m,i,j,10),laim(m,i,j,10),vfc(i,j,10))
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,10),hm(i,j,10),hsd(i,j,10), vfc(i,j,10))
          endif
        enddo
      enddo

      ! convert the rest of sparse veg to crop 15 if present
      do j=1,JMn
        do i=1,IMn
           if( vfc(i,j,N_BARE) > .0 .and. laic(i,j,N_BARE) > .0
     &         .and. vfc(i,j,15) > .0 ) then
!             print *, 'Converting spare to crop/bare',i,j,
!     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), vfc(i,j,15)
             call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &            vfc(i,j,15), laic(i,j,15), laic(i,j,15))
             if (vfc(i,j,N_BARE)<0.0)  print *, 
     &           'After conversion:            ',i,j,
     &            vfc(i,j,N_BARE), laic(i,j,N_BARE), 
     &            vfc(i,j,15), laic(i,j,15)
             do m=1,12
                call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &               vfm(m,i,j,15),laim(m,i,j,15),vfc(i,j,15))
             enddo
             call convert_vfh(
     &            vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &            vfh(i,j,15),hm(i,j,15),hsd(i,j,15), vfc(i,j,15))
          endif
        enddo
      enddo

      ! convert the rest of sparse veg to pft with biggest fraction
      ! (if present)
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) > .0 .and. laic(i,j,N_BARE) > .0 ) then

            maxpft = maxloc( vfc(i,j,1:16), 1 )
            !print *, "max pft is ",maxpft
            if ( vfc(i,j,maxpft) < .0001 ) cycle
            
            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,maxpft), laic(i,j,maxpft), laic(i,j,maxpft))
            do m=1,12
               call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE)
     &           ,vfm(m,i,j,maxpft),laim(m,i,j,maxpft),vfc(i,j,maxpft))
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,maxpft),hm(i,j,maxpft),hsd(i,j,maxpft), 
     &           vfc(i,j,maxpft))

          endif
        enddo
      enddo

      ! convert the rest of sparse veg to arid adapted shrub 10
      do j=1,JMn
         do i=1,IMn
            if( vfc(i,j,N_BARE) > .0 .and. laic(i,j,N_BARE) > .0 ) then

               call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &              vfc(i,j,10), laic(i,j,10), .0 )
               do m=1,12
                  call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &                 vfm(m,i,j,10),laim(m,i,j,10),vfc(i,j,10))
                  if ((vfm(m,i,j,10).lt.0.).or.
     &                 (vfm(m,i,j,N_BARE).lt.0.)) then
                     write(*,*) 'vfm<0: ',i,j,m,vfm(m,i,j,10)
     &                    ,vfm(m,i,j,N_BARE),vfc(i,j,10),vfc(i,j,N_BARE)
                     STOP
                  endif
               enddo

               if (vfc(i,j,10) > 0.) then
                  hm(i,j,10) = 2.0 !Check simard.f Set_shrub_height for value!
                  hsd(i,j,10) = undef
               endif
               call convert_vfh(
     &              vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &              vfh(i,j,10),hm(i,j,10),hsd(i,j,10), vfc(i,j,10))

            endif
         enddo
      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif
      titlefoo = "after split"
      write(91) titlefoo, bs_brightratio
      write(91) titlec(N_BARE-1), vfc(:,:,N_BARE-1)
      write(91) titlec(N_BARE), vfc(:,:,N_BARE)
      
      ! check titles
      write(*,*) 'titlec: N_BARE=',N_BARE
      do k=1,N_BARE
         write(*,*) trim(titlec(k))
      enddo
      write(*,*) 'titlem:'
      do m=1,12
         write(*,*) MONTH(m)
         do k=1,N_BARE
            write(*,*) trim(titlem(m,k))
         enddo
      enddo


      call write_output(IMn,JMn,titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_pure_"//res_out//".ij"
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     ,"max_pure_"//trim(fversion),"   ",res_out)

      write(*,*) "pure"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)
      do m=1,12
         call write_output(IM,JM,titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE,
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &        ,"pure_"//trim(fversion),MONTH(m), res_out)
      enddo


      call write_output_h(IMn, JMn, N_BARE, titleh, hm, hsd,
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout//
     &     "_height_pure_"//trim(fversion),"   ",res_out)
      foo(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
           do k=1,KM
              foo(i,j) = foo(i,j) + vfc(i,j,k)
           enddo
        enddo
      enddo
      titlefoo = "LC pure checksum"
      write(100) titlefoo, foo

      call sum_lai_single(100,IMn,JMn,KM,N_BARE,"MODIS MAX", vfc,laic)
      call sum_lai(100,IMn,JMn,KM,N_BARE,"MODIS PURE", vfm,laim)

      ! convert arid adapted shrub with lai < .15 to bare soil
      ! and restrict lai >= .15 
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,10) > .0 .and. laic(i,j,10) < .15 ) then

            if ((IMn.eq.144).and.(JMn.eq.90)) then !## ONLY FOR 2.5X2 !!!
             if ((i.eq.70).and.(j.eq.11)) then
               print *,'vfc10 before',vfc(i,j,10),vfc(i,j,N_BARE),
     &              laic(i,j,10),laic(i,j,N_BARE),i,j
!               STOP
             endif
            endif
            if (vfc(i,j,10).le.0.0) then
               print *,'a vfc10=',vfc(i,j,10),vfc(i,j,N_BARE),
     &              laic(i,j,10),laic(i,j,N_BARE),i,j
               STOP
            endif
            call convert_vf(vfc(i,j,N_BARE), laic(i,j,N_BARE),
     &           vfc(i,j,10), laic(i,j,10), .15 )
                                ! lai >= .15
            if (vfc(i,j,10).le.0.0) then !ERROR CHECK
               print *,'b vfc10=',i,j,vfc(i,j,10),vfc(i,j,N_BARE),
     &              laic(i,j,10),laic(i,j,N_BARE),i,j
               STOP
            endif
            do m=1,12
               if (vfc(i,j,10).gt.0.0) then  
                  call convert_vfm(vfm(m,i,j,N_BARE),laim(m,i,j,N_BARE),
     &                 vfm(m,i,j,10),laim(m,i,j,10), vfc(i,j,10))
               else
                  write(*,*) 'vfc(i,j,10) le 0.:', m,i,j,N_BARE
                  write(*,*) vfc(i,j,:)
                  write(*,*) vfm(m,i,j,:)
                  write(*,*) laic(i,j,:)
                  write(*,*) laim(m,i,j,:)
               endif
               if ((vfm(m,i,j,10).lt.0.).or.(vfm(m,i,j,N_BARE).lt.0.))
     &              then !CHECK ERROR
                  write(*,*) 'vfm<0:',i,j,m,N_BARE
                  write(*,*) vfc(i,j,:)
                  write(*,*) vfm(m,i,j,:)
                  write(*,*) laic(i,j,:)
                  write(*,*) laim(m,i,j,:)
                  STOP
               endif
            enddo

            call convert_vfh(
     &           vfh(i,j,N_BARE),hm(i,j,N_BARE),hsd(i,j,N_BARE),
     &           vfh(i,j,10),hm(i,j,10),hsd(i,j,10), vfc(i,j,10))

          endif
          s = sum(vfc(i,j,1:N_BARE))
          if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
             write(*,*) 'ERROR trim:  max and monthly lc different'
     &            ,i,j, s,sum(vfm(1,i,j,1:N_BARE))
             write(*,*) vfc(i,j,1:N_BARE)
             write(*,*) vfm(1,i,j,1:N_BARE)
          endif
               
        enddo
      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif

      call write_output(IM,JM,titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_"//res_out//".ij"
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     ,"max_trimmed_"//trim(fversion),"   ",res_out)
      write(*,*) "trimmed"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)

      do m=1,12
         call write_output(IM,JM,titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE,
     &        "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &        ,"trimmed_"//trim(fversion),trim(MONTH(m)), res_out)
      enddo

      call write_output_h(IMn,JMN, N_BARE, titleh, hm, hsd, 
     &     "lc_lai_ent16/V"//res_out_int//
     &     "_"//filepreout//"_height_trimmed_"//trim(fversion)
     &     ,"   ",res_out)
      
      foo(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
           do k=1,KM
              foo(i,j) = foo(i,j) + vfc(i,j,k)
           enddo
        enddo
      enddo
      titlefoo = "LC trimmed checksum"
      write(100) titlefoo, foo
      call sum_lai(100,IMn,JMn,KM,N_BARE,"MODIS trimmed", vfm, laim)


      ! rescale fractions so that they sum to 1
      write(*,*) 'Rescaling...'
      do j=1,JMn
         do i=1,IMn
            s = sum( vfc(i,j,1:N_BARE) )
            if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
               write(*,*) 'ERROR scale0:  max and monthly lc different'
     &              ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
               write(*,*) 'vfc',i,j,vfc(i,j,1:N_BARE)
               write(*,*) 'vfm',i,j,vfm(1,i,j,1:N_BARE)
            endif
            if ( s > 1.00001 ) then
!            if ( abs(s-1.0) > 0.00001 ) then
               vfc(i,j,1:N_BARE) = vfc(i,j,1:N_BARE) / s
               vfm(:,i,j,1:N_BARE) = vfm(:,i,j,1:N_BARE) / s
              !vfh(i,j,1:N_BARE) = vfh(i,j,1:N_BARE) / s !Heights are vertical, no rescale
            endif
            s = sum( vfc(i,j,1:N_BARE) ) 
            if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then !#DEBUG
               write(*,*) 'ERROR scale1:  max and monthly lc different'
     &              ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
               write(*,*) 'vfc',i,j,vfc(i,j,1:N_BARE)
               write(*,*) 'vfm',i,j,vfm(1,i,j,1:N_BARE)
            endif
         enddo
      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif


      call write_output(IM,JM,titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_scaled_"
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     ,"max_trimmed_scaled_"//trim(fversion),"   ",res_out)
      write(*,*) "trimmed scaled"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)

      do m=1,12
         call write_output(IM,JM,titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_lai_trimmed_scaled_"
!     &        //MONTH(m)//"_"//res_out//".ij"
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &        ,"trimmed_scaled_"//trim(fversion),MONTH(m), res_out)
      enddo

      call write_output_h(IMn,JMN, N_BARE, titleh, hm, hsd, 
     &     "lc_lai_ent16/V"//res_out_int//
     &     "_"//filepreout//"_height_trimmed_scaled_"//trim(fversion)
     &     ,"   ",res_out)

      foo(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
           do k=1,KM
              foo(i,j) = foo(i,j) + vfc(i,j,k)
           enddo
        enddo
      enddo
      titlefoo = "LC trimmed scaled checksum"
      write(100) titlefoo, foo
      call sum_lai(100,IMn,JMn,KM,N_BARE,"MODIS trimmed scaled"
     &     , vfm, laim)


      !Generate fill-in crop cover from trimmed_scaled before doing nocrops
      !Herb crop only, since right now zero woody crops.
      laicropext(:,:) = 0.d0 !(i,j) only
      laimcropext(:,:,:) = 0.d0 !(m,i,j)
      hmcropext(:,:) = 0.d0 !(i,j)
      call fill_crops(IMn,JMn,vfc(:,:,15),laic(:,:,15)
     i     ,vfm(:,:,:,15),laim(:,:,:,15),hm(:,:,15),hsd(:,:,15)
     o     ,laicropext(:,:),laimcropext(:,:,:)
     o     ,hmcropext,hsdcropext)
      
      titlefoo = "15 - max crops herb ext"
      call write_output_single(IMn,JMn,titlefoo, trim(ent_names16(15)),
     &     vfc(:,:,15), laicropext,
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     ,"max_trimmed_scaled_crops_ext1_"//trim(fversion)
     &     ,"   ",res_out)
      write(*,*) "laimax trimmed scaled crops ext1"

      do m=1,12
         titlefoo = "15 - crops herb ext"
         call write_output_single(IMn,JMn,titlefoo, ent_names16(15),
     &        vfm(m,:,:,15), laimcropext(m,:,:)
     &        ,"lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &        ,"trimmed_scaled_crops_ext1_"//trim(fversion)
     &        ,MONTH(m), res_out)
      enddo
      write(*,*) "lai month trimmed scaled crops ext1"

      titlefoo = "15 - crops herb ext"
      call write_output_h_single(IMn,JMn
     &     , titleh(:,15), 15, hmcropext, hsdcropext
     &     , "lc_lai_ent16/V"//res_out_int//"_"//filepreout//
     &     "_height_trimmed_scaled_crops_ext1_"//trim(fversion)
     &     ,"     ",res_out)
      write(*,*) "height trimmed scaled crops ext1"


      !Generate nocrops
#ifndef COMBINE_CROPS_C3_C4
      write(*,*) 'STOP:  HAVE TO DO THIS CODE'
      STOP
      ! remove combo crops 15 and rescale fractions so that they sum to 1
!      vfc(:,:,15) = 0.
      laic(:,:,15) = 0.
!      vfm(:,:,:,15) = 0.
      laim(:,:,:,15) = 0.
!!      if ( N_VEG == 17 ) then   ! remove separate C4 crops
!        vfc(:,:,16) = 0.
        laic(:,:,16) = 0.
!        vfm(:,:,:,16) = 0.
        laim(:,:,:,16) = 0.
!!      endif
      flag = 0

#else
      !THIS ONLY WORKS FOR N_VEG=16 !!!
      ! remove combo crops 15 and rescale fractions so that they sum to 1

      !* Set crops LAI to zero, keep cover fraction to add to natural veg type
      !* crops C3+C4, excl woody crops
!      vfc(:,:,15) = 0.
      laic(:,:,15) = 0.  
!      vfm(:,:,:,15) = 0.
      laim(:,:,:,15) = 0.
      !* crops woody 
      laic(:,:,16) = 0.        
      !vfm(:,:,:,16) = 0.
        laim(:,:,:,16) = 0.        
      do k=1,18
         titlefoo = 'vfc before '//titlec(k)
         write(997) titlefoo, vfc(:,:,k)
         titlefoo = 'laic before '//titlec(k)
         write(997) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm repcrops'//titlem(m,4)
         write(997) titlefoo, vfm(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'laim repcrops'//titlem(m,4)
         write(997) titlefoo, laim(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'vfm before'//titlem(m,15)
         write(997) titlefoo, vfm(m,:,:,15) !Just check crops
      enddo
      
      flag = 0
      write(*,*) "N_BARE", N_BARE
      do while (flag.le.1)
         nonaturalcount=0
         LAYER(:,:) = 0.d0 !Map the cells with no natural veg
      do j=1,JMn
        do i=1,IMn
           s = sum( vfc(i,j,1:18) ) 
           if (sum(vfc(i,j,15:16))>0.d0) then !Cell has crops
              !write(*,*) 'Cell has crops',i,j,s
              !Replace with closest non-zero natural cover
              call replace_crops(IMn,JMn,KM,i,j !,s
     &             ,bs_brightratio
     &             ,vfc,vfm,laic,laim,naturalfound)
!??     &             ,vfc,vfm,vfc,vfm,naturalfound)
              if (.not.naturalfound) then
                 nonaturalcount = nonaturalcount+1
                 LAYER(i,j) = 1.
              endif
           endif !else do nothing
        enddo
      enddo
      write(*,*) 'Replaced crops'

      !* DEBUG *
      write(*,*) 'Grid cells with no natural veg: ',nonaturalcount
      titlefoo = 'Grid cells with no natural veg: '
      write(996) titlefoo, LAYER
      do k=1,18
         write(*,*) k, titlec(k)
         titlefoo = 'vfc repcrops '//titlec(k)
         write(998) titlefoo, vfc(:,:,k)
         titlefoo = 'laic repcrops '//titlec(k)
         write(998) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm repcrops'//titlem(m,4)
         write(998) titlefoo, vfm(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'laim repcrops'//titlem(m,4)
         write(998) titlefoo, laim(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'vfm repcrops'//titlem(m,15)
         write(998) titlefoo, vfm(m,:,:,15) !Just check crops
      enddo
#endif
      !Rescale after removing crops
      do i=1,IMn
         do j=1,JMn
           ! Rescale always to account for removal of crops 
           ! and to make coasts sum to 1.
           s = sum( vfc(i,j,1:N_BARE) ) 
           if ( (s.gt.0.).and.(s.ne.1.) ) then
              vfc(i,j,1:N_BARE) = vfc(i,j,1:N_BARE) / s
           endif

           do m=1,12
              s = sum( vfm(m,i,j,1:N_BARE) ) 
              if ( (s.gt.0.).and.(s.ne.1.) ) then
                 vfm(m,i,j,1:N_BARE) = vfm(m,i,j,1:N_BARE) / s
              endif
           enddo
           
           s = sum( vfc(i,j,1:N_BARE) ) !#DEBUG
           if (s.ne.sum(vfm(1,i,j,1:N_BARE))) then 
              write(*,*) 'ERROR nocrops:  max and monthly lc differ'
     &             ,i,j, s,sum( vfm(1,i,j,1:N_BARE) ) 
              !write(*,*) vfc(i,j,1:N_BARE)
              !write(*,*) vfm(1,i,j,1:N_BARE)
              do m=1,12
                 vfm(m,i,j,1:N_BARE) = vfc(i,j,N_BARE)
              enddo
           endif
        enddo
      enddo

      !** DEBUG after rescaling crops
      do k=1,18
         write(*,*) k, titlec(k)
         titlefoo = 'vfc rescalecrops '//titlec(k)
         write(999) titlefoo, vfc(:,:,k)
         titlefoo = 'laic rescalecrops '//titlec(k)
         write(999) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm rescalecrops'//titlem(m,4)
         write(999) titlefoo, vfm(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'laim repcrops'//titlem(m,4)
         write(999) titlefoo, laim(m,:,:,4) !Just check evergrneedle
      enddo
      do m=1,12
         titlefoo = 'vfm rescalecrops'//titlem(m,15)
         write(999) titlefoo, vfm(m,:,:,15) !Just check crops
      enddo

      if (.FALSE.) then !DON'T NEED TO DO THIS ANY MORE, DONE IN replace_crops
         !CAUSES SOME WRONG RESCALING
         if (nonaturalcount.eq.0) then
            write(*,*) 'All crops cells successfully fixed.'
            flag = 2
         else
            flag = flag + 1
            if (flag.eq.1) then
               write(*,*) 'Some all-crop cells. Iterating once...'
     &              ,nonaturalcount
            elseif (flag.gt.1) then
               write(*,*) 'Remaining all-crop cells,',nonaturalcount
            endif
         endif
      else
         flag = 2 !Finish the while loop
      endif

      titlefoo = 'All-crop cells with no near natural veg.'
      write(93) titlefoo, LAYER

      enddo                     !do while flag ------

      !## HACK -nk
!      do m=1,12
!         vfm(m,:,:,:) = vfc(:,:,:)
!      enddo

#ifdef SPLIT_BARE_SOIL
      call split_bare_soil(N_VEG, IMn,JMn,KM,N_BARE
     &     ,bs_brightratio,vfc,laic,vfm,laim,vfh,hm,hsd
     &     ,titlec, titlem, titleh,res_out)
#endif

      !** DEBUG after rescaling crops and split_bare
      do k=1,18
         write(*,*) k, titlec(k)
         titlefoo = 'vfc rescalecrops '//titlec(k)
         write(1000) titlefoo, vfc(:,:,k)
         titlefoo = 'laic rescalecrops '//titlec(k)
         write(1000) titlefoo, laic(:,:,k)
      enddo
      do m=1,12
         titlefoo = 'vfm rescalecrops'//titlem(m,15)
         write(1000) titlefoo, vfm(m,:,:,15) !Just check crops
      enddo
      
      call write_output(IM,JM,titlec, vfc, laic, N_BARE,
!     &     "lc_lai_ent16/EntMM16_lc_laimax_trimmed_scaled_nocrops_"
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     ,"max_trimmed_scaled_nocrops_"//trim(fversion),"   ",res_out)
      write(*,*) "trimmed, scaled, no crops"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)


      do m=1,12
         call write_output(IM,JM,titlem(m,:), vfm(m,:,:,:)
     &        , laim(m,:,:,:), N_BARE, 
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &    ,"trimmed_scaled_nocrops_"//trim(fversion),MONTH(m),res_out)
      enddo

      call write_output_h(IM,JM,N_BARE
     &     ,titleh, hm, hsd
     &     ,"lc_lai_ent16/V"//res_out_int//"_"//filepreout//
     &     "_height_trimmed_scaled_nocrops_"//trim(fversion)
     &     ,"   ",res_out)

      foo(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
           do k=1,KM
              foo(i,j) = foo(i,j) + vfc(i,j,k)
           enddo
        enddo
      enddo
      titlefoo = "LC trimmed scaled nocrops checksum"
      write(100) titlefoo, foo

      ! write the final output
      do k=1,N_BARE
        write(7) titlec(k), vfc(:,:,k)
      enddo
      do k=1,N_BARE
        write(7) titlec(k), laic(:,:,k)
      enddo




      lai_yy(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          lai_yy(i,j) = laic(i,j,N_BARE)*vfc(i,j,N_BARE)
        enddo
      enddo

      write(8) title_yy, lai_yy(:,:)

      lai_yy(:,:) = 0.
      do k=1,KM
        if (k==1 .or. k==19) cycle
         lai_yy(:,:) = lai_yy(:,:) + vfn(:,:,k)*lain(:,:,k)
      enddo

      write(9) title_yy, lai_yy(:,:)

      lai_yy(:,:) = 0.
      do k=1,N_BARE
         lai_yy(:,:) = lai_yy(:,:) + vfc(:,:,k)*laic(:,:,k)
      enddo

      write(9) title_yy, lai_yy(:,:)



      ! find the biggest vf for cells which still have sparse veg
      vf_yy(:,:) = -1e30
      lai_yy(:,:) = 0.
      do j=1,JMn
        do i=1,IMn
          if( vfc(i,j,N_BARE) * laic(i,j,N_BARE) > .0 ) then

            maxpft = maxloc( vfc(i,j,1:16), 1 )

            vf_yy(i,j) = vfc(i,j,maxpft)
            lai_yy(i,j) = maxpft

          endif
        enddo
      enddo
      title_yy = "max fraction"
      write(10) title_yy, vf_yy(:,:)
      title_yy = "max pft"
      write(10) title_yy, lai_yy(:,:)

      !------------------------------------------------------------
      ! convert monthly LAI to nocrops vfc.

      ! Merge crops.ext1 laimax, laimonthly, and height into nocrops
      ! to create laimax_ext1, laimonthly_ext1, and height_ext1 files.

      laic(:,:,15) = laicropext(:,:)
      laim(:,:,:,15) = laimcropext(:,:,:)
      hm(:,:,15) =  hmcropext(:,:)
      hsd(:,:,15) = hsdcropext(:,:)

      call write_output_lai(IMn,JMn,titlec, laic, N_BARE,
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     ,"max_trimmed_scaled_ext1_"//trim(fversion),"   ",res_out)
      write(*,*) "trimmed, scaled, ext1"
      !call check_lc_lai_mismatch(KM,IMn,JMn,vfc,laic,'vfc',titlec)

      do m=1,12
         call write_output_lai(IMn,JMn,titlem(m,:)
     &        , laim(m,:,:,:), N_BARE, 
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     , "trimmed_scaled_ext1_"//trim(fversion),MONTH(m),res_out)
      enddo

      call write_output_lai_monthly_nc(IMn,JMn,N_BARE
     &        , laim(:,:,:,:),
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout
     &     , "trimmed_scaled_ext1_"//trim(fversion)
     &     ,"m^2 leaf/m^2 ground",res_out)
      
      call write_output_h(IMn,JMn,N_BARE,titleh, hm, hsd,
     &     "lc_lai_ent16/V"//res_out_int//"_"//filepreout//
     &     "_height_trimmed_scaled_ext1_"//trim(fversion)
     &     ,"   ",res_out)

      write(*,*) 'Done.'

      close(100) !Checksum foo

      end program convert

