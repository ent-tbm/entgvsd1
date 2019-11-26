      module convertnc
      use netcdf
      implicit none
      include 'netcdf.inc'
      save

      public handle_nf_error
      public my_nf_open, my_nf_close
      public my_nf_inq_varid
      public my_nf_get_var_real32
      public my_nf_get_var_real32_2
      public my_nf_get_var_real32_3
      public my_nf_get_var_real32_4
      public my_nf_inq_get_var_real32
      public my_nf_inq_get_var_real32_2
      public my_nf_inq_get_var_real32_3
      public my_nf_inq_get_var_real32_4
      public my_nf_inq_put_var_int
      public my_nf_inq_put_var_char
      public my_nf_inq_put_var_real32
      public my_nf_inq_put_var_real32_2
      public my_nf_inq_put_var_real32_3
      public my_nf_inq_put_var_real32_4
      public my_nf_inq_def_put_var_real32_2
      public my_nf_inq_def_put_var_real32_3
      public my_int2char
      public calc_lon_lat_2HX2, calc_lon_lat
      public my_nf_create_ij
      public my_nf_inq_put_att_any
!      public GISSbands

!  In EntGVSD_util.f      
!      character*3, parameter :: cmon(12) =
!     &     (/ 'JAN','FEB','MAR','APR','MAY','JUN'
!     &       ,'JUL','AUG','SEP','OCT','NOV','DEC' /)

      real*4, dimension(72), parameter :: lon_4x5 =  (/
     &     -177.5, -172.5, -167.5, -162.5, -157.5, -152.5, -147.5,
     &     -142.5, -137.5, -132.5, -127.5, -122.5, -117.5, -112.5,
     &     -107.5, -102.5, -97.5,  -92.5,   -87.5, -82.5, -77.5,
     &     -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5,
     &     -32.5, -27.5, -22.5, -17.5, -12.5, -7.5, -2.5, 2.5, 7.5,
     &     12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5,
     &     57.5, 62.5, 67.5,  72.5, 77.5, 82.5, 87.5, 92.5, 97.5,
     &     102.5, 107.5, 112.5, 117.5, 122.5,  127.5, 132.5, 137.5,
     &     142.5, 147.5, 152.5, 157.5, 162.5, 167.5, 172.5,  177.5 /)

      real*4, dimension(46), parameter :: lat_4x5 = (/
     &     -89, -86, -82, -78, -74, -70, -66, -62, -58, -54, -50,
     &     -46, -42, -38, -34, -30, -26, -22, -18, -14, -10, -6,
     &     -2, 2, 6, 10, 14, 18, 22, 26, 30,  34, 38, 42, 46, 50,
     &     54, 58, 62, 66, 70, 74, 78, 82, 86, 89 /)

!      character(len=19), dimension(6), parameter :: GISSbands = (/
!     &     'VIS (300-770 nm)   ',
!     &     'NIR1 (770-860 nm)  ',
!     &     'NIR2 (860-1250 nm) ',
!     &     'NIR3 (1250-1500 nm)',
!     &     'NIR4 (1500-2200 nm)',
!     &     'NIR5 (2200-4000 nm)' /)
      
      contains

      subroutine calc_lon_lat_4X5(IM,JM,lon,lat)
      implicit none
      integer,intent(in) :: IM,JM
      real*4,intent(inout) :: lon(IM),lat(JM)
      !---Local---
      real*4 :: XX
      integer :: i,j
      
!* Get and write lat, lon.
C     Calc lats & Longs
C
      XX = -90.000000 + 1.00000 !  (4. * .5 = 1.0)
      lat(1) = XX
      XX = XX + 3.00000
      lat(2) = XX
      DO J=3,(JM-1)
         XX = XX + 4.00000
         lat(J) = XX
c        write(0,*) 'Lat ',J,XX
      END DO
      lat(JM) = 90.000000 - 1.00000
C
      XX = -180.00000 + 5*.5 
      DO I=1,IM
         lon(I) = XX
c        write(0,*) 'Lon ',I,XX
         XX = XX + 5.000000
      END DO
C

      end subroutine calc_lon_lat_4X5

      
      subroutine calc_lon_lat_2HX2(IM,JM,lon,lat)
      implicit none
      integer,intent(in) :: IM,JM
      real*4,intent(inout) :: lon(IM),lat(JM)
      !---Local---
      real*4 :: XX
      integer :: i,j,m

!* Get and write lat, lon.
C     Calc lats & Longs
C
      XX = -90.000000 + 1.00000   !  (2. * .5 = 1.0)
      DO J=1,JM
         lat(J) = XX
c        write(0,*) 'Lat ',J,XX
         XX = XX + 2.00000
      END DO
      lat(JM) = 90.000000 - 1.00000
C
      XX = -178.75000
      DO I=1,IM
         lon(I) = XX
c        write(0,*) 'Lon ',I,XX
         XX = XX + 2.500000
      END DO
C
      end subroutine calc_lon_lat_2HX2

!************************************************************************
cc     subroutine calc_lon_lat(IM,JM,lon,lat)
cc      implicit none
cc      integer,intent(in) :: IM,JM
cc      real*4,intent(inout) :: lon(IM),lat(JM)
cc      integer :: i,j
cc
cc      if (IM.eq.144) then
cc         call calc_lon_lat_2HX2(IM,JM,lon,lat)
cc      else
cc
cc         do i=1,IM
cc            lon(i) = -180. + (360./IM)*i - (360./IM)/2.
cc         end do
cc         do j=1,JM
cc            lat(j) = -90. + (180./JM)*j - (180./JM)/2.
cc         end do
cc         endif
cc      end subroutine calc_lon_lat

       subroutine calc_lon_lat(IM,JM,lon,lat)
      implicit none
      integer,intent(in) :: IM,JM 
      !character(len=*) :: res
      real*4,intent(inout) :: lon(IM),lat(JM)
      !---Local---
      real*4 :: XX
      integer :: i,j
      !real*4 :: di, dj

      !*Assumes lon (-180 to 180) and lat (-90 to 90)

      if (IM.eq.72) then   !if (res.eq.'4X5') then
         call calc_lon_lat_4X5(IM,JM,lon,lat)
         !di = 5.0
         !dj = 4.0
         return
      else if (IM.eq.144) then !if (res.eq.'2HX2') then
         call calc_lon_lat_2HX2(IM,JM,lon,lat)
         !di = 2.5
         !dj = 2.0
         return
      else if (IM.eq.360) then  !if (res.eq.'1x1') then
         !call calc_lon_lat_1x1(longout,latout,lon,lat)
         !di = 1.0
         !dj = 1.0
      else if (IM.eq.720) then !if (res.eq.'HXH') then
         !di = 0.5
         !dj = 0.5
      elseif (IM.eq.1440) then !0.25 degree
         !res = 'QXQ'
      elseif (IM.eq.7200) then
        !res = '6km'
      elseif (IM.eq.43200) then
        !res = '1km'
      elseif (IM.eq.86400) then
        !res = '500m'
      else
         write(*,*) 'Fix resolution specs'
         return
      endif

      do i=1,IM
         !lon(i) = -180.0 + (i - 0.5)*di
         lon(i) = -180. + (360./IM)*i - (360./IM)*0.5
      enddo

      do j=1,JM
         !lat(j) = -90.0 + (j - 0.5)*dj
         lat(j) = -90. + (180./JM)*j - (180./JM)*0.5

      enddo

      end subroutine calc_lon_lat     
!************************************************************************
      subroutine handle_nf_error(status, message)
      use netcdf
      integer, intent(in) :: status
      character*(*) :: message
      
       write(*,*) status,message, ' '
     &     ,trim(nf90_strerror(status))
      if(status /= nf90_NoErr)  then
        write(*,*) 'DONT FORGET TO DELETE THE OLD NETCDF FILES FIRST'
	STOP
      endif
      end subroutine handle_nf_error
      
      
      integer function my_nf_open(filein,NFRW,ncid)
      use netcdf
      implicit none
      include 'netcdf.inc'
      character(len=*),intent(in) :: filein
      integer,intent(in) ::NFRW !0=read, NF_WRITE=write
      integer,intent(inout) :: ncid
      !--- Local ---
      integer :: status
      
      status = nf_open(trim(filein),NFRW,ncid)
      write(0,*) trim(filein),'  ',trim(nf90_strerror(status))
      if (status.ne.0) STOP
      my_nf_open = status
      end function my_nf_open


      integer function my_nf_close(ncid)
      use netcdf
      implicit none
      include 'netcdf.inc'
      integer,intent(inout) :: ncid
      !--- Local ---
      integer :: status
      
      status = nf_close(ncid)
      write(0,*) trim(nf90_strerror(status))
      my_nf_close = status
      end function my_nf_close


      integer function my_nf_inq_varid(ncidin,varname,varid)
        use netcdf
        implicit none
        include 'netcdf.inc'
        integer,intent(in) :: ncidin
        character(len=*),intent(in) :: varname
        integer,intent(inout) :: varid
        !--- Local ----
        integer :: status
        
        status = nf_inq_varid(ncidin,varname,varid)
        write(0,*) 'File Inquiry status',status, varname
        my_nf_inq_varid = status
      end function my_nf_inq_varid

      integer function my_nf_get_var_real32(ncidin,
     &     varid,varreal32)
        use netcdf
        implicit none
        include 'netcdf.inc'
        integer,intent(in) :: ncidin
        integer,intent(inout) :: varid
        real*4, dimension(:) :: varreal32
        !--- Local ----
        integer :: status
        
        status = nf_get_var_real(ncidin,varid,varreal32)
        write(0,*) 'File get var status',status,varid
        my_nf_get_var_real32 = status
        end function my_nf_get_var_real32
      
      integer function my_nf_get_var_real32_2(ncidin,
     &     varid,varreal32)
        use netcdf
        implicit none
        include 'netcdf.inc'
        integer,intent(in) :: ncidin
        integer,intent(inout) :: varid
        real*4, dimension(:,:) :: varreal32
        !--- Local ----
        integer :: status
        
        status = nf_get_var_real(ncidin,varid,varreal32)
        write(0,*) 'File get var status',status,varid
        my_nf_get_var_real32_2 = status
      end function my_nf_get_var_real32_2
      
      integer function my_nf_get_var_real32_3(ncidin,
     &     varid,varreal32)
        use netcdf
        implicit none
        include 'netcdf.inc'
        integer,intent(in) :: ncidin
        integer,intent(inout) :: varid
        real*4, dimension(:,:,:) :: varreal32
        !--- Local ----
        integer :: status
        
        status = nf_get_var_real(ncidin,varid,varreal32)
        write(0,*) 'File get var status',status,varid
        my_nf_get_var_real32_3 = status
      end function my_nf_get_var_real32_3
      
      integer function my_nf_get_var_real32_4(ncidin,
     &     varid,varreal32)
        use netcdf
        implicit none
        include 'netcdf.inc'
        integer,intent(in) :: ncidin
        integer,intent(inout) :: varid
        real*4, dimension(:,:,:,:) :: varreal32
        !--- Local ----
        integer :: status
        
        status = nf_get_var_real(ncidin,varid,varreal32)
        write(0,*) 'File get var status',status,varid
        my_nf_get_var_real32_4 = status
      end function my_nf_get_var_real32_4

      
      integer function my_nf_inq_get_var_real32(ncidin,
     &     varname, varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidin
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:) :: varreal32
      !--- Local ----
      integer :: status

      status = my_nf_inq_varid(ncidin,varname,varid)
      write(*,*) 'inq_var ', status, varname
      status = my_nf_get_var_real32(ncidin,varid,varreal32)
      my_nf_inq_get_var_real32 = status
      end function my_nf_inq_get_var_real32

      integer function my_nf_inq_get_var_real32_2(ncidin,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidin
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:,:) :: varreal32
      !--- Local ----
      integer :: status

      status = my_nf_inq_varid(ncidin,varname,varid)
      write(0,*) 'File get var real32_2 status',status,varname,varid
      status = nf_get_var_real(ncidin,varid,varreal32)
      write(0,*) 'File get var real32_2 status',status,varname,varid
      my_nf_inq_get_var_real32_2 = status
      end function my_nf_inq_get_var_real32_2

      integer function my_nf_inq_get_var_real32_3(ncidin,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidin
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:,:,:) :: varreal32
      !--- Local ----
      integer :: status

      status = my_nf_inq_varid(ncidin,varname,varid)
      status = nf_get_var_real(ncidin,varid,varreal32)
      write(0,*) 'File get var real32_3 status',status,varname,varid
      my_nf_inq_get_var_real32_3 = status
      end function my_nf_inq_get_var_real32_3
      
      integer function my_nf_inq_get_var_real32_4(ncidin,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidin
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:,:,:,:) :: varreal32
      !--- Local ----
      integer :: status

      status = my_nf_inq_varid(ncidin,varname,varid)
      status = nf_get_var_real(ncidin,varid,varreal32)
      write(0,*) 'File get var real32_4 status',status,varname,varid
      my_nf_inq_get_var_real32_4 = status
      end function my_nf_inq_get_var_real32_4

      
      integer function my_nf_inq_put_var_int(ncidout,
     &     varname,varid,varint)
      use netcdf
      implicit none
      integer,intent(in) :: ncidout
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      integer, dimension(:) :: varint
      !--- Local ---
      integer :: status
   
      status = my_nf_inq_varid(ncidout,trim(varname),varid)
      status = nf_put_var_int(ncidout,varid, varint)
      write(0,*) 'File put var status',status,trim(varname),varid
      my_nf_inq_put_var_int = status
      end function my_nf_inq_put_var_int

      integer function my_nf_inq_put_var_char(ncidout,
     &     varname,varid,varchar)
      use netcdf
      implicit none
      integer,intent(in) :: ncidout
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      character(len=*), dimension(:) :: varchar
      !--- Local ---
      integer :: status
   
      status = my_nf_inq_varid(ncidout,trim(varname),varid)
      status = nf_put_var_text(ncidout,varid, varchar)
      write(0,*) 'File put var status',status,trim(varname),varid
      my_nf_inq_put_var_char = status
      end function my_nf_inq_put_var_char

      
      integer function my_nf_inq_put_var_real32(ncidout,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidout
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:) :: varreal32
      !--- Local ---
      integer :: status
   
      status = my_nf_inq_varid(ncidout,trim(varname),varid)
      status = nf_put_var_real(ncidout,varid, varreal32)
      write(0,*) 'File put var status',status,trim(varname),varid
      my_nf_inq_put_var_real32 = status
      end function my_nf_inq_put_var_real32

      integer function my_nf_inq_put_var_real32_2(ncidout,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidout
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:,:) :: varreal32
      !--- Local ---
      integer :: status
   
      status = my_nf_inq_varid(ncidout,varname,varid)
      status = nf_put_var_real(ncidout,varid, varreal32)
      write(0,*) 'File put var status',status,varname,varid
      my_nf_inq_put_var_real32_2 = status
      end function my_nf_inq_put_var_real32_2

      integer function my_nf_inq_put_var_real32_3(ncidout,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidout
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:,:,:) :: varreal32
      !--- Local ---
      integer :: status
   
      status = my_nf_inq_varid(ncidout,varname,varid)
      status = nf_put_var_real(ncidout,varid, varreal32)
      write(0,*) 'File put var status',status,varname,varid
      my_nf_inq_put_var_real32_3 = status
      end function my_nf_inq_put_var_real32_3

      integer function my_nf_inq_put_var_real32_4(ncidout,
     &     varname,varid,varreal32)
      use netcdf
      implicit none
      integer,intent(in) :: ncidout
      character(len=*),intent(in) :: varname
      integer,intent(inout) :: varid
      real*4, dimension(:,:,:,:) :: varreal32
      !--- Local ---
      integer :: status
      
      status = my_nf_inq_varid(ncidout,varname,varid)
      status = nf_put_var_real(ncidout,varid, varreal32)
      write(0,*) 'File put var status',status,varname,varid
      my_nf_inq_put_var_real32_4 = status
      end function my_nf_inq_put_var_real32_4

      character*5 function my_int2char(numint)
        use netcdf
        implicit none
        integer, intent(in) :: numint
        !--- Local ---
        character(len=8) :: fmt ! format descriptor
        character*5 :: x1
        
        fmt = '(I5)'   ! an integer of width 5 with zeros at the left
        write (x1,fmt) numint   ! convert integer to string using a 'internal file'
        my_int2char = x1
      end function my_int2char

!     !------------------------------------------------------------------
      subroutine my_nf_create_ij(filename,IM,JM,
     &   ncid,dimlon,dimlat)
      implicit none
      include 'netcdf.inc'
      character(len=*), intent(in) :: filename
      integer,intent(in) :: IM,JM
      integer,intent(out) :: ncid
      integer,intent(out) :: dimlon, dimlat
      !--- Local ----
      integer :: status, idlon, idlat, varid
      real*4 :: lon(IM), lat(JM)
      !character*80 :: filenc
      character*10 :: res

      if (IM.eq.72) then
        res = '72X46' !5X4
      elseif (IM.eq.144) then
        res = '2HX2' !'144X90'
      elseif (IM.eq.360) then
        res = '1X1'
      elseif (IM.eq.720) then
        res = 'HXH'
      elseif (IM.eq.1440) then
        res = 'QXQ'
      elseif (IM.eq.7200) then
        res = '6km'
      elseif (IM.eq.43200) then
        res = '1km'
      elseif (IM.eq.86400) then
        res = '500m'
      endif

      !filenc =trim(FILEPREFIX)//trim(filename)//'_'//trim(res)//'.nc'
      !filenc = trim(filenc)
      !filenc = trim(filename)
      !Open file if already exists, otherwise make new file
      !if ( nf_open(filename, NF_WRITE, ncid).ne.nf_noerr)
      status=nf_open(filename, NF_WRITE, ncid) !Get ncid if file exists
      write(0,*) 'nf_open status', status, filename
      if (status.ne.nf_noerr)
     &        then
        ! netcdf output file needs to be created
           write(0,*) 'Creating ',filename
           status=nf_create(filename, NF_NOCLOBBER, ncid)
           write(0,*) 'nf_create status',status,filename,ncid
           status=nf_def_dim(ncid, 'lon', IM, dimlon)
           status=nf_def_dim(ncid, 'lat', JM, dimlat)
           !write(0,*) 'dimlon', status, dimlon
           !write(0,*) 'dimlat', status, dimlat
           status=nf_def_var(ncid, 'lon', NF_FLOAT, 1, dimlon, idlon)
           write(0,*) 'nf_def_var lon',status,ncid,dimlon
           status=nf_def_var(ncid, 'lat', NF_FLOAT, 1, dimlat, idlat)
           write(0,*) 'nf_def_var lat',status,ncid,dimlat
           status=nf_put_att_text(ncid, idlon, 'long_name'
     &          ,len('longitude'), 'longitude')
           status=nf_put_att_text(ncid, idlat, 'long_name'
     &          , len('latitude'), 'latitude')
           status=nf_put_att_text(ncid, idlon, 'units'
     &          , len('degrees east'), 'degrees east')
           status=nf_put_att_text(ncid, idlat, 'units'
     &          , len('degrees north'), 'degrees north')
           status=nf_put_att_real(ncid, idlon, '_FillValue', NF_FLOAT
     &          , 1, -1.e30)
           status=nf_put_att_real(ncid, idlat, '_FillValue', NF_FLOAT
     &          , 1, -1.e30)

           status=nf_enddef(ncid)
           call calc_lon_lat(IM,JM,lon,lat)
           !write(*,*) 'lon', lon
           !write(*,*) 'lat', lat
           status=my_nf_inq_put_var_real32(ncid,'lon'
     &          ,varid,lon)
           !status=nf_put_vara_real(ncid,idlon,1,IM,lon)
           write(0,*) 'put lon', status, varid, dimlon, IM
           status=my_nf_inq_put_var_real32(ncid,'lat'
     &          ,varid,lat)
           !status=nf_put_vara_real(ncid,idlat,1,JM,lat)
           write(0,*) 'put lat', status, varid, dimlat, JM
           !status=nf_enddef(ncid)
           status=nf_close(ncid)
           write(0,*) 'Created netcdf ij file'
      else
           dimlon=1
           dimlat=2
           write(0,*) 'Netcdf file was previously created. ',filename
      endif



      end subroutine my_nf_create_ij

!     !------------------------------------------------------------------
      integer function my_nf_inq_def_put_var_int(ncid,
     &     KM,K0,K1,dimk, varname, long_name, units, varint)

      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: ncid
      !character(len=*) :: file
      integer, intent(in) ::  KM,K0,K1,dimk
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: long_name
      character(len=*),intent(in) :: units
      integer, dimension(K0:K1) :: varint
      !--- Local ---
      integer :: status, varid
      integer :: dimgrid(1), start(1), count(1)

      if ( nf_inq_varid(ncid,varname,varid).ne.NF_NOERR) then
         !Need to def var
         status=nf_redef(ncid)
         write(0,*) 'File redef status ', varname,status, ncid
         dimgrid(1) = dimk
         write(0,*) 'dimgrid', dimgrid
         status=nf_def_var(ncid, varname, NF_INT, 1,
     &        dimgrid, varid)
         write(0,*) 'File def var status',status, varname,varid
         status=nf_put_att_text(ncid, varid, 'long_name',
     &        len(long_name), long_name)
         status=nf_put_att_text(ncid, varid, 'units'
     &        , len(units), units)
         write(0,*) 'Put att status ',status, varname,units,varid
!         status=nf_put_att_real(ncid, varid, '_FillValue', NF_FLOAT
!     &          , 1, -1.e30)
!         write(0,*) 'nf_put_att_real status ', status
         status = nf_enddef(ncid)
         write(0,*) 'nf_enddef status ', status
      endif
       
      if (nf_inq_varid(ncid,varname,varid) == NF_NOERR) then
         start(:) = (/K0/)
         count(:) = (/K1-K0+1/)
         status=nf_put_vara_int(ncid, varid,start,count,varint)
         write(0,*) 'File put vara status',status,varname,varid
      endif
      
      my_nf_inq_def_put_var_int = status !return
      
      end function my_nf_inq_def_put_var_int


!     !------------------------------------------------------------------
      integer function my_nf_inq_def_put_var_real32_2(ncid,
     &     IM,JM,I0,I1,J0,J1,dimlon,dimlat,
     &     varname, long_name, units, varreal32)

      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: ncid
      !character(len=*) :: file
      integer, intent(in) ::  IM,JM,I0,I1,J0,J1,dimlon,dimlat
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: long_name
      character(len=*),intent(in) :: units
      real*4, dimension(I0:I1,J0:J1) :: varreal32
      !--- Local ---
      integer :: status, varid
      integer :: dimgrid(2), start(2), count(2)
      !real*4, dimension(I0:I1,J0:J1) :: varreal32

      !varreal32(:,:) = varreal64(:,:)
      if ( nf_inq_varid(ncid,varname,varid).ne.NF_NOERR) then
         !Need to def var
         status=nf_redef(ncid)
         write(0,*) 'File redef status ', varname,status, ncid
         dimgrid(1) = dimlon
         dimgrid(2) = dimlat
         write(0,*) 'dimgrid', dimgrid
         status=nf_def_var(ncid, varname, NF_FLOAT, 2,
     &        dimgrid, varid)
         write(0,*) 'File def var status',status, varname,varid
         status=nf_put_att_text(ncid, varid, 'long_name',
     &        len(long_name), long_name)
         status=nf_put_att_text(ncid, varid, 'units'
     &        , len(units), units)
         write(0,*) 'Put att status ',status, varname,units,varid
         status=nf_put_att_real(ncid, varid, '_FillValue', NF_FLOAT
     &          , 1, -1.e30)
         write(0,*) 'nf_put_att_real status ', status
         status = nf_enddef(ncid)
         write(0,*) 'nf_enddef status ', status
      endif
       
      if (nf_inq_varid(ncid,varname,varid) == NF_NOERR) then
         start(:) = (/I0, J0/)
         count(:) = (/I1-I0+1, J1-J0+1/)
         status=nf_put_vara_real(ncid, varid,start,count,varreal32)
         write(0,*) 'File put vara status',status,varname,varid
      endif
      
      my_nf_inq_def_put_var_real32_2 = status !return
      
      end function my_nf_inq_def_put_var_real32_2


!     !------------------------------------------------------------------
      integer function my_nf_inq_def_put_var_real32_3(ncid,
     &     IM,JM,KM,I0,I1,J0,J1,dimlon,dimlat,dim3,
     &     varname, long_name, units, varreal32)

      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: ncid
      !character(len=*) :: file
      integer, intent(in) ::  IM,JM,KM,I0,I1,J0,J1,dimlon,dimlat,dim3
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: long_name
      character(len=*),intent(in) :: units
      real*4, dimension(1:KM,I0:I1,J0:J1) :: varreal32
      !--- Local ---
      integer :: status, varid
      integer :: dimvar(3), start(3), count(3)
      !real*4, dimension(I0:I1,J0:J1) :: varreal32

      !varreal32(:,:) = varreal64(:,:)

      if ( nf_inq_varid(ncid,varname,varid).ne.NF_NOERR) then
         !Need to def var
         status=nf_redef(ncid)
         write(0,*) 'File redef status ', varname,status, ncid
!         dimvar(1) = dim3
!         dimvar(2) = dimlon
!         dimvar(3) = dimlat
         dimvar(1) = dimlon
         dimvar(2) = dimlat
         dimvar(3) = dim3
         write(0,*) 'dimvar', dimvar
         status=nf_def_var(ncid, varname, NF_FLOAT, 3,
     &        dimvar, varid)
         write(0,*) 'File def var status',status, varname,varid
         status=nf_put_att_text(ncid, varid, 'long_name',
     &        len(long_name), long_name)
         status=nf_put_att_text(ncid, varid, 'units'
     &        , len(units), units)
         write(0,*) 'Put att status ',status, varname,units,varid
!         status=nf_put_att_text(ncid, varid, 'bands'
!     &        , len(descriptiondim3), descriptiondim3)
!         write(0,*) 'Put att status ',status, varname,units,varid
         status=nf_put_att_real(ncid, varid, '_FillValue', NF_FLOAT
     &          , 1, -1.e30)
         write(0,*) 'nf_put_att_real status ', status
         status = nf_enddef(ncid)
         write(0,*) 'nf_enddef status ', status
      endif
       
      if (nf_inq_varid(ncid,varname,varid) == NF_NOERR) then
         start(1) = 1
         count(1) = KM
         start(2:3) = (/I0, J0/)
         count(2:3) = (/I1-I0+1, J1-J0+1/)
         status=nf_put_vara_real(ncid, varid,start,count,varreal32)
         write(0,*) 'File put vara status',status,varname,varid
      endif
      
      my_nf_inq_def_put_var_real32_3 = status !return
      
      end function my_nf_inq_def_put_var_real32_3

!     !------------------------------------------------------------------
      integer function my_nf_inq_def_put_var_real32_3t(ncid,
     &     IM,JM,KM,I0,I1,J0,J1,dimlon,dimlat,dim3,
     &     varname, long_name, units, varreal32)

      !* Same as my_nf_inq_def_put_var_real32_3 but varreal32 indices have K last.
      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: ncid
      !character(len=*) :: file
      integer, intent(in) ::  IM,JM,KM,I0,I1,J0,J1,dimlon,dimlat,dim3
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: long_name
      character(len=*),intent(in) :: units
      real*4, dimension(I0:I1,J0:J1,1:KM) :: varreal32  
      !--- Local ---
      integer :: status, varid
      integer :: dimvar(3), start(3), count(3)
      !real*4, dimension(I0:I1,J0:J1) :: varreal32

      !varreal32(:,:) = varreal64(:,:)

      if ( nf_inq_varid(ncid,varname,varid).ne.NF_NOERR) then
         !Need to def var
         status=nf_redef(ncid)
         write(0,*) 'File redef status ', varname,status, ncid
!         dimvar(1) = dim3
!         dimvar(2) = dimlon
!         dimvar(3) = dimlat
         dimvar(1) = dimlon
         dimvar(2) = dimlat
         dimvar(3) = dim3
         write(0,*) 'dimvar', dimvar
         status=nf_def_var(ncid, varname, NF_FLOAT, 3,
     &        dimvar, varid)
         write(0,*) 'File def var status',status, varname,varid
         status=nf_put_att_text(ncid, varid, 'long_name',
     &        len(long_name), long_name)
         status=nf_put_att_text(ncid, varid, 'units'
     &        , len(units), units)
         write(0,*) 'Put att status ',status, varname,units,varid
!         status=nf_put_att_text(ncid, varid, 'bands'
!     &        , len(descriptiondim3), descriptiondim3)
!         write(0,*) 'Put att status ',status, varname,units,varid
         status=nf_put_att_real(ncid, varid, '_FillValue', NF_FLOAT
     &          , 1, -1.e30)
         write(0,*) 'nf_put_att_real status ', status
         status = nf_enddef(ncid)
         write(0,*) 'nf_enddef status ', status
      endif
       
      if (nf_inq_varid(ncid,varname,varid) == NF_NOERR) then
         start(1) = 1
         start(2:3) = (/I0, J0/)
         count(1:2) = (/I1-I0+1, J1-J0+1/)
         count(3) = KM
         status=nf_put_vara_real(ncid, varid,start,count,varreal32)
         write(0,*) 'File put vara status',status,varname,varid
      endif
      
      my_nf_inq_def_put_var_real32_3t = status !return
      
      end function my_nf_inq_def_put_var_real32_3t

!     !------------------------------------------------------------------
      integer function my_nf_inq_put_att_varid(ncid,
     &     varname, long_name, units)
!     For already created variable, add attributes.

      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: ncid
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: long_name
      character(len=*),intent(in) :: units
      !--- Local ---
      integer :: status, varid

      if (nf_inq_varid(ncid,varname,varid) == NF_NOERR) then
         status=nf_redef(ncid)
         write(0,*) 'File redef status ', varname,status, ncid
         status=nf_put_att_text(ncid, varid, 'long_name',
     &        len(long_name), long_name)
         write(0,*) 'Put long_name status',status,varid
         status=nf_put_att_text(ncid, varid, 'units'
     &        , len(units), units)
         write(0,*) 'Put att status ',status, varname,units,varid
         status=nf_put_att_real(ncid, varid, '_FillValue', NF_FLOAT
     &        , 1, -1.e30)
         status = nf_enddef(ncid)
      endif
      
      end function my_nf_inq_put_att_varid

!     !------------------------------------------------------------------
      integer function my_nf_inq_put_att_any(ncid,
     &     varname, anyatt, anytext)
!     For already created variable, add attributes.
      
      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: ncid
      character(len=*),intent(in) :: varname
      character(len=*),intent(in) :: anyatt
      character(len=*),intent(in) :: anytext
      !--- Local ---
      integer :: status, varid

      if (nf_inq_varid(ncid,varname,varid) == NF_NOERR) then
         status=nf_redef(ncid)
         write(0,*) 'File redef status ', varname,status, ncid
         status=nf_put_att_text(ncid, varid, anyatt,
     &        len(anytext), anytext)
         write(0,*) 'Put att status ',status,varname,varid
         status = nf_enddef(ncid)
      endif
      my_nf_inq_put_att_any = status

      end function my_nf_inq_put_att_any



!     integer function my_nf_inq_put_global(ncid,field,string)
!     !Add global attribute to nc file
!     implicit none
!     include 'netcdf.inc'
!     integer :: ncid
!     character(len=*) :: field
!     character(len=*) :: string
!     !-- Local ----
!     integer :: status
!     
!     status=nf_redef(ncid)
!     write(0,*) 'File redef status ', status, ncid
!     status=nf_put_global_attr(ncid, field, string)
!     write(0,*) 'Global attribute put ', status, ncid
!     status = nf_enddef(ncid)
!     
!     end function my_nf_inq_put_global
      
!     !------------------------------------------------------------------
      
      end module convertnc
      
