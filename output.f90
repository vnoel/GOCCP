!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!! NETCDF RECORDING SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


!----------------------------------------------------------------------------!
! *** CHECK *** This subroutine check the status of the nf90 command         !
!----------------------------------------------------------------------------!
subroutine check(status)
    use netcdf

    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, nf90_strerror(status)
      stop "Stopped"
    end if

end subroutine check  
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** CREATE_PROFNC *** This routine create a netcdf prof file, and its      !
!                       dimensions                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MAP3D file                              !
! dname      : period of MAP3D file (description of ncdf file)               !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of the altitude boxes                                  !
! vlon       : values of model longitude from 180 to -180                    !
! vlat       : values of model latitude from 90 to -90                       !
! vpres      : values of model altitude from 0 to 19.2km or 40.5km           !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
! dim        : dimension id of the MaP3D variables recorded in the ncdf      !
!              files                                                         !
!----------------------------------------------------------------------------!
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
! ndims      : dimension of dim (=4 : lon,lat,alt,time)                      !
! nc_id      : netcdf file id                                                !
! lon_varid  : variable id of longitude                                      !
! lat_varid  : variable id of latitude                                       !
! pres_varid : variable id of altitude                                       !
! time_varid : variable id of time                                           !
! lon_dimid  : dimension id of longitude                                     !
! lat_dimid  : dimension id of latitude                                      !
! pres_dimid : dimension id of altitude                                      !
! time_dimid : dimension id of time                                          !
!----------------------------------------------------------------------------!
!                                                                            !
! ex  create_profnc(file8,file9,lonmod,latmod,altmod,resd,dimidsp,altmax,    !
!                   lonmax,latmax)                                           !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine create_profnc(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtime,dim,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim2(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_Cloud_Fraction_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))

    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))


    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

 
   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_profnc
!----------------------------------------------------------------------------!


subroutine create_temp3d(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,  &
                         vtime,dim,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim2(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_Cloud_Fraction_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))

    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'temp', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))


    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'temp_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the temperature bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'temp_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the temperature bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

 
   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_temp3d
!----------------------------------------------------------------------------!


subroutine create_depolnc3d(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtime,dim,dim2,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2, ncat1=5,ncat2=25
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim3(nv),dim2(5),dim4(2)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2,cat_varid
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid,cat_dimid1,cat_dimid2
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    character(len=25),dimension(ncat1,ncat2)  ::  vcat
    

vcat(1,:)='UNDEFINED'
vcat(2,:)='FALSE LIQ'
vcat(3,:)='FALSE ICE'
vcat(4,:)='Horizontally Oriented'
vcat(5,:)='Unphysical value (NOISE)'

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_CloudFraction_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))  
 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'cat1', ncat1, cat_dimid1))
    call check(nf90_def_dim(nc_id, 'cat2', ncat2, cat_dimid2))

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, alt_dimid, cat_dimid1, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/cat_dimid2, cat_dimid1/)
    
    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))
 
    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim4, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, cat_varid, vcat))
    call check(nf90_put_var(nc_id, time_varid, vtime))



    call check(nf90_close(nc_id))

end subroutine create_depolnc3d
!----------------------------------------------------------------------------!

subroutine create_ind3d(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtime,dim,dim2,alt,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer, parameter ::  ndims = 4, nv=2, ncat1=5,ncat2=25
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat ,dim(ndims), alt,dim3(nv),dim2(5),dim4(2)
    integer  ::  lon_varid,lat_varid,alt_varid, time_varid, nc_id,alt_varid2,cat_varid
    integer  ::  lon_dimid,lat_dimid,alt_dimid,time_dimid, nv_dimid,cat_dimid1,cat_dimid2
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    character(len=25),dimension(ncat1,ncat2)  ::  vcat
    

vcat(1,:)='UNDEFINED'
vcat(2,:)='FALSE LIQ'
vcat(3,:)='FALSE ICE'
vcat(4,:)='Horizontally Oriented'
vcat(5,:)='Unphysical value (NOISE)'

  
    call date_and_time(date,time,zone,value)
    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_Three-dimensionnal_CloudFraction_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))  
 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'cat1', ncat1, cat_dimid1))
    call check(nf90_def_dim(nc_id, 'cat2', ncat2, cat_dimid2))

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/lon_dimid, lat_dimid, alt_dimid, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, alt_dimid, cat_dimid1, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/cat_dimid2, cat_dimid1/)
    
    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))
 
    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim4, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  

    call check(nf90_enddef(nc_id))
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, cat_varid, vcat))
    call check(nf90_put_var(nc_id, time_varid, vtime))



    call check(nf90_close(nc_id))

end subroutine create_ind3d
!----------------------------------------------------------------------------!




!----------------------------------------------------------------------------!
! *** INSTANTSR2NC *** This routine create & record a netcdf instant_SR file,!
!                      and its dimensions                                    !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf instant_SR file                         !
! daynight   : switch which allow to select day or night mode                !
! mod        : name of the grid selected                                     !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of the altitude boxes                                  !
! longi      : Longitude in degrees, has dimension (nprof)                   !
! lati       : Latitude in degrees, has dimension (nprof)                    !
! vprestop   : values of model altitude from 0 to 19.2km or 40.5km           !
! SEi        : Surface_Elevation in kilometers, has dimension (nprof)        !
! timei      : Time converted in UTC fractionned hour of a day               !
! SR         : scattering ratio calculated, dim=(altmax,nprof)               !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
!----------------------------------------------------------------------------!
! fname2     : period of instant_SR file (description of ncdf file)          !
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
! dim        : dimension id of the SR                                        !
! ndims      : dimension of dim (= 2 : alt,nprof)                            !
! nc_id      : netcdf file id                                                !
! varid3     : variable id of longitude                                      !
! varid2     : variable id of latitude                                       !
! alt_varid  : variable id of altitude                                       !
! varid5     : variable id of time                                           !
! varid6     : variable id of SE                                             !
! varid1     : variable id of instant_SR                                     !
! alt_dimid  : dimension id of altitude                                      !
! it_dimid   : dimension id of nprof                                         !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : instantSR2nc(file4,altmod,resd,altmax,switch,gcm,it,lat,lon,SE,temps2,!
!                   srmoy)                                                   !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine instantSR2nc(fname,vprestop_mid,vprestop_bound,vtime,alt,daynight,mod,nprof,lati,    &
                        longi,SEi,timei,SR)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*14
    character  ::  date*8,time*10,zone*5,daynight*5,mod*8
    integer, parameter ::  ndims = 2, nv=2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof,dim2(ndims)
    integer  ::  alt_varid, varid1, varid2, varid3, varid5, varid6, nc_id,alt_varid2
    integer  ::  alt_dimid, it_dimid, nv_dimid, cat_varid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

    real*4,dimension(alt,nprof)  ::  SR
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.

fname2=fname(12:25)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP/GOCCP/instant_SR/temp/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',fname2))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 


    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))

    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))


    dim = (/alt_dimid, it_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'rejected_value',rej))    
    call check(nf90_put_att(nc_id, varid1, 'Surface_value',se))

    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_close(nc_id))

end subroutine instantSR2nc
!----------------------------------------------------------------------------!

subroutine SR_CR_DR_2nc(fname,vprestop_mid,vprestop_bound,vtime,alt,daynight,mod,nprof,lati,    &
                        longi,SEi,timei,SR,CR,DEPOL)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*14
    character  ::  date*8,time*10,zone*5,daynight*5,mod*4
    integer, parameter ::  ndims = 2, nv = 2 
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof, dim2(ndims)
    integer  ::  alt_varid, alt_varid2, varid1, varid2, varid3, varid5, varid6,varid7,varid8, nc_id
    integer  ::  alt_dimid, it_dimid, nv_dimid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real*4,dimension(alt,nprof)  ::  SR,CR,DEPOL
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

fname2=fname(18:39)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP/GOCCP/instant_SR_CR_DR/temp/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_DR_CR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',trim(fname2)))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 

    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/alt_dimid, it_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer')) 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid1, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_CR', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(nc_id, varid7, 'lon_name',                       &
               'instantaneous Color Ratio'))
    call check(nf90_put_att(nc_id, varid7, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid7, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid7, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid7, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_DR', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(nc_id, varid8, 'lon_name',                       &
               'instantaneous Depolarization Ratio'))
    call check(nf90_put_att(nc_id, varid8, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid8, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid8, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid8, 'Surface_value',se))


    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_put_var(nc_id, varid7, CR))
    call check(nf90_put_var(nc_id, varid8, DEPOL))



    call check(nf90_close(nc_id))

end subroutine SR_CR_DR_2nc
!----------------------------------------------------------------------------!

! Same routine as SR_CR_DR_2nc including the record of ATB & ATBmol variables
subroutine SR_CR_DR_ATB_nc(fname,vprestop_mid,vprestop_bound,vtime,alt,      &
                           daynight,mod,nprof,lati,longi,SEi,timei,SR,CR,    &
                           DEPOL,atb,atbm,atbr,atbl,temp,cloudfrac)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*21
    character  ::  date*8,time*10,zone*5,daynight*5,mod*4
    integer, parameter ::  ndims = 2, nv = 2 
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof, dim2(ndims)
    integer  ::  alt_varid, alt_varid2, varid1, varid2, varid3, varid5, varid6,varid7,varid8,varid9,varid10, nc_id,varid11,varid12,varid13, varid14
    integer  ::  alt_dimid, it_dimid, nv_dimid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real,dimension(alt,nprof)  ::  SR,CR,DEPOL,atb,atbm,atbl,atbr,temp,cloudfrac
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound

    fname2=fname(18:39)

    call date_and_time(date,time,zone,value)
    call check(nf90_create(fname, NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_DR_CR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',trim(fname2)))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 

    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    dim = (/alt_dimid, it_dimid/)
    dim2 = (/alt_dimid, nv_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim2, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer')) 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid1, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_Cloud', NF90_FLOAT, dim, varid14))
    call check(nf90_put_att(nc_id, varid14, 'lon_name',                       &
               'instantaneous Cloud Mask'))
    call check(nf90_put_att(nc_id, varid14, 'Legend','1=NaN, 2=Clear, 3=Cloud, 4=Uncertain, 6=Surface, 7=Rejected, 8=Fully Attenuated' ))
    call check(nf90_put_att(nc_id, varid14, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'instant_CR', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(nc_id, varid7, 'lon_name',                       &
               'instantaneous Color Ratio'))
    call check(nf90_put_att(nc_id, varid7, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid7, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid7, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid7, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'instant_DR', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(nc_id, varid8, 'lon_name',                       &
               'instantaneous Depolarization Ratio'))
    call check(nf90_put_att(nc_id, varid8, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid8, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid8, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid8, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB', NF90_FLOAT, dim, varid9))
    call check(nf90_put_att(nc_id, varid9, 'lon_name',                       &
               'Attenuated Total Backscatter'))
    call check(nf90_put_att(nc_id, varid9, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid9, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid9, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid9, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB_mol', NF90_FLOAT, dim, varid10))
    call check(nf90_put_att(nc_id, varid10, 'lon_name',                       &
               'ATB of Molecular'))
    call check(nf90_put_att(nc_id, varid10, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid10, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid10, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid10, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB_per', NF90_FLOAT, dim, varid11))
    call check(nf90_put_att(nc_id, varid11, 'lon_name',                       &
               'Attenuated Perpendicular Backscatter'))
    call check(nf90_put_att(nc_id, varid11, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid11, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid11, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid11, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'ATB_par', NF90_FLOAT, dim, varid12))
    call check(nf90_put_att(nc_id, varid12, 'lon_name',                       &
               'Attenuated Parallel Backscatter'))
    call check(nf90_put_att(nc_id, varid12, 'units','km-1 sr-1'))
    call check(nf90_put_att(nc_id, varid12, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid12, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid12, 'Surface_value',se))

    call check(nf90_def_var(nc_id, 'TEMP', NF90_FLOAT, dim, varid13))
    call check(nf90_put_att(nc_id, varid13, 'lon_name',                       &
               'Temperature'))
    call check(nf90_put_att(nc_id, varid13, 'units','Celcius degree'))
    call check(nf90_put_att(nc_id, varid13, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid13, 'Rejected_value',rej))
    call check(nf90_put_att(nc_id, varid13, 'Surface_value',se))


    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_put_var(nc_id, varid14, cloudfrac))
    call check(nf90_put_var(nc_id, varid7, CR))
    call check(nf90_put_var(nc_id, varid8, DEPOL))
    call check(nf90_put_var(nc_id, varid9, atb))
    call check(nf90_put_var(nc_id, varid10, atbm))
    call check(nf90_put_var(nc_id, varid11, atbr))
    call check(nf90_put_var(nc_id, varid12, atbl))
    call check(nf90_put_var(nc_id, varid13, temp))




    call check(nf90_close(nc_id))

end subroutine SR_CR_DR_ATB_nc
!----------------------------------------------------------------------------!


subroutine instant_phase(fname,nalt,nprof,phase)
use netcdf
implicit none
character(len=*)  ::  fname
integer,dimension(2)  ::  dimm
integer  ::  ncid,varid, alt_dimid,it_dimid,ndims
integer  ::  nalt,nprof
real*4  ::  phase(nalt,nprof)
real,parameter  ::   nan=-9999.


  call check(NF90_OPEN('/bdd/CFMIP/GOCCP/instant_SR_CR_DR/temp/'//trim(fname),NF90_WRITE,ncid))

  call check(NF90_REDEF(ncid))

    call check(nf90_inq_dimid(ncid, 'altitude',  alt_dimid))
    call check(nf90_inq_dimid(ncid, 'it', it_dimid))

    dimm = (/alt_dimid, it_dimid/)

    call check(nf90_def_var(ncid, 'instant_Phase', NF90_FLOAT, dimm, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                       &
               'instantaneous Cloud Phase Mask'))
    call check(nf90_put_att(ncid, varid, 'Legend','1=LIQ, 2=ICE, 3=UNDEFINED, 4=FALSE LIQ, 5=FALSE ICE, 6=Horizontally Oriented, 7=Unphysical value (NOISE)' ))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_put_att(ncid, varid, 'Surface_value','-888'))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid, phase))

    call check(nf90_close(ncid))

end subroutine instant_phase


subroutine SR_DEPOL_2nc(fname,vprestop,vtime,alt,daynight,mod,nprof,lati,    &
                        longi,SEi,timei,SR,DEPOL)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname
    character  ::  fname2*13
    character  ::  date*8,time*10,zone*5,daynight*5,mod*8
    integer, parameter ::  ndims = 2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims), alt, nprof
    integer  ::  alt_varid, varid1, varid2, varid3, varid5, varid6,varid7,varid8, nc_id
    integer  ::  alt_dimid, it_dimid
    real  ::  vtime
    real*8, dimension(nprof) :: timei  
    real, dimension(nprof) :: SEi    
    real, dimension(nprof) :: longi
    real, dimension(nprof) :: lati
    real,dimension(alt)  ::  vprestop
    real*4,dimension(alt,nprof)  ::  SR,DEPOL
    real,parameter  ::   nan=-9999.

fname2=fname(12:26)

    call date_and_time(date,time,zone,value)
    call check(nf90_create('/bdd/CFMIP_TEMP/DEPOL/'//fname,     &
               NF90_CLOBBER, nc_id))

    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description',               &
               'GOCCP_instant_SR_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',fname2))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))

    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'it', nprof, it_dimid))

    dim = (/alt_dimid, it_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, it_dimid, varid3))
    call check(nf90_put_att(nc_id, varid3, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, varid3, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, varid3, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, it_dimid, varid2))
    call check(nf90_put_att(nc_id, varid2, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, varid2, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, varid2, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, it_dimid, varid5))
    call check(nf90_put_att(nc_id, varid5, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, varid5, 'units','Fractionned UTC hour'))
    call check(nf90_put_att(nc_id, varid5, 'axis','T')) 

    call check(nf90_def_var(nc_id, 'SE', NF90_FLOAT, it_dimid, varid6))
    call check(nf90_put_att(nc_id, varid6, 'lon_name','Surface_Elevation'))
    call check(nf90_put_att(nc_id, varid6, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'instant_SR', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(nc_id, varid1, 'lon_name',                       &
               'instantaneous Scattering Ratio'))
    call check(nf90_put_att(nc_id, varid1, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid1, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid1, 'Surface_value','-888'))

    call check(nf90_def_var(nc_id, 'instant_DEPOL', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(nc_id, varid8, 'lon_name',                       &
               'instantaneous Depolarization'))
    call check(nf90_put_att(nc_id, varid8, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, varid8, '_FillValue',nan))
    call check(nf90_put_att(nc_id, varid8, 'Surface_value','-888'))

    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, alt_varid, vprestop))
    call check(nf90_put_var(nc_id, varid2, lati))
    call check(nf90_put_var(nc_id, varid3, longi))
    call check(nf90_put_var(nc_id, varid5, timei))
    call check(nf90_put_var(nc_id, varid6, SEi))

    call check(nf90_put_var(nc_id, varid1, SR))
    call check(nf90_put_var(nc_id, varid8, DEPOL))



    call check(nf90_close(nc_id))

end subroutine SR_DEPOL_2nc
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** PROF_RECVAR2NC *** This routine record the prof variables in the netcdf!
!                        prof file                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MAP3D file                              !
! dim        : dimension id of the SR                                        !
! cloud      : fraction of cloudy point                                      !
! clear      : fraction of clear point                                       !
! sat        : fraction of fully attenuated point                            !
! uncer      : fraction of unclassify point                                  !
! nan        : fraction of nan point                                         !
! se         : fraction of se point                                          !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of the altitude boxes                                  !
!----------------------------------------------------------------------------!
! nanb       : Not a Number value                                            !
! ndims      : dimension of dim (= 4 : lon,lat,alt,time)                     !
! ncid       : netcdf file id                                                !
! varid3     : variable id of cloud fraction                                 !
! varid4     : variable id of clear fraction                                 !
! varid5     : variable id of fully attenuated fraction                      !
! varid6     : variable id of unclassify fraction                            !
! varid7     : variable id of nan fraction                                   !
! varid8     : variable id of se fraction                                    !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : prof_recvar2nc(monthcloudfract,monthclearfract,monthsatfract,         !
!                     monthuncertfract,monthnanfract,monthsefract,dimidsp,   ! 
!                     file8,altmax,lonmax,latmax)                            !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine prof_recvar2nc(cloud,clear,uncer,dim,fname,alt,nlon,nlat)!nan,se,sat
    use netcdf
    implicit none

    integer, parameter ::  ndims=4
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  cloud,clear,uncer!,sat,nan,se

    call check(nf90_open('./'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clcalipso', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Cloud fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clrcalipso', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Clear fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

 

   call check(nf90_def_var(ncid, 'uncalipso', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO 3D Undefined fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','1 fraction'))
   call check(nf90_put_att(ncid, varid6, '_FillValue',nan))


    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, cloud))
    call check(nf90_put_var(ncid, varid4, clear))
 
    call check(nf90_put_var(ncid, varid6, uncer))
   
    call check(nf90_close(ncid))

end subroutine prof_recvar2nc
!----------------------------------------------------------------------------!

subroutine temp_recvar2nc(cloud,liq,ice,phase,dim,fname,alt,nlon,nlat)
    use netcdf
    implicit none

    integer, parameter ::  ndims=4
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  cloud,ice,liq,phase

    call check(nf90_open('./'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cltemp', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Cloud fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltemp_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Liquid fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltemp_ice', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO 3D Ice fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltemp_phase', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'CALIPSO 3D Ratio fraction'))
    call check(nf90_put_att(ncid, varid7, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, cloud))
    call check(nf90_put_var(ncid, varid4, liq))
    call check(nf90_put_var(ncid, varid6, ice))
    call check(nf90_put_var(ncid, varid7, phase))

    call check(nf90_close(ncid))

end subroutine temp_recvar2nc
!----------------------------------------------------------------------------!


subroutine depol_recvar2nc(ice,water,dim,fname,alt,nlon,nlat)
    use netcdf
    implicit none

    integer, parameter ::  ndims=4
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  ice,water

    call check(nf90_open('./'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'ice_cloud', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Ice Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'water_cloud', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Water Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

 
    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, ice))
    call check(nf90_put_var(ncid, varid4, water))


    call check(nf90_close(ncid))

end subroutine depol_recvar2nc
!----------------------------------------------------------------------------!

subroutine depol_recvar2ncocc(ice,water,un,phase,ind,dim,dim2,fname,alt,nlon,nlat)
    use netcdf
    implicit none

    integer, parameter ::  ndims=4, ncat=5
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt,dim2(ncat)
    integer  ::  varid3,varid4,varid5,varid6,varid7,varid8, varid9, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  ice,water,phase,ind
    real*4,dimension(nlon,nlat,alt,ncat)  ::  un  

    call check(nf90_open('./'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clcalipso_ice', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Ice Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Liquid Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_un', NF90_FLOAT, dim2, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO 3D UNCLASS Cloud Phase fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_RPIC', NF90_FLOAT, dim, varid9))
    call check(nf90_put_att(ncid, varid9, 'lon_name',                        &
               'CALIPSO 3D Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid9, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid9, '_FillValue',nan))
    call check(nf90_put_att(ncid, varid9, 'Other_Phase',rej))


    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid3, ice))
    call check(nf90_put_var(ncid, varid4, water))
    call check(nf90_put_var(ncid, varid5, un))
    call check(nf90_put_var(ncid, varid9, phase))

    call check(nf90_close(ncid))

end subroutine depol_recvar2ncocc
!----------------------------------------------------------------------------!
subroutine record_ind3d(cloud,tot,ice,water,un,dim,dim2,fname,alt,nlon,nlat)
    use netcdf
    implicit none

    integer, parameter ::  ndims=4, ncat=5
    real,parameter  ::  nan=-9999. , rej=-777. , se=-888.
    integer  ::  nlon , nlat,dim(ndims),alt,dim2(ncat)
    integer  ::  varid2,varid1,varid3,varid4,varid5,varid6,varid7,varid8, varid9, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt)  ::  tot,cloud,ice,water,phase
    real*4,dimension(nlon,nlat,alt,ncat)  ::  un  

    call check(nf90_open('./'//fname,NF90_WRITE,ncid))
    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clcalipso', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO 3D Cloud'))
    call check(nf90_put_att(ncid, varid2, 'units','1'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'sample', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO 3D Sample'))
    call check(nf90_put_att(ncid, varid1, 'units','1'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_ice', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO 3D Ice Cloud Phase occurrences'))
    call check(nf90_put_att(ncid, varid3, 'units','1'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO 3D Liquid Cloud Phase occurrences'))
    call check(nf90_put_att(ncid, varid4, 'units','1 occurrences'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clcalipso_un', NF90_FLOAT, dim2, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO 3D UNCLASS Cloud Phase occurrences'))
    call check(nf90_put_att(ncid, varid5, 'units','1'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))


    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid2, cloud))
    call check(nf90_put_var(ncid, varid1, tot))
    call check(nf90_put_var(ncid, varid3, ice))
    call check(nf90_put_var(ncid, varid4, water))
    call check(nf90_put_var(ncid, varid5, un))

    call check(nf90_close(ncid))

end subroutine record_ind3d
!----------------------------------------------------------------------------!




!----------------------------------------------------------------------------!
! *** CREATE_MAPNC ***  This routine create a netcdf map file, and its       !
!                       dimensions                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MapLowMidHigh file                      !
! dname      : period of MapLowMidHigh file (description of ncdf file)       !
! grid       : name of the grid selected (cfmip/lmdz/nasa)                   !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! vlon       : values of model longitude from 180 to -180                    !
! vlat       : values of model latitude from 90 to -90                       !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
! dim        : dimension id of the MaPLowMidHigh var recorded in the ncdf    !
!----------------------------------------------------------------------------!
! toplvl     : values of isccp top levels                                    !
!              files                                                         !
! ndims      : dimension of dim (=3 : lon,lat,time)                          !
! nc_id      : netcdf file id                                                !
! lon_varid  : variable id of longitude                                      !
! lat_varid  : variable id of latitude                                       !
! pres_varid : variable id of toplvl                                         !
! time_varid : variable id of time                                           !
! lon_dimid  : dimension id of longitude                                     !
! lat_dimid  : dimension id of latitude                                      !
! pres_dimid : dimension id of toplvl                                        !
! time_dimid : dimension id of time                                          !
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : create_mapnc(file8,file9,lonmod,latmod,resd,dimidsm,gcm,lonmax,latmax)!
!                                                                            !
!----------------------------------------------------------------------------!
subroutine create_mapnc(fname,dname,vlon,vlat,vtime,dim,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*8
    integer, parameter ::  ndims = 3
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(3)  ::  toplvl

if(trim(grid).eq.'CFMIP')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2


elseif(trim(grid).eq.'NASA')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2

!fixme : toplowl -> toplvl(1) etc. (I guess)
! elseif(trim(grid).eq.'LMDZ40')then
! toplowl = 14;       
! topmidl = 18;      
! tophighl = 40  
! elseif(trim(grid).eq.'LMDZ')then
! toplowl = 7;        
! topmidl = 9       
! tophighl = 19   !
endif
  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_Low_Mid_High_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'toplvl', 3, pres_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

   call check(nf90_def_var(nc_id, 'toplvl', NF90_FLOAT,pres_dimid,pres_varid))
   call check(nf90_put_att(nc_id, pres_varid,'lon_name','Top Altitude Level'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, toplvl))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_mapnc
!----------------------------------------------------------------------------!


subroutine create_mapnc_phase(fname,dname,vlon,vlat,vtime,dim,dim2,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*8
    integer, parameter ::  ndims = 3, ndims2 = 4, ncat1=5, ncat2=25
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims), dim2(ndims2),dim3(2)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id,cat_varid
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid, cat_dimid1,cat_dimid2
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(3)  ::  toplvl
    character(len=25),dimension(ncat1,ncat2)  ::  vcat


vcat(1,:)='UNDEFINED'
vcat(2,:)='FALSE LIQ'
vcat(3,:)='FALSE ICE'
vcat(4,:)='Horizontally Oriented'
vcat(5,:)='Unphysical Value (NOISE)'

if(trim(grid).eq.'CFMIP')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2


elseif(trim(grid).eq.'NASA')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2

! fixme
! elseif(trim(grid).eq.'LMDZ40')then
! toplowl = 14;       
! topmidl = 18;      
! tophighl = 40  
! elseif(trim(grid).eq.'LMDZ')then
! toplowl = 7;        
! topmidl = 9       
! tophighl = 19   
endif
  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_Low_Mid_High_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'cat1', ncat1, cat_dimid1))
    call check(nf90_def_dim(nc_id, 'cat2', ncat2, cat_dimid2))

    call check(nf90_def_dim(nc_id, 'toplvl', 3, pres_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)
     dim2 = (/lon_dimid, lat_dimid, cat_dimid1, time_dimid/)
     dim3=(/cat_dimid2,cat_dimid1 /)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim3, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'toplvl', NF90_FLOAT,pres_dimid,pres_varid))
    call check(nf90_put_att(nc_id, pres_varid,'lon_name','Top Altitude Level'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, cat_varid, vcat))
    call check(nf90_put_var(nc_id, pres_varid, toplvl))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_mapnc_phase
!----------------------------------------------------------------------------!



subroutine create_mapnc2(fname,dname,vlon,vlat,vtime,dim,dim2,dim3,dim4,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*5
    integer, parameter ::  ndims = 3, ndims2=4, histmax=10, histmax2=40,histmax3=28
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims), dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid,hist_dimid,hist2_dimid,hist3_dimid
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real,dimension(3)  ::  toplvl

if(trim(grid).eq.'CFMIP')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2


elseif(trim(grid).eq.'NASA')then

toplvl(1)=3.36
toplvl(2)=6.72
toplvl(3)=19.2

! fixme
! elseif(trim(grid).eq.'LMDZ40')then
! toplowl = 14;       
! topmidl = 18;      
! tophighl = 40  
! elseif(trim(grid).eq.'LMDZ')then
! toplowl = 7;        
! topmidl = 9       
! tophighl = 19   !
endif
  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_Low_Mid_High_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'hist', histmax, hist_dimid))
    call check(nf90_def_dim(nc_id, 'hist2', histmax2, hist2_dimid))
    call check(nf90_def_dim(nc_id, 'hist3', histmax3, hist3_dimid))

    call check(nf90_def_dim(nc_id, 'toplvl', 3, pres_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)
     dim2 = (/lon_dimid, lat_dimid, hist_dimid, time_dimid/)
     dim3 = (/lon_dimid, lat_dimid, hist2_dimid, time_dimid/)
     dim4 = (/lon_dimid, lat_dimid, hist3_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

   call check(nf90_def_var(nc_id, 'toplvl', NF90_FLOAT,pres_dimid,pres_varid))
   call check(nf90_put_att(nc_id, pres_varid,'lon_name','Top Altitude Level'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, toplvl))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_mapnc2
!----------------------------------------------------------------------------!




!----------------------------------------------------------------------------!
! *** MAP_RECVAR2NC *** This routine record the map variables in the netcdf  !
!                       prof file                                            !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf MaPLowMidHigh file                      !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! dim        : dimension id of the isccp variables                           !
! low        : isccp low cloud flag                                          !
! mid        : isccp mid cloud flag                                          !
! high       : isccp high cloud flag                                         !
! colcloud   : isccp column cloud flag                                       !
! colclear   : isccp column clear flag                                       !
!----------------------------------------------------------------------------!
! nan        : Not a Number value                                            !
! ndims      : dimension of dim (= 3 : lon,lat,time)                         !
! ncid       : netcdf file id                                                !
! varid1     : variable id of cloud fraction                                 !
! varid2     : variable id of clear fraction                                 !
! varid3     : variable id of fully attenuated fraction                      !
! varid4     : variable id of unclassify fraction                            !
! varid5     : variable id of nan fraction                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : map_recvar2nc(monthisccplow,monthisccpmid,monthisccphigh,monthcolcloud!
!                    monthcolclear,dimidsm,file8,lonmax,latmax)              !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine map_recvar2nc2(low,mid,high,colcloud,colclear, &
                         dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10
    character(LEN=*)  ::  fname

    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, colclear


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

  
    call check(nf90_def_var(ncid, 'cllcalipso', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clccalipso', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Total Clear Fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

 
    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid5, colclear))


    call check(nf90_close(ncid))

end subroutine map_recvar2nc2
!----------------------------------------------------------------------------!

subroutine create_maphighnc(fname,dname,vlon,vlat,vtime,dim,grid,nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5,grid*8
    integer, parameter ::  ndims = 3
    integer,dimension(8)  ::  value
    integer  ::  nlon , nlat,dim(ndims)
    integer  ::  lon_varid,lat_varid,pres_varid, time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,time_dimid
	real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat


  
    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Map_High_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105')) 
    call check(nf90_def_dim(nc_id, 'lon', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'lat', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
 
     dim = (/lon_dimid, lat_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'lon', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'lat', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

     call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
	
    call check(nf90_enddef(nc_id))

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, time_varid, vtime))

    call check(nf90_close(nc_id))

end subroutine create_maphighnc
!----------------------------------------------------------------------------!


subroutine maphigh(high,top,base,dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid3,varid5, varid6, ncid !
    character(LEN=*)  ::  fname

    real*4,dimension(nlon,nlat)  ::  high,top,base


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'clhcalipso', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'topcalipso', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Top Cloud Height'))
    call check(nf90_put_att(ncid, varid5, 'units','Kilometer'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))
 
    call check(nf90_def_var(ncid, 'bascalipso', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Base Cloud Height'))
    call check(nf90_put_att(ncid, varid6, 'units','Kilometer'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid5, top))
    call check(nf90_put_var(ncid, varid6, base))

    call check(nf90_close(ncid))

end subroutine maphigh


subroutine map_recvar2nc2phase(liq,ice,dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10
    character(LEN=*)  ::  fname

    real*4,dimension(nlon,nlat,4)  ::  liq,ice


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

  
    call check(nf90_def_var(ncid, 'cll_liq', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_liq', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_liq', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cll_ice', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_ice', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_ice', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid7, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_ice', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid8, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))


    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid1, liq(:,:,1)))
    call check(nf90_put_var(ncid, varid2, liq(:,:,2)))
    call check(nf90_put_var(ncid, varid3, liq(:,:,3)))
    call check(nf90_put_var(ncid, varid4, liq(:,:,4)))

    call check(nf90_put_var(ncid, varid5, ice(:,:,1)))
    call check(nf90_put_var(ncid, varid6, ice(:,:,2)))
    call check(nf90_put_var(ncid, varid7, ice(:,:,3)))
    call check(nf90_put_var(ncid, varid8, ice(:,:,4)))



    call check(nf90_close(ncid))

end subroutine map_recvar2nc2phase
!----------------------------------------------------------------------------!

subroutine map_recvar2nc2phaseocc(liq,ice,indlow,indmid,indtot,dim,fname,nlon,nlat)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3
    real,parameter  ::  nan=-9999
    integer  ::  nlon , nlat ,dim(ndims)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11
    character(LEN=*)  ::  fname
   
    real*4,dimension(nlon,nlat,4)  ::  liq,ice
    real*4,dimension(nlon,nlat)  ::  indtot,indlow,indmid



    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

  
    call check(nf90_def_var(ncid, 'cll_liq', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_liq', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_liq', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_liq', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cll_ice', NF90_FLOAT, dim, varid5))
    call check(nf90_put_att(ncid, varid5, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid5, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid5, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clm_ice', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid6, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clh_ice', NF90_FLOAT, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid7, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clt_ice', NF90_FLOAT, dim, varid8))
    call check(nf90_put_att(ncid, varid8, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid8, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid8, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'ind_low', NF90_FLOAT, dim, varid9))
    call check(nf90_put_att(ncid, varid9, 'lon_name',                        &
               'indice low'))
    call check(nf90_put_att(ncid, varid9, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid9, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'ind_mid', NF90_FLOAT, dim, varid10))
    call check(nf90_put_att(ncid, varid10, 'lon_name',                        &
               'indice mid'))
    call check(nf90_put_att(ncid, varid10, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid10, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'ind_tot', NF90_FLOAT, dim, varid11))
    call check(nf90_put_att(ncid, varid11, 'lon_name',                        &
               'indice high/tot'))
    call check(nf90_put_att(ncid, varid11, 'units','Occurrence'))
    call check(nf90_put_att(ncid, varid11, '_FillValue',nan))

    call check(nf90_enddef(ncid))


    call check(nf90_put_var(ncid, varid1, liq(:,:,1)))
    call check(nf90_put_var(ncid, varid2, liq(:,:,2)))
    call check(nf90_put_var(ncid, varid3, liq(:,:,3)))
    call check(nf90_put_var(ncid, varid4, liq(:,:,4)))

    call check(nf90_put_var(ncid, varid5, ice(:,:,1)))
    call check(nf90_put_var(ncid, varid6, ice(:,:,2)))
    call check(nf90_put_var(ncid, varid7, ice(:,:,3)))
    call check(nf90_put_var(ncid, varid8, ice(:,:,4)))

    call check(nf90_put_var(ncid, varid9, indlow(:,:)))
    call check(nf90_put_var(ncid, varid10, indmid(:,:)))
    call check(nf90_put_var(ncid, varid11, indtot(:,:)))


    call check(nf90_close(ncid))

end subroutine map_recvar2nc2phaseocc
!----------------------------------------------------------------------------!


subroutine map_recvar2nc2phaseocc2(liq,ice,un2,phase,dim,dim2,fname,nlon,nlat,catmax)

    use netcdf
    implicit none

    integer, parameter ::  ndims=3, ndims2=4
    real,parameter  ::  nan=-9999.
    integer  ::  nlon, nlat, catmax,dim(ndims),dim2(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid21,varid22,varid23,varid24,varid25, varid26,varid27, varid28
    integer  ::  varid8,varid9,varid10,varid11, varid12,varid13,varid14,varid15,varid16
    character(LEN=*)  ::  fname
    real*4,dimension(nlon,nlat,4)  ::  ice,liq,phase
    real*4,dimension(nlon,nlat,4,catmax)  ::  un2


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

 
    call check(nf90_def_var(ncid, 'cllcalipso_liq', NF90_FLOAT, dim, varid21))
    call check(nf90_put_att(ncid, varid21, 'lon_name',                        &
               'Low Level liquid Cloud'))
    call check(nf90_put_att(ncid, varid21, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid21, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_liq', NF90_FLOAT, dim, varid22))
    call check(nf90_put_att(ncid, varid22, 'lon_name',                        &
               'Middle Level liq Cloud'))
    call check(nf90_put_att(ncid, varid22, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid22, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_liq', NF90_FLOAT, dim, varid23))
    call check(nf90_put_att(ncid, varid23, 'lon_name',                        &
               'High Level liq Cloud'))
    call check(nf90_put_att(ncid, varid23, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid23, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_liq', NF90_FLOAT, dim, varid24))
    call check(nf90_put_att(ncid, varid24, 'lon_name',                        &
               'Total Column Level liq Cloud'))
    call check(nf90_put_att(ncid, varid24, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid24, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cllcalipso_ice', NF90_FLOAT, dim, varid25))
    call check(nf90_put_att(ncid, varid25, 'lon_name',                        &
               'Low-level ice Cloud'))
    call check(nf90_put_att(ncid, varid25, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid25, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_ice', NF90_FLOAT, dim, varid26))
    call check(nf90_put_att(ncid, varid26, 'lon_name',                        &
               'CALIPSO Mid-level ice Cloud'))
    call check(nf90_put_att(ncid, varid26, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid26, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_ice', NF90_FLOAT, dim, varid27))
    call check(nf90_put_att(ncid, varid27, 'lon_name',                        &
               'CALIPSO High-level ice Cloud'))
    call check(nf90_put_att(ncid, varid27, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid27, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_ice', NF90_FLOAT, dim, varid28))
    call check(nf90_put_att(ncid, varid28, 'lon_name',                        &
               'CALIPSO Total ice Cloud'))
    call check(nf90_put_att(ncid, varid28, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid28, '_FillValue',nan))

    print *,'step2'

 
    call check(nf90_def_var(ncid, 'cllcalipso_un', NF90_FLOAT, dim2, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'Low Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid1, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_un', NF90_FLOAT, dim2, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'Middle Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid2, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_un', NF90_FLOAT, dim2, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'High Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid3, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_un', NF90_FLOAT, dim2, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'Total Column Level unclassified Cloud'))
    call check(nf90_put_att(ncid, varid4, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cllcalipso_RPIC', NF90_FLOAT, dim, varid13))
    call check(nf90_put_att(ncid, varid13, 'lon_name',                        &
               'CALIPSO Low-level Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid13, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid13, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clmcalipso_RPIC', NF90_FLOAT, dim, varid14))
    call check(nf90_put_att(ncid, varid14, 'lon_name',                        &
               'CALIPSO Mid-level Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid14, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid14, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'clhcalipso_RPIC', NF90_FLOAT, dim, varid15))
    call check(nf90_put_att(ncid, varid15, 'lon_name',                        &
               'CALIPSO High-level Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid15, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid15, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'cltcalipso_RPIC', NF90_FLOAT, dim, varid16))
    call check(nf90_put_att(ncid, varid16, 'lon_name',                        &
               'CALIPSO Total Relative Percentage of Ice in Cloud'))
    call check(nf90_put_att(ncid, varid16, 'units','Fraction'))
    call check(nf90_put_att(ncid, varid16, '_FillValue',nan))


    call check(nf90_enddef(ncid))

    print *,'step3'


    call check(nf90_put_var(ncid, varid21, liq(:,:,1)))
    call check(nf90_put_var(ncid, varid22, liq(:,:,2)))
    call check(nf90_put_var(ncid, varid23, liq(:,:,3)))
    call check(nf90_put_var(ncid, varid24, liq(:,:,4)))

    call check(nf90_put_var(ncid, varid25, ice(:,:,1)))
    call check(nf90_put_var(ncid, varid26, ice(:,:,2)))
    call check(nf90_put_var(ncid, varid27, ice(:,:,3)))
    call check(nf90_put_var(ncid, varid28, ice(:,:,4)))

    call check(nf90_put_var(ncid, varid1, un2(:,:,1,:)))
    call check(nf90_put_var(ncid, varid2, un2(:,:,2,:)))
    call check(nf90_put_var(ncid, varid3, un2(:,:,3,:)))
    call check(nf90_put_var(ncid, varid4, un2(:,:,4,:)))

    call check(nf90_put_var(ncid, varid13, phase(:,:,1)))
    call check(nf90_put_var(ncid, varid14, phase(:,:,2)))
    call check(nf90_put_var(ncid, varid15, phase(:,:,3)))
    call check(nf90_put_var(ncid, varid16, phase(:,:,4)))

 

    call check(nf90_close(ncid))

    print *,'done with this shit'

end subroutine map_recvar2nc2phaseocc2
!----------------------------------------------------------------------------!


subroutine map_recvar2nc3(low,mid,high,colcloud,colclear,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,dim,dim2,dim3,fname,nlon,nlat)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    integer  ::  ilon, ilat
    character(LEN=*)  ::  fname

    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, colclear, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    real*4,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod


histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo

do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif


enddo
enddo


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'Number of orbit passages'))
    call check(nf90_put_att(ncid, varid7, 'units','Arbitrary Unit'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))


    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Cloud Top Height level'))
    call check(nf90_put_att(ncid, varid6, 'units','1-40'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
    call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
    call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
    call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
    call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
    call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

    call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
    call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
    call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
    call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
    call check(nf90_def_var(ncid, 'h_CZ', NF90_FLOAT, dim3, varid17))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))

    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
 call check(nf90_put_var(ncid, varid9, f_CAM)) 
 call check(nf90_put_var(ncid, varid10, f_CAH)) 
 call check(nf90_put_var(ncid, varid11, f_CA)) 
 call check(nf90_put_var(ncid, varid12, f_CZ))
 call check(nf90_put_var(ncid, varid13, hlow ))
 call check(nf90_put_var(ncid, varid14, hmid )) 
 call check(nf90_put_var(ncid, varid15, hhigh)) 
 call check(nf90_put_var(ncid, varid16, hcol)) 
 call check(nf90_put_var(ncid, varid17, hheight))
 

    call check(nf90_close(ncid))

end subroutine map_recvar2nc3
!----------------------------------------------------------------------------!

subroutine map_recvar2nc7(low,mid,high,colcloud,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,dim,dim2,dim3,dim4,fname,nlon,nlat, &
                         lowtemp,midtemp,hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
                         hcoltemp)
! h_CA utile car 31 v!aleur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4, histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    character(LEN=*)  ::  fname

    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH

    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    real*4,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp

    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp
    integer  ::  varid27, varid28, varid29, varid30, varid18, varid19, varid20, varid21
    integer  ::  varid31, varid32, varid33, varid34
    integer  ::  ilon, ilat


histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo


do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif

if (lowtemp(ilon,ilat).gt.0)then
f_CTL(ilon,ilat)=100
else
f_CTL(ilon,ilat)=0
endif

if (midtemp(ilon,ilat).gt.0)then
f_CTM(ilon,ilat)=100
else
f_CTM(ilon,ilat)=0
endif

if (hightemp(ilon,ilat).gt.0)then
f_CTH(ilon,ilat)=100
else
f_CTH(ilon,ilat)=0
endif


if (coltemp(ilon,ilat).gt.0)then
f_CT(ilon,ilat)=100
else
f_CT(ilon,ilat)=0
endif

enddo
enddo


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'Number of orbit passages'))
    call check(nf90_put_att(ncid, varid7, 'units','Arbitrary Unit'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))


    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Cloud Top Height level'))
    call check(nf90_put_att(ncid, varid6, 'units','1-40'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
    call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
    call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
    call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
    call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
    call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

    call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
    call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
    call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
    call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
    call check(nf90_def_var(ncid, 'h_CZ', NF90_FLOAT, dim3, varid17))

    call check(nf90_def_var(ncid, 'a_CTL', NF90_FLOAT, dim, varid27))
    call check(nf90_def_var(ncid, 'a_CTM', NF90_FLOAT, dim, varid28))
    call check(nf90_def_var(ncid, 'a_CTH', NF90_FLOAT, dim, varid29))
    call check(nf90_def_var(ncid, 'a_CT', NF90_FLOAT, dim, varid30))

    call check(nf90_def_var(ncid, 'f_CTL', NF90_INT4, dim, varid18))
    call check(nf90_def_var(ncid, 'f_CTM', NF90_INT4, dim, varid19))
    call check(nf90_def_var(ncid, 'f_CTH', NF90_INT4, dim, varid20))
    call check(nf90_def_var(ncid, 'f_CT', NF90_INT4, dim, varid21))

    call check(nf90_def_var(ncid, 'h_CTL', NF90_INT4, dim4, varid31))
    call check(nf90_def_var(ncid, 'h_CTM', NF90_INT4, dim4, varid32))
    call check(nf90_def_var(ncid, 'h_CTH', NF90_INT4, dim4, varid33))
    call check(nf90_def_var(ncid, 'h_CT', NF90_INT4, dim4, varid34))


    call check(nf90_enddef(ncid))

  
    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
    call check(nf90_put_var(ncid, varid9, f_CAM)) 
    call check(nf90_put_var(ncid, varid10, f_CAH)) 
    call check(nf90_put_var(ncid, varid11, f_CA)) 
    call check(nf90_put_var(ncid, varid12, f_CZ))
    call check(nf90_put_var(ncid, varid13, hlow ))
    call check(nf90_put_var(ncid, varid14, hmid )) 
    call check(nf90_put_var(ncid, varid15, hhigh)) 
    call check(nf90_put_var(ncid, varid16, hcol)) 
    call check(nf90_put_var(ncid, varid17, hheight))
    call check(nf90_put_var(ncid, varid18, f_CTL)) 
    call check(nf90_put_var(ncid, varid19, f_CTM)) 
    call check(nf90_put_var(ncid, varid20, f_CTH)) 
    call check(nf90_put_var(ncid, varid21, f_CT)) 
    call check(nf90_put_var(ncid, varid27, lowtemp))
    call check(nf90_put_var(ncid, varid28, midtemp))
    call check(nf90_put_var(ncid, varid29, hightemp))
    call check(nf90_put_var(ncid, varid30, coltemp))
    call check(nf90_put_var(ncid, varid31, hlowtemp))
    call check(nf90_put_var(ncid, varid32, hmidtemp))
    call check(nf90_put_var(ncid, varid33, hhightemp))
    call check(nf90_put_var(ncid, varid34, hcoltemp))

 

    call check(nf90_close(ncid))

end subroutine map_recvar2nc7



subroutine map_recvar2nc6(low,mid,high,colcloud,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,dim,dim2,dim3,fname,nlon,nlat)!,&

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4,histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    integer  ::  ilon, ilat
    character(LEN=*)  ::  fname
   
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    integer,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod

   
    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp
    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp
    integer,dimension(histmax3)  ::  histmod3





histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo




do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif


enddo
enddo


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_put_att(ncid, varid7, 'lon_name',                        &
               'Number of orbit passages'))
    call check(nf90_put_att(ncid, varid7, 'units','Arbitrary Unit'))
    call check(nf90_put_att(ncid, varid7, '_FillValue',nan))


    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_put_att(ncid, varid1, 'lon_name',                        &
               'CALIPSO Low-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid1, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid1, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                        &
               'CALIPSO Mid-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid2, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                        &
               'CALIPSO High-level Cloud Fraction'))
    call check(nf90_put_att(ncid, varid3, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_put_att(ncid, varid4, 'lon_name',                        &
               'CALIPSO Total Cloud Fraction'))
    call check(nf90_put_att(ncid, varid4, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid4, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
    call check(nf90_put_att(ncid, varid6, 'lon_name',                        &
               'CALIPSO Cloud Top Height level'))
    call check(nf90_put_att(ncid, varid6, 'units','1-40'))
    call check(nf90_put_att(ncid, varid6, '_FillValue',nan))
 
 call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
 call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
 call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
 call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
 call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

 call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
 call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
 call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
 call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
 call check(nf90_def_var(ncid, 'h_CZ', NF90_INT4, dim3, varid17))

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
 call check(nf90_put_var(ncid, varid9, f_CAM)) 
 call check(nf90_put_var(ncid, varid10, f_CAH)) 
 call check(nf90_put_var(ncid, varid11, f_CA)) 
 call check(nf90_put_var(ncid, varid12, f_CZ))
 call check(nf90_put_var(ncid, varid13, hlow ))
 call check(nf90_put_var(ncid, varid14, hmid )) 
 call check(nf90_put_var(ncid, varid15, hhigh)) 
 call check(nf90_put_var(ncid, varid16, hcol)) 
 call check(nf90_put_var(ncid, varid17, hheight))
 

    call check(nf90_close(ncid))

end subroutine map_recvar2nc6
!----------------------------------------------------------------------------!











subroutine map_recvar2nc4(low,mid,high,colcloud,colclear,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,lowtemp,midtemp, &
                         hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
                         hcoltemp,dim,dim2,dim3,dim4,fname,nlon,nlat)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4,histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17,varid18,varid19,varid20,varid21,varid22,varid23,varid24,varid25,varid26,varid27,varid28,varid29,varid30
    integer  ::  ilon, ilat
    character(LEN=*)  ::  fname
    
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, colclear, height
    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp

    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    integer,dimension(nlon,nlat,histmax2)  ::  hheight   
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp

    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod
    integer,dimension(histmax3)  ::  histmod3


histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo



print *, 'f_ begin'

do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif

if (lowtemp(ilon,ilat).gt.0)then
f_CTL(ilon,ilat)=100
else
f_CTL(ilon,ilat)=0
endif

if (midtemp(ilon,ilat).gt.0)then
f_CTM(ilon,ilat)=100
else
f_CTM(ilon,ilat)=0
endif

if (hightemp(ilon,ilat).gt.0)then
f_CTH(ilon,ilat)=100
else
f_CTH(ilon,ilat)=0
endif


if (coltemp(ilon,ilat).gt.0)then
f_CT(ilon,ilat)=100
else
f_CT(ilon,ilat)=0
endif


enddo
enddo

print *, 'f_ ok'


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

     call check(nf90_def_var(ncid, 'a_CTL', NF90_FLOAT, dim, varid27))
     call check(nf90_def_var(ncid, 'a_CTM', NF90_FLOAT, dim, varid28))
     call check(nf90_def_var(ncid, 'a_CTH', NF90_FLOAT, dim, varid29))
     call check(nf90_def_var(ncid, 'a_CT', NF90_FLOAT, dim, varid30))

 call check(nf90_def_var(ncid, 'f_CTL', NF90_INT4, dim, varid18))
 call check(nf90_def_var(ncid, 'f_CTM', NF90_INT4, dim, varid19))
 call check(nf90_def_var(ncid, 'f_CTH', NF90_INT4, dim, varid20))
 call check(nf90_def_var(ncid, 'f_CT', NF90_INT4, dim, varid21))

 call check(nf90_def_var(ncid, 'h_CTL', NF90_INT4, dim4, varid23))
 call check(nf90_def_var(ncid, 'h_CTM', NF90_INT4, dim4, varid24))
 call check(nf90_def_var(ncid, 'h_CTH', NF90_INT4, dim4, varid25))
 call check(nf90_def_var(ncid, 'h_CT', NF90_INT4, dim4, varid26))


print *, 'def var ok'

    call check(nf90_enddef(ncid))
print *, 'f_ ok'

 call check(nf90_put_var(ncid, varid18, f_CTL)) 
 call check(nf90_put_var(ncid, varid19, f_CTM)) 
 call check(nf90_put_var(ncid, varid20, f_CTH)) 
 call check(nf90_put_var(ncid, varid21, f_CT)) 
print *, 'f_ ok'

 call check(nf90_put_var(ncid, varid23, hlowtemp ))
 call check(nf90_put_var(ncid, varid24, hmidtemp )) 
 call check(nf90_put_var(ncid, varid25, hhightemp)) 
 call check(nf90_put_var(ncid, varid26, hcoltemp)) 
print *, 'f_ ok'

    call check(nf90_put_var(ncid, varid27, lowtemp))
    call check(nf90_put_var(ncid, varid28, midtemp))
    call check(nf90_put_var(ncid, varid29, hightemp))
    call check(nf90_put_var(ncid, varid30, coltemp))



print *, 'f_ ok'

 

    call check(nf90_close(ncid))

end subroutine map_recvar2nc4
!----------------------------------------------------------------------------!


subroutine map_recvar2nc5(low,mid,high,colcloud,height,indtot,&
                         hlow,hmid,hhigh,hcol,hheight,lowtemp,midtemp, &
                         hightemp,coltemp,hlowtemp,hmidtemp,hhightemp, &
                         hcoltemp,dim,dim2,dim3,dim4,fname,nlon,nlat)

! h_CA utile car 31 valeur au max...
    use netcdf
    implicit none

    integer, parameter ::  ndims=3,histmax=10,histmax2=40,ndims2=4,histmax3=28
    real,parameter  ::  nan=-999
    integer  ::  nlon , nlat ,dim(ndims),ihist,dim2(ndims2),dim3(ndims2),dim4(ndims2)
    integer  ::  varid1,varid2,varid3,varid4,varid5, varid6,varid7, ncid !
    integer  ::  varid8,varid9,varid10,varid11,varid12,varid13,varid14
    integer  ::  varid15,varid16,varid17
    integer  ::  ilon, ilat
    character(LEN=*)  ::  fname
   
    real*4,dimension(nlon,nlat)  ::  low, mid, high
    real*4,dimension(nlon,nlat)  ::  colcloud, height
    integer,dimension(nlon,nlat)  ::  indtot,f_CA,f_CAL,f_CAM,f_CAH,f_CZ
    integer,dimension(nlon,nlat,histmax)  ::  hlow,hmid,hhigh,hcol
    integer,dimension(nlon,nlat,histmax2)  ::  hheight    
    integer,dimension(histmax2)  ::  histmod2
    real,dimension(histmax)  ::  histmod

    integer  ::  varid18,varid19,varid20,varid21,varid22,varid23,varid24,varid25,varid26,varid27,varid28,varid29,varid30

    real*4,dimension(nlon,nlat)  ::  lowtemp,midtemp,hightemp,coltemp

    integer,dimension(nlon,nlat)  ::  f_CT,f_CTL,f_CTM,f_CTH
    integer,dimension(nlon,nlat,histmax3)  ::  hlowtemp,hmidtemp,hhightemp,hcoltemp

histmod(:)=0;
do ihist=1,histmax
   histmod(ihist+1)=histmod(ihist)+0.1
enddo

histmod2(:)=0;
do ihist=1,histmax2
   histmod2(ihist+1)=histmod2(ihist)+1
enddo




do ilon=1,nlon
  do ilat=1,nlat

if (colcloud(ilon,ilat).eq.0)then
f_CA(ilon,ilat)=0
elseif(colcloud(ilon,ilat).eq.-999)then
f_CA(ilon,ilat)=-999
else
f_CA(ilon,ilat)=100
endif

if (low(ilon,ilat).eq.0)then
f_CAL(ilon,ilat)=0
elseif(low(ilon,ilat).eq.-999)then
f_CAL(ilon,ilat)=-999
else
f_CAL(ilon,ilat)=100
endif

if (mid(ilon,ilat).eq.0)then
f_CAM(ilon,ilat)=0
elseif(mid(ilon,ilat).eq.-999)then
f_CAM(ilon,ilat)=-999
else
f_CAM(ilon,ilat)=100
endif

if (high(ilon,ilat).eq.0)then
f_CAH(ilon,ilat)=0
elseif(high(ilon,ilat).eq.-999)then
f_CAH(ilon,ilat)=-999
else
f_CAH(ilon,ilat)=100
endif

if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif


if (height(ilon,ilat).eq.0)then
f_CZ(ilon,ilat)=0
elseif(height(ilon,ilat).eq.-999)then
f_CZ(ilon,ilat)=-999
else
f_CZ(ilon,ilat)=100
endif

if (lowtemp(ilon,ilat).gt.0)then
f_CTL(ilon,ilat)=100
else
f_CTL(ilon,ilat)=0
endif

if (midtemp(ilon,ilat).gt.0)then
f_CTM(ilon,ilat)=100
else
f_CTM(ilon,ilat)=0
endif

if (hightemp(ilon,ilat).gt.0)then
f_CTH(ilon,ilat)=100
else
f_CTH(ilon,ilat)=0
endif


if (coltemp(ilon,ilat).gt.0)then
f_CT(ilon,ilat)=100
else
f_CT(ilon,ilat)=0
endif


enddo
enddo

print *, 'creation fichier'

    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'n_tot', NF90_INT4, dim, varid7))
    call check(nf90_def_var(ncid, 'a_CAL', NF90_FLOAT, dim, varid1))
    call check(nf90_def_var(ncid, 'a_CAM', NF90_FLOAT, dim, varid2))
    call check(nf90_def_var(ncid, 'a_CAH', NF90_FLOAT, dim, varid3))
    call check(nf90_def_var(ncid, 'a_CA', NF90_FLOAT, dim, varid4))
    call check(nf90_def_var(ncid, 'a_CZ', NF90_FLOAT, dim, varid6))
 
 call check(nf90_def_var(ncid, 'f_CAL', NF90_INT4, dim, varid8))
 call check(nf90_def_var(ncid, 'f_CAM', NF90_INT4, dim, varid9))
 call check(nf90_def_var(ncid, 'f_CAH', NF90_INT4, dim, varid10))
 call check(nf90_def_var(ncid, 'f_CA', NF90_INT4, dim, varid11))
 call check(nf90_def_var(ncid, 'f_CZ', NF90_INT4, dim, varid12))

 call check(nf90_def_var(ncid, 'h_CAL', NF90_INT4, dim2, varid13))
 call check(nf90_def_var(ncid, 'h_CAM', NF90_INT4, dim2, varid14))
 call check(nf90_def_var(ncid, 'h_CAH', NF90_INT4, dim2, varid15))
 call check(nf90_def_var(ncid, 'h_CA', NF90_INT4, dim2, varid16))
 call check(nf90_def_var(ncid, 'h_CZ', NF90_INT4, dim3, varid17))


    call check(nf90_enddef(ncid))
print *, 'end def'

    call check(nf90_put_var(ncid, varid1, low))
    call check(nf90_put_var(ncid, varid2, mid))
    call check(nf90_put_var(ncid, varid3, high))
print *, 'file gewex1'

    call check(nf90_put_var(ncid, varid4, colcloud))
    call check(nf90_put_var(ncid, varid6, height)) 
    call check(nf90_put_var(ncid, varid7, indtot)) 
    call check(nf90_put_var(ncid, varid8, f_CAL)) 
print *, 'file gewex1'

 call check(nf90_put_var(ncid, varid9, f_CAM)) 
 call check(nf90_put_var(ncid, varid10, f_CAH)) 
 call check(nf90_put_var(ncid, varid11, f_CA)) 
 call check(nf90_put_var(ncid, varid12, f_CZ))
print *, 'file gewex1'

 call check(nf90_put_var(ncid, varid13, hlow ))
 call check(nf90_put_var(ncid, varid14, hmid )) 
 call check(nf90_put_var(ncid, varid15, hhigh)) 
 call check(nf90_put_var(ncid, varid16, hcol)) 
print *, 'file gewex1'

 call check(nf90_put_var(ncid, varid17, hheight))
 
print *, 'file gewex1'


    call check(nf90_close(ncid))

end subroutine map_recvar2nc5
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** CREATE_MAPNC ***  This routine create a netcdf diag file, and its      !
!                       dimensions                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf diagSR file                             !
! dname      : period of diagSR file (description of ncdf file)              !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! vlon       : values of model longitude from 180 to -180                    !
! vlat       : values of model latitude from 90 to -90                       !
! vprestop   : values of model altitude                                      !
! vsrmod     : values of model diag boxes                                    !
! alt        : values of model altitude                                      !
! dim        : dimension id of the diagSR var recorded in the ncdf           !
!              files                                                         !
! vtime      : number of days since 2000/01/01 for the trimonthly perdiod,   !
!              in day.                                                       !
!----------------------------------------------------------------------------!
! date       : date from the real-system clock and has form yyyymmdd.        !
! time       : time from the real-system clock and has form hhmmss.sss.      !
! zone       : represente the difference with respect to Coordinated         !
!              Universal Time (UTC), and has form (+-)hhmm.                  !
! value      : 8 dimension value which contains the year,month,day,hour,     !
!              minute, seconds and milliseconds of the real-time.            !
! diagmax    : values of SR diagbox                                          !
! ndims      : dimension of dim (=5 : lon,lat,alt,diag,time)                 !
! nc_id      : netcdf file id                                                !
! lon_varid  : variable id of longitude                                      !
! lat_varid  : variable id of latitude                                       !
! pres_varid : variable id of altitude                                       !
! srmod_varid: variable id of diagsr                                         !
! time_varid : variable id of time                                           !
! lon_dimid  : dimension id of longitude                                     !
! lat_dimid  : dimension id of latitude                                      !
! pres_dimid : dimension id of altitude                                      !
! srmod_dimid: dimension id of srbox                                         !
! time_dimid : dimension id of time                                          !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : create_diagnc(file8,file9,lonmod,latmod,altmod,srmod,resd,dimidsd,    !
!                    altmax,lonmax,latmax)                                   !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine create_diagnc(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vsrmod,vtime,dim,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon, nlat, iz
    integer, parameter :: diagmax = 19, ndims = 5, nv=2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),dim2(ndims),alt,dim3(nv),dim4(nv),dim5(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, srmod_varid,time_varid,srmod_varid2,srmod_varid3,srmod_varid4,nc_id,alt_varid2,srmod_varid33,srmod_varid44
    integer  ::  lon_dimid,lat_dimid,alt_dimid,srmod_dimid,time_dimid,srmod_dimid2,srmod_dimid3,srmod_dimid4, nv_dimid,sr_dimid,sr_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    real, dimension(diagmax-1)  ::  vsrmod2
    real, dimension(diagmax-1,2)  ::  vsrmod_bound


vsrmod2(:)=0;
vsrmod_bound(:,:)=0;

   
vsrmod2(1)=-888
vsrmod2(2)=-777

      do iz=3,diagmax-1
            vsrmod2(iz) = (vsrmod(iz)+vsrmod(iz+1))/2
       enddo

vsrmod_bound(1,:)=-888
vsrmod_bound(2,:)=-777
vsrmod_bound(3,1)=-776
vsrmod_bound(3,2)=0

      do iz=4,diagmax-1
         vsrmod_bound(iz,1)=vsrmod(iz)
         vsrmod_bound(iz,2)=vsrmod(iz+1)
      enddo



    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Scattering_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    call check(nf90_def_dim(nc_id, 'box', diagmax-4, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'box2', 3, srmod_dimid4 ))   

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, alt_dimid, srmod_dimid2, time_dimid/)
    dim2 = (/lon_dimid, lat_dimid, alt_dimid, srmod_dimid4, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/srmod_dimid2, nv_dimid/)
    dim5 = (/srmod_dimid4, nv_dimid/)



    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'srbox_mid', NF90_FLOAT, srmod_dimid2,         &
               srmod_varid2))
    call check(nf90_put_att(nc_id, srmod_varid2, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'axis','I')) 
  
    call check(nf90_def_var(nc_id, 'srbox_mid2', NF90_FLOAT, srmod_dimid4,         &
               srmod_varid4))
    call check(nf90_put_att(nc_id, srmod_varid4, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid4, 'units','Arbitrary unit'))
 
    call check(nf90_def_var(nc_id, 'srbox_bound', NF90_FLOAT, dim4,         &
               srmod_varid33))
    call check(nf90_put_att(nc_id, srmod_varid33, 'lon_name',                  &
               'Boundarie Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid33, 'units','Arbitrary unit'))

    call check(nf90_def_var(nc_id, 'srbox_bound2', NF90_FLOAT, dim5,         &
               srmod_varid44))
    call check(nf90_put_att(nc_id, srmod_varid44, 'lon_name',                  &
               'Boundarie Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid44, 'units','Arbitrary unit'))
 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  



    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))

    call check(nf90_put_var(nc_id, srmod_varid2, vsrmod2(4:18)))
    call check(nf90_put_var(nc_id, srmod_varid4, vsrmod2(1:3)))
    call check(nf90_put_var(nc_id, srmod_varid33, vsrmod_bound(4:18,:)))
    call check(nf90_put_var(nc_id, srmod_varid44, vsrmod_bound(1:3,:)))



    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagnc
!----------------------------------------------------------------------------!

subroutine create_diagncpha(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vsrmod,vtime,dim,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon, nlat, iz 
    integer, parameter :: diagmax = 19, ndims = 5, nv=2
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),dim2(ndims),alt,dim3(nv),dim4(nv),dim5(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, srmod_varid,time_varid,srmod_varid2,srmod_varid3,srmod_varid4,nc_id,alt_varid2,srmod_varid33,srmod_varid44
    integer  ::  lon_dimid,lat_dimid,alt_dimid,srmod_dimid,time_dimid,srmod_dimid2,srmod_dimid3,srmod_dimid4, nv_dimid,sr_dimid,sr_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    real, dimension(diagmax-8)  ::  vsrmod2
    real, dimension(diagmax-8,2)  ::  vsrmod_bound


vsrmod2(:)=0;
vsrmod_bound(:,:)=0;

      do iz=8,diagmax-1
            vsrmod2(iz-7) = (vsrmod(iz)+vsrmod(iz+1))/2
       enddo

      do iz=8,diagmax-1
         vsrmod_bound(iz-7,1)=vsrmod(iz)
         vsrmod_bound(iz-7,2)=vsrmod(iz+1)
      enddo

print *, "srmod creer"

    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Phase_Histogram_of_Scattering_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    call check(nf90_def_dim(nc_id, 'box', diagmax-8, srmod_dimid2 ))   

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, alt_dimid, srmod_dimid2, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/srmod_dimid2, nv_dimid/)



    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'srbox_mid', NF90_FLOAT, srmod_dimid2,         &
               srmod_varid2))
    call check(nf90_put_att(nc_id, srmod_varid2, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid2, 'axis','I')) 
  
    call check(nf90_def_var(nc_id, 'srbox_bound', NF90_FLOAT, dim4,         &
               srmod_varid33))
    call check(nf90_put_att(nc_id, srmod_varid33, 'lon_name',                  &
               'Boundarie Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, srmod_varid33, 'units','Arbitrary unit'))

 

    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  



    call check(nf90_enddef(nc_id))
	

    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))


   call check(nf90_put_var(nc_id, srmod_varid2, vsrmod2))
    call check(nf90_put_var(nc_id, srmod_varid33, vsrmod_bound))



    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagncpha
!----------------------------------------------------------------------------!

subroutine create_diagnc2(fname,dname,vlon,vlat,vprestop,vsrmod,vtime,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 19, ndims = 4
    integer,dimension(8)  ::  value
    integer  ::  dim2(ndims),alt
    integer  ::  lon_varid,lat_varid,pres_varid, srmod_varid,time_varid,nc_id
    integer  ::  lon_dimid,lat_dimid,pres_dimid,srmod_dimid,time_dimid,srmod_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real,dimension(alt)  ::  vprestop
    real, dimension(diagmax-1)  ::  vsrmod2


vsrmod2=vsrmod(2:19)
    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Scattering_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, pres_dimid))
    call check(nf90_def_dim(nc_id, 'box2', diagmax-1, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'box', diagmax, srmod_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
  
    dim2 = (/lon_dimid, lat_dimid, pres_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, pres_dimid, pres_varid))
    call check(nf90_put_att(nc_id, pres_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

    call check(nf90_def_var(nc_id, 'srbox', NF90_FLOAT, srmod_dimid,         &
               srmod_varid))
    call check(nf90_put_att(nc_id, srmod_varid, 'lon_name',                  &
               'Scattering Ratio Value'))
    call check(nf90_put_att(nc_id, srmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid, 'axis','I'))  

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, vprestop))
    call check(nf90_put_var(nc_id, srmod_varid, vsrmod))
    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagnc2
!----------------------------------------------------------------------------!


subroutine create_depolnc2(fname,dname,vlon,vlat,vprestop,vsrmod,vdepolmod,vtime,dim2,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 17, ndims = 5, depolmax = 21
    integer,dimension(8)  ::  value
    integer  ::  dim2(ndims),alt
    integer  ::  lon_varid,lat_varid,pres_varid, srmod_varid,time_varid,nc_id, depolmod_varid
    integer  ::  lon_dimid,lat_dimid,pres_dimid,srmod_dimid,time_dimid,srmod_dimid2, depolmod_dimid, depolmod_dimid2
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real, dimension(depolmax)  ::  vdepolmod
    real,dimension(alt)  ::  vprestop
    real, dimension(diagmax-1)  ::  vsrmod2
    real, dimension(depolmax-1)  ::  vdepolmod2


vsrmod2=vsrmod(2:17)
vdepolmod2=vdepolmod(2:21)

    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Depolarization_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, pres_dimid))
    call check(nf90_def_dim(nc_id, 'srbox2', diagmax-1, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'srbox', diagmax, srmod_dimid))
    call check(nf90_def_dim(nc_id, 'depolbox2', depolmax-1, depolmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'depolbox', depolmax, depolmod_dimid))


    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
  
    dim2 = (/lon_dimid, lat_dimid, pres_dimid, srmod_dimid, time_dimid/)

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, pres_dimid, pres_varid))
    call check(nf90_put_att(nc_id, pres_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

    call check(nf90_def_var(nc_id, 'srbox', NF90_FLOAT, srmod_dimid,         &
               srmod_varid))
    call check(nf90_put_att(nc_id, srmod_varid, 'lon_name',                  &
               'Scattering Ratio Value'))
    call check(nf90_put_att(nc_id, srmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid, 'axis','I'))  

    call check(nf90_def_var(nc_id, 'depolbox', NF90_FLOAT, depolmod_dimid,         &
               depolmod_varid))
    call check(nf90_put_att(nc_id, depolmod_varid, 'lon_name',                  &
               'Depolarization Ratio Value'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'axis','I'))  


   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, vprestop))
    call check(nf90_put_var(nc_id, srmod_varid, vsrmod))
    call check(nf90_put_var(nc_id, depolmod_varid, vdepolmod))

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_depolnc2
!----------------------------------------------------------------------------!


subroutine create_depolnc(fname,dname,vlon,vlat,vprestop,vsrmod,vdepolmod,vtime,dim,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat 
    integer, parameter :: diagmax = 17, ndims = 6, depolmax = 21
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),dim2(ndims-1),alt
    integer  ::  lon_varid,lat_varid,pres_varid, srmod_varid,time_varid,nc_id, depolmod_varid
    integer  ::  lon_dimid,lat_dimid,pres_dimid,srmod_dimid,time_dimid,srmod_dimid2, depolmod_dimid2, depolmod_dimid
    real  ::  vtime  
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(diagmax)  ::  vsrmod
    real, dimension(depolmax)  ::  vdepolmod
    real,dimension(alt)  ::  vprestop
    real, dimension(diagmax-1)  ::  vsrmod2
    real, dimension(depolmax-1)  ::  vdepolmod2



vsrmod2=vsrmod(2:17)
vdepolmod2=vdepolmod(2:21)

    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Depolarization_Ratio'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, pres_dimid))
    call check(nf90_def_dim(nc_id, 'srbox2', diagmax-1, srmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'srbox', diagmax, srmod_dimid))
    call check(nf90_def_dim(nc_id, 'depolbox2', depolmax-1, depolmod_dimid2 ))   
    call check(nf90_def_dim(nc_id, 'depolbox', depolmax, depolmod_dimid))
    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, pres_dimid, srmod_dimid2, depolmod_dimid2, time_dimid/)
    
  

    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'altitude', NF90_FLOAT, pres_dimid, pres_varid))
    call check(nf90_put_att(nc_id, pres_varid, 'lon_name','Altitude'))
    call check(nf90_put_att(nc_id, pres_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, pres_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, pres_varid, 'axis','Z'))   

    call check(nf90_def_var(nc_id, 'srbox', NF90_FLOAT, srmod_dimid,         &
               srmod_varid))
    call check(nf90_put_att(nc_id, srmod_varid, 'lon_name',                  &
               'Scattering Ratio Value'))
    call check(nf90_put_att(nc_id, srmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, srmod_varid, 'axis','I'))  

    call check(nf90_def_var(nc_id, 'depolbox', NF90_FLOAT, depolmod_dimid,         &
               depolmod_varid))
    call check(nf90_put_att(nc_id, depolmod_varid, 'lon_name',                  &
               'Depolarization Ratio Value'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, depolmod_varid, 'axis','I'))  

   call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  
    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, pres_varid, vprestop))
    call check(nf90_put_var(nc_id, srmod_varid, vsrmod))
    call check(nf90_put_var(nc_id, depolmod_varid, vdepolmod))

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_depolnc
!----------------------------------------------------------------------------!

subroutine create_diagPHAnc(fname,dname,vlon,vlat,vprestop_mid,vprestop_bound,vtempmod,vtime,dim,alt,&
                         nlon,nlat)

    use netcdf
    implicit none

    character(LEN=*)   ::  fname,dname
    character  ::  date*8,time*10,zone*5
    integer :: nlon , nlat, iz 
    integer, parameter :: tempmax = 35, ndims = 6, nv=2, ncat=3
    integer,dimension(8)  ::  value
    integer  ::  dim(ndims),alt,dim3(nv),dim4(nv),dim5(nv)
    integer  ::  lon_varid,lat_varid,alt_varid, tempmod_varid,time_varid,tempmod_varid2,srmod_varid3,srmod_varid4,nc_id,alt_varid2,srmod_varid33,srmod_varid44,cat_varid
    integer  ::  lon_dimid,lat_dimid,alt_dimid,tempmod_dimid,time_dimid,tempmod_dimid2,srmod_dimid3,srmod_dimid4, nv_dimid,sr_dimid,sr_dimid2, cat_dimid
    real  ::  vtime
    real*4, dimension(nlon) :: vlon
    real*4, dimension(nlat) :: vlat
    real, dimension(tempmax-1)  ::  vtempmid
    real, dimension(tempmax)  ::  vtempmod
    real, dimension(tempmax-1,2)  ::  vtempmod2
    real,dimension(alt)  ::  vprestop_mid
    real,dimension(alt,2)  ::  vprestop_bound
    character(len=3),dimension(ncat,ncat) :: vcat

vcat(1,:)='LIQ'
vcat(2,:)='ICE'
vcat(3,:)='VAP'



      do iz=1,tempmax-1
         vtempmod2(iz,1)=vtempmod(iz)
         vtempmod2(iz,2)=vtempmod(iz+1)
         vtempmid(iz)=(vtempmod(iz)+vtempmod(iz+1))/2
      enddo



    call date_and_time(date,time,zone,value)

    call check(nf90_create('./'//fname, NF90_CLOBBER,  &
               nc_id))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Description','GOCCP_Histogram_of_Phase_file'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Date',dname))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Version','Prog_version'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Author',                    &
               'Gregory CESANA, Helene CHEPFER, LMD/IPSL'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Scientific_contact',        &
               'helene.chepfer@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Technical_support',         &
               'gregory.cesana@lmd.polytechnique.fr'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Creationdate',date))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'Website',                &
               'http://climserv.ipsl.polytechnique.fr/cfmip-obs.html'))
    call check(nf90_put_att(nc_id, NF90_GLOBAL, 'References',                &
               'Chepfer et al, 2010, The GCM Oriented CALIPSO Cloud Product (CALIPSO-GOCCP), JGR, 105'))
    call check(nf90_def_dim(nc_id, 'longitude', nlon, lon_dimid))
    call check(nf90_def_dim(nc_id, 'latitude', nlat, lat_dimid))
    call check(nf90_def_dim(nc_id, 'altitude', alt, alt_dimid))
    call check(nf90_def_dim(nc_id, 'tempbox', tempmax-1, tempmod_dimid )) 
    call check(nf90_def_dim(nc_id, 'cat', ncat, cat_dimid))
    call check(nf90_def_dim(nc_id, 'nv', nv, nv_dimid))

    call check(nf90_def_dim(nc_id, 'time', NF90_UNLIMITED, time_dimid))
  
    dim = (/lon_dimid, lat_dimid, alt_dimid, tempmod_dimid, cat_dimid, time_dimid/)
    dim3 = (/alt_dimid, nv_dimid/)
    dim4 = (/tempmod_dimid, nv_dimid/)
    dim5 = (/cat_dimid, cat_dimid /)


    call check(nf90_def_var(nc_id, 'longitude', NF90_FLOAT, lon_dimid, lon_varid))
    call check(nf90_put_att(nc_id, lon_varid, 'lon_name','Longitude'))
    call check(nf90_put_att(nc_id, lon_varid, 'units','degrees_east'))
    call check(nf90_put_att(nc_id, lon_varid, 'axis','X'))

    call check(nf90_def_var(nc_id, 'latitude', NF90_FLOAT, lat_dimid, lat_varid))
    call check(nf90_put_att(nc_id, lat_varid, 'lon_name','Latitude'))
    call check(nf90_put_att(nc_id, lat_varid, 'units','degrees_north'))
    call check(nf90_put_att(nc_id, lat_varid, 'axis','Y'))

    call check(nf90_def_var(nc_id, 'alt_mid', NF90_FLOAT, alt_dimid, alt_varid))
    call check(nf90_put_att(nc_id, alt_varid, 'lon_name','Middle of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid, 'units','kilometer'))
    call check(nf90_put_att(nc_id, alt_varid, 'positive','up'))
    call check(nf90_put_att(nc_id, alt_varid, 'axis','Z'))    

    call check(nf90_def_var(nc_id, 'alt_bound', NF90_FLOAT, dim3, alt_varid2))
    call check(nf90_put_att(nc_id, alt_varid2, 'lon_name','Boundaries of the altitude bin'))
    call check(nf90_put_att(nc_id, alt_varid2, 'units','kilometer'))

    call check(nf90_def_var(nc_id, 'temp_mid', NF90_FLOAT, tempmod_dimid,         &
               tempmod_varid))
    call check(nf90_put_att(nc_id, tempmod_varid, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, tempmod_varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, tempmod_varid, 'axis','I')) 

    call check(nf90_def_var(nc_id, 'temp_bound', NF90_FLOAT, dim4,         &
               tempmod_varid2))
    call check(nf90_put_att(nc_id, tempmod_varid2, 'lon_name',                  &
               'Middle Value of Scattering Ratio Boxes'))
    call check(nf90_put_att(nc_id, tempmod_varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(nc_id, tempmod_varid2, 'axis','I')) 
  
    call check(nf90_def_var(nc_id, 'category', NF90_CHAR, dim5, cat_varid))
    call check(nf90_put_att(nc_id, cat_varid, 'lon_name','Category'))
    call check(nf90_put_att(nc_id, cat_varid, 'units','Arbitrary unit'))


    call check(nf90_def_var(nc_id, 'time', NF90_FLOAT, time_dimid, time_varid))
    call check(nf90_put_att(nc_id, time_varid, 'lon_name','Time'))
    call check(nf90_put_att(nc_id, time_varid, 'units',                      &
               'days since 2000-01-01 00:00:00'))
    call check(nf90_put_att(nc_id, time_varid, 'axis','T'))  
    call check(nf90_put_att(nc_id, time_varid, 'comment',                    &
               'monthly means: date is set to the 15th of the month'))  



    call check(nf90_enddef(nc_id))
	
    call check(nf90_put_var(nc_id, lon_varid, vlon))
    call check(nf90_put_var(nc_id, lat_varid, vlat))
    call check(nf90_put_var(nc_id, alt_varid, vprestop_mid))
    call check(nf90_put_var(nc_id, alt_varid2, vprestop_bound))

    call check(nf90_put_var(nc_id, tempmod_varid, vtempmid))
    call check(nf90_put_var(nc_id, tempmod_varid2, vtempmod2))

    call check(nf90_put_var(nc_id, cat_varid, vcat))

    call check(nf90_put_var(nc_id, time_varid, vtime))
  
    call check(nf90_close(nc_id))

end subroutine create_diagPHAnc
!----------------------------------------------------------------------------!

subroutine diagPHA_recvar2nc3(diagpha,dim,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: tempmax = 34, ndims = 6
    integer  ::  dim(ndims),alt
    integer  ::  varid, varid2, ncid
    character(LEN=*)   ::  fname
    real,dimension(nlon,nlat,alt,tempmax,3)  ::  diagpha


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarPHASE', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar ICE/LIQUID/VAPOR Phase occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))


    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diagpha(:,:,:,:,:)))

    call check(nf90_close(ncid))
    
end subroutine diagPHA_recvar2nc3

!----------------------------------------------------------------------------!
! *** DIAG_RECVAR2NC *** This routine record the diag variables in the netcdf!
!                        prof file                                           !
!----------------------------------------------------------------------------!
! fname      : name of output netcdf diagSR file                             !
! diag       : diagSR variable                                               !
! dim        : dimension id of the isccp variables                           !
! nlon       : number of the longitude boxes                                 !
! nlat       : number of the latitude boxes                                  !
! alt        : number of model altitude                                      !
!----------------------------------------------------------------------------!
! nan        : Not a Number value                                            !
! diagmax    : number of SR diagbox                                          !
! ndims      : dimension of dim (= 5 : lon,lat,alt,diagbox,time)             !
! ncid       : netcdf file id                                                !
! varid      : variable id of diagSR                                         !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : diag_recvar2nc(monthdiagSR,dimidsd,file8,altmax,lonmax,latmax)        !
!                                                                            !
!----------------------------------------------------------------------------!
subroutine diag_recvar2nc(diag,diag1,dim,dim2,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 15, ndims = 5
    integer  ::  dim(ndims), dim2(ndims-1),alt
    integer  ::  varid,varid2, ncid
    character(LEN=*)   ::  fname
    real*4,dimension(nlon,nlat,alt,diagmax)  ::  diag
    real*4,dimension(nlon,nlat,alt)  ::  diag1


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Number of Scattering Ratio occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))

    call check(nf90_def_var(ncid, 'negdiagSR_Occ', NF90_FLOAT, dim2, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                         &
             'Number of Scattering Ratio occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))
    call check(nf90_put_var(ncid, varid2, diag1))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc
!----------------------------------------------------------------------------!



subroutine diag_recvar2nc2(diag,dim,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 16, ndims = 5
    integer  ::  dim(ndims), alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    real*8,dimension(nlon,nlat,alt,diagmax)  ::  diag


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Frac', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
           'Lidar Scattering Ratio CFAD (532nm) Fraction of occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','1 fraction'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc2
!----------------------------------------------------------------------------!

subroutine diag_recvar2nc3(diag,dim,dim2,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 18, ndims = 5
    integer  ::  dim(ndims),dim2(ndims),alt
    integer  ::  varid, varid2, ncid
    character(LEN=*)   ::  fname
    real,dimension(nlon,nlat,alt,diagmax)  ::  diag


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ2', NF90_FLOAT, dim2, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))



    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag(:,:,:,4:18)))
    call check(nf90_put_var(ncid, varid2, diag(:,:,:,1:3)))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc3
!----------------------------------------------------------------------------!

subroutine diag_recvar2nc3pha(diag,dim,dim2,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999.
    integer, parameter :: diagmax = 11, ndims = 5
    integer  ::  dim(ndims),dim2(ndims),alt
    integer  ::  varid, varid2, varid3, ncid
    character(LEN=*)   ::  fname
    real,dimension(nlon,nlat,alt,diagmax,3)  ::  diag


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_liq', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cfad_lidarsr532_ice', NF90_FLOAT, dim, varid2))
    call check(nf90_put_att(ncid, varid2, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid2, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid2, '_FillValue',nan))

   call check(nf90_def_var(ncid, 'cfad_lidarsr532_un', NF90_FLOAT, dim, varid3))
    call check(nf90_put_att(ncid, varid3, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid3, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid3, '_FillValue',nan))



    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag(:,:,:,:,1)))
    call check(nf90_put_var(ncid, varid2, diag(:,:,:,:,2)))
    call check(nf90_put_var(ncid, varid3, diag(:,:,:,:,3)))

    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc3pha
!----------------------------------------------------------------------------!

subroutine depol_recvar2nc3(diag,dim,fname,alt,nlon,nlat)

    use netcdf
    implicit none

    integer :: nlon , nlat 
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 17, ndims = 6, depolmax = 21
    integer  ::  dim(ndims),alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    real*8,dimension(nlon,nlat,alt,diagmax-1,depolmax-1)  ::  diag


    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidardepol532_Occ', NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Depolarization Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))


    call check(nf90_close(ncid))
    
end subroutine depol_recvar2nc3
!----------------------------------------------------------------------------!


subroutine diag_recvar2nc4(diag,dim,fname,alt,nlon,nlat,ndiag)

    use netcdf
    implicit none

    integer :: nlon , nlat
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 18, ndims = 4
    integer  ::  dim(ndims),alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    character  ::  ndiag*2
    integer,dimension(nlon,nlat,alt)  ::  diag

print *, 'monthdiagSR_Occ'// trim(ADJUSTL(ndiag))
    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidarsr532_Occ'//trim(ADJUSTL(ndiag)), NF90_INT4, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Scattering Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))


    call check(nf90_close(ncid))
    
end subroutine diag_recvar2nc4
!----------------------------------------------------------------------------!


subroutine depol_recvar2nc4(diag,dim,fname,alt,nlon,nlat,ndiag)

    use netcdf
    implicit none

    integer :: nlon , nlat
    real,parameter  ::   nan=-9999
    integer, parameter :: diagmax = 16, ndims = 5
    integer  ::  dim(ndims),alt
    integer  ::  varid, ncid
    character(LEN=*)   ::  fname
    character  ::  ndiag*2
    real*8,dimension(nlon,nlat,alt,diagmax)  ::  diag

print *, 'monthdepolSR_Occ'// trim(ADJUSTL(ndiag))
    call check(nf90_open('./'//fname,NF90_WRITE,ncid))

    call check(nf90_redef(ncid))

    call check(nf90_def_var(ncid, 'cfad_lidardepol532_Occ'//trim(ADJUSTL(ndiag)), NF90_FLOAT, dim, varid))
    call check(nf90_put_att(ncid, varid, 'lon_name',                         &
             'Lidar Depolarization Ratio CFAD (532nm) occurence accumulated over a month'))
    call check(nf90_put_att(ncid, varid, 'units','Arbitrary unit'))
    call check(nf90_put_att(ncid, varid, '_FillValue',nan))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, diag))


    call check(nf90_close(ncid))
    
end subroutine depol_recvar2nc4
!----------------------------------------------------------------------------!







subroutine rdnc3(fname,var,alt,nlon,nlat,nvar)
  use netcdf
  implicit none
 character(len=*)  ::  nvar, fname
integer  ::  ncid,varid
integer :: nlon , nlat, alt
integer,dimension(nlon,nlat,alt) :: var

print *, nvar

 call check(NF90_OPEN(fname,NF90_NOWRITE,ncid))
 call check(NF90_inq_varid(ncid,nvar,varid))
 call check(NF90_get_var(ncid,varid,var))
 call check(NF90_CLOSE(ncid))


end subroutine rdnc3

subroutine rdnc4(fname,var,alt,nlon,nlat,nvar)
  use netcdf
  implicit none
 character(len=*)  ::  nvar, fname
integer  ::  ncid,varid
integer, parameter  :: diagmax = 17
integer :: nlon , nlat, alt
real*8,dimension(nlon,nlat,alt,diagmax-1) :: var

print *, nvar

 call check(NF90_OPEN(fname,NF90_NOWRITE,ncid))
 call check(NF90_inq_varid(ncid,nvar,varid))
 call check(NF90_get_var(ncid,varid,var))
 call check(NF90_CLOSE(ncid))


end subroutine rdnc4



!****************************************************************************!
