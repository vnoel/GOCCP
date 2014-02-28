! surface.f90

! subroutines dealing with surface elevation

!----------------------------------------------------------------------------!
! *** SE_KM_2_PRES *** This subroutine convert the SE from km to hPa and     !
!                      allocate -888 to the variable below this threshold    !
!----------------------------------------------------------------------------!
! nprofs : number of profil                                                  !
! i      : loop index on the profil                                          !
! alt    : number of the altitude boxes                                      !
! var1   : Surface_elevation (nprofs)                                        !
! var3   : altitude of lidar lvl in kilometer, (583)                         !
! var4   : Pressure interpolated from 33lvl to 583lvl in hPa (583,nprofs)    !
! var5   : pressure of model boxes                                           !
!----------------------------------------------------------------------------!
! var2   : Surface_elevation in hPa  (nprofs)                                !
! var6   : add the surface_elevation to this variable                        !
!----------------------------------------------------------------------------!
! iz     : loop index on the altitude grid                                   !
! j      : loop index on Surface Elevation level                             !
! ilid   : loop index of the lidar altitude (583)                            !
! alt2   : number of level of the lidar variables                            !
! a,b    : coefficient for the interpolation equation                        !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call SE_km_2_Pres(SE,SEp,altl,pres2,prestop,srmoy,altmax,it,i)        !
!----------------------------------------------------------------------------!
subroutine SE_km_2_Pres(var1,var2,var3,var4,var5,var6,alt,nprofs,i)

implicit none
integer  ::  ilid, iz, alt, nprofs,i,j
integer,parameter  ::  altitude=583
real*4  ::  a,b
real*4,dimension(altitude,nprofs)  ::  var4
real*4,dimension(alt)  ::  var5
real,dimension(nprofs)  ::  var1,var2
real*4,dimension(alt,nprofs)  ::  var6
real*4,dimension(altitude)  ::  var3
! var1 = SE, var2 = SEp, var3 = altl, var4 = pres2, var5 = prestop, var6 = srmoy
do ilid=1,altitude
  ! print *, var1(1,i), var3(ilid)
      if( var1(i).eq.var3(ilid) )then
         var2(i)=var4(ilid,i)
      else
        if( var1(i).eq.var3(ilid+1) )then
             var2(i)=var4(ilid+1,i)
        else
          if((var1(i).lt.var3(ilid)).and.(var1(i).gt.var3(ilid+1)))then
             if( var4(ilid,i).ne.-9999 )then
               a=(var4(ilid,i)-var4(ilid+1,i))/(var3(ilid)-var3(ilid+1))
               b=var4(ilid,i)-a*var3(ilid)
               var2(i)=a*var1(i)+b
             else
               var2(i)=-9999
             endif
          endif
        endif
      endif
   enddo

   do iz=alt,2,-1 
         if ( (var2(i).gt.var5(iz)).and.(var2(i).lt.var5(iz-1)) )then
            do j=1,iz
               var6(j,i)=-888
            enddo
         endif
   enddo

end subroutine SE_km_2_Pres
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SE_KM_2_PRES *** This subroutine convert the SE from km to hPa and     !
!                      allocate -9999 to the variable below this threshold   !
!----------------------------------------------------------------------------!
! nprofs : number of profil                                                  !
! i      : loop index on the profil                                          !
! alt    : number of the altitude boxes                                      !
! var1   : Surface_elevation (nprofs)                                        !
! var3   : altitude of lidar lvl in kilometer, (583)                         !
! var4   : Pressure interpolated from 33lvl to 583lvl in hPa (583,nprofs)    !
! var5   : pressure of model boxes                                           !
!----------------------------------------------------------------------------!
! var2   : Surface_elevation in hPa  (nprofs)                                !
! var6   : add the surface_elevation to this variable                        !
!----------------------------------------------------------------------------!
! iz     : loop index on the altitude grid                                   !
! j      : loop index on Surface Elevation level                             !
! ilid   : loop index of the lidar altitude (583)                            !
! alt2   : number of level of the lidar variables                            !
! a,b    : coefficient for the interpolation equation                        !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call SE_km_2_Pres(SE,SEp,altl,pres2,prestop,srmoy,altmax,it,i)        !
!----------------------------------------------------------------------------!
subroutine SE_km_2_Pres2(var1,var2,var3,var4,var5,var6,alt,nprofs,i)

implicit none
integer  ::  ilid, iz, alt, nprofs,i,j
integer,parameter  ::  altitude=583
real*4  ::  a,b
real*4,dimension(altitude,nprofs)  ::  var4
real*4,dimension(alt)  ::  var5
real,dimension(nprofs)  ::  var1,var2
real*4,dimension(alt,nprofs)  ::  var6
real*4,dimension(altitude)  ::  var3
! var1 = SE, var2 = SEp, var3 = altl, var4 = pres2, var5 = prestop, var6 = srmoy
do ilid=1,altitude
  ! print *, var1(1,i), var3(ilid)
      if( var1(i).eq.var3(ilid) )then
         var2(i)=var4(ilid,i)
      else
        if( var1(i).eq.var3(ilid+1) )then
             var2(i)=var4(ilid+1,i)
        else
          if((var1(i).lt.var3(ilid)).and.(var1(i).gt.var3(ilid+1)))then
             if( var4(ilid,i).ne.-9999 )then
               a=(var4(ilid,i)-var4(ilid+1,i))/(var3(ilid)-var3(ilid+1))
               b=var4(ilid,i)-a*var3(ilid)
               var2(i)=a*var1(i)+b
             else
               var2(i)=-9999
             endif
          endif
        endif
      endif
   enddo

   do iz=alt,2,-1 
         if ( (var2(i).gt.var5(iz)).and.(var2(i).lt.var5(iz-1)) )then
            do j=1,iz
               var6(j,i)=-9999
            enddo
         endif
   enddo

end subroutine SE_km_2_Pres2
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SURF_DETECT *** This subroutine detect the surface elevation and       !
!                     allocate -888 to the variable below this threshold     !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! j      : loop index on Surface Elevation level                             !
! iz     : loop index on the altitude grid                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call Surf_detect(SE,altmod,pr2moy,altmax,it,i)                       !
!----------------------------------------------------------------------------!
subroutine Surf_detect(var1,var2,var3,alt,nprofs,i)

  implicit none

  integer  ::  iz, alt, nprofs, i,j
  real*4,dimension(alt+1)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3

  B1: do iz=1,alt-1 
         if ( (var1(i).gt.var2(iz)).and.(var1(i).le.var2(iz+1)) )then
          B2:   do j=1,iz
                   var3(j,i)=-888
                enddo  B2
         exit B1            
         endif 
  
      enddo B1              
   
end subroutine Surf_detect
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SURF_DETECT2 *** This subroutine detect the surface elevation and      !
!                      allocate -9999 to the variable below this threshold   !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! j      : loop index on Surface Elevation level                             !
! iz     : loop index on the altitude grid                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call Surf_detect2(SE,altmod,pr2moy,altmax,it,i)                      !
!----------------------------------------------------------------------------!
subroutine Surf_detect2(var1,var2,var3,alt,nprofs,i,mod)

  implicit none

  integer  ::  iz, alt, nprofs, i,j
  character  ::  mod*8
  real*4,dimension(alt+1)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3

if(mod=='altitude')then
 B11:  do iz=1,alt-1 
         if ( (var1(i).gt.var2(iz)).and.(var1(i).le.var2(iz+1)) )then
            B22: do j=1,iz
                   var3(j,i)=-888
                enddo B22
                exit B11
         endif
      enddo B11             
elseif(mod=='pressure')then
  B12:  do iz=1,alt-1 
         if ( (var1(i).lt.var2(iz)).and.(var1(i).gt.var2(iz+1)) )then
         B21:    do j=1,iz
                   var3(j,i)=-888
                enddo B21
                exit B12
         endif  
      enddo  B12 
endif
  
end subroutine Surf_detect2
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** SE_ALT_CHIM *** This subroutine detect the surface elevation and       !
!                     allocate -9999 to the variable below this threshold for!
!                     the CHIMERE grid                                       !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! l      : loop index on Surface Elevation level                             !
! ilid   : loop index on the lidar altitude (583)                            !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call SE_alt_atb(SE,altl,mol3,altitude,it,i)                         !
!----------------------------------------------------------------------------!
subroutine SE_alt_atb(var1,var2,var3,alt,nprofs,i)

  implicit none

  integer  ::  ilid, i, nprofs,ilid2,k,iSEdown,iSEup
  integer(kind=2)  ::  alt
  real  ::  echothres
  
  real*4,dimension(alt)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3
              
iSEdown=0

! Control of SE value
! if SE=NaN ==> Forced SE to 0
if(var1(i).eq.-9999.)then
var1(i)=0
endif

! Put values below the Surface Elevation threshold to -9999.
! Find the lidar level called "iSEdown" matching with the Surface elevation
   b2 :do ilid=2,alt
         if ( (var1(i).ge.var2(ilid)) )then 
             iSEdown=ilid
              b1:  do ilid2=583,ilid,-1
                 var3(ilid2,i)=-9999.
                   enddo b1
                   exit b2
         endif
   enddo b2

echothres=1.

if(var1(i).gt.0.)then
  b4 :do ilid=2,iSEdown-33
         if(var3(ilid,i).gt.0.15)then
         echothres=0.4
         exit b4
         endif
      enddo b4
endif

iSEup=-9999

! Find the highest lidar ATB value of the lidar ground echo called "iSEup"
if(any(var3(iSEdown-30:iSEdown,i).ge.echothres))then
  b3: do ilid=iSEdown-30,iSEdown
     iSEup=ilid
        if(var3(ilid,i).eq.maxval(var3(iSEdown-30:iSEdown,i)))exit b3
      enddo  b3        
endif

if(iSEup.ne.-9999)then
! Put values from the ground to 3 levels above the iSEup level to -9999. 
   do ilid=iSEup-3,iSEdown
      var3(ilid,i)=-9999.
   enddo
endif



end subroutine SE_alt_atb
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
! *** SE_ALT_CHIM *** This subroutine detect the surface elevation and       !
!                     allocate -9999 to the variable below this threshold for!
!                     the CHIMERE grid                                       !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1   : Surface Elevation                                                 !
! var2   : altitude of model boxes                                           !
!----------------------------------------------------------------------------!
! var3   : varically averaged variable                                       !
!----------------------------------------------------------------------------!
! l      : loop index on Surface Elevation level                             !
! ilid   : loop index on the lidar altitude (583)                            !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call SE_alt_atb(SE,altl,mol3,altitude,it,i)                         !
!----------------------------------------------------------------------------!
subroutine SE_alt_mol(var1,var2,var3,alt,nprofs,i)

  implicit none

  integer  ::  ilid, i, nprofs,l
  integer(kind=2)  ::  alt
  real*4,dimension(alt)  ::  var2
  real,dimension(nprofs)  ::  var1
  real*4,dimension(alt,nprofs)  ::  var3
              
   b2 :do ilid=2,alt
         if ( (var1(i).gt.var2(ilid)) )then 
              b1:  do l=583,(ilid+1),-1
                 var3(l,i)=-9999
                   enddo b1
               exit b2
         endif 
   enddo b2

end subroutine SE_alt_mol
!----------------------------------------------------------------------------!
