! subroutines relative to ATB calculations


!----------------------------------------------------------------------------!
! ***** ATB_MOL_INTERP ***** This subroutine extrapolate the atbmol from the !
!                            last atbmol known value to the Surface Elevation!
!                            in the limit of the 5 first kilometers          !
!----------------------------------------------------------------------------!
! i        : loop index on the profil                                        !
! nprofs   : number of profil                                                !
! alt      : altitude of lidar lvl in kilometer = altl                       ! 
!----------------------------------------------------------------------------!
! var3     : molecular interpolated in km-1 sr-1                             !
!----------------------------------------------------------------------------!
! l        : loop index on alt3 altitude levels                              !
! alt3     : selection of value of alt the nearest to 0 1 2 3 4 & 5km        !
! a,b      : coefficient for the interpolation equation                      !
! ilid     : loop index of the lidar altitude (583)                          !
! SeuilMol1km : treshold to add to the molecular to get the value of         !
!               molecular 1kilometer lower (this treshold is available under !
!               5km altitude)                                                !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call atb_mol_interp(mol3,altl,i,it)                                  !
!                                                                            !
!----------------------------------------------------------------------------!
!  Purpose : When the molecular value is missing close to the ground (e.g.   !
!            under 5km), this routine extrapolate this value with the tresh- !
!            old "SeuilMol1km" applied to the last known value of molecular  !
!            at 5 4 3 2 or 1 km. Then we have 2 values of molecular distant  !
!            from 1km. At this moment the routine perform the interpolation  !
!            between this 2 values at the other altitude. This process       !
!            don't stop until reaching the surface elevation.                !
!----------------------------------------------------------------------------!
subroutine atb_mol_interp(var3,alt,i,nprofs,seuil,SE)
  
  implicit none
  integer  ::  ilid, l, i, n, k

  real*4  :: seuil 
  real*8  ::  a,b
  integer  ::  nprofs  
  integer,dimension(10)  ::  alt3
  real,dimension(583,nprofs)  ::  var3    
  real,dimension(583) :: alt,SE



! altitude of 5 level from 0 to 5km each km (made from altl)
alt3(1)=275; alt3(2)=296; alt3(3)=329; alt3(4)=362; alt3(5)=395;  
alt3(6)=429; alt3(7)=462; alt3(8)=495; alt3(9)=529; alt3(10)=562;


if( sum(var3(1:275,i)).lt.0. )then
! Exclude entiere NaN profile
var3(:,i)=-777.

else

do ilid=275,563 

   ! Exclude the nan value
   if( (var3(ilid,i).eq.(-9999.)) .or. (var3(ilid,i).eq.(-777.)))then
     do l=1,9

      if( (ilid.gt.alt3(l)).and.(ilid.le.alt3(l+1)) ) then
           ! calculation of extrapolated value 
           var3(alt3(l+1),i)=var3(alt3(l)-1,i)+seuil
      ! Calculation of the parameter to interpolate the value between the
      ! extrapolated value and the last known value        
      a=(alt(alt3(l+1))-alt(alt3(l)-1))/(var3(alt3(l+1),i)-var3(alt3(l)-1,i))
      b=alt(alt3(l+1))-a*var3(alt3(l+1),i)
         
         ! Calculation of the new interpolated variable  
         do k=alt3(l),alt3(l+1)
            var3(k,i)=(alt(k)-b)/a
         enddo
      endif

     enddo
   endif

enddo

endif
end subroutine atb_mol_interp
!----------------------------------------------------------------------------!

subroutine atb_temp_interp(var3,alt,i,nprofs,seuil,SE)
  
  implicit none
  integer  ::  ilid, l, i, n, k

  real*4  :: seuil 
  real*8  ::  a,b
  integer  ::  nprofs  
  integer,dimension(10)  ::  alt3
  real,dimension(583,nprofs)  ::  var3    
  real,dimension(583) :: alt,SE



! altitude of 5 level from 0 to 5km each km (made from altl)
alt3(1)=275; alt3(2)=296; alt3(3)=329; alt3(4)=362; alt3(5)=395;  
alt3(6)=429; alt3(7)=462; alt3(8)=495; alt3(9)=529; alt3(10)=562;


if ( ((var3(150,i).eq.(-9999.)).and.(var3(200,i).eq.(-9999.)).and.       &
     (var3(250,i).eq.(-9999.)).and.(var3(329,i).eq.(-9999.)).and.        &
     (var3(429,i).eq.(-9999.)).and.(var3(529,i).eq.(-9999.))) .or.       &
     ((var3(150,i).eq.(-777.)).and.(var3(200,i).eq.(-777.)).and.         &
     (var3(250,i).eq.(-777.)).and.(var3(329,i).eq.(-777.)).and.          &
     (var3(429,i).eq.(-777.)).and.(var3(529,i).eq.(-777.))) )then


continue
else

do ilid=275,563 

   ! Exclude the nan value
   if( (var3(ilid,i).eq.(-9999.)) .or. (var3(ilid,i).eq.(-777.)))then
     do l=1,9

      if( (ilid.gt.alt3(l)).and.(ilid.le.alt3(l+1)) ) then
           ! calculation of extrapolated value 
           var3(alt3(l+1),i)=var3(alt3(l)-1,i)+seuil
      ! Calculation of the parameter to interpolate the value between the
      ! extrapolated value and the last known value        
      a=(alt(alt3(l+1))-alt(alt3(l)-1))/(var3(alt3(l+1),i)-var3(alt3(l)-1,i))
      b=alt(alt3(l+1))-a*var3(alt3(l+1),i)
         
         ! Calculation of the new interpolated variable  
         do k=alt3(l),alt3(l+1)
            var3(k,i)=(alt(k)-b)/a
         enddo
      endif

     enddo
   endif

enddo

endif
end subroutine atb_temp_interp
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
! ***** ATB_MOL ***** This routine calculate the atbmol from 33lvl to 583lvl !
!                     including the calculation of the rapport               !
!----------------------------------------------------------------------------!
! alt1     : top calculation altitude level of the normalized ratio          !
! alt2     : toplow calculation altitude level of the normalized ratio       !
! nprofs   : number of profil                                                !
! i        : loop index on the profil                                        !
! var      : attenuated backscatter (=atb(583,nprof))                        !
! var2     : molecular from 33lvl to 583lvl in CN (=mol2(583,nprof))         !
!----------------------------------------------------------------------------!
! var3     : molecular in km-1 sr-1 from mol2 (=mol3(583,nprof))             !
!----------------------------------------------------------------------------!
! n        : loop index of the sliding average                               !
! ilid     : loop index of the lidar altitude (583)                          !
! matb     : sum of the atb                                                  !
! mmol     : sum of the molecular                                            !
! rapport3 : normalized ratio for each profil dim=(nprof)                    !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call atb_mol(atb,mol2,mol3,i,it,62,92)                               !
!                                                                            !
!----------------------------------------------------------------------------!
!  Purpose : In order to caculate the scattering ratio, we need to estimate  !
!            the molecular in sr-1 km-1 because :                            ! 
!            SR = atb(km-1 sr-1) / mol(km-1 sr-1)                            !
!            and we already know the atb value.                              !
!            The calculation is a sliding average of 66 profils on 20km      !
!            (1profil every 330m, 66*330m = 20km). This operation allows to  !
!            refine the molecular.                                           !
!----------------------------------------------------------------------------!
subroutine atb_mol(var,var2,var3,i,nprofs,alt1,alt2)

  implicit none
  integer  ::  ilid, i , n ,alt1, alt2, it     
  real  ::  matb, mmol 
  integer  ::  nprofs                    
  real,dimension(583,nprofs)  ::  var, var2, var3    
  real,dimension(nprofs)  ::  rapport3
  
  matb=0
  mmol=0
  rapport3(i)=0

  ! change vnoel 20140224 - "it" was not defined in the subroutine. But it 
  ! appears in the example comments above in the place of "nprofs" arguments.
  it = nprofs

     ! Average for the first profils
     if(i.lt.34)then  !!!!! loop on profil
        matb=0
        mmol=0 

        do n=0,65
           ! average between 22-25km or 20-25km
           do ilid=alt1,alt2   
               ! exclude the nan values
               if ((var(ilid,i).ne.(-9999)).and.(var(ilid,i).lt.1).and.      &
                   (var2(ilid,i).ne.-9999))then
                matb=matb+var(ilid,1+n)
                mmol=mmol+var2(ilid,1+n) 
              ! else
               ! rcompt=rcompt+1    !! check this process & the value rcompt
               endif
            enddo
         enddo
         ! calculation of the ratio
         rapport3(i)=matb/mmol

    ! Average for other profils
    else 
       if( (i.gt.32).and.(i.lt.(it-32)) )then
           matb=0
           mmol=0 

          do n=0,65
             ! average between 22-25km or 20-25km
             do ilid=alt1,alt2  
                ! exclude the nan values
                if ((var(ilid,i).ne.(-9999)).and.(var(ilid,i).lt.1).and.     &
                    (var2(ilid,i).ne.-9999))then
                   matb=matb+var(ilid,i-32+n)
                   mmol=mmol+var2(ilid,i-32+n) 
            !    else
              !     rcompt=rcompt+1   !! check this process & the value rcompt
                endif
             enddo
          enddo
          rapport3(i)=matb/mmol
          
    ! Average for the last profils      
          else
             matb=0
             mmol=0 
          
           do n=0,65
              ! average between 22-25km or 20-25km
              do ilid=alt1,alt2 
                 ! exclude the nan values
                 if ((var(ilid,i).ne.(-9999)).and.(var(ilid,i).lt.1).and.    &
                    (var2(ilid,i).ne.-9999))then
                   matb=matb+var(ilid,it-65+n);
                   mmol=mmol+var2(ilid,it-65+n) ;
                 endif
              enddo
           enddo
          rapport3(i)=matb/mmol

       endif

    endif !!!!! END loop on profil

! Calculation of the new molecular in km-1 sr-1   
do ilid=1,583
   if((var2(ilid,i).ne.(-9999)).and.(rapport3(i).ge.0.25e-28).and.(rapport3(i).le.0.97e-28))then
      var3(ilid,i)=rapport3(i)*var2(ilid,i)
   elseif(var2(ilid,i).eq.(-9999.))then
      var3(ilid,i)=-9999. ! allocation of nan values 
   else
    
     var3(ilid,i)=-777.
   endif

enddo
 
end subroutine atb_mol
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! ***** INTERP ***** This routine proceed a linear interpolation of variables!
!                    entered                                                 !
!----------------------------------------------------------------------------!
! i        : loop index on the profil                                        !
! var      : non interpolated variable dim=33,nprofs                         !
! fvar     : altitude of meteo lvl in kilometer, = altm                      !
! fvar2    : altitude of lidar lvl in kilometer = altl                       !
! nprofs   : number of profil                                                !
!----------------------------------------------------------------------------!
! var2     : interpolated variable dim=583,nprofs                            !
!----------------------------------------------------------------------------!
! a,b      : coefficient for the interpolation equation                      !
! imol     : loop index of the meteo altitude (33)                           !
! ilid     : loop index of the lidar altitude (583)                          !
! altitude : number of level for the lidar variables                         !
! altitude2: number of level for the meteo variables                         !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call interp(pres,pres2,altm,altl,i,it)                                !
!----------------------------------------------------------------------------!
subroutine interp(var,var2,fvar,fvar2,i,nprofs)
  implicit none

!     Indexes & parameters 
      real  ::  a,b                   
      integer  ::  i
      integer  ::  ilid                       
      integer  ::  imol                      
      integer(kind=2),parameter :: altitude2 = 33 
      integer(kind=2),parameter :: altitude = 583  
      integer  ::nprofs 
    
!     META variables
      real*4,dimension(583)  ::  fvar2  
      real*4,dimension(33)  ::   fvar 
      real,dimension(33,nprofs)  ::  var 
      real,dimension(583,nprofs)  ::  var2


    do imol=2,altitude2 
      ! exclude the nan value
      if( ( (var(imol,i).lt.1E+16).and.(var(imol,i).gt.1E+13).or.            &
         (var(imol,i).eq.-9999.) ).or.( (var(imol-1,i).lt.1E+16).and.         &
         (var(imol-1,i).gt.1E+13).or.(var(imol-1,i).eq.-9999.) ) )then

        do ilid=1,altitude
        if ((fvar2(ilid).ge.fvar(imol)).and.(fvar2(ilid).lt.fvar(imol-1)))then
            var2(ilid,i)=-9999.
        endif
        enddo

      else
      ! calculation of coefficient
         a=(var(imol,i)-var(imol-1,i))/(fvar(imol)-fvar(imol-1))
         b=var(imol,i)-a*fvar(imol)
        
        ! calculation of new interpolated variable
        do ilid=1,altitude
        if ((fvar2(ilid).ge.fvar(imol)).and.(fvar2(ilid).lt.fvar(imol-1)))then
           var2(ilid,i)=a*fvar2(ilid)+b
        endif
        enddo

      endif
    enddo

end subroutine interp
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! **** ZERO_DETECT **** This routine detect the 0 values and substitute them !
!                        by -777                                            !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! iz     : loop index on the altitude grid                                   !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
!----------------------------------------------------------------------------!
! var    : verticaly averaged observed variables                             !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call zero_detect(pr2moy,i,iz,it,altmax)                              !
!----------------------------------------------------------------------------!
subroutine zero_detect(var,i,iz,nprofs,alt)

  implicit none
  integer ::  i,iz,nprofs,alt
  real*4,dimension(alt,nprofs)  ::  var
  
 if(var(iz,i).eq.0)then       ! if no value is detected then allocate -9999
            var(iz,i)=-9999.
 endif

end subroutine zero_detect
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! **** ZERO_DETECT_CHIM **** This routine detect the 0 values and substitute !
!                              them by -9999                                 !
!----------------------------------------------------------------------------!
! j      : loop index on latitude grid                                       !
! k      : loop index on longitude grid                                      ! 
! lat    : number of the latitude boxes                                      !
! lon    : number of the longitude boxes                                     !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
!----------------------------------------------------------------------------!
! var    : verticaly averaged observed variables                             !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call zero_detect_chim(pr2moy,iz,altmax,j,k,latmax,lonmax)            !
!----------------------------------------------------------------------------!
subroutine zero_detect_chim(var,iz,alt,j,k,lat,lon)

  implicit none
  integer ::  iz,alt,j,k,lat,lon
  real*4,dimension(lat,lon,alt)  ::  var
  
 if(var(j,k,iz).eq.0)then       ! if no value is detected then allocate -9999
            var(j,k,iz)=-9999
 endif

end subroutine zero_detect_chim
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** SR_CR_DEPOL_chim *** This routine do the SR CR & DEPOL average except  !
!                          for the -9999 values                              !
!----------------------------------------------------------------------------!
! j      : loop index on latitude grid                                       !
! k      : loop index on longitude grid                                      ! 
! lat    : number of the latitude boxes                                      !
! lon    : number of the longitude boxes                                     !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
! num    : index to select the method of calculation                         !
! var1, var2 : verticaly averaged observed variables                         !
! var4   : scattering ratio vertically averaged                              ! 
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            !
!----------------------------------------------------------------------------!
! var3   : result of var1/var2                                               !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call SR_CR_DEPOL_chim(pr2moy,molmoy,srmoy,indice,indicem,iz,altmax,  !
!                             j,k,latmax,lonmax,1,srmoy)                     !
!----------------------------------------------------------------------------!
subroutine SR_CR_DEPOL_chim(var1,var2,var3,ind1,ind2,iz,alt,j,k,lat,lon,num, &
                            var4)

  implicit none
  integer  ::  iz,alt,j,k,lat,lon,num
  real*4,dimension(lat,lon,alt)  ::  var1,var2,var3,var4,ind1,ind2

if(num.eq.1)then

   ! calculation if var3 = srmoy
   if ((var1(j,k,iz).eq.(-9999)).or.(var2(j,k,iz).EQ.(-9999))) then
      var3(j,k,iz) = -9999
   else
      var3(j,k,iz) = (var1(j,k,iz)*ind2(j,k,iz))/(var2(j,k,iz)*ind1(j,k,iz))
   endif

else 

    ! calculation if var3 is different than srmoy  
   if(var4(j,k,iz).eq.(-9999))then
      var3(j,k,iz) = -9999
   else
      if ((var1(j,k,iz).eq.(-9999)).or.(var2(j,k,iz).EQ.(-9999))) then    
       var3(j,k,iz) = -9999
      else
       var3(j,k,iz) = (var1(j,k,iz)*ind2(j,k,iz))/(var2(j,k,iz)*ind1(j,k,iz))
      endif
   endif

endif

end subroutine SR_CR_DEPOL_chim
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** SR_CR_DEPOL_mean *** This routine do the SR CR & DEPOL average except  !
!                          for the -9999 values                              !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! iz     : loop index on the altitude grid                                   !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1, var2 : verticaly averaged observed variables                         !
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            !
!----------------------------------------------------------------------------!
! var3   : result of var1/var2                                               ! 
!----------------------------------------------------------------------------!
!                                                                            !
! ex :   call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   !
!                              altmax)                                       !
!----------------------------------------------------------------------------!
subroutine SR_CR_DEPOL_mean(var1,var2,var3,ind1,ind2,i,iz,nprofs,alt)
 
  implicit none
  integer  ::  i,iz,nprofs,alt
  real,dimension(alt,nprofs)  ::  var1,var2,var3
  real,dimension(alt,nprofs)  ::  ind1,ind2

    if ((var1(iz,i).eq.(-9999.)).or.(var2(iz,i).EQ.(-9999.))) then    
     var3(iz,i) = -9999.
    elseif ((var1(iz,i).eq.(-888.)).or.(var2(iz,i).EQ.(-888.))) then 
     var3(iz,i) = -888.
    elseif ((var1(iz,i).eq.(-777.)).or.(var2(iz,i).EQ.(-777.))) then 
     var3(iz,i) = -777.  
    else
     var3(iz,i) = (var1(iz,i)*ind2(iz,i)) / (var2(iz,i)*ind1(iz,i))
    endif

end subroutine SR_CR_DEPOL_mean
!----------------------------------------------------------------------------!



!----------------------------------------------------------------------------!
! *** FILTRE_2LVL *** This routine delete the cloud over 21km and do the     !
!                     the SR average execpt for the nan values.              !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! iz     : loop index on the altitude grid                                   !
! nprofs : number of profil                                                  !
! alt    : number of the altitude boxes                                      !
! var1, var2 : verticaly averaged observed variables                         !
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            !
!----------------------------------------------------------------------------!
! var3   : result of var1/var2                                               ! 
!----------------------------------------------------------------------------!
!                                                                            !
! ex :   call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   !
!                              altmax)                                       !
!----------------------------------------------------------------------------!
!  Purpose : noise could be interprated as cloudy point in the program so,   !
!            we have built this routine to filter this phenomena over 21km   !
!            because cloud couldn't appear in this area execpt during the    !
!            polar winter.                                                   !
!----------------------------------------------------------------------------!
! change 20140224 : "gcm" argument was named "grid"
subroutine filtre_2lvl(var1,var2,var3,ind1,ind2,i,iz,nprofs,alt,gcm)

 implicit none
 integer  ::  i,iz,nprofs,alt,lvl
 real*4,dimension(alt,nprofs)  ::  var1,var2,var3
 real*4,dimension(alt,nprofs)  ::  ind1,ind2
 character  ::  gcm*8

if(gcm.eq.'LMDZ')then
   lvl=15   ! level to go past in order to apply the filter = 20.8km
elseif(gcm.eq.'LMDZ40')then
   lvl=30   !  level to go past in order to apply the filter = 21.1km
elseif(gcm.eq.'WRF')then
   lvl=46   !  level to go past in order to apply the filter = 21.7km
endif



  if (iz.ge.lvl) then

  ! filter the cloud over 21km + nan values
    if((((var1(iz,i)*ind2(iz,i))/(var2(iz,i)*ind1(iz,i))).ge.5.) .or.&
           (var1(iz,i).eq.(-9999)).or.(var2(iz,i).EQ.(-9999)).or.&
           ((var1(iz,i).eq.(-777)).or.(var2(iz,i).EQ.(-777))))then
       var3(iz,i) = -777       
    else
       ! calculation of the average
       var3(iz,i) = (var1(iz,i)*ind2(iz,i))/(var2(iz,i)*ind1(iz,i))
    endif

  elseif((var1(iz,i).eq.(-9999)).or.(var2(iz,i).EQ.(-9999))) then 
! filter the nan values       
     var3(iz,i) = -9999
  elseif ((var1(iz,i).eq.(-888)).or.(var2(iz,i).EQ.(-888))) then 
     var3(iz,i) = -888
  elseif ((var1(iz,i).eq.(-777)).or.(var2(iz,i).EQ.(-777))) then 
     var3(iz,i) = -777       
  else
        ! calculation of the average
     var3(iz,i) = (var1(iz,i)*ind2(iz,i))/(var2(iz,i)*ind1(iz,i))
  endif

   
end subroutine filtre_2lvl
!----------------------------------------------------------------------------!
