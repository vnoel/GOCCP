! subgrid.f90

! subroutines relative to subgrid processing


!----------------------------------------------------------------------------!
! *** FRACTION_SUBGRID2 *** This routine calculate the cloudy,clear,         !
!                           saturated... fraction flag(0/1)                  !
!----------------------------------------------------------------------------!
! nprofs   : number of profil                                                !
! alt      : number of the altitude boxes                                    !
! i        : loop index on the profil                                        !
! switch       : switch to chose cloudy or sat mode                          !
! var          : vertically averaged SR                                      !
! var1, var2   : verticaly averaged observed variables (atb & mol moy)       !
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            ! 
!----------------------------------------------------------------------------!
! frac1  : fully attenuated fraction                                         !
! frac2  : cloudy fraction                                                   !
! frac3  : clear fraction                                                    !
! frac4  : uncertain fraction                                                !
! frac5  : NaN fraction                                                      !
! frac6  : Surface_elevation fraction                                        !
!----------------------------------------------------------------------------!
! iz       : loop index on the altitude grid                                 !
! fracttot : sum of all the fraction (=1)                                    !
! toplvlcloud  : first cloudy point encountered from the top to the ground   !
! toplvlsat    : first fully attenuated point encountered from the top to    !
!                the ground                                                  ! 
! SeuilClearSr : clear treshold detection = 1.2                              !
! SeuilSatSr   : fully attenuated treshold detection = 0.01                  !
! SeuilSrCloud : cloudy treshold detection = 5                               !
! SeuilDeltAtb : delta atb treshold = 1.4e-03 km-1 sr-1                      !
! delta        : atb - mol                                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call fraction_subgrid2(srmoy,pr2moy,indice,molmoy,indicem,satfraction,! 
!                             cloudfraction,clearfraction,uncertfraction,    ! 
!                             nanfraction,sefraction,i,altmax,it,switch2)    !
!----------------------------------------------------------------------------!
 
subroutine fraction_subgrid2(var,var1,ind1,var2,ind2,frac1,frac2,frac3,frac4,&
                             frac5,frac6,frac7,i,alt,nprofs,switch)!frac7,&
                          !   frac8,i,alt,nprofs)!frac1bis,frac2bis

  implicit none
  integer  ::  i,iz,nprofs,alt
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.    
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac1,frac2,frac3,frac4,frac5,frac6,frac7!,frac8!,frac2bis,frac1bis


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2



  do iz=1,alt
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! pts nuageux
               endif
elseif(trim(switch).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).gt.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
             (delta.gt.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
  enddo   



 do  iz=1,alt  
 fracttot=frac1(iz,i)+frac2(iz,i)+frac3(iz,i)+frac4(iz,i)+frac5(iz,i)+       &
          frac6(iz,i)+frac7(iz,i) !check = 1
  if (fracttot.ne.1) then
     print *, "fraction error"
     print *, 'instant',fracttot,frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i)
  endif
enddo
end subroutine fraction_subgrid2
!----------------------------------------------------------------------------!

subroutine fraction_subgrid2_8km(seuilsnrlow,seuilsnrhigh,var,var1,ind1,var2,ind2,frac1,frac2,frac3,frac4,&
                             frac5,frac6,frac7,i,alt,nprofs,toplow,topmid,switch,switch2)!frac7,&
                          !   frac8,i,alt,nprofs)!frac1bis,frac2bis

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid,seuilsnrhigh,seuilsnrlow
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac1,frac2,frac3,frac4,frac5,frac6,frac7!,frac8!,frac2bis,frac1bis


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'


do iz=seuilsnrhigh,alt

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif

          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif 

          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif          
enddo   

do iz=1,seuilsnrlow
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif

           if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=seuilsnrlow+1,seuilsnrhigh-1
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif

           if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac3(iz,i)=frac3(iz,i)+1              
          endif

          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -777 (surface elevation) 
          endif

          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


endif

!print *, frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
!                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i)

 do  iz=1,alt  
 fracttot=frac1(iz,i)+frac2(iz,i)+frac3(iz,i)+frac4(iz,i)+frac5(iz,i)+       &
          frac6(iz,i)+frac7(iz,i) !check = 1
  if (fracttot.ne.1) then
     print *, "fraction error"
     print *, 'instant',fracttot,frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i), i,iz,var(iz,i)
!     stop
  endif
enddo


end subroutine fraction_subgrid2_8km
!----------------------------------------------------------------------------!


subroutine fraction_subgrid2_8km_delta(var,var1,ind1,var2,ind2,frac1,frac2,frac3,frac4,&
                             frac5,frac6,frac7,i,alt,nprofs,toplow,topmid,switch,switch2)!frac7,&
                          !   frac8,i,alt,nprofs)!frac1bis,frac2bis

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac1,frac2,frac3,frac4,frac5,frac6,frac7!,frac8!,frac2bis,frac1bis


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'


do iz=18,alt

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3))then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then!.or. &
         !   ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
         !   (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
         !   (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then!.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif          
enddo   

do iz=1,5
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then!.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=6,17
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac5(iz,i)=frac5(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac6(iz,i)=frac6(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac7(iz,i)=frac7(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac2(iz,i)=frac2(iz,i)+1   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac2(iz,i)=frac2(iz,i)+1    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac1(iz,i)=frac1(iz,i)+1  !
               elseif(iz.gt.toplvlsat)then
                  frac3(iz,i)=frac3(iz,i)+1    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac3(iz,i)=frac3(iz,i)+1   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac4(iz,i)=frac4(iz,i)+1   ! indice nb de points incertain
           endif
enddo   


endif

!print *, frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
!                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i)

 do  iz=1,alt  
 fracttot=frac1(iz,i)+frac2(iz,i)+frac3(iz,i)+frac4(iz,i)+frac5(iz,i)+       &
          frac6(iz,i)+frac7(iz,i) !check = 1
  if (fracttot.ne.1) then
     print *, "fraction error"
     print *, 'instant',fracttot,frac1(iz,i),frac2(iz,i),frac3(iz,i),        &
                        frac4(iz,i),frac5(iz,i),frac6(iz,i),frac7(iz,i), i,iz,var(iz,i)
!     stop
  endif
enddo


end subroutine fraction_subgrid2_8km_delta
!----------------------------------------------------------------------------!


subroutine fraction_subgrid3_8km(seuilsnrlow,seuilsnrhigh,altvar,var,var1,ind1,var2,ind2,frac,i,alt,nprofs,switch,switch2)

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid,seuilsnrhigh,seuilsnrlow
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt+1)  ::  altvar
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8


delta = 0
toplvlcloud=0; toplvlsat=0



B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8

do iz=seuilsnrhigh,alt   !=18 for cfmip2

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif          
enddo   


do iz=1,seuilsnrlow  ! =5 for cfmip2
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=seuilsnrlow+1,seuilsnrhigh-1
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -777 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3).and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ))then
                  frac(iz,i)=frac(iz,i)+2   
          endif
          if((var(iz,i).ge.SeuilSrCloud3).and.(delta.lt.SeuilDeltAtb)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -777 (surface elevation)            
          endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
             (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


endif


end subroutine fraction_subgrid3_8km
!----------------------------------------------------------------------------!


subroutine fraction_subgrid3_8km_delta(var,var1,ind1,var2,ind2,frac,i,alt,nprofs,switch,switch2)

  implicit none
  integer  ::  i,iz,nprofs,alt,toplow,topmid
  real*4  ::  fracttot
  real  ::  toplvlcloud,toplvlsat
  real, parameter  ::  SeuilStrat = 30.  
  real, parameter  ::  SeuilClearSr = 1.2  
  real, parameter  ::   SeuilSatSr = 0.01  
  real, parameter  ::    SeuilSrCloud = 5.
  real, parameter  ::    SeuilSrCloud2 = 15
  real  ::  SeuilSrCloud3
  real,parameter  ::    SeuilDeltAtb = 2.5e-03   
  real*4  ::  delta    ! delta atb = atb-atbmol
  character  ::  switch*5,switch2*6
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8


delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2


!print *, 'switch2',trim(switch2)


if(trim(switch).eq.'day')then

SeuilSrCloud3 = SeuilSrCloud

!print *, 'test'

! uncert+4 nan+1 surf+6 rejec+7 clear+2 cld+3 sat+8

do iz=18,alt

     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif          
enddo   


do iz=1,5
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif     
enddo   


B3:   do iz=1,toplow-1
      if(var(iz,i).gt.SeuilStrat)then
      SeuilSrCloud3 = SeuilSrCloud2 
      exit B3
      endif
   enddo B3

do iz=6,17
 
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then !.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then !.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


!
elseif(trim(switch).eq.'night')then


do iz=1,alt

SeuilSrCloud3 = SeuilSrCloud
     
     delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+7      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud3) )then !.and.(delta.ge.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch2).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8    ! pts fully attenuated
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts fully attenuated over the cloud ==> clear
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! 1st pt fully attenuated ==> cloud
               endif
elseif(trim(switch2).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).ge.(SeuilSatSr)) )then!.or. &
          !  ((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr).and.&
          !  (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud3).and.   &
          !  (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud3).and.(var(iz,i).ge.SeuilClearSr) )then!.and.&
          !   (delta.ge.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif
enddo   


endif


end subroutine fraction_subgrid3_8km_delta
!----------------------------------------------------------------------------!








!----------------------------------------------------------------------------!
! *** FRACTION_SUBGRID3 *** This routine calculate the cloudy,clear,         !
!                           saturated... fraction in order to be used by     !
!                           matlab (flag 1/2/3/4/5/6)                        !
!----------------------------------------------------------------------------!
! i        : loop index on the profil                                        !
! nprofs   : number of profil                                                !
! alt      : number of the altitude boxes                                    !
! var1, var2   : verticaly averaged observed variables (atb & mol moy)       !
! var    : vertically averaged SR                                            ! 
! ind1   : number of occurence of obs variables contained in var1            !
! ind2   : number of occurence of obs variables contained in var2            ! 
!----------------------------------------------------------------------------!
! frac   : fraction of fully, cloudy, clear, uncert, nan, se point           !
!          cloud=6, clear=2, fully_att=1, uncert=3, nan=4, SE=5              !
!----------------------------------------------------------------------------!
! iz       : loop index on the altitude grid                                 !
! toplvlcloud  : first cloudy point encountered from the top to the ground   !
! toplvlsat    : first fully attenuated point encountered from the top to    !
!                the ground                                                  ! 
! SeuilClearSr : clear treshold detection = 1.2                              !
! SeuilSatSr   : fully attenuated treshold detection = 0.01                  !
! SeuilSrCloud : cloudy treshold detection = 3                               !
! SeuilDeltAtb : delta atb treshold = 1.4e-03 km-1 sr-1                      !
! delta        : atb - mol                                                   !
!----------------------------------------------------------------------------!
!                                                                            !
! ex : call fraction_subgrid3(srmoy,pr2moy,indice,molmoy,indicem,fractot,i,  !
!                             altmax,it)                                     !
!----------------------------------------------------------------------------!
subroutine fraction_subgrid3(var,var1,ind1,var2,ind2,frac,i,alt,nprofs,switch)

  implicit none
  integer  ::  i,iz,nprofs,alt
  character  ::  switch*6

  real, parameter  ::  SeuilClearSr = 1.2      ! Seuil de saturation en atb en attendant flag
  real, parameter  ::   SeuilSatSr = 0.01     ! seuil de saturation en ATB=valeur ATBmol à 30km
  real, parameter  ::    SeuilSrCloud = 5.    ! seuil detection nuageuse sr
  real,parameter  ::    SeuilCrCloud = 0.   ! seuil detection nuageuse  cr>0.6, ne filtre pas les gros aerosols (atbd calipso)
  real,parameter  ::    SeuilDeltAtb = 2.5e-03    ! seuil détection unclassify
  real*4  ::  delta    ! delta atb = atb-atbmol
  integer  ::  toplvlcloud,toplvlsat
  real*4,dimension(alt,nprofs)  ::  var
  real*4,dimension(alt,nprofs)  ::  var1,var2
  real*4,dimension(alt,nprofs)  ::  ind1,ind2
  real*4,dimension(alt,nprofs)  ::  frac!,frac7

!flag 1= satfraction, 2=clear, 3=uncert, 5=nan, 4=SE, 6=cloud



delta = 0

toplvlcloud=0; toplvlsat=0

B1 : do iz=alt,1,-1
   if((var1(iz,i).eq.(-9999)).or.(var(iz,i).eq.(-888)).or.(var(iz,i).eq.(-777))) exit B1
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
     if ( (var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb) ) then
        toplvlcloud=iz     ! lvl from first cloudy point
        exit B1
     endif
  enddo B1

 B2 : do iz=alt,1,-1
    if(iz.lt.toplvlcloud)then
   if((var1(iz,i).eq.(-9999)).or.(var(iz,i).eq.(-888)).or.(var(iz,i).eq.(-777))) exit B2

       if ((var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.             &
          (var(iz,i).lt.(SeuilSatSr))) then 
        toplvlsat=iz    ! lvl from first saturated point
        exit B2
       endif
    endif
  enddo B2



  do iz=1,alt
   if((var1(iz,i).eq.(-9999)).or.(var(iz,i).eq.(-888)).or.(var(iz,i).eq.(-777)))then
   delta=0
   else
   delta= (var1(iz,i)/ind1(iz,i)) - (var2(iz,i)/ind2(iz,i))
   endif
           if ( (var(iz,i).eq.(-9999) )  )then! .or.((var(1,iz,i).ne.(-888)).and.(var(1,iz,i).lt.(SeuilSatSr2))) )  then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -9999 
           endif
           if (var(iz,i).eq.(-888)) then
              frac(iz,i)=frac(iz,i)+6      ! indice nb de -888 (surface elevation) 
           endif
           if (var(iz,i).eq.(-777)) then
              frac(iz,i)=frac(iz,i)+1      ! indice nb de -888 (surface elevation) 
           endif
           
           if ((var(iz,i).ge.SeuilSrCloud).and.(delta.gt.SeuilDeltAtb)) then
                  frac(iz,i)=frac(iz,i)+3   ! indice nb de points nuageux ds boite
           endif
           if ( (var(iz,i).ne.(-9999)).and.(var(iz,i).ne.(-888)).and.(var(iz,i).ne.(-777)).and.        &
                (var(iz,i).lt.(SeuilSatSr)) ) then 
if(trim(switch).eq.'cloudy')then
               if(iz.lt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
                  elseif(iz.eq.toplvlsat)then
                     frac(iz,i)=frac(iz,i)+3    ! pts nuageux
               endif
elseif(trim(switch).eq.'sat')then
               if(iz.le.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+8  !
               elseif(iz.gt.toplvlsat)then
                  frac(iz,i)=frac(iz,i)+2    ! pts clair
               endif
endif
           endif
          if((var(iz,i).lt.SeuilClearSr).and.(var(iz,i).gt.(SeuilSatSr)).or. &
            ((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
            (delta.lt.SeuilDeltAtb) ).or.((var(iz,i).ge.SeuilSrCloud).and.   &
            (delta.lt.SeuilDeltAtb))) then
                  frac(iz,i)=frac(iz,i)+2   
           
           endif
          if((var(iz,i).lt.SeuilSrCloud).and.(var(iz,i).gt.SeuilClearSr).and.&
             (delta.gt.SeuilDeltAtb))then
                  frac(iz,i)=frac(iz,i)+4   ! indice nb de points incertain
           endif

!!$if((i.ge.35618).and.(i.le.35629))then
!!$print *,i
!!$if(iz.lt.4)then
!!$    print *, iz,i, delta,toplvlcloud,toplvlsat
!!$    print *, var(iz,i),frac(iz,i)
!!$endif
!!$
!!$endif


  enddo   
         

end subroutine fraction_subgrid3
!----------------------------------------------------------------------------!
