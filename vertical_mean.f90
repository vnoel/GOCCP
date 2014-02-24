! vertical_mean.f90

!----------------------------------------------------------------------------!
! *** VERTICAL_MEAN_CHIM *** This routine calculate the vertical mean for    !
!                             the CHIMERE grid                               !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! ilid   : loop index of the lidar altitude (583)                            !
! nprofs : number of profil                                                  !
! j      : loop index on latitude grid                                       !
! k      : loop index on longitude grid                                      ! 
! lat    : number of the latitude boxes                                      !
! lon    : number of the longitude boxes                                     !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
! nvar   : switch to chose cloudy or sat mode                                !
! var3   : attenuated backscatter                                            !
! var1   : observed variable to average                                      !
! ind    : number of occurence of obs variables contained in var1            !
!----------------------------------------------------------------------------!
! var2   : averaged variable                                                 !
!----------------------------------------------------------------------------!
!                                                                            !
! ex  call vertical_mean_chim(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,     !
!                             altitude,j,k,latmax,lonmax,1)                  !
!----------------------------------------------------------------------------!
subroutine vertical_mean_chim(var1,var2,var3,ind,i,iz,ilid,nprofs,alt,alt2,j,&
                              k,lat,lon,nvar) 

  implicit none
  integer  ::  nvar
  integer  ::  i,nprofs,j,k,iz,alt,ilid,lat,lon
  integer(kind=2)  ::  alt2
  real*4,dimension(alt2,nprofs)  ::  var1,var3
  real*4,dimension(lat,lon,alt)  ::  var2
  real*4,dimension(lat,lon,alt)  ::  ind
  

  if (nvar==1) then
     ! Vertical average for atb
     if ( (var1(ilid,i).ne.(-9999)).and.(var1(ilid,i).lt.1) ) then
        var2(j,k,iz)=var2(j,k,iz)+var1(ilid,i)
        ind(j,k,iz)=ind(j,k,iz)+1 
     endif

  else
      ! Vertical average for perp
     if(nvar==2) then
        if ( (var1(ilid,i).ne.(-9999)).and.(var3(ilid,i).ne.(-9999)) .and.   &
             (var3(ilid,i).lt.1) )then 
                var2(j,k,iz)=var2(j,k,iz)+(var3(ilid,i)-var1(ilid,i))
                ind(j,k,iz)=ind(j,k,iz)+1 
        endif
     else
      ! Other vertical average
        if ( var1(ilid,i).ne.(-9999) ) then                            
           var2(j,k,iz)=var2(j,k,iz)+var1(ilid,i)
           ind(j,k,iz)=ind(j,k,iz)+1 
        endif
     endif

  endif
                    
end subroutine vertical_mean_chim
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! *** VERTICAL_MEAN *** This routine calculate the vertical mean for other   !
!                       grid                                                 !
!----------------------------------------------------------------------------!
! i      : loop index on the profil                                          !
! ilid   : loop index of the lidar altitude (583)                            !
! nprofs : number of profil                                                  !
! iz     : loop index on the altitude grid                                   !
! alt    : number of the altitude boxes                                      !
! alt2   : number of level of the lidar variables                            !
! nvar   : switch to chose cloudy or sat mode                                !
! var3   : attenuated backscatter                                            !
! var1   : observed variable to average                                      !
! ind    : number of occurence of obs variables contained in var1            ! 
!----------------------------------------------------------------------------!
! var2   : averaged variable                                                 !
!----------------------------------------------------------------------------!
!                                                                            !
! ex :  call vertical_mean(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,        !
!                          altitude,1)                                       !
!----------------------------------------------------------------------------!
subroutine vertical_mean(var1,var2,var3,ind,i,iz,ilid,nprofs,alt,alt2,nvar)
  
  implicit none
  integer  ::  nvar
  integer  ::  i,iz,nprofs,alt,ilid
  integer(kind=2)  ::  alt2
  real*4,dimension(alt2,nprofs)  ::  var1,var3
  real*4,dimension(alt,nprofs)  ::  var2
  real*4,dimension(alt,nprofs)  ::  ind
  

  if (nvar==1) then ! Vertical average for atb

     if ( (var1(ilid,i).ne.(-9999.)).and.(var1(ilid,i).lt.1).and.(var1(ilid,i).ne.(-888.)).and.(var1(ilid,i).ne.(-777.))) then         
        var2(iz,i)=var2(iz,i)+var1(ilid,i)
        ind(iz,i)=ind(iz,i)+1 
     endif

  elseif(nvar==2) then ! Vertical average for perp

     if ( (var1(ilid,i).ne.(-9999.)).and.(var3(ilid,i).ne.(-9999.)) .and.   &
          (var3(ilid,i).lt.1).and.(var1(ilid,i).ne.(-888.)) .and.  &
          (var1(ilid,i).ne.(-777.)) )then 
         var2(iz,i)=var2(iz,i)+(var3(ilid,i)-var1(ilid,i))
         ind(iz,i)=ind(iz,i)+1 
     endif

  elseif(nvar==3) then ! Other vertical average 
   
     if (( var1(ilid,i).ne.(-9999.).and.(var1(ilid,i).ne.(-888.)) ) .and.   &
        (var3(ilid,i).lt.1).and.(var3(ilid,i).ne.(-9999.)).and.  &
        (var1(ilid,i).ne.(-777.)) ) then                            
           var2(iz,i)=var2(iz,i)+var1(ilid,i)
           ind(iz,i)=ind(iz,i)+1 
     endif

  elseif(nvar==4) then ! TEMP vertical average 
  
     if (( var1(ilid,i).ne.(-9999.).and.(var1(ilid,i).ne.(-888.)) ) .and.   &
        (var1(ilid,i).ne.(-777.)) ) then                            
           var2(iz,i)=var2(iz,i)+var1(ilid,i)
           ind(iz,i)=ind(iz,i)+1 
     endif


  endif

                    
end subroutine vertical_mean
!----------------------------------------------------------------------------!


subroutine vertical_mean_hori(var1,var2,var3,ind,i,iz,ilid,nprofs,alt,alt2,nvar,profavg,profmax)
  
  implicit none
  integer  ::  nvar, n
  integer  ::  i,iz,nprofs,alt,ilid,profavg,profmax,it
  integer(kind=2)  ::  alt2
  real*4  ::  var2m, indm
  real*4,dimension(alt2,nprofs)  ::  var1,var3
  real*4,dimension(alt,nprofs)  ::  var2
  real*4,dimension(alt,nprofs)  ::  ind
  
  ! change vnoel 20140224 - same as in atb_mol
  it = nprofs


  if (nvar==1) then
  
    if (i.le.profavg)then
       var2m=0
       indm=0
      do n = 0,profmax   
       ! Vertical average for atb
         if ( (var1(ilid,1+n).ne.(-9999)).and.(var1(ilid,1+n).lt.1) ) then         
            var2m=var2m+var1(ilid,1+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif


    elseif ( (i.gt.profavg).and.(i.le.(it-(profavg+1))) ) then   
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,i-profavg+n).ne.(-9999)).and.(var1(ilid,i-profavg+n).lt.1) ) then         
            var2m=var2m+var1(ilid,i-profavg+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    else
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,it-profmax+n).ne.(-9999)).and.(var1(ilid,it-profmax+n).lt.1) ) then         
            var2m=var2m+var1(ilid,it-profmax+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    endif
  

  else

     ! Vertical average for perp
     if(nvar==2) then
        if ( (var1(ilid,i).ne.(-9999)).and.(var3(ilid,i).ne.(-9999)) .and.   &
             (var3(ilid,i).lt.1) )then 
               var2(iz,i)=var2(iz,i)+(var3(ilid,i)-var1(ilid,i))
               ind(iz,i)=ind(iz,i)+1 
        endif
     else

    if (i.le.profavg)then
       var2m=0
       indm=0
      do n = 0,profmax   
       ! Vertical average for atb
         if ( (var1(ilid,1+n).ne.(-9999)) ) then         
            var2m=var2m+var1(ilid,1+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif


    elseif ( (i.gt.profavg).and.(i.le.(it-(profavg+1))) ) then   
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,i-profavg+n).ne.(-9999)) ) then         
            var2m=var2m+var1(ilid,i-profavg+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    else
       var2m=0
       indm=0
      do n = 0,profmax   
         ! Vertical average for atb
         if ( (var1(ilid,it-profmax+n).ne.(-9999)) ) then         
            var2m=var2m+var1(ilid,it-profmax+n)
            indm=indm+1 
         endif
      enddo
      if (var2m.gt.0) then
      var2(iz,i)=var2(iz,i)+(var2m/indm)
      ind(iz,i)=ind(iz,i)+1
      endif

    endif
  



     endif

  endif
                  
 

                    
end subroutine vertical_mean_hori
!----------------------------------------------------------------------------!
