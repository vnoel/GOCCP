! It would be nice to move all HDF reading routines here, e.g. sdsread
! It might not be doable since variables are allocated within these subroutines


!----------------------------------------------------------------------------!
! ***** NPROF ***** This routine find the hdf file number of profiles        !
!----------------------------------------------------------------------------!
subroutine nprof(filename,nb,nprofs)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: ret
	integer 	:: sd_id, sds_id
        integer         :: npts, nprofs, nb
	integer*4	:: rank, istat, attributes, num_type !index
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
 
     
	sd_id = sfstart (filename,	DFACC_READ)
  
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)
        
        ! Rtrieving the dimesions
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)
print *, 'le nombre de profil du fichier est :', nprofs

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)

end subroutine nprof
!----------------------------------------------------------------------------!
