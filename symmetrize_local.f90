!
subroutine symmetrize_local(ng, data, ns, inds, mstar, phase, &
     work, gsinv)
  !
  !     enforces correct symmetry, and also that data(G)=conjg(data(-G))
  !
  !     1995 by Bernd Pfrommer, while at UCB
  !
  use constants
  implicit none    ! implicit? Just say no!
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: &
       ng, &
       ns, &                   ! number of stars
       gsinv(ng), &            ! inversion information
       inds(ng), &             ! star index for each gvector
       mstar(ns)               ! size of stars
  complex(dp), intent(in) :: &
       phase(ng)               ! phase factors for mixing gvectors
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  complex(dp), intent(inout) :: &
       data(ng)                ! the data to be symmetrized
  !
  !     WORK ARRAYS:
  !     -----------
  !
  complex(dp) :: work(ns)
  !  
  !     ------------------ local variables ------------------------
  !
  integer :: i, ii  
  complex(dp) :: dcmp  
  !
  !     ------------------ first symmetrize with symmetry operations -----
  !
  do i = 1, ns  
     work(i) = zzero
  end do
  do i = 1, ng  
     work(inds(i)) = work(inds(i)) + data(i) * phase(i)  
     !         if(inds(i).eq.42) print*,i,inds(i),
     !     $        data(i),phase(i),data(i)*phase(i)
  end do
  do i = 1, ns  
     work(i) = work(i) / real(mstar(i), dp)  
  end do
  !      write(9,*) 'symmetrization:'
  do i = 1, ng  
     !         dcmp=data(i)
     data(i) = work(inds(i)) * conjg(phase(i))  
     !         if(inds(i).eq.42) write(9,*)  i,inds(i), dcmp,data(i),
     !     $        work(inds(i)), phase(i)
  end do
  !
  !     now symmetrize data(G)=conjg(data(-G))
  !
  do i = 1, ng  
     ii = gsinv(i)  
     if (ii > ng .or. ii <= 0) then  
        write(9, *) 'symmetrize_local: inversion array is incorrect'  
        call mystop  
     end if
     dcmp = dhalf * (data(i) + conjg(data(ii)))  
     data(i) = cmplx(real(dcmp, dp), aimag(dcmp), dp)  
     data(ii) = conjg(data(i))  
  end do

  return  

end subroutine symmetrize_local
