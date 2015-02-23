!     @process extchk
!
subroutine take_nonloc_deriv_kq_k(dir,gs, rk, psi, nanl,nloc_deriv_kq,&
     nloc_deriv_k,dpsi, crys) 

  use all_to_all_module

  include 'use.h'
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  include 'all_to_all.h'
  include 'flibcalls.ph'


  !
  !     INPUT:
  !     -----

  type(parallel_gspace) :: gs 

  type(crystal) :: crys

  integer :: &
       nanl, &              ! number of nonlocal projectors
       dir                  ! direction into which to take the derivative

  real(dp) :: rk(3)

  complex(dp) :: &
       psi(gs%length),&     ! the wavefn to which dVnl/dk is applied
       nloc_deriv_kq(gs%length+nanl,nanl,2,3),&!derivative of nonlocal, kq 
       nloc_deriv_k(gs%length+nanl,nanl,2,3)  ! derivative of nonlocal, k
  ! 
  !     OUTPUT:
  !     ------
  !
  complex(dp) :: dpsi(gs%length)  ! the dVnl/dk * |psi>

  !
  !     DESCRIPTION:
  !     -----------
  !
  !     * applies the nonlocal potential derivative operator to the wave
  !       function 
  !
  !     * computes derivative wrt spatial direction: dir
  !       we compute data(g)*(k+G), so set the k-point of the gspace
  !       to zero if you just want the gradient!
  !
  !     (1996) Bernd Pfrommer
  !

  !
  !     ----------- local variables -----------------------------------
  !
  complex(dp) dcone, dczero, dcmone
  parameter(dcone=(1.d0,0.d0), dczero=(0.d0,0.d0), dcmone=(-1.d0,0.d0))
  complex(dp), allocatable :: work(:)
  integer :: k

  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4),fftn(4), ffth(4),igv3, irod,iord,igs

  real(dp) :: gv(3)
  !
  !
  if((dir.lt.1).or.(dir.gt.3)) then
     write(9,*) ' derivative_gspace: illegal direction: ',dir
     call mystop
  endif
  !
  !     --------- apply the nonlocal part of the hamiltonian ----------
  !

  dpsi=0

  !      goto 200                  ! use this to disable nonlocal part...

  if(nanl.gt.0) then
     allocate(work(nanl))
     !
     !        the +q part
     !
     call mzgemm('C', 'N', nanl, 1, gs%length,dcone,nloc_deriv_k(1,1,2,dir),&
          gs%length+nanl,psi(1),gs%length,dczero, work(1),nanl)

     do k=1,nanl
        work(k)=work(k)*dble(nloc_deriv_k(gs%length+k,1,1,1))
     end do

     call all_sum_all(work,nanl)
     !
     !        linearly combine the projectors to get +q part
     !     
     call mzgemm('N','N',gs%length,1,nanl,dcone,nloc_deriv_kq(1,1,2,dir),&
          gs%length+nanl, work(1), nanl,dczero, dpsi(1),gs%length)
     !
     !        the -q part
     !
     call mzgemm('C', 'N', nanl, 1, gs%length, dcone, nloc_deriv_k(1,1,1,dir),&
          gs%length+nanl, psi(1),gs%length,dczero,work(1),nanl)

     do k=1,nanl
        work(k)=work(k)*dble(nloc_deriv_k(gs%length+k,1,1,1))
     end do


     call all_sum_all(work,nanl)
     !
     !        linearly combine the projectors to get -q part
     !     
     call mzgemm('N','N',gs%length,1,nanl,dcmone,nloc_deriv_kq(1,1,1,dir),&
          gs%length+nanl, work(1), nanl,dcone, dpsi(1),gs%length); 

     deallocate(work)
  endif

200 continue
  !
  !     --------- take the gradient ----------
  !
  fftn(1:3)=gs%fftsize(1:3)
  fftn(4)=gs%fftsize(3)
  ffth(:)=fftn(:)/2
  igs = 0 
  do iord=1, gs%lorder    ! loop through x/y gspace 
     irod=gs%order(1,iord)
     igv(1)= irod/gs%fftsize(2)
     igv(2)= mod(irod,gs%fftsize(2))
     igv(3)= gs%order(2,iord)
     igv(4)= gs%order(3,iord)
     igv(:)= mod(igv(:)+ffth(:),fftn(:))-ffth(:) 
     gv(1:2)=dble(igv(1:2))+rk(1:2)
     do igv3=igv(3),igv(4)     ! loop over z axis
        gv(3)=dble(igv3)+rk(3)
        igs = igs + 1
        dpsi(igs)=dpsi(igs)+psi(igs)*(crys%bvec(dir,1)*gv(1)+&
             crys%bvec(dir,2)*gv(2)+crys%bvec(dir,3)*gv(3))
     end do
  end do
  return

end subroutine take_nonloc_deriv_kq_k
