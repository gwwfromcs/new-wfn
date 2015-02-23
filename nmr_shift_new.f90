subroutine nmr_shift_new(pw_params,crys,gs,chi_g, &
     dia_corr_shift,para_corr_shift)

  use all_to_all_module
  include 'use.h'
  implicit none             ! implicit? Just say no!
  include 'interface.h'
  include 'all_to_all.h'
  !
  !     computes NMR shift at the nuclei
  !
  !
  !     INPUT:
  !     -----

  type (crystal)         :: crys       ! for the lattice vectors
  type (parallel_gspace) :: gs         ! potential gspace
  type (pw_parameter)    :: pw_params  ! plane-wave parameters

  complex(dp) :: chi_g(gs%length,3,3)  ! chi in gspace

  integer, parameter :: current=3, magnet=3

  real(dp)  :: dia_corr_shift(magnet,current,crys%ntype,crys%mxdatm) 
  real(dp)  :: para_corr_shift(magnet,current,crys%ntype,crys%mxdatm) 


  !
  !
  !     1996 Bernd Pfrommer, UCB
  !
  !

  !     ---------------- local variables ----------------------------------

  real(dp) :: ft(3,3)

  integer :: i,iat,nt,n

  real(dp)              :: r(3),t0, gimmetime
  real(dp), allocatable :: ftbare(:,:,:)
  external gimmetime
  !
  !     ----- variables for the gspace loop ---------
  !
  integer :: igv(4),fftn(4), ffth(4),igv3, irod,iord,igs

  real(dp) :: gv(3)  ! if required

  t0=gimmetime()

  allocate(ftbare(3,3,crys%mxdatm*crys%ntype))

  write(9,50)

  fftn(1:3)=gs%fftsize(1:3)
  fftn(4)=gs%fftsize(3)
  ffth(:)=fftn(:)/2

  n=0
  do nt=1,crys%ntype
     do iat=1,crys%natom(nt)
        n=n+1
        ft=0.d0
        r=crys%rat(:,iat,nt)
        igs = 0 
        do iord=1, gs%lorder ! loop through x/y gspace 
           irod=gs%order(1,iord)
           igv(1)= irod/gs%fftsize(2)
           igv(2)= mod(irod,gs%fftsize(2))
           igv(3)= gs%order(2,iord)
           igv(4)= gs%order(3,iord)
           igv(:)= mod(igv(:)+ffth(:),fftn(:))-ffth(:) 
           gv(1:2)=dble(igv(1:2))+gs%rk(1:2)
           do igv3=igv(3),igv(4) ! loop over z axis
              gv(3)=dble(igv3)+gs%rk(3)
              igs = igs + 1
              ft=ft+real(chi_g(igs,:,:)*exp(cmplx(0.d0,r(1)*gv(1)+r(2)*gv(2)+&
                   r(3)*gv(3),kind=8)),dp)
           end do
        end do
        call all_sum_all(ft,9) ! 3x3 matrix sum across all processors
        ftbare(:,:,n)=ft
        ft=-1.d0*transpose(ft)/(137.036**2.d0)*1.d6
        write(9,100) crys%nameat(nt),iat,crys%rat(:,iat,nt)/pi2,&
             (ft(1,1)+ft(2,2)+ft(3,3))/3.d0, ft
     end do
  end do

  if(iand(pw_params%output(1),8).eq.8) write(9,940) gimmetime()-t0

  write(9,57)
  n=0
  do nt=1,crys%ntype
     do iat=1,crys%natom(nt)
        n=n+1
        ft=-transpose(ftbare(:,:,n))/(137.036**2.d0)*1.d6
        ft=ft+dia_corr_shift(:,:,nt,iat)
        write(9,101) crys%nameat(nt),iat,crys%rat(:,iat,nt)/pi2,&
             (ft(1,1)+ft(2,2)+ft(3,3))/3.d0,ft
     end do
  end do

  write(9,58)
  n=0
  do nt=1,crys%ntype
     do iat=1,crys%natom(nt)
        n=n+1
        ft = -ftbare(:,:,n)/(137.036**2.d0)*1.d6
        ft = transpose(ft)
        ft = ft+para_corr_shift(:,:,nt,iat)
        ft = ft+dia_corr_shift(:,:,nt,iat)
        write(9,102) crys%nameat(nt),iat,crys%rat(:,iat,nt)/pi2,&
             (ft(1,1)+ft(2,2)+ft(3,3))/3.d0,ft
     end do
  end do

  call myflush(9)

  deallocate(ftbare)

  return

50 format(/' -------- BARE MAGNETIC SHIFTS -----------'//,&
        ' ATOM  NR.           POSITION',21x,'SHIFT [ppm] '/)
57 format(/' ---- DIA-CORRECTED MAGNETIC SHIFTS ------'//,&
        ' ATOM  NR.           POSITION',21x,'SHIFT [ppm] '/)
58 format(/' ---DIA&PARA-CORRECTED MAGNETIC SHIFTS ----'//,&
        ' ATOM  NR.           POSITION',21x,'SHIFT [ppm] '/)

100 format('BARE ',1x,a2,3x,i3,3f12.6,3x,f12.4,//3(9x,3f12.4/))
101 format('DIA  ',1x,a2,3x,i3,3f12.6,3x,f12.4,//3(9x,3f12.4/))
102 format('PARA ',1x,a2,3x,i3,3f12.6,3x,f12.4,//3(9x,3f12.4/))


940 format(' TIME FOR NMR FOURIER SUM:',f12.3)

end subroutine nmr_shift_new




