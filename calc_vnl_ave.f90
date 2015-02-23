!*
subroutine calc_vnl_ave(vnl_ave,ffts,gs,pspot,crys,vnl_atom,nt ) 
  use all_to_all_module  
  include 'use.h'  
  implicit none               ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  

  !
  !     INPUT:  
  !     -----
  !
  integer nt
  type(parallel_gspace), intent(in) :: gs   
  type(fft_struc), intent(in) :: ffts  
  type(pseudo_potential), intent(in) :: pspot
  type(crystal), intent(in):: crys  
  !
  !     OUTPUT:
  !     ------
  !
  type(complex_gspace_array), intent(inout) :: &
   vnl_ave  ! average nonlocal in real space for TF mixing

  !  local variables

  integer i,n,j

  real(dp) a,b,vqj,qj,xn,delql,q2vn,q2vp,q2vm

  real(dp) :: vnl_atom(pspot%mxdlqp ) ,temp(2*gs%r_size),s,dasum,dzasum

  delql=pspot%delqnl(nt)
        
      do j = 1, gs%length 
        !
        !          interpolate vnl_atom
        !
        qj = sqrt(gs%ekin(j))  
        if (qj > dzero) then  
           xn = qj / delql + dtwo  
           n = xn + dhalf  

           if (n < pspot%nqnl(nt) ) then  
              if (n <= 3) n = 4  
              xn = xn - real(n, dp)  
             vqj = vnl_atom(n) * (done + xn) * (done - xn) + dhalf * (vnl_atom(n + 1) * &
                (done + xn) - vnl_atom(n - 1) * (done - xn)) * xn
  
            vqj=vqj/crys%vcell
           !
           !             sum up the charge density
           !

              vnl_ave%data(j,1,1) = vnl_ave%data(j,1,1) + vqj * &
                   gs%struc(nt, j)
           end if
        end if
      end do

  return

end subroutine calc_vnl_ave
