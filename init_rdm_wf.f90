  subroutine init_rdm_wf(psi,iranm,neig,ffts,gspace,crys)
  !
  ! calculates random wavefunction for 1 kpt 
  !
  use all_to_all_module  
  include 'use.h' 
  implicit NONE
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'
  !
  !     INPUT:
  !     -----
  !
  type(fft_struc), intent(inout) :: ffts             ! information for the FFT
  type(crystal), intent(in) :: crys                  ! crystal structure
  type(parallel_gspace), intent(in) :: gspace
  integer, intent(in) :: iranm,neig
  !
  !     OUTPUT:
  !     -------
  !
  complex(dp) psi(*)  
  !
  !     --------------- local variables -------------------------
  !
  real(dp)  x1,x2,ran1,qk(3),qcar(3)
  integer len,iwfn,i,m,n1,n2,n3,j,k,ir,rlen
  complex(dp) temp(gspace%r_size)
  !
  !     --------------------------------------------------------------
  !
  len=gspace%length
  rlen=gspace%r_size
  n1=ffts%fftsize(1)
  n2=ffts%fftsize(2)
  n3=ffts%fftsize(3)

  do iwfn=0,neig-1
!    do i=1,n1 
!     do j=1,n2
!      do k=1,n3
     do i=1,rlen
      x1=ran1(iranm)
      x2=ran1(iranm)
!      ffts%rspacebuf(i+(j-1)*n1+(k-1)*n1*n2)=cmplx(x1-0.5d0,x2-0.5d0,dp)
      temp(i)=cmplx(x1-0.5,x2-0.5,dp)
    enddo
!    end do
!   end do

   call fourier_transform(1,ffts,gspace,psi(iwfn*len+1),temp(1),1 ) 

  enddo
 
!*************************************************
!**** end generate the initial wavefunction from random
!*************************************************
  do m=1,neig
    call orth_comp(psi(1+(m-1)*len),psi,m-1,len,neig,crys%vcell,1)
  enddo

  return
  end subroutine init_rdm_wf
  !
  !
  !  Gram-Scmidt orthogonalization  *******************************
  !
  !
  subroutine orth_comp(psi_mt,psi,mt,len,neig,vol,isign)

  use all_to_all_module  
  include 'use.h' 
  implicit NONE
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: len,neig,mt,isign
  complex(dp), intent(in) :: psi(len,neig)
  real(dp), intent(in) :: vol 
  !
  !     INPUT/OUTPUT:
  !     -------------
  !
  complex(dp) psi_mt(len)
  !
  !     --------------- local variables -------------------------
  !
  complex(dp) cc,ccc,cc1,cc2
  complex(dp) sumdumc(mt)
  integer mi,m,i      
  real(dp) s

  if(isign.eq.3.or.mt.le.1) then  

    mi=mt

    if(mt .gt. 0) then  !  le.0) goto 25

      do m=mi,mt
  
        cc=zzero

        do i=1,len
         cc=cc+psi_mt(i)*conjg(psi(i,m))
        enddo

        call all_sum_all(cc)
        cc=cc*vol

        do i=1,len
         psi_mt(i)=psi_mt(i)-cc*psi(i,m)
        enddo

      end do

    end if

  else
 
    call mzgemv('c',len,mt,zone,psi(1,1),len,psi_mt(1),1,zzero,sumdumc(1),1)

    call all_sum_all(sumdumc,mt)

    do i = 1,mt
      sumdumc(i) = sumdumc(i)*vol
    enddo

    call mzgemv('n',len,mt,zmone,psi(1,1),len,sumdumc(1),1,zone,psi_mt(1),1)

  endif 

!**************************************************
!****  renormalize
!**************************************************

  if(isign.eq.1) then

    s=dzero

    do i=1,len
      s=s+abs(psi_mt(i))**2
    enddo

    call all_sum_all(s)

    if(s.lt.1.D-200) then
      write(6,*) "test, warning, s=", s
    endif

    s=dsqrt(done/(s*vol))

    do i=1,len
      psi_mt(i)=s*psi_mt(i)
    enddo

  endif

  return
  end subroutine orth_comp
  !
  !  random fnumber genrator
  !
   FUNCTION RAN1(IDUM)
   include 'use.h'
   implicit double precision (a-h,o-z)
   DIMENSION R(97)
   PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
   PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
   PARAMETER (M3=243000,IA3=4561,IC3=51349)
   real(dp) ran1
   save ix1,ix2,ix3
   save R,IFF
   DATA IFF /0/

   IF (IDUM.LT.0.OR.IFF.EQ.0) THEN

     IFF=1
     IX1=MOD(IC1-IDUM,M1)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX2=MOD(IX1,M2)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX3=MOD(IX1,M3)

     do J=1,97
       IX1=MOD(IA1*IX1+IC1,M1)
       IX2=MOD(IA2*IX2+IC2,M2)
       R(J)=(real(IX1,dp)+real(IX2,dp)*RM2)*RM1
     end do
     IDUM=1

   ENDIF

   IX1=MOD(IA1*IX1+IC1,M1)
   IX2=MOD(IA2*IX2+IC2,M2)
   IX3=MOD(IA3*IX3+IC3,M3)
   J=1+(97*IX3)/M3
   IF(J.GT.97.OR.J.LT.1) call mystop( ' error in init_rdm_wf.f90' )
   RAN1=R(J)
   R(J)=(real(IX1,dp)+real(IX2,dp)*RM2)*RM1

   RETURN
   END




