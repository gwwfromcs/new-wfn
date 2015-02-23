!
subroutine apply_ham_wavefn (time,imonitor, hpsi, ham, psi, nwfn, work, &
     nbandsfft)

  use all_to_all_module  
  include 'use.h'  
  implicit none              ! implicit? Just say no!
  include 'interface.h'  
  include 'all_to_all.h'  
  include 'flibcalls.ph'  
  !
  !     INPUT:
  !     -----
  !
  type(hamiltonian), intent(inout) :: ham      ! hamiltonian at a given k-point
  complex(dp), intent(in) :: psi(*)         ! input wave functions
  integer, intent(in) :: imonitor(4), &        ! monitor flags
       nwfn, &                              ! number of wave functions
       nbandsfft
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  real(dp), intent(inout):: time(num_time)
  !
  !     OUTPUT:
  !     ------
  !
  complex(dp), intent(out) :: hpsi(*)  
  !
  !  applies the hamiltonian ham to the first nwfn wave functions contained
  !  in psi, and returns the answer in the first nwfn wave functions in
  !
  !
  !     WORK:
  !     -----
  ! work(nanl*nwfn)
  complex(dp), intent(inout) :: work(*) 
  !
  !     ------------- local variables --------------------------------
  !
  integer :: i, len, nanl, j, dnx, nrodst1, nrodst2,rlen,nmatom_tot,ico1,&
             ico2,isum_cnt,ia,nref,isum_cnt_st,iwfn,jj,ifn,n1,n2,n3,k,ir,&
             iatom,nfn, iwfn2
  real(dp) :: emod, ekinj, ekinmod (3), t0, t1, t2, n, opcount, t3
  !
  real(dp), external :: myderf, gimmetime  
  complex(dp) sumdum(ham%vnloc%nvecs),y,s
  integer natom_tot,icnt
  complex(dp), allocatable :: workr2_psi(:) 
  !
  !     apply local potential in realspace
  !
  ! length of gspace
  len = ham%gspace%length
  rlen=ham%gspace%r_size
  natom_tot=ham%pspot%natom_tot
  nanl = ham%vnloc%nvecs 

  n1=ham%fftstruc%fftsize(1)
  n2=ham%fftstruc%fftsize(2)
  n3=ham%fftstruc%fftsize(3)

  if (iand(imonitor(1), 8) == 8) t3 = gimmetime() 
  t0=gimmetime()

  if( ham%pspot%NL_rspace(1)) then

   allocate(workr2_psi(nbandsfft*rlen))

   do iwfn = 0, nwfn - 1, nbandsfft    

     nfn=min(nbandsfft, nwfn - iwfn)
     call fourier_transform(-1, ham%fftstruc, ham%gspace, psi(iwfn*len+1), &
              ham%fftstruc%rspacebuf(1), nfn)
  
     if (iand(imonitor(1), 8) == 8) call get_timing(time(1),t3)

     call mzcopy(nfn*rlen,ham%fftstruc%rspacebuf(1),1,workr2_psi(1),1)

     do iwfn2=0,nfn-1

     ham%fftstruc%rspacebuf(iwfn2*rlen+1:iwfn2*rlen+rlen)= &
        ham%fftstruc%rspacebuf(iwfn2*rlen+1:iwfn2*rlen + rlen)* &
              real(ham%vloc(1:rlen))

     sumdum = zzero 
     ico1=0
     ico2=1
     isum_cnt_st=1

     do ia=1,natom_tot

      nref=ham%pspot%numref(ia)
      
      do i=1,ham%pspot%nmap(ia)
        ico1=ico1+1
        y=ham%cphase(ico1)*workr2_psi(ham%pspot%indm(ico1)+iwfn2*rlen)  

        isum_cnt=isum_cnt_st

        select case(nref)

         case(8)
         sumdum(isum_cnt)=sumdum(isum_cnt)+ham%pspot%kb_rsp(ico2) *y
         sumdum(isum_cnt+1)=sumdum(isum_cnt+1)+ham%pspot%kb_rsp(ico2+1)*y      
         sumdum(isum_cnt+2)=sumdum(isum_cnt+2)+ham%pspot%kb_rsp(ico2+2)*y
         sumdum(isum_cnt+3)=sumdum(isum_cnt+3)+ham%pspot%kb_rsp(ico2+3)*y
         sumdum(isum_cnt+4)=sumdum(isum_cnt+4)+ham%pspot%kb_rsp(ico2+4)*y
         sumdum(isum_cnt+5)=sumdum(isum_cnt+5)+ham%pspot%kb_rsp(ico2+5)*y
         sumdum(isum_cnt+6)=sumdum(isum_cnt+6)+ham%pspot%kb_rsp(ico2+6)*y
         sumdum(isum_cnt+7)=sumdum(isum_cnt+7)+ham%pspot%kb_rsp(ico2+7)*y
         case(6)
         sumdum(isum_cnt)=sumdum(isum_cnt)+ham%pspot%kb_rsp(ico2)*y
         sumdum(isum_cnt+1)=sumdum(isum_cnt+1)+ham%pspot%kb_rsp(ico2+1)*y      
         sumdum(isum_cnt+2)=sumdum(isum_cnt+2)+ham%pspot%kb_rsp(ico2+2)*y
         sumdum(isum_cnt+3)=sumdum(isum_cnt+3)+ham%pspot%kb_rsp(ico2+3)*y
         sumdum(isum_cnt+4)=sumdum(isum_cnt+4)+ham%pspot%kb_rsp(ico2+4)*y
         sumdum(isum_cnt+5)=sumdum(isum_cnt+5)+ham%pspot%kb_rsp(ico2+5)*y
         case(5)
         sumdum(isum_cnt)=sumdum(isum_cnt)+ham%pspot%kb_rsp(ico2)*y
         sumdum(isum_cnt+1)=sumdum(isum_cnt+1)+ham%pspot%kb_rsp(ico2+1)*y      
         sumdum(isum_cnt+2)=sumdum(isum_cnt+2)+ham%pspot%kb_rsp(ico2+2)*y
         sumdum(isum_cnt+3)=sumdum(isum_cnt+3)+ham%pspot%kb_rsp(ico2+3)*y
         sumdum(isum_cnt+4)=sumdum(isum_cnt+4)+ham%pspot%kb_rsp(ico2+4)*y
         case(4)
         sumdum(isum_cnt)=sumdum(isum_cnt)+ham%pspot%kb_rsp(ico2)*y
         sumdum(isum_cnt+1)=sumdum(isum_cnt+1)+ham%pspot%kb_rsp(ico2+1)*y      
         sumdum(isum_cnt+2)=sumdum(isum_cnt+2)+ham%pspot%kb_rsp(ico2+2)*y
         sumdum(isum_cnt+3)=sumdum(isum_cnt+3)+ham%pspot%kb_rsp(ico2+3)*y
         case(3)
         sumdum(isum_cnt)=sumdum(isum_cnt)+ham%pspot%kb_rsp(ico2)*y
         sumdum(isum_cnt+1)=sumdum(isum_cnt+1)+ham%pspot%kb_rsp(ico2+1)*y      
         sumdum(isum_cnt+2)=sumdum(isum_cnt+2)+ham%pspot%kb_rsp(ico2+2)*y
         case(1)
         sumdum(isum_cnt)=sumdum(isum_cnt)+ham%pspot%kb_rsp(ico2)*y
         case default
          write(0, *) 'ERROR # of projectors is wrong', nref
         call mystop
        end select
         ico2=ico2+nref
      enddo  

      isum_cnt_st=isum_cnt_st + nref
     enddo

     call fast_all_sum_all_alloc_dc(nanl)
     call fast_all_sum_all_complex(sumdum(1), nanl) 

     ico1=0
     ico2=1
     isum_cnt_st=1
 
     do ia=1,natom_tot
       nref=ham%pspot%numref(ia)
       isum_cnt=isum_cnt_st
 
       select case(nref)
       case(8)
         sumdum(isum_cnt) = sumdum(isum_cnt) *ham%rsp_norm(isum_cnt)  
         sumdum(isum_cnt+1) = sumdum(isum_cnt+1) *ham%rsp_norm(isum_cnt+1)
         sumdum(isum_cnt+2) = sumdum(isum_cnt+2) *ham%rsp_norm(isum_cnt+2)
         sumdum(isum_cnt+3) = sumdum(isum_cnt+3) *ham%rsp_norm(isum_cnt+3)
         sumdum(isum_cnt+4) = sumdum(isum_cnt+4) *ham%rsp_norm(isum_cnt+4)
         sumdum(isum_cnt+5) = sumdum(isum_cnt+5) *ham%rsp_norm(isum_cnt+5)
         sumdum(isum_cnt+6) = sumdum(isum_cnt+6) *ham%rsp_norm(isum_cnt+6)
         sumdum(isum_cnt+7) = sumdum(isum_cnt+7) *ham%rsp_norm(isum_cnt+7)
       case(6)  
         sumdum(isum_cnt) = sumdum(isum_cnt) *ham%rsp_norm(isum_cnt)  
         sumdum(isum_cnt+1) = sumdum(isum_cnt+1) *ham%rsp_norm(isum_cnt+1)
         sumdum(isum_cnt+2) = sumdum(isum_cnt+2) *ham%rsp_norm(isum_cnt+2)
         sumdum(isum_cnt+3) = sumdum(isum_cnt+3) *ham%rsp_norm(isum_cnt+3)
         sumdum(isum_cnt+4) = sumdum(isum_cnt+4) *ham%rsp_norm(isum_cnt+4)
         sumdum(isum_cnt+5) = sumdum(isum_cnt+5) *ham%rsp_norm(isum_cnt+5)
       case(5)
         sumdum(isum_cnt) = sumdum(isum_cnt) *ham%rsp_norm(isum_cnt)  
         sumdum(isum_cnt+1) = sumdum(isum_cnt+1) *ham%rsp_norm(isum_cnt+1)
         sumdum(isum_cnt+2) = sumdum(isum_cnt+2) *ham%rsp_norm(isum_cnt+2)
         sumdum(isum_cnt+3) = sumdum(isum_cnt+3) *ham%rsp_norm(isum_cnt+3)
         sumdum(isum_cnt+4) = sumdum(isum_cnt+4) *ham%rsp_norm(isum_cnt+4)
       case(4)
         sumdum(isum_cnt) = sumdum(isum_cnt) *ham%rsp_norm(isum_cnt)  
         sumdum(isum_cnt+1) = sumdum(isum_cnt+1) *ham%rsp_norm(isum_cnt+1)
         sumdum(isum_cnt+2) = sumdum(isum_cnt+2) *ham%rsp_norm(isum_cnt+2)
         sumdum(isum_cnt+3) = sumdum(isum_cnt+3) *ham%rsp_norm(isum_cnt+3)
       case(3)
         sumdum(isum_cnt) = sumdum(isum_cnt) *ham%rsp_norm(isum_cnt)  
         sumdum(isum_cnt+1) = sumdum(isum_cnt+1) *ham%rsp_norm(isum_cnt+1)
         sumdum(isum_cnt+2) = sumdum(isum_cnt+2) *ham%rsp_norm(isum_cnt+2)
       case(1)
         sumdum(isum_cnt) = sumdum(isum_cnt) *ham%rsp_norm(isum_cnt)  
       case default
        write(0, *) 'ERROR # of projectors is wrong', nref
        call mystop
       end select

       do i=1,ham%pspot%nmap(ia)
         ico1=ico1+1
         s=zzero
         isum_cnt=isum_cnt_st

         select case(nref)

         case(8)
         s=s+sumdum(isum_cnt)*ham%pspot%kb_rsp(ico2) +&
             sumdum(isum_cnt+1)*ham%pspot%kb_rsp(ico2+1) +& 
             sumdum(isum_cnt+2)*ham%pspot%kb_rsp(ico2+2) +& 
             sumdum(isum_cnt+3)*ham%pspot%kb_rsp(ico2+3) +& 
             sumdum(isum_cnt+4)*ham%pspot%kb_rsp(ico2+4) +& 
             sumdum(isum_cnt+5)*ham%pspot%kb_rsp(ico2+5) +& 
             sumdum(isum_cnt+6)*ham%pspot%kb_rsp(ico2+6) +& 
             sumdum(isum_cnt+7)*ham%pspot%kb_rsp(ico2+7) 
         case(6)
         s=s+sumdum(isum_cnt)*ham%pspot%kb_rsp(ico2) +&
             sumdum(isum_cnt+1)*ham%pspot%kb_rsp(ico2+1) +& 
             sumdum(isum_cnt+2)*ham%pspot%kb_rsp(ico2+2) +& 
             sumdum(isum_cnt+3)*ham%pspot%kb_rsp(ico2+3) +& 
             sumdum(isum_cnt+4)*ham%pspot%kb_rsp(ico2+4) +& 
             sumdum(isum_cnt+5)*ham%pspot%kb_rsp(ico2+5) 
         case(5)
         s=s+sumdum(isum_cnt)*ham%pspot%kb_rsp(ico2) +&
             sumdum(isum_cnt+1)*ham%pspot%kb_rsp(ico2+1) +& 
             sumdum(isum_cnt+2)*ham%pspot%kb_rsp(ico2+2) +& 
             sumdum(isum_cnt+3)*ham%pspot%kb_rsp(ico2+3) +& 
             sumdum(isum_cnt+4)*ham%pspot%kb_rsp(ico2+4)  
         case(4)
         s=s+sumdum(isum_cnt)*ham%pspot%kb_rsp(ico2) +&
             sumdum(isum_cnt+1)*ham%pspot%kb_rsp(ico2+1) +& 
             sumdum(isum_cnt+2)*ham%pspot%kb_rsp(ico2+2) +& 
             sumdum(isum_cnt+3)*ham%pspot%kb_rsp(ico2+3)  
         case(3)
         s=s+sumdum(isum_cnt)*ham%pspot%kb_rsp(ico2) +&
             sumdum(isum_cnt+1)*ham%pspot%kb_rsp(ico2+1) +& 
             sumdum(isum_cnt+2)*ham%pspot%kb_rsp(ico2+2) 
         case(1)
         s=s+sumdum(isum_cnt)*ham%pspot%kb_rsp(ico2) 
         case default
          write(0, *) 'ERROR # of projectors is wrong', nref
          call mystop
         end select
         ico2=ico2+nref

         ham%fftstruc%rspacebuf(ham%pspot%indm(ico1)+iwfn2*rlen)=&
                ham%fftstruc%rspacebuf(ham%pspot%indm(ico1)+iwfn2*rlen)&
                     +  s*conjg(ham%cphase(ico1))
       enddo
     isum_cnt_st=isum_cnt_st + nref

     enddo

     end do

     if (iand(imonitor(1), 8) == 8) call get_timing(time(2),t3)

     call fourier_transform(1, ham%fftstruc, ham%gspace, hpsi(iwfn*len+1), &
             ham%fftstruc%rspacebuf(1),nfn ) 
     if (iand(imonitor(1), 8) == 8) call get_timing(time(1),t3)
   end do 

    deallocate(workr2_psi)

  else
  

  if (iand(imonitor(1), 1) == 1) t0 = gimmetime()  
 
  do i = 0, nwfn - 1, nbandsfft  
     call fft_convolute(ham%gspace, ham%fftstruc, hpsi(i * len + 1), &
          ham%vloc(1), psi(i * len + 1), min(nbandsfft, nwfn - i))
  end do
  if (iand(imonitor(1), 8) == 8) time(1)=time(1)+ gimmetime() - t3

   if (iand(imonitor(1), 1) == 1) then  
     t2 = gimmetime() - t0    
     dnx = mod(ham%gspace%fftsize(4) - ham%gspace%fftsize(5) + &
          ham%gspace%fftsize(1), ham%gspace%fftsize(1)) + 1
     nrodst1 = (dnx * ham%gspace%fftsize(3)) / ham%gspace%nproc  
     if (mod(dnx * ham%gspace%fftsize(3), ham%gspace%nproc) > &
          ham%gspace%myproc) nrodst1 = nrodst1 + 1
     nrodst2 = (ham%gspace%fftsize(2) * ham%gspace%fftsize(3)) / &
          ham%gspace%nproc
  
     if (mod(ham%gspace%fftsize(2) * ham%gspace%fftsize(3), &
          ham%gspace%nproc) > ham%gspace%myproc) nrodst2 = nrodst2 + 1
     !
     !        assume base 2 cmplx-cmplx: flops = 5n log_2 n-6n+6
     !
     n = real(ham%gspace%fftsize(3), dp)  
     opcount = real(ham%gspace%lorder, dp) * (dfive*n*log(n) / log(dtwo) - &
          dsix * n + dsix)
     n = real(ham%gspace%fftsize(2), dp)  
     opcount = opcount + real(nrodst1, dp) * (dfive*n*log(n) / log(dtwo) - &
          dsix * n + dsix)
     n = real(ham%gspace%fftsize(1), dp)  
     opcount = opcount + real(nrodst2, dp) * (dfive*n*log(n) / log(dtwo) - &
          dsix * n + dsix)
     ! fft-1 and convolution
     opcount = opcount * dtwo + real(ham%gspace%r_size, dp) * dsix  
     opcount = opcount * nwfn       ! for each wave function
     if (t2 > dzero) then  
        write (9, '(a,i4,a12,f12.6,a,f12.6)') 'TIME FOR ', nwfn, ' FFTS:', &
             & t2, ', MFLOPS=', (opcount / t2 * 1.0d-6)
        call myflush (9)  
     end if
  
   endif
  end if
  !
  !     kinetic energy term: multiply complex vector with real matrix
  !
  ekinmod = ham%ekinmod  
  do i = 0, nwfn - 1  
     do j = 1, len  
        ekinj = ham%gspace%ekin(j)  
        !
        !           modified kinetic energy expression
        !
        emod = ekinj - ham%shift + ekinmod(1) * (done + myderf((ekinj - &
             ekinmod(2)) / ekinmod(3)))
        hpsi(j + i * len) = hpsi(j + i * len) + psi(j + i * len) * emod
     end do
  end do
  !
  !     apply the nonlocal part of the hamiltonian
  !
  if( .not. ham%pspot%NL_rspace(1)) then

    if (iand(imonitor(1), 8) == 8) t3 = gimmetime() 
    if (iand(imonitor(1), 1) == 1) then  
      call parallel_barrier()  
      t0 = gimmetime()  
    end if
    if (nanl > 0) then                 ! do matrix-vector multiply
      if (nwfn == 1) then  
        call mzgemv('C', len, nanl, zone, ham%vnloc%data(1, 1, 1), &
             len, psi(1), 1, zzero, work(1), 1)
      else  
        call mzgemm('C', 'N', nanl, nwfn, len, zone, ham%vnloc%data(1, 1, 1), &
             len, psi(1), len, zzero, work(1), nanl)
      end if
    end if
    !
    !     scale with xnorm - This one scales poorly as it thrashes memory
    !

!    do i = 1, nanl  
!      call mzdscal(nwfn, ham%xnorm(i), work(i), nanl)
!    end do
 
    icnt=0
     do j=1,nwfn
     do  i = 1, nanl 
       icnt=icnt+1
       work(icnt)=work(icnt)* ham%xnorm(i)
     end do
    end do

    call fast_all_sum_all_alloc_dc(nanl * nwfn)
    call fast_all_sum_all_complex(work(1), nanl * nwfn)  
    !
    !     linearly combine the projectors to get H_nloc*psi
    !

    if (nanl > 0) then                 ! do matrix-vector multiply
      if (nwfn == 1) then  
        call mzgemv('N', len, nanl, zone, ham%vnloc%data(1, 1, 1), &
             len, work(1), 1, zone, hpsi(1), 1)
      else  
        call mzgemm('N', 'N', len, nwfn, nanl, zone, ham%vnloc%data(1, 1, 1), &
             len, work(1), nanl, zone, hpsi(1), len)
      end if
    end if

    if (iand(imonitor(1), 1) == 1) then  
      opcount = 1.6d1 * real(nanl * nwfn * len, dp)
      t2 = gimmetime() - t0  
      if (t2 > dzero) write (9, '(a,i4,a12,f12.6,a,f12.6)') 'TIME FOR ', &
          nanl, ' PROJECTORS:', t2, ', MFLOPS=', opcount / t2 * 1.0d-6
      call myflush(9)
      end if

    if (iand(imonitor(1), 8) == 8)  time(2)=time(2)+ gimmetime() - t3

  end if

  return  

end subroutine apply_ham_wavefn

