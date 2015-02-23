!     @process extchk
!
subroutine neighbor(crys, pw_params)  
  !
  !     find the nearest neighbor distances and other information
  !
  !
  !     1996 by Bernd Pfrommer, while at UCB
  include 'use.h'  
  implicit none             ! implicit? Just say no!
  include 'interface.h'  
  !
  !     INPUT:
  !     -----
  !
  type(crystal), intent(in) :: crys  
  type(pw_parameter), intent(in) :: pw_params  
  !
  !     ---------- local variables ---------------------
  !
  real(dp) :: rc(3), rc1(3), drc(3), d, d2, del, dr(3), &
       dum, alld(2, crys%mxdatm * crys%ntype), angl(crys%mxdatm * crys%ntype), &
       alldr(3, crys%mxdatm * crys%ntype, crys%mxdatm * crys%ntype),&
       dnn(crys%mxdatm * crys%ntype)
  integer :: indx(crys%mxdatm * crys%ntype), nt, ia, nt1, ia1, i, &
       j, k, i1, i2, i3, iv (3), nn, nn0, nn1, inn, i_t
  character (len=255) :: fmtstr  
  character (len=5),allocatable :: label(:)
  if (iand(pw_params%output(1), 64) /= 64 .and. &
       iand(pw_params%output(1), 8192) /= 8192) return
  !
  !     count number of atoms
  !
  nn0 = 0  
  do nt = 1, crys%ntype  
     do ia = 1, crys%natom(nt)  
        nn0 = nn0 + 1  
     end do
  end do

  allocate (label(nn0))

  if (iand(pw_params%output(1), 64) == 64) then  
     write(9, 10)  
!     write(9, 20) ((crys%nameat(nt), j, j = 1, crys%natom(nt)), &
!          nt = 1, crys%ntype)
     write(9, 21)
!     write(fmtstr, 24) nn0  
!     write(9, fmtstr) (j, j = 1, nn0)  
  end if

  alldr = 0  
  nn0 = 0  
  do nt = 1, crys%ntype  
     do ia = 1, crys%natom(nt)  
        nn0 = nn0 + 1  
        write (label(nn0),1000)crys%nameat(nt),ia
1000 Format (A2,I3) 
!        write (9,*) label(nn0)
        rc = crys%rat(:, ia, nt) / pi2  
        !
        !           make sure it is in first unit cell
        !
        rc(:) = mod(rc(:) + done, done)  
        rc = matmul(crys%avec, rc)  
        !
        !           now check for neighbors
        !
        nn = 0  
        nn1 = 0  
        do nt1 = 1, crys%ntype  
           do ia1 = 1, crys%natom(nt1)  
              nn1 = nn1 + 1  
              nn = nn + 1  
              d = 1.0d6  
              d2 = 1.0d6  
              do i1 = -2, 2
                 do i2 = -2, 2
                    do i3 = -2, 2  
                       iv(1) = i1
                       iv(2) = i2
                       iv(3) = i3  
                       rc1 = crys%rat(:, ia1, nt1) / pi2  
                       rc1(:) = mod(rc1(:) + done, done)  
                       rc1 = matmul(crys%avec, real(iv, dp) + rc1)  
                       !     write(9,*) 'rc test=',rc1
                       drc = rc1 - rc  
                       !     write(9,*) 'deltarc=',drc
                       del = sqrt(drc(1) * drc(1) + drc(2) * drc(2) + &
                            drc(3) * drc(3))
                       if (del < d - 1.0d-6 .and. del > 1.0d-6) then  
                          d2 = d  
                          ! found nearest neighbor
                          d = del  
                          dr = drc  
                       else if (del < d2 .and. del > 1.0d-6 .and. &
                            del > d + 1.0d-6) then
                          ! found 2nd nearest neighbor
                          d2 = del  
                       end if
                    end do
                 end do
              end do  
              alld(1, nn) = d  
              alld(2, nn) = d2  
              alldr(:, nn, nn0) = dr  
           end do
        end do
        if (iand(pw_params%output(1), 64) == 64) then  
!           write(9, 100) crys%nameat(nt), ia, crys%rat(:, ia, nt) / pi2, &
!                (alld(1, i), i = 1, nn)
!           write(9, 110) nn0, (alld(2, i), i = 1, nn)  
           !
           !              compute ordered list of neighbors
           !
           do i = 1, nn  
              indx(i) = i  
           end do
           call sort(nn, alld(1, :), indx)  
           dnn(1)=alld(1, indx(1))

           inn=1
           do i=2,nn
             if (alld(1, indx(i)) .lt. 1.2*dnn(1)) then
               inn=inn+1
             else
               exit
             end if
           end do
           i_t=i 

           write(9, 115) crys%nameat(nt), nn0, crys%rat(:, ia, nt) / pi2,&
                            (indx(i), alld(1, indx(i)), i = 1, inn)  
!           call sort(nn, alld(2, :), indx) 

           dnn(1)=alld(1, indx(i_t))

!           inn=1
           do i=i_t,nn
             if (alld(1, indx(i)) .lt. 1.1*dnn(1)) then
               inn=inn+1
!               dnn(inn)=alld(2, indx(i))
             end if
           end do

          write(9, 120) (indx(i), alld(1, indx(i)), i = i_t, inn) 
 !         write(9, 120) (dnn(i), i = 1, inn) 

        end if
     end do
  end do
  !
  !     -------------- compute and print the angles ------------
  !
  if (iand(pw_params%output(1), 8192) == 8192) then  

     write(9, 12)  
     do i = 1, nn0  
        write(9, 30) label(i)  
        write(9,45)
        do j = 1, nn0  
           do k = 1, j-1  
              dr = alldr(:, j, i)  
              del = sqrt(dot_product(dr, dr))  
              drc = alldr(:, k, i)  
              d = sqrt(dot_product(drc, drc))  
              if (del > 1.0d-6 .and. d > 1.0d-6) then  
                 dum = dot_product(dr, drc) / del / d  
                 if (dum > done) dum = done
                 if (dum.lt. - 1.d0) dum = - 1.d0  
                 angl(k) = 180.d0 / pi * acos(dum)  
                 if (del < pw_params%bd_lgth .and. d < pw_params%bd_lgth) then
                    write(9,50) label(j),del,label(k),d,angl(k)
                 endif
              else  
                 angl(k) = dzero  
              end if
           end do
!           write(9, 40) j, (angl(k), k = 1, j)  
        end do
     end do
  end if

  call myflush(9)
  deallocate (label)
10 format(/1x,20('-'),' NEIGHBOR INFORMATION ',20('-'))  
12 format(/1x,20('-'),'  ANGLE   INFORMATION ',20('-'))  
20 format(/2x,'ATOM',12x,'COORD',9x,1000(5x,a4,i3))  
21 format(/2x,'ATOM',8x,'COORD',8x,'neighbor list')  
22 format(2x,'(NR)',12x,'     ',9x,1000(6x,'(',i3,')',1x)/)  
24 format('(2x,''(NR)'',12x,''     '',9x,', &
       &     i5,'(6x,''('',i3,'')'',1x)/)')
30 format(/' %%%%%%% BOND CENTER: ATOM ',A5,' %%%%%%%%%')  
40 format(i3,1000f7.2)
45 format(1x,' Atom ',' distance ',' Atom ',' distance ','  Angle  ')
50 format(1x,2(A5,f10.5,1x),f7.2) 
100 format(1x,a2,i3,3f9.4,'|', 500f12.6)  
110 format(1x,'(',i3,')',27x,'|', 500f12.6)  
115 format(1x,a2,i3,3f6.2,'|', 500(i4,'(',f4.2,')'))  
120 format(24x,'|', 500(i4,'(',f4.2,')'))   

end subroutine neighbor
