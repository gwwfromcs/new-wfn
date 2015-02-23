!*
subroutine warn(i, xarg, iarg, rout)
  !
  !      HANDLES THE WARNINGS
  !      AND KEEPS TRACK OF THEIR NUMBER
  !
  use constants
  implicit none
  !
  integer, intent(in) :: i, iarg
  real(dp), intent(in) :: xarg
  character(len=6), intent(in) :: rout
  integer, save :: nwarn = 0
  integer :: i10, ir
  !
  if (i == 0) then  
     if (nwarn >= 1) write(9, 1001) nwarn
  else  
     nwarn = nwarn + 1  
     i10 = i / 10  
     ir = i - i10 * 10  
     write(9, 1000) i, rout
     !
     if (i10 == 13) then  
        if (ir == 2) then  
           write(9, 1132)  
        else if (ir == 3) then  
           write(9, 1133)  
        else if (ir == 4) then  
           write(9, 1134) xarg
        else if (ir == 5) then  
           write(9, 1135) xarg  
        end if
     else if (i10 == 17) then  
        if (ir == 0) then  
           write(9, 1170)  
        end if
     else if(i10 == 22) then  
        if (ir == 0) then  
           write(9, 1220) iarg 
        else if(ir == 1) then  
           write(9, 1221) iarg 
        end if
     else if (i10 == 25) then  
        if (ir == 0) then  
           write(9, 1250) xarg  
        else if (ir == 2) then  
           write(9, 1252) iarg, xarg  
        end if
     else if (i10 == 27) then  
        if (ir == 0) then  
           write(9, 1270) iarg  
        else if (ir == 1) then  
           write(9, 1271) iarg, xarg  
        else if (ir == 2) then  
           write(9, 1272) xarg  
        end if
     else if (i10 == 28) then  
        if (ir == 0) then  
           write(9, 1280)  
        end if
     else if (i10 == 29) then  
        if (ir == 0) then  
           write(9, 1290) xarg  
        else if (ir == 1) then  
           write(9, 1291) xarg  
        end if
     else if (i10 == 30) then  
        if (ir == 0) then  
           write(9, 1300) iarg  
        end if
     else if (i10 == 34) then  
        if (ir == 0) then  
           write(9, 1340)  
        end if
     else if (i10 == 36) then  
        if (ir == 0) then  
           write(9, 1360)  
        end if
     else if (i10 == 47) then  
        if (ir == 0) then  
           write(9, 1470) xarg  
        end if
     end if
     !
  end if

  return
  
1000 FORMAT ('  *** WARNING *** ',I6,' IN ',A6)  
1001 FORMAT (/,2X,I3,' WARNINGS UP TO NOW')  
1132 FORMAT ('  CHEMICAL SYMBOLS DO NOT MATCH')  
1133 FORMAT ('  CORRELATION POTENTIAL DOES NOT MATCH')  
1134 FORMAT ('  TAPE CONTAINS INFORMATION ABOUT', &
       &        ' LOCAL POTENTIAL UP TO GMAX = ',E10.3)
1135 FORMAT ('  TAPE CONTAINS INFORMATION ABOUT', &
       &        ' NON-LOCAL POTENTIAL UP TO GMAX = ',E10.3)
1170 FORMAT ('  ERROR READING SCREENING TAPE ', &
       &        'DEFAULT SCREENING WAS USED')
1220 FORMAT ('  ERROR ',I5,' DIAGONALIZING HAMILTONIAN ')  
1221 FORMAT ('  ERROR ',I5,' CALCULATING EIGENVECTORS')  
1250 FORMAT ('  NON HERMITIAN DIIS MATRIX ',E10.3)  
1252 FORMAT ('  INACCURATE DIIS EIGENVECTOR ',I4,E10.3)  
1270 FORMAT ('  N(FFT) REDUCED IERR = ',I4)  
1271 FORMAT ('  COMPLEX CHARGE DENSITY IN', &
       &        ' REAL SPACE - IERR = ',I5,' MAXERR = ',E10.3)
1272 FORMAT ('  NEGATIVE CHARGE DENSITY IN', &
       &        ' REAL SPACE, RHO-MIN = ',E10.3)
1280 FORMAT ('  UNKNOWN CORRELATION, EXC ASSUMED ZERO')  
1290 FORMAT ('  WRONG TOTAL CHARGE DENSITY    DIFF = ',E12.6)  
1291 FORMAT ('  COMPLEX TOTAL CHARGE DENSITY DENI(1) = ',E10.3)  
1300 FORMAT ('  TRIVIAL UPDATE MATRIX, MXDUPD = ',I5)  
1340 FORMAT ('  PROBLEMS IN BASIS SET REPLACEMENT')  
1360 FORMAT ('  SOME SUBROUTINE CALLS IN MAIN PROGRAM', &
       &        ' MAY BE INCOMPATIBLE WITH LIBRARY')
1470 FORMAT ('  FLAT BAND AT E = ',E10.3)  

end subroutine warn
!*
!*
subroutine fatal(i, xarg, iarg, rout)  
  !
  !      HANDLES THE FATAL ERRORS
  !      AND STOPS THE EXECUTION OF THE PROGRAM
  !
  use constants
  implicit none
  !
  integer, intent(in) :: i, iarg
  real(dp), intent(in) :: xarg
  character(len=6), intent(in) :: rout  
  integer :: i10, ir
  !
  write(9, 1000) rout
  i10 = i / 10  
  ir = i - i10 * 10  
  if (i10 == 5) then  
     !        CRSTL
     if (ir == 0) then  
        write(9, 1050) xarg  
     else if (ir == 1) then  
        write(9, 1051) iarg  
     else if (ir == 2) then  
        write(9, 1052) iarg  
     end if
  else if (i10 == 9) then  
     if (ir == 0) then  
        write(9, 1090) xarg, iarg  
     else if (ir == 1) then  
        write(9, 1091) xarg, iarg  
     else if (ir == 2) then  
        write(9, 1092) xarg, iarg  
     else if (ir == 3) then  
        write(9, 1093) iarg
     end if
  else if (i10 == 11) then  
     if (ir == 0) then  
        write(9, 1110) xarg  
     end if
  else if (i10 == 13) then  
     if (ir == 0) then  
        write(9, 1130) iarg  
     else if (ir == 1) then  
        write(9, 1131) iarg  
     end if
  else if (i10 == 19) then  
     if (ir == 0) then  
        write(9, 1190) iarg  
     else if (ir == 1) then  
        write(9, 1191) iarg  
     else if (ir == 2) then  
        write(9, 1192) iarg  
     end if
  else if (i10 == 21) then  
     if (ir == 0) then  
        write(9, 1210) iarg  
     end if
  else if (i10 == 22) then  
     if (ir == 0) then  
        write(9, 1220) iarg  
     end if
  else if (i10 == 23) then  
     if (ir == 0) then  
        write(9, 1230)  
     end if
  else if (i10 == 25) then  
     if (ir == 1) then  
        write(9, 1251) iarg, xarg  
     end if
  else if (i10 == 27) then  
     if (ir == 0) then  
        write(9, 1270) iarg  
     end if
  else if (i10 == 33) then  
     if (ir == 0) then  
        write(9, 1330) iarg  
     end if
  else if (i10 == 34) then  
     if (ir == 0) then  
        write(9, 1340) iarg  
     else if (ir == 1) then  
        write(9, 1341) iarg  
     else if (ir == 2) then  
        write(9, 1342) iarg  
     else if (ir == 3) then  
        write(9, 1343) iarg  
     end if
  else if (i10 == 35) then  
     if (ir == 0) then  
        write(9, 1350) iarg  
     else if (ir == 1) then  
        write(9, 1351) iarg  
     end if
  else if (i10 == 36) then  
     if (ir == 0) then  
        write(9, 1360)  
     end if
  else if (i10 == 37) then  
     if (ir == 0) then  
        write(9, 1370)  
     end if
  end if
  call mystop

1000 FORMAT ('  ***FATAL ERROR*** IN ',A6)  
1050 FORMAT ('  CELL VOLUME = ',E12.5)  
1051 FORMAT ('  TYPES OF ATOMS = ',I5)  
1052 FORMAT ('  NUMBER OF ATOMS = ',I5)  
1090 FORMAT ('  GMOD = ',E10.3,'  TRY  NO. OF G-VEC.= ',I6)  
1091 FORMAT ('  GMAX = ',E10.3,'  NO. OF G-VEC.= ',I6)  
1092 FORMAT ('  GMAX**2 = ',E10.3, &
       &        '  FOR STAR NO.= ',I4)
1093 FORMAT ('  CHANGE MXDCUB TO AT LEAST  ',I5)  
1110 FORMAT ('  PHASE AND STRUCTURE FACTORS', &
       &        ' DO NOT AGREE, ERROR = ',E10.3)
1130 FORMAT ('  DIMENSIONS FOR LOCAL PSEUDOPO', &
       &        'TENTIAL ARRAY TOO SMALL   MXD = ',I5)
1131 FORMAT ('  DIMENSIONS FOR NON-LOCAL PSEUDOPOTENTIAL', &
       &        ' ARRAY TOO SMALL   MXD = ',I5)
1190 FORMAT ('  CHANGE MXDBND TO AT LEAST  ',I5)  
1191 FORMAT ('  DIMENSIONS FOR UNSYMMETRIZED INTEGRATION', &
       &        ' K-POINTS TOO SMALL   N=',I6)
1192 FORMAT ('  DIMENSIONS FOR NUMBER OF INTEGRATION K-POINTS', &
       &        'TOO SMALL  N=',I5)
1210 FORMAT ('  CHANGE MXDDIM TO AT LEAST  ',I5)  
1220 FORMAT ('  CHANGE MXDANL TO AT LEAST  ',I5)  
1230 FORMAT ('  UNRECOGNIZED MATRIX STORAGE')  
1251 FORMAT ('  WRONG DIIS EIGENVECTOR ',I4,E10.3)  
1270 FORMAT ('  SIZE OF FFT ARRAY TOO SMALL  NFFT=',I7)  
1330 FORMAT ('  UNABLE TO FIND EQUIVALENT ATOM FOR K= ',I5)  
1340 FORMAT ('  CHANGE MXDSML TO AT LEAST  ',I5)  
1341 FORMAT ('  INCREASE MXDSML MAYBE TO A VALUE OF  ',I5)  
1342 FORMAT ('  INCREASE SMALL MATRIX SIZE, MTXDS, TO  ',I5)  
1343 FORMAT ('  PROBLEMS WITH GRAM-SCHMIDT, NSIZE = ',I5)  
1350 FORMAT ('  CHANGE MXDFFT TO AT LEAST  ',I5, &
       &        '  OR USE ANOTHER WORK ARRAY')
1351 FORMAT ('  THE FFT MESH IS NOT SUFFICIENTLY LARGE', &
       &        '  TRY MXDFFT TWICE OF ',I7,/, &
       &        '  OR INCREASE KMSCR (KMAX)')
1360 FORMAT ('  INCOMPATIBLE VERSION OF DRIVER PROGRAM')  
1370 FORMAT ('  INVALID CHOICE OF PSEUDOPOTENTIAL')  

end subroutine fatal
!
!
!
subroutine alccheck(wst, isize, istat)

  implicit none  

  integer, intent(in) :: isize, istat  
  character(len=*), intent(in) :: wst

  if (istat /= 0) then  
     write(9, '(a,a)') ' *** MEMORY ALLOCATION FAILED ON ARRAY ', wst
     write(9, '(a,i12)') '     ARRAY SIZE: ', isize  
     call myflush(9)  
     call mystop  
     !      else
     !         write(9,'(a,a)')
     !     $        ' *** MEMORY ALLOCATION SUCCEEDED ON ARRAY ',wst
     !         write(9,'(a,i12)') '     ARRAY SIZE: ',isize
     !         call myflush(9)
  end if

end subroutine alccheck
