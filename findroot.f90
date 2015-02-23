!/*
!**                 determines the solution of
!**
!**                 a*x^3 + b*x^2 + c*x +d = 0
!**
!**       which has   3*a*x^2 + 2*b*x + c >0
!**
!**       and is the minimum of
!**        (1/4) a*x^4 + (1/3) b*x^3 +(1/2) c*x^2 + d*x
!**
!*/

function myasinh(x)

  use constants  
  implicit none

  real(dp) :: myasinh  
  real(dp), intent(in) :: x

  myasinh = log(x + sqrt(x * x + done))  

  return  

end function myasinh

function myacosh(x)

  use constants
  implicit none

  real(dp) :: myacosh
  real(dp), intent(in) :: x

  myacosh = log(x + sqrt(x * x - done))

  return  

end function myacosh

subroutine findroot(ret, ainp, binp, cinp, dinp)  

  use constants
  implicit none  

  real(dp), intent(out) :: ret       ! the return value
  real(dp), intent(in) :: ainp, &    ! a
       binp, &                       ! b
       cinp, &                       ! c
       dinp                          ! d
  !
  !     ------------- local variables ---------------------------
  !
  real(dp), parameter :: find_root_eps = 1.0d-40, &
       pi3rd = pi * dthird
  real(dp) :: myacosh, myasinh
  real(dp) :: asq, &
       acubed, bsq, bcubed, p, q, r, phi3rd, a, b, c, d, &
       x, x1, x2, x3, x1p, x2p, x3p, x1f, x2f, x3f
  integer, parameter :: imode = 0  
  !
  !     --------------------------------------------------------
  !
  !      write(9,*) 'a=',ainp
  !      write(9,*) 'b=',binp
  !      write(9,*) 'c=',cinp
  !      write(9,*) 'd=',dinp

  a = ainp ; b = binp ; c = cinp ; d = dinp

  asq = a * a ; acubed = asq * a ; bsq = b * b ; bcubed = bsq * b

  if (abs(a) <= find_root_eps) then      ! it is a polynomial of second order
     !
     !     this situation should rarely occur, and could be due to a
     !     failure of the conjugate gradient method. In case anything
     !     suspicious is detected, findroot returns a negative value to
     !     initiate  a reset of the conjugate gradient algorithm
     !
     p = c * c - dfour * b * d
     if (p < dzero) then  
        write(9, *) 'findroot: 2nd order polynomial has no real roots.'  
        write(9, *) '      resetting algorithm'  
        ret = -0.001d0  
        return  
     end if
     if (abs(b) < find_root_eps) then  
        write(9, *) 'findroot: encountered linear polynomial.'  
        write(9, *) '          resetting algorithm.'  
        ret = -0.0011d0  
        return  
     end if
     !        only the positive root is of interest, since the other one
     !        corresponds to a maximum
     !
     x = (-c + sqrt(p)) / (dtwo * b)  
     if (x > dfive) then  
        write(9, *) &
             'findroot: 2nd order polynomial shows runaway solution.'
        write(9, *) '          resetting algorithm.'  
        ret = -0.0012d0  
        return  
     end if
     ret = x  
     return  
  end if
  !     /* calculate p and q */
  !      q = (2.0*bcubed-9.0*a*b*c + 27.0*asq*d)  / (54.0 * acubed);
  !
  q = (b / a) * (b / a) * (b / a) / 27.0d0 - &
       (b / a) * (c / a) / dsix + d / (dtwo * a)

  !     /*  p = (3.0*a*c-bsq)/(9.0 *asq); */
  p = c / (a * dthree) - (b / a) * (b / a) / 9.0d0  
  !     /*  printf("p=%.5e, q=%.5e\n",p,q); */
  !     /* Find real roots depending on p, q */
  ! easy real solution   ??
  if (abs(p) < find_root_eps) then  
     if (q < dzero) then  
        x = (dtwo * abs(q))**dthird - b / (dthree * a)  
     else  
        x = -(dtwo * abs(q))**dthird - b / (dthree * a)  
     end if
     if ((dthree * x * x + dtwo * b * x + c) > dzero .and. x >= dzero) then  
        ret = x
        return  
     else  
        ret = -0.002d0
        return  
        !            write(9,100);   call mystop
     end if
     ! no, then have a closer look
  else  
     if (q >= dzero) then  
        r = sqrt(abs(p))  
     else  
        r = -sqrt(abs(p))  
     end if
     if (p > dzero) then                                   ! one real solution
        x = dmtwo * r * sinh(myasinh(q / (r * r * r) ) * dthird) - &
             b / (dthree * a)
        !            write(9,*) 3.d0*a*x*x + 2.d0*b*x + c
        !            write(9,*) x
        if ((dthree * a * x * x + dtwo * b * x + c) > dzero .and. &
             x >= dzero) then
           ret = x
           return  
        else  
           ret = -0.003d0
           return  
        end if
        
     else                                           ! look at the discriminant
        if ((q * q + p * p * p) <= dzero) then          ! three real solutions
           phi3rd = acos(q / (r * r * r)) * dthird 
           x1 = dmtwo * r * cos(phi3rd) - b / (dthree * a)
           x2 = dtwo * r * cos(pi3rd-phi3rd) - b / (dthree * a)
           x3 = dtwo * r * cos(pi3rd+phi3rd) - b / (dthree * a)
           x1p = dthree * a * x1 * x1 + dtwo * b * x1 + c
           x2p = dthree * a * x2 * x2 + dtwo * b * x2 + c
           x3p = dthree * a * x3 * x3 + dtwo * b * x3 + c
           x1f = dqtr * a * (x1 * x1 * x1 * x1) + &
                dthird * b * (x1 * x1 * x1) + dhalf * c * x1 * x1 + d * x1
           x2f = dqtr * a * (x2 * x2 * x2 * x2) + &
                dthird * b * (x2 * x2 * x2) + dhalf * c * x2 * x2 + d * x2
           x3f = dqtr * a * (x3 * x3 * x3) + &
                dthird * b * (x3 * x3 * x3) + dhalf * c * x3 * x3 + d * x3
           !              sort them according to size
           x = -1.0d20
           ! take the one that lowers the energy the most
           if (imode == 1) then  
              if (x3 > dzero .and. x3p > dzero) then  
                 if (x2 > dzero .and. x2p > dzero) then  
                    if (x1 > dzero .and. x1p > dzero) then  
                       if (x1f <= x2f .and. x1f <= x3f) x = x1
                       if (x2f <= x1f .and. x2f <= x3f) x = x2
                       if (x3f <= x1f .and. x3f <= x2f) x = x3
                    else                                    ! only x2, x3 left
                       if (x2f <= x3f) then  
                          x = x2
                       else  
                          x = x3
                       end if
                    end if
                 else                                       ! only x1, x3 left
                    if (x1 > dzero .and. x1p > dzero) then  
                       if (x1f <= x3f) then  
                          x = x1  
                       else  
                          x = x3 
                       end if
                    else  
                       x = x3
                    end if
                 end if
              else                                          ! only x1, x2 left
                 if (x2 > dzero .and. x2p > dzero) then  
                    if (x1 > dzero .and. x1p > dzero) then  
                       if (x1f <= x2f) then  
                          x = x1  
                       else  
                          x = x2  
                       end if
                    else  
                       x = x2
                    end if
                 else  
                    if (x1 > dzero .and. x1p > dzero) x = x1
                 end if
              end if
           else                  ! take smallest root which represents minimum
              if (x3 > dzero .and. x3p > dzero) then  
                 if (x2 > dzero .and. x2p > dzero) then  
                    if (x1 > dzero .and. x1p > dzero) then  
                       if (x1 <= x2 .and. x1 <= x3) x = x1
                       if (x2 <= x1 .and. x2 <= x3) x = x2
                       if (x3 <= x1 .and. x3 <= x2) x = x3
                    else                                   !  only x2, x3 left
                       if (x2 <= x3) then  
                          x = x2
                       else  
                          x = x3
                       end if
                    end if
                 else                                       ! only x1, x3 left
                    if (x1 > dzero .and. x1p > dzero) then  
                       if (x1 <= x3) then  
                          x = x1
                       else  
                          x = x3
                       end if
                    else  
                       x = x3
                    end if
                 end if
              else                                          ! only x1, x2 left
                 if (x2 > dzero .and. x2p > dzero) then  
                    if (x1 > dzero .and. x1p > dzero) then  
                       if (x1 <= x2) then  
                          x = x1
                       else  
                          x = x2
                       end if
                    else  
                       x = x2
                    end if
                 else  
                    if (x1 > dzero .and. x1p > dzero) x = x1
                 end if
              end if
           end if
           if (x > dzero) then  
              ret = x
              return
           else  
              ret = -0.004d0
              return  
              !                  write(9,100); call mystop
           end if
        else                                              !  one real solution
           x = dmtwo * r * cosh(myacosh(q / (r * r * r)) * dthird) - &
                b / (dthree * a)
           if ((dthree * a * x * x + dtwo * b * x + c) > dzero .and. &
                x >= dzero) then
              ret = x
              return
           else  
              ret = -0.005d0
              return  
              !                  write(9,100); call mystop
           end if
        end if                                      
     end if                                         ! look at the discriminant
  end if                                             ! then have a closer look

  ret = -1.0d20  

  return  

100 format('*** MATRIX DIAGONALIZATION BOMBOUT!', &
       &     /'-->  CHECK INPUT CHARGE DENSITY,' &
       &     /'-->  MAYBE REDUCE MIXING,' &
       &     /'-->  INCREASE DIAGSAFETY')

end subroutine findroot
