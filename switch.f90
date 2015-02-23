!
function switch(x, c)  
  !
  !  Converts input values according to the control parameter c
  !
  !  Last modification: MOD_DATE
  !
  use constants
  implicit none
  !
  real(dp) :: switch
  !
  !     .. Scalar Arguments ..
  real(dp), intent(inout) :: x  
  character(len=1), intent(in) :: c
  !     ..
  !     .. Local Scalars ..
  real(dp) :: sign  
  !     ..
  !     .. External Functions ..
  logical, external :: leqi  
  !     ..
  sign = done  
  if (x < dzero) sign = -sign  
  x = abs(x)  
  if (leqi(c, 's')) x = sqrt(x)  
  if (leqi(c, 'c')) x = x**dthird
  if (leqi(c, 't')) x = x * dthird
  if (leqi(c, 'd')) x = x * dhalf
  if (leqi(c, 'r')) x = x / drt2
  if (leqi(c, 'h')) x = x / drt3
  if (leqi(c, 'i')) x = done / x  
  if (leqi(c, 'a')) x = sqrt(x * dthird)  
  if (leqi(c, 'b')) x = x / (dsix * 7.5d0)  
  !
  switch = sign * x  
  !
  return  
  !
end function switch
