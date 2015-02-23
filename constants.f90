! Module containing constants for Paratec code
!
! Peter Haynes March 2000

module constants

  implicit none

  integer, parameter :: dp = kind(1.0d0),  num_time=14 !itmaxpul=30
  integer, parameter :: llmax= ( (2+1)**2+1) ! maximum lm quantum number

  real(dp), parameter :: dzero = 0.0d0, &
       done = 1.0d0, &
       dmone = -done, &
       dtwo = 2.0d0, &
       dmtwo = -dtwo, &
       dthree = 3.0d0, &
       dfour = 4.0d0, &
       dmfour = -dfour, &
       dfive = 5.0d0, &
       dsix = 6.0d0

  real(dp), parameter :: dhalf = 0.5d0, &
       dmhalf = -dhalf, &
       dtrhalf = 1.5d0, &
       dthird = 1.0d0 / 3.0d0, &
       dmthird = -dthird, &
       dttrd = 2.0d0 / 3.0d0, &
       dmttrd = -dttrd, &
       dftrd = 4.0d0 / 3.0d0, &
       dqtr = 0.25d0, &
       dsixth = 1.0d0 / 6.0d0

  real(dp), parameter :: drt2 = 1.41421356237309504880168872421d0, &
       drt3 = 1.73205080756887729352744634151d0, &
       drt5 = 2.23606797749978969640917366873d0, &
       drt15 = drt3 * drt5, &
       oort2 = done / drt2, &
       oort3 = done / drt3

  real(dp), parameter :: pi = 3.14159265358979323846264338328d0, &
       pi2 = dtwo*pi, &
       pi4 = dfour*pi, &
       pi43 = pi4 / dthree, &
       ootwopi = done / pi2

  real(dp), parameter :: ryd = 13.60580d0

  real(dp), parameter :: &
       B_to_M = 0.529177d-10, &  
       M_to_KG = 1.6605402d-27, &
       eV_to_J = 1.60217733d-19, &
       Kb_eV = 8.6217d-5, &
       Kb_Ryd = Kb_eV/ryd

  complex(dp), parameter :: zzero = (0.0d0,0.0d0), &
       zone = (1.0d0,0.0d0), &
       zmone = -zone, &
       ztwo = (2.0d0,0.0d0), &
       zmtwo = -ztwo, &
       zthree = (3.0d0,0.0d0), &
       zmthree = -zthree, &
       znine = (9.0d0,0.0d0), &
       zhalf = (0.5d0,0.0d0), &
       zmhalf = -zhalf, &
       ztrhalf = (1.5d0,0.0d0)

  complex(dp), parameter :: zi = (0.0d0,1.0d0)

end module constants
