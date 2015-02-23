Module bessfn_mod

  Use constants

  Implicit none

Contains

  !========================================================
  !=  Function returns value of spherical bessel function =
  !=                                                      =
  !=               Chris J. Pickard Jan. 1999             =
  !========================================================

  Function bessfn(x,l)

    Integer,  Intent(in) :: l
    Real(dp), Intent(in) :: x

    Real(dp) :: bessfn

    ! Local

    Integer             :: i,m,isum
    Integer,  Parameter :: IACC=40
    Real(dp)            :: j,jm,jp,oox,summ,sn,cs
    Real(dp), Parameter :: xsplit = 0.01d0

    If(x<10.d0*Tiny(x)) Then
       If(l.Eq.0) Then
          bessfn=1.d0
       Else
          bessfn=0.d0
       Endif
    Else If(x>xsplit*l) Then
       oox = 1.d0/x
       jm=Sin(x)*oox
       j=(jm-Cos(x))*oox
       Do i=1,l-1
          jp=(2*i+1)*oox*j-jm
          jm=j
          j=jp
       End Do
       If(l.Eq.0) Then
          bessfn=jm
       Else
          bessfn=j
       End If
    Else
       oox = 1.d0/x
       m=2*((l+Int(Sqrt(Dble(IACC*l))))/2)
       bessfn=0.d0
       isum=0
       summ=0.d0
       jp=0.d0
       j=1.d0
       Do i=m,1,-1
          jm=(2*i+1)*oox*j-jp
          jp=j
          j=jm
          If(isum/=0) summ=summ+j
          isum=1-isum
          If(i==l) bessfn=jp
       Enddo
       summ=2.d0*summ-j
       bessfn=bessfn/summ
    Endif

  End Function bessfn

End Module bessfn_mod
