!  
!  E l e c t r o n i c   S t r u c t u r e   D a t a   F o r m a t
!  ---------------------------------------------------------------
!  
!                             E S D F
!                             =======
!  
!  Author: Chris J. Pickard (c)
!  Email : cp@min.uni-kiel.de
!  Place : Kiel, Germany
!  Date  : 5/6th August 1999
!  
!  Summary
!  -------
!  
!  This module is designed to simplify and enhance the input of data into
!  electronic structure codes (for example, CASTEP). It works from a
!  highly flexible input file. Data is input in a "label <value>"
!  fashion, and input is independent of the ordering of the input
!  file. An important feature is the requirement that most inputs require
!  default settings to be supplied within the main program calling
!  ESDF. This means that rarely used variables will not clutter everyday
!  input files, and, even more usefully, "intelligence" may be built into
!  the main code as the defaults may be dependent of other set
!  variables. Block data may also be read in. Another important feature
!  is the ability to define "physical" values. This means that the input
!  files need not depend on the internal physical units used by the main
!  program.
!  
!  
!  History
!  -------
!  
!  ESDF has been written from scratch in F90, but is heavily based
!  (especially for the concept) on the FDF package developed by Alberto
!  Garcia and Jose Soler. It is not as "flexible" as FDF - there is no
!  provision for branching to further input files. This simplifies the
!  code, and I hope that it is still useful without this feature. Also,
!  the input and defaults are not dumped to a output file currently. I've
!  not found this a hindrance as of now.
!  
!  
!  Future
!  ------ 
!  
!  My intention is to make this release available to Alberto Garcia and
!  Jose Soler for their comments. It might be a good idea to use this as
!  a base for fully converting the FDF package to F90. Or it may remain
!  as a cut down version of FDF. I certainly hope that a package of the
!  FDF sort becomes widely used in the electronic structure community. My
!  experience has been very positive.
!  
!  
!  Usage
!  -----
!  
!  First, "Use esdf" wherever you wish to make use of its features. In
!  the main program call the initialisation routine: Call
!  esdf_init('input.esdf'). "input.esdf" is the name of the input file -
!  it could be anything. This routine opens the input file, and reads
!  into a dynamically allocated storage array. The comments and blank
!  lines are stripped out. You are now ready to use the
!  esdf_functions. For example, if you want to read in the number of
!  atoms in your calculation, you would use: natom =
!  esdf_integer('NumberOfAtoms',1), where 'NumberOfAtoms' is the label to
!  search for in the input file, and '1' is the default. Call esdf_close to
!  deallocate the data arrays. You may then open another input file using
!  esdf_init. It is not currently possible to open more that on input 
!  file at one time.
!  
!  
!  Syntax
!  ------
!  
!  The input file can contain comments. These are defined as anything to
!  the right of, and including, '#', ';', or '!'. It is straightforward
!  to modify the routine to accept further characters. Blank lines are
!  ignored -- use comments and blank lines to make you input file
!  readable.
!  
!  The "labels" are case insensitive (e.g. unitCell is equivalent to
!  UnItceLL) and punctuation insensitive (unit.cell is equivalent to
!  unit_cell is equivalent to unitcell). Punctuation characters are '.'
!  and '-' at the moment. Again - use this feature to improve
!  readability.
!  
!  The following are equivalent ways of defining a physical quantity:
!  
!  "AgeOfUniverse = 24.d0 s" or "AgeOfUniverse : 24.d0 S" or
!  "AgeOfUniverse 24.d0 S"
!   
!  It would be read in by the main program in the following way:
!  
!  aou = esdf_physical('ageofuniverse',77.d0,'ns')
!  
!  "aou" is the double precision variable, 77.d0 is the default number of
!  "ns" or nanoseconds. 24s will be converted automatically to its
!  equivalent number of nanoseconds.
!  
!  Block data should be placed in the input file as follows:
!  
!  begin cellvectors 
!  1.0 1.0 0.0 
!  0.0 1.0 1.0 
!  1.0 0.0 1.0 
!  end cellvectors
!  
!  And it may be read:
!  
!    if(esdf_block('CellVectors',nlines))
!      if(nlines.ne.3) then (... break out here if the incorrect number
!  of lines)
!      do i=1,nlines
!        read(block_data(i),*) x,y,z
!      end do
!    endif
!  
!  
!  List of functions
!  -----------------
!  
!  Self explanatory:
!  
!  esdf_string(label,default)
!  esdf_integer(label,default)
!  esdf_single(label,default)
!  esdf_double(label,default)
!  esdf_physical(label,default,unit)
!  
!  A little more explanation:
!  
!  esdf_defined(label) is true if "label" found, false otherwise
!  
!  esdf_boolean(label,default) is true if "label yes/true/t (case/punct.insens)
!                              is false if"label no/false/f (case/punct.insens)
!  
!  The Help feature
!  ----------------
!  
!  The routine "esdf_help(helpword,searchword)" can be used to access the
!  information contained within the "esdf_key_mod" module.
!  
!  If "helpword" is "search" (case insensitive), then all variables whose
!  description contains "searchword" will be output.
!  
!  If "helpword" is "basic", "inter", "expert" or "dummy" the varibles of
!  that type will be displayed.
!  
!  If "helpword" is one of the valid labels, then a description of this
!  label will be output.
!  
!  
!  Finishing off
!  -------------
!  
!  Two routines, "esdf_warnout" and "esdf_close", can be used to finish
!  the use of ESDF. "esdf_warnout" outputs ESDF warnings to screen, and
!  "esdf_close" deallocates the allocated ESDF arrays.
!  
!  Contact the Author
!  ------------------
!  
!  This code is under development, and the author would be very happy to
!  receive comments by email. Any use in a commercial software package is 
!  forbidden without prior arrangement with the author (Chris J. Pickard).

Module esdf

  Use esdf_key

  Implicit None

  ! Kind parameters 

  Integer, Private, Parameter :: I4B = Kind(1)
  Integer, Private, Parameter :: DP  = Kind(1.d0)
  Integer, Private, Parameter :: SP  = Kind(1.0)

  ! Set the length of the lines

  Integer(I4B), Public, Parameter :: llength=80
  Integer(I4B), Private, Parameter ::  nphys = 54
  Integer(I4B), Private :: nrecords,nwarns
  Character(llength), Private, Dimension(:), Allocatable :: llist,warns
  Character(llength), Private, Dimension(:,:), Allocatable :: tlist 

  ! The public block data array

  Character(llength), Public, Dimension(:), Allocatable :: block_data 

  ! Set the physical units database
  
  Character(10), Private :: phy_d(nphys),phy_n(nphys) ! d - dimension n - name
  Real(DP), Private      :: phy_u(nphys)              ! u - unit

  !
  !     We allow case variations in the units. This could be dangerous
  !     (meV --> MeV !!) in real life, but not in this restricted field.
  !     
  ! m - mass l - length t - time e - energy f - force p - pressure c - charge
  ! d - dipole mom - mom inert ef - efield
  !

  Data phy_d(1) /'m'/;Data phy_n(1) /'kg'/;Data phy_u(1) /1.0_dp/
  Data phy_d(2) /'m'/;Data phy_n(2) /'g'/;Data phy_u(2) /1.e-3_dp/
  Data phy_d(3) /'m'/;Data phy_n(3) /'amu'/;Data phy_u(3) /1.66054e-27_dp/
  Data phy_d(4) /'l'/;Data phy_n(4) /'m'/;Data phy_u(4) /1.0_dp/
  Data phy_d(5) /'l'/;Data phy_n(5) /'nm'/;Data phy_u(5) /1.e-9_dp/
  Data phy_d(6) /'l'/;Data phy_n(6) /'ang'/;Data phy_u(6) /1.e-10_dp/
  Data phy_d(7) /'l'/;Data phy_n(7) /'bohr'/;Data phy_u(7) /0.529177e-10_dp/
  Data phy_d(8) /'t'/;Data phy_n(8) /'s'/;Data phy_u(8) /1.0_dp/
  Data phy_d(9) /'t'/;Data phy_n(9) /'ns'/;Data phy_u(9) /1.e-9_dp/
  Data phy_d(10) /'t'/;Data phy_n(10) /'ps'/;Data phy_u(10) /1.e-12_dp/
  Data phy_d(11) /'t'/;Data phy_n(11) /'fs'/;Data phy_u(11) /1.e-15_dp/
  Data phy_d(12) /'e'/;Data phy_n(12) /'j'/;Data phy_u(12) /1.0_dp/
  Data phy_d(13) /'e'/;Data phy_n(13) /'erg'/;Data phy_u(13) /1.e-7_dp/
  Data phy_d(14) /'e'/;Data phy_n(14) /'ev'/;Data phy_u(14) /1.60219e-19_dp/
  Data phy_d(15) /'e'/;Data phy_n(15) /'mev'/;Data phy_u(15) /1.60219e-22_dp/
  Data phy_d(16) /'e'/;Data phy_n(16) /'ry'/;Data phy_u(16) /2.17991e-18_dp/
  Data phy_d(17) /'e'/;Data phy_n(17) /'mry'/;Data phy_u(17) /2.17991e-21_dp/
  Data phy_d(18) /'e'/;Data phy_n(18) /'hartree'/;Data phy_u(18) /4.35982e-18_dp/
  Data phy_d(19) /'e'/;Data phy_n(19) /'kcal/mol'/;Data phy_u(19) /6.94780e-21_dp/
  Data phy_d(20) /'e'/;Data phy_n(20) /'mhartree'/;Data phy_u(20) /4.35982e-21_dp/
  Data phy_d(21) /'e'/;Data phy_n(21) /'kj/mol'/;Data phy_u(21) /1.6606e-21_dp/
  Data phy_d(22) /'e'/;Data phy_n(22) /'hz'/;Data phy_u(22) /6.6262e-34_dp/
  Data phy_d(23) /'e'/;Data phy_n(23) /'thz'/;Data phy_u(23) /6.6262e-22_dp/
  Data phy_d(24) /'e'/;Data phy_n(24) /'cm-1'/;Data phy_u(24) /1.986e-23_dp/
  Data phy_d(25) /'e'/;Data phy_n(25) /'cm^-1'/;Data phy_u(25) /1.986e-23_dp/
  Data phy_d(26) /'e'/;Data phy_n(26) /'cm**-1'/;Data phy_u(26) /1.986e-23_dp/
  Data phy_d(27) /'f'/;Data phy_n(27) /'N'/;Data phy_u(27) /1.0_dp/
  Data phy_d(28) /'f'/;Data phy_n(28) /'ev/ang'/;Data phy_u(28) /1.60219e-9_dp/
  Data phy_d(29) /'f'/;Data phy_n(29) /'ry/bohr'/;Data phy_u(29) /4.11943e-8_dp/
  Data phy_d(30) /'l'/;Data phy_n(30) /'cm'/;Data phy_u(30) /1.e-2_dp/
  Data phy_d(31) /'p'/;Data phy_n(31) /'pa'/;Data phy_u(31) /1.0_dp/
  Data phy_d(32) /'p'/;Data phy_n(32) /'mpa'/;Data phy_u(32) /1.e6_dp/
  Data phy_d(33) /'p'/;Data phy_n(33) /'gpa'/;Data phy_u(33) /1.e9_dp/
  Data phy_d(34) /'p'/;Data phy_n(34) /'atm'/;Data phy_u(34) /1.01325e5_dp/
  Data phy_d(35) /'p'/;Data phy_n(35) /'bar'/;Data phy_u(35) /1.e5_dp/
  Data phy_d(36) /'p'/;Data phy_n(36) /'mbar'/;Data phy_u(36) /1.e11_dp/
  Data phy_d(37) /'p'/;Data phy_n(37) /'ry/bohr**3'/;Data phy_u(37) /1.47108e13_dp/
  Data phy_d(38) /'p'/;Data phy_n(38) /'ev/ang**3'/;Data phy_u(38) /1.60219e11_dp/
  Data phy_d(39) /'c'/;Data phy_n(39) /'c'/;Data phy_u(39) /1.0_dp/
  Data phy_d(40) /'c'/;Data phy_n(40) /'e'/;Data phy_u(40) /1.602177e-19_dp/
  Data phy_d(41) /'d'/;Data phy_n(41) /'C*m'/;Data phy_u(41) /1.0_dp/
  Data phy_d(42) /'d'/;Data phy_n(42) /'D'/;Data phy_u(42) /3.33564e-30_dp/
  Data phy_d(43) /'d'/;Data phy_n(43) /'debye'/;Data phy_u(43) /3.33564e-30_dp/
  Data phy_d(44) /'d'/;Data phy_n(44) /'e*bohr'/;Data phy_u(44) /8.47835e-30_dp/
  Data phy_d(45) /'d'/;Data phy_n(45) /'e*ang'/;Data phy_u(45) /1.602177e-29_dp/
  Data phy_d(46) /'mom'/;Data phy_n(46) /'kg*m**2'/;Data phy_u(46) /1.0_dp/
  Data phy_d(47) /'mom'/;Data phy_n(47) /'ry*fs**2'/;Data phy_u(47) /2.1799e-48_dp/
  Data phy_d(48) /'ef'/;Data phy_n(48) /'v/m'/;Data phy_u(48) /1.0_dp/
  Data phy_d(49) /'ef'/;Data phy_n(49) /'v/nm'/;Data phy_u(49) /1.e9_dp/
  Data phy_d(50) /'ef'/;Data phy_n(50) /'v/ang'/;Data phy_u(50) /1.e10_dp/
  Data phy_d(51) /'ef'/;Data phy_n(51) /'v/bohr'/;Data phy_u(51) /1.8897268e10_dp/
  Data phy_d(52) /'ef'/;Data phy_n(52) /'ry/bohr/e'/;Data phy_u(52) /2.5711273e11_dp/
  Data phy_d(53) /'ef'/;Data phy_n(53) /'har/bohr/e'/;Data phy_u(53) /5.1422546e11_dp/
  Data phy_d(54) /'e'/;Data phy_n(54) /'k'/;Data phy_u(54) /1.38066e-23_dp/

Contains

  Subroutine esdf_init(filename)

    Character(*), Intent(in) :: filename

    ! Local

    Integer(I4B), Parameter :: ncomm=3,ndiv=3
    Integer(I4B) :: unit,ierr,i,j,ic,nt,ndef,nread,itemp,itemp2
    Character(llength) :: cjunk,ctemp
    Character(1) :: comment(ncomm),divide(ndiv)
    Logical :: inblock

!    Integer :: temp1!, Index !  DJR:  Had to change this to get
                              !        it to compile under linux.

    ! Define comment characters

    Data comment /'#',';','!'/
    Data divide /' ','=',':'/

    ! "reduce" the keyword list for comparison

    Do i = 1,numkw
       ctemp = kw_label(i)
       kw_label(i) = esdf_reduce(ctemp)
    End Do

    ! Open the esdf file

    Call esdf_file(unit,filename,ierr)
    cjunk = 'Unable to open main input file "'//Trim(filename)//'"'

    If(ierr.Eq.1) Then
       Write(*,*) 'ESDF WARNING: '//Trim(cjunk)//' - using defaults'
       nread = 0
    Else
       nread = Huge(1)
    Endif

    ! Count the number of records (excluding blank lines and commented lines)

    nrecords = 0

    Do i=1,nread
       Read(unit,'(a)',End=100) cjunk
       Do j=1,ncomm
          ic=Index(cjunk,comment(j))
          If(ic.Gt.0) cjunk(ic:) = ' '
       End Do
       If(Len_trim(cjunk).Gt.0) Then
          nrecords = nrecords + 1
       Endif
    End Do
100 Rewind(unit)

    ! Allocate the array to hold the records and tokens

    Allocate(llist(nrecords),block_data(nrecords),tlist(llength,nrecords),&
         warns(nrecords))

    ! Set the number of warnings to zero

    nwarns = 0 ; warns = ' '

    ! Read in the records

    nrecords = 0
    Do i=1,nread
       Read(unit,'(a)',End=101) cjunk
       Do j=1,ncomm
          ic=Index(cjunk,comment(j))
          If(ic.Gt.0) cjunk(ic:) = ' '
       End Do

       If(Len_trim(cjunk).Gt.0) Then
          nrecords=nrecords+1
          llist(nrecords) = Adjustl(cjunk)
       Endif
    End Do
101 Close(unit)

    ! Now read in the tokens from llist

    tlist = ' '

    Do i=1,nrecords
       ctemp = llist(i)
       nt=0
       Do While(Len_trim(ctemp).Gt.0)
 !  apparently this a hack to make it compile under linux 
!          temp1=Index(ctemp,divide(1))
!          ic = Minval(temp1)
!          ic = Minval(Index(ctemp,divide),mask=Index(ctemp,divide)>0)
          ic = len_trim(ctemp)+1
          do itemp=1,size(divide)
             itemp2 = Index(ctemp,divide(itemp))
             if(itemp2.eq.0) itemp2=len_trim(ctemp)+1
             if(itemp2.lt.ic) ic=itemp2
          end do
          If(ic.Gt.1) Then
             nt=nt+1
             tlist(nt,i) = Adjustl(ctemp(:ic-1))
          Endif
          ctemp = Adjustl(ctemp(ic+1:))
       End Do
    End Do

    ! Check if any of the "labels" in the input file are unrecognised

    inblock=.False.
    Do i=1,nrecords
       ! Check if we are in a block
       If(esdf_reduce(tlist(1,i)).Eq.'begin') Then
          inblock = .True.
          ! Check if block label is recognised
          If((Count(esdf_reduce(tlist(2,i)).Eq.kw_label).Eq.0)) Then
             ctemp='Label "'//Trim(esdf_reduce(tlist(2,i)))//&
                  &'" not in keyword list'
             If(Count(ctemp.Eq.warns).Eq.0) Call esdf_warn(ctemp) 
          Endif
          ! Check if "label" is multiply defined in the input file
          ndef=0
          Do j=1,nrecords
             If(esdf_reduce(tlist(2,i)).Eq.esdf_reduce(tlist(2,j))) ndef=ndef+1
          End Do
          ctemp='Label "'//Trim(esdf_reduce(tlist(2,i)))//&
               &'" is multiply defined in the input file. '
          If((ndef.Gt.2).And.(Count(ctemp.Eq.warns).Eq.0))&
               Call esdf_warn(ctemp)
       Endif
       ! Check it is in the list of keywords
       If((Count(esdf_reduce(tlist(1,i)).Eq.kw_label).Eq.0)&
            .And.(.Not.inblock)) Then
          ctemp='Label "'//Trim(esdf_reduce(tlist(1,i)))//&
               &'" not in keyword list'
          If(Count(ctemp.Eq.warns).Eq.0) Call esdf_warn(ctemp)
       Endif
       If(.Not.inblock) Then
          ! Check if "label" is multiply defined in the input file
          ndef=0
          Do j=1,nrecords
             If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(tlist(1,j))) ndef=ndef+1
          End Do
          ctemp='Label "'//Trim(esdf_reduce(tlist(1,i)))//&
               &'" is multiply defined in the input file. '
          If((ndef.Gt.1).And.(Count(ctemp.Eq.warns).Eq.0)) &
               Call esdf_warn(ctemp)
       Endif
       ! Check if we have left a block
       If(esdf_reduce(tlist(1,i)).Eq.'end') inblock= .False.

    End Do

  End Subroutine esdf_init

  !  
  ! Return the string attached to the "label"
  !

  Function esdf_string(label,default)

    Character(*), Intent(in) :: label,default
    Character(llength) :: esdf_string

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp

    ! Check "label" is defined

    Call esdf_lblchk(label,'T')

    ! Set to default

    esdf_string = default

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          esdf_string = llist(i)(Index(llist(i),Trim(tlist(2,i))):)
          Exit
       Endif

    End Do

    Return

  End Function esdf_string

  !  
  ! Return the integer attached to the "label"
  !

  Function esdf_integer(label,default)

    Integer(I4B), Intent(in) :: default
    Character(*), Intent(in) :: label
    Integer(I4B) :: esdf_integer

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp

    ! Check "label" is defined

    Call esdf_lblchk(label,'I')

    ! Set to default

    esdf_integer = default

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,err=100) esdf_integer
          Exit
       Endif

    End Do

    Return


100 ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_integer'
    Call esdf_die(ctemp)

  End Function esdf_integer

  !  
  ! Return the single precisioned value attached to the "label"
  !

  Function esdf_single(label,default)

    Real(SP), Intent(in) :: default
    Character(*), Intent(in) :: label
    Real(SP) :: esdf_single

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp

    ! Check "label" is defined

    Call esdf_lblchk(label,'S')

    ! Set to default

    esdf_single = default

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,err=100) esdf_single
          Exit
       Endif

    End Do

    Return

100 ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_single'
    Call esdf_die(ctemp)

  End Function esdf_single

  !  
  ! Return the double precisioned value attached to the "label"
  !

  Function esdf_double(label,default)

    Real(DP), Intent(in) :: default
    Character(*), Intent(in) :: label
    Real(DP) :: esdf_double

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp

    ! Check "label" is defined

    Call esdf_lblchk(label,'D')
    
    ! Set to default

    esdf_double = default

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,err=100) esdf_double
          Exit
       Endif

    End Do  
 
    Return

100 esdf_double = default
    ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_double'
    Call esdf_die(ctemp)

  End Function esdf_double

  !  
  ! Return the double precisioned physical value attached to the "label"
  ! Units converted to "dunit"
  !

  Function esdf_physical(label,default,dunit)

    Real(DP), Intent(in) :: default
    Character(*), Intent(in) :: label,dunit
    Real(DP) :: esdf_physical

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp,iunit

    ! Check "label" is defined

    Call esdf_lblchk(label,'P')

    ! Set to default

    esdf_physical = default

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,err=100) esdf_physical
          iunit = dunit//repeat(' ',llength-len(dunit))
          Read(tlist(3,i),*,err=100,end=103) iunit
          esdf_physical = esdf_convfac(iunit,dunit)*esdf_physical
          Exit
       Endif

    End Do
    
    Return

103 iunit = dunit//repeat(' ',llength-len(dunit))
    esdf_physical = esdf_convfac(iunit,dunit)*esdf_physical

    Return

100 esdf_physical = default
    ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_physical'
    Call esdf_die(ctemp)

  End Function esdf_physical

  !
  ! Is the "label" defined in the input file
  !

  Function esdf_defined(label)

    Character(*), Intent(in) :: label
    Logical :: esdf_defined

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp

    ! Check "label" is defined

    Call esdf_lblchk(label,'E')

    ! Set to default

    esdf_defined = .False.

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          esdf_defined = .True.
          Exit
       Endif

    End Do

    Return

  End Function esdf_defined
  !
  ! Is the "label" defined in the input file
  !

  Function esdf_boolean(label,default)

    Character(*), Intent(in) :: label
    Logical, Intent(in) :: default
    Logical :: esdf_boolean

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp,positive(3),negative(3)

    Data positive /'yes','true','t'/
    Data negative /'no','false','f'/

    ! Check "label" is defined

    Call esdf_lblchk(label,'L')

    ! Set to default

    esdf_boolean = default

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          If(Len_trim(tlist(2,i)).Eq.0) Then
             esdf_boolean = .True.
             Exit
          Endif
          If(Any(Index(positive,esdf_reduce(tlist(2,i))).Gt.0)) Then
             esdf_boolean = .True.
             Exit
          Endif
          If(Any(Index(negative,esdf_reduce(tlist(2,i))).Gt.0)) Then
             esdf_boolean = .False.
             Exit
          Endif
          Call esdf_die('Unable to parse boolean value')

       Endif

    End Do

    Return

  End Function esdf_boolean

  Function esdf_block(label,nlines)

    Character(*), Intent(in) :: label
    Integer(I4B), Intent(out) :: nlines
    Logical :: esdf_block

    ! Local

    Integer(I4B) :: i,j
    Character(llength) :: ctemp

    ! Check "label" is defined

    Call esdf_lblchk(label,'B')

    ctemp ='Block "'//Trim(esdf_reduce(label))//'" not closed correctly '
    
    esdf_block=.False.

    nlines = 0

    Do i=1,nrecords
       If((esdf_reduce(tlist(1,i)).Eq.esdf_reduce('begin'))&
            .And.(esdf_reduce(tlist(2,i)).Eq.esdf_reduce(label))) Then
          esdf_block = .True.
          Do While(esdf_reduce(tlist(1,i+nlines+1))&
               .Ne.esdf_reduce('end'))
             nlines=nlines+1
             If(nlines+i.Gt.nrecords) Call esdf_die(ctemp)
             block_data(nlines)=llist(i+nlines)
          End Do

          If(esdf_reduce(tlist(2,i+nlines+1)).Ne.esdf_reduce(label))&
               Call esdf_die(ctemp)
          Exit
       Endif
    End Do

    Return

  End Function esdf_block

  !
  ! Reduce the string to lower case and remove punctuation
  !

  Function esdf_reduce(string)

    Character(*), Intent(in) :: string
    Character(llength) :: esdf_reduce

    ! Local

    Integer(I4B), Parameter :: npunct=2
    Integer(I4B) :: iA,iZ,ishift,ic,i,ln
    Character(llength) :: ctemp
    Character(1) :: punct(npunct)
    Logical :: keep_going

    ! Define the punctuation to be removed

    Data punct /'.','-'/

    ! Initialise system dependant bounds in collating sequence

    iA = Ichar('A');iZ = Ichar('Z') 
    ishift = Ichar('a')-iA 

    ! Initialise output

    ln = len(string)

    if(ln.lt.llength) then
       esdf_reduce(1:ln) = string(1:ln) ; esdf_reduce(ln+1:)=' '
    else
       esdf_reduce(1:llength) = string(1:llength)
    end if

    ! Drop all upper case characters to lower case

    Do i=1,llength
       ic = Ichar(esdf_reduce(i:i))
       If((ic.Ge.iA).And.(ic.Le.iZ)) esdf_reduce(i:i) = Char(ishift+ic) 
    Enddo

    ! Now remove punctuation

    Do i=1,npunct

       keep_going=.True.

       Do While(keep_going)
          ic = Index(esdf_reduce,punct(i))
          If(ic.Gt.0) Then
             ctemp = esdf_reduce
             esdf_reduce(ic:)=ctemp(ic+1:)
          Else
             keep_going=.False.
          End If
       End Do

    End Do

    esdf_reduce = Trim(Adjustl(esdf_reduce))

  End Function esdf_reduce

  !
  ! Find the conversion factor between physical units
  !

  Function esdf_convfac(from,to)

    Character(*), Intent(in) :: from,to
    Real(DP) :: esdf_convfac

    ! Local

    Integer(I4B) :: i,ifrom,ito
    Character(llength) :: ctemp

    ! Find the index numbers of the from and to units

    ifrom = 0 ; ito = 0
    Do i=1,nphys
       If(esdf_reduce(from).Eq.esdf_reduce(phy_n(i))) ifrom = i
       If(esdf_reduce(to).Eq.esdf_reduce(phy_n(i))) ito = i
    End Do

    ! Check that the units were recognised

    If(ifrom.Eq.0) Then
       ctemp = 'Units not recognised in input file : '//Trim(esdf_reduce(from))
       Call esdf_die(ctemp)
    Endif

    If(ito.Eq.0) Then
       ctemp = 'Units not recognised in Program : '//Trim(esdf_reduce(to))
       Call esdf_die(ctemp)
    Endif

    ! Check that from and to are of the same dimensions

    If(phy_d(ifrom).Ne.phy_d(ito)) Then
       ctemp = 'Dimensions Do not match : '//Trim(esdf_reduce(from))&
            & //' vs '//Trim(esdf_reduce(to))
       Call esdf_die(ctemp)
    Endif

    ! Set the conversion factor

    esdf_convfac = phy_u(ifrom)/phy_u(ito)

  End Function esdf_convfac

  ! 
  ! Find an unused i/o unit
  !

  Function esdf_unit(ierr)
    Integer(I4B), Intent(out) :: ierr
    Integer(I4B) :: esdf_unit
    ! Local 
    Logical :: op
    ierr=0
    Do esdf_unit=10,99
       Inquire(unit=esdf_unit,opened=op,err=100)
       If(.Not.op) Return
    End Do
    Call esdf_warn('Unable to find a free i/o unit using esdf_unit')
    ierr = 1
    Return
100 Call esdf_die('Error opening files by esdf_unit')
  End Function esdf_unit

  !
  ! Open an old file
  !

  Subroutine esdf_file(unit,filename,ierr)
    Character(*), Intent(in) :: filename
    Integer(I4B), Intent(out) :: unit,ierr
    Logical :: ex
    unit = esdf_unit(ierr)
    If(ierr.Gt.0) Return
    Inquire(file=Trim(filename),exist=ex,err=100)
    If(.Not.ex) Goto 100
    Open(unit=unit,file=Trim(filename),form='formatted',status='old',err=100)
    Return
100 ierr=1
    Return
  End Subroutine esdf_file

  ! Open a new file

  Subroutine esdf_newfile(unit,filename,ierr)
    Character(*), Intent(in) :: filename
    Integer(I4B), Intent(out) :: unit,ierr
    unit = esdf_unit(ierr)
    If(ierr.Gt.0) Return
    Open(unit=unit,file=Trim(filename),form='formatted',status='replace',err=100)
    Return
100 ierr=1
    Return
  End Subroutine esdf_newfile

  !
  ! Check that the label is known, and used correctly
  !

  Subroutine esdf_lblchk(string,typ)
    Character(*), Intent(in) :: string
    Character(1), Intent(in) :: typ
    ! Local
    Character(llength) :: ctemp
    Character(1) :: tp
    Integer(I4B) :: i
    ! Check if label is recognised
    i=Count(esdf_reduce(string).Eq.kw_label)
    ctemp = 'Label "'//Trim(esdf_reduce(string))//'" not recognised in&
         & keyword list'
    If(i.Eq.0) Call esdf_die(ctemp)
    ctemp = 'Label "'//Trim(esdf_reduce(string))//'" is multiply defined'
    If(i.Gt.1) Call esdf_die(ctemp)
    ctemp = 'Label "'//Trim(esdf_reduce(string))//'" has been used with the wrong type'
    tp = ' '
    i=0
    Do While(tp.Eq.' ')
       i=i+1
       If(esdf_reduce(string).Eq.kw_label(i)) tp=kw_typ(i)
    End Do
    If(typ.Ne.tp) Call esdf_die(ctemp)
  End Subroutine esdf_lblchk


  Subroutine esdf_help(helpword,searchword)

    Implicit None

    Character(*) :: helpword,searchword

    ! Local

    Integer(I4B)  :: i,indx,indx2,ln
    Character(20) :: ctyp,clev
    Character(60) :: title,fm
    Character(80) :: ctemp
    Character(1)  :: cl

    helpword = esdf_reduce(helpword)
    searchword = esdf_reduce(searchword)

    If(esdf_reduce(helpword).Eq.'search') Then

       If(Len_trim(searchword).Lt.1) Call esdf_die('help: "searchword" is empty')

       ! Search for useful keywords

       Do i=1,numkw
          If((Index(kw_label(i),Trim(searchword)).Gt.0).Or.&
               (Index(kw_dscrpt(i),Trim(searchword)).Gt.0)) Then 
             indx=Index(kw_dscrpt(i),'!*')-1
             If(indx.Eq.-1) &
                  Call esdf_die('help: keyword description incorrectly formatted')
             title = kw_dscrpt(i)(1:indx)
             ln=Len(Trim(title))
             If(ln.Gt.80) Call esdf_die('help: keyword title too long')

             Write (*,*) kw_label(i),Trim(title)

          End If
       End Do
       Call esdf_die('help: search done')


    Endif

    ! All keywords, short description

    If('all'.Eq.helpword) Then
       Do i=1,numkw
          If(Len_trim(kw_label(i)).Gt.0) Then 
             indx=Index(kw_dscrpt(i),'!*')-1
             If(indx.Eq.-1) &
                  Call esdf_die('help: keyword description incorrectly formatted')
             title = kw_dscrpt(i)(1:indx)
             ln=Len(Trim(title))
             If(ln.Gt.80) Call esdf_die('help: keyword title too long')

             Write (*,*) kw_label(i),Trim(title)

          End If
       End Do
       Call esdf_die('help: all done')
    End If

    ! All specific levels of keywords

    If(Any((/'basic ','inter ','expert','dummy '/).Eq.helpword)) Then

       Select Case(helpword)
       Case('basic')  ; cl = 'B'
       Case('inter')  ; cl = 'I'
       Case('expert') ; cl = 'E'
       Case('dummy')  ; cl = 'D'
       End Select

       Do i=1,numkw
          If(kw_typ(i)(3:3).Eq.cl) Then 
             indx=Index(kw_dscrpt(i),'!*')-1
             If(indx.Eq.-1) &
                  Call esdf_die('help: keyword description incorrectly formatted')
             title = kw_dscrpt(i)(1:indx)
             ln=Len(Trim(title))
             If(ln.Gt.80) Call esdf_die('help: keyword title too long')

             Write (*,*) kw_label(i),Trim(title)

          End If
       End Do
       Call esdf_die('help: level done')
    End If

    ! More information about a specific keyword

    If(.Not.Any(kw_label.Eq.helpword)) &
         Call esdf_die('help: keyword not recognised')
    If(Count(kw_label.Eq.helpword).Gt.1) &
         Call esdf_die('help: keyword entry duplicated')
    Do i=1,numkw
       If(kw_label(i).Eq.helpword) Then 
          indx=Index(kw_dscrpt(i),'!*')+1
          If(indx.Eq.1) &
               Call esdf_die('help: keyword description incorrectly formatted')
          title = kw_dscrpt(i)(1:indx)
          ln=Len(Trim(title))
          If(ln.Gt.80) &
               Call esdf_die('help: keyword title too long')
          If(ln.Le.9) Write(fm,'("(",i2,"x,a",i1,")")') 40-ln/2,ln
          If(ln.Gt.9) Write(fm,'("(",i2,"x,a",i2,")")') 40-ln/2,ln
          Write (*,fm) Trim(title)
          Write (*,*)
          Select Case(kw_typ(i)(1:1))
          Case('I') ; ctyp ='Integer'
          Case('S') ; ctyp ='Single Precision'
          Case('D') ; ctyp ='Double Precision'
          Case('P') ; ctyp ='Physical'
          Case('T') ; ctyp ='String'
          Case('E') ; ctyp ='Defined'
          Case('B') ; ctyp ='Block'
          Case('L') ; ctyp ='Boolean'   
          End Select
          Select Case(kw_typ(i)(3:3))
          Case('B') ; clev ='Basic'
          Case('I') ; clev ='Intermediate'
          Case('E') ; clev ='Expert'
          Case('D') ; clev ='Dummy'
          End Select
          Write (fm,'(a,i2,a)') '("Type: ",a,',&
               78-(7+Len_trim(clev))-(6+Len_trim(ctyp)),'x," Level: ",a)'
          Write (ctemp,fm) Trim(ctyp),Trim(clev)
          Write (*,'(a)') Trim(ctemp)
          Write (*,*)
          indx=indx+1
          ln = Len(Trim(kw_dscrpt(i)))
          Do While (indx.Lt.ln)
             ctemp = kw_dscrpt(i)(indx:Min(indx+80,ln))
             indx2=Index(ctemp,' ',back=.True.)
             Write (ctemp,'(a)') Adjustl(ctemp(:indx2))
             Write (*, '(A)') Trim(ctemp)
             indx=indx+Len(ctemp(:indx2))
          End Do

       End If
    End Do

    Call esdf_die('help: done')

  End Subroutine esdf_help

  ! 
  ! Stop execution due to an error cause by esdf
  !

  Subroutine esdf_die(string)
    Character(*), Intent(in) :: string
    Write (*,'(a,a)') ' ESDF ERROR: ',Trim(string)
    Write (*,'(a)') ' Stopping now'
    Call mystop    
  End Subroutine esdf_die

  ! 
  ! Warning due to an error cause by esdf
  !

  Subroutine esdf_warn(string)
    Character(*), Intent(in) :: string
    nwarns=nwarns+1
    warns(nwarns) = string
  End Subroutine esdf_warn

  !
  ! Dump the warnings to screen
  !

  Subroutine esdf_warnout
    Integer(I4B) :: i
    Do i=1,nwarns
       Write (*,*) 'ESDF WARNING: '//Trim(warns(i))
    End Do
  End Subroutine esdf_warnout

  !
  ! Deallocate the data arrays --- call this before re-initialising
  !

  Subroutine esdf_close
    Deallocate(llist,tlist,block_data)
  End Subroutine esdf_close

End Module esdf
 
