subroutine mol_dyn(plen,MD,x,g,pwcount,pw_params,t0, p,&
                     myproc,  relaxstress, relaxforce, &
                     iwfnreuse_c, adjustpressure,altkpoints,&
                     bands,altbands,syms,energs,crys)
  !
  !     ---------------------------------------------------------------
  !    
  include 'use.h'
  IMPLICIT NONE
  include 'flibcalls.ph'
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  type(pw_parameter), intent(inout) :: pw_params ! plane wave parameters
  type(kpoint), intent(inout) :: altkpoints  ! alternative BZ integration kpts
  type(band), intent(inout) :: bands         ! eigenvalues, occ. numbers etc
  type(band), intent(inout) :: altbands      ! another set of eigenvalues etc..
  type(symmetry), intent(inout) :: syms      ! symmetry operations
  type(energy), intent(inout) :: energs      ! energies and other information
  type(crystal), intent(inout) :: crys       ! crystal structure
  type(molecular_dynamics), intent(inout) :: MD  ! molecular dynamics structure
  !
  !     INPUT:
  !     ------
  !
  integer, intent(in) :: &
        plen, &
       myproc
  real(dp), intent(in) :: & 
            relaxstress, relaxforce  
  !
  !     INPUT/OUTPUT:
  !     ------------
  !
  integer, intent(inout) ::  iwfnreuse_c,adjustpressure,pwcount
  real(dp), intent(inout) :: p,t0
  real(dp), intent(inout) :: x(plen)
  real(dp), intent(inout) :: &
       g(plen)  
  !
  !     -------------------------- local variables -----------------------
  !
  real(dp) t_step_fs ,KE_scal_unit,&
              F_scal_unit,temp,ext_entropy,  ext_kin, eta_new,eta_old,&
             zeta_new,zeta_old, matinvert, tmp(3) , etot
  real(dp) g_temp(plen), g_temp_old(plen),gold(plen),xold(plen)
  integer ip,iter,i,j,k,n,extrap
  real(dp) ekin_ion, &
                del_r(plen,4) 
  real(dp) feps(3,3),lvect_inv(3,3),atam1(3,3),f1,average,funeval
  real(dp) in_pr11,in_pr12,in_pr13,in_pr22,in_pr23,in_pr33, inv_fac,Temp_curr

  t_step_fs =MD%time_step*1d-15 

  KE_scal_unit = ( m_to_kg/(13.6058*ev_to_j) ) * (b_to_m/t_step_fs)**2

  F_scal_unit=(MD%time_step**2)*((13.6058*ev_to_j)/m_to_kg) *(1d-30/(b_to_m**2))

  gold=dzero
  g_temp=dzero
  !
  ! set mass of heat bath for canonical MD
  !
  if (MD%ensemble .eq. 2 )  then 
   if (MD%Q_mass .eq. 0d0) then
    MD%Q_mass= &
         Kb_Ryd * MD%n_deg_freedom*MD%temp/10
   end if
  end if
  write(34,*) MD%Q_mass, '=  Mass of heat bath in Ryd'

  ! convert to energy*time_step^2
  MD%Q_mass= MD%Q_mass * MD%time_step**2

  zeta_new=dzero

  do iter = 1, MD%itmax - 1  

    pw_params%ilog = iter


    lvect_inv = transpose(crys%lvec)
    temp = matinvert(lvect_inv(1, 1))
    feps =  crys%lvec
    temp = matinvert(feps(1, 1))
    call matmul(atam1,lvect_inv,feps )

    call mdgemm('N','N',3,crys%mxdatm*crys%ntype+3,3,done,atam1(1,1),3,g(1),&
               3,dzero,g_temp_old(1),3)  

    xold=x
    ip=9
    do i=1,crys%ntype
      do j=1,crys%natom(i)
        do k=1,3
          ip=ip+1
          g(ip)=-g_temp_old(ip)*F_scal_unit /(2*crys%mass(i))
        end do
      end do
    end do

    if (iter .eq. 1 ) then

      gold=g
      
      if (MD%Q_mass>0) then
        zeta_new=dzero 
      end if
        
      eta_new= dhalf * zeta_new       

      ip=9
      do i=1,crys%ntype
       do j=1,crys%natom(i)
        do k=1,3
          ip=ip+1
          g_temp(ip)=done / (done + zeta_new/dtwo) * (MD%velocity(ip) + g(ip) )
        end do
       end do
      end do

    else

      call EKIN(0,plen,ekin_ion,crys%ntype,crys%natom,crys%mass,&
                         MD%time_step,crys%lvec,g_temp)

      if (MD%Q_mass>0) then
        zeta_new=zeta_old + done/MD%Q_mass * &
             ( 2*ekin_ion-MD%n_deg_freedom*Kb_Ryd*MD%temp)
        eta_new=eta_old + dhalf*(zeta_new+zeta_old)
      else
        zeta_new=dzero
        eta_new=dzero
      end if

      ip=9
      do i=1,crys%ntype
        do j=1,crys%natom(i)
          do k=1,3 
            ip=ip+1

            MD%velocity(ip)=g_temp(ip)+g(ip)-dhalf*zeta_new*g_temp(ip)

            g_temp(ip)=done/(done + zeta_new/dtwo) * (MD%velocity(ip) + g(ip) )
            gold(ip)=g(ip)
             
          end do
        end do
      end do

    end if

    ext_kin=dhalf*MD%Q_mass*(zeta_new)**2
    ext_entropy =MD%n_deg_freedom*Kb_Ryd*MD%temp*(eta_new)

    zeta_old=zeta_new
    eta_old=eta_new

    ip=9
    do i=1,crys%ntype
      do j=1,crys%natom(i)
        do k=1,3 
          ip=ip+1
          x(ip)=x(ip)+g_temp(ip)
        end do
      end do
    end do

    call EKIN(1,plen,ekin_ion,crys%ntype,crys%natom,crys%mass,&
                         MD%time_step,crys%lvec,MD%velocity)
    !
    ! OUTPUT TEMP and TOTAL ENERGY
    !
    Temp_curr=dtwo*ekin_ion/(real(MD%n_deg_freedom,dp)*Kb_Ryd)

    write(34,*) Temp_curr, 'TEMP'
    write(34,33) ekin_ion + energs%total+ ext_entropy+ ext_kin,&
               energs%total, ekin_ion ,ext_kin,ext_entropy
    call myflush(34)
    ! 
    ! extrapolate positions and get alpha,beta for extrapolating wavefunctions
    !
    if (MD%extrapolation .ge. 0) then 

      extrap=mod( MD%extrapolation , 10)

      if (iter .eq. 1) then
        del_r(:,1)=x(:)-xold(:)
        pw_params%MD_alpha=dzero
        pw_params%MD_beta=dzero         
      else if (iter .eq. 2 .or. extrap .eq. 0) then
        del_r(:,2)= del_r(:,1)
        del_r(:,1)=x(:)-xold(:)
        pw_params%MD_alpha= done
        pw_params%MD_beta=dzero
      else 
        del_r(:,3)= del_r(:,2)
        del_r(:,2)= del_r(:,1)
        del_r(:,1)=x(:)-xold(:)

        in_pr11 = mddot(plen-9, del_r(10,1),1, del_r(10,1), 1)
        in_pr12 = mddot(plen-9, del_r(10,1),1, del_r(10,2), 1)
        in_pr13 = mddot(plen-9, del_r(10,1),1, del_r(10,3), 1)
        in_pr22 = mddot(plen-9, del_r(10,2),1, del_r(10,2), 1)
        in_pr23 = mddot(plen-9, del_r(10,2),1, del_r(10,3), 1)
        in_pr33 = mddot(plen-9, del_r(10,3),1, del_r(10,3), 1)

        inv_fac = 1/( in_pr22* in_pr33 - in_pr23* in_pr23)

        if (mod(iter,extrap) .eq. 0) then
          pw_params%MD_alpha= inv_fac* (in_pr33*in_pr12 - in_pr23*in_pr13)
          pw_params%MD_beta= inv_fac* (-in_pr23*in_pr12 + in_pr22*in_pr13)
        else
          pw_params%MD_alpha=done
          pw_params%MD_beta=dzero
        end if

      end if 

    else

      pw_params%MD_alpha=dzero
      pw_params%MD_beta=dzero 
    end if

    write(34,*) pw_params%MD_beta, pw_params%MD_alpha, ' beta  alpha'
    call myflush(34)
    gold=g 

    f1 = funeval(t0, plen, x, p, g, myproc, relaxstress, relaxforce, &
           iwfnreuse_c, adjustpressure,pw_params,altkpoints,&
          bands,altbands,syms,energs,crys)
          pwcount = pwcount + 1
        
  end do  ! MD steps

33 FORMAT('TOT_E= ',F12.6,2x,'KIN_ION=',F12.6,2x,'POT_E= ',F12.6,2x,'KIN_EXT=',F12.6,2x,&
            'POT_EXT',F12.6)

  return
  end subroutine mol_dyn
  !
  !     ---------------------------------------------------------------
  ! 
  subroutine initialize_velocity(MD,plen,mxdatm,ntyp,natom,avec,mass,natom_tot,lvec)
  use constants
  use molecular_dynamics_module
  implicit none
  include 'flibcalls.ph' 

  integer mxdatm,ntyp,plen,natom_tot
  real(dp) avec(3,3)
  real(dp) mass(ntyp),lvec(3,3)
  type(molecular_dynamics) :: MD  ! molecular dynamics structure

  integer natom(ntyp)

  real(dp) rdm(2*mxdatm*ntyp)
  real(dp) velocity_tmp(plen)
  real(dp) phi,theta,energy,v_length(ntyp),tmp(3),average
  integer i,j,icnt_rdm,i_vel,n

  real(dp) unit_time,fact,ekin_ion,scale
 
  unit_time =MD%time_step*1d-15  

  fact= (M_to_KG/eV_to_J)*(B_to_M/unit_time)**2

! energy= 3/2 KbT 
  energy = (real(MD%n_deg_freedom,dp)/(real(natom_tot,dp)))*Kb_eV * MD%temp   ! * D-5 Joules
!

  do i=1,ntyp
! 1/2mv2 = 3/2kbT=energy - >  solve for v in m^2/s^2
    v_length(i) = sqrt (energy/(mass(i)*fact)) ! 1d2
  end do

  call random_number(rdm)

  velocity_tmp=dzero
  icnt_rdm=1
  i_vel=9
  do i=1,ntyp
    do j=1,natom(i)
      phi=rdm(icnt_rdm)*360
      theta=rdm(icnt_rdm+1)*180

      velocity_tmp(i_vel+1)=v_length(i)*cos(phi)*sin(theta)
      velocity_tmp(i_vel+2)=v_length(i)*sin(phi)*sin(theta)
      velocity_tmp(i_vel+3)=v_length(i)*cos(theta)

      icnt_rdm=icnt_rdm+2

      i_vel=i_vel+3
    end do
  end do

  call mdgemm('N','N',3,mxdatm*ntyp+3,3,done,avec(1,1),3,velocity_tmp(1),&
               3,dzero,MD%velocity(1),3)

  tmp = dzero
  average = dzero
  i_vel=9
  do i= 1,ntyp
    do n=1,natom(i)
      average = average + mass(i)
      do j=1,3
         i_vel=i_vel+1
         tmp(j) = tmp(j)+MD%velocity(i_vel)*mass(i)
      end do
    end do
  end do

  do j=1,3
    tmp(j)=-tmp(j)/average
  end do

  write(34,*) tmp(1),tmp(2),tmp(3), 'drift'
  call myflush(34)

  i_vel=9
  do i= 1,ntyp
    do n=1,natom(i)
      do j=1,3
        i_vel=i_vel+1
        MD%velocity(i_vel) = MD%velocity(i_vel) + tmp(j)
      end do
    end do
  end do

  call EKIN(0,plen,ekin_ion,ntyp,natom,mass,MD%time_step,lvec,MD%velocity)

  scale=sqrt(MD%temp*MD%n_deg_freedom*Kb_Ryd/(2*ekin_ion))

  i_vel=9
  do i= 1,ntyp
    do n=1,natom(i)
      do j=1,3
        i_vel=i_vel+1
        MD%velocity(i_vel) = MD%velocity(i_vel)*scale
      end do
    end do
  end do
  
  return
  end subroutine  initialize_velocity
  !
  !     ---------------------------------------------------------------
  ! 
  subroutine EKIN(ipr,plen,ekin_ion,ntyp,natom,mass,time_step,lvec,velocity)
  use constants
  implicit none
  include 'flibcalls.ph' 
  real(dp), intent(out) ::  ekin_ion
  integer, intent(in) :: ntyp,natom(ntyp),plen,ipr
  real(dp), intent(in) :: mass(ntyp),lvec(3,3),velocity(plen),time_step

  real(dp) velocity_tmp(plen)
  integer ip,i,j,n

  real(dp) t_step_fs ,KE_scal_unit, average,tmp(3)


  t_step_fs =time_step*1d-15 

  KE_scal_unit = ( m_to_kg/(ryd*ev_to_j) ) * (b_to_m/t_step_fs)**2

  call mdgemm('T','N',3,plen/3,3,done,lvec(1,1),3,velocity(1),&
               3,dzero,velocity_tmp(1),3)  

  ekin_ion=dzero
  ip=10
  do i=1,ntyp
    do j=1,natom(i)
       ekin_ion = ekin_ion + dhalf*mass(i) * &
         (velocity_tmp(ip)**2 + velocity_tmp(ip+1)**2 + &
                velocity_tmp(ip+2)**2)* KE_scal_unit   
       ip=ip+3
    end do
  end do

  tmp = dzero
  average = dzero
  ip=9
  do i= 1,ntyp
    do n=1,natom(i)
      average = average + mass(i)
      do j=1,3
        ip=ip+1
        tmp(j) = tmp(j)+velocity_tmp(ip)*mass(i)
      end do
    end do
  end do

  do j=1,3
    tmp(j)=-tmp(j)/average
  end do

  if (ipr > 0) write(34,*) tmp(1),tmp(2),tmp(3), 'drift'

  return
  end  subroutine EKIN
