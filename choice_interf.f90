module choice_interf
  !
  ! Module programmer: Tomas Grönstedt
  !
  ! v1.0 17.04.01  Date written.
  !
  ! Choice interface. Used by both standalone and gestpan integrated version.
  !
  use choice_physics
  use sim_precision
  use units_and_constants
  use choice_aux
  !
  implicit none
  !
  private 
  !
  public :: set_performance_choice, set_weight_choice, set_configuration_choice, & 
    & compute_noise_sources, interpolate_to_t_source, compute_flight_effects, set_rotational_speeds_choice, set_noise_choice, & 
    & compute_certification_data, get_string, get_inputNoise_ci, set_trajectory_subr, save_noise_points, deallocate_module_arrays, & 
    & certificationLimits, get_MTOW, get_no_engines, get_version_num
  !
  integer :: n_modules, nt
  character(len=cstr_len), dimension(:), allocatable :: modules ! the list of modules in the particular engine being analyzed
  logical, private :: output_file_open = .false.
  !
  integer, private :: output_unit = 0
  !
contains
  subroutine reset_choice_interf

    if( allocated(modules) ) deallocate(modules)

  end subroutine reset_choice_interf
  !
  subroutine get_version_num(fname)
    character(len=*), intent(in) :: fname
    if( .not. output_file_open ) call open_output_file(fname,output_unit)
    write(output_unit,'(2A)') '! Choice version: ',version_num
    !
    contains 
     subroutine open_output_file(fname,unit)
       character(len=*), intent(in) :: fname
       integer, intent(out) :: unit
     
       unit = get_unit_number() 
       open(unit,file=fname,status='unknown')
       
       output_file_open = .true.
       
     end subroutine open_output_file
  end  subroutine get_version_num
  !
  !
  subroutine save_noise_points(fname,opPoint)
    character(len=*), intent(in) :: opPoint  
    character(len=*), intent(in) :: fname
    real(kind=rp) :: kernel

    write(output_unit,'(A)') 'Operating point is '//opPoint
    write(output_unit,'(A,F16.4)') 'Fan inlet EPNL is ',EPNL_inlet_fan
    write(output_unit,'(A,F16.4)') 'Fan discharge EPNL is ',EPNL_discharge_fan
    write(output_unit,'(A,F16.4)') 'Inlet LPC EPNL is ',EPNL_inlet_lpc
    write(output_unit,'(A,F16.4)') 'LPT EPNL is ',EPNL_lpt
    write(output_unit,'(A,F16.4)') 'Comb EPNL is ',EPNL_comb
    write(output_unit,'(A,F16.4)') 'Caj EPNL is ',EPNL_caj
    write(output_unit,'(A,F16.4)') 'Airframe EPNL is ',EPNL_airfrm
    
    kernel = 10.0**(EPNL_inlet_fan/10.0_rp) + 10.0**(EPNL_discharge_fan/10.0_rp) + 10.0**(EPNL_inlet_lpc/10) + & 
        & 10.0**(EPNL_lpt/10.0_rp) + 10.0**(EPNL_comb/10.0_rp) + 10.0**(EPNL_caj/10.0_rp)
    
    write(output_unit,'(A,F16.4)') 'EPNL_tot is ',10.*log10(kernel);


  end subroutine save_noise_points
  !
  !
  !
  subroutine compute_certification_data()
  
  ! create total SPLp
  
     call compute_PNL()
     
     call compute_PNLT()
     
     call compute_EPNL()
  
  end subroutine compute_certification_data
  !
  ! Nopred validation data: EPNL_fanIn =  95.277362, EPNL_fanAft = 94.969237, EPNL_LPCin = 52.6732949
  !
  subroutine compute_EPNL()
     integer :: i 
     
     do i = 1, n_modules ! evaluate all noise sources in this engine model
        if( trim(modules(i)).eq.'Fan' ) then
            EPNL_inlet_fan = getEPNL(n_times,PNLT_inlet_fan)
            EPNL_discharge_fan = getEPNL(n_times,PNLT_discharge_fan)
        elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
            EPNL_inlet_lpc = getEPNL(n_times,PNLT_inlet_lpc)
        elseif( trim(modules(i)).eq.'Lpt' ) then
            EPNL_lpt = getEPNL(n_times,PNLT_lpt)
        elseif( trim(modules(i)).eq.'Comb' ) then
            EPNL_comb = getEPNL(n_times,PNLT_Comb)
        elseif( trim(modules(i)).eq.'cold_nozzle' ) then
            EPNL_caj = getEPNL(n_times,PNLT_caj)
        else
            ! no further modules implemented
        end if
     end do
     
     EPNL_airfrm = getEPNL(n_times,PNLT_airfrm)
     
  end subroutine compute_EPNL  
  !
  !
  !
  subroutine compute_PNL()
     real(kind=rp), dimension(:,:), allocatable :: wrk ! work array
     real(kind=rp), dimension(:),   allocatable :: vec1, vec2 ! work vectors
     integer :: i 
     
     call allocate_variables()
     
     do i = 1, n_modules ! evaluate all noise sources in this engine model
        if( trim(modules(i)).eq.'Fan' ) then
            PNL_inlet_fan = getPNL(nfreq,n_times,tmic,fobs,SPLp_inlet_fan,vec1,vec2,wrk)
            PNL_discharge_fan = getPNL(nfreq,n_times,tmic,fobs,SPLp_discharge_fan,vec1,vec2,wrk)
        elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
            PNL_inlet_lpc = getPNL(nfreq,n_times,tmic,fobs,SPLp_inlet_lpc,vec1,vec2,wrk)
        elseif( trim(modules(i)).eq.'Lpt' ) then
            PNL_lpt = getPNL(nfreq,n_times,tmic,fobs,SPLp_lpt,vec1,vec2,wrk)
        elseif( trim(modules(i)).eq.'Comb' ) then
            PNL_comb = getPNL(nfreq,n_times,tmic,fobs,SPLp_comb,vec1,vec2,wrk)
        elseif( trim(modules(i)).eq.'cold_nozzle' ) then
            PNL_caj = getPNL(nfreq,n_times,tmic,fobs,SPLp_caj,vec1,vec2,wrk)
        else
            ! no further modules implemented
        end if
     end do    
     PNL_airfrm = getPNL(nfreq,n_times,tmic,fobs,SPLp_airfrm,vec1,vec2,wrk)

     deallocate(wrk) ; deallocate(vec1) ; deallocate(vec2)
     
  contains
     subroutine allocate_variables()
       integer :: i 
       
       allocate(wrk(nfreq,n_times)) ; allocate(vec1(n_times)); allocate(vec2(n_times))
     
       do i = 1, n_modules ! evaluate all noise sources in this engine model
          if( trim(modules(i)).eq.'Fan' ) then
              allocate(PNL_inlet_fan(n_times)) ; allocate(PNL_discharge_fan(n_times))
          elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
              allocate(PNL_inlet_lpc(n_times)) 
          elseif( trim(modules(i)).eq.'Lpt' ) then
              allocate(PNL_lpt(n_times)) 
          elseif( trim(modules(i)).eq.'Comb' ) then
              allocate(PNL_comb(n_times)) 
          elseif( trim(modules(i)).eq.'cold_nozzle' ) then
              allocate(PNL_caj(n_times)) 
          else
              ! no further modules implemented
          end if
       end do
       
       allocate(PNL_airfrm(n_times)) 
     
     end subroutine allocate_variables
  end subroutine compute_PNL
  !
  !
  !
  function get_no_engines()
    integer :: get_no_engines 
    
    get_no_engines = no_engines
  
  end function get_no_engines
  !
  !
  !
  function get_MTOW()
    real(kind=rp) :: get_MTOW ! get maximum take-off weight for aircraft in tonnes
    
    get_MTOW = total_weight_airfrm/1000.0_rp
  
  end function get_MTOW
  !
  !
  !
  subroutine compute_PNLT()
     real(kind=rp), dimension(:), allocatable :: vec1 ! work vector
     real(kind=rp), dimension(:,:), allocatable :: wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7 ! work arrays
     integer :: i 
     
     call allocate_variables()
     
     do i = 1, n_modules ! evaluate all noise sources in this engine model
        if( trim(modules(i)).eq.'Fan' ) then
            PNLT_inlet_fan = getPNLT(nfreq,n_times,fband,PNL_inlet_fan,SPLp_inlet_fan,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
            PNLT_discharge_fan = getPNLT(nfreq,n_times,fband,PNL_discharge_fan,SPLp_discharge_fan,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
        elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
            PNLT_inlet_lpc = getPNLT(nfreq,n_times,fband,PNL_inlet_lpc,SPLp_inlet_lpc,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
        elseif( trim(modules(i)).eq.'Lpt' ) then
            PNLT_lpt = getPNLT(nfreq,n_times,fband,PNL_lpt,SPLp_lpt,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
        elseif( trim(modules(i)).eq.'Comb' ) then
            PNLT_Comb = getPNLT(nfreq,n_times,fband,PNL_comb,SPLp_comb,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
        elseif( trim(modules(i)).eq.'cold_nozzle' ) then
            PNLT_caj = getPNLT(nfreq,n_times,fband,PNL_caj,SPLp_caj,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
        else
            ! no further modules implemented
        end if
     end do
     
     PNLT_airfrm = getPNLT(nfreq,n_times,fband,PNL_airfrm,SPLp_airfrm,vec1,wrk1,wrk2,wrk3,wrk4plus,wrk5,wrk6,wrk7)
     
     open (unit = 29, file = "C:\Users\fatbah\Documents\Approach_comb_performance\PNLTout.txt")
    
    do i=1, n_times
       write (29,5) i, PNLT_inlet_fan(i), PNLT_discharge_fan(i), PNLT_inlet_lpc(i),PNLT_lpt(i),PNLT_Comb(i),PNLT_caj(i),PNLT_airfrm(i)
    end do
5   format (i3,f16.10,f16.10,f16.10,f16.10,f16.10,f16.10,f16.10)
    close(29)
     
     deallocate(vec1)
     deallocate(wrk1) ; deallocate(wrk2) ; deallocate(wrk3)
     deallocate(wrk4plus) ; deallocate(wrk5) ; deallocate(wrk6) ; deallocate(wrk7)
     
  contains
     subroutine allocate_variables()
       integer :: i 
       
       allocate(vec1(n_times)) 
       allocate(wrk1(nfreq,n_times)) ; allocate(wrk2(nfreq,n_times)) ; allocate(wrk3(nfreq,n_times))
       allocate(wrk4plus(nfreq+1,n_times)) ; allocate(wrk5(nfreq,n_times)) ; allocate(wrk6(nfreq,n_times)) ; allocate(wrk7(nfreq,n_times))
      
       do i = 1, n_modules ! evaluate all noise sources in this engine model
          if( trim(modules(i)).eq.'Fan' ) then
              allocate(PNLT_inlet_fan(n_times)) ; allocate(PNLT_discharge_fan(n_times))
          elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
              allocate(PNLT_inlet_lpc(n_times)) 
          elseif( trim(modules(i)).eq.'Lpt' ) then
              allocate(PNLT_lpt(n_times)) 
          elseif( trim(modules(i)).eq.'Comb' ) then
              allocate(PNLT_comb(n_times)) 
          elseif( trim(modules(i)).eq.'cold_nozzle' ) then
              allocate(PNLT_caj(n_times)) 
          else
              ! no further modules implemented
          end if
       end do
       
       allocate(PNLT_airfrm(n_times)) 
     
     end subroutine allocate_variables
  end subroutine compute_PNLT
  !
  !
  !
  subroutine compute_noise_sources(gen_noise_source_matr,operatingPoint)
     logical, intent(in) :: gen_noise_source_matr
     character(len=*), intent(in) :: operatingPoint
     !
     integer, parameter :: pnt = 1 ! pointer along flight path
     integer :: i, j
 
     call allocate_module_arrays(n_traj_pts) ! allocate compooutput_unitnent storage arrays 

     call set_vector(zero,180.0_rp,5.0_rp,nthet,theta)     
     call set_frequencies(nb,nfreq,real(fmin,rp),real(fmax,rp),fband,f,freq)
     

! evaluate component noise models
     do i = 1, n_traj_pts ! for all points along trajectory - evaluate noise sources
       do j = 1, n_modules ! evaluate all noise sources in this engine model
          if( trim(modules(j)).eq.'Fan' ) then
            call calcFanandCompressor(i,'Fan',operatingPoint,MtipD_fan,N_rotors_fan,N_stators_fan,rss_fan,Mtip_fan(i),Mu_fan(i),dt_fan(i),xnl_fan(i),g1_fan(i), trim(modules(j)),nfreq,nthet,theta,fband,f,freq, & 
               &  prms_inlet_tone_fan(:,:,i), prms_discharge_tone_fan(:,:,i), prms_inlet_broadband_fan(:,:,i), prms_discharge_broadband_fan(:,:,i), prms_inlet_combination_fan(:,:,i), &  
               &  prms_inlet_fan(:,:,i),prms_discharge_fan(:,:,i)) ! create prms_inlet_fan, prms_discharge_fan
          elseif( trim(modules(j)).eq.'Ipc' .or. trim(modules(j)).eq.'Lpc' ) then
            call calcFanandCompressor(i,trim(modules(j)),operatingPoint,MtipD_lpc,N_rotors_lpc,N_stators_lpc,rss_lpc,Mtip_lpc(i),Mu_lpc(i),dt_lpc(i),xnl_lpc(i),g1_lpc(i),trim(modules(j)),nfreq,nthet,theta,fband,f,freq, & 
               & prms_inlet_tone_Lpc(:,:,i), prms_discharge_tone_Lpc(:,:,i), prms_inlet_broadband_Lpc(:,:,i), prms_discharge_broadband_Lpc(:,:,i), prms_inlet_combination_Lpc(:,:,i), &  
               & prms_inlet_Lpc(:,:,i),prms_discharge_Lpc(:,:,i)) ! create prms_inlet_Lpc, prms_discharge_Lpc
          elseif( trim(modules(j)).eq.'Lpt' ) then
            call calcTurbine(i,NOPRED_test,N_rotors_lpt,n_stages_lpt,V_ref_lpt,c_ref_lpt,Vtr_lpt(i),Texit_lpt(i),xnl_lpt(i),mcore_lpt(i),Cax_lpt(i),Ma(i),Ta(i),clgr(i),alpha(i),SRS_lpt,K_lpt,& 
               & nfreq,nthet,theta,fband,f,freq,prms_Lpt(:,:,i)) ! 
          elseif( trim(modules(j)).eq.'Comb' ) then
            call calcComb(i,NOPRED_test,type_comb, Nf_comb, R0_comb, Aec_comb, De_comb, Dh_comb, Lc_comb, h_comb, Nf_comb, Nfmax_comb, pattern_comb, get_p_ambient(y(i)), & 
               & P3_comb(i),P4_comb(i),P7_comb(i),ta(i),T3_comb(i),T4_comb(i),T5_comb(i),W3_comb(i),nfreq,nthet,theta,fband,f,freq,prms_Comb(:,:,i)) 
          elseif( trim(modules(j)).eq.'cold_nozzle' ) then
            call calcCaj(i,dmdt_1_caj(i),dmdt_2_caj(i),v_1_caj(i),v_2_caj(i),T_1_caj(i),T_2_caj(i),gamma_gas,gamma_air,Ta(i),get_p_ambient(y(i)),A_core_caj, A_bypass_caj, r_const_caj, & 
               & nfreq,nthet,fband,theta,prms_Caj(:,:,i))
          else
              ! at the moment only the noise sources above are implemented
          end if
       end do
       
       ! there must always be an aircraft with an engine....
       call calcAirfrm(i,operatingPoint,NOPRED_test,NoAil_airfrm,NoInb_airfrm,NoOut_airfrm,NoSlat_airfrm,N_wheelsm_airfrm,N_wheelsn_airfrm,& 
          & width_wheelm_airfrm,width_wheeln_airfrm,d_wheelm_airfrm,d_wheeln_airfrm,spann_aileron_airfrm,lscale_aileron_airfrm,& 
          & chord_outb_flap_airfrm,circulation_flap_airfrm,v_spanwise_flap_airfrm,r_const_airfrm,nlsm_airfrm,length_strutm_airfrm,&
          & nlsn_airfrm,length_strutn_airfrm, d_strutm_airfrm, d_strutn_airfrm, geometric_strutm_airfrm, geometric_strutn_airfrm, total_weight_airfrm,& 
          & Ta(i),get_p_ambient(y(i)),Ma(i),Va(i),psi_airfrm(i), Cl_outb_flap_airfrm(i), defl_flap_airfrm(i), aoa_ail_airfrm(i), aoa_slat_airfrm(i), & 
          & nfreq,nthet,theta,fband,f,freq,prms_Airfrm(:,:,i)) 
       
       !
     end do 
     deallocate(prms_discharge_Lpc) ! not used for LPC

!     call save_3D_matrix(nfreq,nthet,n_traj_pts,prms_Airfrm,'.m') ; stop        
!     call save_matrix(nfreq,nthet,prms_inlet_Lpc(:,:,1),'.m') ; stop        

! include multiple engines
     do i = 1, n_modules 
        if( trim(modules(i)).eq.'Fan' ) then
           call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_inlet_fan)
           call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_discharge_fan)
        elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
           call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_inlet_Lpc)
        elseif( trim(modules(i)).eq.'Lpt' ) then
           call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_Lpt)
        elseif( trim(modules(i)).eq.'Comb' ) then
           call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_Comb)
        elseif( trim(modules(i)).eq.'cold_nozzle' ) then
           call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_caj)
        else
            ! at the moment only the fan module is implemented
        end if
     end do
     !
     call include_multiple_engines(no_engines,nfreq,nthet,n_traj_pts,prms_Airfrm)
     
     if( gen_noise_source_matr ) call gen_noise_source_matr_subr(operatingPoint,nfreq,nthet,n_traj_pts,p0,prms_inlet_fan,& 
        & prms_inlet_tone_fan, prms_discharge_tone_fan, prms_inlet_broadband_fan, prms_discharge_broadband_fan, prms_inlet_combination_fan, & 
        & prms_inlet_tone_Lpc, prms_inlet_broadband_Lpc, prms_inlet_combination_Lpc, & 
        & prms_discharge_fan,prms_inlet_Lpc,prms_Lpt,prms_Comb,prms_caj,prms_Airfrm)
     
     ! run tests after completion of calculation
     if(testfunc) then
        call test_functions()    
     endif
     !
  contains
! nfreq and nthet are integer parameters set in choice_physics. Only ntp needs to be variable since it depends on the trajectory which is user input
    subroutine allocate_module_arrays(ntp)
      integer, intent(in) :: ntp   ! number of trajectory points. 

      allocate(prms_inlet_fan(nfreq,nthet,ntp)) ; allocate(prms_discharge_fan(nfreq,nthet,ntp)) 
      allocate(prms_inlet_lpc(nfreq,nthet,ntp)) ; allocate(prms_discharge_lpc(nfreq,nthet,ntp)) 

      allocate(prms_inlet_tone_fan(nfreq,nthet,ntp)) ; allocate(prms_discharge_tone_fan(nfreq,nthet,ntp)) 
      allocate(prms_inlet_broadband_fan(nfreq,nthet,ntp)) ; allocate(prms_discharge_broadband_fan(nfreq,nthet,ntp)) 
      allocate(prms_inlet_combination_fan(nfreq,nthet,ntp)) 

      allocate(prms_inlet_tone_Lpc(nfreq,nthet,ntp)) ; allocate(prms_discharge_tone_Lpc(nfreq,nthet,ntp)) 
      allocate(prms_inlet_broadband_Lpc(nfreq,nthet,ntp)) ; allocate(prms_discharge_broadband_Lpc(nfreq,nthet,ntp)) 
      allocate(prms_inlet_combination_Lpc(nfreq,nthet,ntp)) 

      allocate(prms_Lpt(nfreq,nthet,ntp)) 
      allocate(prms_Comb(nfreq,nthet,ntp)) 
      allocate(prms_Airfrm(nfreq,nthet,ntp)) 
      allocate(prms_Caj(nfreq,nthet,ntp)) 
      
    end subroutine allocate_module_arrays
  end subroutine compute_noise_sources
  !
  !
  !
  subroutine deallocate_module_arrays()

      deallocate(prms_inlet_fan) 
      deallocate(prms_discharge_fan) 
      
      deallocate(prms_inlet_tone_fan) ; deallocate(prms_discharge_tone_fan) 
      deallocate(prms_inlet_broadband_fan) ; deallocate(prms_discharge_broadband_fan) 
      deallocate(prms_inlet_combination_fan) 
      
      deallocate(prms_inlet_tone_Lpc) ; deallocate(prms_discharge_tone_Lpc) 
      deallocate(prms_inlet_broadband_Lpc) ; deallocate(prms_discharge_broadband_Lpc) 
      deallocate(prms_inlet_combination_Lpc) 

      deallocate(prms_inlet_lpc) 
      if( allocated(prms_discharge_lpc) ) deallocate(prms_discharge_lpc) 
      deallocate(prms_Lpt) 
      deallocate(prms_Comb) 
      deallocate(prms_Airfrm) 
      deallocate(prms_Caj) 
      
!
      deallocate(prmsi_inlet_fan) 
      deallocate(prmsi_discharge_fan) 
      deallocate(SPLi_inlet_fan) 
      deallocate(SPLi_discharge_fan) 

      deallocate(prmsi_inlet_lpc) 
      deallocate(SPLi_inlet_lpc) 
      
      deallocate(prmsi_lpt) 
      deallocate(SPLi_lpt) 

      deallocate(prmsi_comb) 
      deallocate(SPLi_comb) 

      deallocate(prmsi_caj) 
      deallocate(SPLi_caj) 
      deallocate(prmsi_airfrm)  
      deallocate(SPLi_airfrm) 
      
      deallocate(xsii_alpha)
      deallocate(fobs) 
      deallocate(atm_absorption) 

      
      deallocate(SPLp_inlet_fan)
      deallocate(prmsp_inlet_fan) 
      deallocate(SPLp_discharge_fan)
      deallocate(prmsp_discharge_fan)  
      deallocate(SPLp_inlet_lpc)
      deallocate(prmsp_inlet_lpc) 
      deallocate(SPLp_lpt)
      deallocate(prmsp_lpt) 
      deallocate(SPLp_comb)  
      deallocate(prmsp_comb) 
      deallocate(SPLp_caj)    
      deallocate(prmsp_caj) 
    
      deallocate(SPLp_airfrm)
      deallocate(prmsp_airfrm) 

      deallocate(PNL_inlet_fan) 
      deallocate(PNL_discharge_fan)
      deallocate(PNL_inlet_lpc) 
      deallocate(PNL_lpt) 
      deallocate(PNL_comb) 
      deallocate(PNL_caj) 
      deallocate(PNL_airfrm) 
      
      deallocate(PNLT_inlet_fan)
      deallocate(PNLT_discharge_fan)
      deallocate(PNLT_inlet_lpc) 
      deallocate(PNLT_lpt) 
      deallocate(PNLT_comb) 
      deallocate(PNLT_caj) 
      deallocate(PNLT_airfrm) 
      
  end subroutine deallocate_module_arrays
  !
  !
  !
  subroutine include_multiple_engines(no_engines,nfreq,nthet,ntp,prms)
    integer, intent(in) :: no_engines
    integer, intent(in) :: nthet, nfreq, ntp
    real(kind=rp), dimension(nfreq,nthet,ntp) :: prms
    integer :: i, j, k
    
    do i = 1, nfreq
      do j = 1, nthet
         do k = 1, ntp
            prms(i,j,k) = sqrt(real(no_engines,rp)*prms(i,j,k)**2)
         end do 
      end do
    end do
      
  end subroutine include_multiple_engines
  !
  !
  !
  subroutine compute_flight_effects(iptr,use_ground_refl,spherical_spr,atm_atten,operatingPoint)
    integer, intent(in) :: iptr
    logical, intent(in) :: use_ground_refl, spherical_spr, atm_atten
    character(len=*), intent(in) :: operatingPoint
    logical, parameter :: save_SPLp_data = .false.
    integer :: i 
    
!    call save_3D_matrix(nfreq,nthet,n_times,SPLi_inlet_fan,'.m') ; stop
    call allocate_arrays()
    
    xsii_alpha = get_theta(n_times,x_source,xsii) ! add aircraft alpha to xsii
    
    fobs(1:nfreq,1:n_times) = getDopplerShift(nfreq,n_times,fband,xsii,Mai)

    do i = 1, n_times 
      atm_absorption(i,1:nfreq) = get_atm_abs(nfreq,Tai(i),RH,fobs(:,i))
    end do 
    
    do i = 1, n_modules
      if( trim(modules(i)).eq.'Fan' ) then
        call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption, & 
           & SPLi_inlet_fan,SPLp_inlet_fan,prmsp_inlet_fan) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
        call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption, & 
           & SPLi_discharge_fan,SPLp_discharge_fan,prmsp_discharge_fan) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
      elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
        call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption, & 
           & SPLi_inlet_lpc,SPLp_inlet_lpc,prmsp_inlet_lpc) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
      elseif( trim(modules(i)).eq.'Lpt' ) then
        call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption, & 
           & SPLi_lpt,SPLp_lpt,prmsp_lpt) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
      elseif( trim(modules(i)).eq.'Comb' ) then
        call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption,& 
           & SPLi_comb,SPLp_comb,prmsp_comb) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
      elseif( trim(modules(i)).eq.'cold_nozzle' ) then
        call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption, & 
           & SPLi_caj,SPLp_caj,prmsp_caj) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
      end if 
    end do 

    call flightEffects(iptr,nfreq,nthet,n_times,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x_source,y_source,r1,Tai,fobs,atm_absorption, & 
           & SPLi_airfrm,SPLp_airfrm,prmsp_airfrm) ! propagate to microphone (select directivity, propagation & doppler, ground reflection)
    
    if( save_SPLp_data ) call gen_SPLp_matr_subr(operatingPoint,nfreq,n_times,prmsp_inlet_fan,& 
        & prmsp_discharge_fan,prmsp_inlet_lpc,prmsp_lpt,prmsp_comb,prmsp_caj,prmsp_airfrm)

  
  contains
    subroutine allocate_arrays()
      integer :: i 
      
      allocate(xsii_alpha(n_times))
      allocate(fobs(nfreq,n_times)) 
      allocate(atm_absorption(n_times,nfreq)) 

      do i = 1, n_modules
        if( trim(modules(i)).eq.'Fan' ) then
          allocate(SPLp_inlet_fan(nfreq,n_times))     ; allocate(prmsp_inlet_fan(nfreq,n_times)) 
          allocate(SPLp_discharge_fan(nfreq,n_times)) ; allocate(prmsp_discharge_fan(nfreq,n_times))  
      elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
          allocate(SPLp_inlet_lpc(nfreq,n_times))     ; allocate(prmsp_inlet_lpc(nfreq,n_times)) 
      elseif( trim(modules(i)).eq.'Lpt' ) then
          allocate(SPLp_lpt(nfreq,n_times))     ; allocate(prmsp_lpt(nfreq,n_times)) 
      elseif( trim(modules(i)).eq.'Comb' ) then
          allocate(SPLp_comb(nfreq,n_times))     ; allocate(prmsp_comb(nfreq,n_times)) 
      elseif( trim(modules(i)).eq.'cold_nozzle' ) then
          allocate(SPLp_caj(nfreq,n_times))     ; allocate(prmsp_caj(nfreq,n_times)) 
      end if 
    end do 
    
    allocate(SPLp_airfrm(nfreq,n_times))     ; allocate(prmsp_airfrm(nfreq,n_times)) 
    
    end subroutine allocate_arrays
  end subroutine compute_flight_effects
  !
  !
  !
  subroutine gen_SPLp_matr_subr(operatingPoint,nfreq,n_times,inlet_fan,discharge_fan,inlet_Lpc,Lpt,Comb,Caj,Airfrm)
    character(len=*), intent(in) :: operatingPoint
    integer, intent(in) :: nfreq, n_times
    real(kind=rp), dimension(nfreq,n_times), intent(in) :: inlet_fan,discharge_fan,inlet_Lpc,Lpt,Comb,Caj,Airfrm

    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_fanInlet',nfreq,n_times,inlet_fan)
    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_fanDischarge',nfreq,n_times,discharge_fan)
    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_lpcInlet',nfreq,n_times,inlet_Lpc)
    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_Lpt',nfreq,n_times,Lpt)
    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_Comb',nfreq,n_times,Comb)
    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_Caj',nfreq,n_times,Caj) 
    call save_matrix('choiceOutput/'//trim(operatingPoint)//'_Airfrm',nfreq,n_times,Airfrm)
  
  end subroutine gen_SPLp_matr_subr 
  !
  !
  !
  function getDopplerShift(nfr,nti,fband,xsii,Mai)
    integer, intent(in) :: nfr,nti
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), intent(in), dimension(nti) :: xsii, Mai
    real(kind=rp), dimension(nfr,nti) :: getDopplerShift
    !
    integer :: i, j
    
! Doppler shift the frequency vector 
    do i = 1, nfr
       do j = 1, nti
          getDopplerShift(i,j) = fband(i)/(one-Mai(j)*cos(xsii(j)))
       end do
    end do
    
  end function getDopplerShift
  !
  !
  !
  subroutine interpolate_to_t_source()
     integer :: i, j, k, l

     ! [xsi_plane_inter_TO]=interpolation(xsi_plane,t_flight_step,t_source_TO,1);
     
     ! interpolation is right, thinking is right but xsi_plane is wrong. It needs to contain angles between observer and flightpath.
     
     if( allocated(xsii) )     deallocate(xsii) ;    allocate(xsii(n_times)) 
     if( allocated(Mai)  )      deallocate(Mai) ;     allocate(Mai(n_times)) 
     if( allocated(Tai)  )      deallocate(Tai) ;     allocate(Tai(n_times)) 
     if( allocated(alphai) ) deallocate(alphai) ;  allocate(alphai(n_times)) 
     
     xsii(1) = xsi(1) ; Mai(1) = Ma(1) ; Tai(1) = Ta(1) ; alphai(1)  = alpha(1) 
     do i = 2, n_times 
       xsii(i) = get_vector_interpolation(t_source(i),nfreq,nthet,n_traj_pts,time,xsi) ! create xi for the t_source times rather than times
        Mai(i) = get_vector_interpolation(t_source(i),nfreq,nthet,n_traj_pts,time, Ma) ! create Ma for the t_source times rather than times
        Tai(i) = get_vector_interpolation(t_source(i),nfreq,nthet,n_traj_pts,time, Ta) ! create Ta for the t_source times rather than times
        alphai(i) = get_vector_interpolation(t_source(i),nfreq,nthet,n_traj_pts,time,alpha) ! create Ta for the t_source times rather than times
     end do
     
     ! interpolate 
     do i = 1, n_modules
       if( trim(modules(i)).eq.'Fan' ) then 
           ! 
           allocate(prmsi_inlet_fan(nfreq,nthet,n_times)) ; allocate(prmsi_discharge_fan(nfreq,nthet,n_times)) 
           allocate(SPLi_inlet_fan(nfreq,nthet,n_times)) ; SPLi_inlet_fan = zero
           allocate(SPLi_discharge_fan(nfreq,nthet,n_times)) ; SPLi_discharge_fan = zero 
           !
           call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_inlet_fan,prmsi_inlet_fan)
           call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_discharge_fan,prmsi_discharge_fan)
           !
           call compute_SPLi(nfreq,nthet,n_times,prmsi_inlet_fan,SPLi_inlet_fan)
           call compute_SPLi(nfreq,nthet,n_times,prmsi_discharge_fan,SPLi_discharge_fan)
           !
       elseif( trim(modules(i)).eq.'Ipc' .or. trim(modules(i)).eq.'Lpc' ) then
           !
           allocate(prmsi_inlet_lpc(nfreq,nthet,n_times)) 
           allocate(SPLi_inlet_lpc(nfreq,nthet,n_times)) 
           !
           call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_inlet_lpc,prmsi_inlet_lpc)
           !
           call compute_SPLi(nfreq,nthet,n_times,prmsi_inlet_lpc,SPLi_inlet_lpc)
           !
       elseif( trim(modules(i)).eq.'Lpt' ) then
           !
           allocate(prmsi_lpt(nfreq,nthet,n_times)) 
           allocate(SPLi_lpt(nfreq,nthet,n_times)) 
           !
           call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_lpt,prmsi_lpt)
           !
           call compute_SPLi(nfreq,nthet,n_times,prmsi_lpt,SPLi_lpt)
           !
       elseif( trim(modules(i)).eq.'Comb' ) then
           !
           allocate(prmsi_comb(nfreq,nthet,n_times)) 
           allocate(SPLi_comb(nfreq,nthet,n_times)) 
           !
           call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_comb,prmsi_comb)
           !
           call compute_SPLi(nfreq,nthet,n_times,prmsi_comb,SPLi_comb)
           !
       elseif( trim(modules(i)).eq.'cold_nozzle' ) then
           !
           allocate(prmsi_caj(nfreq,nthet,n_times)) 
           allocate(SPLi_caj(nfreq,nthet,n_times)) 
           !
           call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_caj,prmsi_caj)
           !
           call compute_SPLi(nfreq,nthet,n_times,prmsi_caj,SPLi_caj)
           !
       else
           ! more modules needed....
       end if
     end do 
     
     allocate(prmsi_airfrm(nfreq,nthet,n_times)) ; allocate(SPLi_airfrm(nfreq,nthet,n_times)) 
     call source_interpolation(nfreq,nthet,n_traj_pts,n_times,prms_airfrm,prmsi_airfrm)
     !
     call compute_SPLi(nfreq,nthet,n_times,prmsi_airfrm,SPLi_airfrm)
     
!    call save_3D_matrix(nfreq,nthet,n_times,SPLi_lpt,'.m') ; stop
!     call save_matrix(nfreq,nthet,SPLi_inlet_fan(:,:,1),'.m') ; stop
     
  end subroutine interpolate_to_t_source
  !
  !
  !
  subroutine compute_SPLi(nfr,nth,nti,prmsi,SPLi)
    integer, intent(in) :: nfr, nth, nti
    real(kind=rp), intent(in), dimension(nfr,nth,nti) :: prmsi
    real(kind=rp), intent(out), dimension(nfr,nth,nti) :: SPLi
    !
    integer :: i, j, k

    do i = 1, nfr 
       do j = 1, nth
          do k = 1, nti
             !
             SPLi(i,j,k) = 20.0_rp*log10(prmsi(i,j,k)/p0) 
             !
          end do 
       end do
    end do
  
  end subroutine compute_SPLi
  !
  !
  !
  subroutine source_interpolation(nfr,nth,ntp,nti,prms,prmsi)
    integer, intent(in) :: nfr,nth,ntp,nti
    real(kind=rp), intent(in), dimension(nfr,nth,ntp) :: prms 
    real(kind=rp), intent(out), dimension(nfr,nth,nti) :: prmsi
    !
    integer :: i 

! first observer and source pair is set to be the same
    prmsi(1:nfr,1:nth,1) = prms(1:nfreq,1:nthet,1)

    !
    do i = 2, nti
      prmsi(1:nfreq,1:nthet,i) = get_matrix_interpolation(t_source(i),nfr,nth,ntp,time,prms)
    end do

!    call save_3D_matrix(nfreq,nthet,ntp,prms,'.m') ; stop
    
  end subroutine source_interpolation
  !
  ! set_performance_choice: transfer data to choice_physics, variable by variable.
  ! The whole input file is transferred as argument. 
  !
  subroutine set_performance_choice(ptr,trajPerf,point)
    integer, intent(in) :: ptr
    logical, intent(in) :: trajPerf
    character(len=*), intent(in) :: point ! Take-off, top-of-climb, cruise, approach, cutback etc...
    !
    integer, parameter :: maxNoColumns = 8 ! largest number of columns any module has in its performance file.
    character(len=32), dimension(maxNoColumns) :: strVec
    real(kind=rp), allocatable, dimension(:,:) :: storageMat
    character(len=32), allocatable, dimension(:,:) :: charStorageMat ! used to pre-processed "starred" input coming from mission analysis
    logical :: openRotor
    integer :: i, slen
    
    openRotor = getIsOpenRotor()
    if( openRotor ) stop 'Model is open rotor - noise modelling not implemented yet'
    
    if( trajPerf ) allocate(storageMat(n_traj_pts,maxNoColumns))
    
    do i = 1, n_modules
       if( trim(modules(i)).eq.'Fan' ) call setFan()
       if( trim(modules(i)).eq.'Ipc' .or.  trim(modules(i)).eq.'Lpc' ) call setLpc(openRotor)
       if( trim(modules(i)).eq.'Lpt' ) call setLpt()
       if( trim(modules(i)).eq.'Comb' ) call setComb()
       if( trim(modules(i)).eq.'cold_nozzle' ) call coAxialJet()
    end do 
    
    call setAirfrm(ptr)
    
    if( trajPerf ) deallocate(storageMat)
  
  contains 
     function getIsOpenRotor()
       logical :: getIsOpenRotor
       integer :: i, slen

       getIsOpenRotor = .false.
       
       do i = 1, n_modules
         slen = len(trim(modules(i))) 
         if( trim(modules(i)(1:slen)).eq.'Prop' ) then
             getIsOpenRotor = .true.
             exit
         end if
       end do 
       
     end function getIsOpenRotor
     !
     subroutine coAxialJet()
       real(kind=rp) :: W8, W18, M8, M18
       real(kind=rp) :: T8, T18, C8, C18
       integer :: i
       !
       if( allocated(dmdt_1_caj) ) deallocate(dmdt_1_caj) ; allocate(dmdt_1_caj(n_traj_pts)) 
       if( allocated(dmdt_2_caj) ) deallocate(dmdt_2_caj) ; allocate(dmdt_2_caj(n_traj_pts)) 
       if( allocated(v_1_caj) ) deallocate(v_1_caj) ; allocate(v_1_caj(n_traj_pts)) 
       if( allocated(v_2_caj) ) deallocate(v_2_caj) ; allocate(v_2_caj(n_traj_pts)) 
       if( allocated(T_1_caj) ) deallocate(T_1_caj) ; allocate(T_1_caj(n_traj_pts)) 
       if( allocated(T_2_caj) ) deallocate(T_2_caj) ; allocate(T_2_caj(n_traj_pts)) 

       if( trajPerf ) then
           call loadStorageMat(trim(point)//'_coAxialJet_performance.txt',n_traj_pts,6,strVec,storageMat)
           do i = 1, n_traj_pts
               dmdt_1_caj(i)  = storageMat(i,1) ; dmdt_2_caj(i) = storageMat(i,2) ; v_1_caj(i) = storageMat(i,3)
               v_2_caj(i) = storageMat(i,4) ; T_1_caj(i) = storageMat(i,5) ; T_2_caj(i) = storageMat(i,6) 
           end do
       else
          ! 2 is cold. 
          W8 = get_variable(point,'W8:') ; W18 = get_variable(point,'W18:') 
          M8 = get_variable(point,'M8:') ; M18 = get_variable(point,'M18:') 
          T8 = get_variable(point,'T8:') ; T18 = get_variable(point,'T18:') 
          !
          C18 = getVel(gamma_air,T18,M18)
          C8 = getVel(gamma_gas,T8,M8)
          !
          do i = 1, n_traj_pts
             !
             dmdt_1_caj(i)  = W8 ; dmdt_2_caj(i) = W18 ; v_1_caj(i) = C8
             v_2_caj(i) = C18 ; T_1_caj(i) = T8 ; T_2_caj(i) = T18 
             !
          end do
      end if
     
     end subroutine coAxialJet
     !
     subroutine setAirfrm(ptr)
       integer, intent(in) :: ptr
       integer :: i
       !
       if( allocated(psi_airfrm) ) deallocate(psi_airfrm) ; allocate(psi_airfrm(n_traj_pts)) 
       if( allocated(Cl_outb_flap_airfrm) ) deallocate(Cl_outb_flap_airfrm) ; allocate(Cl_outb_flap_airfrm(n_traj_pts)) 
       if( allocated(defl_flap_airfrm) ) deallocate(defl_flap_airfrm) ; allocate(defl_flap_airfrm(n_traj_pts)) 
       if( allocated(aoa_ail_airfrm) ) deallocate(aoa_ail_airfrm) ; allocate(aoa_ail_airfrm(n_traj_pts)) 
       if( allocated(aoa_slat_airfrm) ) deallocate(aoa_slat_airfrm) ; allocate(aoa_slat_airfrm(n_traj_pts)) 

       if( trajPerf ) then
           call loadStorageMat(trim(point)//'_airfrm_performance.txt',n_traj_pts,5,strVec,storageMat)
           do i = 1, n_traj_pts
               psi_airfrm(i)  = storageMat(i,1) ; Cl_outb_flap_airfrm(i) = storageMat(i,2) ; defl_flap_airfrm(i) = storageMat(i,3)
               aoa_ail_airfrm(i) = storageMat(i,4) ; aoa_slat_airfrm(i) = storageMat(i,5) 
           end do
       else
          ! 
          do i = 1, n_traj_pts
             !
             psi_airfrm(i) = psi_airfrm_vec(ptr) 
             Cl_outb_flap_airfrm(i) = Cl_outb_flap_airfrm_vec(ptr) 
             defl_flap_airfrm(i) = defl_flap_airfrm_vec(ptr) 
             aoa_ail_airfrm(i) = aoa_ail_airfrm_vec(ptr) 
             aoa_slat_airfrm(i) = aoa_slat_airfrm_vec(ptr) 
             !
          end do
          !
      end if
       
     end subroutine setAirfrm
     !
     subroutine setLpt()
       integer :: i
       real(kind=rp) :: xnl, t5, w5, Cax, p5

       if( allocated(Vtr_lpt) )     deallocate(Vtr_lpt) ;   allocate(Vtr_lpt(n_traj_pts)) 
       if( allocated(Texit_lpt) ) deallocate(Texit_lpt) ; allocate(Texit_lpt(n_traj_pts)) 
       if( allocated(xnl_lpt) ) deallocate(xnl_lpt) ; allocate(xnl_lpt(n_traj_pts)) ! new for every operating point 
       if( allocated(mcore_lpt) ) deallocate(mcore_lpt) ; allocate(mcore_lpt(n_traj_pts)) ! new for every operating point 
       if( allocated(Cax_lpt) ) deallocate(Cax_lpt) ; allocate(Cax_lpt(n_traj_pts)) ! new for every operating point 

       if( trajPerf ) then
           if( containsStars(trim(point)//'_lpt_performance.txt') ) then
               write(*,*) 'Discovered raw trajectory mission file (lpt) - computing *** and ** from other data and re-writing file (this is done only once)'
               call preProcessLptFile(trim(point)//'_lpt_performance.txt',De_lpt,Ae_lpt)
           end if
           !
           call loadStorageMat(trim(point)//'_lpt_performance.txt',n_traj_pts,5,strVec,storageMat)
           do i = 1, n_traj_pts
               Vtr_lpt(i) = storageMat(i,1) ; Texit_lpt(i) = storageMat(i,2) ; xnl_lpt(i) = storageMat(i,3) ; mcore_lpt(i) = storageMat(i,4) ; Cax_lpt(i) = storageMat(i,5) 
           end do
       else
          ! use performanceResults.txt file and additional data to estimate the trajactory variables.
          xnl = get_variable(point,'NL:') ; t5 = get_variable(point,'T5:') ; w5 = get_variable(point,'W5:') ; p5 = get_variable(point,'P5:') 
          Cax = get_Cax(gamma_gas,p5,t5,w5,Ae_lpt)
          !
          do i = 1, n_traj_pts
             !
             Vtr_lpt(i) = (pi_num*De_lpt)*xnl
             Texit_lpt(i) = t5
             xnl_lpt(i) = xnl
             mcore_lpt(i) = w5
             Cax_lpt(i) = Cax
             !
          end do
          !
      end if
       
       
     end subroutine setLpt
     !
     subroutine setComb()
       real(kind=rp) :: P3, P4, P7, T3, T4, T5, W3
       integer :: i
       
       if( allocated(P3_comb) ) deallocate(P3_comb) ; allocate(P3_comb(n_traj_pts)) 
       if( allocated(P4_comb) ) deallocate(P4_comb) ; allocate(P4_comb(n_traj_pts)) 
       if( allocated(P7_comb) ) deallocate(P7_comb) ; allocate(P7_comb(n_traj_pts)) 
       if( allocated(T3_comb) ) deallocate(T3_comb) ; allocate(T3_comb(n_traj_pts)) 
       if( allocated(T4_comb) ) deallocate(T4_comb) ; allocate(T4_comb(n_traj_pts)) 
       if( allocated(T5_comb) ) deallocate(T5_comb) ; allocate(T5_comb(n_traj_pts)) 
       if( allocated(W3_comb) ) deallocate(W3_comb) ; allocate(W3_comb(n_traj_pts)) 

       if( trajPerf ) then
           call loadStorageMat(trim(point)//'_comb_performance.txt',n_traj_pts,7,strVec,storageMat)
           do i = 1, n_traj_pts
               P3_comb(i) = storageMat(i,1) ; P4_comb(i) = storageMat(i,2) ; P7_comb(i) = storageMat(i,3) ; T3_comb(i) = storageMat(i,4)
               T4_comb(i) = storageMat(i,5) ; T5_comb(i) = storageMat(i,6) ; W3_comb(i) = storageMat(i,7) 
           end do
       else
          ! use performanceResults.txt file and additional data to estimate the trajactory variables.
          P3 = get_variable(point,'P3:') ; P4= get_variable(point,'P4:') ;  P7= get_variable(point,'P5:')
          T3 = get_variable(point,'T3:') ; T4= get_variable(point,'T4:') ;  T5= get_variable(point,'T5:') ;  W3= get_variable(point,'W3:')  
          !
          do i = 1, n_traj_pts
             P3_comb(i) = P3 ; P4_comb(i) = P4 ; P7_comb(i) = P7
             T3_comb(i) = T3 ; T4_comb(i) = T4 ; T5_comb(i) = T5 ; W3_comb(i) = W3
          end do
          !
      end if
       
     end subroutine setComb
     !
     subroutine setFan()
       integer :: i
       real(kind=rp) :: xnl, g1, t1, dt, Mtip, Mu, p1, Umid
       
! allocate fan variables (re-allocate for every certification point)
       if( allocated(Mtip_fan) ) deallocate(Mtip_fan) ; allocate(Mtip_fan(n_traj_pts)) ! new for every operating point 
       if( allocated(Mu_fan) ) deallocate(Mu_fan) ; allocate(Mu_fan(n_traj_pts)) ! new for every operating point 
       if( allocated(dt_fan) ) deallocate(dt_fan) ; allocate(dt_fan(n_traj_pts)) ! new for every operating point 
       if( allocated(xnl_fan) ) deallocate(xnl_fan) ; allocate(xnl_fan(n_traj_pts)) ! new for every operating point 
       if( allocated(g1_fan) ) deallocate(g1_fan) ; allocate(g1_fan(n_traj_pts)) ! new for every operating point 

       if( trajPerf ) then
           if( containsStars(trim(point)//'_fan_performance.txt') ) then
               write(*,*) 'Discovered raw trajectory mission file (fan) - computing Mtip and Mu from other data and re-writing file (this is done only once)'
               call preProcessFanFile(trim(point)//'_fan_performance.txt',D1_fan,A2_fan)
           end if
           !
           call loadStorageMat(trim(point)//'_fan_performance.txt',n_traj_pts,5,strVec,storageMat)
           do i = 1, n_traj_pts
               Mtip_fan(i) = storageMat(i,1) ; Mu_fan(i) = storageMat(i,2) ; xnl_fan(i) = storageMat(i,3) ; dt_fan(i) = storageMat(i,4) ; g1_fan(i) = storageMat(i,5)
           end do
       else
          ! use performanceResults.txt file and additional data to estimate the trajactory variables.
          xnl = get_variable(point,'NL:') ! derived from the WEICO top-of-climb value using performance factors => it is rps
          g1 = get_variable(point,'W2:') ; t1 = get_variable(point,'T2:') 
          dt = get_variable(point,'T13:') - t1 ; p1 = get_variable(point,'P2:') 
          call setMachNumbers(p1,t1,g1,A2_fan,D1_fan,xnl,gamma_air,Mtip,Mu,Umid)
          !
          do i = 1, n_traj_pts
             xnl_fan(i) = xnl ; dt_fan(i) = dt ; g1_fan(i) = g1
             Mtip_fan(i) = Mtip ; Mu_fan(i) = Mu 
          end do
      end if
                    
     end subroutine setFan
     !
     subroutine setLpc(or)
       logical, intent(in) :: or
       real(kind=rp) :: xnl, xni, g1, t1, dt, Mtip, Mu, p1, Umid, A, A2_lpc
       integer :: i
       
! allocate fan variables (re-allocate for every certification point)
       if( allocated(Mtip_lpc) ) deallocate(Mtip_lpc) ; allocate(Mtip_lpc(n_traj_pts)) ! new for every operating point 
       if( allocated(Mu_lpc) ) deallocate(Mu_lpc) ; allocate(Mu_lpc(n_traj_pts)) ! new for every operating point 
       if( allocated(dt_lpc) ) deallocate(dt_lpc) ; allocate(dt_lpc(n_traj_pts)) ! new for every operating point 
       if( allocated(xnl_lpc) ) deallocate(xnl_lpc) ; allocate(xnl_lpc(n_traj_pts)) ! new for every operating point 
       if( allocated(g1_lpc) ) deallocate(g1_lpc) ; allocate(g1_lpc(n_traj_pts)) ! new for every operating point 

       if( trajPerf ) then
           if( containsStars(trim(point)//'_lpc_performance.txt') ) then
               write(*,*) 'Discovered raw trajectory mission file (lpc) - computing Mtip and Mu from other data and re-writing file (this is done only once)'
               A2_lpc = pi_num*((D1_lpc/two)**2 - (Dh1_lpc/two)**2) 
               call preProcessFanFile(trim(point)//'_lpc_performance.txt',D1_lpc,A2_lpc)
           end if
           call loadStorageMat(trim(point)//'_lpc_performance.txt',n_traj_pts,5,strVec,storageMat)
           do i = 1, n_traj_pts
               Mtip_lpc(i) = storageMat(i,1) ; Mu_lpc(i) = storageMat(i,2) ; xnl_lpc(i) = storageMat(i,3) ; dt_lpc(i) = storageMat(i,4) ;  g1_lpc(i) = storageMat(i,5)
           end do
       else
          ! use performanceResults.txt file and additional data to estimate the trajactory variables.
          if( .not. or ) then
             xnl = get_variable(point,'NL:')            
             xni = xnl*gbx_ratio
          else
             stop 'no open rotor method yet - lpc'
          end if
          
          g1 = get_variable(point,'W23:') ; t1 = get_variable(point,'T23:')  ; p1 = get_variable(point,'P23:') 
          A = (pi_num/four)*(D1_lpc**2-Dh1_lpc**2) 
          call setMachNumbers(p1,t1,g1,A,D1_lpc,xni,gamma_air,Mtip,Mu,Umid)
          !
          if( gbx_ratio.eq.one ) then
             dt = (0.6_rp*Umid**2)/(cp_air*n_stages_LPC)
          else 
             dt = (0.45_rp*Umid**2)/(cp_air*n_stages_LPC)
          end if
          !
          do i = 1, n_traj_pts             
             xnl_lpc(i) = xni ; dt_lpc(i) = dt ; g1_lpc(i) = g1
             Mtip_lpc(i) = Mtip ; Mu_lpc(i) = Mu 
          end do
          ! use performanceResults.txt file and additional data to estimate the trajactory variables.
        end if
                    
     end subroutine setLpc
  end subroutine set_performance_choice
  !
  !
  ! set_rotational_speeds_choice: run prior to handling performance. All rotational speeds are 
  ! re-computed in two steps:
  !            
  !     1. Use rotational speed from the weico aero design point. 
  !     2. Take out all performance rotational speeds and recompute these to ratios, by dividing all rotational speed with the 
  !        aero design point of the performance model. Multiply these with the value from point 1 above. 
  !
  subroutine set_rotational_speeds_choice()
     !
     real(kind=rp) :: N_to, N_sl, N_toc, N_cr, N_appr, N_cutb
 
     if( mpd%n_lines_to.gt.0 )       N_to   = get_variable('Take-off'    ,'NL:') 
     if( mpd%n_lines_sl.gt.0 )       N_sl   = get_variable('Sideline'    ,'NL:') 
     if( mpd%n_lines_toc.gt.0 )      N_toc  = get_variable('Top-of-climb','NL:') 
     if( mpd%n_lines_cr.gt.0 )       N_cr   = get_variable('Cruise'      ,'NL:') 
     if( mpd%n_lines_approach.gt.0 ) N_appr = get_variable('Approach'    ,'NL:') 
     if( mpd%n_lines_cutback.gt.0 )  N_cutb = get_variable('Cutback'     ,'NL:')
     
     if( mpd%n_lines_sl.gt.0 ) N_to = N_sl ! use sideline data instead of take-off if sideline exists.
     
     ! form ratios and use the WEICO value
     N_to = xnlD_fan*(N_to/N_toc) ; N_cr = xnlD_fan*(N_cr/N_toc) ; N_appr = xnlD_fan*(N_appr/N_toc) ; N_cutb = xnlD_fan*(N_cutb/N_toc) ; N_sl = xnlD_fan*(N_sl/N_toc) 
     N_toc = xnlD_fan*one  
 
     ! set back into the mpd structure 
     if( mpd%n_lines_to.gt.0 )       call set_variable('Take-off'    ,'NL:',N_to) 
     if( mpd%n_lines_toc.gt.0 )      call set_variable('Top-of-climb','NL:',N_toc) 
     if( mpd%n_lines_sl.gt.0 )       call set_variable('Sideline'    ,'NL:',N_sl) 
     if( mpd%n_lines_cr.gt.0 )       call set_variable('Cruise'      ,'NL:',N_cr) 
     if( mpd%n_lines_approach.gt.0 ) call set_variable('Approach'    ,'NL:',N_appr) 
     if( mpd%n_lines_cutback.gt.0 )  call set_variable('Cutback'     ,'NL:',N_cutb) 
      
  end subroutine set_rotational_speeds_choice
  !
  !
  !
  function get_variable(point,name)
    character(len=*), intent(in) :: point, name
    real(kind=rp) :: get_variable 
    !
    integer :: slen
    
    slen = len(point)
    
    if( point(1:slen).eq.'Take-off' ) then
      get_variable = get_value(name,mpd%n_lines_to,mpd%perf_file_to(1:mpd%n_lines_to))
    elseif( point(1:slen).eq.'Sideline' ) then
      get_variable = get_value(name,mpd%n_lines_sl,mpd%perf_file_sl(1:mpd%n_lines_to))
    elseif( point(1:slen).eq.'Cutback' ) then
      get_variable = get_value(name,mpd%n_lines_cutback,mpd%perf_file_cutback(1:mpd%n_lines_cutback))
    elseif( point(1:slen).eq.'Approach' ) then
      get_variable = get_value(name,mpd%n_lines_approach,mpd%perf_file_approach(1:mpd%n_lines_approach))
    elseif( point(1:slen).eq.'Cruise' ) then
      get_variable = get_value(name,mpd%n_lines_cr,mpd%perf_file_cr(1:mpd%n_lines_cr))
    elseif( point(1:slen).eq.'Top-of-climb' ) then
      get_variable = get_value(name,mpd%n_lines_toc,mpd%perf_file_toc(1:mpd%n_lines_toc))
    end if
    
  end function get_variable
  !
  !
  !
  function get_string(name,nl,file)
    integer, intent(in) :: nl
    character(len=*), intent(in) :: name
    character(len=cline_len), intent(in), dimension(nl) :: file
    character(len=32) :: get_string
    !
    integer :: i, cp
    logical :: found
    
    found = .false.
    get_string = ' '
    
    do i = 1, nl
        cp = index(file(i),':')
        if( cp.eq.0 ) cycle
        if( file(i)(1:cp).eq.name ) then
            get_string = trim(adjustl(file(i)(cp+1:cline_len)))
            found = .true.
            exit
        end if
    end do 
    
    if( .not.found ) call report_error('name '//trim(name)//' not found in performance file', & 
       & 'get_value','choice_interf')
  
  end function get_string
  !
  !
  !
  function get_value(name,nl,file)
    integer, intent(in) :: nl
    character(len=*), intent(in) :: name
    character(len=cline_len), intent(in), dimension(nl) :: file
    real(kind=rp) :: get_value
    !
    integer :: i, cp
    logical :: found
    
    found = .false.
    
    do i = 1, nl
        cp = index(file(i),':')
        if( cp.eq.0 ) cycle
        if( file(i)(1:cp).eq.name ) then
            get_value = string_number_cast(file(i)(cp+1:cline_len))
            found = .true.
            exit
        end if
    end do 
    
    if( .not.found ) call report_error('name '//trim(name)//' not found in performance file', & 
       & 'get_value','choice_interf')
  
  end function get_value
  !
  !
  !
  subroutine set_variable(point,name,value)
    character(len=*), intent(in) :: point, name
    real(kind=rp), intent(in) :: value
    !
    integer :: slen
    
    slen = len(point)
    
    if( point(1:slen).eq.'Take-off' ) then
      call set_value(name,mpd%n_lines_to,mpd%perf_file_to(1:mpd%n_lines_to),value)
    elseif( point(1:slen).eq.'Sideline' ) then
      call set_value(name,mpd%n_lines_sl,mpd%perf_file_sl(1:mpd%n_lines_sl),value)
    elseif( point(1:slen).eq.'Cutback' ) then
      call set_value(name,mpd%n_lines_cutback,mpd%perf_file_cutback(1:mpd%n_lines_cutback),value)
    elseif( point(1:slen).eq.'Approach' ) then
      call set_value(name,mpd%n_lines_approach,mpd%perf_file_approach(1:mpd%n_lines_approach),value)
    elseif( point(1:slen).eq.'Cruise' ) then
      call set_value(name,mpd%n_lines_cr,mpd%perf_file_cr(1:mpd%n_lines_cr),value)
    elseif( point(1:slen).eq.'Top-of-climb' ) then
      call set_value(name,mpd%n_lines_toc,mpd%perf_file_toc(1:mpd%n_lines_toc),value)
    end if
    
  end subroutine set_variable
  !
  !
  !
  subroutine set_value(name,nl,file,val)
    integer, intent(in) :: nl
    character(len=*), intent(in) :: name
    character(len=cline_len), intent(inout), dimension(nl) :: file
    real(kind=rp), intent(in) :: val
    !
    integer :: i, cp
    logical :: found
    
    found = .false.
    
    do i = 1, nl
        cp = index(file(i),':')
        if( cp.eq.0 ) cycle
        if( file(i)(1:cp).eq.name ) then
            file(i) = name//number_string_cast(val)
            found = .true.
            exit
        end if
    end do 
    
    if( .not.found ) call report_error('name '//trim(name)//' not found in performance file', & 
       & 'set_variable','choice_interf')
  
  end subroutine set_value
  !
  !
  !
  subroutine set_trajectory_subr(ipnt,opPnt)
    integer, intent(in) :: ipnt
    character(len=*), intent(in) :: opPnt
    
! read x, y, Va, alpha from trajectory file (Cutback.txt, Take-off.txt or Approach.txt).
    call setTrajectory(ipnt,opPnt) ! Compute r, clgr and xsi, time
       
    call setVelocity() ! compute ta from y(i) and dtIsa, a(i) from ta(i) and Mach number from Va(i) and a(i)    
    call setTimeVector(ipnt) ! compute "retarded" vectors at source (see manual) 

    
  
  contains 
     subroutine setTimeVector(ipnt)
       integer, intent(in) :: ipnt
       !
       integer :: i, idump,k
       real(kind=rp), parameter :: dt = 0.5_rp
       real(kind=rp) :: dr, dy, temp, z_source
       !
       integer :: ntab_data
       real(kind=rp), dimension(:), allocatable :: clgr_table, va_table, a_table
       
       ntab_data = 2+2*n_traj_pts ! correct for 1D table
       allocate(clgr_table(ntab_data)) ; allocate(va_table(ntab_data)) ; allocate(a_table(ntab_data))  
       call set_table(n_traj_pts,x,clgr,clgr_table)
       call set_table(n_traj_pts,x,  va,  va_table)
       call set_table(n_traj_pts,x,   a,   a_table)
       
! create tmic (starting from when the first noise reaches the microphone from xstart (i.e. r(1)/a(1)) to when the noise reaches the microphone from the last point)
        call init_vector(r(1)/a(1),time(n_traj_pts)+r(n_traj_pts)/a(n_traj_pts),0.5_rp,n_times,tmic)
        
        if( allocated(t_source) ) deallocate(t_source)
        allocate(t_source(n_times)) 

        if( allocated(x_source) ) deallocate(x_source)
        allocate(x_source(n_times)) 
        
        if( allocated(y_source) ) deallocate(y_source)
        allocate(y_source(n_times)) 

        if( allocated(Va_source) ) deallocate(Va_source)
        allocate(Va_source(n_times)) 

        if( allocated(clgr_source) ) deallocate(clgr_source)
        allocate(clgr_source(n_times)) 
        
        if( allocated(a_source) ) deallocate(a_source)
        allocate(a_source(n_times)) 
        
        if( allocated(r1) ) deallocate(r1)
        allocate(r1(n_times)) 
        
        t_source(1) = zero ; 
        x_source(1) = x(1) - xmic(ipnt) ! use mic-related coordinate system
        y_source(1) = y(1) - ymic(ipnt) ! use mic-related coordinate system
        va_source(1) = va(1) ; a_source(1) = a(1) ; clgr_source(1) = clgr(1)
        z_source = zmic(ipnt)
        r1(1) = sqrt(x_source(1)**2 + y_source(1)**2 + z_source**2) ! distance between source and microphone. 
        do i = 2, n_times
           if( va_source(i-1).eq.zero ) then
               dx = zero
           else               
               dx = get_dx(x_source(i-1),y_source(i-1),z_source,dt,clgr_source(i-1),a_source(i-1),va_source(i-1),r1(i-1)/a_source(i-1)+dt) 
           end if
           ! update using dx
           x_source(i) = x_source(i-1) + dx 
           dy = dx*tand(clgr_source(i-1))
           y_source(i) = y_source(i-1) + dy
           dr = sqrt(dx**2 + dy**2)
           r1(i) = sqrt(x_source(i)**2 + y_source(i)**2 + z_source**2) 
           if( dr.eq.zero ) then
               t_source(i) = t_source(i-1) + dt
           else
               t_source(i) = t_source(i-1) + dr/va_source(i-1)
           end if
           ! interpolate
           call linear_interp(clgr_table,1,x_source(i)+xmic(ipnt),zero,clgr_source(i),idump)
           call linear_interp(  va_table,1,x_source(i)+xmic(ipnt),zero,  va_source(i),idump)
           call linear_interp(   a_table,1,x_source(i)+xmic(ipnt),zero,   a_source(i),idump)
           !
        end do 
        
        !mess first was here
        !!! my messs
        
        
    !if(ipnt.EQ.2) then
        
    open (unit = 30, file = "C:\Users\fatbah\Documents\Approach_comb_performance\tsource_Cutback.txt")
    
    do i=1, n_times
       write (30,9) i,t_source(i),x_source(i),y_source(i)
    end do

    close(30)
    
9   format (i3,f16.4,f16.4,f16.4)
    !end if
    
    
    
    !!!!!!!!!!
        
        
        
        deallocate(clgr_table) ; deallocate(va_table) ; deallocate(a_table)  

!  get_trajectory_variable(x_source(i)+xMic_cutback,va) <============ PREPARED FOR GENERALIZATION TO TRAJECTORY
        
     end subroutine setTimeVector
     !
     subroutine setVelocity() ! define velocity along flight trajectory
       real(kind=rp) :: R, gamma
       integer :: i, slen
       
      
       if( allocated(ta) ) deallocate(ta) ; if( allocated(a) ) deallocate(a) ; ; if( allocated(Ma) ) deallocate(Ma)
       allocate(ta(n_traj_pts))  ; allocate(a(n_traj_pts)) ; allocate(Ma(n_traj_pts))
       
       R = get_R(zero)
       do i = 1, n_traj_pts
         ta(i) = get_ambient_temperature(y(i),dtIsa)
         a(i) = sqrt(gamma_air*R*ta(i))
         Ma(i) = Va(i)/a(i)
       end do
       
     end subroutine setVelocity
     !
     subroutine setTrajectory(ipnt,opPnt)
       integer, intent(in) :: ipnt
       character(len=*), intent(in) :: opPnt
       logical :: cutback ! removeMe 
       integer :: i_last_climb ! removeMe
       integer :: slen
       real(kind=rp) :: temp, V_avg
       integer :: i,m 
       
       call read_trajectory_input(opPnt,n_traj_pts,x,y,Va,alpha)
       
      
    
       if( opPnt(1:8).eq.'Approach' ) x = parse_xVecApproach(n_traj_pts,x) 
                            
       if( allocated(r) ) deallocate(r) ; allocate(r(n_traj_pts)) 
       do i = 1, n_traj_pts
          r(i) = sqrt( (x(i)-xmic(ipnt))**2 + (y(i)-ymic(ipnt))**2 + zmic(ipnt)**2 ) ! distance along trajectory to microphone 
       end do 
       !
       if( allocated(clgr) ) deallocate(clgr) ; allocate(clgr(n_traj_pts)) 
       do i = 1, n_traj_pts-1
          if( y(i+1).eq.zero ) then
             clgr(i) = zero 
          else
             clgr(i) = (180.0_rp/pi_num)*atan( (y(i+1)-y(i)) / (x(i+1)-x(i)) ) ! clgr along trajectory
          end if
          if( ISNAN(clgr(i)) ) stop 'not a number (NaN) computed for clgr'
       end do 
       clgr(n_traj_pts) = clgr(n_traj_pts-1)
!       call save_x_versus_y('clgr_versus_x',n_traj_pts,n_traj_pts,x,clgr)
       !
       if( allocated(xsi) ) deallocate(xsi) ; allocate(xsi(n_traj_pts))         
       if( allocated(xsid) ) deallocate(xsid) ; allocate(xsid(n_traj_pts))         
       do i = 1, n_traj_pts
          xsi(i) = get_xsi(clgr(i),x(i),y(i),xmic(ipnt),ymic(ipnt)) ! Xsi angle along trajectory 
          xsid(i) = rad2deg(xsi(i))
          if( ISNAN(xsi(i)) ) stop 'not a number computed for xsi'
       end do 
!       call save_x_versus_y('xsi_versus_x',n_traj_pts,n_traj_pts,x,xsid)

! time is computed for 2D trajectory
       if( allocated(time) ) deallocate(time) ; allocate(time(n_traj_pts)) 
       time(1) = zero
       do i = 2, n_traj_pts
         if( sqrt( (y(i)-y(i-1))**2 + (x(i)-x(i-1))**2 ).eq.zero  ) then
           time(i) = time(i-1) + one ! aircraft is not moving - impossible to compute time => set a fictitios time
         else
           ! use average velocity. GESTPAN will occasionally register a small velocity value
           ! prior to any distance covered => jf this velocity is used to compute time elapsed it may become very large.
           V_avg = (Va(i)+Va(i-1))/two 
           time(i) = time(i-1) + sqrt( (y(i)-y(i-1))**2 + (x(i)-x(i-1))**2 )/V_avg
           if( ISNAN(time(i)) ) stop 'not a number computed for time'
         end if   
       end do 

!       call save_x_versus_y('time_versus_x',n_traj_pts,n_traj_pts,x,time)

 !!! my messs
        
    !if(ipnt.eq.2) then
        open (unit = 31, file = "C:\Users\fatbah\Documents\Approach_comb_performance\trajec_Cutback.txt")
    do m=1, n_traj_pts
       write(31,7) m,time(m),x(m)-xmic(ipnt),y(m)-ymic(ipnt)
    end do    
    close(31)
7   format (i3,f16.4,f16.4,f16.4)
    !end if
    
       
    end subroutine setTrajectory  
  end subroutine set_trajectory_subr
  !
  !
  !
  subroutine set_table(npts,x,y,table)  
    integer, intent(in) :: npts
    real(kind=rp), intent(in), dimension(npts) :: x, y
    real(kind=rp), intent(out), dimension(*) :: table
    !
    integer :: istart, iend
    integer :: i, ntab_data
    ! 
    if( x(1).eq.x(2) ) then ! we have constant values in the beginning
      do i = 1, npts
          if( x(i+1).ne.x(i) ) then
              istart = i
              exit
          end if
      end do
    else
      istart = 1
    end if
    
    if( x(npts).eq.x(npts-1) ) then ! we have constant values in the end
      do i = npts,1,-1 
          if( x(i-1).ne.x(i) ) then
              iend = i
              exit
          end if
      end do
    else
      iend = npts
    end if
    
    ntab_data = 2+2*(iend-istart+1)
    table(1:ntab_data) = get_gestpan_1d_table(iend-istart+1,x(istart:iend),y(istart:iend))
  
  end subroutine set_table
  !
  !
  !
  subroutine set_noise_choice(opPnt,nops,trajectory_performance,use_trajectory_preparser,gen_noise_source_matr, & 
       & use_ground_reflection,use_spherical_spreading,use_atmospheric_attenuation,nlines,noiseFile)
    character(len=*), intent(out), dimension(3) :: opPnt
    integer, intent(out) :: nops
    logical, intent(out) :: trajectory_performance, gen_noise_source_matr, use_ground_reflection, & 
       & use_spherical_spreading, use_atmospheric_attenuation, use_trajectory_preparser
    integer, intent(in) :: nlines
    character(len=32), dimension(3) :: str
    character(len=cline_len), dimension(*), intent(in) :: noiseFile
    !
    integer :: nvals
    
    opPnt(1) = get_string('point:',nlines,noiseFile) ! x-lift
    nops = get_n_tokens_in_string(opPnt(1))
    if( nops.gt.maxNoPts ) stop 'number of noise points greater than allowed'
    
    trajectory_performance = .false.
    if( index( get_string('trajectory_performance:',nlines,noiseFile),'true' ).gt.0 ) trajectory_performance = .true.
    
    use_trajectory_preparser = .false.
    if( index( get_string('use_trajectory_preparser:',nlines,noiseFile),'true' ).gt.0 ) use_trajectory_preparser = .true.

    gen_noise_source_matr = .false.
    if( index( get_string('generate_noise_sources:',nlines,noiseFile),'true' ).gt.0 ) gen_noise_source_matr = .true.

    use_ground_reflection = .true. ! default is true
    if( index( get_string('use_ground_reflection:',nlines,noiseFile),'false' ).gt.0 ) use_ground_reflection = .false.

    use_spherical_spreading = .true. ! default is true
    if( index( get_string('use_spherical_spreading:',nlines,noiseFile),'false' ).gt.0 ) use_spherical_spreading = .false.

    use_atmospheric_attenuation = .true. ! default is true
    if( index( get_string('use_atmospheric_attenuation:',nlines,noiseFile),'false' ).gt.0 ) use_atmospheric_attenuation = .false.
    
    !
    if( nops.eq.1 ) then
       ! 
       ! opPnt already ok from above
      xmic(1) = get_value('xmic:',nlines,noiseFile)
      ymic(1) = get_value('ymic:',nlines,noiseFile)
      zmic(1) = get_value('zmic:',nlines,noiseFile)
      !
    elseif( nops.gt.1 ) then ! if more than 1 we assume three points
       !
       if( nops.gt.maxNoPts ) stop 'number of noise points greater than allowed'
       opPnt(1:nops) = tokenize_string(opPnt(1),nops)
       ! 
       call set_microphone_positions() ! set xmic(1:3), ymic(1:3), zmic(1:3)
       !
    end if
    
    no_engines = get_value('no_engines:',nlines,noiseFile)
    dtIsa = get_value('dtIsa:',nlines,noiseFile) 
    total_weight_airfrm = get_value('aircraft_mass:',nlines,noiseFile) 
    width_wheelm_airfrm = get_value('width_wheelm:',nlines,noiseFile)
    width_wheeln_airfrm = get_value('width_wheeln:',nlines,noiseFile)
    d_wheelm_airfrm= get_value('d_wheelm:',nlines,noiseFile)
    d_wheeln_airfrm= get_value('d_wheeln:',nlines,noiseFile)
 

!    spann_aileron_airfrm= get_value('semi_span_airframe:',nlines,noiseFile)
!    lscale_aileron_airfrm= get_value('lscale_aileron:',nlines,noiseFile)
!    chord_outb_flap_airfrm = get_value('chord_outb_flap:',nlines,noiseFil
! Airframe Noise Sub-Component Definition and Model , Rahul Sen et al, NASA/CR-2004-213255 - Chapter 3
    spann_aileron_airfrm  = 17.05_rp ! hard coded until parameter significance has been understood
    lscale_aileron_airfrm = 2.7_rp ! hard coded until parameter significance has been understood
    chord_outb_flap_airfrm = 0.5_rp ! hard coded until parameter significance has been understood
    
    
    NoAil_airfrm = get_value('noAil:',nlines,noiseFile)
    NoInb_airfrm = get_value('noInb_flaps:',nlines,noiseFile)
    NoOut_airfrm = get_value('noOut_flaps:',nlines,noiseFile)
    NoSlat_airfrm = get_value('noSlats:',nlines,noiseFile)
    
    Wref_airfrm = get_value('Wref_Boeing737:',nlines,noiseFile)
    Nref_airfrm = get_value('Nwheels_Boeing737:',nlines,noiseFile)
    Lref_airfrm = get_value('Lstruts_Boeing737:',nlines,noiseFile)
    
    N_wheelsm_airfrm = get_value('Nwheelsm:',nlines,noiseFile)
    N_wheelsn_airfrm = get_value('Nwheelsn:',nlines,noiseFile)
         
    
    nlsm_airfrm = get_value('Nostrutsm:',nlines,noiseFile)  
    
        
    call setVectorData('length_strutm:',nlines,noiseFile,nvals,length_strutm_airfrm)
    call setVectorData('diameter_strutm:',nlines,noiseFile,nvals,d_strutm_airfrm)
    call setVectorData('depth_strutm:',nlines,noiseFile,nvals,geometric_strutm_airfrm)
    
    nlsn_airfrm = get_value('Nostrutsn:',nlines,noiseFile) 
    
    call setVectorData('length_strutn:',nlines,noiseFile,nvals,length_strutn_airfrm)
    call setVectorData('diameter_strutn:',nlines,noiseFile,nvals,d_strutn_airfrm)
    call setVectorData('depth_strutn:',nlines,noiseFile,nvals,geometric_strutn_airfrm)    
    
    call setVectorData('psi:',nlines,noiseFile,nvals,psi_airfrm_vec)    
    call setVectorData('defl_flap:',nlines,noiseFile,nvals,defl_flap_airfrm_vec)    

! 
!    call setVectorData('cl_outb_flap:',nlines,noiseFile,nvals,Cl_outb_flap_airfrm_vec)    
!    call setVectorData('aoa_ail:',nlines,noiseFile,nvals,aoa_ail_airfrm_vec)    
!    call setVectorData('aoa_slat:',nlines,noiseFile,nvals,aoa_slat_airfrm_vec )    

    allocate(Cl_outb_flap_airfrm_vec(1:nvals)) ; allocate(aoa_ail_airfrm_vec(1:nvals)) ; allocate(aoa_slat_airfrm_vec(1:nvals))  
    Cl_outb_flap_airfrm_vec(1:nvals) = (/1.432065_rp,1.432065_rp,1.378180_rp/)
    aoa_ail_airfrm_vec(1:nvals) = (/0.000000_rp,0.244346_rp,0.436332_rp/)
    aoa_slat_airfrm_vec(1:nvals) = (/0.174533_rp,0.174533_rp,0.174533_rp/)

       
  contains 
     subroutine set_microphone_positions() 
        character(len=cline_len) :: temp
        integer :: i, ndata
! xmic 
        temp = get_string('xmic:',nlines,noiseFile) ; ndata = get_n_tokens_in_string(temp)
        if( ndata.ne.nops ) stop 'xmic data not consistent with number of points'
        str(1:ndata) = tokenize_string(temp,ndata)
        do i = 1, ndata
          xmic(i) = string_number_cast(str(i))
       end do
    
! ymic 
        temp = get_string('ymic:',nlines,noiseFile) ; ndata = get_n_tokens_in_string(temp)
        if( ndata.ne.nops ) stop 'ymic data not consistent with number of points'
        str(1:ndata) = tokenize_string(temp,ndata)
        do i = 1, ndata
          ymic(i) = string_number_cast(str(i))
       end do
! xmic 
        temp = get_string('zmic:',nlines,noiseFile) ; ndata = get_n_tokens_in_string(temp)
        if( ndata.ne.nops ) stop 'zmic data not consistent with number of points'
        str(1:ndata) = tokenize_string(temp,ndata)
        do i = 1, ndata
          zmic(i) = string_number_cast(str(i))
        end do

     end subroutine set_microphone_positions
  end subroutine set_noise_choice
  !
  !
  !
  subroutine setVectorData(str,nl,noiseFile,nvals,vec)
    character(len=*), intent(in) :: str
    integer, intent(in) :: nl
    character(len=cline_len), dimension(*), intent(in) :: noiseFile
    integer, intent(out) :: nvals 
    real(kind=rp), allocatable, dimension(:), intent(out) :: vec
    !
    logical :: found
    integer :: i, j, cp
    character(len=32), dimension(:), allocatable :: strvec
    
    found = .false.
    
    
    do i = 1, nl
        cp = index(noiseFile(i),':')
        if( cp.eq.0 ) cycle
        if( noiseFile(i)(1:cp).eq.str ) then
            found = .true.
            !
            nvals = get_n_tokens_in_string(noiseFile(i)(cp+1:cline_len))
            allocate(vec(nvals)) ; allocate(strvec(nvals))
            vec = zero ; strvec = ' '
            strvec = tokenize_string(noiseFile(i)(cp+1:cline_len),nvals)
            do j = 1, nvals
                vec(j) = string_number_cast(strvec(j))
            end do
            !
            exit
        end if
    end do 
    
    deallocate(strvec)
    if( .not.found ) call report_error('name '//trim(str)//' not found in performance file','setVectorData','choice_interf')
  
  end subroutine setVectorData
  !
  ! get_inputNoise_ci: replaces functionality of choice_read_and_write.f90 which is used when choice is run in stand alone mode.
  !
  function get_inputNoise_ci(nchars,nlns,fname)
    integer, intent(in) :: nchars
    integer, intent(in) :: nlns
    character(len=*), intent(in) :: fname
    character(len=nchars), dimension(nlns) :: get_inputNoise_ci
    !
    integer :: wa_unit ! weight aircraft unit
     
    wa_unit = get_unit_number()
    open(wa_unit,file=fname,status='old')
  
    call load_file(wa_unit, nlns, get_inputNoise_ci)    
    
    close(wa_unit)    
    
  end function get_inputNoise_ci
  !
  !
  !
  subroutine read_trajectory_input(opPnt,n_traj_pts,x,y,Va,alpha)
    character(len=*), intent(in) :: opPnt
    integer, intent(out) :: n_traj_pts 
    real(kind=rp), dimension(:), allocatable, intent(out) :: x, y, Va, alpha
    !
    integer :: unit, i
    character(len=cline_len) :: temp    
    character(len=32), dimension(4) :: str
    !
    unit = get_unit_number()
    open(unit,file=trim(opPnt)//'.txt',status='old') 
    
    n_traj_pts = get_n_lines(unit)
  
    allocate(x(n_traj_pts)) ; x = zero
    allocate(y(n_traj_pts)) ; y = zero
    allocate(Va(n_traj_pts)) ; Va = zero
    allocate(alpha(n_traj_pts)) ; alpha = zero

    i = 1
    do 
       read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
       if( is_a_comment_line(temp) ) cycle
       ! data line
       !
       str(1:4) = tokenize_string(temp,4)
       x(i) = string_number_cast(str(1)) ; y(i) = string_number_cast(str(2)) ; Va(i) = string_number_cast(str(3)) ; alpha(i) = string_number_cast(str(4))
       !
       if( i.eq.n_traj_pts ) exit
       !
       i = i + 1 
       !
    end do

100 close(unit)
  
  end subroutine read_trajectory_input
  !
  ! set_weights_choice: transfer data to choice_physics, variable by variable.
  ! The whole input file (weightAircraft.txt) is transferred as argument. 
  !
  subroutine set_weight_choice(nlines,weightFile)
    integer, intent(in) :: nlines
    character(len=cline_len), dimension(*), intent(in) :: weightFile
    !
    integer :: i 
    
    do i = 1, n_modules
      if( trim(modules(i)).eq.'Fan' ) call setFan()
      if( trim(modules(i)).eq.'Lpc' .or. trim(modules(i)).eq.'Ipc' ) call setLpc()
      if( trim(modules(i)).eq.'Lpt' ) call setLpt()
      if( trim(modules(i)).eq.'Comb' ) call setComb()
      if( trim(modules(i)).eq.'cold_nozzle' ) call setCaj()
    end do 
  
  contains 
     subroutine setFan()
      
        MtipD_fan = get_value('MtipFan:',nlines,weightFile) 
        if( MtipD_fan.gt.3.0_rp .or. MtipD_fan.lt.0.7_rp ) call report_error('unexpected value on MtipD_fan', & 
           & 'setFan','set_weight_choice')
        
        xnlD_fan = get_value('xnl:',nlines,weightFile)
        if( xnlD_fan.gt.600 ) call report_error('xnl from weightFile very large. Unit should always be rps!', & 
           & 'setFan','set_weight_choice')
        
        rss_fan = get_value('FanRss:',nlines,weightFile) 
        if( rss_fan.gt.300.0_rp .or. rss_fan.lt.5.0_rp ) call report_error('unexpected value on rss_fan', & 
           & 'setFan','set_weight_choice')

        N_rotors_fan  = get_value('FanR1BNb:',nlines,weightFile) 
        N_stators_fan = get_value('FanVsOgvNb:',nlines,weightFile) 
        
        A2_fan = get_value('FanA2:',nlines,weightFile) 
        if( A2_fan.gt.50.0_rp .or. A2_fan.lt.0.1_rp ) call report_error('unexpected value on A2_fan', & 
           & 'setFan','set_weight_choice')

        D1_fan = get_value('FanR1BDiaOuter:',nlines,weightFile) 
        if( D1_fan.lt.0.1_rp .or. D1_fan.gt.10.0_rp ) call report_error('unexpected value on D1_fan', & 
           & 'setFan','set_weight_choice')


     end subroutine setFan
     !
     subroutine setLpc()
     
        n_stages_LPC = get_value('stages_LPC:',nlines,weightFile) 

        GBX_ratio = get_value('GBX_ratio:',nlines,weightFile) 
        MtipD_lpc = get_value('MtipLpc:',nlines,weightFile) 
        xnlD_lpc = xnlD_fan*GBX_ratio
        rss_lpc = get_value('RSS_compr:',nlines,weightFile) 

        N_rotors_lpc = get_value('N_compr:',nlines,weightFile) 
        N_stators_lpc = get_value('S_compr:',nlines,weightFile) 

        D1_lpc  = two*get_value('r_compr:',nlines,weightFile) 
        Dh1_lpc = two*get_value('rh_compr:',nlines,weightFile) 
        
     end subroutine setLpc
     !
     subroutine setLpt()
     
        N_rotors_lpt = get_value('LptStgLastBNb:',nlines,weightFile) 
        SRS_lpt = get_value('SRS:',nlines,weightFile) 
        n_stages_lpt = get_value('stages_LPT:',nlines,weightFile) 
        De_lpt = get_value('LptStgLastDiaOuter:',nlines,weightFile) 
        Ae_lpt = get_value('LptStgLastAExit:',nlines,weightFile) 

     end subroutine setLpt
     !
     subroutine setComb()
     
        type_comb = get_string('CombType:',nlines,weightFile)

     end subroutine setComb
     !
     subroutine setCaj()
     
        A_core_caj = get_value('A_core:',nlines,weightFile)
        A_bypass_caj = get_value('A_bypass:',nlines,weightFile)
     
     end subroutine setCaj 
  end subroutine set_weight_choice
  !
  !
  !
  function is_in_model(str)
    character(len=*), intent(in) :: str
    logical :: is_in_model
    integer :: i
    
    is_in_model = .false.
    
    do i = 1, n_modules
      if( trim(modules(i)).eq.str ) then
        is_in_model = .true.
      end if
    end do
    
  end function is_in_model
  !
  !
  !
  subroutine set_configuration_choice(nmods,mods)
    integer, intent(in) :: nmods ! n_modules
    character(len=cstr_len), intent(in), dimension(nmods) :: mods
    
    n_modules = nmods
    
    if( .not.allocated(modules) ) allocate(modules(n_modules))
    modules(1:n_modules) = mods(1:n_modules)
  
  end subroutine set_configuration_choice 
  !
  !
  !
  function get_dim_variable(spos,epos,dimFile,str)
    integer, intent(in) :: spos, epos
    character(len=cline_len), dimension(*), intent(in) :: dimFile
    character(len=*), intent(in) :: str
    real(kind=rp) :: get_dim_variable
    integer :: cp
    integer :: i 
    !
    do i = spos+1, epos-1
       if( index(trim(dimFile(i)),str) .gt. 0 ) then
         cp = index(dimFile(i),':')
         get_dim_variable = string_number_cast(dimFile(i)(cp+1:cline_len))  
         exit
       end if
    end do
  
  end function get_dim_variable
  !
  !
  !
  subroutine set_spos_and_epos(nlines,dimFile,strs,stre,spos,epos)
    integer, intent(in) :: nlines
    character(len=cline_len), dimension(nlines), intent(in) :: dimFile
    character(len=*), intent(in) :: strs ! string start
    character(len=*), intent(in) :: stre ! string end
    integer, intent(out) :: spos, epos
    !
    character(len=cline_len) :: temp
    logical :: endModuleFound , ModuleFound 
    integer :: i, slen
    
    spos = 0 ; epos = 0 ; endModuleFound = .false. ; ModuleFound = .false.
    slen = len(stre)
    
    do i = 1, nlines
       temp = dimFile(i)
       if( temp(1:slen).eq.stre ) then
         endModuleFound = .true.
         epos = i
         exit
       end if
   end do 
    
    if( .not.endModuleFound ) call report_error('no end module tag found','set_spos_and_epos','choice_interf')

    slen = len(strs)
    do i = epos-1,1,-1
       temp = dimFile(i)
       if( temp(1:slen).eq.strs ) then
         spos = i
         ModuleFound = .true.
         exit
       end if
    end do 
    if( .not.ModuleFound ) call report_error('no end module tag found','set_spos_and_epos','choice_interf')
    
  end subroutine set_spos_and_epos
  !
  !
  !
  subroutine certificationLimits(noEngines,MTOW)
    integer, intent(in) :: noEngines
    real(kind=rp) :: MTOW ! maximum take-off weight in tonnes
    !
    real(kind=rp) :: EPNL_lateral,EPNL_cutback,EPNL_approach
    real(kind=rp) :: EPNL_cum

    call chapter3(noEngines,MTOW,EPNL_lateral,EPNL_cutback,EPNL_approach)
    
    EPNL_cum = EPNL_lateral + EPNL_cutback + EPNL_approach
    
  end subroutine certificationLimits
end module choice_interf
