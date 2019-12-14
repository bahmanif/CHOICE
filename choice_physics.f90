 module choice_physics
  !
  ! Module programmer: Tomas Grönstedt
  !
  ! v1.0 17.04.01  Date written.
  !
  ! Physics for the noise prediction broken down on modules
  !
  use choice_data
  use choice_aux
  use sim_precision
  use units_and_constants
  use interpolation
  use choice_aux
  use error_handling
  !
  implicit none
  !
  public :: test_functions, get_ambient_temperature, get_trajectory_variable, flightEffects, get_vector_interpolation, & 
     & get_theta, get_closest_ptr, getPNL, setMachNumbers, getVel
  public :: calcFanandCompressor, calcTurbine, calcComb, calcCaj, get_xsi, get_atm_abs
  !
  real(kind=rp), public :: dtIsa ! delta from ISA 
  integer :: nfrec ! number of frequencies
  !
  integer, public :: no_engines
  integer, parameter, public :: fmin = 50 
  integer, parameter, public :: fmax = 10000 
  integer, parameter, public :: nb = 1  ! number of bands
  integer, parameter, public :: nfreq = nint(three*nb*log10(float(fmax)/float(fmin))/log10(two)) + 1 ! for fmin = 50 and fmax = 10000 this computes to 24
  !
  real(kind=rp), parameter :: theta_step = 5.0_rp
  integer, parameter, public :: nthet = nint(180.0_rp/theta_step) + 1  ! number of thetas
  real(kind=rp), dimension(nthet), public :: theta
  real(kind=rp), dimension(:), allocatable, public :: source_theta
  integer, dimension(:), allocatable, public :: closest_value_ptr
  !
  real(kind=rp), dimension(nfreq), public :: fband ! for standard settings of noPred,  (/50.0000, 62.9961,....,10159.0)
  real(kind=rp), dimension(nfreq+1), public :: f
  real(kind=rp), dimension(fmax-fmin+1), public :: freq
  !
  real(kind=rp), parameter, public :: RH = 70.0_rp
  real(kind=rp), public, parameter :: pa2psia = 6.9E-3_rp 
  !
  real(kind=rp), parameter, public :: V_ref_lpt = 0.305_rp ! Reference velocity
  real(kind=rp), parameter, public :: c_ref_lpt = 340.3_rp ! Reference speed of sound
  real(kind=rp), parameter, public :: K_lpt = -10.0_rp ! Correction factor depending on nozzel exit plane, if not panel symmetric set K=-10 
  real(kind=rp), public :: xMic_cutback = 6500.0_rp
  !
  !
  !
  !----------------------------------------
  ! FAN DATA 
  logical, parameter, private :: fan_IGV = .false.
  logical, parameter, private :: fan_distortion = .true.
  integer, parameter, private :: no_fan_stages = 1
  !
  real(kind=rp), public :: N_rotors_fan ! N_fan => N_rotors_fan
  real(kind=rp), public :: N_stators_fan ! N_fan => N_stators_fan
  real(kind=rp), public :: MtipD_fan ! M_trD => MtipD_fan
  real(kind=rp), public :: xnlD_fan  ! top of climb rotational speed (must be in rps)
  real(kind=rp), public :: D1_fan, A2_fan ! fan face area
  real(kind=rp), public :: rss_fan ! fan rotor stator spacing
  !
  real(kind=rp), public, allocatable, dimension(:) :: Mtip_fan ! Relative tip Mach number 
  real(kind=rp), public, allocatable, dimension(:) :: Mu_fan ! Blade Mach number 
  real(kind=rp), public, allocatable, dimension(:) :: dt_fan ! Stage temperature rise over fan rotor.
  real(kind=rp), public, allocatable, dimension(:) :: xnl_fan ! Rotational speed (in rps) 
  real(kind=rp), public, allocatable, dimension(:) :: g1_fan ! Mass flow at fan face. 
  !
  real(kind=rp), public, dimension(nfreq,nthet) :: prms_inlet_broadband, prms_discharge_broadband
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_inlet_fan, prms_discharge_fan
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_inlet_tone_fan, prms_discharge_tone_fan, &  
      & prms_inlet_broadband_fan, prms_discharge_broadband_fan,prms_inlet_combination_fan
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_inlet_tone_Lpc, prms_discharge_tone_Lpc, &  
      & prms_inlet_broadband_Lpc, prms_discharge_broadband_Lpc, prms_inlet_combination_Lpc
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prmsi_inlet_fan, prmsi_discharge_fan ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:,:), allocatable :: SPLi_inlet_fan, SPLi_discharge_fan ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:), allocatable :: SPLp_inlet_fan, SPLp_discharge_fan ! propageted onto observer. Only one directivity propagates.
  real(kind=rp), public, dimension(:,:), allocatable :: prmsp_inlet_fan, prmsp_discharge_fan ! 
  real(kind=rp), public, dimension(:), allocatable :: PNL_inlet_fan, PNL_discharge_fan
  real(kind=rp), public, dimension(:), allocatable :: PNLT_inlet_fan, PNLT_discharge_fan
  real(kind=rp), public :: EPNL_inlet_fan, EPNL_discharge_fan
  ! LPC
  integer, public :: n_stages_LPC
  real(kind=rp), public :: gbx_ratio ! gear box ratio (used to compute xnld_lpc from xnld_fan) 
  real(kind=rp), public :: MtipD_lpc  
  real(kind=rp), public :: rss_lpc
  real(kind=rp), public :: xnld_lpc
  real(kind=rp), public :: N_stators_lpc,N_rotors_lpc
  !
  real(kind=rp), public :: D1_lpc, Dh1_lpc
  real(kind=rp), public, allocatable, dimension(:) :: Mtip_lpc ! Relative tip Mach number 
  real(kind=rp), public, allocatable, dimension(:) :: Mu_lpc ! Blade Mach number 
  real(kind=rp), public, allocatable, dimension(:) :: dt_lpc ! Stage temperature rise over fan rotor.
  real(kind=rp), public, allocatable, dimension(:) :: xnl_lpc ! Rotational speed (in rps) 
  real(kind=rp), public, allocatable, dimension(:) :: g1_lpc ! Mass flow at fan face. 
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_inlet_lpc, prms_discharge_lpc
  real(kind=rp), public, dimension(:,:,:), allocatable :: prmsi_inlet_lpc ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:,:), allocatable :: SPLi_inlet_lpc ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:), allocatable :: SPLp_inlet_lpc ! propageted onto observer. Only one directivity propagates.
  real(kind=rp), public, dimension(:,:), allocatable :: prmsp_inlet_lpc ! 
  real(kind=rp), public, dimension(:), allocatable :: PNL_inlet_lpc 
  real(kind=rp), public, dimension(:), allocatable :: PNLT_inlet_lpc
  real(kind=rp), public :: EPNL_inlet_lpc
  ! LPT
  real(kind=rp), public :: N_rotors_lpt, SRS_lpt, De_lpt, n_stages_lpt, Ae_lpt
  real(kind=rp), public, allocatable, dimension(:) :: Vtr_lpt ! Relative tip speed
  real(kind=rp), public, allocatable, dimension(:) :: Texit_lpt ! Exit temperature 
  real(kind=rp), public, allocatable, dimension(:) :: xnl_lpt ! Exit temperature 
  real(kind=rp), public, allocatable, dimension(:) :: mcore_lpt ! Mass flow
  real(kind=rp), public, allocatable, dimension(:) :: Cax_lpt ! Axial velocity
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_lpt
  real(kind=rp), public, dimension(:,:,:), allocatable :: prmsi_lpt ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:,:), allocatable :: SPLi_lpt ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:), allocatable :: SPLp_lpt ! propageted onto observer. Only one directivity propagates.
  real(kind=rp), public, dimension(:,:), allocatable :: prmsp_lpt ! 
  real(kind=rp), public, dimension(:), allocatable :: PNL_lpt 
  real(kind=rp), public, dimension(:), allocatable :: PNLT_lpt
  real(kind=rp), public :: EPNL_lpt
  ! Combustor
  character(len=32), public :: type_comb ! SAC or DAC
  real(kind=rp), public, parameter :: Nf_comb = 12.0_rp      ! number of ignited fuel nozzles for annular combustor
  real(kind=rp), public, parameter :: R0_comb = 1.0_rp       ! observation radius (ft) 
  real(kind=rp), public, parameter :: Aec_comb = 5.0_rp      ! combustor exit area (ft^2)
  real(kind=rp), public, parameter :: De_comb = 0.4_rp       ! exhaust nozzle exit plane effective diameter (ft)
  real(kind=rp), public, parameter :: Dh_comb =  5.0_rp      ! exhaust nozzle exit plane hydraulic diameter (ft)  
  real(kind=rp), public, parameter :: Lc_comb = 1.0_rp       ! combustor nominal length (ft) 
  real(kind=rp), public, parameter :: h_comb = 0.4_rp        ! annulus height at combustor exit (ft)
  real(kind=rp), public, parameter :: Nfmax_comb = 12.0_rp   ! total number of dac fuel nozzles
  real(kind=rp), public, parameter :: pattern_comb = 40.0_rp ! firing pattern alternative: 20+20, 20+10, 20+2.25, 20+0.0
  !
  real(kind=rp), public, allocatable, dimension(:) :: P3_comb
  real(kind=rp), public, allocatable, dimension(:) :: P4_comb
  real(kind=rp), public, allocatable, dimension(:) :: P7_comb
  real(kind=rp), public, allocatable, dimension(:) :: T3_comb
  real(kind=rp), public, allocatable, dimension(:) :: T4_comb
  real(kind=rp), public, allocatable, dimension(:) :: T5_comb
  real(kind=rp), public, allocatable, dimension(:) :: W3_comb
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_comb
  real(kind=rp), public, dimension(:,:,:), allocatable :: prmsi_comb ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:,:), allocatable :: SPLi_comb ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:), allocatable :: SPLp_comb ! propageted onto observer. Only one directivity propagates.
  real(kind=rp), public, dimension(:,:), allocatable :: prmsp_comb ! 
  real(kind=rp), public, dimension(:), allocatable :: PNL_comb
  real(kind=rp), public, dimension(:), allocatable :: PNLT_comb
  real(kind=rp), public :: EPNL_comb
  ! Airframe
  real(kind=rp), public :: total_weight_airfrm
  real(kind=rp), public :: width_wheelm_airfrm 
  real(kind=rp), public :: d_wheelm_airfrm ! Wheel diameter
  real(kind=rp), public :: d_wheeln_airfrm  ! Nose: Wheel diameter 
  real(kind=rp), public :: width_wheeln_airfrm  ! Nose: Wheel width
  real(kind=rp), public :: spann_aileron_airfrm     ! Semi-span
  real(kind=rp), public :: lscale_aileron_airfrm   ! Length scale aileron
  real(kind=rp), public :: chord_outb_flap_airfrm                  ! Outboard flap edge chord
  real(kind=rp), public, parameter :: circulation_flap_airfrm = one                  ! Circulation of flap
  real(kind=rp), public, parameter :: v_spanwise_flap_airfrm  = 35.0_rp              ! Spanwise velocity
  real(kind=rp), public, parameter :: r_const_airfrm          = one                  ! Far-field spectra distance
  !
  integer, public :: NoAil_airfrm 
  integer, public :: NoInb_airfrm  ! number of inboard flaps
  integer, public :: NoOut_airfrm   ! number of outboard flaps
  integer, public :: NoSlat_airfrm  ! number of slats
  !
  ! Geometric reference values
  real(kind=rp) :: Wref_airfrm  ! Weight of Boeing 737 [lb]
  real(kind=rp) :: Nref_airfrm  ! Number of wheels of Boeing 737
  real(kind=rp) :: Lref_airfrm  ! Length of struts of Boeing 737 [inch]
  !
  ! landing gear parameters
  integer :: nlsm_airfrm
  real(kind=rp), public, allocatable, dimension(:) :: length_strutm_airfrm ! OK-Boeing777: Respective length of the struts
  real(kind=rp), public, allocatable, dimension(:) :: d_strutm_airfrm  ! OK-Boeing777  MAIN: Diameter of strut/Width of strut
  real(kind=rp), public, allocatable, dimension(:) :: geometric_strutm_airfrm  ! OK-Boeing777  MAIN: depth of strut if rectangular. OBS zero if strut is circular
  integer :: nlsn_airfrm 
  real(kind=rp), public, allocatable, dimension(:) :: length_strutn_airfrm  ! OK-Boeing777 NOSE: Respective length of the struts
  real(kind=rp), public, allocatable, dimension(:) :: d_strutn_airfrm  ! OK-Boeing777  NOSE: Diameter of strut/Width of strut
  real(kind=rp), public, allocatable, dimension(:) :: geometric_strutn_airfrm   ! OK-Boeing777  NOSE: depth of strut if rectangular. OBS zero if strut is circular
  ! performance data to be used for airframes
  integer, parameter, public :: maxNoPts = 3
  real(kind=rp), public, allocatable, dimension(:) :: psi_airfrm_vec, Cl_outb_flap_airfrm_vec, defl_flap_airfrm_vec, aoa_ail_airfrm_vec, aoa_slat_airfrm_vec 
  !
  integer, public :: N_wheelsm_airfrm   ! number of wheels in the landing gear assembly 
  integer, public :: N_wheelsn_airfrm ! number of wheels in the nose landing gear assembly   
  !
  real(kind=rp), public, dimension(:), allocatable :: psi_airfrm          ! "alpha + gamma"-angle, i.e. angle between wheel track and direction of the flow
  real(kind=rp), public, dimension(:), allocatable :: Cl_outb_flap_airfrm ! lift coefficient for outboard flap (1st flap element)
  real(kind=rp), public, dimension(:), allocatable :: defl_flap_airfrm    ! flap deflection angle
  real(kind=rp), public, dimension(:), allocatable :: aoa_ail_airfrm      ! angle of attack of aileron           
  real(kind=rp), public, dimension(:), allocatable :: aoa_slat_airfrm     ! angle_attack_slat
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_airfrm
  real(kind=rp), public, dimension(:,:,:), allocatable :: prmsi_airfrm ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:,:), allocatable :: SPLi_airfrm ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:), allocatable :: SPLp_airfrm ! propageted onto observer. Only one directivity propagates.
  real(kind=rp), public, dimension(:,:), allocatable :: prmsp_airfrm ! 
  real(kind=rp), public, dimension(:), allocatable :: PNL_airfrm
  real(kind=rp), public, dimension(:), allocatable :: PNLT_airfrm
  real(kind=rp), public :: EPNL_airfrm
  ! co-axial jet
  real(kind=rp), public :: A_core_caj, A_bypass_caj
  real(kind=rp), public :: r_const_caj = 4.0_rp
  !
  real(kind=rp), public, dimension(:), allocatable :: dmdt_1_caj
  real(kind=rp), public, dimension(:), allocatable :: dmdt_2_caj
  real(kind=rp), public, dimension(:), allocatable :: v_1_caj
  real(kind=rp), public, dimension(:), allocatable :: v_2_caj
  real(kind=rp), public, dimension(:), allocatable :: T_1_caj
  real(kind=rp), public, dimension(:), allocatable :: T_2_caj
  !
  real(kind=rp), public, dimension(:,:,:), allocatable :: prms_caj
  real(kind=rp), public, dimension(:,:,:), allocatable :: prmsi_caj ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:,:), allocatable :: SPLi_caj ! interpolated onto observer discretization
  real(kind=rp), public, dimension(:,:), allocatable :: SPLp_caj ! propageted onto observer. Only one directivity propagates.
  real(kind=rp), public, dimension(:,:), allocatable :: prmsp_caj ! 
  real(kind=rp), public, dimension(:), allocatable :: PNL_caj
  real(kind=rp), public, dimension(:), allocatable :: PNLT_caj
  real(kind=rp) :: EPNL_caj
  !
  real(kind=rp), public, dimension(3) :: xmic, ymic, zmic
  real(kind=rp), public :: dx      ! step on the flight envelope
  real(kind=rp), public :: mic_lateral ! position of mic from extended centerline
  real(kind=rp), public :: Va_takeoff, Va_cutback   
  ! raw trajectory data. Could also be read from file. Now established from a few variables and clgr as well as fixed Va.
  integer, public :: n_traj_pts
  real(kind=rp), allocatable, dimension(:), public :: xsi, xsid, r1, xsii, xsii_alpha ! trajectory variables, angle 
  real(kind=rp), allocatable, dimension(:), public :: Ma  ! Mach number on original trajectory
  real(kind=rp), allocatable, dimension(:), public :: Mai ! Mach number on "retarded observer time" interpolated 
  real(kind=rp), allocatable, dimension(:), public :: Tai ! Temperature on "retarded observer time" interpolated 
  real(kind=rp), allocatable, dimension(:,:), public :: fobs ! doppler shifted frequencies. Each point on the trajectory will have a different doppler shift. 
  real(kind=rp), allocatable, dimension(:,:), public :: atm_absorption
  real(kind=rp), allocatable, dimension(:), public :: x, y, clgr ! trajectory variables, x-position, y-position. 
  real(kind=rp), allocatable, dimension(:), public :: xi, yi ! interpolated onto trajectory variables, x-position, y-position. 
  real(kind=rp), allocatable, dimension(:), public :: r ! distance
  real(kind=rp), allocatable, dimension(:), public :: ta ! trajectory variable, ambient temperature
  real(kind=rp), allocatable, dimension(:), public :: a ! trajectory variable, local speed of sound
  real(kind=rp), allocatable, dimension(:), public :: Va ! trajectory variable, flight velocity
  real(kind=rp), allocatable, dimension(:), public :: alpha ! trajectory variable, flight velocity
  real(kind=rp), allocatable, dimension(:), public :: alphai ! trajectory variable, flight velocity
  real(kind=rp), allocatable, dimension(:), public :: time ! trajectory variable, time
  !
  integer :: n_times
  real(kind=rp), allocatable, dimension(:), public :: tmic ! tmic cutback (time from first sound wave reach microphone counted from when aircraft point for start of simulation to aircraft passes end of simulation)
  real(kind=rp), allocatable, dimension(:), public :: t_source  ! t at source, i.e. the aircraft (noise pairs with tmic)
  real(kind=rp), allocatable, dimension(:), public :: clgr_source !
  real(kind=rp), allocatable, dimension(:), public :: x_source  ! aircraft x-coordinates at noise pairs
  real(kind=rp), allocatable, dimension(:), public :: y_source  ! aircraft y-coordinates at noise pairs  
  real(kind=rp), allocatable, dimension(:), public ::  a_source ! aircraft local speed of sound 
  real(kind=rp), allocatable, dimension(:), public :: Va_source ! aircraft velocity at noise pairs  
  !
contains
  subroutine reset_choice_physics
     
     
  end subroutine reset_choice_physics
  !
  !
  !
  function get_ambient_temperature(alt,dtIsa)
    real(kind=rp), intent(in) :: alt, dtIsa
    real(kind=rp) :: get_ambient_temperature
    
    if ( alt .lt. 11000.0_rp  ) then
       get_ambient_temperature = 288.15_rp-6.5E-3_rp*alt
    elseif ( alt.ge.11000.0_rp .and. alt.lt. 25000.0_rp) then
       get_ambient_temperature = 216.65_rp
    endif
    
    get_ambient_temperature  = get_ambient_temperature  + dtIsa
  
  end function get_ambient_temperature
  !
  ! get_trajectory_variable: interpolate in x along trajectory to get other trajectory variable. 
  !                          having this routine prepares for more general trajectory definitions
  !
  function get_trajectory_variable(xval,a)
    real(kind=rp), intent(in) :: xval 
    real(kind=rp), intent(in), dimension(*) :: a 
    real(kind=rp) :: get_trajectory_variable
    !
    real(kind=rp) :: frac
    integer :: i , ptr
    !
    do i = 1, n_traj_pts-1
      if( xval.gt.x(n_traj_pts) ) then
         frac = one
         ptr = n_traj_pts-1
         exit         
      end if
      if( x(i+1).gt.xval ) then
         frac = (xval-x(i))/(x(i+1)-x(i))
         ptr = i
         exit
      end if
    end do
    
    get_trajectory_variable = a(ptr) + frac*(a(ptr+1)-a(ptr))

  end function get_trajectory_variable
  !
  !
  !
  function getVel(gam,T,M)
    real(kind=rp), intent(in) :: gam, T, M
    real(kind=rp) :: getVel
    real(kind=rp) :: ts
    
    ts = T/(one + ((gam-one)/two)*M**2)
    getVel = sqrt(gam*R_air*ts)*M
    
  end function getVel
  !
  !
  !
  subroutine calcAirfrm(cnt,operatingPoint,testCalcTurbine,NoAil,NoInb,NoOut,NoSlat,N_wheelsm,N_wheelsn,width_wheelm, & 
          & width_wheeln,d_wheelm,d_wheeln,b,lscale_aileron,chord_outb_flap,circ_flap,v_spanwise_flap,r_const,nlsm, & 
          & length_strutm,nlsn,length_strutn,d_strutm,d_strutn,geometric_strutm, geometric_strutn,total_weight, & 
          & Ta,Pa,Ma,Va,psi,Cl_outb_flap,defl_flap,aoa_ail,aoa_slat,nfrec,nt,theta,fband,f,freq,prms)
    integer, intent(in) :: cnt ! loop counter - useful for controlling debugging printouts 
    character(len=*), intent(in) :: operatingPoint
    logical, intent(in) :: testCalcTurbine ! if true => load hard coded data to allow comparison with NoPred. False => run on usual input
    integer, intent(in) :: NoAil,NoInb,NoOut,NoSlat
    integer, intent(in) :: N_wheelsm, N_wheelsn
    integer, intent(in) :: nlsm, nlsn 
    real(kind=rp), intent(in) :: width_wheelm, width_wheeln, d_wheelm, d_wheeln, b, lscale_aileron, & 
      & chord_outb_flap, circ_flap, v_spanwise_flap, r_const, Ta, Pa, Ma, Va
    real(kind=rp), intent(in), dimension(*) :: length_strutm, length_strutn, d_strutm,d_strutn,geometric_strutm, geometric_strutn
    real(kind=rp), intent(in) :: psi, Cl_outb_flap, defl_flap, aoa_ail, aoa_slat, total_weight
    integer, intent(in) :: nfrec ! number of frequencies
    integer, intent(in) :: nt ! number of theta 
    real(kind=rp), intent(in), dimension(*) :: theta ! directivity
    real(kind=rp), intent(in), dimension(nfreq) :: fband
    real(kind=rp), intent(in), dimension(nfreq+1) :: f
    real(kind=rp), intent(in), dimension(*) :: freq
    real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms
    real(kind=rp), dimension(nfreq,nthet) :: SPL_tot, Prms_landinggear_main, Prms_landinggear_nose
    !
    real(kind=rp) ::  Z, to_withdraw, rho0
    real(kind=rp),dimension(nfreq) :: AA
    real(kind=rp), dimension(nfreq) :: aileron, inboard, outboard, slat, ps_tot
    integer :: i, j, slen
    
    Z = b/47.4_rp
    
    aileron = get_aileron(nfreq,Z,b,Ma,Cl_outb_flap,defl_flap,aoa_ail,lscale_aileron,chord_outb_flap,Va,fband)
    ps_tot=real(NoAil,rp)*aileron(:)**2
    
    do i = 1, noInb
       inboard = get_inboard_flap(nfreq,Ma,Cl_outb_flap,defl_flap,v_spanwise_flap,lscale_aileron,chord_outb_flap,Va,circ_flap,fband);
       ps_tot=inboard(:)**2+ps_tot
    end do
    
    do i = 1, NoOut
       outboard = get_outboard_flap(nfreq,Z,Ma,Cl_outb_flap,defl_flap,v_spanwise_flap,lscale_aileron,chord_outb_flap,Va,circ_flap,fband);
       ps_tot=outboard(:)**2+ps_tot
    end do
    
    do i = 1, NoSlat
       slat = get_slat(nfreq,Z,Ma,Cl_outb_flap,aoa_slat,lscale_aileron,chord_outb_flap,Va,fband)
       ps_tot=slat(:)**2+ps_tot    
    end do
    
    SPL_tot(1:nfreq,1:nthet) = get_directivity(nthet,nfreq,Va,fband,theta,20.0_rp*log10(sqrt(ps_tot(:))/P0)) ! this is implemented as the transpose of the original noPred code. 
    !    call save_matrix(nfreq,nthet,SPL_tot,'.m') ; stop

    Prms=P0*10.0_rp**(SPL_tot/20.0_rp)

    to_withdraw = 394.0_rp*0.305_rp - r_const

! Atmospheric Attenuation & Spherical Spreading
    AA =  get_atm_abs(nfreq,Ta,RH,fband) ! dB/100 m
    do i = 1, nfreq
        do j = 1, nthet
          SPL_tot(i,j)=20.0_rp*log10((Prms(i,j)*to_withdraw)/P0)+(to_withdraw/100.0_rp)*AA(i)   
        end do
    end do
    Prms = P0*10.0_rp**(SPL_tot/20.0_rp) ! everything apart from landing gear noise

    slen = len(operatingPoint)
    if( operatingPoint(1:slen).eq.'Approach' ) then
        ! 
        rho0 = Pa/(R_air*Ta)
        Prms_landinggear_main = getLandingGear(nfreq,nthet,N_wheelsm,nlsm,length_strutm,Ma,sqrt(gamma_air*R_air*Ta),& 
            & d_strutm,geometric_strutm,width_wheelm,d_wheelm,psi,total_weight,rho0,r_const,fband,theta,AA)
        !
        Prms_landinggear_nose = getLandingGear(nfreq,nthet,N_wheelsn,nlsn,length_strutn,Ma,sqrt(gamma_air*R_air*Ta),& 
            & d_strutn,geometric_strutn,width_wheeln,d_wheeln,psi,total_weight,rho0,r_const,fband,theta,AA)
        !
    else
        Prms_landinggear_main = zero ; Prms_landinggear_nose = zero
    end if
    
    ! form total PRMS
    Prms = sqrt(Prms**2 + Prms_landinggear_main**2 + Prms_landinggear_nose**2)
   
!    call save_matrix(nfreq,nthet,Prms,'.m') ; stop
  
  end subroutine calcAirfrm
  !
  !        
  !
  function getLandingGear(nfr,nth,Nw,nls,Lj,M,c0,d_strut,geom_strutn,width_wheel,d_wheel,psi,W,rho0,r_const,fband,theta,AA)
    integer, intent(in) :: nfr, nth
    integer, intent(in) :: nW ! number of wheels
    integer, intent(in) :: nls
    real(kind=rp), intent(in), dimension(nls) :: Lj, d_strut, geom_strutn
    real(kind=rp), intent(in) :: M, c0, width_wheel, d_wheel, psi, W, rho0, r_const
    real(kind=rp), intent(in), dimension(nfr) :: fband, AA 
    real(kind=rp), intent(in), dimension(nth) :: theta ! directivity
    real(kind=rp), dimension(nfr,nth) :: getLandingGear
    !
    integer, parameter :: np = 3
    real(kind=rp), dimension(np), parameter :: St0 = (/one, 0.3_rp, 0.1_rp/)   ! The three spectra peaks
    real(kind=rp), dimension(np), parameter :: sigma_e = (/4.0_rp, 3.0_rp, 2.0_rp/) ! Empirical values
    real(kind=rp), dimension(np), parameter :: q = (/2.6_rp, 4.2_rp, 4.2_rp/) ! Empirical values
    real(kind=rp), dimension(np), parameter :: my = (/2.5_rp, 1.5_rp, 1.1_rp/) ! Empirical values
    real(kind=rp), dimension(np), parameter :: A = ((q*my/sigma_e)**q)*(St0**(q*my-sigma_e)) ! Analytical values based on empirical values
    real(kind=rp), dimension(np), parameter :: B = (my*q/sigma_e-one)*St0**my ! Analytical value based on empirical values
    real(kind=rp), dimension(np), parameter :: h = (/0.2_rp, 0.6_rp, 1.0_rp/)
    real(kind=rp), dimension(np), parameter :: beta = (/4.5E-8_rp,1.5E-8_rp,3.2E-5_rp/) ! empirical amplitudes of the three landing gear noise comp.
    !
    real(kind=rp), dimension(nthet) :: Do 
    real(kind=rp) :: L, U, Sl, Sm, am, eta
    real(kind=rp), dimension(max(nlsm_airfrm,nlsn_airfrm)) :: sj 
    real(kind=rp), dimension(np) :: S
    real(kind=rp), dimension(nfreq,np) :: St, F
    real(kind=rp), dimension(nthet,np) :: D
    real(kind=rp), dimension(nthet,nfreq,3) :: P    
    real(kind=rp), dimension(nthet,nfreq) :: prms_check    
    real(kind=rp) :: fact1, fact2, fact3
    integer :: i, j
    !
    L = sum(Lj(1:nls))
    U = M*c0
    
! create perimeter of the cross section of the j.th strut.
    do i = 1, nls
        if( geom_strutn(i).eq.zero ) then
          sj(i)=pi_num*d_strut(i)            
        else
          sj(i)=two*(d_strut(i) + geom_strutn(i))
        end if
    end do

! Installation Effects
    Do=1.2_rp*(one-0.9_rp*(cos(theta*pi_num/180.0_rp))**2) ! Takes care of the refraction/reflection/diffraction from the airplane (Empirical)

! Low frequency        
    Sl = pi_num*real(Nw,rp)*width_wheel*d_wheel

! Medium frequency
    Sm = sum(sj(1:nls)*Lj(1:nls))
    am = Sm/(pi_num*L)
    
! High frequency
    eta = (one + 0.028_rp*((Nw*L*(W/0.45359237_rp))/(Nref_airfrm*Lref_airfrm*Wref_airfrm)-one))*(one + two*((Nw-two)/Nw)*sin(two*psi)) ! Complexity value of the small details
    l = 0.15*am ! length-scale of the small details
    
    S=(/pi_num*Nw*width_wheel*d_wheel, sum(sj(1:nls)*Lj(1:nls)), eta*l**2/)
!    
    St(:,1)= fband(:)*(d_wheel/U)
    St(:,2)= fband(:)*(am/U)
    St(:,3)= fband(:)*(l/U)
!
    F(:,1) = (A(1)*((St(:,1)**sigma_e(1))/((B(1)+St(:,1)**my(1))**q(1))))
    F(:,2) = (A(2)*((St(:,2)**sigma_e(2))/((B(2)+St(:,2)**my(2))**q(2))))
    F(:,3) = (A(3)*((St(:,3)**sigma_e(3))/((B(3)+St(:,3)**my(3))**q(3))))
!
    D(:,1)=((one+h(1)*(cos(theta*pi_num/180.0_rp))**2)**2)
    D(:,2)=((one+h(2)*(cos(theta*pi_num/180.0_rp))**2)**2)
    D(:,3)=((one+h(3)*(cos(theta*pi_num/180.0_rp))**2)**2)
!
    do i = 1, nth
        do j = 1, nfr
          P(i,j,1)=beta(1)*S(1)*D(i,1)*F(j,1)
          P(i,j,2)=beta(2)*S(2)*D(i,2)*F(j,2)
          P(i,j,3)=beta(3)*S(3)*D(i,3)*F(j,3)
        end do
    end do
    
!    call save_3D_matrix(nth,nfr,3,P,'.m') 
    
    getLandingGear = zero
    do i = 1, nthet
        do j = 1, nfreq
           fact1 = ((((rho0*c0**2)**2)*exp(-(AA(j)/100.0_rp)*r_const)*(M**6))*Do(i))
           fact2 = (r_const**2)*(one-M*cos(theta(i)*pi_num/180.0_rp))**4
           fact3 = (P(i,j,1)+P(i,j,2)+P(i,j,3))
           getLandingGear(j,i)=(fact1/fact2)*fact3 ! do the transpose here to avoid having SPL:s defined on non-standard form
!           prms_check(i,j) = (fact1/fact2)*fact3 
        end do
    end do
!    call save_matrix(nthet,nfreq,prms_check,'.m') 
    
    getLandingGear = sqrt(getLandingGear)        
    
  end function getLandingGear
  !        
  !
  !
  function get_directivity(nthet,nfreq,Va,fband,theta,SPL)
    integer, intent(in) :: nthet, nfreq
    real(kind=rp), intent(in) :: Va
    real(kind=rp), intent(in), dimension(nfreq) :: fband
    real(kind=rp), intent(in), dimension(nthet) :: theta
    real(kind=rp), intent(in), dimension(nfreq) :: SPL
    real(kind=rp), dimension(nfreq,nthet) :: get_directivity
    !
    real(kind=rp) :: D
    real(kind=rp), dimension(nfreq) :: St, X
    integer, dimension(nfreq) :: maxPos
    integer, parameter :: ndata = 11
    real(kind=rp), dimension(ndata) :: angles=(/25.0_rp,32.0_rp,39.0_rp,48.0_rp,57.0_rp,68.0_rp,80.0_rp,90.0_rp,106.0_rp,120.0_rp,133.0_rp/)
    real(kind=rp), dimension(ndata,3) :: constants_data=(/-0.3897_rp,-0.4547_rp,-0.4792_rp,-0.5964_rp,-0.5097_rp,-0.4392_rp,-0.3839_rp,0.0000_rp,0.1552_rp,0.1565_rp,-0.0426_rp, & 
       & -0.5585_rp,-0.4823_rp,-0.2639_rp,-0.0682_rp,-0.0378_rp,-0.0494_rp,-0.0588_rp,0.0000_rp,-0.0665_rp,-0.1438_rp,-0.7227_rp, & 
       & -2.13_rp,-1.16_rp,-0.40_rp,0.080_rp,0.470_rp,0.400_rp,-0.34_rp,0.0000_rp,-0.85_rp,-1.39_rp,-2.01_rp/)
    real(kind=rp), dimension(nthet,3) :: const
    real(kind=rp), dimension(nfreq,3) :: variable
    integer :: i, j, iptr1, iptr2, idump
   
    maxPos = maxloc(SPL)
    D=Va/fband(maxPos(1))
    St=fband*D/Va
    X=log10(St)

    iptr1 = findFirstElementGtThanBnd(angles(1),nthet,theta)
    iptr2 = findFirstElementGtThanBnd(angles(ndata),nthet,theta)-1
    
    do i = 1,3
       do j = iptr1, iptr2
          call linear_interp(get_gestpan_1d_table(ndata,angles,constants_data(:,i)),1,theta(j),zero,const(j,i),idump)
       end do
       const(1:iptr1-1,i) = constants_data(1,i)
       const(iptr2+1:nthet,i) = constants_data(ndata,i)
    end do
        
    variable(1:nfreq,1) = X(:)**2
    variable(1:nfreq,2) = X(:)
    variable(1:nfreq,3) = one
    
    do i = 1, nfreq
        do j = 1, nthet
            get_directivity(i,j)=SPL(i) + sum(const(j,:)*variable(i,:))
        end do
    end do
      
  end function get_directivity
  !
  !
  !
  function get_slat(nfr,Z,M,C1,alfa,xsi,t,Va,fband)
     integer, intent(in) :: nfr
     real(kind=rp), intent(in) :: Z, M,C1,alfa,xsi,t,Va
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), dimension(nfreq) :: get_slat
    !
    real(kind=rp), dimension(3) :: c = (/7.6_rp,-4.1355_rp,1.7486_rp/)
    real(kind=rp) :: n, OASPL 
    real(kind=rp), dimension(nfreq) :: X, Xones
    real(kind=rp), dimension(nfreq,7) :: lmat
    real(kind=rp), dimension(7) :: cvec = (/-12.7289_rp,4.05419_rp,-4.3004_rp,-1.8453_rp,1.9039E-1_rp,2.1658E-1_rp,-3.9728E-1_rp/)
    real(kind=rp), dimension(nfreq,7) :: constants
    real(kind=rp), dimension(nfreq) :: SPL
    integer :: i 
    !    
    n = 10.0_rp*(c(3)*Z**2 + c(2)*Z + c(1))
    !
    OASPL = 117.4_rp + n*log10(M) + 10.3_rp*log10(sin(alfa)) + 44.5_rp*log10(C1) + 11.5_rp*log10(xsi*t)
    
! lmat
    X=log10(fband(:)*xsi*t/Va)    
    lmat(1:nfreq,1) = one
    do i = 2, 7
        lmat(:,i) = X(:)**(i-1)
    end do    
! constants
    do i = 1, nfreq
       constants(i,:) = cvec(:)
    end do
!    
    SPL=sum(lmat*constants,2)+OASPL
    get_slat = P0*10.0_rp**(SPL/20.0_rp)
     
  end function get_slat
  !
  !
  !
  function get_outboard_flap(nfr,Z,M,C1,delta,Vy,xsi,t,Va,Circ,fband)
    integer, intent(in) :: nfr
    real(kind=rp), intent(in) :: Z, M,C1,delta,Vy,xsi,t,Va,Circ
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), dimension(nfreq) :: get_outboard_flap
    !
    real(kind=rp), dimension(nfreq) :: X, Xones
    real(kind=rp), dimension(nfreq,7) :: lmat
    real(kind=rp) :: OASPL
    real(kind=rp), dimension(4) :: bvec = (/4.1204_rp,-6.220_rp,5.5040_rp,-1.6245_rp/)
    real(kind=rp), dimension(7) :: cvec = (/-14.5121_rp,9.1991_rp,-4.7859_rp,-3.8899_rp,1.3238_rp,4.1053E-1_rp,-5.8664E-1_rp/)
    real(kind=rp), dimension(nfreq,7) :: constants
    real(kind=rp), dimension(nfreq) :: SPL
    !
    real(kind=rp) :: msmall
    integer :: i 

    msmall=10.0_rp*(bvec(4)*Z*3+bvec(3)*Z*2+bvec(2)*Z+bvec(1))
    
    OASPL = 93.3_rp +53.0_rp*log10(M) + 20.15_rp*log10(sin(delta)) + 0.5_rp*log10(C1) + msmall*log10(xsi*t) + 0.11_rp*log10(Vy) + 0.12_rp*log10(Circ)

! lmat
    X=log10(fband(:)*xsi*t/Va)    
    lmat(1:nfreq,1) = one
    do i = 2, 7
        lmat(:,i) = X(:)**(i-1)
    end do    
! constants
    do i = 1, nfreq
       constants(i,:) = cvec(:)
    end do
!    
    SPL=sum(lmat*constants,2)+OASPL
    get_outboard_flap = P0*10.0_rp**(SPL/20.0_rp)

  end function get_outboard_flap
  !
  !
  !
  function get_inboard_flap(nfr,M,C1,delta,Vy,xsi,t,Va,Circ,fband)
    integer, intent(in) :: nfr
    real(kind=rp), intent(in) :: M,C1,delta,Vy,xsi,t,Va,Circ
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), dimension(nfreq) :: get_inboard_flap
    !
    real(kind=rp), dimension(nfreq) :: X, Xones
    real(kind=rp), dimension(nfreq,7) :: lmat
    real(kind=rp) :: OASPL
    real(kind=rp), dimension(7) :: cvec = (/-13.1814_rp,10.0319_rp,-9.4077_rp,-6.6166_rp,4.0970_rp,1.4162_rp,-1.0740_rp/)
    real(kind=rp), dimension(nfreq,7) :: constants
    real(kind=rp), dimension(nfreq) :: SPL
    !
    integer :: i 

    
!---    
    OASPL = 92.1_rp +53.0_rp*log10(M) + 20.15_rp*log10(sin(delta)) + 0.5_rp*log10(C1) + 17.8_rp*log10(xsi*t) + 0.11_rp*log10(Vy) + 0.12_rp*log10(Circ)

! lmat
    X=log10(fband(:)*xsi*t/Va)    
    lmat(1:nfreq,1) = one
    do i = 2, 7
        lmat(:,i) = X(:)**(i-1)
    end do
!
    do i = 1, nfreq
       constants(i,:) = cvec(:)
    end do
!    
    SPL=sum(lmat*constants,2)+OASPL
    get_inboard_flap = P0*10.0_rp**(SPL/20.0_rp)

  end function get_inboard_flap 
  !
  !
  !
  function get_aileron(nfr,Z,b,M,C1,delta,alfa,xsi,t,Va,fband)
    integer, intent(in) :: nfr
    real(kind=rp), intent(in) :: Z,b,M,C1,delta,alfa,xsi,t,Va
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), dimension(nfreq) :: get_aileron
    !
    real(kind=rp) :: n, OASPL
    !
    real(kind=rp), dimension(nfreq) :: X, Xones
    real(kind=rp), dimension(nfreq,7) :: lmat
    real(kind=rp), dimension(7) :: cvec = (/-14.0467_rp,9.0187_rp,-6.4007_rp,-4.3267_rp,2.1149_rp,5.725E-1_rp,-7.062E-1_rp/)
    real(kind=rp), dimension(nfreq,7) :: constants
    real(kind=rp), dimension(nfreq) :: SPL
    integer :: i 
    
    n = 52.0_rp - 5.8565_rp*(Z-one)
    
    !OASPL for 90 Degrees
    OASPL = 109.7_rp + n*log10(M) + 3.36_rp*log10(C1*sin(delta)) + 3.46_rp*log10(sin(alfa)) + 3.46_rp*log10(xsi*t)
  
    X=log10(fband(:)*xsi*t/Va)    
    lmat(1:nfreq,1) = one
    do i = 2, 7
        lmat(:,i) = X(:)**(i-1)
    end do
!
    do i = 1, nfreq
       constants(i,:) = cvec(:)
    end do
!
    
    SPL=sum(lmat*constants,2)+OASPL
    get_aileron=P0*10.0_rp**(SPL/20.0_rp)
    
  end function get_aileron
  !
  ! calcCaj: Statino 1 generally refers to the core flow, whereas station 2 refers to the bypass flow.
  !
  ! Based on NASA CR-3786, Russel, J. W., "An empirical Method for Predicting
  ! the mixing Noise Levels of Subsonic Circular and Coaxial Jets", 1984
  !
  ! Parameters:
  !           A_e - Nozzle exit flow area of single equivalent jet [m^2]
  !           A_ref - Reference area used in computing normalized overall power level [m^2]
  !           A_1 - Nozzle exit flow area of inner stream or circular jet [m^2]
  !           A_2 - Nozzle exit flow area of outer stream [m^2]
  !           ALT - Flight altitude [m]
  !           c_amb - Ambient speed of sound [m/s]
  !           c_ISA - Speed of sound at ISA SLS conditions [m/s]
  !           D1...D36 - Derivative terms
  !           D(theta_c) - Directivity index at coordinate point [dB]
  !           D_e - Nozzle exit flow diameter of single equivalent jet [m]
  !           D_j(theta_c) - Derivative value for directivity index at coordinate point [-] 
  !           dmdt_e - mass flow rate of single equivalent jet [kg/s] 
  !           dmdt_1 - mass flow rate of inner stream or circular jet [kg/s]
  !           dmdt_2 - mass flow rate of outer stream [kg/s]
  !           eta_c - Normalized frequency parameter at coordinate point [-]
  !           eta_p - Normalized frequency parameter at predicetd point [-]
  !           F(eta_c) - 1/3 octave band normalized power spectrum at coordinate point [dB] 
  !           F(eta_p) - 1/3 octave band normalized power spectrum at predicted point [dB]
  !           F_j - Derivative value of 1/3 ocatave band normalized power
  !                 spectrum at coordinate point [-]
  !           f_p - 1/3 octave band predicted center frequency [Hz]
  !           gamma_e - ratio of specific heats [-] 
  !           gamma_1 - ratio of spec. heats of inner stream or circular jet [-]
  !           gamma_2 - ratio of spec. heats of outer stream [-]
  !           OAPWL_norm - Normalized overall acoustic power level [dB]
  !           OASPL(theta_c) - Overall sound pressure level at directivity coordinate point [dB]
  !           OASPL(theta_p) - Overall sound pressure level at predicted directivity angle [dB]
  !           PWL_j - Derivative values for OAPWL [-]
  !           R - gas constant [J/kg*K]
  !           r - Radial distance from nozzle exit to observer [m]
  !           rho_amb - Ambient air density [kg/m^3]
  !           rho_e - Nozzle exit flow density of single equivalent jet [kg/m^3]
  !           rho_ISA - Density of air under ISA SLS conditions [kg/m^3]
  !           RSL(theta_c,eta_c) - Normalized relative spectrum level at coordinate point [dB]
  !           RSL(theta_p,eta_p) - Normalized relative spectrum level at
  !                                predicted directivity angle and predicted normalized frequency prarameter [dB]
  !           RSL_j - Derivative values of normalized relative spectrum level
  !                   at coordinate point [dB]
  !           SPL(theta_p,f_p) - SOund pressure level at the preicted directivity angle and predicted frequency [dB]
  !           v_e - Nozzle exit equivalent flow velocity [m/s]
  !           t_amb - Ambient static temperature [K]
  !           T_e - Nozzle exit flow equivalent jet total temp [K]
  !           T_1 - Nozzle exit flow total temp. of inner stream or circular jet [K]
  !           T_2 - Nozzle exit flow total temp. of outer stream [K]
  !           theta_c - Directivity angle at coordinate point [deg]
  !           theta_p - Directivity angle at predicted point [deg]
  !           v_1 - Nozzle exit flow velocity of inner stream or circular jet [m/s]
  !           v_2 - Nozzle exit flow velocity of outer stream [m/s]
  !     x_1...x_5 - Equivalent flow parameters [-]  
  !
  ! Limitations:
  !           0.3 <= v_e/c_amb <= 2.0
  !           0.7 <= T_e/t_amb <= 4.5
  !          0.02 <= v_2/v_1 <= 2.5
  !           0.2 <= T_2/T_1 <= 4.0
  !           0.5 <= A_2/A_1 <= 10.0
  !
  subroutine calcCaj(cnt,dmdt_1,dmdt_2,v_1,v_2,T_1,T_2,gamma_1, gamma_2, Ta, Pa, A_1, A_2, r_const, nfr,nth,fband,theta,prms)
    integer, intent(in) :: cnt ! loop counter - useful for controlling debugging printouts 
    integer, intent(in) :: nfr, nth
    real(kind=rp), intent(in) :: dmdt_1, dmdt_2, v_1, v_2, T_1, T_2, gamma_1, gamma_2, Ta, Pa, A_1, A_2, r_const
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), intent(in), dimension(nth) :: theta
    real(kind=rp), intent(out), dimension(nfr,nth) :: prms
    !
    real(kind=rp) :: dmdt_e, v_e, gamma_e_gamma_e_1, c_amb, c_isa, gamma_e, rho_amb, rho_e, rho_isa, & 
    & t_amb, t_e, a_e, d_e, x_1, x_2, x_3, x_4, x_5,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,D21,& 
         & D22,D23,D24,D25,D26,D27,D28,D29,D30,D31,D32,D33,D34,D35,D36,OAPWL_norm, A_ref
    real(kind=rp), dimension(36) :: derivMultipliers
    integer, parameter :: n_thetas = 7
    real(kind=rp), dimension(n_thetas) :: theta_c = (/0.0_rp,30.0_rp,60.0_rp,90.0_rp,120.0_rp,150.0_rp,180.0_rp/)
    real(kind=rp), dimension(n_thetas) :: eta_c = (/-1.5_rp,-1.0_rp,-0.5_rp,0.0_rp,0.5_rp,1.0_rp,1.5_rp/)
    real(kind=rp), dimension(n_thetas) :: D, F, OASPL, derivative
    real(kind=rp), dimension(n_thetas,n_thetas) :: RSL 
    real(kind=rp), dimension(nthet) :: OASPL_p, OASPL_p_der
    real(kind=rp), dimension(nfreq) :: eta_p, F_p, F_p_der
    real(kind=rp), dimension(nfreq,nthet) :: SPL, SPL_temp
    integer, parameter :: lwk = 2*n_thetas
    real(kind=rp),  dimension(lwk) :: wk
    logical :: spline, skip
    !
    integer :: i, j, k, ierr, switch, incfd, idump
    integer, parameter :: IBEG = 1
    integer, parameter :: IEND = 1
    integer, dimension(2), parameter :: IC = (/IBEG,IEND/)
    real(kind=rp), dimension(2), parameter :: VC = (/zero,zero/)
    real(kind=rp), dimension(2+n_thetas+n_thetas+n_thetas*n_thetas) :: gpvec
    real(kind=rp), dimension(nfreq,nthet) :: RSL_p
    real(kind=rp) :: zval

    !    IC(1) = IBEG, desired condition at beginning of data.
!    IC(2) = IEND, desired condition at end of data.
    
    !
    dmdt_e = dmdt_1 + dmdt_2 ! equivalent mass flow
    v_e = (dmdt_1*v_1 + dmdt_2*v_2)/dmdt_e ! equivalent velocity
    T_e = (dmdt_1*((gamma_1)/(gamma_1 - one))*T_1 + dmdt_2*((gamma_2)/(gamma_2 - one))*T_2)/(dmdt_1*((gamma_1)/(gamma_1 - one)) + dmdt_2*((gamma_2)/(gamma_2 - one))) ! Compute equivalent temperature
    gamma_e_gamma_e_1 = (dmdt_1*((gamma_1)/(gamma_1 - one)) + dmdt_2*((gamma_2)/(gamma_2 - 1.0)))/(dmdt_1 + dmdt_2) ! Compute the equivalent ratio of spec. heat:
    gamma_e = gamma_e_gamma_e_1/(gamma_e_gamma_e_1-one)
    
    rho_amb = pa/(R_air*Ta)
    c_amb = sqrt(gamma_2*R_air*Ta)
    t_amb = Ta
    rho_ISA = 101325.0_rp/(R_air*288.15_rp)
    c_ISA = sqrt(gamma_2*R_air*288.15_rp)

    rho_e = rho_amb*(T_e/t_amb - ((gamma_e-one)/two) * (v_e/c_amb)**two)**(-one)
    
    A_e = dmdt_e/(rho_e*v_e)     ! compute the equivalent jet area
    D_e = sqrt(four*A_e/pi_num)  ! equivalent diameter
    
! turbofan case
    x_1 = log10((v_e/c_amb)/one) 
    x_2 = log10((T_e/t_amb)/two)
    x_3 = log10((v_2/v_1)/one)
    x_4 = log10((T_2/T_1)/one)
    x_5 = log10((A_2/A_1)/one)
    
    if (v_e/c_amb.lt.0.30_rp .or. v_e/c_amb.gt. 2.0_rp ) call error_proc(200,'v_e/c_amb is not within valid range (0.3 - 2.0)',    'call from calcCaj')        
    if (T_e/t_amb.lt.0.70_rp .or. T_e/t_amb.gt. 4.5_rp ) call error_proc(200,'T_e/t_amb is not within the valid range (0.7 - 4.5)','call from calcCaj')        
    if (v_2/v_1  .lt.0.02_rp .or.  v_2/v_1 .gt. 2.5_rp ) call error_proc(200,'v_2/v_1 is not within the valid range (0.02 - 2.5)', 'call from calcCaj')        
    if (T_2/T_1  .lt.0.20_rp .or.  T_2/T_1 .gt. 4.0_rp ) call error_proc(200,'T_2/T_1 is not within the valid range (0.2 - 4.0)',  'call from calcCaj')        
    if (A_2/A_1  .lt.0.50_rp .or.  A_2/A_1 .gt.10.0_rp ) call error_proc(200,'A_2/A_1 is not within the valid range (0.5 - 10.0)', 'call from calcCaj')
    
! STEP 3 - Compute Derivative Multipliers 
! Circular and coannular jet
    D1 = one
    D2 = x_1
    D3 = x_2
    D4 = (x_1**two)/two
    D5 = x_1*x_2
    D6 = (x_2**two)/two
    D7 = (x_1**two*x_2)/two
    D8 = (x_1*x_2**two)/two
! Coannular jet only:
    D9 = x_3
    D10 = x_4
    D11 = x_1*x_3
    D12 = x_1*x_4
    D13 = x_1*x_5
    D14 = x_2*x_3
    D15 = x_2*x_4
    D16 = x_2*x_5
    D17 = (x_3**two)/two
    D18 = x_3*x_4
    D19 = x_3*x_5
    D20 = (x_4**two)/two
    D21 = x_4*x_5
    D22 = (x_1**two*x_3)/two
    D23 = (x_1**two*x_5)/two
    D24 = (x_1*x_3**two)/two
    D25 = x_1*x_3*x_4
    D26 = x_1*x_3*x_5
    D27 = x_1*x_4*x_5
    D28 = (x_1*x_5**two)/two
    D29 = (x_3**3.0)/6.0
    D30 = (x_3**two*x_4)/two
    D31 = (x_3**two*x_5)/two
    D32 = (x_3*x_4**two)/two
    D33 = x_3*x_4*x_5
    D34 = (x_3*x_5**two)/two
    D35 = (x_4**two*x_5)/two
    D36 = (x_4*x_5**two)/two

    derivMultipliers = (/D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,D21,& 
         & D22,D23,D24,D25,D26,D27,D28,D29,D30,D31,D32,D33,D34,D35,D36/)
    
    !
    OAPWL_norm = zero
    do i = 1, 36
       OAPWL_norm = OAPWL_norm + Table_IV(i,1)*DerivMultipliers(i)
    end do
    
    D = zero ; F = zero
    do j = 1, n_thetas
       do i = 1, 36
          D(j) = D(j) + Table_IV(i,j+1)*DerivMultipliers(i)
          F(j) = F(j) + Table_V(i,j)*DerivMultipliers(i)
       end do 
    end do
    
    RSL = zero
    do k = 1, n_thetas
       do j = 1, n_thetas
           do i = 1, 36
              RSL(j,k) = RSL(j,k) + Table_VI(i,j,k)*DerivMultipliers(i) 
           end do
       end do
    end do   
!    call save_matrix(n_thetas,n_thetas,RSL) - RSL was verified. Very accurate. 
    
    A_ref = dmdt_e/(rho_amb*c_amb)
    
    do i = 1, n_thetas
       OASPL(i) = OAPWL_norm + D(i) + 20.0_rp*log10((rho_amb*c_amb**two)/(rho_ISA*c_ISA**2.0)) + & 
           &  10.0_rp*log10(A_ref/(4.0_rp*pi_num*r_const**two)) + 197.0_rp
    end do

! STEP 9 - Compute the predicted OASPL 
! gives about a 0.3% error from noPred interpolation
    switch = 0 ; incfd = 1 ; skip = .false.
    call pchic ( ic, vc, zero, n_thetas, theta_c, OASPL, derivative, incfd, wk, lwk, ierr )
    call pchfe ( n_thetas, theta_c, OASPL, derivative, incfd, skip, nthet, theta, OASPL_p, ierr )
    
!-------------------- FULINTERPOLATION  
!    spline = .false.
!    call pchez ( n_thetas, theta_c, OASPL, derivative, spline, wk, lwk, ierr )
!    call pchev ( n_thetas, theta_c, OASPL, derivative, nthet, theta, OASPL_p, OASPL_p_der, ierr )
!    
! 
    eta_p = log10((fband*D_e)/v_e)
    gpvec = get_gestpan_2d_table(n_thetas,theta_c,n_thetas,eta_c,transpose(RSL))
    !
    do i = 1, nfreq
       do j = 1, nthet
          call linear_interp(gpvec,1,theta(j),eta_p(i),RSL_p(i,j),idump)
       end do
    end do
    
! STEP 12 - Compute the norm frequency params 
!    switch = 0 ; incfd = 1 ; skip = .false.
!    call pchic ( ic, vc, zero, n_thetas, eta_c, F, derivative, incfd, wk, lwk, ierr )
!    call pchfe ( n_thetas, eta_c, F, derivative, incfd, skip, nfreq, eta_p, F_p, ierr )

! interpolate F_p
    spline = .true.     ! spline shape preserving interpolation - ala MATLAB
    call pchez ( n_thetas, eta_c, F, derivative, spline, wk, lwk, ierr )
    call pchev ( n_thetas, eta_c, F, derivative, nfreq, eta_p, F_p, F_p_der, ierr )

!-------------------- FULINTERPOLATION  
!    spline = .false.
!    call pchez ( n_thetas, eta_c, F, derivative, spline, wk, lwk, ierr )
!    call pchev ( n_thetas, eta_c, F, derivative, nfreq, eta_p, F_p, F_p_der, ierr )

!---- STEP 13 - Compute the norm frequency params 
    do i = 1, nthet
       do j = 1, nfreq
          SPL_temp(j,i) = RSL_p(j,i) + OASPL_p(i) 
       end do
    end do
    do i = 1, nfreq
       do j = 1, nthet
          SPL(i,j) = SPL_temp(i,j) + F_p(i) 
       end do
    end do
    Prms = P0*10.0_rp**(SPL/20.0_rp)
    
  end subroutine calcCaj
  !
  !-------------------- code use to generate parameter arrays
!   call parseTable('Table_IV.m')
!   call parseTable('Table_V.m')
!   call parseTable('Table_VII.m')
!   call parseTable('Table_VI2.m')
!   call parseTable('Table_VI3.m')
!   call parseTable('Table_VI4.m')
!   call parseTable('Table_VI5.m')
!   call parseTable('Table_VI6.m')
!   call parseTable('Table_VI7.m')
!   stop
  !
  !
  subroutine calcComb(cnt,testCalcTurbine,cmbtype,Nf_comb,R0,Aec,De,Dh,Lc,h,Nf,Nfmax,pattern,pa,p3,p4,p7,ta,t3,t4,t5,w3,nfrec,nt,theta,fband,f,freq,prms)
    integer, intent(in) :: cnt ! loop counter - useful for controlling debugging printouts 
    logical, intent(in) :: testCalcTurbine ! if true => load hard coded data to allow comparison with NoPred. False => run on usual input
    character(len=*), intent(in) :: cmbtype
    real(kind=rp), intent(in) :: Nf_comb,R0,Aec,De,Dh,Lc,h,Nf,Nfmax,pattern
    real(kind=rp), intent(in) :: pa,p3,p4,p7
    real(kind=rp), intent(in) :: ta,t3,t4,t5
    real(kind=rp), intent(in) :: w3
    integer, intent(in) :: nfrec ! number of frequencies
    integer, intent(in) :: nt ! number of theta 
    real(kind=rp), intent(in), dimension(*) :: theta ! directivity
    real(kind=rp), intent(in), dimension(nfreq) :: fband
    real(kind=rp), intent(in), dimension(nfreq+1) :: f
    real(kind=rp), intent(in), dimension(*) :: freq
    real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms
    !
    real(kind=rp) :: pa_psia, p3_psia, p4_psia, p7_psia, t3_fahr, t4_fahr, t5_fahr
    real(kind=rp), dimension(nfreq) :: AA
    real(kind=rp), dimension(nfreq,nthet) :: SPL
    !
    integer :: slen, i 
    !
    pa_psia = pa*pa2psia ; p3_psia = p3*pa2psia ; p4_psia = p4*pa2psia ; p7_psia = p7*pa2psia   ! pressures in psia
    t3_fahr = 1.8_rp*t3 + 32.0_rp ; t4_fahr = 1.8_rp*t4 + 32.0_rp ; t5_fahr = 1.8_rp*t5 + 32.0_rp ! temperatures in Fahrenheit
    
    AA =  get_atm_abs(nfreq,Ta,RH,fband) ! dB/100 m
    AA(1:nfreq) = (AA(1:nfreq)/100.0_rp)*0.3048_rp*1000.0_rp

    slen = len(cmbtype) 
    if( cmbtype(1:slen).eq.'SAC' ) then
        SPL(1:nfreq,1:nthet) = getSAC(nfreq,nthet,R0,Aec,De,Dh,Lc,h,Nf,Nfmax,Pa_psia,P3_psia,P4_psia,P7_psia,T3_fahr,T4_fahr,T5_fahr,W3,AA,theta,fband)
    elseif( cmbtype(1:slen).eq.'DAC' ) then
       SPL(1:nfreq,1:nthet) = getDAC(nfreq,nthet,R0,Aec,De,Dh,Lc,h,Nf,Nfmax,pattern,Pa_psia,P3_psia,P4_psia,P7_psia,T3_fahr,T4_fahr,T5_fahr,W3,AA,theta,fband)
    else
        stop 'combtype not recognised in calcComb'
    end if
    
    ! Fix SPL to radius 1 meter
    do i=1, nthet
      SPL(:,i)=SPL(:,i)-AA(:)*(one/(0.3048_rp*1000.0_rp))  
    end do 
    
    Prms = P0*10.0_rp**(SPL/20.0_rp)
  
  end subroutine calcComb
  !
  ! getSAC: calculate the noise levels from a Single Annular Combustor
  !
  function getSAC(nfr,nth,R0,Aec,De,Dh,Lc,h,Nf,Nfmax,Pa,P3,P4,P7,T3,T4,T5,W3,AA,theta,fband)
     integer, intent(in) :: nfr, nth
     real(kind=rp), intent(in) :: R0,Aec,De,Dh,Lc,h,Nf,Nfmax,Pa,P3,P4,P7,T3,T4,T5,W3
     real(kind=rp), intent(in), dimension(nfr) :: AA,fband
     real(kind=rp), intent(in), dimension(nth) :: theta
     real(kind=rp), dimension(nfr,nth) :: getSAC
     integer, parameter :: nang = 3
     real(kind=rp), dimension(nang) :: fp =  (/63.0_rp,160.0_rp,630.0_rp/) ! step 1
     real(kind=rp), dimension(nang) :: ap = (/150.0_rp,130.0_rp,130.0_rp/) ! step 1
     !
     real(kind=rp), dimension(nfreq,3) :: fN, SPLN
     real(kind=rp), dimension(nthet,3) :: aN, OASPLN, OASPL
     real(kind=rp), dimension(nfreq,nthet,3) :: SPL, prms
     real(kind=rp), dimension(nfreq,nthet) :: prms_total
     real(kind=rp) :: CP, Fc, Ft, Tl, SPL_Fc, SPL_Tl
     real(kind=rp), dimension(nang) :: Hcp
     integer :: iptr1,iptr2,iptr3     
     integer :: i, j
     
! STEP 2 - normalized frequency
     fN(1:nfr,1) =fband(1:nfr)/fP(1) ; fN(1:nfr,2) =fband(1:nfr)/fP(2) ;  fN(1:nfr,3) =fband(1:nfr)/fP(3) 
     aN(1:nth,1) =theta(1:nth)/aP(1) ; aN(1:nth,2) =theta(1:nth)/aP(2) ;  aN(1:nth,3) =theta(1:nth)/aP(3) 

! STEP 3
     OASPLN(1:nth,1) =   -67.8_rp*(aN(1:nth,1))**2 +  141.7_rp*aN(1:nth,1) - 66.84_rp
     OASPLN(1:nth,2) = -26.019_rp*(aN(1:nth,2))**3 - 5.2974_rp*aN(1:nth,2)**2 + 93.43_rp*aN(1:nth,2) - 61.75_rp 
     OASPLN(1:nth,3) =  -156.5_rp*(aN(1:nth,3))**2 + 322.34_rp*aN(1:nth,3) - 164.89_rp

!     call save_matrix(nthet,3,OASPLN,'.m') ; stop
     
! STEP 5
     CP = (W3*sqrt(T3)/P3)*((T4-T3)/T4)*(P3/Pa)*((Dh/De)**0.5_rp)
     Hcp = (/76.45_rp+14.256_rp*log10(CP),108.5_rp+3.31_rp*log10(CP),106.38_rp+6.938_rp*log10(CP)/)

     OASPL=zero

     iptr1 = findFirstElementGeThanBnd(aP(1),nth,theta)
     OASPL(iptr1,1) = -20.0_rp*log(R0) + Hcp(1)*((30.0_rp/Nf)**(-0.225_rp))

     iptr2 = findFirstElementGeThanBnd(aP(2),nth,theta)
     OASPL(iptr2,2) = -20.0_rp*log(R0) + Hcp(2)*((30.0_rp/Nf)**(0.05_rp))

     iptr3 = findFirstElementGeThanBnd(aP(3),nth,theta)
     OASPL(iptr3,3) = -20.0_rp*log(R0) + Hcp(3)*((30.0_rp/Nf)**(0.02_rp))

!     call save_3D_matrix(nfreq,nthet,3,SPL,'.m') ; stop
     
     
! STEP 4
! Loss parameter
     Fc=(W3*sqrt(T4-T3)/(P3*Aec**2*sqrt(Nf)))
     Ft=(P4/P7)*sqrt(T5/T4)
     Tl=((one+Ft)**2)/(4.0_rp*Lc*Ft/(pi_num*h))

     SPL_Fc=20.0_rp*log10(Fc)
     SPL_Tl=20.0_rp*log10(Tl)

     OASPL(:,1) = OASPLN(:,1) + OASPL(iptr1,1) + 0.4_rp*(SPL_Fc-SPL_Tl)
     OASPL(:,2) = OASPLN(:,2) + OASPL(iptr2,2) + 0.1_rp*(SPL_Fc-SPL_Tl)
     OASPL(:,3) = OASPLN(:,3) + OASPL(iptr3,3) + 0.3_rp*(SPL_Fc-SPL_Tl)

!     call save_matrix(nthet,3,OASPL,'.m') ; stop
     
! STEP 6
     SPLN(:,1) = -152.70_rp + 295.46_rp*fN(:,1) - 145.61_rp*fN(:,1)**2 
     SPLN(:,2) = -170.07_rp + 331.33_rp*fN(:,2) - 163.34_rp*fN(:,2)**2 
     SPLN(:,3) = -147.50_rp + 286.40_rp*fN(:,3) - 142.31_rp*fN(:,3)**2

     do i = 1, nthet
         do j = 1, nfreq
           SPL(j,i,1) = OASPL(i,1) + SPLN(j,1) + AA(j)*(R0/1000.0_rp)
           SPL(j,i,2) = OASPL(i,2) + SPLN(j,2) + AA(j)*(R0/1000.0_rp)
           SPL(j,i,3) = OASPL(i,3) + SPLN(j,3) + AA(j)*(R0/1000.0_rp)
         end do
     end do
!     
!
     Prms = P0*10.0_rp**(SPL/20.0_rp)
     Prms_total(:,:) = sqrt(Prms(:,:,1)**2 + Prms(:,:,2)**2 + Prms(:,:,3)**2);
     
     WHERE( Prms_total.lt.P0 ) Prms_total=P0

     getSAC = 20.0_rp*log10(Prms_total/P0)
     
  end function getSAC
  !
  ! GETDAC calculates the noise from a dual annular combustor
  !
  function getDAC(nfr,nth,R0,Aec,De,Dh,Lc,h,Nf,Nfmax,pattern,Pa,P3,P4,P7,T3,T4,T5,W3,AA,theta,fband)
     integer, intent(in) :: nfr, nth
     real(kind=rp), intent(in) :: R0,Aec,De,Dh,Lc,h,Nf,Nfmax,pattern,Pa,P3,P4,P7,T3,T4,T5,W3
     real(kind=rp), intent(in), dimension(nfr) :: AA,fband
     real(kind=rp), intent(in), dimension(nth) :: theta
     real(kind=rp), dimension(nfr,nth) :: getDAC
     !
     integer, parameter :: n_fp = 2, n_ap = 1
     real(kind=rp), dimension(n_ap) :: ap = (/130.0_rp/) ! step 1
     real(kind=rp), dimension(n_fp) :: fp = (/160.0_rp,500.0_rp/) ! step 1
     real(kind=rp), dimension(2) :: Mf = (/0.020_rp,0.180_rp/)
     real(kind=rp), dimension(4,2) :: Knf = (/1.2_rp,0.98_rp,0.98_rp,1.10_rp,1.0_rp,0.90_rp,0.90_rp,0.98_rp/) ! fortran stores column by column
     real(kind=rp), dimension(4,2) :: Xk = (/0.25_rp,0.25_rp,0.2_rp,0.0_rp,0.25_rp,0.25_rp,0.2_rp,0.0_rp/) ! fortran stores column by column
     real(kind=rp), dimension(nfreq,n_fp) :: fN, SPLN
     real(kind=rp), dimension(nthet) :: aN
     real(kind=rp), dimension(nthet,2) :: OASPLN, OASPL
     real(kind=rp), dimension(2) :: CP
     real(kind=rp), dimension(4) :: Hcp
     real(kind=rp), dimension(nfreq,nthet,2) :: SPL, prms
     real(kind=rp), dimension(nfreq,nthet) :: prms_total
     integer :: ipat, iptr
     real(kind=rp) :: Fc, Ft, Tl, SPL_Fc, SPL_Tl
     integer :: i, j
     
!STEP 2 - normalized frequency
     fN(:,1)=fband(:)/fP(1)   
     fN(:,2)=fband(:)/fP(2)   
     aN(:)=theta(:)/aP(1)

! STEP 3
     OASPLN(:,1) = -116.95_rp*aN(:)**2 + 235.23_rp*aN(:) - 120.65_rp
     OASPLN(:,2) = -137.59_rp*aN(:)**2 + 283.40_rp*aN(:) - 147.73_rp
     
     !
     
! STEP 5
     CP(1:2)  = (/(W3*sqrt(T3)/P3)*((T4-T3)/T4)*(P3/Pa)*((Dh/De)**2),(W3*sqrt(T3)/P3)*((T4-T3)/T4)*(P3/Pa)*((Dh/De)**1.2)/)
     Hcp(1:4) = (/76.45_rp+14.256_rp*log10(CP(1)),76.45_rp+14.256_rp*log10(CP(2)),110.62_rp+2.997_rp*log10(CP(1)),110.62_rp+2.997_rp*log10(CP(2))/)
     !
     OASPL = zero
     !
     ipat = get_ipattern(pattern)
     
     iptr = findFirstElementGeThanBnd(aP(1),nth,theta)
     OASPL(iptr,1) = Knf(ipat,1)*(-20.0_rp*log(R0) + Hcp(1)*(((20.0_rp + Nf)/Nfmax)**(-Xk(ipat,1)))*((30.0_rp/(20.0_rp + Nf))**Mf(1))) ! eq 227
     OASPL(iptr,2) = Knf(ipat,2)*(-20.0_rp*log(R0) + Hcp(2)*(((20.0_rp + Nf)/Nfmax)**(-Xk(ipat,2)))*((30.0_rp/(20.0_rp + Nf))**Mf(2))) ! eq 227
     
     ! STEP 4
     ! Loss paramterer
     Fc=(W3*sqrt(T4-T3)/(P3*Aec**2*sqrt(20.0_rp+Nf)))
     Ft=(P4/P7)*sqrt(T5/T4)
     Tl=((one+Ft)**2)/(4.0_rp*Lc*Ft/(pi_num*h))
     
     !
     SPL_Fc=20.0_rp*log10(Fc)
     SPL_Tl=20.0_rp*log10(Tl)

     OASPL(:,1) = OASPLN(:,1) + OASPL(iptr,1) + 0.45_rp*(SPL_Fc-SPL_Tl)
     OASPL(:,2) = OASPLN(:,2) + OASPL(iptr,2) - 0.10_rp*(SPL_Fc-SPL_Tl)

! STEP 6
     SPLN(:,1) = -170.07_rp + 331.33_rp*fN(:,1) - 163.34_rp*fN(:,1)**2
     SPLN(:,2) = -137.21_rp + 268.99_rp*fN(:,2) - 135.81_rp*fN(:,2)**2

     !
     
     do i = 1, nthet
         do j = 1, nfreq
            SPL(j,i,1) = OASPL(i,1) + SPLN(j,1) + AA(j)*(R0/10000_rp)
            SPL(j,i,2) = OASPL(i,2) + SPLN(j,2) + AA(j)*(R0/10000_rp)
         end do
     end do
!     
!
     Prms = P0*10.0_rp**(SPL/20.0_rp)
     Prms_total(:,:) = sqrt(Prms(:,:,1)**2 + Prms(:,:,2)**2)
     
     WHERE( Prms_total.lt.P0 ) Prms_total=P0

     getDAC = 20.0_rp*log10(Prms_total/P0)
     
     
  end function getDAC
  !
  !
  !
  function get_ipattern(pat)
    real(kind=rp) :: pat
    integer :: get_ipattern
    
    if( pat .eq. 40_rp ) then 
      get_ipattern = 1
    elseif( pat .eq. 30_rp ) then 
      get_ipattern = 2
    elseif( pat .eq. 22.5_rp ) then 
      get_ipattern = 3
    else
      get_ipattern = 4
    end if
    
  end function get_ipattern
  !
  !
  !
  subroutine calcTurbine(cnt,testCalcTurbine,N_rotors,n_stages,Vref,cref,Vtr,Texit,xnl,mcore,Cax,Ma,ta,clgr,alpha,SRS,K,nfrec,nt,theta,fband,f,freq,prms)
    integer, intent(in) :: cnt ! loop counter - useful for controlling debugging printouts 
    logical, intent(in) :: testCalcTurbine ! if true => load hard coded data to allow comparison with NoPred. False => run on usual input
    real(kind=rp), intent(in) :: N_rotors,n_stages,Vref,cref,Vtr,Texit,xnl,mcore, Cax, Ma, ta, clgr, alpha, SRS, K
    integer, intent(in) :: nt ! number of theta 
    integer, intent(in) :: nfrec ! number of frequencies
    real(kind=rp), intent(in), dimension(*) :: theta ! directivity
    real(kind=rp), intent(in), dimension(nfreq) :: fband
    real(kind=rp), intent(in), dimension(nfreq+1) :: f
    real(kind=rp), intent(in), dimension(*) :: freq
    real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms
    !
    real(kind=rp) :: t_static, c_exit, fac1, fac3, Prms_temp, Prms_46m_t, Prms_b, Prms_t, Prms_tot, SPL_tot, SPL_tot_46m
    !
    real(kind=rp), dimension(nthet) :: SPL_peak_b ! peak SPL broadband 
    real(kind=rp), dimension(nthet) :: SPL_peak_t ! peak SPL tone
    real(kind=rp), dimension(nthet) :: BPF, F1, F1der, F3, F3der, derivative, BPF_max
    integer, parameter :: noData = 12
    real(kind=rp), dimension(noData) :: broadband_corr = (/-37.0_rp,-29.0_rp,-21.0_rp,-13.0_rp,-7.4_rp,-4.0_rp,-1.2_rp,0.0_rp,-1.2_rp,-9.15_rp,-19.0_rp,-29.0_rp/)
    real(kind=rp), dimension(noData) :: angle = (/0.0_rp,20.0_rp,40.0_rp,60.0_rp,80.0_rp,90.0_rp,100.0_rp,110.0_rp,120.0_rp,140.0_rp,160.0_rp,180.0_rp/) ! empirical F1 curve
    real(kind=rp), dimension(noData) :: tone_corr=(/-47.0_rp,-37.0_rp,-27.0_rp,-18.2_rp,-10.0_rp,-6.0_rp,-2.5_rp,0.0_rp,-2.5_rp,-14.8_rp,-26.0_rp,-37.0_rp/)
    
    integer, parameter :: lwk = 2*noData
    real(kind=rp),  dimension(lwk) :: wk
    real(kind=rp), dimension(nthet,nfreq) :: SPL_46m_b, SPL_46m_t
    real(kind=rp), dimension(nfreq,nthet) :: SPL_1m
    real(kind=rp), dimension(:,:), allocatable :: BPF_harmonics
    real(kind=rp), dimension(nfreq) :: atm_abs
    !
    logical :: spline
    integer :: iptr, ierr, n_harmonics
    integer :: i, j
    
    t_static = Texit - Cax**2/(two*Cp_gas)
    c_exit = sqrt(gamma_gas*get_R(zero)*t_static) ! incorrect value for gas constant. Gamma_gas is only approximative.
    
    do i = 1, nthet
       BPF(i) = nint(N_rotors*xnl/(60.0_rp*(one-Ma*cos(abs(alpha+clgr - theta(i))*(pi_num/180.0_rp)))))
    end do
    
!---- broadband noise
    spline = .true.     ! spline shape preserving interpolation - ala MATLAB
    call pchez ( noData, angle, broadband_corr, derivative, spline, wk, lwk, ierr )
    call pchev ( noData, angle, broadband_corr, derivative, nthet, theta, F1, F1der, ierr )
   
    do i = 1, nthet
      fac1 = (Vtr*cref)/(Vref*c_exit)
      fac3 = (one-Ma*cos(abs(alpha+clgr - theta(i))*(pi_num/180.0_rp))) 
      SPL_peak_b(i)=10.0_rp * log10( fac1**3 * (mcore/mref)* fac3**(-4) ) + F1(i) - 10.0_rp
    end do
    
    do i = 1, nthet
        do j = 1, nfreq
            if( fband(j)/BPF(i) .le. one ) then
              SPL_46m_b(i,j)=SPL_peak_b(i)+10.0_rp*log10(fband(j)/BPF(i))
            else
              SPL_46m_b(i,j)=SPL_peak_b(i)-20.0_rp*log10(fband(j)/BPF(i))
            end if
        end do 
    end do 
    
!---- tone noise
    spline = .true.     ! spline shape preserving interpolation - ala MATLAB
    call pchez ( noData, angle, tone_corr, derivative, spline, wk, lwk, ierr )
    call pchev ( noData, angle, tone_corr, derivative, nthet, theta, F3, F3der, ierr )
   
    do i = 1, nthet
        fac1 = (Vtr/Vref)**0.6_rp
        fac3 = (one-Ma*cos(abs(alpha+clgr - theta(i))*(pi_num/180.0_rp))) 
        SPL_peak_t(i) = 10.0_rp * log10( fac1*(cref/c_exit)**3 * (mcore/mref)*(one/srs)*fac3**(-4) ) + F3(i) + 56.0_rp + K
        BPF_max(i)=fband(nfreq)/BPF(i)
    end do
    n_harmonics = floor(maxval(BPF_max))
    allocate(BPF_harmonics(nthet,n_harmonics)) ; BPF_harmonics = zero
    !
    do i = 1, nthet
        do j = 1, floor(BPF_max(i))
          BPF_harmonics(i,j)=BPF(i)*j
        end do
    end do
    
    SPL_46m_t = zero
    do i = 1, nthet
        do j = 1, floor(BPF_max(i))
            iptr = findFirstElementGtThanBnd(j*BPF(i),nfreq,f)-1
            if( SPL_46m_t(i,iptr).eq.0 ) then
               SPL_46m_t(i,iptr)=SPL_peak_t(i)-(j-1)*10.0_rp
            else
               Prms_temp=P0*10.0_rp**(SPL_46m_t(i,iptr)/20.0_rp)
               Prms_46m_t=P0*10.0_rp**((SPL_peak_t(i)-(j-1)*10.0_rp)/20.0_rp)
               Prms_46m_t=sqrt(Prms_temp**2+Prms_46m_t**2)
               SPL_46m_t(i,iptr)=20*log10(Prms_46m_t/P0)
            end if
            if( SPL_46m_t(i,iptr).lt.zero ) SPL_46m_t(i,iptr) = zero
        end do
    end do   
!    call save_matrix(nthet,nfreq,SPL_46m_t,'.m') ; stop
    !
    atm_abs = get_atm_abs(nfreq,ta,RH,fband)
    do i = 1, nthet
       do j = 1, nfreq
           Prms_b = p0*10.0_rp**(SPL_46m_b(i,j)/20.0_rp)
           Prms_t = p0*10.0_rp**(SPL_46m_t(i,j)/20.0_rp)
           Prms_tot = sqrt(Prms_b**2+Prms_t**2)
           SPL_tot = 20.0_rp*log10(Prms_tot/P0)
           SPL_tot_46m = SPL_tot + 10.0_rp*log10(n_stages)
           !
           SPL_1m(j,i) = SPL_tot_46m + 33.2_rp + atm_abs(j)*0.457_rp
       end do
    end do

    do i = 1, nfreq
        do j = 1, nthet
            prms(i,j) = p0*10.0_rp**(SPL_1m(i,j)/20.0_rp)
        end do
    end do

    deallocate(BPF_harmonics)
    
  contains
    subroutine setTestDataTurbine(xpos) ! load hard coded noPred test data (from file)
      real(kind=rp), intent(in) :: xpos 

      
    end subroutine setTestDataTurbine
  end subroutine calcTurbine
  !
  !
  !
  subroutine calcFanandCompressor(cnt,component,operatingPoint,MtipD,N_rotors,N_stators,rss,Mtip,Mu,dt,xnl,g1,comp,nfrec,nt,theta,fband,f,freq,& 
        & prms_inlet_tone, prms_discharge_tone, prms_inlet_broadband, prms_discharge_broadband, prms_inlet_combination, &  
        & prms_inlet,prms_discharge)
    character(len=3), intent(in) :: component
    integer, intent(in) :: cnt ! loop counter - useful for controlling debugging printouts 
    character(len=*), intent(in) :: operatingPoint
    character(len=*), intent(in) :: comp
    real(kind=rp), intent(in) :: MtipD,N_rotors,N_stators,rss
    real(kind=rp), intent(in) :: Mtip,Mu,dt,xnl,g1
    integer, intent(in) :: nt ! number of theta 
    integer, intent(in) :: nfrec ! number of frequencies
    real(kind=rp), intent(in), dimension(*) :: theta ! directivity
    real(kind=rp), intent(in), dimension(nfreq) :: fband
    real(kind=rp), intent(in), dimension(nfreq+1) :: f
    real(kind=rp), intent(in), dimension(*) :: freq
    real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms_inlet, prms_discharge
    real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms_inlet_tone, prms_discharge_tone, prms_inlet_combination, prms_inlet_broadband, prms_discharge_broadband
    !
    real(kind=rp) :: omega ! angular speed
    integer :: fb ! blade passing frequency
    integer :: acc ! blade passing frequency
    real(kind=rp) :: delta_cutoff, mt
    integer :: i, j
 
!    if( testCalcFan ) call setTestDataFan(ctype,xpos) 
    
    omega = xnl*two*pi_num
    fb    = nint(xnl*N_rotors)
    
    acc=nint(f(nfreq+1)/fb)   
!    if( .not.testCalcFan ) call set_mach_numbers(Mtip,Mu) 

    delta_cutoff = abs(Mu/(one-N_stators/N_rotors)) ! cut-off factor 
! get broadband noise
    call get_fanCompr_broadband(MtipD,Mtip,real(fb,rp),dt,g1,rss,fband,prms_inlet_broadband,prms_discharge_broadband) ! Prms_inlet(1,1)=0.0032; Prms_inlet(20,20)=3.3804
        
! get tone noise
    call get_fanCompr_tone(fan_distortion,acc,MtipD,Mtip,Mu,dt,g1,rss,delta_cutoff,prms_inlet_tone,prms_discharge_tone)    
    
! get combination noise
    if( Mtip .gt. one ) then 
        call get_fanCompr_comb_tone(MtipD,Mtip,dt,g1,fband,f,prms_inlet_combination)
    else
        prms_inlet_combination = zero
    end if
!    call save_matrix(nfrec,nthet,prms_inlet_broadband,'.m') ; stop
    
    if( component(1:3).ne.'Fan' ) prms_discharge_tone = zero ! discharge zero for Lpc/Ipc 
    if( component(1:3).ne.'Fan' ) prms_discharge_broadband = zero ! discharge zero for Lpc/Ipc 
    do i = 1, nfrec
       do j = 1, nt
! create inlet combination (broadband + combination + tone)
          prms_inlet(i,j) = sqrt( prms_inlet_broadband(i,j)**2 + prms_inlet_tone(i,j)**2 + prms_inlet_combination(i,j)**2 )  
! create outlet combintion (tone + broadband)
          prms_discharge(i,j) = sqrt(prms_discharge_tone(i,j)**2 + prms_discharge_broadband(i,j)**2)
!          prms_discharge(i,j) = sqrt(prms_discharge_tone(i,j)**2 )
!          prms_discharge(i,j) = sqrt(prms_discharge_broadband(i,j)**2)
          !
          if( isNaN(prms_discharge(i,j)) .or. isNaN(prms_inlet(i,j)) ) then
              write(*,*) 'something wrong'
          end if
          !
       end do
    end do
    
    

    
!    if( cnt .eq. 2 ) then
!      call save_matrix(nfrec,nthet,prms_discharge_tone,'.m') 
!      stop
!    end if    
    
  contains
    subroutine get_fanCompr_comb_tone(MtrD,Mtr,dt,g1,fband,f,Prms_total)
      real(kind=rp), intent(in) :: MtrD,Mtr,dt,g1
      real(kind=rp), intent(in), dimension(nfreq) :: fband
      real(kind=rp), intent(in), dimension(nfreq+1) :: f
      real(kind=rp), dimension(nfreq,nthet), intent(out) :: Prms_total
      !      
      real(kind=rp), dimension(n_dir_data) :: f2_data = (/-9.50_rp,-9.00_rp,-8.50_rp,-7.75_rp,-7.00_rp,-6.00_rp,-5.00_rp,-3.50_rp,-2.00_rp,-1.00_rp, & 
         & 0.00_rp, 0.00_rp, 0.00_rp,-1.75_rp,-3.50_rp,-5.50_rp,-7.50_rp,-8.25_rp,-9.00_rp,-9.25_rp,-9.50_rp,-9.75_rp,-10.00_rp,-10.25_rp,-10.50_rp,-10.75_rp,& 
         & -11.00_rp,-11.25_rp,-11.50_rp,-11.75_rp,-12.00_rp,-12.25_rp,-12.50_rp,-12.75_rp,-13.00_rp,-13.25_rp,-13.50_rp /)
      real(kind=rp), dimension(3,nthet) :: peak
      real(kind=rp), dimension(nthet) :: Lc_combination
      real(kind=rp), dimension(nfreq,nthet,3) :: SPL_combination, Prms_combination
      real(kind=rp), dimension(3,nfreq) :: F3
      real(kind=rp) :: C
      real(kind=rp), dimension(3) :: F1, peakf
      integer :: i, k, idump

! Creating F1 according to fig 15(a) F1=[f/fb=1/2;f/fb=1/4;f/fb=1/8]
      call linear_interp(get_gestpan_1d_table(3,(/1.0_rp,1.14_rp,2.0_rp/),(/30.0_rp,72.5_rp, 4.4_rp/)),1,Mtr,zero,F1(1),idump)      
      call linear_interp(get_gestpan_1d_table(3,(/1.0_rp,1.25_rp,2.0_rp/),(/30.0_rp,68.6_rp,10.5_rp/)),1,Mtr,zero,F1(2),idump)      
      call linear_interp(get_gestpan_1d_table(3,(/1.0_rp,1.61_rp,2.0_rp/),(/36.0_rp,60.6_rp,56.5_rp/)),1,Mtr,zero,F1(3),idump)            
      
      if( no_fan_stages.eq.1 .and. fan_IGV ) then
         C = 5.0_rp
      else
         C = zero
      end if
      
      do i = 1, nt
         Lc_combination(i) = 20.0_rp*log10(dT/dTref) + 10.0_rp*log10(g1/mref) + f2_data(i) + C ! eq 8     
      end do     
      
! set peak
      do i = 1, 3
         do j = 1, nt
            peak(i,j) = Lc_combination(j) + F1(i)
         end do
      end do
      
!
      peakf(1) = findFirstElementGtThanBnd(0.500_rp*fb,nfrec,f) 
      peakf(2) = findFirstElementGtThanBnd(0.250_rp*fb,nfrec,f) 
      peakf(3) = findFirstElementGtThanBnd(0.125_rp*fb,nfrec,f)
      

! Figure 14
      do i = 1, nfrec ! first row (half) - top curve 
          if( i.lt.peakf(1) ) then
             F3(1,i) = 30.0_rp*log10(two*fband(i)/fb)  
          else
             F3(1,i) = -30.0_rp*log10(two*fband(i)/fb)  
          end if
      end do
!      call save_matrix(3,nt,peak(:,:),'.m')
      
      
      do i = 1, nfrec ! second row (quarter) - mid curve 
          if( i.lt.peakf(2) ) then
             F3(2,i) = 50.0_rp*log10(4.0_rp*fband(i)/fb)  
          else
             F3(2,i) = -50.0_rp*log10(4.0_rp*fband(i)/fb)  
          end if
      end do
      do i = 1, nfrec ! third row (eighth) - bottom curve
          if( i.lt.peakf(3) ) then
             F3(3,i) = 50.0_rp*log10(8.0_rp*fband(i)/fb)  
          else
             F3(3,i) = -30.0_rp*log10(8.0_rp*fband(i)/fb) 
          end if
      end do
      
      do i = 1, nfrec
          do j = 1, nt
            SPL_combination(i,j,1) = peak(1,j) + F3(1,i)
            SPL_combination(i,j,2) = peak(2,j) + F3(2,i) 
            SPL_combination(i,j,3) = peak(3,j) + F3(3,i)
          end do
      end do
!      call save_matrix(nfrec,nt,SPL_combination(:,:,1),'.m')

      do i = 1, nfrec
          do j = 1, nt
              do k = 1, 3
                 Prms_combination(i,j,k) = P0*10**(SPL_combination(i,j,k)/20.0_rp)
              end do
              Prms_total(i,j) = sqrt( Prms_combination(i,j,1)**2 + Prms_combination(i,j,2)**2 + Prms_combination(i,j,3)**2 )    
          end do
      end do
!      call save_matrix(nfrec,nt,Prms_total(:,:),'.m')
    
    end subroutine get_fanCompr_comb_tone
    !
    subroutine get_fanCompr_tone(distortion,acc,MtrD,Mtr,Mu,dt,g1,rss,delta_cutoff,prms_inlet,prms_discharge)
      logical, intent(in) :: distortion
      integer, intent(in) :: acc
      real(kind=rp), intent(in) :: MtrD,Mtr,Mu,dt,g1,rss,delta_cutoff
      real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms_inlet, prms_discharge
      !
      integer, dimension(:), allocatable :: k ! size depends on acc which is computed based on user input
      real(kind=rp), dimension(:), allocatable :: F4a, F4b, F5 ! size depends on acc which is computed based on user input
      !
      real(kind=rp), dimension(nthet) :: Lc_inlet, Lc_discharge
      real(kind=rp), dimension(nthet) :: BPF1, BPF2
      real(kind=rp), dimension(nfreq,nthet) :: SPL_inlet, SPL_discharge
      real(kind=rp), dimension(n_dir_data) :: f3a_data = &
        (/-3.00_rp,-2.25_rp,-1.50_rp,-0.75_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,-0.60_rp,-1.20_rp,-2.35_rp,-3.50_rp,-5.15_rp, & 
          & -6.80_rp, -8.65_rp,-10.50_rp,-12.50_rp,-14.50_rp,-16.75_rp,-19.00_rp, -21.25_rp,-23.50_rp,-25.75_rp,-28.00_rp,-30.25_rp,-32.50_rp,-34.75_rp, & 
          & -37.00_rp,-39.25_rp,-41.50_rp,-43.75_rp,-46.00_rp,-48.25_rp,-50.50_rp,-52.75_rp,-55.00_rp/) ! figure 13(a) FOR INLET (extrapolateded from 100 degrees to 180)
      real(kind=rp), dimension(n_dir_data) :: f3b_data = (/-39.0_rp,-37.0_rp,-35.0_rp,-33.0_rp,-31.0_rp,-29.0_rp,-27.0_rp, & 
        & -25.0_rp,-23.0_rp,-21.0_rp,-19.0_rp,-17.0_rp,-15.0_rp,-13.0_rp,-11.0_rp, -9.5_rp, -8.0_rp, -6.5_rp, -5.0_rp, -4.0_rp, -3.0_rp, & 
        & -2.0_rp, -1.0_rp, -0.5_rp,  0.0_rp,  0.0_rp,  0.0_rp, -1.0_rp,-2.0_rp,-3.75_rp, -5.5_rp,-7.25_rp, -9.0_rp,-11.0_rp,-13.0_rp,-15.5_rp,-18.0_rp/)
      integer, parameter :: np = 16
      real(kind=rp), dimension(2*np+2) :: gtab
      real(kind=rp), dimension(np) :: angles = (/10.0_rp,20.0_rp,30.0_rp,40.0_rp,50.0_rp,60.0_rp,70.0_rp,80.0_rp,90.0_rp,100.0_rp,110.0_rp,120.0_rp,130.0_rp,140.0_rp,150.0_rp,160.0_rp/)
      
! GE "Flight Cleanup" -- TCS Suppression - Table 4.3. Idea is to make sure that
! the installation effects in take-off and approach are removed.
      real(kind=rp), dimension(np) ::  approach_BPF = (/5.6_rp, 5.8_rp, 4.7_rp, 4.6_rp, 4.9_rp, 5.1_rp, 2.9_rp, 3.2_rp, 1.6_rp, 1.6_rp, 1.8_rp, 2.1_rp, 2.4_rp, 2.2_rp, 2.0_rp, 2.8_rp/)
      real(kind=rp), dimension(np) :: approach_2BPF = (/5.4_rp, 4.3_rp, 3.4_rp, 4.1_rp, 2.0_rp, 2.9_rp, 1.6_rp, 1.3_rp, 1.5_rp, 1.1_rp, 1.4_rp, 1.5_rp, 1.0_rp, 1.8_rp, 1.6_rp, 1.6_rp/)
      real(kind=rp), dimension(np) ::   takeoff_BPF = (/4.8_rp, 5.5_rp, 5.5_rp, 5.3_rp, 5.3_rp, 5.1_rp, 4.4_rp, 3.9_rp, 2.6_rp, 2.3_rp, 1.8_rp, 2.1_rp, 1.7_rp, 1.7_rp, 2.6_rp, 3.5_rp/)
      real(kind=rp), dimension(np) ::  takeoff_2BPF = (/5.8_rp, 3.8_rp, 5.3_rp, 6.4_rp, 3.5_rp, 3.0_rp, 2.1_rp, 2.1_rp, 1.1_rp, 1.4_rp, 0.9_rp, 0.7_rp, 0.7_rp, 0.4_rp, 0.6_rp, 0.8_rp/)
      real(kind=rp) :: F1a, F1b, F2a, F2b, Prms_temp, Prms_tot
      real(kind=rp) :: temp1, temp2, temp3
      integer :: i, j, spos, nk, pos, idump
      integer :: na, slen ! number of angles
      !
      !Double layer Liner Data is from NASA CR-2002-211672 Report
      !Advanced Nacelle Acoustic lining Concepts Development
      !Linearly extrapolated from 0 - 90 deg and form 160-180
      !The data was obtained only for brodband noise (we need to find data for tonal noise)
      real(kind=rp), dimension(n_dir_data) :: liner_data_app = (/&
         &  zero, zero, zero, zero, zero, zero, -0.31149048974951_rp, -0.75394155342520_rp, -1.19639261710089_rp, -1.63884368077659_rp, -2.08129474445228_rp, -2.52374580812797_rp, &
         & -2.96619687180366_rp, -3.40864793547935_rp, -3.85109899915504_rp, -4.29355006283073_rp, -4.73600112650642_rp, -5.17845219018211_rp, -5.62090325385780_rp, &
         & -6.06335431753349_rp, -6.50580538120918_rp, -6.55856451491287_rp, -6.60963444765545_rp, -6.01819645866627_rp, -5.40364544803953_rp, -4.88024833985105_rp, &
         & -4.36166961630789_rp, -3.91871616323577_rp, -3.47659239938625_rp, -3.08872017651543_rp, -2.70434534445140_rp, -2.53311635918184_rp, -2.36917034195449_rp, &
         & -2.20522432472714_rp, -2.04127830749979_rp, -1.87733229027244_rp, -1.71338627304509_rp/)
      real(kind=rp), dimension(n_dir_data) :: liner_data_cut = (/&
         & -1.13306098626510_rp, -1.45911467702874_rp, -1.78516836779238_rp, -2.11122205855601_rp, -2.43727574931965_rp, -2.76332944008329_rp, -3.08938313084693_rp, &
         & -3.41543682161056_rp, -3.74149051237420_rp, -4.06754420313784_rp, -4.39359789390148_rp, -4.71965158466511_rp, -5.04570527542875_rp, -5.37175896619239_rp, &
         & -5.69781265695602_rp, -6.02386634771966_rp, -6.34992003848330_rp, -6.67597372924694_rp, -7.00202742001057_rp, -7.32808111077421_rp, -7.65413480153785_rp, &
         & -7.30873649087542_rp, -6.95315039191579_rp, -7.17795829812488_rp, -7.41074393672352_rp, -7.01984807951155_rp, -6.61979710153490_rp, -5.66478064174521_rp, &
         & -4.68471569886730_rp, -3.93953462123672_rp, -3.20789049050141_rp, -2.75748299134357_rp, -2.31773313596766_rp, -1.87798328059175_rp, -1.43823342521584_rp, &
         & -0.99848356983992_rp, -0.55873371446401_rp/)
      real(kind=rp), dimension(n_dir_data) :: liner_data_sid = (/&
         & -5.54049164663243_rp, -5.69865908856118_rp, -5.85682653048993_rp, -6.01499397241868_rp, -6.17316141434743_rp, -6.33132885627618_rp, -6.48949629820493_rp, &
         & -6.64766374013368_rp, -6.80583118206242_rp, -6.96399862399117_rp, -7.12216606591992_rp, -7.28033350784867_rp, -7.43850094977742_rp, -7.59666839170617_rp, &
         & -7.75483583363492_rp, -7.91300327556367_rp, -8.07117071749241_rp, -8.22933815942116_rp, -8.38750560134991_rp, -8.54567304327866_rp, -8.68994505933804_rp, &
         & -8.32871946157918_rp, -7.96749386382032_rp, -7.48047644940415_rp, -6.99409867992745_rp, -6.86526861357872_rp, -6.73643854723000_rp, -6.60760848088128_rp, &
         & -6.47877841453256_rp, -6.34994834818384_rp, -6.22111828183511_rp, -6.09228821548639_rp, -5.96345814913767_rp, -5.83462808278895_rp, -5.70579801644022_rp, &
         & -5.57696795009150_rp, -5.44813788374278_rp/)
      real(kind=rp) :: liner_data
      logical :: liners_tone = .true.
      
      
      F1a = get_Fig10a(MtrD,Mtr) ! Creating F1a according to figure 10(a) FOR INLET
      F1b = get_Fig10b(MtrD,Mtr) ! Creating F1b according to figure 10(b) FOR DISCHARGE
    
! Creating F2a according to figure 12 FOR INLET
      F2a = zero;

! Creating F2b according to figure 12 FOR DISCHARGE
      F2b=-10.0_rp*log10(rss/300.0_rp)
      
      call init_ivector(1,acc+1,1,nk,k)
      
! re-interpolate GE-data onto theta angles (GE "Flight Cleanup" -- TCS Suppression - Table 4.3. Idea is to make sure that)      
      slen = len(operatingPoint) 
      if( operatingPoint(1:slen).eq.'Take-off' .or. operatingPoint(1:slen).eq.'Cutback' .or. operatingPoint(1:slen).eq.'Sideline' ) then
         gtab(1:2*np+2) = get_gestpan_1d_table(np,angles, takeoff_BPF)
         do i = 1, nt
            call linear_interp(gtab,1,theta(i),zero,BPF1(i),idump)
         end do 
         gtab(1:2*np+2) = get_gestpan_1d_table(np,angles, takeoff_2BPF)
         do i = 1, nt
            call linear_interp(gtab,1,theta(i),zero,BPF2(i),idump)
         end do           
      elseif( operatingPoint(1:slen).eq.'Approach' ) then
         gtab(1:2*np+2) = get_gestpan_1d_table(np,angles, approach_BPF)
         do i = 1, nt
            call linear_interp(gtab,1,theta(i),zero,BPF1(i),idump)
         end do 
         gtab(1:2*np+2) = get_gestpan_1d_table(np,angles, approach_2BPF)
         do i = 1, nt
            call linear_interp(gtab,1,theta(i),zero,BPF2(i),idump)
         end do           
      else
          BPF1 = zero ; BPF2 = zero
      end if
! --------------- BPF1 & BPF2 interpolation complete
      
! set inlet data
      allocate(F5(1:nk)) 
      F5(1:nk) = get_Fig9(no_fan_stages,nk,k)
    
      do i = 1, nt ! nt is in scope becaouse this is an internal subroutine
         Lc_inlet(i)= 20.0_rp*log10(dT/dTref) + 10.0_rp*log10(g1/mref) + F1a + F2a + F3a_data(i) ! eq 6
      end do


      allocate(F4a(1:nk))  
      F4a(1:nk) = get_Fig8_inlet(Mu,j,nk,k,delta_cutoff)
      
      SPL_inlet = zero ; prms_inlet = zero 
      do i = 1, nt
          do j = 1, acc
              !
              pos = findFirstElementGtThanBnd(real(fb*j,rp),nfrec,f) -1
              if( pos.eq.-1 ) cycle
              !
              if( SPL_inlet(pos,i).eq.zero ) then
                if( j.eq.1 ) then ! first tone
                  SPL_inlet(pos,i) = Lc_inlet(i) + 10.0_rp*log10( 10.0_rp**(0.10_rp*F4a(1)) + 10.0_rp**(-0.10_rp*BPF1(i)) ) 
                elseif( j.eq.2 ) then ! second tone
                  SPL_inlet(pos,i) = Lc_inlet(i) + 10.0_rp*log10( 10.0_rp**(0.10_rp*F4a(2)) + 10.0_rp**(-0.10_rp*BPF2(i)) )
                else
                  SPL_inlet(pos,i) = Lc_inlet(i) + 10.0_rp*log10( 10.0_rp**(0.10_rp*F4a(j)) + 10.0_rp**(0.10_rp*F5(j)) )
                end if
              else ! coinciding tones
                  Prms_temp = P0*10**(SPL_inlet(pos,i)/20.0_rp)
                  if( j.eq.1 ) then ! first tone
                    Prms_inlet=       P0*10.0_rp**((Lc_inlet(i)+10.0_rp*log10( 10.0_rp**(0.10_rp*F4a(1)) + 10.0_rp**(-0.10_rp**BPF1(i))))/20.0_rp)
                  elseif( j.eq.2 ) then ! second tone
                    Prms_inlet=       P0*10.0_rp**((Lc_inlet(i)+10.0_rp*log10( 10.0_rp**(0.10_rp*F4a(2)) + 10.0_rp**(-0.10_rp**BPF2(i))))/20.0_rp)
                  else
                    Prms_inlet(pos,i)=P0*10.0_rp**((Lc_inlet(i)+10.0_rp*log10( 10.0_rp**(0.10_rp*F4a(j)) + 10.0_rp**(0.10_rp*F5(j))))/20.0_rp)
                  end if
                  Prms_tot = sqrt(Prms_temp**2 + Prms_inlet(pos,i)**2)
                  SPL_inlet(pos,i)=20.0_rp*log10(Prms_tot/P0)
              end if
          end do
      end do
      
!     call save_matrix(nfreq,nthet,SPL_inlet,'.m') ; stop

! set discharge data
      allocate(F4b(1:nk)) ; F4b = zero 
      F4b(1:nk) = get_Fig8_discharge(fan_IGV,j,nk,k,delta_cutoff)
      
      Lc_discharge = zero
      do i=1, nt
          if (liners_tone) then 
              if (comp(1:3) .eq. 'Fan') then
                  if (operatingPoint(1:slen).eq.'Sideline') then
                      liner_data = liner_data_sid(i) 
                  elseif (operatingPoint(1:slen).eq.'Cutback') then
                      liner_data = liner_data_cut(i)
                  elseif (operatingPoint(1:slen).eq.'Approach') then
                      liner_data = liner_data_app(i)
                  endif 
              else
                  liner_data = zero
              endif 
          endif 
        Lc_discharge(i) = 20.0_rp*log10(dT/dTref) + 10.0_rp*log10(g1/mref) + F1b + F2b + F3b_data(i) + get_C() + liner_data ! eq 6
      end do      
    
      SPL_discharge = zero
      do i = 1, nt
          do j = 1, acc
              pos = findFirstElementGtThanBnd(real(fb*j,rp),nfrec,f) -1
              if( pos.eq.-1 ) cycle
              if( SPL_discharge(pos,i).eq.zero ) then
                 SPL_discharge(pos,i) = Lc_discharge(i) + F4b(j)                 
              elseif( SPL_discharge(pos,i).ne.zero ) then
                 Prms_temp      = P0*10**(SPL_discharge(pos,i)/20.0_rp)
                 Prms_discharge(pos,i) = P0*10**((Lc_discharge(i)+F4b(j))/20.0_rp)
                 Prms_tot = sqrt(Prms_temp**2 + Prms_discharge(pos,i)**2)
                 SPL_discharge(pos,i)=20.0_rp*log10(Prms_tot/P0)
              end if
          end do
      end do
      
      do i=1, nt
        do j = 1, nfrec
           prms_inlet(j,i) = p0*10**(SPL_inlet(j,i)/20.0_rp)
           prms_discharge(j,i) = p0*10**(SPL_discharge(j,i)/20.0_rp)
        end do 
      end do      
      
      deallocate(k) ; deallocate(F5) ; deallocate(F4a) ; deallocate(F4b)
               
    end subroutine get_fanCompr_tone
    !
    subroutine get_fanCompr_broadband(MtrD,Mtr,fb,dt,g1,rss,fband,prms_inlet,prms_discharge)
      real(kind=rp), intent(in) :: MtrD, Mtr, fb, dt, g1, rss
      real(kind=rp), intent(in), dimension(nfreq) :: fband
      real(kind=rp), dimension(nfreq,nthet), intent(out) :: prms_inlet, prms_discharge
      !
      real(kind=rp) :: F1a, F1b
      real(kind=rp) :: F2a, F2b
      real(kind=rp), dimension(:), allocatable :: F3b, F4b
      real(kind=rp), dimension(nfreq,nthet) :: Lc_inlet, Lc_discharge, SPL_inlet, SPL_discharge
      real(kind=rp), dimension(n_dir_data) :: f3a_data = (/-2.0_rp,-1.5_rp,-1.0_rp,-0.5_rp,zero,zero,zero,zero,zero,-1.0_rp,-2.0_rp,-3.25_rp,-4.5_rp,  &
         &-6.0_rp,-7.5_rp,-9.25_rp,-11.0_rp,-13.0_rp,-15.0_rp,-17.0_rp,-19.0_rp,-22.0_rp,-25.0_rp,-28.0_rp,-31.0_rp,-34.0_rp,  &
         & -37.0_rp,-40.0_rp,-43.0_rp,-46.0_rp,-49.0_rp,-52.0_rp,-55.0_rp,-58.0_rp,-61.0_rp,-64.0_rp,-67.0_rp/)
      real(kind=rp), dimension(n_dir_data) :: f3b_data = (/&
         & -41.600_rp,-39.450_rp,-37.300_rp,-35.150_rp,-33.000_rp,-30.850_rp,-28.700_rp,-26.550_rp,-24.400_rp,-22.250_rp,-20.100_rp,-17.950_rp,-15.8000_rp, & 
         & -13.650_rp,-11.500_rp,-9.750_rp,-8.000_rp,-6.500_rp,-5.000_rp,-3.850_rp,-2.700_rp,-1.950_rp,-1.200_rp,-0.750_rp,-0.300_rp,-0.1500_rp, & 
         &  0.000_rp,-1.000_rp,-2.000_rp,-4.000_rp,-6.000_rp,-8.000_rp,-10.000_rp,-12.500_rp,-15.000_rp,-17.500_rp,-20.000_rp/)

      integer :: i, j, slen
      !
      !Double layer Liner Data is from NASA CR-2002-211672 Report
      !Advanced Nacelle Acoustic lining Concepts Development
      !Linearly extrapolated from 0 - 90 deg and form 160-180
      real(kind=rp), dimension(n_dir_data) :: liner_data_app = (/&
         &  zero, zero, zero, zero, zero, zero, -0.31149048974951_rp, -0.75394155342520_rp, -1.19639261710089_rp, -1.63884368077659_rp, -2.08129474445228_rp, -2.52374580812797_rp, &
         & -2.96619687180366_rp, -3.40864793547935_rp, -3.85109899915504_rp, -4.29355006283073_rp, -4.73600112650642_rp, -5.17845219018211_rp, -5.62090325385780_rp, &
         & -6.06335431753349_rp, -6.50580538120918_rp, -6.55856451491287_rp, -6.60963444765545_rp, -6.01819645866627_rp, -5.40364544803953_rp, -4.88024833985105_rp, &
         & -4.36166961630789_rp, -3.91871616323577_rp, -3.47659239938625_rp, -3.08872017651543_rp, -2.70434534445140_rp, -2.53311635918184_rp, -2.36917034195449_rp, &
         & -2.20522432472714_rp, -2.04127830749979_rp, -1.87733229027244_rp, -1.71338627304509_rp/)
      real(kind=rp), dimension(n_dir_data) :: liner_data_cut = (/&
         & -1.13306098626510_rp, -1.45911467702874_rp, -1.78516836779238_rp, -2.11122205855601_rp, -2.43727574931965_rp, -2.76332944008329_rp, -3.08938313084693_rp, &
         & -3.41543682161056_rp, -3.74149051237420_rp, -4.06754420313784_rp, -4.39359789390148_rp, -4.71965158466511_rp, -5.04570527542875_rp, -5.37175896619239_rp, &
         & -5.69781265695602_rp, -6.02386634771966_rp, -6.34992003848330_rp, -6.67597372924694_rp, -7.00202742001057_rp, -7.32808111077421_rp, -7.65413480153785_rp, &
         & -7.30873649087542_rp, -6.95315039191579_rp, -7.17795829812488_rp, -7.41074393672352_rp, -7.01984807951155_rp, -6.61979710153490_rp, -5.66478064174521_rp, &
         & -4.68471569886730_rp, -3.93953462123672_rp, -3.20789049050141_rp, -2.75748299134357_rp, -2.31773313596766_rp, -1.87798328059175_rp, -1.43823342521584_rp, &
         & -0.99848356983992_rp, -0.55873371446401_rp/)
      real(kind=rp), dimension(n_dir_data) :: liner_data_sid = (/&
         & -5.54049164663243_rp, -5.69865908856118_rp, -5.85682653048993_rp, -6.01499397241868_rp, -6.17316141434743_rp, -6.33132885627618_rp, -6.48949629820493_rp, &
         & -6.64766374013368_rp, -6.80583118206242_rp, -6.96399862399117_rp, -7.12216606591992_rp, -7.28033350784867_rp, -7.43850094977742_rp, -7.59666839170617_rp, &
         & -7.75483583363492_rp, -7.91300327556367_rp, -8.07117071749241_rp, -8.22933815942116_rp, -8.38750560134991_rp, -8.54567304327866_rp, -8.68994505933804_rp, &
         & -8.32871946157918_rp, -7.96749386382032_rp, -7.48047644940415_rp, -6.99409867992745_rp, -6.86526861357872_rp, -6.73643854723000_rp, -6.60760848088128_rp, &
         & -6.47877841453256_rp, -6.34994834818384_rp, -6.22111828183511_rp, -6.09228821548639_rp, -5.96345814913767_rp, -5.83462808278895_rp, -5.70579801644022_rp, &
         & -5.57696795009150_rp, -5.44813788374278_rp/)
      logical :: liners
      real(kind=rp) :: liner_data
      liners = .true.
      !
      !
      F1a = get_Fig4a(MtrD,Mtr) 
      F1b = get_Fig4b(MtrD,Mtr) 
       
      ! Creating F2a according to Fig 6(a) FOR INLET
      F2a=zero
      ! Creating F2b according to Fig 6(b) FOR DISCHARGE
      F2b=-5.0_rp*log10(rss/300.0_rp)
    
      do i=1, nt
        do j = 1, nfrec
           Lc_inlet(j,i) = 20.0_rp*log10(dT/dTref) + 10.0_rp*log10(g1/mref) + F1a + F2a + f3a_data(i)! eq 4
           SPL_inlet(j,i) = Lc_inlet(j,i) + get_Fig3a(fband(j)/fb) ! eq 5
        end do 
      end do
 
    slen = len(operatingPoint)       
      do i=1, nt
        do j = 1, nfrec
            if (liners) then
                if (comp(1:3) .eq. 'Fan') then
                    if (operatingPoint(1:slen).eq.'Sideline') then
                        liner_data = liner_data_sid(i) 
                    elseif (operatingPoint(1:slen).eq.'Cutback') then
                        liner_data = liner_data_cut(i)
                    elseif (operatingPoint(1:slen).eq.'Approach') then
                        liner_data = liner_data_app(i)
                    endif
                else
                    liner_data = zero
                endif
            endif 
           Lc_discharge(j,i) = 20.0_rp*log10(dT/dTref) + 10.0_rp*log10(g1/mref) + F1b + F2b + f3b_data(i) + liner_data ! eq 4
           SPL_discharge(j,i) = Lc_discharge(j,i) + get_Fig3a(fband(j)/fb) ! eq 5
        end do 
      end do
            
      do i=1, nt
        do j = 1, nfrec
           prms_inlet(j,i) = p0*10**(SPL_inlet(j,i)/20.0_rp)
           prms_discharge(j,i) = p0*10**(SPL_discharge(j,i)/20.0_rp)
        end do 
      end do
      
    end subroutine get_fanCompr_broadband
  end subroutine calcFanandCompressor
  !
  !
  !
  function get_C()
    real(kind=rp) :: get_C
    
    if( fan_IGV .or. no_fan_stages.eq.2 ) then
       get_C = 6.0_rp
    else
       get_C = zero
    end if  
  
  end function get_C
  !
  ! test cases: get_F4b(50.0_rp/2111.0_rp)=  -75.8259 ;   get_F4b(8063.5_rp/2111.0_rp) =  -0.6276 
  !
  function get_Fig3a(f_over_fb)
    real(kind=rp), intent(in) :: f_over_fb
    real(kind=rp) :: get_Fig3a
    !
    real(kind=rp), parameter :: sigmae = 2.2_rp
  
    ! Creating F4a according to equation 2 FOR INLET
    get_Fig3a = 10.0_rp*log10(  one / (exp( (half) * ( log(f_over_fb/2.5_rp) / log(sigmae)) **2 )))   
 
  end function get_Fig3a
  !
  ! Get peak broadband sound pressure levels from Fig 4(a) for inlet according to Ref. 1 (see above)
  !
  function get_Fig4a(M_trd,M_tr)
     real(kind=rp), intent(in) :: M_trd, M_tr
     real(kind=rp) :: get_Fig4a 
     
     if ( M_trD.le.one .and. M_tr.le.0.9_rp ) then
       get_Fig4a = 58.5_rp
     elseif( M_trD.le.one .and. M_tr.gt.0.9_rp ) then
       get_Fig4a = 58.5_rp - 20.0_rp*log10(M_tr/0.9_rp)
     elseif( M_trD.gt.one .and. M_tr.le.0.9_rp ) then
       get_Fig4a = 58.5_rp + 20.0_rp*log10(M_trD)
     elseif( M_trD.gt.one .and. M_tr.gt.0.9_rp ) then
!      get_Fig4a = 58.5_rp + 20.0_rp*log10(M_trD) - 20.0_rp*log10(M_tr/0.9_rp) ! original Heidmann method
       get_Fig4a = 58.5_rp + 20.0_rp*log10(M_trD) - 50.0_rp*log10(M_tr/0.9_rp) ! ANOPP update (1996)
     else
        call report_error('unexpected combination of tip Mach numbers','get_Fig4a','choice_physics')
     endif
  end function get_Fig4a
  !
  !
  !
  function get_Fig8_inlet(Mu,j,nk,k,delta)
    real(kind=rp), intent(in) :: Mu
    integer, intent(in) :: j 
    integer, intent(in) :: nk ! number of elements in k
    integer, intent(in), dimension(nk) ::  k
    real(kind=rp), intent(in) :: delta
    real(kind=rp), dimension(nk) :: get_Fig8_inlet
    !
    integer :: i 
    
    if( Mu.lt.1.15_rp ) then
        if( delta .gt. 1.05_rp ) then
            do i = 1, nk
              get_Fig8_inlet(i) = 6.0_rp - 6.0_rp*k(i) 
            end do
        else
            get_Fig8_inlet(1)  = -8.0_rp  
            do i = 2, nk
              get_Fig8_inlet(i) =  6.0_rp - 6.0_rp*k(i)
            end do
        end if
    else
        if ( delta .gt. 1.05_rp ) then
            do i = 1, nk
              get_Fig8_inlet(i) = 9.0_rp-9.0_rp*k(i)
            end do
        else
           get_Fig8_inlet(1)  = -8.0_rp 
           do i = 1, nk
             get_Fig8_inlet(i) = -9.0_rp - 9.0_rp*k(i)
           end do
        end if
    end if
    
  end function get_Fig8_inlet
  !
  !
  !
  function get_Fig8_discharge(igv,j,nk,k,delta)
    logical, intent(in) :: igv
    integer, intent(in) :: j 
    integer, intent(in) :: nk ! number of elements in k
    integer, intent(in), dimension(nk) ::  k
    real(kind=rp), intent(in) :: delta
    real(kind=rp), dimension(nk) :: get_Fig8_discharge
    !
    integer :: i 
    
    if( .not.IGV ) then
        if( delta .gt. 1.05_rp ) then
            do i = 1, nk
              get_Fig8_discharge(i) = 3.0_rp - 3.0_rp*k(i) 
            end do
        else
            get_Fig8_discharge(1:2)  = -13.0_rp + 5.0_rp*k(1:2)
            get_Fig8_discharge(3:nk) =   3.0_rp - 3.0_rp*k(3:nk)
        end if
    else
        if ( delta .gt. 1.05_rp ) then
           get_Fig8_discharge(1:2)=9.0_rp-9.0_rp*k(1:2)
           get_Fig8_discharge(3:nk)=-3.0_rp-3.0_rp*k(3:nk)
        else
           get_Fig8_discharge(1:2)  = -7.0_rp - k(1:2)
           get_Fig8_discharge(3:nk) = -3.0_rp - 3.0_rp*k(3:nk)
        end if
    end if
    
  end function get_Fig8_discharge
  !
  !
  !
  function get_Fig9(no_fan_stages,nk,k)
    integer, intent(in) :: no_fan_stages
    integer, intent(in) :: nk ! number of elements in k
    integer, intent(in), dimension(nk) ::  k
    real(kind=rp), dimension(nk) :: get_Fig9
    !
    integer :: i 
    
    if( no_fan_stages.eq.1 ) then
       do i = 1, nk
         get_Fig9(i) = 10.0_rp - 10.0_rp*k(i) 
       end do
     elseif( no_fan_stages.eq.2 ) then
       do i = 1, nk
         get_Fig9(i) = zero  
       end do
     else
       call report_error('Unexpected number of stages','get_Fig9','choice_physics')
    end if
    
  end function get_Fig9
  !
  ! Get get_Fig10a according to figure 10(a) FOR INLET
  !
  function get_Fig10a(M_trd,M_tr)
     real(kind=rp), intent(in) :: M_trd, M_tr
     real(kind=rp) :: get_Fig10a
     !
     real(kind=rp) :: F1a1, F1a2
     
     !Creating get_Fig10a according to figure 10(a) FOR INLET
     if( M_tr.le.0.72_rp  ) then
        if( M_trd.le.one ) then 
           get_Fig10a = 60.5_rp
        elseif( M_trd.gt.one ) then  
           get_Fig10a = 60.5_rp + 20._rp*log10(M_trD)
        end if
     elseif( M_tr.gt.0.72 ) then 
       F1a1 = 60.5_rp + 20.0_rp*log10(M_trD) + 50.0_rp*log10(M_tr/0.72_rp)
       F1a2 = 64.5_rp + 80.0_rp*log10(M_trD/M_tr)
       get_Fig10a = min(F1a1,F1a2)
    else
      call report_error('unexpected combination of tip Mach numbers','get_Fig10a','choice_physics')
    endif
    
  end function get_Fig10a
  !
  ! Get get_Fig10b according to figure 10(b) FOR INLET
  !
  function get_Fig10b(M_trd,M_tr)
     real(kind=rp), intent(in) :: M_trd, M_tr
     real(kind=rp) :: get_Fig10b
     !
     real(kind=rp) :: F1b
     !
     !Creating F1b according to figure 10(b) FOR DISCHARGE
     if ( M_tr.le.one ) then
        if( M_trD .le. one ) then
          get_Fig10b = 63.0_rp
        elseif ( M_trD .gt. one ) then
          get_Fig10b = 63.0_rp + 20.0_rp*log10(M_trD);
        end if
     elseif( M_tr.gt.one ) then 
       get_Fig10b = 63.0_rp + 20.0_rp*log10(M_trD) - 20.0_rp*log10(M_tr)
     else
       call report_error('unexpected combination of tip Mach numbers','get_Fig10b','choice_physics')
     endif
    
  end function get_Fig10b
  !
  ! Get peak broadband sound pressure levels from Fig 4(b) for discharge according to Ref. 1 (see above)
  !
  function get_Fig4b(M_trd,M_tr)
     real(kind=rp), intent(in) :: M_trd, M_tr
     real(kind=rp) :: get_Fig4b
     
     if ( M_trD.le.one .and. M_tr.le.one ) then
       get_Fig4b = 60.0_rp
     elseif( M_trD.gt.one .and. M_tr.le.one ) then
       ! get_Fig4b = 60.0_rp + 20.0_rp*log10(M_trD) ! original Heidmann correlation
       get_Fig4b = 63.0_rp + 20.0_rp*log10(M_trD) ! ANOPP revised code
     elseif( M_trD.gt.one .and. M_tr.gt.one ) then
       ! get_Fig4b = 60.0_rp + 20.0_rp*log10(M_trD) - 20.0_rp*log10(M_tr) ! original Heidmann correlation
       get_Fig4b = 63.0_rp + 20.0_rp*log10(M_trD) - 30.0_rp*log10(M_tr) ! ANOPP revised code
     else
        call report_error('unexpected combination of tip Mach numbers','get_Fig4b','choice_physics')
     endif
  
  end function get_Fig4b
  !
  !
  !
  function get_dx(x1,y1,zmic,dt,gamma,c0,Va,tnext)
    real(kind=rp), intent(in) :: x1
    real(kind=rp), intent(in) :: y1
    real(kind=rp), intent(in) :: zmic
    real(kind=rp), intent(in) :: dt
    real(kind=rp), intent(in) :: gamma
    real(kind=rp), intent(in) :: c0
    real(kind=rp), intent(in) :: Va
    real(kind=rp), intent(in) :: tnext
    real(kind=rp) :: a, b, c
    real(kind=rp) :: get_dx
    !
    real(kind=rp) :: x2, y2, dr, tg, r1, sol1, sol2, attempt1, attempt2, res
    logical :: solution_found
    
    solution_found = .false.
    
    tg = tand(gamma)
    r1 = sqrt(x1**2 + y1**2 + zmic**2)

    !
    a = ((one+tg**2)*(c0**2-Va**2))/Va**2
    b = -two*(r1+c0*dt)*c0*sqrt(one+tg**2)/Va - two*x1 - two*y1*tg
    c = r1**2 + (c0**2)*(dt**2) + 2*r1*c0*dt - x1**2 - y1**2-zmic**2         !error fixed
       
    sol1 = (-b+sqrt(b**2-4.0_rp*a*c))/(two*a)
    sol2 = (-b-sqrt(b**2-4.0_rp*a*c))/(two*a)

! pick solution 
    x2 = x1+sol1 ; y2 = y1+sol1*tg ; dr = sqrt(sol1**2+(sol1*tg)**2)
    attempt1 = sqrt( x2**2 + y2**2+zmic**2 )/c0 + dr/Va               !error fixed
    res = abs(attempt1-tnext)/max(one,tnext)
    if( res.lt.10.0_rp*epsilon(one) ) then 
      get_dx = sol1
      solution_found = .true.
    end if
    
    x2 = x1+sol2 ; y2 = y1+sol2*tg ; dr = sqrt(sol2**2+(sol2*tg)**2)
    attempt2 = sqrt( x2**2 + y2**2+zmic**2 )/c0 + dr/Va             !error fixed
    res = abs(attempt2-tnext)/max(one,tnext)
    if( res.lt.10.0_rp*epsilon(one) ) then
      get_dx = sol2
      solution_found = .true.
    end if
    
    if( .not.solution_found ) call report_error('failed to establish solution','get_dx','choice_physics')
       
  end function get_dx
  !
  !
  !
  function get_theta(npt,xs,xsi)
    integer, intent(in) :: npt
    real(kind=rp), intent(in), dimension(npt) :: xs ! x_source
    real(kind=rp), intent(in), dimension(npt) :: xsi
    !
    real(kind=rp), dimension(npt) :: get_theta
    integer :: i 
    
    do i = 1, npt
      get_theta(i) = xsi(i) + (alphai(i)*pi_num)/180.0_rp 
    end do 
    
    do i = 1, npt
       if( get_theta(i) .gt. pi_num ) then
           get_theta(i) = two*pi_num - get_theta(i) 
       end if
    end do
    
  end function get_theta
  !
  !
  !
  function get_xsi(clgr,x,y,xmic,ymic)
    real(kind=rp), intent(in) :: clgr
    real(kind=rp) :: x,y,xmic,ymic
    real(kind=rp) :: get_xsi
    !
    real(kind=rp) :: fi, fid, psi
    integer :: slen
    !
    fi = atan( (xmic-x)/(y-ymic) )
    if( y.le.ymic ) fi = pi_num/two ! in reality the airframe and engines are always higher than the microphone
    psi = pi_num/two - fi 
    !
    get_xsi = psi + (pi_num/180.0_rp)*clgr  
    !    
  end function get_xsi
  !
  ! test code below:
  !
  subroutine test_functions
    integer :: uno
    !                                              also test 0.5
    integer, parameter :: nfa1 = 5
    integer, parameter :: maxPts = 20
    real(kind=rp), parameter, dimension(nfa1) :: fa1_vec = (/one,1.2_rp,1.4_rp,1.6_rp,1.8_rp/) ! various MRD:s
    real(kind=rp) :: M_trd, M_tr
    integer :: i, j
    
    call plot_trajectory()

! decomment tests below 
   call plot_Fig4a 
   call plot_Fig4b 

   call plot_Fig10a    
   call plot_Fig10b    
   
!   write(*,*) get_Fig3a(50.0_rp/2111.0_rp),  get_Fig3a(8063.5_rp/2111.0_rp)     
   call plot_Fig3a()
   
   call test_equation_A1()
       
  contains 
     subroutine plot_trajectory()
     
        uno = get_unit_number()
        open(uno,file='trajectory.m')
   
! test get_Fig4a (from M_tr 0.2 up to M_tr = 0.9)
       write(uno,'(A)') 'clf'
       write(uno,'(A)') 'x = ['
       do i = 1, n_traj_pts
          write(uno,'(F12.3,$)') x(i)
          if( i.eq.n_traj_pts ) write(uno,'(A)') ']'
       end do 
       write(uno,'(A)') 'y = ['
       do i = 1, n_traj_pts
         write(uno,'(F12.3,$)') y(i)
         if( i.eq.n_traj_pts ) write(uno,'(A)') ']'
       end do 
       !
       write(uno,'(A)') 'plot(x,y)'
       
       close(uno)
     
     end subroutine plot_trajectory
     !
     subroutine plot_Fig3a()
        real(kind=rp) :: frat
        real(kind=rp), dimension(15) :: xvec
     
       uno = get_unit_number()
       open(uno,file='testgetFig3a.m')
   
! test get_Fig4a (from M_tr 0.2 up to M_tr = 0.9)
    write(uno,'(A)') 'clf'
    write(uno,'(A)') 'x = ['
    j = 1 
    do i = 1, 8
      xvec(j) =  0.1_rp*real(2**i,rp)
      j = j + 2
    end do 
    j = 1 
    do i = 1, 7
      xvec(j+1) =  (xvec(j)+xvec(j+2))/2.0_rp
      j = j + 2
    end do 
    do i = 1, 15
        write(uno,'(F12.3,$)') xvec(i)
        if( i.eq.15 ) write(uno,'(A)') ']'
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, 15
        write(uno,'(F12.3,$)') get_Fig3a(xvec(i))
        if( i.eq.15 ) write(uno,'(A)') ']'
    end do     
    !
    write(uno,'(A)') 'semilogx(x,y)'

       close(uno)
    
    
     end subroutine plot_Fig3a
     !
     subroutine plot_Fig4a

    uno = get_unit_number()
  
    open(uno,file='testgetFig4a.m')
    
! test get_Fig4a (from M_tr 0.2 up to M_tr = 0.9)
    write(uno,'(A)') 'clf'
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (0.9_rp-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (0.9_rp-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig4a(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
       if( i.eq.1 ) write(uno,'(A)') 'hold on'
    end do 
    
! more testF4a (from M_tr 0.9 up to M_tr = M_trd)
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.9_rp + (M_trd-0.9_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.9_rp + (M_trd-0.9_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig4a(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
    end do 

    
    close(uno)
     
     
     end subroutine plot_Fig4a
     !
     ! 
     !
     subroutine plot_Fig4b

    uno = get_unit_number()
  
    open(uno,file='testgetFig4b.m')
    
! test get_F4b (from M_tr 0.2 up to M_tr = 1.0)
    write(uno,'(A)') 'clf'
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (1.0_rp-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (1.0_rp-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig4b(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
       if( i.eq.1 ) write(uno,'(A)') 'hold on'
    end do 
    
! more testF4b (from M_tr 1.0 up to M_tr = M_trd)
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 1.0_rp + (M_trd-1.0_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 1.0_rp + (M_trd-1.0_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig4b(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
    end do 

    
    close(uno)
     
     
     end subroutine plot_Fig4b
     !
     !
     !
     subroutine plot_Fig10a

    uno = get_unit_number()
  
    open(uno,file='testgetFig10a.m')
    
! test get_Fig10a (from M_tr 0.2 up to M_tr = 0.72)
    write(uno,'(A)') 'clf'
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (0.72_rp-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (0.72_rp-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig10a(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
       if( i.eq.1 ) write(uno,'(A)') 'hold on'
    end do 
    
! more testF10a (from M_tr 0.72 up to M_tr = 1.0)
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.72_rp + (one-0.72_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.72_rp + (one-0.72_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig10a(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
    end do 

    
! more testF10a (from M_tr = 1.0 up to M_tr = M_trd)
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = one + (M_trd-one)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = one + (M_trd-one)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig10a(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
    end do 
    
    close(uno)
     
     
     end subroutine plot_Fig10a
     !
     !
     !
     subroutine plot_Fig10b

    uno = get_unit_number()
  
    open(uno,file='testgetFig10b.m')
    
! test get_Fig10b (from M_tr = 0.2 up to M_tr = 1.0)
    write(uno,'(A)') 'clf'
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (one-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = 0.2_rp + (one-0.2_rp)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig10b(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
       if( i.eq.1 ) write(uno,'(A)') 'hold on'
    end do 
    
! more testF10a (from M_tr one up to M_tr = M_trd)
    write(uno,'(A)') 'x = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = one + (M_trd-one)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') M_tr
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do 
    write(uno,'(A)') 'y = ['
    do i = 1, nfa1
        M_trd = fa1_vec(i)
        do j = 1, maxPts
          M_tr = one + (M_trd-one)*float(j-1)/float(maxPts-1)
          write(uno,'(F12.3,$)') get_Fig10b(M_trd,M_tr)
        end do 
        if( i.eq.nfa1 ) then
          write(uno,'(A)') ']'
        else
          write(uno,'(A)') ' '
        end if
    end do     
    !
    do i = 1, nfa1
       write(uno,'(A)') 'semilogx(x('//trim(adjustl(inumber_string_cast(i)))//',:),y('//trim(adjustl(inumber_string_cast(i)))//',:))'
    end do 
    
    close(uno)
          
     end subroutine plot_Fig10b
     !
     !
     ! 
     subroutine test_equation_A1()
       real(kind=rp), parameter :: alpha = zero  ! level flight
       real(kind=rp), parameter :: c0 = 340.0_rp ! speed of sound
       real(kind=rp), parameter :: Va = -100.0_rp ! flight velocity 
       real(kind=rp), parameter :: dt = 0.5 ! sampling specified 
       !
       real(kind=rp) :: ta ! tangens alpha = ta
       real(kind=rp) :: r1, y1, x1, sol1, sol2
       real(kind=rp) :: a, b, c
       real(kind=rp) :: t2a,t2b, r2
       
       ta = tand(alpha) 
       
! level flight. Starting 4c0 away at 45 degrees angle 
       r1=4.0_rp*c0
       x1=r1*cosd(45.0_rp)
       y1=r1*sind(45.0_rp)
       r1 = sqrt(x1**2 + y1**2)
       
       a = ((one+ta**2)*(c0**2-Va**2))/Va**2
       b = -two*(r1+c0*dt)*c0*sqrt(one+ta**2)/Va - two*x1 - two*y1*ta
       c= r1**2 + (c0**2)*(dt**2) + 2*r1*c0*dt - x1**2 - y1**2
       
       sol1 = (-b+sqrt(b**2-4.0_rp*a*c))/(two*a)
       sol2 = (-b-sqrt(b**2-4.0_rp*a*c))/(two*a)
       
       r2 = sqrt( (x1+sol1)**2 + y1**2 )
       t2a = r2/c0 + abs(sol1)/Va

       r2 = sqrt( (x1+sol2)**2 + y1**2 )
       t2b = r2/c0 + abs(sol1)/Va
       
     end subroutine test_equation_A1
     !
     !
     !
  end subroutine test_functions
  !
  !
  !
  function get_matrix_interpolation(tobs,nfreq,nt,n_traj_pts,time,matrix)
    real(kind=rp), intent(in) :: tobs ! time onto which the new matrix should be interpolated
    integer, intent(in) :: nfreq, nt, n_traj_pts
    real(kind=rp), intent(in), dimension(n_traj_pts) :: time ! times for which aircraft data are sampled (they are stored in n_traj_pts number of matrices)
    real(kind=rp), intent(in), dimension(nfreq,nt,n_traj_pts) :: matrix
    real(kind=rp), dimension(nfreq,nt)  :: get_matrix_interpolation
    !
    integer :: i, j, k
    real(kind=rp) :: frac
    
    if( tobs.le.time(1) ) get_matrix_interpolation(1:nfreq,1:nt) = matrix(1:nfreq,1:nt,1)
    if( tobs.gt.time(n_traj_pts) ) get_matrix_interpolation(1:nfreq,1:nt) = matrix(1:nfreq,1:nt,n_traj_pts)
    
    do i = 2, n_traj_pts
        ! set interpolation pointers....
        ! otherwise
        if( tobs.gt.time(i-1) .and. tobs.le.time(i) ) then
           frac = (time(i-1)-tobs)/(time(i-1)-time(i))
           do j = 1, nfreq
             do k = 1, nt
                get_matrix_interpolation(j,k) = matrix(j,k,i)*frac + matrix(j,k,i-1)*(one-frac) 
             end do
           end do
           exit
        end if
        !
    end do
  
  end function get_matrix_interpolation
  !
  !
  !
  function get_vector_interpolation(tmic,nfreq,nt,n_traj_pts,time,vector)
    real(kind=rp), intent(in) :: tmic ! time onto which the new matrix should be interpolated
    integer, intent(in) :: nfreq, nt, n_traj_pts
    real(kind=rp), intent(in), dimension(n_traj_pts) :: time ! times for which aircraft data are sampled (they are stored in n_traj_pts number of matrices)
    real(kind=rp), intent(in), dimension(n_traj_pts) :: vector
    real(kind=rp) :: get_vector_interpolation
    !
    integer :: i, j, k
    real(kind=rp) :: frac
    
    if( tmic.le.time(1) ) get_vector_interpolation = vector(1)
    if( tmic.gt.time(n_traj_pts) ) get_vector_interpolation = vector(n_traj_pts)
    
    do i = 1, n_traj_pts-1
        if( tmic.gt.time(i) .and. tmic.le.time(i+1) ) then
           frac = (tmic-time(i))/(time(i+1)-time(i))
           do j = 1, nfreq
             do k = 1, nt
                get_vector_interpolation = (one-frac)*vector(i) + frac*vector(i+1) 
             end do
           end do
           exit
        end if
        !
    end do
  
  end function get_vector_interpolation
  !
  !
  !
  subroutine flightEffects(iptr,nfr,nth,nti,use_ground_refl,spherical_spr,atm_atten,theta,xsii_alpha,x,y,r1,tai,fobs,atm_absorption,SPLi,SPLp,prmsp)
    integer, intent(in) :: iptr
    integer, intent(in) :: nfr,nth,nti
    logical, intent(in) :: use_ground_refl,spherical_spr,atm_atten
    real(kind=rp), intent(in), dimension(nth) :: theta
    real(kind=rp), intent(in), dimension(nti) :: xsii_alpha
    real(kind=rp), intent(in), dimension(nti) :: x,y, tai
    real(kind=rp), intent(in), dimension(nti) :: r1
    real(kind=rp), intent(in), dimension(nfr,nti) :: fobs
    real(kind=rp), intent(in), dimension(nti,nfr) :: atm_absorption
    real(kind=rp), intent(in), dimension(nfr,nth,nti) :: SPLi
    real(kind=rp), intent(out), dimension(nfr,nti) :: SPLp, prmsp
    integer :: n_freqs
    integer, parameter :: N_b = 5, max_n_freqs = 150, max_ns = 15000
    real(kind=rp), dimension(max_n_freqs) :: fband,f
    real(kind=rp), dimension(max_ns) :: freq
    real(kind=rp), dimension(:,:), allocatable :: G13damp, G1mat, sub_band
    real(kind=rp) :: meanVal, term1, term2
    integer :: i, j, n_start, n_end, ptr, n_freqs_last 
    
!     call save_3D_matrix(nfr,nth,nti,SPLi,'.m') ; stop
    
! directivity. Only one directivity reaches the microphone. Here we choose the element in the directivity vector closes to the computed angle.
    allocate(source_theta(nti)) ; allocate(closest_value_ptr(nti)) 
    do i = 1, nti
       closest_value_ptr(i) = get_closest_ptr(degree(xsii_alpha(i)),nth,theta)
       source_theta(i) = theta(closest_value_ptr(i))
    end do
    

! extract the directivity elements from the SPL matrix
    do i = 1, nfr
      do j = 1, nti
         SPLp(i,j) = SPLi(i,closest_value_ptr(j),j)
         prmsp(i,j) = p0*10.0_rp**(SPLp(i,j)/20.0_rp)
      end do
    end do

!    call save_matrix(nfr,nti,SPLp,'.m') ; stop        


! atmospheric attenuation & spherical spreading
    do i = 1, nfr
       do j = 1, nti
          if( spherical_spr .and. atm_atten ) then
            SPLp(i,j)=20.0_rp*log10(prmsp(i,j)/(r1(j)*P0)) - (r1(j)/100.0_rp)*atm_absorption(j,i)
          elseif( spherical_spr ) then
            SPLp(i,j)=20.0_rp*log10(prmsp(i,j)/(r1(j)*P0))            
          elseif( atm_atten ) then
            SPLp(i,j)= SPLp(i,j) - (r1(j)/100.0_rp)*atm_absorption(j,i)            
          end if
          if( SPLp(i,j).lt.zero ) SPLp(i,j) = zero
          Prmsp(i,j)=P0*10.0_rp**(SPLp(i,j)/20.0_rp)
       end do
    end do    
!    call save_matrix(nfr,nti,Prmsp,'.m') ; stop        
    
! Groundreflection
    allocate(G1mat(max_n_freqs,nti)) ; allocate(sub_band(max_n_freqs,nti)) 
    do i = 1, nti
      !
      n_freqs = nint(three*N_b*log10(fobs(nfr,i)/fobs(1,i))/log10(two)) + 1 
      if( n_freqs.gt.max_n_freqs ) stop 'max_n_freqs insufficient in flight effects'
      !
      call set_frequencies(N_b,n_freqs,fobs(1,i),fobs(nfr,i),fband,f,freq)
      sub_band(1:n_freqs,i) = fband(1:n_freqs)
      !
      call ground_reflection(iptr,N_b,y(i),x(i),r1(i),tai(i),n_freqs,fband,G1mat(1:n_freqs,i))
      !
    end do
    
    ptr = 2
    allocate(G13damp(nfr,nti))
    do i = N_b+1, n_freqs-(N_b-1)/2, N_b
        n_start = i-(N_b-1)/2 ; n_end = i+(N_b-1)/2
        do j = 1, nti
           G13damp(ptr,j) = sum(G1mat(n_start:n_end,j))/real(n_end-n_start+1,rp)
        end do        
        ptr = ptr + 1
    end do
    ptr = ptr - 1 ! remove last unused update
    !
    G13damp(1,1:nti) = sum(G1mat(1:(N_b-1)/2+1,1))/real((N_b-1)/2+1,rp)  ! first row special
    !
    n_freqs_last = nint(three*N_b*log10(fobs(nfr,nti)/fobs(1,nti))/log10(two)) + 1 ! make sure you have the right size. Never comfirmed that n_freqs is constant in a batch 
    G13damp(nfr,1:nti) =sum(G1mat(n_freqs_last-(N_b-1)/2:n_freqs_last,1))/real((N_b-1)/2+1,rp)    ! last row special
    !
    do i = 1, nfr
       do j = 1, nti
          if( use_ground_refl ) Prmsp(i,j)=sqrt(G13damp(i,j)*Prmsp(i,j)**2)
          if( Prmsp(i,j).lt.P0 ) Prmsp(i,j) = P0
          !
          SPLp(i,j) =20.0_rp*log10(Prmsp(i,j)/P0)
          !
       end do
    end do 
    !
    deallocate(source_theta) ; deallocate(closest_value_ptr) ; deallocate(G1mat)
    deallocate(sub_band) ; deallocate(G13damp)  
          
  end subroutine flightEffects
  !
  !
  !
  subroutine ground_reflection(iptr,N_b,y_plane,x1,r1,ta,nf,fband,G1)
    integer, intent(in) :: iptr, N_b
    real(kind=rp), intent(in) :: y_plane,x1,r1,ta
    integer, intent(in) :: nf
    real(kind=rp), dimension(nf), intent(in) :: fband
    real(kind=rp), dimension(max_n_freqs), intent(out) :: G1
    !
    real(kind=rp), dimension(max_n_freqs) :: k,eta, Rvec, alfa, G, arrow, Ccap, Rvect
    complex(kind=rp), dimension(max_n_freqs) :: ny, tau, Fcap, W, gam, Amat
    !
    real(kind=rp), dimension(201) :: t
    complex(kind=rp), dimension(201) :: vec
    complex(kind=rp) :: temp
    real(kind=rp), parameter :: ainc=0.01_rp ! incoherence coefficent
    real(kind=rp), parameter :: sigma = 225000_rp ! specific flow resistance of the ground
    real(kind=rp), parameter :: R=287.05_rp 
    real(kind=rp) :: r2,y1,c,x,theta,dr,nyMod,nyArg, U, rho, Kcap, sub
    integer :: no_pts, fptr
    integer, parameter :: hundred = 100
    integer :: i, j, l
    logical :: firstPos
    !
    y1 = y_plane + ymic(iptr)
    c=sqrt(R*gamma_air*ta)
    x=abs(x1)
        
    ! microphone angle
    theta = pi_num/two - atan((y1+ymic(iptr))/x) ! likely bug. Should be "pi_num/two - atan((y_plane-y_mic)/x)"
    dr=2.0_rp*ymic(iptr)*(y1-ymic(iptr))/r1
    r2=r1+dr
    
    do i = 1, nf 
        k(i) = two*pi_num*fband(i)/c ! create the wave number k   
        eta(i) = (two*pi_num*rhoisa/sigma)*fband(i) ! dimensionless frequency
        nyMod = one+(6.86_rp*eta(i))**(-0.75_rp)
        nyArg = (4.36_rp*eta(i))**(-0.73_rp)
        ny(i)=one/CMPLX(nyMod,nyArg) ! the complex specific ground admittance
        temp = (k(i)*r2)/(two*cmplx(zero,one))
        tau(i) = sqrt(temp)*(cos(theta)+ny(i))
    end do
!    tau=((k.*r2)./(2*i)).^(0.5).*(cos(theta)+ny);

    j = 1
    do i = -hundred,hundred
       t(j) = i
       j = j + 1 
    end do
    
    do i = 1, nf
      if( abs(tau(i)).gt.10.0_rp ) then
         if( real(tau(i)).lt.zero ) then
           U=one ! U is a step function
         elseif( real(tau(i)).eq.zero ) then
           U=half   
         elseif( real(tau(i)).gt.zero ) then
           U=zero
         else 
             stop 'impossible code location'
         end if
         Fcap(i)=-two*sqrt(pi_num)*U*tau(i)*exp(tau(i)**2)+1/(2*tau(i)**2)-3/((two*tau(i)**2)**2);
      else
         do j = 1,2*hundred+1
           vec(j) = exp(-t(j)**2)/(cmplx(zero,one)*tau(i)-t(j))
         end do
         W(i)=(cmplx(zero,one)/pi_num)*sum(vec) ! Imag(z)>0 complex error function  
         Fcap(i)=one-sqrt(pi_num)*tau(i)*W(i)
      endif 
    end do
    
    !
    do i = 1, nf
       gam(i)=(cos(theta)-ny(i))/(cos(theta)+ny(i)) !  complex plane-wave reflection coefficient 
       Amat(i)=Gam(i)+(one-Gam(i))*Fcap(i) ! a complex spherical-wave reflection coefficent
       Rvec(i)=abs(Amat(i)) ! R=magnitude of the complex spherical-wave reflection coefficent  
       alfa(i)=atan2(aimag(Amat(i)),real(Amat(i))) 
       Ccap(i) = exp(-(ainc*k(i)*dr)**2)  ! coherence coefficient
    end do

    Kcap=two**(one/(6.0_rp*N_b))

    !
    do i = 1, nf
      G(i)  = one + Rvec(i)**2 + two*Rvec(i)*Ccap(i)*cos(alfa(i)+k(i)*dr)*sin((Kcap-one)*k(i)*dr)/((Kcap-one)*k(i)*dr)
      G1(i) = one + Rvec(i)**2 + two*Rvec(i)*Ccap(i)*cos(alfa(i)+k(i)*dr)*sin((Kcap-one)*k(i)*dr)/((Kcap-one)*k(i)*dr)
    end do
    
    if( y1.gt.ymic(iptr) ) then
        no_pts = 1 ; firstPos = .true.
        do i = 1, nf
           if( fband(i).gt.4000.0_rp ) then
              if( firstPos ) then
                  firstPos = .false. ; fptr = i
              end if
              Rvect(no_pts) = G(i)
              Arrow(no_pts) = Rvect(1)
              no_pts = no_pts + 1
           end if
        end do
        no_pts = no_pts - 1 
        !
        Arrow(no_pts) = 1
        sub = (arrow(1)-arrow(no_pts))/real(no_pts,rp);
        !
        j = 1
        do i = 2, no_pts-1
          arrow(i)=arrow(i)-sub*j
          j = j + 1
        end do
        !
        j = 1
        do i = fptr, nf
          G1(i)=arrow(j)  
          j = j + 1
        end do
        !
    end if
    
    if( x.eq.zero .or. theta.gt.(89.0_rp/180.0_rp)*pi_num ) then
        G(1:nf) = one ; G1(1:nf) = one
    end if
  
  end subroutine ground_reflection
  !
  !
  !
  function get_closest_ptr(val,nth,theta)
    real(kind=rp), intent(in) :: val
    integer, intent(in) :: nth
    real(kind=rp), dimension(nth) :: theta
    integer :: get_closest_ptr
    !
    real(kind=rp) :: min_distance 
    integer :: i 
    
    min_distance = huge(one)
    
    do i = 1, nth
        if( abs(theta(i)-val) .lt. min_distance ) then
            get_closest_ptr = i
            min_distance = abs(theta(i)-val)
        end if
    end do
  
  end function get_closest_ptr
  !
  ! alfa_class & alfa_mol = [dB/100 m]
  ! temp - Temperature in the atmosphere (degrees Kelvin)
  ! RH - Relative Humidity (Percent)
  !
  function get_atm_abs(nfr,t_k,relHum,fband)
    integer, intent(In) :: nfr
    real(kind=rp), intent(in) :: t_k, relHum
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), dimension(nfr) :: get_atm_abs
    ! Standard Classical Atmospheric Coefficient
    real(kind=rp), dimension(4) :: b= (/1.328924_rp, -3.179768E-2_rp, 2.173716E-4_rp, -1.7496E-6_rp/)
    ! Normalized parameters, quadratic interpolation is recommended
    integer, parameter :: no_data = 29
    real(kind=rp), parameter, dimension(no_data) :: h_norm_data=(/0.0_rp,0.25_rp,0.5_rp,0.6_rp,0.7_rp,0.8_rp,0.9_rp,1.0_rp,1.1_rp,1.2_rp,1.3_rp,1.5_rp,1.7_rp, & 
        & 2.0_rp,2.3_rp,2.5_rp,2.8_rp,3.0_rp,3.3_rp,3.6_rp,4.15_rp,4.45_rp,4.80_rp,5.25_rp,5.7_rp,6.05_rp,6.5_rp,7.0_rp,10.0_rp/)
    real(kind=rp), parameter, dimension(no_data) :: alfa_norm_data=(/0.0_rp,0.315_rp,0.7_rp,0.84_rp,0.93_rp,0.975_rp,0.996_rp,1.0_rp,0.97_rp,0.9_rp,0.840_rp,0.750_rp,& 
       & 0.670_rp,0.570_rp,0.495_rp,0.45_rp,0.4_rp,0.37_rp,0.33_rp,0.30_rp,0.26_rp,0.245_rp,0.230_rp,0.220_rp,0.210_rp,0.205_rp,0.2_rp,0.2_rp,0.2_rp/)    
    real(kind=rp), dimension(no_data) :: derivative
    integer, parameter :: lwk = 2*no_data
    real(kind=rp),  dimension(lwk) :: wk
    integer :: ierr
    !
    real(kind=rp), dimension(nfreq) :: alfa_class, alfa_molmax, alfa_mol, h_molmax, h_norm, dval, alfa_norm

    real(kind=rp) :: t_C, Bpar, ha
    logical :: spline
    integer :: i 
    !
    t_C = t_k - 273.15_rp
    !
    do i = 1, nfr
      alfa_class(i)=10.0_rp**(2.05_rp*log10(fband(i)/1000_rp) + 1.1394E-3_rp*t_C - 1.916984_rp)
      alfa_molmax(i)=10.0_rp**(log10(fband(i))+(8.42994E-3_rp)*t_C-2.755624_rp)
      Bpar=b(1)+b(2)*t_C+b(3)*t_C**2+b(4)*t_C**3
      ha = 10.0_rp**(log10(relHum)-Bpar)
      h_molmax(i)=sqrt(fband(i)/1010.0_rp)
      h_norm(i)=ha/h_molmax(i)
!      if( h_norm(i).gt.h_norm_data(no_data) )  h_norm(i) = h_norm_data(no_data) 
    end do
    ! cubic shape preserving interpolation - ala MATLAB
    spline = .false.
    call pchez ( no_data, h_norm_data, alfa_norm_data, derivative, spline, wk, lwk, ierr )
    call pchev ( no_data, h_norm_data, alfa_norm_data, derivative, nfr, h_norm, alfa_norm, dval, ierr )

    !
    do i = 1, nfr
      alfa_mol(i)=alfa_molmax(i)*alfa_norm(i)
      get_atm_abs(i)=alfa_class(i)+alfa_mol(i)    
    end do
    
  end function get_atm_abs
  !
  !
  !
  function getPNLT(nfr,nti,fband,PNL,SPL,C,s,encircled,SPLnew,snew,sbar,SPLfinal,F)
    integer, intent(in) :: nfr, nti
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), intent(in), dimension(nti) :: PNL
    real(kind=rp), intent(in), dimension(nfr,nti) :: SPL
    real(kind=rp), dimension(nti) :: C ! work vector
    real(kind=rp), dimension(nfr,nti) :: s, encircled,SPLnew,sbar,SPLfinal,F ! work arrays allocated in calling routine to avoid allocate/reallocate many times
    real(kind=rp), dimension(nfr+1,nti) :: snew ! work array allocated in calling routine to avoid allocate/reallocate many times
    real(kind=rp), dimension(nti) :: getPNLT  
    !
    integer i, k
    
    s = zero ; encircled = zero ; snew = zero ; sbar = zero ; SPLfinal = zero
    
    do i = 4, nfr
      s(i,:) = SPL(i,:)-SPL(i-1,:)
    end do
    
    do i = 2, nfr
        do k = 1, nti
            if( abs(s(i,k)-s(i-1,k)).gt.5.0_rp ) then
                if( s(i,k).gt.zero .and. s(i,k).gt.s(i-1,k) ) then
                    encircled(i,k) = one
                elseif( s(i,k).le.zero .and. s(i-1,k).gt.zero ) then
                    encircled(i-1,k) = one
                else
                    encircled(i,k) = zero
                end if
            end if
        end do
    end do
    

    SPLnew(1,:) = SPL(1,:)
    do i = 1, nfr
        do k = 1, nti
            if( encircled(i,k).eq.zero ) then
                SPLnew(i,k) = SPL(i,k)
            elseif( i.lt.24 .and. encircled(i,k).ne.zero ) then
                SPLnew(i,k) = half*(SPL(i-1,k)+SPL(i+1,k))                
            elseif( i.eq.24 .and. encircled(i,k).ne.zero ) then
                SPLnew(24,k) = SPL(23,k)+s(23,k)
            end if      
        end do
    end do

!   call save_matrix(nfr,nti,SPLnew,'.m') ; stop
! new slopes
! snew
    do i=4,nfr
      snew(i,:)=SPLnew(i,:)-SPLnew(i-1,:)
    end do
    snew(3,:)=snew(4,:)
    snew(nfr+1,:)=snew(nfr,:)
!
! sbar 
    do i = 3, nfr-1
      sbar(i,:)=(one/three)*(snew(i,:)+snew(i+1,:)+snew(i+2,:))
    end do
    
! SPLfinal    
    SPLfinal(3,:)=SPL(3,:)
    do i = 4, nfr
       SPLfinal(i,:)=SPLfinal(i-1,:)+sbar(i-1,:);
    end do

! difference
    F=SPL-SPLfinal
    !
    where( F.lt.3.0_rp ) F=zero
    F(1:2,:)=zero
    !
    do i = 3, nfr
      do k = 1, nti
        if( fband(i).lt.500.0_rp .and. fband(i).ge.50.0_rp ) then
            if( F(i,k).ge.3.0_rp .and. F(i,k).lt.20.0_rp ) then
                F(i,k)=F(i,k)/6.0_rp
            elseif( F(i,k).ge.20.0_rp ) then
                F(i,k)=10.0_rp/3.0_rp
            end if
        elseif( fband(i).ge.500.0_rp .and. fband(i).le.5000.0_rp ) then
            if( F(i,k).ge.3.0_rp .and.  F(i,k).lt.20.0_rp ) then
                F(i,k)=F(i,k)/3.0_rp
            elseif( F(i,k).ge.20.0_rp ) then
                F(i,k)=20.0_rp/3.0_rp
            end if
        elseif( fband(i).gt.5000.0_rp ) then 
            if( F(i,k).ge.3.0_rp .and. F(i,k).lt.20.0_rp ) then
                F(i,k)=F(i,k)/6.0_rp
            elseif( F(i,k).ge.20.0_rp ) then
                F(i,k)=10.0_rp/3.0_rp
            end if
        end if
      end do 
    end do


    do i = 1, nti
       C(i) = maxval(F(:,i))
       getPNLT(i) = PNL(i) + C(i)
    end do
    
  end function getPNLT 
  !
  !
  !
  function getEPNL(nti,PNLT)
    integer, intent(in) :: nti
    real(kind=rp), intent(in), dimension(nti) :: PNLT
    real(kind=rp) :: getEPNL
    !
    real(kind=rp) :: PNLTsum, PNLTM, D 
    integer :: i, istart, istop
    
    PNLTM = maxval(PNLT)
    
    if( PNLTM.gt.100.0_rp ) then
       istart=findFirstElementLtThanBnd(10.0_rp,nti,PNLTM-PNLT)
       istop =istart + findFirstElementGtThanBnd(10.0_rp,nti-istart+1,PNLTM-PNLT(istart:nti))
       PNLTsum=zero
       do i = istart, istop
          PNLTsum=PNLTsum+(10.0_rp**(PNLT(i)/10.0_rp)) 
       end do
       D = 10.0_rp*log10(PNLTsum)-PNLTM-13.0_rp;
    elseif( PNLTM.lt.90_rp ) then
       D=zero
    else
       istart=findFirstElementGtThanBnd(90.0_rp,nti,PNLT)
       istop =istart + findFirstElementLtThanBnd(90.0_rp,nti-istart+1,PNLT(istart:nti))
       !
       PNLTsum=zero
       do i = istart, istop
          PNLTsum=PNLTsum+(10.0_rp**(PNLT(i)/10.0_rp)) 
       end do
       D = 10.0_rp*log10(PNLTsum)-PNLTM-13.0_rp;
    end if
    
    if ( (D.LE.0.0_rp).AND.(abs(D).GT.abs(PNLTM-90.0_rp))) then
    D=90.0_rp-PNLTM
    end if
  
    getEPNL = D+PNLTM
    
    
  
  end function getEPNL
  !
  !
  !
  function getPNL(nfr,nti,tmic,fdoppler,SPLp,vec1,vec2,wrk)
    integer, intent(in) :: nfr ! number of frequencies
    integer, intent(in) :: nti ! number of times (for interpolated grid)
    real(kind=rp), intent(in), dimension(nti) :: tmic ! nti
    real(kind=rp), intent(in), dimension(nfr,nti) :: fdoppler ! fdoppler is fobs(1:nfreq,1:n_times) = getDopplerShift(nfreq,n_times,fband,xsii,Mai)
    real(kind=rp), intent(in), dimension(nfr,nti) :: SPLp
    real(kind=rp), dimension(nti) :: vec1, vec2 ! work arrays allocated in calling routine to avoid allocate/reallocate many times
    real(kind=rp), dimension(nfr,nti) :: wrk ! work array allocated in calling routine to avoid allocate/reallocate many times
    real(kind=rp), dimension(nti) :: getPNL 
    !
    real(kind=rp), parameter :: Ffactor = 0.15_rp 
    integer :: i 
    
    do i = 1, nti
       wrk(:,i) = getnoys(nfr,tmic(i),fdoppler(:,i),SPLp(:,i))
       vec1(i) = maxval(wrk(:,i))
    end do    
!    call save_matrix(nfr,nti,wrk,'.m') ; stop    

    do i=1,nti
      vec2(i)=(one-Ffactor)*vec1(i) + Ffactor*sum(wrk(:,i))
      getPNL(i) = 40.0_rp + (10.0_rp/(log10(two)))*log10(vec2(i))
    end do
     
  end function getPNL 
  !
  !
  !
  function getNoys(nfr,time,fband,SPL_vec) 
    integer, intent(in) :: nfr
    real(kind=rp), intent(in) :: time
    real(kind=rp), intent(in), dimension(nfr) :: fband
    real(kind=rp), intent(in), dimension(nfr) :: SPL_vec
    !
    integer, parameter :: noPars = 24
    real(kind=rp), parameter, dimension(noPars) :: Mb = (/0.043478_rp,0.040570_rp,0.036831_rp,0.036831_rp,0.035336_rp,0.033333_rp,0.033333_rp,0.032051_rp,0.030675_rp, & 
        & zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, 0.042285_rp,0.042285_rp/)
    real(kind=rp), parameter, dimension(noPars) :: Mc = (/0.030103_rp,0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, & 
        & 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, 0.030103_rp, & 
        & 0.029960_rp,0.029960_rp,0.029960_rp,0.029960_rp,0.029960_rp,0.029960_rp,0.029960_rp,0.029960_rp,0.029960_rp/)
    real(kind=rp), parameter, dimension(noPars) :: SPLa = (/91.0_rp,85.9_rp,87.3_rp,79.9_rp,79.8_rp,76.0_rp,74.0_rp,74.9_rp,94.6_rp, & 
       & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,44.3_rp,50.7_rp/)
    real(kind=rp), parameter, dimension(noPars) :: SPLb = (/64.0_rp,60.0_rp,56.0_rp,53.0_rp,51.0_rp,48.0_rp,46.0_rp,44.0_rp,42.0_rp, & 
       & zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,37.0_rp,41.0_rp/)
    real(kind=rp), parameter, dimension(noPars) :: SPLc = (/52.0_rp,51.0_rp,49.0_rp,47.0_rp,46.0_rp,45.0_rp,43.0_rp,42.0_rp,41.0_rp,& 
       & 40.0_rp,40.0_rp,40.0_rp,40.0_rp,40.0_rp,38.0_rp,34.0_rp,32.0_rp,30.0_rp,29.0_rp,29.0_rp,30.0_rp,31.0_rp,34.0_rp,37.0_rp/)
    real(kind=rp), dimension(nfr) :: getNoys
    !
    integer :: i
    !
    do i = 1, nfr
       if( fband(i).lt.400.0_rp .or. fband(i).gt.6300.0_rp) then
          if( SPL_vec(i).lt.SPLa(i) ) then
              getNoys(i)=10**(Mb(i)*(SPL_vec(i)-SPLb(i)))
          elseif( SPL_vec(i).ge.SPLa(i) ) then
              getNoys(i)=10**(Mc(i)*(SPL_vec(i)-SPLc(i)))
          end if
       elseif( fband(i).le.6300.0_rp .and. fband(i).ge.400.0_rp ) then
          getNoys(i)=10**(Mc(i)*(SPL_vec(i)-SPLc(i)))
       end if
    end do

   end function getNoys
end module choice_physics
