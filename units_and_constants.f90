module units_and_constants
  !
  use sim_precision
  !
  implicit none
  !
  real(kind=rp), parameter :: pi_num = 3.14159265358979_rp
  !
  real(kind=rp), parameter :: zero = 0.0_rp
  real(kind=rp), parameter :: half  = 0.5_rp
  real(kind=rp), parameter ::  one = 1.0_rp
  real(kind=rp), parameter ::  two = 2.0_rp
  real(kind=rp), parameter :: three = 3.0_rp
  real(kind=rp), parameter :: threepointfive = 3.5_rp
  real(kind=rp), parameter :: four = 4.0_rp
  !
  real(kind=rp), parameter, public :: pbar = 100000.0_rp
  integer, parameter, public :: max_parameter = 100
  !
  real(kind=rp), parameter :: gamma_air = 1.400_rp
  real(kind=rp), parameter :: gamma_gas = 1.333_rp
  !
  real(kind=rp), parameter :: cp_air = 1005.0_rp
  real(kind=rp), parameter :: cp_gas = 1148.0_rp
  !
  ! unit conversion
  real(kind=rp), parameter :: inch = 0.0254_rp
  real(kind=rp), parameter :: lbft_to_NM = 1.355818181_rp
  real(kind=rp), parameter :: lbm2kg = 0.45359237_rp
  real(kind=rp), parameter :: f2m = 0.3048_rp 
  real(kind=rp), parameter :: psia2pascal = 6894.757_rp
  real(kind=rp), parameter :: hp2watt = 745.7_rp
  real(kind=rp), parameter :: rankine2kelvin = 1.0_rp/1.8_rp
  real(kind=rp), parameter :: lbf2N = 4.448_rp
  real(kind=rp), parameter :: btu2J = 1055.0_rp
  real(kind=rp), parameter :: g = 9.81_rp
  real(kind=rp), parameter :: R_air = 287.0121_rp
  real(kind=rp), parameter :: Ru = 8314.5_rp !universal gas constant
  real(kind=rp), parameter :: square_inches2square_metres = (inch)**2
  real(kind=rp), parameter :: pph_to_kg_p_s = lbm2kg/3600.0_rp
  real(kind=rp), parameter :: nautical_to_meter = 1852.0_rp
  ! hot structure scaling parameter
  real(kind=rp), parameter :: hsfact = 3.36_rp ! thickness and volume factor variation parameter for tcf and tec
  !
  real(kind=rp), parameter, public :: fhval = 43200000.0_rp ! obtained from CEA at 288.15 K (Liquid Jet A)
  !
  real(kind=rp), parameter, public :: P_ref = 101325.0_rp
  real(kind=rp), parameter, public :: T_ref = 288.15_rp
  real(kind=rp), parameter :: mcorr_star = 70.0_rp
  !
  logical, parameter :: plot_discs = .true.
  logical, parameter :: plot_texts = .false.
  !
contains
  function square_inches(square_metres)
    real(kind=rp), intent(in) :: square_metres
    real(kind=rp) :: square_inches
    !
    square_inches = square_metres/square_inches2square_metres
    !
  end function square_inches
  !
  function square_metres(square_inches)
    real(kind=rp), intent(in) :: square_inches
    real(kind=rp) :: square_metres
    !
    square_metres = square_inches*square_inches2square_metres
    !
  end function square_metres
  !
  function btu(J)
    real(kind=rp), intent(in) :: J
    real(kind=rp) :: btu
    !
    btu = J/btu2J
    !
  end function btu
  !
  function ms(kts)
    real(kind=rp), intent(in) :: kts
    real(kind=rp) :: ms
    !
    ms = (1.852_rp/3.6_rp)*kts
    !
  end function ms
  !
  function J(btu)
    real(kind=rp), intent(in) :: btu
    real(kind=rp) :: J
    !
    J = btu*btu2J
    !
  end function J
  !
  function lbf(N)
    real(kind=rp), intent(in) :: N
    real(kind=rp) :: lbf
    !
    lbf = N/lbf2N
    !
  end function lbf
  !
  function N(lbf)
    real(kind=rp), intent(in) :: lbf
    real(kind=rp) :: N
    !
    N = lbf*lbf2N
    !
  end function N
  !
  function lbm(kg)
    real(kind=rp), intent(in) :: kg
    real(kind=rp) :: lbm
    !
    lbm = kg/lbm2kg
    !
  end function lbm
  !
  function kg(lbm)
    real(kind=rp), intent(in) :: lbm
    real(kind=rp) :: kg
    !
    kg = lbm*lbm2kg
    !
  end function kg
  !
  function rankine(kelvin)
    real(kind=rp), intent(in) :: kelvin
    real(kind=rp) :: rankine
    !
    rankine = kelvin/rankine2kelvin
    !
  end function rankine
  !
  function kelvin(rankine)
    real(kind=rp), intent(in) :: rankine
    real(kind=rp) :: kelvin
    !
    kelvin = rankine*rankine2kelvin
    !
  end function kelvin
  !
  function kelvin2(farenheit)
    real(kind=rp), intent(in) :: farenheit
    real(kind=rp) :: kelvin2
    !
    kelvin2 = (farenheit-32.0_rp)/1.8_rp + 273.15_rp
    !
  end function kelvin2
  !
  function kelvin3(celcius)
    real(kind=rp), intent(in) :: celcius
    real(kind=rp) :: kelvin3
    !
    kelvin3 = celcius + 273.15_rp
    !
  end function kelvin3  
  !
  function farenheit(kelvin)
    real(kind=rp), intent(in) :: kelvin
    real(kind=rp) :: farenheit
    !
    farenheit = rankine(kelvin) - (273.15_rp*1.8_rp) + 32.0_rp 
    !
  end function farenheit
  !
  function hp(watt)
    real(kind=rp), intent(in) :: watt
    real(kind=rp) :: hp
    !
    hp = watt/hp2watt
    !
  end function hp
  !
  function watt(hp)
    real(kind=rp), intent(in) :: hp
    real(kind=rp) :: watt
    !
    watt = hp*hp2watt
    !
  end function watt
  !
  function psia(pascal)
    real(kind=rp), intent(in) :: pascal
    real(kind=rp) :: psia
    !
    psia = pascal/psia2pascal
    !
  end function psia
  !
  function pascal(psia)
    real(kind=rp), intent(in) :: psia
    real(kind=rp) :: pascal
    !
    pascal = psia*psia2pascal
    !
  end function pascal
  !
  function feet(meter)
    real(kind=rp), intent(in) :: meter
    real(kind=rp) :: feet
    !
    feet = meter/f2m
    !
  end function feet
  !
  function meter(feet)
    real(kind=rp), intent(in) :: feet
    real(kind=rp) :: meter
    !
    meter = feet*f2m
    !
  end function meter
  !
  function btu_per_lbm(mj_per_kg)
    real(kind=rp), intent(in) :: mj_per_kg
    real(kind=rp) :: btu_per_lbm
    !
    btu_per_lbm = (mj_per_kg*1.0E6_rp/btu2J)*lbm2kg
    !
  end function btu_per_lbm
  !
  function kg_per_s(lbm_per_h)
    real(kind=rp), intent(in) :: lbm_per_h
    real(kind=rp) :: kg_per_s
    !
    kg_per_s = lbm_per_h * pph_to_kg_p_s
    !
  end function kg_per_s
  !
  function lbm_per_h(kg_per_s)
    real(kind=rp), intent(in) :: kg_per_s
    real(kind=rp) :: lbm_per_h
    !
    lbm_per_h = kg_per_s/pph_to_kg_p_s
    !
  end function lbm_per_h
  !
  function JPerKg(btu_per_lbm)
    real(kind=rp), intent(in) :: btu_per_lbm
    real(kind=rp) :: JPerKg
    !
    JPerKg = btu2J*(btu_per_lbm/lbm2kg)
    !
  end function JPerKg
  !
  function btuPerLbm(JPerKg)
    real(kind=rp), intent(in) :: JPerKg
    real(kind=rp) :: btuPerLbm 
    !
    btuPerLbm = (JPerKg*lbm2kg)/btu2J
    !
  end function btuPerLbm
  !
  !
  !
  function mgPerNPerS(lbm_per_hour_per_lbf)
     real(kind=rp), intent(in) :: lbm_per_hour_per_lbf
     real(kind=rp) :: mgPerNPerS
     !
     mgPerNPerS = (kg(lbm_per_hour_per_lbf)/3600.0_rp/lbf2N)*1.0E6_rp
     
  end function mgPerNPerS  
  !
  !
  !
  function farenheit_to_kelvin(fahrenheit)
     real(kind=rp), intent(in) :: fahrenheit
	 real(kind=rp) :: farenheit_to_kelvin

     farenheit_to_kelvin = ((fahrenheit-32.0_rp)*(5.0_rp/9.0_rp))+273.15_rp

  end function farenheit_to_kelvin
  !
  !
  !
  function kelvin_to_farenheit(kelvin)
     real(kind=rp), intent(in) :: kelvin
	 real(kind=rp) :: kelvin_to_farenheit

     kelvin_to_farenheit = (kelvin-273.15_rp)*(9.0_rp/5.0_rp)+32.0_rp

  end function kelvin_to_farenheit
  !
  !
  !
  function rotsp_rpm(omega)
     real(kind=rp), intent(in) :: omega
	 real(kind=rp) :: rotsp_rpm

	 rotsp_rpm = omega/(two*pi_num)*60.0_rp

  end function rotsp_rpm
  !
  !
  !
  function kfact(kilounit)
    real(kind=rp), intent(in) :: kilounit
	real(kind=rp) :: kfact

    kfact = kilounit*1000.0_rp

  end function kfact
  !
  !
  !
  function rad(d)
    real(kind=rp), intent(in) :: d
    real(kind=rp) :: rad
      
    rad = (pi_num/180.0_rp)*d
      
  end function rad
  !
  !
  !
  function degr(r)
    real(kind=rp), intent(in) :: r
    real(kind=rp) :: degr
      
    degr = (180.0_rp/pi_num)*r
      
  end function degr
    end module units_and_constants

    
    
!        spline = .false.     ! spline shape preserving interpolation - ala MATLAB
!    call pchez ( n_thetas, theta_c, OASPL, derivative, spline, wk, lwk, ierr )
!    call pchev ( n_thetas, theta_c, OASPL, derivative, nthet, theta, OASPL_p, OASPL_p_der, ierr )

!    call save_vector(nthet,OASPL_p,'.m') ; stop
    
    !    IC(1) = IBEG, desired condition at beginning of data.
!    IC(2) = IEND, desired condition at end of data.
