module sim_precision
  ! module programmer: Tomas Groenstedt
  !
  ! v1.0.0 98.09.10  release version
  !
  ! v2.0.0 99.06.07  * sim_precision_info added to report module name
  !                  and module version number
  !
  ! Description:
  ! The module contains an integer parameter which decides the
  ! number of significant decimals to be used in all real
  ! variables and parameters.
  ! 
  integer, parameter :: rp=selected_real_kind(12)
  character(len=15), private, parameter :: version_c='v2.0.0 99.06.07'  
  character(len=13), private, parameter :: module_name_c='sim_precision'
  !
contains
  subroutine sim_precision_info(version,name)
    character(len=*), intent(out) :: version,name
    !
    version = version_c 
    name = module_name_c 
    !
  end subroutine sim_precision_info
end module sim_precision


