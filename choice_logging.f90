module choice_logging
  !
  ! Module programmer: Tomas Grönstedt
  !
  ! v1.0 17.04.01  Date written.
  !
  ! choice logging routines. 
  !
  use sim_precision
  use units_and_constants
  !
  implicit none
  !
  public :: report_error, get_unit_number
  !
  private 
  !
  logical, parameter, public :: release_version = .true.
  !
  logical, private :: error_unit_open 
  integer, private :: n_error_invocations
  integer, private :: error_unit = 0    
  !
contains
  subroutine reset_choice_logging

     n_error_invocations = 0
     error_unit_open = .false.

  end subroutine reset_choice_logging
  !
  !
  !
  recursive subroutine report_error(mess,rout,mod)
    character(len=*), intent(in) :: mess ! error message
    character(len=*), intent(in) :: mod ! module where error occurred
    character(len=*), intent(in) :: rout ! routine where error occurred
    
    if( release_version ) return

    n_error_invocations = n_error_invocations + 1   
    
    if( error_unit_open ) then 	
      write(error_unit,fmt='(A)') mess
      write(error_unit,fmt='(A)') 'occurred in module '//mod 
      write(error_unit,fmt='(A)') 'in routine '//rout
    else
      error_unit_open = .true.
      error_unit = get_unit_number()
      open(error_unit,file='weicoErrorReport.txt',status='unknown')
      write(error_unit,fmt='(A)') mess
      write(error_unit,fmt='(A)') 'occurred in module '//mod 
      write(error_unit,fmt='(A)') 'in routine '//rout
    end if

  end subroutine report_error
  !
  !
  !
  function get_unit_number()
    integer :: get_unit_number
    ! locals
    logical :: opnd
    !
    get_unit_number = 242
    !
    do 
       inquire(get_unit_number,OPENED=opnd)
       if(opnd) then
          get_unit_number = get_unit_number + 10
       else 
          exit
       end if
    end do
    !
  end function get_unit_number
end module choice_logging
