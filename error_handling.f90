module error_handling
  !
  ! Programmer: Anders Lundbladh
  !
  ! v1.0.0 96.09.20 Release version
  !
  ! v2.0.0 99.11.09 error_stop_thr and error_print_thr declared 
  !                 as public
  !                 character length for error_type_l and
  !                 error_prog_l increased to 150
  !                 error_handling_info added to report module name and
  !                 module version
  !
  ! v2.1.0 02.04.20 Stop statement in error_proc routine removed. Superflous
  !                 implicit none statements removed (module implicit none 
  !                 statement sets implicit none in all module internal 
  !                 subroutines).
  !
  ! v2.2.0 02.04.25 Pause statement removed.
  !
  ! Description:
  ! This is an error handling module with a global error record.
  ! 
  implicit none
  !
  private
  !
  public :: error_proc,error_report,error_pop,error_push,error_handling_info
  public :: error_stop_thr,error_print_thr
  public :: error_count_max,error_unit,error_count
  public :: error_proc_info
  public :: set_error_occurred, get_error_occurred
  !
  ! recorded error information :
  ! error_count : the total number of errors during this execution
  ! error_ls : error list size
  ! error_type_l : list of encountered error types
  ! error_prog_l : list of encountered error occurance locations
  ! error_sev_l : list of encountered error severities
  ! error control parameters :
  ! error_stop_thr : the threshold of error severity above
  !                  which the program will stop
  ! error_print_thr : the threshold of error severity above
  !                  which the errors will be reported immeadately
  ! error_count_max : the max number of errors at which the program 
  !                   will continue
  ! error_unit : unit number of the error log file
  ! error_count_stack  : pushed error_count values
  ! error_count_sp  : pointer to the most recently pushed error count
  !
  integer           :: error_count = 0
  integer           :: error_stop_thr = 2000
  integer           :: error_unit = 6 !300
  integer           :: error_count_max = 1000
  integer           :: error_print_thr = 50
  integer,parameter :: error_ls = 1000
  character(150)    :: error_type_l(error_ls+2)
  character(150)    :: error_prog_l(error_ls+2)
  integer           :: error_sev_l(error_ls+2)
  integer,parameter :: error_count_ss = 1000
  integer           :: error_count_stack(error_count_ss)
  integer           :: error_count_sp = 0
  !
  character(len=15), parameter :: version_c='v2.2.0 02.04.25'
  character(len=14), parameter :: module_name_c='error_handling'    
  !
  logical, private :: error_occurred = .false.
  !
contains 
  !
  subroutine error_handling_info(version,name)
    character(len=*), intent(out) :: version, name
    !
    version = version_c 
    name = module_name_c
    !
  end subroutine error_handling_info  
  !
  !
  !
  subroutine set_error_occurred(eo)
    logical, intent(in) :: eo

    error_occurred = eo

  end subroutine set_error_occurred
  !
  !
  !
  function get_error_occurred()
    logical :: get_error_occurred

    get_error_occurred = error_occurred

  end function get_error_occurred
  !
  !
  !
  recursive subroutine error_proc(severity,error_type,prog_unit)
    !
    ! this is the error processor
    ! to be called when a part of the program encounters an error
    !
    ! severity    :   the severity of the error
    ! error_type  :   description of the error
    ! prog_unit   :   the program unit which encountered the error
    integer, intent(in)         :: severity
    character(len=*), intent(in) :: error_type,prog_unit
    !
    !  add new error to list
    if(error_count.eq.error_ls+2) error_count=error_ls
    error_count=error_count+1
    error_type_l(error_count)=error_type
    error_prog_l(error_count)=prog_unit
    error_sev_l(error_count)=severity
    !
    if(severity.gt.error_stop_thr) then
      write(error_unit,'(a)') 
      write(error_unit,'(a)') 'Severe error encountered'
      write(error_unit,'(a)') 
      call error_report('full')
      stop
    endif
    if(severity.gt.error_print_thr) then
      write(error_unit,*) trim(error_type)
      write(error_unit,*) 'occurred in ',trim(prog_unit)
      write(error_unit,*) 'severity : ',severity
    end if
    if(error_count.gt.error_count_max) then
      write(error_unit,'(a)') 
      write(error_unit,'(a)') 'Error count limit exceeded.'
      write(error_unit,'(a)') 
      call error_report('full')
      ! pause
    endif
    if(error_count.eq.error_ls+1) then
      call error_proc(1000000,'Error list overflow','subroutine error_proc in&
           & module error_handling')
    end if
    return
  end subroutine error_proc
  !
  !
  !
  subroutine error_proc_info(mess,whereabout)
    character(len=*), intent(in) :: mess ! error message
    character(len=*), intent(in) :: whereabout ! module where error occurred
	
    write(*,'(A)') mess
    write(*,'(A)') 'occurred in '//whereabout 

  end subroutine error_proc_info
  !
  !
  !
  subroutine error_push
    !
    ! pushes the present error_count onto the error_count stack
    !
    if(error_count_sp.ge.error_count_ss) call error_proc(1000000&
         & ,'error-stack overflow','subroutine error_push in module&
         & error_handling')
    error_count_sp=error_count_sp+1
    error_count_stack(min(error_count_sp,error_count_ss))=error_count
    return
  end subroutine error_push
  !
  subroutine error_pop
    !
    ! pops the error_count from the error_count stack
    !
    if(error_count_sp.lt.1) then
      call error_proc(1000000,'error-stack underflow','subroutine&
           & error_pop in module error_handling')
    else
      error_count=error_count_stack(error_count_sp)
    endif
    error_count_sp=error_count_sp-1
    return
  end subroutine error_pop
  !
  subroutine error_report(report_type)
    !
    ! produces an error report to the error log file
    !
    ! report_type     the type of report desired, 'full' to generate a list
    !                 of all errors encountered, 'short' to get first
    !                 and worst error
    character(len=*), intent(in) :: report_type
    ! local variables
    ! sev_max   : the most severe of errors encountered
    ! sev_max_i : the index of the same
    integer i,sev_max,sev_max_i
    !
    if(error_count.gt.0) then
      ! find the most severe error      
      sev_max=0
      do i=1,error_count
        if(error_sev_l(i).gt.sev_max) then
          sev_max=error_sev_l(i)
          sev_max_i=i
        end if
      end do
      ! 
      write(error_unit,*) 
      write(error_unit,'(a)') 'Error report'
      write(error_unit,*) 
      write(error_unit,'(a,i10)') '  Total number of errors : '&
           & ,error_count
      if(report_type.eq.'short') then
        write(error_unit,'(2a)') '  The first error occurred in :&
             & ',trim(error_prog_l(1))
        write(error_unit,'(2a)') '    Error type : '&
             & ,trim(error_type_l(1))
        write(error_unit,'(a,i10)') '    Error severity : '&
             & ,error_sev_l(1)
      endif
      write(error_unit,'(2a)') '  The most severe error occurred in :&
           & ',trim(error_prog_l(sev_max_i))
      write(error_unit,'(2a)') '    Error type : '&
           & ,trim(error_type_l(sev_max_i))
      write(error_unit,'(a,i10)') '    Error severity : '&
           & ,error_sev_l(sev_max_i)
      if(report_type.eq.'full') then
        write(error_unit,'(2a)') '  Error list : '
        write(error_unit,*) 
        do i=1,error_count
          write(error_unit,*) '    ', trim(error_type_l(i)), &
               &' occurred in '&
               & ,trim(error_prog_l(i)), &
               &', severity : ',error_sev_l(i)
        end do
      endif
      write(error_unit,*) 
      write(error_unit,'(a)') 'End of error report'
      write(error_unit,*) 
    endif
    return
  end subroutine error_report
end module error_handling
