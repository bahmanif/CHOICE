module choice_read_and_write
  !
  ! Module programmer: Tomas Grönstedt
  !
  ! v1.0 17.04.01  Date written.
  !
  ! v1.1 17.04.02  Parsing code included.
  !
  ! Used to read external files. Replaces GESTPAN set methods. Only used in 
  ! stand-alone version.
  !
  !
  !
  use choice_interf
  use choice_aux
  use sim_precision
  use units_and_constants
  !
  implicit none
  !
  public :: read_external_files_craw, get_mpd_craw
  public :: get_n_dim_lines_craw, get_n_weight_lines_craw, get_n_noise_lines_craw, get_dimFile_craw,  & 
    & get_configuration_craw, get_n_modules_craw, get_weightFile_craw, get_noiseFile_craw, preparse_trajectories
  ! 
  private
  ! 
  character(len=cline_len), dimension(:), allocatable :: noiseFile
  character(len=cline_len), dimension(:), allocatable :: weightFile
  character(len=cline_len), dimension(:), allocatable :: dimFile
  character(len=cline_len), dimension(:), allocatable :: perfFile  
  character(len=22), parameter ::    perf_file = 'performanceResults.txt'
  character(len=20), parameter ::     dim_file =   'dimensionsWeight.txt'
  character(len=18), parameter ::  weight_file =     'weightAircraft.txt'
  character(len=14), parameter ::   noise_file =         'inputNoise.txt'
  !
  integer :: perf_unit ! performance unit
  integer :: dim_unit ! dimensions unit
  integer :: weight_unit ! weight unit
  integer :: noise_unit
  !
  logical :: first_time = .true.
  logical :: mpde = .true. ! multiple design point data exists
  integer :: n_modules ! counted from dimensionsNoise file
  !
  integer, private :: n_lines_perf ! number of data lines in perf. file  
  integer, private :: n_lines_dim ! number of data lines in dimensions file
  integer, private :: n_lines_weight ! number of data lines in weightAircraft file
  integer, private :: n_lines_noise ! number of data lines in noise file
  character(len=7), public :: perf_type 
  !
  character(len=cline_len), dimension(:), allocatable :: modules
  !
contains      
  subroutine reset_craw()
  
    first_time = .true.
    mpde= .true.
    n_modules = 0
    modules = ' '
  
  end subroutine reset_craw
  !
  !
  !
  function get_noiseFile_craw(nld)
    integer, intent(in) :: nld 
    character(len=cline_len), dimension(nld) :: get_NoiseFile_craw
    
    get_NoiseFile_craw(1:nld) = noiseFile(1:nld)
  
  end function get_NoiseFile_craw
  !
  !
  !
  function get_weightFile_craw(nld)
    integer, intent(in) :: nld 
    character(len=cline_len), dimension(nld) :: get_weightFile_craw
    
    get_weightFile_craw(1:nld) = weightFile(1:nld)
  
  end function get_weightFile_craw
  !
  !
  !
  function get_dimFile_craw(nld)
    integer, intent(in) :: nld 
    character(len=cline_len), dimension(nld) :: get_dimFile_craw
    
    get_dimFile_craw(1:nld) = dimFile(1:nld)
  
  end function get_dimFile_craw
  !
  !
  !
  function get_n_modules_craw()
    integer :: get_n_modules_craw
    
    get_n_modules_craw = n_modules
  
  end function get_n_modules_craw
  !
  !
  !
  function get_configuration_craw(n_modules)
    integer, intent(in) :: n_modules
    character(len=cstr_len), dimension(n_modules) :: get_configuration_craw ! extraced from dim file
    
    get_configuration_craw(1:n_modules) = modules(1:n_modules)
  
  end function get_configuration_craw
  !
  !
  !
  function get_n_weight_lines_craw()
    integer :: get_n_weight_lines_craw
    
    get_n_weight_lines_craw = n_lines_weight
  
  end function get_n_weight_lines_craw
  !
  !
  !
  function get_n_dim_lines_craw()
    integer :: get_n_dim_lines_craw
    
    get_n_dim_lines_craw = n_lines_dim
  
  end function get_n_dim_lines_craw
  !
  !
  !
  subroutine preparse_trajectories(log,traj_perf,no_pts,opPnt)
     logical, intent(in) :: log
     logical, intent(in) :: traj_perf
     integer, intent(in) :: no_pts
     character(len=*), dimension(no_pts), intent(in) :: opPnt     
     logical, parameter :: logParsing = .true.
     !
     integer :: i, nBegSkip, nEndSkip
     
     write(*,'(A)') 'Preparsing routines was requested...'
     
     do i = 1, no_pts
        if( logParsing ) write(*,'(A)') '   parsing '//trim(opPnt(i))//'.txt'
        call parseTrajectoryFile(logParsing,trim(opPnt(i)),nBegSkip,nEndSkip)
        !
        if( traj_perf ) then
           if( logParsing ) write(*,'(A)') 'Trajectory performance files being parsed'
           call parsePerformanceFiles(log,trim(opPnt(i)),nBegSkip,nEndSkip)
        end if
        ! 
     end do
     
     stop 'stopped after preparsing - turn preparsing false to continue to computations'
  contains
     subroutine parseTrajectoryFile(log,opPnt,nBeg,nEnd)
        logical, intent(in) :: log
        character(len=*), intent(in) :: opPnt
        integer, intent(out) :: nBeg,nEnd
        !
        integer :: i, j, unit, npts, ntoks, xpos, Vapos, parseStart
        character(len=180) :: temp
        character(len=180), dimension(:), allocatable :: file
        character(len=32), dimension(20) :: str
        character(len=32) :: str1
        integer :: nchars 
        
        unit = get_unit_number()
        open(unit,file=trim(opPnt)//'.txt',status='old') 
    
        npts = get_n_raw_lines(unit) ! also counts comment lines
        allocate(file(npts))
        str(1:10) = ' '
   
        i = 1
        do j = 1, npts
           read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
           if( is_a_comment_line(temp) .or. is_a_data_line(trim(temp)) ) then
               file(i) = temp
               i = i + 1
               cycle
           else
               nchars = len(trim(temp))
               if( nchars.eq.0 ) write(*,'(A)') '      skipping empty line'
           end if
        end do
100     close(unit)
        npts = i - 1 ! take the new npts counted on data lines and comment lines only

! count leading and trailing line skips
        if( is_a_comment_line(file(1)) ) then ! if the file starts with a comment line try to figure out the position of x and Va (method is based on units)
           ntoks = get_n_tokens_in_string(file(1))
           str(1:ntoks) = tokenize_string(file(1),ntoks)
           xpos = 0 ; j = 0
           do i = 1, ntoks
               str1 = ' ' ; str1 = adjustl(str(1))
               if( str1(1:1).eq.'!' ) cycle
               if( index(str(i),'(m)').gt.0 ) cycle
               if( index(str(i),'(m/s)').gt.0 ) cycle
               if( index(str(i),'(s)').gt.0 ) cycle
               if( index(str(i),'(deg)').gt.0 ) cycle
               !
               j = j + 1 
               if( index(str(i),'xpos').gt.0 ) xpos = j
               if( index(str(i),'Va').gt.0 ) vapos = j
               !
           end do
        else ! if no clues from comments on first line put default positions
           xpos = 1
           vapos = 3
        end if
        !
        nBeg = 0
        do i = 1, npts 
           if( is_a_comment_line(file(i)) ) cycle
           str = ' ' ; ntoks = get_n_tokens_in_string(file(i))
           str(1:ntoks) = tokenize_string(file(i),ntoks)
           if( string_number_cast(str(xpos)).eq.zero .OR. string_number_cast(str(vapos)).eq.zero ) then
              nBeg = nBeg + 1 
              if( log ) write(*,'(A,I3,A)') 'Detected standstill point in line ',i, & 
                  & ' => removing line. CHOICE will remove corresponding lines in performance files if trajectory performance is true'
           else
              exit
           end if
        end do
! add code to estimate nEnd here  
        nEnd = 0
! 
        

! re-open unit and dump parsed file
        unit = get_unit_number()
        open(unit,file=trim(opPnt)//'.txt',status='old')         
        
        j = 0 
        do i = 1, npts
           if( is_a_comment_line(file(i)) ) then  
              write(unit,'(A)') trim(file(i)) ! format specified must match cline_len to be completely safe
           else
              j = j + 1
              if( j.le.nBeg ) then
                cycle
              else
                write(unit,'(A)') trim(file(i)) ! format specified must match cline_len to be completely safe
              end if
           end if
        end do
        
        deallocate(file)     
     
     end subroutine parseTrajectoryFile
     !
     subroutine parsePerformanceFiles(log,opPnt,nBeg,nEnd)
        logical, intent(in) :: log
        character(len=*), intent(in) :: opPnt
        integer, intent(in) :: nBeg, nEnd 
        !
        integer :: i 
        !
        do i = 1, n_modules
          if( trim(modules(i)).eq.'Fan' ) call parsePerfFile(log,opPnt,'Fan',nBeg,nEnd)
          if( trim(modules(i)).eq.'Ipc' .or.  trim(modules(i)).eq.'Lpc' ) call parsePerfFile(log,opPnt,'Lpc',nBeg,nEnd)
          if( trim(modules(i)).eq.'Lpt' ) call parsePerfFile(log,opPnt,'Lpt',nBeg,nEnd)
          if( trim(modules(i)).eq.'Comb' ) call parsePerfFile(log,opPnt,'Comb',nBeg,nEnd)
          if( trim(modules(i)).eq.'cold_nozzle' ) call parsePerfFile(log,opPnt,'cold_nozzle',nBeg,nEnd)
        end do 
        !
        call parsePerfFile(log,opPnt,'airfrm',nBeg,nEnd)
        !
     end subroutine parsePerformanceFiles
     !
     subroutine parsePerfFile(log,opPnt,modstr,nBeg,nEnd)
       logical, intent(in) :: log
       character(len=*), intent(in) :: opPnt
       character(len=*), intent(in) :: modstr
       integer, intent(in) :: nBeg,nEnd
       !
       integer :: slen, unit
       integer :: i, j, npts, ntoks, nchars
       character(len=180) :: temp
       character(len=180), dimension(:), allocatable :: file
       character(len=32), dimension(20) :: str
       character(len=32) :: str1
       character(len=80) :: fileName
       
       unit = get_unit_number() ; fileName = ' '

       slen = len(modstr)
       if( modstr(1:slen).eq.'Fan' ) fileName = trim(opPnt)//'_fan_performance.txt'
       if( modstr(1:slen).eq.'Lpc' ) fileName = trim(opPnt)//'_lpc_performance.txt' 
       if( modstr(1:slen).eq.'cold_nozzle' ) fileName = trim(opPnt)//'_coAxialJet_performance.txt' 
       if( modstr(1:slen).eq.'airfrm' ) fileName = trim(opPnt)//'_airfrm_performance.txt' 
       if( modstr(1:slen).eq.'Lpt' ) fileName = trim(opPnt)//'_lpt_performance.txt'
       if( modstr(1:slen).eq.'Comb' ) fileName = trim(opPnt)//'_comb_performance.txt' 
       
       if( log ) write(*,'(A)') '    opening file '//trim(fileName)
       open(unit,file=trim(adjustl(fileName)),status='old') 

       npts = get_n_raw_lines(unit) ! also counts comment lines
       allocate(file(npts))
       str(1:10) = ' '

! load only valid lines into file(:)
       i = 1
       do j = 1, npts
          read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
          if( is_a_comment_line(temp) .or. is_a_data_line(trim(temp)) ) then
              file(i) = temp
              i = i + 1
              cycle
           else
               nchars = len(trim(temp))
               if( nchars.eq.0 .and. log ) write(*,'(A)') '      skipping empty line'          
          end if
       end do
100    close(unit)
       npts = i - 1 ! take the new npts counted on data lines and comment lines only       
       !      
      
! re-open unit and dump parsed file
       unit = get_unit_number()
       open(unit,file=trim(adjustl(fileName)),status='old')         
        
        j = 0 
        do i = 1, npts
           if( is_a_comment_line(file(i)) ) then  
              write(unit,'(A)') trim(file(i)) ! format specified must match cline_len to be completely safe
           else
              j = j + 1
              if( j.le.nBeg ) then
                 if( log ) write(*,'(A,I3,A)') '      skipping line number ',i,' to avoid standstill operation'          
                 cycle
              else
                write(unit,'(A)') trim(file(i)) ! format specified must match cline_len to be completely safe
              end if
           end if
        end do
        
        deallocate(file)     
       
     end subroutine parsePerfFile
  end subroutine preparse_trajectories
  !
  !
  !
  function get_n_noise_lines_craw()
    integer :: get_n_noise_lines_craw
    
    get_n_noise_lines_craw = n_lines_noise
  
  end function get_n_noise_lines_craw
  !
  !
  !
  subroutine read_external_files_craw

    ! read external files (dim_file, weight_file, noise_file, perf_type are hard coded names)
    call handle_files(perf_file, dim_file, weight_file, noise_file, perf_type) ! determines performance type (perf_type)

    ! parse dim file and set n_modules & modules
    call set_modules(n_lines_dim,dimFile)
    
  end subroutine read_external_files_craw
  !
  !
  !
  function get_mpd_craw()
    type(multpt_perf_data) :: get_mpd_craw
    
    get_mpd_craw = mpd
  
  end function get_mpd_craw
  !
  !
  !
  subroutine set_modules(nld,dimF)
    integer, intent(in) :: nld 
    character(len=cline_len), dimension(nld), intent(in) :: dimF
    !
    integer :: i, j 

! first count number of modules    
    do i = 1, nld
      if( dimF(i)(1:10) .eq.'end module' ) n_modules = n_modules + 1
    end do 
    
    if( .not.allocated(modules) ) allocate(modules(n_modules))

! then set the modules
    j = 1
    do i = 1, nld
      if( dimF(i)(1:10) .eq.'end module' ) then
         modules(j) = trim(adjustl(dimF(i)(11:cline_len)))
         j = j + 1
      end if
    end do 
    
  end subroutine set_modules
  !
  !
  !
   subroutine handle_files(performance_file, dimensions_file, weight_file, noise_file, perf_type)
      character(len=*), intent(in) :: performance_file, dimensions_file, weight_file, noise_file
      character(len=7), intent(out) :: perf_type

      call retrieve_dimensions_file(dimensions_file)
      call retrieve_weight_file(weight_file)
      call retrieve_noise_file(noise_file)
      
      if( multipoint_performance_data_exist() ) then
         call open_performance_file(perf_file,perf_type) ! set perfFile
         call parse_units()
         call split_performance_file(n_lines_perf,perfFile)
         perf_type = 'multopt'
      else
          write(*,*) 'Multipoint perf. file not found. Case must be traj. perf. type.'
      end if

   end subroutine handle_files
   !
   !
   !
   subroutine parse_units()
      integer :: i 

      ! if( index(perfFile(i),'NL:').gt.0 ) call rps2rpm('NL:',perfFile(i))   <=== obsolescent
      
   end subroutine parse_units
   !
   !
   !
  subroutine retrieve_dimensions_file(dimensions_file)
    character(len=*), intent(in) :: dimensions_file
    !

    call open_file(dimensions_file,dim_unit)

    n_lines_dim  = get_n_lines(dim_unit)
    if( allocated(dimFile) ) deallocate(dimFile)
    allocate(dimFile(n_lines_dim)) ; dimFile(1:n_lines_dim) = ' '

    call load_file(dim_unit, n_lines_dim , dimFile)

    close(dim_unit)

  end subroutine retrieve_dimensions_file
   !
   !
   !
  subroutine retrieve_weight_file(weight_file)
    character(len=*), intent(in) :: weight_file
    !

    call open_file(weight_file,weight_unit)

    n_lines_weight = get_n_lines(weight_unit)
    if( allocated(weightFile) ) deallocate(weightFile)
    allocate(weightFile(n_lines_weight)) ; weightFile(1:n_lines_weight) = ' '

    call load_file(weight_unit, n_lines_weight , weightFile)

    close(weight_unit)

  end subroutine retrieve_weight_file
  !
  ! COMMENT: these three (noise, weight, dimensions) could be merged into a sigle retrieve routine
  !
  subroutine retrieve_noise_file(noise_file)
    character(len=*), intent(in) :: noise_file
    !

    call open_file(noise_file,noise_unit)

    n_lines_noise = get_n_lines(noise_unit)
    if( allocated(noiseFile) ) deallocate(noiseFile)
    allocate(noiseFile(n_lines_noise)) ; noiseFile(1:n_lines_noise) = ' '

    call load_file(noise_unit, n_lines_noise, noiseFile)

    close(noise_unit)

  end subroutine retrieve_noise_file
  !
  !
  !
  function multipoint_performance_data_exist()
    logical :: multipoint_performance_data_exist
    !
    if( file_exist(perf_file) ) then
      mpde = .true.
    else ! performanceNoise (only single 
      mpde = .false.
    end if
    
    multipoint_performance_data_exist = mpde 
    
  end function multipoint_performance_data_exist
  !
  !
  !
  function file_exist(fn)
    character(len=*), intent(in) :: fn
    logical :: file_exist
    
    inquire(file=fn,exist=file_exist) 

  end function file_exist  
  !
  !
  !
    subroutine split_performance_file(nlp,pf)
      integer, intent(in) :: nlp
      character(len=cline_len), intent(in), dimension(nlp) :: pf
      !
      integer :: i, spos, epos
      logical :: tag_found
      !
      do i = 1, nlp
         if( index(pf(i),'Point Name:').gt.0 .or. index(pf(i),'Point name:').gt.0 ) then
            tag_found = .false.
            if( index(pf(i),'ISA').gt.0 .or. index(pf(i),'SLS Hot Day').gt.0 .or. index(pf(i),'ICAO').gt.0 .or. & 
               & index(pf(i),'SLS Hot Day').gt.0 .or. index(pf(i),'design').gt.0 ) cycle
		    if( index(pf(i),'Cruise').gt.0 .or. index(pf(i),'MID-CRUISE').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_cr,mpd%perf_file_cr)
            end if
		    if( index(pf(i),'Approach').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_approach,mpd%perf_file_approach)
            end if
		    if( index(pf(i),'Cutback').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_cutback,mpd%perf_file_cutback)
            end if
		    if( index(pf(i),'Cruise').gt.0 .or. index(pf(i),'MID-CRUISE').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_cr,mpd%perf_file_cr)
            end if
            if( index(pf(i),'Top of Climb').gt.0 .or. index(pf(i),'TOP_OF_CLIMB').gt.0 .or. index(pf(i),'Top-of-climb').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_toc,mpd%perf_file_toc)
            end if
            if( index(pf(i),'Sideline').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_sl,mpd%perf_file_sl)
            end if
            if( index(pf(i),'Take-Off').gt.0 .or. index(pf(i),'Take-off').gt.0 ) then
               spos = i+3 ; tag_found = .true.
               call load_array(spos,nlp,pf,mpd%n_lines_to,mpd%perf_file_to)
            end if
            if( .not. tag_found ) call report_error('Point name specifier likely to have been misspelt. You wrote '//& 
               & trim(pf(i))//'. Try Take-off, Top-of-climb or Cruise', 'split_performance_file','read_and_write_files') 
         else
            cycle
         end if
      end do
      
      
    end subroutine split_performance_file
    !
    !
    !
    function get_perf_type()
     character(len=7) :: get_perf_type 
     !
     character(len=cline_len) :: temp
     integer :: i 

    get_perf_type = 'weight '

    do i = 1, n_lines_perf
       temp = perfFile(i)
       if( (index(temp,' W ').gt.0) .and. & 
            & (index(temp,' T ').gt.0)  .and. & 
            & (index(temp,' P ').gt.0)  .and. & 
            & (index(temp,'WRstd').gt.0) ) then 
          get_perf_type = 'gasturb'
          exit 
       end if
    end do

      end function get_perf_type
!
      !
  
!
    
  subroutine open_performance_file(performance_file,perf_type)
    character(len=*), intent(in) :: performance_file
    character(len=7), intent(out) :: perf_type

    call open_file(performance_file,perf_unit)

    n_lines_perf = get_n_lines(perf_unit)
	if( allocated(perfFile) ) deallocate(perfFile)
    allocate(perfFile(n_lines_perf)) ; perfFile(1:n_lines_perf) = ' '

    call load_file(perf_unit,n_lines_perf,perfFile)

    perf_type = get_perf_type()

    close(perf_unit)

  end subroutine open_performance_file
  !
  !
  !
  function get_var_value(mod,nam)
    character(len=*), intent(in) :: mod, nam
    real(kind=rp) :: get_var_value
    
    get_var_value = zero
  
  end function get_var_value
  !
  !
  !
    subroutine load_array(spos,nlp,pf,nl,farr)
      integer, intent(in)  :: spos
      integer, intent(in)  :: nlp
      character(len=*), intent(in), dimension(nlp) :: pf
      integer, intent(out)  :: nl
      character(len=cline_len), intent(out), dimension(:), allocatable :: farr
      !
      integer :: i, j, epos
      
      nl = 1
      do i = spos,nlp
         if( index(pf(i),':').gt.0 ) then 
            nl = nl + 1
         elseif( index(pf(i),'Isight').gt.0 ) then ! fix to handle DREAM additional variables
            nl = nl + 1
         else
            epos = i - 1
            exit
         end if
         if( i.eq.nlp ) then
!            epos = i + 1   ! line is believed to generate a problem. epos will go one step beyond the last variable
            epos = i 
            exit
         end if
      end do
      
      nl = epos-spos+1
      allocate(farr(nl)) ; farr(1:nl) = ' '
      
      j = 1
      do i = spos,epos
         farr(j) = pf(i)
         j = j + 1
      end do


    end subroutine load_array
end module choice_read_and_write
