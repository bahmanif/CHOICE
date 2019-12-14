module choice_aux
  !
  ! Module programmer: Tomas Grönstedt
  !
  ! v1.0 17.04.01  Date written.
  !
  ! choice auxiliary routines such as error handling and "non-physics" data 
  ! such as string lengths, common lengths of non-allocated arrays etc. 
  !
  use sim_precision
  use units_and_constants
  use interpolation
  use error_handling
  !
  implicit none
  !
  public :: reset_choice_aux, get_n_lines, get_n_raw_lines, load_file, open_file, report_error
  public :: string_number_cast, number_string_cast_fmt, multpt_perf_data, number_string_cast
  public :: set_frequencies, get_last_freq, get_unit_number, string_inumber_cast, inumber_string_cast
  public :: get_gamma, get_M, get_R, init_vector, init_ivector
  public :: save_3D_matrix, save_matrix, save_vector
  public :: findFirstElementGtThanBnd, findFirstElementGeThanBnd, findFirstElementLtThanBnd
  public :: get_gestpan_1d_table, get_gestpan_2d_table
  public :: set_vector
  public :: degree, get_p_ambient, loadStorageMat2
  public :: is_a_comment_line, is_a_data_line, tokenize_string
  public :: loadStorageMat, parseTable
  !
  public :: containsStars, preProcessFanFile, preProcessLptFile
  public :: setMachNumbers, get_Cax
  public :: gen_noise_source_matr_subr
  !
  public :: chapter3
  !
  public :: get_n_tokens_in_string
  public :: parse_xVecApproach
  public :: rad2deg
  public :: save_x_versus_y
  !
  private 
  !
  integer, parameter, public ::  cstr_len = 32 ! max string length for characters in module
  integer, parameter, public :: cline_len = 180 ! max number of characters in line 
  !
  character(len=3), parameter, public :: version_num = '1.0'
  logical, parameter, public :: release_version = .true.
  logical, parameter, public :: NOPRED_test = .false. ! run noPred test case to allow comparison / debugging
  logical, parameter, public :: testfunc = .true. ! plot component correlations to compare with reports
  !
  logical, private :: error_unit_open 
  integer, private :: n_error_invocations
  !
  integer, public :: max_n_freqs = 150 ! parameter used for local allocations in ground reflection scheme to avoid large number of repeated allocations
  !
  type multpt_perf_data
    integer :: n_lines_to, n_lines_sl, n_lines_toc, n_lines_cr, n_lines_approach, n_lines_cutback 
    character(len=cstr_len), allocatable, dimension(:) :: perf_file_to, perf_file_sl, perf_file_toc, perf_file_cr, perf_file_cutback , perf_file_approach
  end type multpt_perf_data
  type(multpt_perf_data), public :: mpd 
  !
contains
  subroutine reset_choice_aux

     ! add init code here

  end subroutine reset_choice_aux
  !
  !
  !
  subroutine loadStorageMat(fname,ndata,nvars,str,Mat)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ndata ! number of data loaded along the trajectory 
    integer, intent(in) :: nvars ! number of variables stored for particular module
    character(len=32), dimension(nvars) :: str
    real(kind=rp), intent(out), dimension(ndata,nvars) :: Mat
    !
    integer :: unit, i, j
    character(len=cline_len) :: temp    
    !
    unit = get_unit_number()
    open(unit,file=fname,status='old') 

    if( get_n_lines(unit).ne.ndata ) then
        write(*,*) get_n_lines(unit)
        stop 'confusion in number of data - ndata different from number of lines in file for file '
    end if
  
    i = 1
    do 
       read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
       if( is_a_comment_line(temp) ) cycle
       if( get_n_tokens_in_string(temp).eq.0 ) cycle ! skip blank lines 
       ! data line
       !
       str(1:nvars) = tokenize_string(temp,nvars)
       do j = 1, nvars
         Mat(i,j) = string_number_cast(str(j)) 
       end do       
       !
       if( i.eq.ndata ) exit
       i = i + 1 
       !
    end do

100 close(unit)
  
  end subroutine loadStorageMat
  !
  ! this routine is only called if file has stars
  !
  subroutine preProcessFanFile(fname,d1,a2)
    character(len=*), intent(in) :: fname
    real(kind=rp), intent(in) :: d1,a2
    !
    integer :: unit, i, j, nvars, nlines
    integer, parameter :: maxNCols = 8 ! valid for fan files with stars
    character(len=32), dimension(maxNCols) :: str
    character(len=cline_len) :: temp    
    character(len=cline_len), allocatable, dimension(:) :: newFile
    character(len=3) :: char3
    character(len=6) :: char6

    unit = get_unit_number()
    open(unit,file=fname,status='old') 
    
    nlines = get_n_raw_lines(unit)
    allocate(newFile(nlines))
    
    do i = 1, nlines
       read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
       !
       if( is_a_comment_line(temp) ) newFile(i) = ' '
       if( is_a_comment_line(temp) ) newFile(i)(1:90) = temp(1:90)
       if( is_a_comment_line(temp) ) cycle
       !
       ! data line
       nvars = get_n_tokens_in_string(temp)
       if( nvars.eq.0 ) then
           cycle 
       elseif( nvars.ne.8 ) then 
           stop 'failure in parsing logic - expected 8 number of columns in fan file' ! consistency check to not forget that this logic is specific - avoid future simulation errors
       end if
       !
       str(1:8) = tokenize_string(temp,8)
       !
       newFile(i) = parseFanLine(8,str,d1,a2)
       !
    end do 

100 close(unit)
    
! dump newFile onto file
    open(unit,file=fname,status='unknown') 
    do i = 1, nlines
       char3 = ' ' ; char3 =  adjustl(inumber_string_cast(len(trim(newFile(i)))))
       char6 = ' ' ; char6 = '(A'//adjustl(char3)//')'
       write(unit,char6) trim(newFile(i)) 
    end do 
    close(unit)
    
  contains 
    function parseFanLine(n,s,d1,a2)
       real(kind=rp), intent(in) :: d1,a2
       character(len=cline_len) :: parseFanLine
       integer, intent(in) :: n
       integer, parameter :: maxNchars = 18
       integer, parameter :: nchars = 32
       character(len=nchars), dimension(n) :: s
       character(len=32) :: varStr
       !
       integer :: i 
       !
       real(kind=rp) :: mdot, t, p, xnl
       real(kind=rp) :: Mtip, Mu, Umid
       
       xnl  = string_number_cast(s(3))
       mdot = string_number_cast(s(5))
       t    = string_number_cast(s(6))
       p    = string_number_cast(s(8))
              
       call setMachNumbers(p,t,mdot,a2,d1,xnl,gamma_air,Mtip,Mu,Umid)

       parseFanLine = ' '

       varStr = number_string_cast(Mtip)
       parseFanLine(1:18) = varStr(1:16)//'  '
       !
       varStr = number_string_cast(Mu)
       parseFanLine(19:36) = varStr(1:16)//'  '
       
       do i = 3, 5
           parseFanLine(1+(i-1)*maxNchars:1+i*maxNchars-1) = ' '//s(i)(1:maxNchars-1)
       end do 
    
    end function parseFanLine
  end subroutine preProcessFanFile
  !
  ! this routine is only called if file has stars
  !
  subroutine preProcessLptFile(fname,de,ae)
    character(len=*), intent(in) :: fname
    real(kind=rp), intent(in) :: de,ae
    !
    integer :: unit, i, j, nvars, nlines
    integer, parameter :: maxNCols = 8 ! valid for fan files with stars
    character(len=32), dimension(maxNCols) :: str
    character(len=cline_len) :: temp    
    character(len=cline_len), allocatable, dimension(:) :: newFile
    character(len=3) :: char3
    character(len=6) :: char6

    unit = get_unit_number()
    open(unit,file=fname,status='old') 
    
    nlines = get_n_raw_lines(unit)
    allocate(newFile(nlines))
    
    do i = 1, nlines
       read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
       !
       if( is_a_comment_line(temp) ) newFile(i) = ' '
       if( is_a_comment_line(temp) ) newFile(i)(1:90) = temp(1:90)
       if( is_a_comment_line(temp) ) cycle
       !
       ! data line
       nvars = get_n_tokens_in_string(temp)
       if( nvars.eq.0 ) then
           cycle 
       elseif( nvars.ne.8 ) then 
          stop 'failure in parsing logic - expected 8 number of columns in fan file' ! consistency check to not forget that this logic is specific - avoid future simulation errors
       end if
       !
       str(1:8) = tokenize_string(temp,8)
       !
       newFile(i) = parseLptLine(8,str,de,ae)
       !
    end do 

100 close(unit)
    
! dump newFile onto file
    open(unit,file=fname,status='unknown') 
    do i = 1, nlines
       char3 = ' ' ; char3 =  adjustl(inumber_string_cast(len(trim(newFile(i)))))
       char6 = ' ' ; char6 = '(A'//adjustl(char3)//')'
       write(unit,char6) trim(newFile(i)) 
    end do 
    close(unit)
    
  contains 
    function parseLptLine(n,s,de,ae)
       real(kind=rp), intent(in) :: de,ae
       character(len=cline_len) :: parseLptLine
       integer, intent(in) :: n
       integer, parameter :: maxNchars = 18
       integer, parameter :: nchars = 32
       character(len=nchars), dimension(n) :: s
       character(len=32) :: varStr
       !
       integer :: i 
       !
       real(kind=rp) :: mdot, t, p, xnl
       real(kind=rp) :: Cax, Vtr
       
       xnl  = string_number_cast(s(3))
       mdot = string_number_cast(s(4))
       t    = string_number_cast(s(6))
       p    = string_number_cast(s(7))
       
       Cax = get_Cax(gamma_gas,p,t,mdot,Ae)
       Vtr = (pi_num*de)*xnl       
       
       parseLptLine = ' '

       varStr = number_string_cast(Vtr)
       parseLptLine(1:18) = varStr(1:16)//'  '
       !       
       do i = 2, 4
           parseLptLine(1+(i-1)*maxNchars:1+i*maxNchars-1) = ' '//s(i)(1:maxNchars-1)
       end do 
       
       i = 5
       varStr = number_string_cast(Cax)
       parseLptLine(1+(i-1)*maxNchars:1+i*maxNchars-1) = varStr(1:16)//'  '
       
    
    end function parseLptLine
  end subroutine preProcessLptFile
  !
  ! checks if there are '**********' data in file (raw file from gestpan mission run)
  !
  function containsStars(fname)
    character(len=*), intent(in) :: fname
    logical :: containsStars
    !
    integer :: unit, i, j, nvars
    integer, parameter :: maxNCols = 50
    character(len=32), dimension(maxNCols) :: str
    character(len=cline_len) :: temp    

    unit = get_unit_number()
    open(unit,file=fname,status='old') 
    
    containsStars = .false.

outer: do 
       read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
       if( is_a_comment_line(temp) ) cycle
       ! data line
       !
       nvars = get_n_tokens_in_string(temp)
       !
       str(1:nvars) = tokenize_string(temp,nvars)
       
       !
       do j = 1, nvars
           if( index(str(j),'*****') ) then
               containsStars = .true.
               exit outer
           end if
       end do
       
       !
    end do outer

100 close(unit)
 
  
  end function containsStars
  !
  ! loadStorageMat2: keeps strings rather than to convert to floating point (see loadStorageMat)
  !
  subroutine loadStorageMat2(fname,ndata,nvars,str,Mat)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ndata ! number of data loaded along the trajectory 
    integer, intent(in) :: nvars ! number of variables stored for particular module
    character(len=32), dimension(nvars) :: str
    character(len=32), intent(out), dimension(ndata,nvars) :: Mat
    !
    integer :: unit, i, j
    character(len=cline_len) :: temp    
    !
    unit = get_unit_number()
    open(unit,file=fname,status='old') 

    if( get_n_lines(unit).ne.ndata ) then
        write(*,*) 'confusion reading '//fname
        stop 'confusion in number of data - ndata different from number of lines in file for file '
    end if
  
    i = 1
    do 
       read(unit,'(A180)',err=100,end=100) temp ! format specified must match cline_len to be completely safe
       if( is_a_comment_line(temp) ) cycle
       ! data line
       !
       str(1:nvars) = tokenize_string(temp,nvars)
       do j = 1, nvars
         Mat(i,j) = str(j) 
       end do       
       !
       if( i.eq.ndata ) exit
       i = i + 1 
       !
    end do

100 close(unit)
  
  end subroutine loadStorageMat2
  !
  ! set_frequencies: set third octave band frequencies
  !
  subroutine set_frequencies(nb,nfreq,fmin,fmax,fband,f,freq)
    integer, intent(in) :: nb ! number of bands
    integer, intent(in) :: nfreq ! number of frequencies
    real(kind=rp), intent(in) :: fmin, fmax
    real(kind=rp), dimension(nfreq), intent(out) :: fband
    real(kind=rp), dimension(nfreq+1), intent(out) :: f
    real(kind=rp), dimension(ceiling(fmax)-ceiling(fmin)+1), intent(out) :: freq
    !
    integer :: i, j
    real(kind=rp) :: pot
    
    f(1) = zero ; fband(1) = fmin ; pot = one/(three*float(nb)) 
    do i = 2, nfreq
      fband(i)=fband(i-1)*2**pot
      f(i)=0.5_rp*(fband(i-1)+fband(i))
    end do
    f(nfreq+1)=fband(nfreq)+half*(fband(nfreq)-fband(nfreq-1));

    j = 1
    do i = fmin,fmax
        freq(j) = i
        j = j + 1
    end do 
    
  end subroutine set_frequencies
  !
  !
  !
  function get_last_freq(nb,fmin,fmax)
    integer, intent(in) :: nb, fmin, fmax
    integer :: get_last_freq
    
    get_last_freq = nint(three*nb*log10(float(fmax)/float(fmin))/log10(two)) + 1 
  
  end function get_last_freq
  !
  !
  !
  subroutine open_file(fn,unit)
    character(len=*), intent(in) :: fn
    integer, intent(out) :: unit
    !
    logical :: file_exists
    logical, parameter :: log_open_file = .false.

    unit = get_unit_number()

    !
    inquire(file=fn,exist=file_exists)
    if( file_exists ) then 
       if( log_open_file ) write(*,*) 'Opening file ',fn
       open(unit,file=fn,status='unknown')
    else
       call report_error(& 
            & 'file '//trim(adjustl(fn))//' does not exist', & 
            & 'open_file','read_and_write_files')
    end if

  end subroutine open_file
  !
  !
  !
  subroutine load_file(unit,nl,file)
    integer, intent(in) :: unit, nl
    character(len=cline_len), intent(out), dimension(nl) :: file
    !
    character(len=cline_len) :: temp
    !
    character(len=6) :: str 
    integer :: i
    
    if( cline_len.lt.100 .or. cline_len.gt.999 ) call report_error('parameter variable string length reading not within programmed range', 'load_file','choice_aux')
    str(1:2) = '(A' ; str(3:5) = trim(inumber_string_cast(cline_len)) ; str(6:6) = ')'

    if( .not. nl.eq.0 ) then
       i = 1
       load_loop: do 
          read(unit,str) temp
          if( is_a_comment_line(temp) ) then 
             cycle load_loop
          else
             file(i) = temp
          end if
          !
          i = i + 1
          !
          if( i.eq.nl+1 ) exit load_loop
          !
       end do load_loop
    end if

  end subroutine load_file
  !
  !
  !
  function inumber_string_cast(inumber)
    integer, intent(in) :: inumber
    character(len=32) :: inumber_string_cast
    !
    inumber_string_cast = ' '
    write(inumber_string_cast(1:32),*) inumber
    !
  end function inumber_string_cast
  !
  !
  !
  function get_n_lines(unit)
    integer, intent(in) :: unit
    integer :: get_n_lines
    character(len=cline_len) :: temp
    !
    integer :: j

    j = 0
    !
    do 
       !
       read(unit,'(A150)',err=100,end=100) temp
       !
       if( is_a_comment_line(temp) ) cycle
       j = j + 1 
       !
    end do

100 rewind(unit)
    get_n_lines = j
    
  end function get_n_lines
  !
  ! like get_n_lines but counts also comment lines
  !
  function get_n_raw_lines(unit)
    integer, intent(in) :: unit
    integer :: get_n_raw_lines
    character(len=cline_len) :: temp
    !
    integer :: j

    j = 0
    !
    do 
       !
       read(unit,'(A150)',err=100,end=100) temp
       j = j + 1 
       !
    end do

100 rewind(unit)
    get_n_raw_lines = j
    
  end function get_n_raw_lines
  !
  !
  !
  function is_a_comment_line(str)
    character(len=cline_len), intent(in) :: str
    logical :: is_a_comment_line
    !
    integer :: i 

    is_a_comment_line = .false.

    do i = 1, cline_len
       if( is_a_white_space(str(i:i)) ) cycle
       if( str(i:i) .eq. '!' ) then
          is_a_comment_line = .true.
          exit
       end if
       if( is_a_character(str(i:i)) ) exit 
    end do

  end function is_a_comment_line
  !
  !
  !
  function is_a_digit(str)
     character(len=1), intent(in) :: str
     logical :: is_a_digit
     !
     integer :: i
     !
     is_a_digit = .false.
     ! 
     do i = 48,57
        if(achar(i).eq.str) is_a_digit = .true.
     enddo
	 
  end function is_a_digit
  !
  !
  !
  function is_a_data_line(str)
    character(len=*), intent(in) :: str
    logical :: is_a_data_line
    integer :: i, slen
    
    slen = len(str)
    is_a_data_line = .false.
    
    do i = 1, slen
       if( part_of_a_token(str(i:i)) ) then
           is_a_data_line = .true.
           exit
       end if
    end do
  
  end function is_a_data_line
  !
  !
  !
  function part_of_a_token(str)
    character(len=1), intent(in) :: str
    logical :: part_of_a_token
    ! locals
    integer :: i
    !
    part_of_a_token = .false.
    !
    ! check digits
    if( is_a_digit(str) ) part_of_a_token = .true.
    !
    ! check upper case and lower case
    if( is_a_character(str) ) part_of_a_token = .true.
    !
    ! check parenthesis
    if(achar(40).eq.str) part_of_a_token = .true.
    if(achar(41).eq.str) part_of_a_token = .true.
    !
    ! check underscore
    if(achar(95).eq.str) part_of_a_token = .true.
    ! 
    ! check signs
    if(achar(43).eq.str) part_of_a_token = .true. ! plus sign
    if(achar(45).eq.str) part_of_a_token = .true. ! minus sign
	!
	! check less then
	if(achar(60).eq.str) part_of_a_token = .true. ! <
	!
	! check greater then
    if(achar(62).eq.str) part_of_a_token = .true. ! >
	!
	! check equal then
	if(achar(61).eq.str) part_of_a_token = .true. ! =
	!
	! check equal then
	if(achar(58).eq.str) part_of_a_token = .true. ! :
	!
	! check equal then
	if(achar(47).eq.str) part_of_a_token = .true. ! /
    !
	! check star then
	if(achar(42).eq.str) part_of_a_token = .true. ! *
    !
  end function part_of_a_token
  !
  !
  !
  function is_a_character(str)
    character(len=1), intent(in) :: str
    logical :: is_a_character
    !
    integer :: i
    !
    is_a_character = .false.
    ! check upper case
    do i = 65,90
       if(achar(i).eq.str) is_a_character = .true.
    enddo
    ! check lower case
    do i = 97,122
       if(achar(i).eq.str) is_a_character = .true.
    enddo
    
  end function is_a_character
  !
  !
  !
  function is_a_white_space(string)
    ! true if string is only whitespace
    character(len=*) :: string
    logical :: is_a_white_space
    ! locals
    integer :: i

    is_a_white_space = .true.
    !
    do i = 1,len(string)
       if( .not.white_space(string(i:i))) then
          is_a_white_space = .false.
          exit
       end if
    end do
    !
  end function is_a_white_space
  !
  !
  !
  function white_space(str)
    character(len=1), intent(in) :: str
    logical :: white_space
    !
    white_space = .false.
    ! check if tab
    if(str.eq.achar(9)) white_space = .true.
    ! check if space
    if(str.eq.achar(32)) white_space = .true.
    !
  end function white_space
  !
  !
  !
  function string_number_cast(str)
    character(len=*), intent(in) :: str
    real(kind=rp) :: string_number_cast
    ! locals
    integer :: i, j
    integer :: np, sp, ep
    !
    integer :: ifact
    logical :: coeff_set 
    real(kind=rp) :: coeff

    np = len(str)
    coeff_set = .false.

    set_start_and_end: do i = 1, np
       if( part_of_a_digit(str(i:i)) ) then
          sp = i

          do j = i, np
             if(  .not. part_of_a_digit(str(j:j)) ) then
                if(  .not. (str(j:j).eq.'e' .or. str(j:j).eq.'E' ) ) then
                   ep = j - 1
                   exit set_start_and_end
                else
                   coeff_set = .true.
                   ep = j - 1
                   if( str(j+1:j+1) .eq. '+' ) then
                      read(str(j+2:j+4),'(I3)')  ifact
                      coeff = 10.0_rp**real(ifact)
                   elseif( str(j+1:j+1) .eq. '-' ) then
                      read(str(j+2:j+4),'(I3)') ifact
                      coeff = 10.0_rp**real(-ifact)
                   end if
                   !
                   exit set_start_and_end
                   !
                end if
             end if
          end do
       end if
    end do set_start_and_end
    ! transfer variable

    if( coeff_set ) then
       read(str(sp:ep),*) string_number_cast
       string_number_cast = string_number_cast*coeff
    else
       read(str(sp:ep),*) string_number_cast
    end if

  end function string_number_cast
  !
  !
  !
  function part_of_a_digit(str)
    character(len=1), intent(in) :: str
    logical :: part_of_a_digit
    !
    part_of_a_digit = .false.
    !
    if(str(1:1).eq.'0') part_of_a_digit = .true.
    if(str(1:1).eq.'1') part_of_a_digit = .true.
    if(str(1:1).eq.'2') part_of_a_digit = .true.
    if(str(1:1).eq.'3') part_of_a_digit = .true.
    if(str(1:1).eq.'4') part_of_a_digit = .true.
    if(str(1:1).eq.'5') part_of_a_digit = .true.
    if(str(1:1).eq.'6') part_of_a_digit = .true.
    if(str(1:1).eq.'7') part_of_a_digit = .true.
    if(str(1:1).eq.'8') part_of_a_digit = .true.
    if(str(1:1).eq.'9') part_of_a_digit = .true.
    if(str(1:1).eq.'E') part_of_a_digit = .true.
    if(str(1:1).eq.'+') part_of_a_digit = .true.
    if(str(1:1).eq.'-') part_of_a_digit = .true.
    if(str(1:1).eq.'.') part_of_a_digit = .true.
    !
  end function part_of_a_digit
  !
  !
  !
  function number_string_cast(number)
    real(kind=rp), intent(in) :: number
    character(len=32) :: number_string_cast
    !
    number_string_cast = ' '
    write(number_string_cast(1:32),*) number
    !
  end function number_string_cast
  !
  !
  !
  function number_string_cast_fmt(number,fmt)
    real(kind=rp), intent(in) :: number
    integer, intent(In) :: fmt
    character(len=fmt) :: number_string_cast_fmt
    !
    number_string_cast_fmt(1:fmt) = ' '
    write(number_string_cast_fmt(1:fmt),'(A18)') number
    !
  end function number_string_cast_fmt
  !
  !
  !
  function string_inumber_cast(str)
    character(len=*), intent(in) :: str
    integer :: string_inumber_cast
    ! locals
    integer :: i, j
    integer :: np, sp, ep

    ! transfer variable
    read(str(1:2),*) string_inumber_cast

  end function string_inumber_cast
  !
  !
  !
  !
  !  Derived from zeroin function downloaded from http://www.netlib.org
  !  Routine outputs approximate root approximating a zero of  f  in the interval ax,bx
  !
  !  It is assumed  that f(ax) and f(bx) have  opposite signs.
  !  This is checked, and an error message is printed if this is not
  !  satisfied.   zeroin  returns a zero x in the given interval
  !  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
  !  the  relative machine precision defined as the smallest representable
  !  number such that  1.+macheps .gt. 1.
  !	
  !  this function subprogram is a slightly  modified  translation  of
  !  the algol 60 procedure  zero  given in Richard Brent, algorithms for
  !  minimization without derivatives, prentice-hall, inc. (1973).
  !
  recursive function zeroin(ax,bx,f,tol,rpar,bracketed)
    real(kind=rp), intent(in) :: ax ! left endpoint of initial interval
    real(kind=rp), intent(in) :: bx ! right endpoint of initial interval
    real(kind=rp) :: f ! function subprogram which evaluates f(x) for any x in the interval  ax,bx
    real(kind=rp), intent(in) :: tol ! desired length of the interval of uncertainty of the final result (.ge.0.)
    real(kind=rp), intent(in), optional, dimension(*) :: rpar
    real(kind=rp) :: zeroin
    real(kind=rp) :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
    logical, intent(out) :: bracketed 
    !
    external f
    !
    bracketed = .true.
10  eps = epsilon(1.0_rp)
    tol1 = eps+1.0d0
    !
    a=ax
    b=bx
    !
    if( present(rpar) ) then
       fa=f(a,rpar)
    else
       fa=f(a)
    end if
    !
    if( present(rpar) ) then
       fb=f(b,rpar)
    else
       fb=f(b)
    end if
    !
    ! check that f(ax) and f(bx) have different signs
    if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
    if (fa * (fb/abs(fb)) .le. 0.0d0) go to 20
     bracketed = .false.
20  c=a
    fc=fa
    d=b-a
    e=d
30  if (abs(fc).ge.abs(fb)) go to 40
    a=b
    b=c
    c=a
    fa=fb
    fb=fc
    fc=fa
40  tol1=2.0d0*eps*abs(b)+0.5d0*tol
    xm = 0.5d0*(c-b)
    if ((abs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
    !
    ! see if a bisection is forced
    !
    if ((abs(e).ge.tol1).and.(abs(fa).gt.abs(fb))) go to 50
    d=xm
    e=d
    go to 110
50  s=fb/fa
    if (a.ne.c) go to 60
    !
    ! linear interpolation
    !
    p=2.0d0*xm*s
    q=1.0d0-s
    go to 70
    !
    ! inverse quadratic interpolation
    !
60  q=fa/fc
    r=fb/fc
    p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
    q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
70  if (p.le.0.0d0) go to 80
    q=-q
    go to 90
80  p=-p
90  s=e
    e=d
    if (((2.0d0*p).ge.(3.0d0*xm*q-abs(tol1*q))).or.(p.ge.abs(0.5d0*s*q))) go to 100
    d=p/q
    go to 110
100 d=xm
    e=d
110 a=b
    fa=fb
    if (abs(d).le.tol1) go to 120
    b=b+d
    go to 140
120 if (xm.le.0.0d0) go to 130
    b=b+tol1
    go to 140
130 b=b-tol1
140 if( present(rpar) ) then 
      fb=f(b,rpar)
	else
	  fb=f(b)
	end if
    if ((fb*(fc/abs(fc))).gt.0.0d0) go to 20
    go to 30
150 zeroin=b
    return
  end function zeroin
  !
  !
  !
  function get_R(fa)
    real(kind=rp), intent(in) :: fa
    real(kind=rp) :: get_R
    !
    get_R = R_air/(fa + one)
    !    
  end function get_R
  !
  !
  !
  function get_M(gam,xfunc0)
    real(kind=rp), intent(in) :: gam
    real(kind=rp), intent(in) :: xfunc0
    !
    real(kind=rp) :: get_M
    real(kind=rp), dimension(2) :: rpar
    logical :: bracketed = .true.
    
    rpar(1) = gam
    rpar(2) = xfunc0
    
    get_M = zeroin(zero,one,xfunc_err,epsilon(1.0_rp),rpar,bracketed)
    if(.not.bracketed) then 
	  call report_error('failed to bracket M','get_M','physics')
    end if
    
  end function get_M
  !
  !
  !
  function xfunc_err(M,rpar)
    real(kind=rp), intent(in) :: M
	real(kind=rp), intent(in), dimension(2) :: rpar
    !
    real(kind=rp) :: xfunc_err,gam,xfunc0
    !
    xfunc_err = xfunc(M,rpar(1)) - rpar(2)
    !
    return
  end function xfunc_err
  !
  !
  !
  function xfunc(M,gam)
    real(kind=rp), intent(in) :: gam, M
    real(kind=rp) :: xfunc 
    !
    xfunc = sqrt(gam)*M/((one+(gam-one)*M*M/two)&
         & **((gam+one)/(two*(gam-one))))
    !
  end function xfunc  
  !
  !
  !
  function get_gamma(arg,fa)
    real(kind=rp), intent(in) :: arg, fa
    real(kind=rp) :: get_gamma 

    if( fa.eq.zero ) then 
       get_gamma = get_gamma_air25_nd(arg) ! nd = non-dissociated
    else
       get_gamma = get_gamma25_nd(arg,fa) ! nd = non-dissociated
    end if
    
  end function get_gamma
  !
  !
  !
  function get_gamma_air25_nd(x)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: get_gamma_air25_nd
    ! smoothness parameter s =           0.00000001000 was used
 ! smoothness parameter s =           0.00000050000 was used
 integer, parameter :: n_c3 =   10
 real(kind=rp), dimension(n_c3), parameter :: c3 = (/ & 
      &   0.140101004522E+01_rp,  0.140187077892E+01_rp,  0.139909670543E+01_rp,  &
      &   0.138211328907E+01_rp,  0.135816186142E+01_rp,  0.132502254773E+01_rp,  &
      &   0.130732269934E+01_rp,  0.129264508483E+01_rp,  0.128777797052E+01_rp,  &
      &   0.128483089906E+01_rp/)
 integer, parameter :: n_t3 =  14
 real(kind=rp), dimension(n_t3),  parameter :: t3 = (/ & 
      &   0.200000000000E+03_rp,  0.200000000000E+03_rp,  0.200000000000E+03_rp,  &
      &   0.200000000000E+03_rp,  0.350000000000E+03_rp,  0.550000000000E+03_rp,  &
      &   0.750000000000E+03_rp,  0.950000000000E+03_rp,  0.160000000000E+04_rp,  &
      &   0.210000000000E+04_rp,  0.300000000000E+04_rp,  0.300000000000E+04_rp,  &
      &   0.300000000000E+04_rp,  0.300000000000E+04_rp/)
    ! local variables 
    integer, parameter :: lwrk =  27
    real(kind=rp) :: wrk(lwrk),z
    integer :: idump
    ! 
    call splev(t3,n_t3,c3,3,x,z,wrk,idump)
    ! 
    get_gamma_air25_nd = z
    !
    return
  end function get_gamma_air25_nd  
  !
  !
  !
  function get_gamma25_nd(x,y)
    real(kind=rp), intent(in) :: x,y
    real(kind=rp) :: get_gamma25_nd
    !
 ! smoothness parameter s =           0.00004000000 was used
 integer, parameter :: n_c3 =   99
 real(kind=rp), dimension(n_c3), parameter :: c3 = (/ & 
      &   0.140100559970E+01_rp,  0.139660650376E+01_rp,  0.139191939788E+01_rp,  &
      &   0.138709494998E+01_rp,  0.138556621199E+01_rp,  0.138453349197E+01_rp,  &
      &   0.138355638473E+01_rp,  0.138312584583E+01_rp,  0.138288226211E+01_rp,  &
      &   0.140190242493E+01_rp,  0.139578548894E+01_rp,  0.138938458824E+01_rp,  &
      &   0.138286675556E+01_rp,  0.138081775125E+01_rp,  0.137944747535E+01_rp,  &
      &   0.137810596285E+01_rp,  0.137750371460E+01_rp,  0.137715917676E+01_rp,  &
      &   0.139910004870E+01_rp,  0.139005537288E+01_rp,  0.138080335618E+01_rp,  &
      &   0.137150603102E+01_rp,  0.136861157755E+01_rp,  0.136669997218E+01_rp,  &
      &   0.136474922235E+01_rp,  0.136385451262E+01_rp,  0.136333630592E+01_rp,  &
      &   0.138206463493E+01_rp,  0.137130985369E+01_rp,  0.136054704095E+01_rp,  &
      &   0.134985699050E+01_rp,  0.134655856132E+01_rp,  0.134440554670E+01_rp,  &
      &   0.134212227012E+01_rp,  0.134105529659E+01_rp,  0.134043083789E+01_rp,  &
      &   0.135824346779E+01_rp,  0.134724190829E+01_rp,  0.133635962316E+01_rp,  &
      &   0.132560451401E+01_rp,  0.132230333975E+01_rp,  0.132016370536E+01_rp,  &
      &   0.131783830765E+01_rp,  0.131673931333E+01_rp,  0.131609215960E+01_rp,  &
      &   0.132495719850E+01_rp,  0.131362373143E+01_rp,  0.130256033140E+01_rp,  &
      &   0.129169534669E+01_rp,  0.128837959653E+01_rp,  0.128624646268E+01_rp,  &
      &   0.128387036280E+01_rp,  0.128273504330E+01_rp,  0.128206257872E+01_rp,  &
      &   0.130740134630E+01_rp,  0.129526178726E+01_rp,  0.128356778218E+01_rp,  &
      &   0.127224678769E+01_rp,  0.126878718125E+01_rp,  0.126655265080E+01_rp,  &
      &   0.126414529057E+01_rp,  0.126301137667E+01_rp,  0.126234472124E+01_rp,  &
      &   0.129581328243E+01_rp,  0.128351167904E+01_rp,  0.127168644419E+01_rp,  &
      &   0.126037308501E+01_rp,  0.125688628576E+01_rp,  0.125460061284E+01_rp,  &
      &   0.125232428302E+01_rp,  0.125129138000E+01_rp,  0.125069642531E+01_rp,  &
      &   0.128944095545E+01_rp,  0.127715922059E+01_rp,  0.126545608369E+01_rp,  &
      &   0.125418087983E+01_rp,  0.125074941864E+01_rp,  0.124854528226E+01_rp,  &
      &   0.124612824478E+01_rp,  0.124498058881E+01_rp,  0.124430292772E+01_rp,  &
      &   0.128621694999E+01_rp,  0.127409783555E+01_rp,  0.126248007654E+01_rp,  &
      &   0.125137827120E+01_rp,  0.124796053634E+01_rp,  0.124572335438E+01_rp,  &
      &   0.124348405346E+01_rp,  0.124246531999E+01_rp,  0.124187765768E+01_rp,  &
      &   0.128489608970E+01_rp,  0.127279636487E+01_rp,  0.126132385380E+01_rp,  &
      &   0.125019226957E+01_rp,  0.124683835484E+01_rp,  0.124471954011E+01_rp,  &
      &   0.124221196276E+01_rp,  0.124098304195E+01_rp,  0.124024548003E+01_rp/)
 integer, parameter :: n_tx3 =  15
 real(kind=rp), dimension(n_tx3),  parameter :: tx3 = (/ & 
      &   0.200000000000E+03_rp,  0.200000000000E+03_rp,  0.200000000000E+03_rp,  &
      &   0.200000000000E+03_rp,  0.350000000000E+03_rp,  0.550000000000E+03_rp,  &
      &   0.750000000000E+03_rp,  0.950000000000E+03_rp,  0.160000000000E+04_rp,  &
      &   0.210000000000E+04_rp,  0.255000000000E+04_rp,  0.300000000000E+04_rp,  &
      &   0.300000000000E+04_rp,  0.300000000000E+04_rp,  0.300000000000E+04_rp/)
 integer, parameter :: n_ty3 =  13
 real(kind=rp), dimension(n_ty3),  parameter :: ty3 = (/ & 
      &   0.000000000000E+00_rp,  0.000000000000E+00_rp,  0.000000000000E+00_rp,  &
      &   0.000000000000E+00_rp,  0.450000000000E-01_rp,  0.525000000000E-01_rp,  &
      &   0.575000000000E-01_rp,  0.640000000000E-01_rp,  0.660000000000E-01_rp,  &
      &   0.682000000000E-01_rp,  0.682000000000E-01_rp,  0.682000000000E-01_rp,  &
      &   0.682000000000E-01_rp/)
    ! local variables 
    integer, parameter :: lwrk = 12
    real(kind=rp) :: wrk(lwrk),z
    integer, parameter :: kwrk = 2
    integer :: iwrk(kwrk),idump
    ! 
    call bispev(tx3,n_tx3,ty3,n_ty3,c3,3,3,x,y,z,wrk,lwrk,iwrk,lwrk,idump)
    ! 
    get_gamma25_nd = z
    !
    return
  end function get_gamma25_nd
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
  !
  !
  !
  function rad2deg(rad)
    real(kind=rp), intent(in) :: rad
    real(kind=rp) :: rad2deg
    
    rad2deg = rad*(180.0_rp/pi_num)
  
  end function rad2deg
  !
  ! vector initializer fixed alloc. style
  !
  subroutine set_vector(xs,xe,st,nd,vect)
    real(kind=rp), intent(in) :: xs ! start x
    real(kind=rp), intent(in) :: xe ! end x
    real(kind=rp), intent(in) :: st ! x step
    integer, intent(in) :: nd ! number of data
    real(kind=rp), dimension(nd) :: vect
    !
    integer :: i
    
    do i = 1, nd
       vect(i) = xs + float(i-1)*st
    end do
  
  end subroutine set_vector   
  !
  ! vector initializer "MATLAB" style
  !
  subroutine init_vector(xs,xe,st,nd,vect)
    real(kind=rp), intent(in) :: xs ! start x
    real(kind=rp), intent(in) :: xe ! end x
    real(kind=rp), intent(in) :: st ! x step
    integer, intent(out) :: nd ! number of data
    real(kind=rp), dimension(:), allocatable, intent(out) :: vect
    !
    integer :: i
    
    nd = floor((xe-xs)/st) + 1
    allocate(vect(nd)) 

    do i = 1, nd
       vect(i) = xs + float(i-1)*st
    end do
  
  end subroutine init_vector 
  !
  !
  !
  function parse_xVecApproach(npts,x)
    integer, intent(in) :: npts 
    real(kind=rp), dimension(npts) :: x
    real(kind=rp), dimension(npts) :: parse_xVecApproach
    !
    integer :: i 
    !
    
    do i = 1, npts
       parse_xVecApproach(i) = x(i) - x(npts)
    end do 
    
  end function parse_xVecApproach
  !
  ! vector initializer "MATLAB" style
  !
  subroutine init_ivector(xs,xe,st,nd,vect)
    integer, intent(in) :: xs ! start x
    integer, intent(in) :: xe ! end x
    integer, intent(in) :: st ! x step
    integer, intent(out) :: nd ! number of data
    integer, dimension(:), allocatable :: vect
    !
    integer :: i
    
    nd = (xe-xs)/st + 1
    allocate(vect(nd)) 

    do i = 1, nd
       vect(i) = xs + (i-1)*st
    end do
  
  end subroutine init_ivector 
  !
  !
  !
  subroutine save_vector(n_rows,vec,ext)
       integer, intent(in) :: n_rows ! number of rows
       real(kind=rp), intent(in), dimension(n_rows) :: vec
       character(len=*), intent(in), optional :: ext
       !
       integer :: unit, i
              
       unit = get_unit_number()
       if( present(ext) ) then
         open(unit,file='choiceVector'//ext,status='unknown')
       else
         open(unit,file='choiceVector',status='unknown')
       end if
       
     
       write(unit,'(A,$)') 'chVec = ['
       
       do i = 1, n_rows
          write(unit,'(F20.10,$)') vec(i)
       end do
       
       write(unit,'(A)') ']'
     
       close(unit)
       
  end subroutine save_vector
  !
  !
  !
  function get_gestpan_1d_table(n,xval,yval)
    integer, intent(in) :: n
    real(kind=rp), intent(in), dimension(n) :: xval, yval
    real(kind=rp), dimension(2+2*n) :: get_gestpan_1d_table
    
    get_gestpan_1d_table(1) = n
    get_gestpan_1d_table(2) = 0
    
    get_gestpan_1d_table(3:2+n) = xval(1:n) 
    get_gestpan_1d_table(3+n:2+2*n) = yval(1:n) 
    
  end function get_gestpan_1d_table
  !
  !
  !
  function get_gestpan_2d_table(nx,xval,ny,yval,zval)
    integer, intent(in) :: nx
    real(kind=rp), intent(in), dimension(nx) :: xval
    integer, intent(in) :: ny
    real(kind=rp), intent(in), dimension(ny) :: yval
    real(kind=rp), intent(in), dimension(nx,ny) :: zval
    real(kind=rp), dimension(2+nx+ny+nx*ny) :: get_gestpan_2d_table
    !
    integer :: i 
    
    get_gestpan_2d_table(1) = nx
    get_gestpan_2d_table(2) = ny
    
    get_gestpan_2d_table(3:2+nx) = xval(1:nx) 
    get_gestpan_2d_table(3+nx:2+nx+ny) = yval(1:ny) 
    
    do i = 1, nx
       get_gestpan_2d_table(3+nx+ny+ny*(i-1):2+nx+ny+ny*i) = zval(i,:)
    end do
    
  end function get_gestpan_2d_table
  !
  !
  !
  function findFirstElementLtThanBnd(val,n,vec)
    real(kind=rp), intent(in) :: val
    integer, intent(in) :: n
    real(kind=rp), intent(in), dimension(n+1) :: vec
    integer :: findFirstElementLtThanBnd
    !
    integer :: i 
    
    do i = 1, n+1
       if( vec(i).lt.val ) then
          findFirstElementLtThanBnd = i
          exit
       end if
    end do
  
  end function findFirstElementLtThanBnd
  !
  !
  !
  function findFirstElementGtThanBnd(val,n,vec)
    real(kind=rp), intent(in) :: val
    integer, intent(in) :: n
    real(kind=rp), intent(in), dimension(n+1) :: vec
    integer :: findFirstElementGtThanBnd
    !
    integer :: i 
    findFirstElementGtThanBnd = 0 ! means no value found
    
    do i = 1, n+1
       if( vec(i).gt.val ) then
          findFirstElementGtThanBnd = i
          exit
       end if
    end do
  
  end function findFirstElementGtThanBnd
  !
  !
  !
  !
  !
  !
  function findFirstElementGeThanBnd(val,n,vec)
    real(kind=rp), intent(in) :: val
    integer, intent(in) :: n
    real(kind=rp), intent(in), dimension(n+1) :: vec
    integer :: findFirstElementGeThanBnd
    !
    integer :: i 
    findFirstElementGeThanBnd = 0 ! means no value found
    
    do i = 1, n+1
       if( vec(i).ge.val ) then
          findFirstElementGeThanBnd = i
          exit
       end if
    end do
  
  end function findFirstElementGeThanBnd
  !
  !
  !
  subroutine gen_noise_source_matr_subr(operatingPoint,nfreq,nthet,n_traj_pts,p0,inlet_fan,discharge_fan,inlet_Lpc, & 
    & inlet_tone_fan, discharge_tone_fan, inlet_broadband_fan, discharge_broadband_fan, inlet_combination_fan, & 
    & inlet_tone_Lpc, inlet_broadband_Lpc, inlet_combination_Lpc, & 
    & Lpt,Comb,Caj,Airfrm)
    !
    character(len=*), intent(in) :: operatingPoint
    integer, intent(in) :: nfreq, nthet, n_traj_pts
    real(kind=rp), intent(in) :: p0
    real(kind=rp), dimension(nfreq,nthet,n_traj_pts), intent(in) :: inlet_tone_fan, discharge_tone_fan, & 
      & inlet_broadband_fan, discharge_broadband_fan, inlet_combination_fan, inlet_tone_Lpc, inlet_broadband_Lpc, inlet_combination_Lpc
    real(kind=rp), dimension(nfreq,nthet,n_traj_pts), intent(in) :: inlet_fan,discharge_fan,inlet_Lpc,Lpt,Comb,Caj,Airfrm
    !
    real(kind=rp), allocatable, dimension(:,:,:) :: temp3DMatrix
    
    allocate(temp3DMatrix(nfreq,nthet,n_traj_pts))
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanInlet',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
!    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,discharge_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanDischarge',nfreq,nthet,n_traj_pts,temp3DMatrix)    
    
!    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_Lpc,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_lpcInlet',nfreq,nthet,n_traj_pts,temp3DMatrix)

!    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,Lpt,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_Lpt',nfreq,nthet,n_traj_pts,temp3DMatrix)

!    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,Comb,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_Comb',nfreq,nthet,n_traj_pts,temp3DMatrix)

!    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,Caj,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_Caj',nfreq,nthet,n_traj_pts,temp3DMatrix) 

!    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,Airfrm,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_Airfrm',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
! implementation for URT
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_tone_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanInletTone',nfreq,nthet,n_traj_pts,temp3DMatrix)

    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,discharge_tone_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanDischargeTone',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_broadband_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanInletBroadband',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,discharge_broadband_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanDischargeBroadband',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_combination_fan,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_fanInletCombination',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_tone_Lpc,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_LpcInletTone',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_broadband_Lpc,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_LpcInletBroadband',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    temp3DMatrix = zero ; call castPrms2SPL(nfreq,nthet,n_traj_pts,p0,inlet_combination_Lpc,temp3DMatrix)
    call save_3D_matrix('choiceOutput/'//trim(operatingPoint)//'_LpcInletCombination',nfreq,nthet,n_traj_pts,temp3DMatrix)
    
    
    deallocate(temp3DMatrix) 
      
  end subroutine gen_noise_source_matr_subr
  !
  !
  !
  subroutine castPrms2SPL(nfr,nth,nti,p0,prms,castSPL)
    integer, intent(in) :: nfr, nth, nti
    real(kind=rp), intent(in) :: p0 
    real(kind=rp), intent(in), dimension(nfr,nth,nti) :: prms
    real(kind=rp), intent(out), dimension(nfr,nth,nti) :: castSPL
    !
    integer :: i, j, k

    do i = 1, nfr 
       do j = 1, nth
           do k = 1, nti
             !
             castSPL(i,j,k) = 20.0_rp * log10(prms(i,j,k)/p0) 
             !
           end do 
       end do
    end do
  
  end subroutine castPrms2SPL
  !
  !
  !
  subroutine save_3D_matrix(fname,n_rows,n_cols,n_2d_mat,mat,ext)
       character(len=*), intent(in) :: fname
       integer, intent(in) :: n_rows ! number of rows
       integer, intent(in) :: n_cols ! number of columns
       integer, intent(in) :: n_2d_mat ! number of 2D submatrices
       real(kind=rp), intent(in), dimension(n_rows,n_cols,n_2d_mat) :: mat
       character(len=*), intent(in), optional :: ext
       !
       integer :: unit, imat, j, k
                     
       unit = get_unit_number()
       if( present(ext) ) then
         open(unit,file=fname//'choice3DMatrix'//ext,status='unknown')
         write(*,*) 'Creating file '//fname//'choice3DMatrix'//ext
       else
         write(*,*) 'Creating file '//fname//'choice3DMatrix.m'
         open(unit,file=fname//'choice3DMatrix.m',status='unknown')
       end if
       
       do imat = 1, n_2d_mat
         write(unit,'(A,$)') 'chMat(:,:,'//trim(adjustl(inumber_string_cast(imat)))//') = ['
         !
         do j = 1, n_rows
           do k = 1, n_cols
             write(unit,'(F20.10,$)') mat(j,k,imat)
           end do
           write(unit,*) ' '
         end do     
         write(unit,'(A)') '];'
       end do
       
       close(unit)
     
  end subroutine save_3D_matrix
  !
  !
  !
  subroutine save_matrix(fname,n_rows,n_cols,mat,ext)
       character(len=*), intent(in) :: fname
       integer, intent(in) :: n_rows ! number of rows
       integer, intent(in) :: n_cols ! number of columns
       real(kind=rp), intent(in), dimension(n_rows,n_cols) :: mat
       character(len=*), intent(in), optional :: ext
       !
       integer :: unit, i, j
              
       unit = get_unit_number()
       if( present(ext) ) then
         open(unit,file=fname//'choiceMatrix'//ext,status='unknown')
       else
         open(unit,file=fname//'choiceMatrix.m',status='unknown')
       end if
       
     
       write(unit,'(A,$)') 'chMat = ['
       
       do i = 1, n_rows
         do j = 1, n_cols
           write(unit,'(F20.10,$)') mat(i,j)
         end do
         write(unit,*) ' '
       end do
       
       write(unit,'(A)') '];'
       
       close(unit)
     
  end subroutine save_matrix
  !
  !
  !
  subroutine save_x_versus_y(fname,nx,ny,x,y,ext)
       character(len=*), intent(in) :: fname
       integer, intent(in) :: nx ! number of x
       integer, intent(in) :: ny ! number of y
       real(kind=rp), intent(in), dimension(nx) :: x
       real(kind=rp), intent(in), dimension(ny) :: y
       character(len=*), intent(in), optional :: ext
       !
       integer :: unit, i, j
              
       unit = get_unit_number()
       if( present(ext) ) then
         open(unit,file=fname//ext,status='unknown')
       else ! default is MATLAB
         open(unit,file=fname//'.m',status='unknown')
       end if
            
       write(unit,'(A)') 'clear all ; close all'
       write(unit,'(A,$)') 'x = ['
       
       do i = 1, nx
          write(unit,'(F20.10,$)') x(i)
       end do
       write(unit,'(A)') '];'
       
       write(unit,'(A,$)') 'y = ['
       
       do i = 1, ny
          write(unit,'(F20.10,$)') y(i)
       end do
       write(unit,'(A)') '];'
       
       write(unit,'(A)') 'plot(x,y)'
       
       close(unit)
     
  end subroutine save_x_versus_y
  !
  !
  !
  function get_p_ambient(altitude)
    real(kind=rp), intent(in) :: altitude
    real(kind=rp) :: get_p_ambient
  
    if ( altitude .lt. 11000.0_rp  ) then
       get_p_ambient = 101325.0_rp*(one-2.2557E-5_rp*altitude)**5.2561_rp
    elseif ( altitude.ge.11000.0_rp .and. altitude.lt. 25000.0_rp) then
       get_p_ambient = 22632.0_rp/exp(1.5769e-4_rp*(altitude-11000.0_rp))
    elseif ( altitude.gt. 25000.0_rp ) then
       get_p_ambient = 22632.0_rp/exp(1.5769e-4_rp*(altitude-11000.0_rp))
    endif

  end function get_p_ambient
  !
  !
  !
  function degree(rad)
    real(kind=rp), intent(in) :: rad
    real(kind=rp) :: degree
    
    degree = rad*180.0_rp/pi_num
  
  end function degree
  !
  ! get_n_tokens_in_string: 
  !
  !   counts the number of substrings in a long string. 
  !   tokenizes based on whitespace
  !
  function get_n_tokens_in_string(str)
    character(len=*), intent(in) :: str
    integer :: get_n_tokens_in_string
    !
    ! locals
    integer :: i, j, k, isave, lwidth
    !
    lwidth = len(str)
    !
    get_n_tokens_in_string = 0
    !
    i = 1
    outer: do 
       if(part_of_a_token(str(i:i))) then
          k = 1
          isave = i
          inner: do j = i,lwidth
             if(white_space(str(j:j))) then
                !
                i = i + 1
                !
                get_n_tokens_in_string = get_n_tokens_in_string + 1
                !
                if(i.ge.lwidth) exit outer
                exit inner
             else
                i = i + 1	
                k = k + 1	     
                if(i.ge.lwidth) exit outer
             endif
          end do inner
       else
          i = i + 1
       endif
       !
       if(i.ge.lwidth) exit outer
    end do outer
    !
    if ( part_of_a_token(str(lwidth:lwidth)) ) then 
       ! no trailing blank => last digit has not been written yet
       !
       get_n_tokens_in_string = get_n_tokens_in_string + 1
       !
    end if
    !
  end function get_n_tokens_in_string
  !  
  ! tokenize_string: 
  !
  ! tokenizes a character string into separate substrings. 
  ! splitting based on whitespace.
  !
  function tokenize_string(str,n_ss)
    character(len=*), intent(in) :: str
    integer, intent(in) :: n_ss ! number of substrings
    !
    integer, parameter :: max_len = 32
    character(len=max_len), dimension(n_ss) :: tokenize_string
    !
    ! locals
    integer :: i, j, k, isave, lwidth, counter
    !
    lwidth = len(str)
    !
    i = 1
    counter = 1
    !
    outer: do 
       if(part_of_a_token(str(i:i))) then
          k = 1
          isave = i
          inner: do j = i,lwidth
             if(white_space(str(j:j))) then
                !
				tokenize_string(counter) = ' '
                tokenize_string(counter) = str(isave:isave+k-1)  
                !
                if( k+1 .gt. max_len ) then 
				  call error_proc(1000,& 
                     & 'tokenize_string too short. max. no. characters is' & 
                     & //inumber_string_cast(max_len),'tokenize_string in service_module')
				end if
                !
                i = i + 1
                !
                if(i.ge.lwidth) exit outer
				counter = counter + 1
				!
                exit inner
             else
                i = i + 1	
                k = k + 1	     
                if(i.ge.lwidth) exit outer
             endif
          end do inner
       else
          i = i + 1
       endif
       !
       if(i.ge.lwidth) exit outer
    end do outer
    !
    if ( part_of_a_token(str(lwidth:lwidth)) ) then 
       ! no trailing blank => last digit has not been written yet
       tokenize_string(counter) = str(isave:isave+k-1)  
    end if
	!
  end function tokenize_string
  !
  ! parseTable: use to develop the code (conversion of MATLAB-tables
  ! to internal gestpan parameter statements).
  !
  subroutine parseTable(fname)
    character(len=*), intent(in) :: fname
    !
    integer, parameter :: maxNoColumns = 7 ! largest number of columns any module has in its performance file.
    character(len=32), allocatable, dimension(:,:) :: storageMat
    integer :: i, slen, unit
    integer, parameter :: nrows = 36
    character(len=32), dimension(maxNoColumns) :: strVec
    
    allocate(storageMat(nrows,maxNoColumns))
    !
    call loadStorageMat2(fname,nrows,maxNoColumns,strVec,storageMat)

    unit = get_unit_number()
    open(unit,file='output.txt',status='unknown')

    do i = 1, maxNoColumns
       call writeColumn(nrows,unit,storageMat(:,i))
    end do
    
    close(unit)
        
  end subroutine parseTable  
  !
  !
  !
  subroutine writeColumn(nd,unit,svec)
    integer, intent(in) :: nd
    integer, intent(in) :: unit
    character(len=32), intent(in), dimension(nd) :: svec
    character(len=32) :: str
    integer :: i, j, ptr 
    
    ptr = 1
    do i = 1, 4
        write(unit,'(A,$)') '   & '
        do j = 1, 9
           str = trim(svec(ptr))//'_rp, '
           write(unit,'(A,$)') trim(str)
           ptr = ptr + 1
        end do
        write(unit,'(A)') '   & '
    end do
    
  
      
    
  end subroutine writeColumn
  !
  !
  !
  subroutine setMachNumbers(p1,t1,g1,Ain,D1,xnl,gamma,Mtip,Mu,Umid)
    real(kind=rp), intent(in) :: p1,t1,g1,Ain,D1,xnl,gamma
    real(kind=rp), intent(out) :: Mtip,Mu
    !
    real(kind=rp) :: rt, rh, R, Max, xfunc0, ts, Utip, Umid
    
    rt = D1 / two
    rh = sqrt(rt**2 - Ain/pi_num) 
    
    R = get_R(zero)

    xfunc0 = g1*sqrt(R*t1)/(p1*Ain)
    Max = get_M(gamma,xfunc0)  
    
    ts = t1/(one+half*(gamma-one)*Max**two)
       
    Utip = two*pi_num*rt*xnl
    Umid = two*pi_num*((rt+rh)/two)*xnl
    Mu = Utip/sqrt(gamma*R*ts) ! blade Mach number
    
    Mtip = sqrt(Mu**2 + Max**2) ! use Mach number triangle, assume zero swirl
        
  end subroutine setMachNumbers
  !
  !
  !
  function get_Cax(gam,p,t,w,A)
    real(kind=rp), intent(in) :: gam,p,t,w,A
    real(kind=rp) :: get_Cax 
    !
    real(kind=rp) :: R, xfunc0, Max, ts 
  
    R = get_R(zero)

    xfunc0 = w*sqrt(R*t)/(p*A)
    Max = get_M(gam,xfunc0)  
    
    ts = t/(one+half*(gam-one)*Max**two)
    get_Cax = Max*sqrt(gam*R*ts)
       
  end function get_Cax
  !
  !
  !
  function get_p_amb(altitude)
    real(kind=rp), intent(in) :: altitude
    real(kind=rp) :: get_p_amb
    
    if ( altitude .lt. 11000.0_rp  ) then
      get_p_amb = 101325.0_rp*(one-2.2557E-5_rp*altitude)**5.2561_rp
    elseif ( altitude.ge.11000.0_rp .and. altitude.lt. 25000.0_rp) then
      get_p_amb = 22632.0_rp/exp(1.5769e-4_rp*(altitude-11000.0_rp))
    elseif ( altitude.gt. 25000.0_rp ) then    
      get_p_amb = 22632.0_rp/exp(1.5769e-4_rp*(altitude-11000.0_rp))
    endif
    
  end function get_p_amb
  !
  !
  !
  subroutine chapter3(noEngines,MTOW,EPNL_lateral,EPNL_cutback,EPNL_approach)
    integer, intent(in) :: noEngines
    real(kind=rp) :: MTOW ! maximum take-off weight in tonnes
    real(kind=rp), intent(out) :: EPNL_lateral,EPNL_cutback,EPNL_approach
    
! lateral
    EPNL_lateral = get_EPNL_lateral(MTOW)

! cutback
    EPNL_cutback = get_EPNL_cutback(noEngines,MTOW)
    
! approach
    EPNL_approach = get_EPNL_approach(MTOW)
  
  end subroutine chapter3
  !
  !
  !
  function get_EPNL_approach(MTOW)
    real(kind=rp), intent(in) :: MTOW ! maximum take-off weight in tonnes
    real(kind=rp) :: get_EPNL_approach
  
    if( MTOW.lt.35.0_rp) then
       get_EPNL_approach = 86.03_rp + 7.75_rp*log10(35.0_rp) ! 97.9965
    elseif( MTOW.ge.35.0_rp .and. MTOW.lt.400.0_rp ) then
       get_EPNL_approach = 86.03_rp + 7.75_rp*log10(MTOW)
    else
       get_EPNL_approach = 86.03_rp + 7.75_rp*log10(400.0_rp) ! 104.995
    end if

  end function get_EPNL_approach
  !
  !
  !
  function get_EPNL_cutback(noEng,MTOW)
    integer, intent(in) :: noEng
    real(kind=rp), intent(in) :: MTOW ! maximum take-off weight in tonnes
    real(kind=rp) :: get_EPNL_cutback

    if( noEng.eq.2 ) then
       if( MTOW.lt.48.1_rp) then
         get_EPNL_cutback = 66.65_rp + 13.29_rp*log10(48.1_rp) ! 89.0057
       elseif( MTOW.ge.48.1_rp .and. MTOW.lt.385.0_rp ) then
         get_EPNL_cutback = 66.65_rp + 13.29_rp*log10(MTOW)
       else
         get_EPNL_cutback = 66.65_rp + 13.29_rp*log10(385.0_rp) ! 101.01077
       end if
    elseif( noEng.eq.3 ) then
       if( MTOW.lt.28.6_rp ) then
         get_EPNL_cutback = 69.65_rp + 13.29_rp*log10(28.6_rp) ! 89.0051
       elseif( MTOW.ge.48.1_rp .and. MTOW.lt.385.0_rp ) then
         get_EPNL_cutback = 69.65_rp + 13.29_rp*log10(MTOW)
       else
         get_EPNL_cutback = 69.65_rp + 13.29_rp*log10(385.0_rp) ! 104.01077
       end if
    elseif( noEng.eq.4 ) then
       if( MTOW.lt.20.2_rp ) then
         get_EPNL_cutback = 71.65_rp + 13.29_rp*log10(20.2_rp) ! 88.998
       elseif( MTOW.ge.48.1_rp .and. MTOW.lt.385.0_rp ) then
         get_EPNL_cutback = 71.65_rp + 13.29_rp*log10(MTOW)
       else
         get_EPNL_cutback = 71.65_rp + 13.29_rp*log10(385.0_rp) ! 106.01077
       end if
    end if
        

  end function get_EPNL_cutback
  !
  !
  !
  function get_EPNL_lateral(MTOW)
    real(kind=rp), intent(in) :: MTOW ! maximum take-off weight in tonnes
    real(kind=rp) :: get_EPNL_lateral
  
    if( MTOW.lt.35.0_rp) then
       get_EPNL_lateral = 80.57_rp + 8.51_rp*log10(35.0_rp) ! 93.71
    elseif( MTOW.ge.35.0_rp .and. MTOW.lt.400.0_rp ) then
       get_EPNL_lateral = 80.57_rp + 8.51_rp*log10(MTOW)
    else
       get_EPNL_lateral = 80.57_rp + 8.51_rp*log10(400.0_rp) ! 102.71
    end if

  end function get_EPNL_lateral
end module choice_aux
