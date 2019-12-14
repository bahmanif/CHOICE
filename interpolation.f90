module interpolation
  !
  ! Module programmer : Tomas Groenstedt. Code is based on the work of
  !                     P. Dierckx. See reference below. 
  !
  !            
  !        Dierckx p. : Curve and surface fitting with splines,
  !                     monographs on numerical analysis, oxford
  !                     university press, 1993. 
  !                     
  !                     Author address: 
  !
  !                     P Dierckx
  !                     dept. computer science, k.u. leuven
  !                     celestijnenlaan 200a, b-3001 heverlee, belgium.
  !                     e-mail : Paul.Dierckx@cs.kuleuven.ac.be
  !
  use sim_precision
  use units_and_constants
  use error_handling
  !
  implicit none
  !
  private
  !
  public :: splev
  public :: splev2
  public :: curfit
  public :: bispev
  public :: interpolate ! not used by GESTPAN 
  public :: linear_interp
  public :: get_spline_coeffs
  public :: set_silent_interp
  !
  logical, private :: silent = .false.
  !
contains
  subroutine set_silent_interp(sl)
    logical, intent(in) :: sl

	silent = sl

  end subroutine set_silent_interp
  !
  !
  !
  subroutine linear_interp(f,ik,x,y,zout,idump)
    !
    ! f     = input vector containing table to be interpolated in
    ! ik    = is the starting position of the relevant table in the
    !         vector f
    ! x     = input value of x-coordinate
    ! y     = input value of y-coordinate
    ! zout  = output value of z-coordinate
    ! idump = information flag (0 = table interpolation ok. 1 = 
    !         table interpolation outside table range)
    !
    ! local variables
    !
    ! xin   = local copy of x input (modified if outside table range)
    ! yin   = local copy of x input (modified if outside table range)
    ! ifstx = location of first x value (left point)
    ! ifsty = location of first y value
    ! ixint = location of lower xvalue for interpolation
    ! nx    = number of x values in table
    ! ny    = number of y values in table
    ! ilstx = location of last x value
    ! ilsty = location of last y value
    ! il    = left side index
    ! ir    = right side index
    !   
    real(kind=rp), intent(in) :: f(*)
    integer, intent(in) :: ik
    real(kind=rp), intent(in) :: x,y
    real(kind=rp), intent(out) :: zout
    integer, intent(out) :: idump
    !
    integer ifstx,nx,ilstx,ny,il,ir,m,indx,ixint,iaint,ifsty,ilsty&
         & ,iyint,i,j
    real(kind=rp) xin,yin,x1,x2,xb1,xb2,x12,y1,yb1,yb2,y12,y2,z,p(2),q(2)
    !
    xin = x   ! xin and yin might be changed below if
    yin = y   ! interpolation input shows to be outside table limit
    !
    idump=0   ! ok so far
    zout = zero
    !
    ifstx = ik + 2
    !
    ny = f (ifstx - 1) + .0001
    nx = f (ifstx - 2) + .0001
    !
    ilstx = (ifstx - 1) + nx
    !
    if (xin.le.f(ifstx)) then  ! low end of xrange
      if(xin.eq.f(ifstx)) then
        ixint = ifstx
      endif
      if(xin.lt.f(ifstx)) then
        ixint = ifstx
        !
        xin=f(ifstx)
        idump=1
      endif
    elseif (xin.ge.f(ilstx)) then  ! high end of xrange
      if(xin.eq.f(ilstx)) then
        ixint = ilstx - 1 
      end if
      if(xin.gt.f(ilstx)) then
        ixint = ilstx - 1 
        !
        xin=f(ilstx)
        idump=1
      endif
    else
      ! x within range of table. bracket the x-value
      !
      il = 1
      ir = nx
      !
      do 

        m = (il + ir) / 2
        indx = (ifstx - 1) + m
        if (xin-f(indx).eq.zero) then 
          ixint =  m + (ifstx - 1) - 1
        endif
        if (xin-f(indx).lt.zero) then
          ir = m
        else
          il = m
        endif
        if (ir-il-1.eq.0) exit
      enddo
      !
      ir = ir + (ifstx - 1)
      ixint = ir - 1
      if (ixint .lt. ifstx) ixint = ifstx
      !
    endif 
    !
    ! x1 = x-value of lower row
    ! x2 = x-value of higher row
    !
    x1 = f (ixint)
    x2 = f (ixint + 1 )
    xb1 = xin - x1
    xb2 = xin - x2
    x12 = x1  - x2
    !
    p (1) =  xb2 / x12   ! non-dimensional distances above
    p (2) = -xb1 / x12   ! non-dimensional distance from below    
    !
    if (ny.eq.0) then 
       iaint = ixint + nx
       ! interpolate
       zout = p (1) * f (iaint) + p (2) * f (iaint + 1)
       return 
    endif
    !
    !
    ifsty = ifstx + nx
    ilsty = ifsty + ny - 1
    !
    !--- find y interpolation coefficients
    !
    if (yin.le.f(ifsty)) then 
       if(yin.eq.f(ifsty)) then
          iyint = ifsty
       endif
       if(yin.lt.f(ifsty)) then
          iyint = ifsty
          !
          yin=f(ifsty)
          idump = 1
       endif
       elseif (yin.ge.f(ilsty)) then
       if(yin.eq.f(ilsty)) then
          iyint = ilsty - 1
       endif
       if(yin.gt.f(ilsty)) then
          iyint = ilsty - 1
          !
          yin=f(ilsty)
          idump= 1
       endif
    else
       il = 1
       ir = ny
       do 
          m = (il + ir) / 2
          indx = (ifsty -1) + m
          if  (yin.eq.f(indx)) then
             iyint = m + (ifsty -1) - 1
          endif
          if  (yin.le.f(indx)) then
             ir = m
          else
             il = m
          endif
          !
          if (ir-il-1.eq.0 ) exit
       enddo
       ir = ir + (ifsty -1)
       !
       iyint = ir - 1
       if (iyint.lt.ifsty) iyint = ifsty
    endif
    !
    y1 = f (iyint)
    y2 = f (iyint + 1 )
    yb1 = yin - y1
    yb2 = yin - y2
    y12 = y1  - y2
    !    
    q (1) =  yb2 / y12
    q (2) = -yb1 / y12
    !    
    iaint = iyint + ny * (ixint - ifstx) - 1
    do i=1,2
      iaint = iaint + ny
      z = zero
      do j=1,2
        indx = iaint + j
        z = z + q(j) * f(indx)
      enddo
      zout = zout + z * p (i)
    enddo
    return
  end subroutine linear_interp
  !
  !
  !
  subroutine interpolate(sm_lev,ik,tvect,t,ikn,tx,ikx,ty,iky,c,wrk&
       & ,lwrk,iwrk,mod_name,x,y,z)
    !
    ! sm_lev = degree of smoothness 
    ! tvect = vector containing the table
    ! ik = index pointing at relevant table in tab_vect
    ! tx = knots for spline interpolation (x-distribution)
    ! ikx = number of x-knots
    ! ty = knots for spline interpolation (y-distribution)
    ! iky =  number of y-knots
    ! c = spline coefficients
    ! iky =  number of spline coefficients
    ! wrk = work vector used by the spline routines
    ! lwrk = number of elements allocated for the work and iwork arrays
    ! iwork = integer vector of length lwrk
    ! mod_name = name of module which is calling
    ! x = x-value to be interpolated for
    ! y = y-value to be interpolated for (zero if 1d interpolation)
    !
    real(kind=rp), intent(in) :: tvect(*),t(*),tx(*),ty(*),c(*),x,y,wrk(*)
    integer, intent(in) :: sm_lev,ik,ikn,ikx,iky,lwrk,iwrk(*)
    real(kind=rp), intent(out) :: z
    !
    character(len=*) :: mod_name
    ! locals 
    integer :: idump,i
    !
    select case( sm_lev )
    case ( 0 )
       call linear_interp(tvect,ik,x,y,z,idump)
       if(idump.eq.1) call error_proc(200,& 
	    & 'linear interpolation out of range','call from '//trim(mod_name))
    case ( 1 )
       call spline_fit(1)
       if(idump.eq.1 .and. .not.silent) call error_proc(200,& 
	    & 'linear spline fit out of range','call from '//trim(mod_name))
    case ( 3 )
       call spline_fit(3)
       if(idump.eq.1 .and. .not.silent) call error_proc(200,& 
	    & 'cubic spline fit out of range',& 
	    & 'call from '//trim(mod_name))
    case ( 5 )
       call spline_fit(5)
       if(idump.eq.1 .and. .not.silent) call error_proc(200,'fifth order spline fit out of range',& 
	    & 'call from '//trim(mod_name))
    case ( 6 )
       ! pre stored cubic spline routine
    end select
    !
  contains
    subroutine spline_fit(degree)
      integer, intent(in) :: degree
      !
      if(int(tvect(ik+1)).eq.0) then
	 ! 1d splines
	 call splev(t,ikn,c,degree,x,z,wrk,idump)
      else
	 ! 2d splines
	 call bispev(tx,ikx,ty,iky,c,degree,degree,x,y,z,wrk,lwrk&
	      & ,iwrk,lwrk,idump)
      endif
      !
    end subroutine spline_fit
  end subroutine interpolate
  !
  !
  !
  subroutine get_spline_coeffs(degree,f,ik,w,s,c,ikc,t,tx,ty,ikn,iknx&
       & ,ikny,wrk,lwrk,iwrk)
    !
    ! * routine interfacing the spline package FITPACK written
    ! by Paul Dierckx (NETLIB code) and the internal gestpan
    ! table format. The routine should be called once to
    ! determine spline coefficients and later spline evaluations
    ! should be done by directly calling splev or bispev
    !
    ! degree = degree of spline fit
    ! f =      vector where functions of module are stored
    ! ik =     starting position in vector f where relevant table starts
    ! c =      vector containing spline coefficients
    ! ikc =    starting position in vector c where relevant spline 
    !          coefficients are going to be stored. Incremented
    !          with n-degree on call, giving the first position
    !          for the spline coefficients for the next table. 
    ! ikn =    starting position in vector c where relevant spline 
    !          coefficients are going to be stored. Incremented with
    !          a value n+1 on call, giving the first position for the
    !          knots of the next table.
    ! w =      vector containing weighting coeffs (in this implementation)
    !          no active use of weighting is carried out (all are weighted
    !          with a factor 1.0). This could easily be activated by
    !          storing weights in the w array mimicking the pointer
    !          algorithm of the f vector (x-position could be used to store
    !          weights and the rest should be initialized to -1.0 to 
    !          increase chance of discovering any programming errors)
    !          
    ! s =      smoothness parameter
    !
    ! nest =   overestimate of storage requirement for knots (a cowards
    !          approach of nest=nx+degree+1 is taken since the arrays
    !          are generally very small for gestpan. if large tables
    !          are going to be used for some analysis nest = m/2 should
    !          be tried
    ! lwrk =   size of f vector. has to be of at least size 
    !          (nx*(degree+1)+nest*(7+3*degree))
    ! wrk  =   work array used by curfit of at least size lwrk
    !      
    !
    ! local variables
    !
    ! iopt =   in this implementation iopt is always set zero (iopt=0)
    !          since the knots and the coefficients are determined
    !          during the init mission of the modules and thereafter
    !          they are only used in splev calls
    ! ier =    error flag
    ! nk =     number of knots
    ! nxk =    number of x-knots (2d problem) 
    ! nyk =    number of y-knots (2d problem) 
    !
       integer, intent(in) :: degree,ik,lwrk
    real(kind=rp), intent(in)  :: f(*),w(*),s
    real(kind=rp), intent(out) :: c(*),t(*),tx(*),ty(*),wrk(*)
    integer, intent(out) :: ikc,ikn,iknx,ikny,iwrk(*)
    ! local variables
    real(kind=rp) :: fp
    integer :: nk,nxk,nyk,nx,ny,fx,ftv
    integer, parameter :: iopt = 0
    integer :: ierr
    !    
    nx = int(f(ik))
    ny = int(f(ik+1))
    !
    if(ny.gt.0) then ! 2d problem
       !
       call regrid(iopt,nx,f(ik+2:ik+1+nx),ny,f(ik+2+nx:ik+1+nx&
            & +ny),f(ik+2+nx+ny:ik+1+nx+ny+nx*ny),f(ik+2),f(ik+1+nx)&
            & ,f(ik+2+nx),f(ik+1+nx+ny),degree,degree,s,nx+degree+1&
            & ,ny+degree+1,nxk,tx,nyk,ty,c,fp,wrk,lwrk,iwrk,lwrk,ierr)
       !
       write(*,*) 'ierr =',ierr,'fp =',fp,'s =',s
       ! increment positional information
       ikc =  nx*ny + 1
       iknx = nxk+1
       ikny = nyk+1
!       write(*,*) 'c constants'
!       do i = 1,ikc
!          write(*,*) i,c(i)
!       enddo
!       write(*,*) 'tx constants'
!       do i = 1,iknx - 1
!          write(*,*) i,tx(i)
!       enddo
!       write(*,*) 'ty constants'
!       do i = 1,ikny - 1
!          write(*,*) i,ty(i)
!       enddo
!       stop 
       !
    elseif(ny.eq.0) then ! 1d problem
       ! determine spline coefficients and knots
       ! nx                    = m 
       ! f(ik+2:ik+1+nx)       = x 
       ! f(ik+2+nx:ik+1+nx+ny) = y 
       ! f(ik+2) = xb
       ! f(ik+1+nx) = xe
       ! nest = nx+degree+1
       !
       call curfit(iopt,nx,f(ik+2:ik+1+nx),f(ik+2+nx:ik+1+nx+nx),w&
            & ,f(ik+2),f(ik+1+nx),degree,s,nx+degree+1,nk,t,c,fp,wrk&
            & ,lwrk,iwrk,ierr)

       ! increment positional information
       ikc = nk-degree
       ikn = nk+1
       !
       write(*,*) 'ierr =',ierr,'fp =',fp,'s =',s
       !
    endif
    !
  end subroutine get_spline_coeffs
  !
  !
  !
  subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
    !
    !  given the set of data points (x(i),y(i)) and the set of positive
    !  numbers w(i),i=1,2,...,m, subroutine curfit determines a smooth 
    !  spline approximation of degree k on the interval xb <= x <= xe.
    !  if iopt = -1 curfit calculates the weighted least-squares spline
    !  according to a given set of knots.
    !  if iopt>= 0 the number of knots of the spline s(x) and the position
    !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
    !  ness of s(x) is then achieved by minimalizing the discontinuity
    !  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,
    !  k+3,...,n-k-1. the amount of smoothness is determined by the 
    !  condition that f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, 
    !  with s a given non-negative constant, called the smoothing factor.
    !  the fit s(x) is given in the b-spline representation (b-spline coef-
    !  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
    !  subroutine splev.
    !
    !  calling sequence:
    !     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
    !    * lwrk,iwrk,ier)
    !
    !  parameters:
    !   iopt  : integer flag. on entry iopt must specify whether a weighted
    !           least-squares spline (iopt=-1) or a smoothing spline (iopt=
    !           0 or 1) must be determined. if iopt=0 the routine will 
    !           start with an initial set of knots t(i)=xb, t(i+k+1)=xe, 
    !           i=1,2,... k+1. if iopt=1 the routine will continue with 
    !           the knots found at the last call of the routine.
    !           attention: a call with iopt=1 must always be immediately
    !           preceded by another call with iopt=1 or iopt=0.
    !           unchanged on exit.
    !   m     : integer. on entry m must specify the number of data points.
    !           m > k. unchanged on exit.
    !   x     : real array of dimension at least (m). before entry, x(i)
    !           must be set to the i-th value of the independent variable 
    !           x, for i=1,2,...,m. these values must be supplied in 
    !           strictly ascending order. unchanged on exit.
    !   y     : real array of dimension at least (m). before entry, y(i)
    !           must be set to the i-th value of the dependent variable y,
    !           for i=1,2,...,m. unchanged on exit.
    !   w     : real array of dimension at least (m). before entry, w(i)
    !           must be set to the i-th value in the set of weights. the
    !           w(i) must be strictly positive. unchanged on exit.
    !           see also further comments.
    !   xb,xe : real values. on entry xb and xe must specify the boundaries
    !           of the approximation interval. xb<=x(1), xe>=x(m).
    !           unchanged on exit.
    !   k     : integer. on entry k must specify the degree of the spline.
    !           1<=k<=5. it is recommended to use cubic splines (k=3).
    !           the user is strongly dissuaded from choosing k even,
    !           together with a small s-value. unchanged on exit.
    !   s     : real.on entry (in case iopt>=0) s must specify the 
    !           smoothing factor. s >=0. unchanged on exit.
    !           for advice on the choice of s see further comments.
    !   nest  : integer. on entry nest must contain an over-estimate of the
    !           total number of knots of the spline returned, to indicate
    !           the storage space available to the routine. nest >=2*k+2.
    !           in most practical situation nest=m/2 will be sufficient.
    !           always large enough is  nest=m+k+1, the number of knots
    !           needed for interpolation (s=0). unchanged on exit.
    !   n     : integer.
    !           unless ier =10 (in case iopt >=0), n will contain the
    !           total number of knots of the spline approximation returned.
    !           if the computation mode iopt=1 is used this value of n
    !           should be left unchanged between subsequent calls.
    !           in case iopt=-1, the value of n must be specified on entry.
    !   t     : real array of dimension at least (nest).
    !           on succesful exit, this array will contain the knots of the
    !           spline,i.e. the position of the interior knots t(k+2),
    !           t(k+3)...,t(n-k-1) as well as the position of the 
    !           additional knots t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=
    !           t(n)=xe needed for the b-spline representation.
    !           if the computation mode iopt=1 is used, the values of t(1),
    !           t(2),...,t(n) should be left unchanged between subsequent
    !           calls. if the computation mode iopt=-1 is used, the values
    !           t(k+2),...,t(n-k-1) must be supplied by the user, before
    !           entry. see also the restrictions (ier=10).
    !   c     : real array of dimension at least (nest).
    !           on succesful exit, this array will contain the coefficients
    !           c(1),c(2),..,c(n-k-1) in the b-spline representation of 
    !           s(x)
    !   fp    : real. unless ier=10, fp contains the weighted sum of
    !           squared residuals of the spline approximation returned.
    !   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
    !           used as working space. if the computation mode iopt=1 is
    !           used, the values wrk(1),...,wrk(n) should be left unchanged
    !           between subsequent calls.
    !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
    !           the array wrk as declared in the calling (sub)program.lwrk
    !           must not be too small (see wrk). unchanged on exit.
    !   iwrk  : integer array of dimension at least (nest).
    !           used as working space. if the computation mode iopt=1 is
    !           used,the values iwrk(1),...,iwrk(n) should be left 
    !           unchanged between subsequent calls.
    !   ier   : integer. unless the routine detects an error, ier 
    !           contains a
    !           non-positive value on exit, i.e.
    !    ier=0  : normal return. the spline returned has a residual sum of
    !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
    !             ive tolerance set to 0.001 by the program.
    !    ier=-1 : normal return. the spline returned is an interpolating
    !             spline (fp=0).
    !    ier=-2 : normal return. the spline returned is the weighted least-
    !             squares polynomial of degree k. in this extreme case fp
    !             gives the upper bound fp0 for the smoothing factor s.
    !    ier=1  : error. the required storage space exceeds the available
    !             storage space, as specified by the parameter nest.
    !             probably causes : nest too small. if nest is already
    !             large (say nest > m/2), it may also indicate that s is
    !             too small
    !             the approximation returned is the weighted least-squares
    !             spline according to the knots t(1),t(2),...,t(n). 
    !             (n=nest)
    !             the parameter fp gives the corresponding weighted sum of
    !             squared residuals (fp>s).
    !    ier=2  : error. a theoretically impossible result was found during
    !             the iteration proces for finding a smoothing spline with
    !             fp = s. probably causes : s too small.
    !             there is an approximation returned but the corresponding
    !             weighted sum of squared residuals does not satisfy the
    !             condition abs(fp-s)/s < tol.
    !    ier=3  : error. the maximal number of iterations maxit (set to 20
    !             by the program) allowed for finding a smoothing spline
    !             with fp=s has been reached. probably causes : s too small
    !             there is an approximation returned but the corresponding
    !             weighted sum of squared residuals does not satisfy the
    !             condition abs(fp-s)/s < tol.
    !    ier=10 : error. on entry, the input data are controlled on 
    !             validity the following restrictions must be satisfied.
    !             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
    !             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
    !             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
    !                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
    !                       the schoenberg-whitney conditions, i.e. there
    !                       must be a subset of data points xx(j) such that
    !                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
    !             if iopt>=0: s>=0
    !                         if s=0 : nest >= m+k+1
    !             if one of these conditions is found to be violated,
    !             control is immediately repassed to the calling program. 
    !             in that case there is no approximation returned.
    !
    !  further comments:
    !   by means of the parameter s, the user can control the tradeoff
    !   between closeness of fit and smoothness of fit of the 
    !   approximation. if s is too large, the spline will be too
    !   smooth and signal will belost ; if s is too small the spline
    !   will pick up too much noise. in the extreme cases the program
    !   will return an interpolating spline if s=0 and the weighted
    !   least-squares polynomial of degree k if s is very large.
    !   between these extremes, a properly chosen s will result
    !   in a good compromise between closeness of fit and smoothness
    !   of fit. to decide whether an approximation, corresponding to
    !   a certain s is satisfactory the user is highly recommended to
    !   inspect the fits graphically. recommended values for s depend
    !   on the weights w(i). if these are taken as 1/d(i) with d(i)
    !   an estimate of the standard deviation of y(i), a good s-value
    !   should be found in the range (m-sqrt(2*m),m+sqrt(2*m)). if
    !   nothing is known about the statistical error in y(i) each
    !   w(i) can be set equal to one and s determined by trial and
    !   error, taking account of the comments above. the best is then
    !   to start with a very large value of s ( to determine the
    !   least-squares polynomial and the corresponding upper bound
    !   fp0 for s) and then to progressively decrease the value of s
    !   (say by a factor 10 in the beginning, i.e. s=fp0/10, fp0/100
    !   ,...and more carefully as the approximation shows more detail)
    !   to obtain closer fits. to economize the search for a good s
    !   -value the program provides with different modes of
    !   computation. at the first call of the routine, or whenever he
    !   wants to restart with the initial set of knots the user must
    !   set iopt=0. if iopt=1 the program will continue with the set
    !   of knots found at the last call of the routine. this will save
    !   a lot of computation time if curfit is called repeatedly for
    !   different values of s. the number of knots of the spline
    !   returned and their location will depend on the value of s and
    !   on the complexity of the shape of the function underlying the
    !   data. but, if the computation mode iopt=1 is used, the knots
    !   returned may also depend on the s-values at previous calls (if
    !   these were smaller). therefore, if after a number of trials
    !   with different s-values and iopt=1, the user can finally
    !   accept a fit as satisfactory, it may be worthwhile for him to
    !   call curfit once more with the selected value for s but now
    !   with iopt=0. indeed, curfit may then return an approximation
    !   of the same quality of fit but with fewer knots and therefore
    !   better if data reduction is also an important objective for
    !   the user.
    !
    !  other subroutines required:
    !    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati
    !    ,fprota
    !
    !  references:
    !   dierckx p. : an algorithm for smoothing, differentiation and
    !                integration of experimental data using spline
    !                functions,j.comp.appl.maths 1 (1975) 165-184.
    !   dierckx p. : a fast algorithm for smoothing data on a
    !                rectangular grid while using spline functions,
    !                siam j.numer.anal. 19 (1982) 1286-1304.
    !   dierckx p. : an improved algorithm for curve fitting with
    !                spline functions, report tw54, dept. computer
    !                science,k.u. leuven, 1981.
    !   dierckx p. : curve and surface fitting with splines,
    !                monographs on numerical analysis, oxford
    !                university press, 1993.
    !
    !  author:
    !    p.dierckx
    !    dept. computer science, k.u. leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  creation date : may 1979
    !  latest update : march 1987
    !
    !  ..
    !  ..scalar arguments..
    real(kind=rp) xb,xe,s,fp
    integer iopt,m,k,nest,n,lwrk,ier
    !  ..array arguments..
    real(kind=rp) x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
    integer iwrk(nest)
    !  ..local scalars..
    real(kind=rp) tol
    integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,maxit,nmin
    !  ..
    !  we set up the parameters tol and maxit
    maxit = 20
    tol = 0.1e-03_rp
    ! before starting computations a data check is made. if the
    ! input data are invalid, control is immediately repassed to the
    ! calling program.
    ier = 10
    if(k.le.0 .or. k.gt.5) go to 50
    k1 = k+1
    k2 = k1+1
    if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
    nmin = 2*k1
    if(m.lt.k1 .or. nest.lt.nmin) go to 50
    lwest = m*k1+nest*(7+3*k)
    if(lwrk.lt.lwest) go to 50
    if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
    do i=2,m
       if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
    enddo
    !
    if(iopt.ge.0) go to 30
    if(n.lt.nmin .or. n.gt.nest) go to 50
    j = n
    do i=1,k1
       t(i) = xb
       t(j) = xe
       j = j-1
    enddo
    call fpchec(x,m,t,n,k,ier)
    if(ier) 50,40,50
30  if(s.lt.0.) go to 50
    if(s.eq.0. .and. nest.lt.(m+k1)) go to 50
    ier = 0
    ! we partition the working space and determine the spline
    ! approximation.
40  ifp = 1
    iz = ifp+nest
    ia = iz+nest
    ib = ia+nest*k1
    ig = ib+nest*k2
    iq = ig+nest*k2
    call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,&
         &wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
50  return
  end subroutine curfit
  !
  !
  !
  subroutine fpback(a,z,n,k,c,nest)
    !  subroutine fpback calculates the solution of the system of
    !  equations a*c = z with a a n x n upper triangular matrix
    !  of bandwidth k.
    !  ..
    !  ..scalar arguments..
    integer n,k,nest
    !  ..array arguments..
    real(kind=rp) a(nest,k),z(n),c(n)
    !  ..local scalars..
    real(kind=rp) store
    integer i,i1,j,k1,l,m
    !  ..
    k1 = k-1
    c(n) = z(n)/a(n,1)
    i = n-1
    if(i.eq.0) go to 30
    do j=2,n
       store = z(i)
       i1 = k1
       if(j.le.k1) i1 = j-1
       m = i
       do l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
       enddo
       c(i) = store/a(i,1)
       i = i-1
    enddo
30  return
  end subroutine fpback
  !
  !
  !
  subroutine fpbspl(t,n,k,x,l,h)
    !  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
    !  degree k at t(l) <= x < t(l+1) using the stable recurrence
    !  relation of de boor and cox.
    !  ..
    !  ..scalar arguments..
    real(kind=rp) x
    integer n,k,l
    !  ..array arguments..
    real(kind=rp) t(n),h(6)
    !  ..local scalars..
    real(kind=rp) f
    integer i,j,li,lj
    !  ..local arrays..
    real(kind=rp) hh(5)
    !  ..
    h(1) = one
    do 20 j=1,k
       do i=1,j
          hh(i) = h(i)
       enddo
       h(1) = zero
       do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
20  continue
    return
  end subroutine fpbspl
  !
  !
  !
  subroutine fpchec(x,m,t,n,k,ier)
    !  subroutine fpchec verifies the number and the position of the
    !  knots t(j),j=1,2,...,n of a spline of degree k, in relation to
    !  the number and the position of the data points x(i),i=1,2,...
    !  ,m. if all of the following conditions are fulfilled, the
    !  error parameter ier is set to zero. if one of the conditions
    !  is violated ier is set to ten.
    !      1) k+1 <= n-k-1 <= m
    !      2) t(1) <= t(2) <= ... <= t(k+1)
    !         t(n-k) <= t(n-k+1) <= ... <= t(n)
    !      3) t(k+1) < t(k+2) < ... < t(n-k)
    !      4) t(k+1) <= x(i) <= t(n-k)
    !      5) the conditions specified by schoenberg and whitney must
    !         hold for at least one subset of data points, i.e. there
    !         must be a subset of data points y(j) such that
    !         t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
    !  ..
    !  ..scalar arguments..
    integer m,n,k,ier
    !  ..array arguments..
    real(kind=rp) x(m),t(n)
    !  ..local scalars..
    integer i,j,k1,k2,l,nk1,nk2,nk3
    real(kind=rp) tj,tl
    !  ..
    k1 = k+1
    k2 = k1+1
    nk1 = n-k1
    nk2 = nk1+1
    ier = 10
    !  check condition no 1
    if(nk1.lt.k1 .or. nk1.gt.m) go to 80
    !  check condition no 2
    j = n
    do i=1,k
       if(t(i).gt.t(i+1)) go to 80
       if(t(j).lt.t(j-1)) go to 80
       j = j-1
    enddo
    !  check condition no 3
    do i=k2,nk2
       if(t(i).le.t(i-1)) go to 80
    enddo
    !  check condition no 4
    if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
    !  check condition no 5
    if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
    i = 1
    l = k2
    nk3 = nk1-1
    if(nk3.lt.2) go to 70
    do 60 j=2,nk3
       tj = t(j)
       l = l+1
       tl = t(l)
40     i = i+1
       if(i.ge.m) go to 80
       if(x(i).le.tj) go to 40
       if(x(i).ge.tl) go to 80
60  continue
70  ier = 0
80  return
  end subroutine fpchec
  !
  !
  !
  subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,&
       & n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
    !  ..
    !  ..scalar arguments..
    real(kind=rp) xb,xe,s,tol,fp
    integer iopt,m,k,nest,maxit,k1,k2,n,ier
    !  ..array arguments..
    real(kind=rp) x(m),y(m),w(m),t(nest),c(nest),fpint(nest),&
         & z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
    integer nrdata(nest)
    !  ..local scalars..
    real(kind=rp) acc,con1,con4,con9,cos,half,fpart,fpms,fpold,fp0,f1,f2,f3,&
         & one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
    integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0,&
         & mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
    !  ..local arrays..
    real(kind=rp) h(7)
    !  ..function references
    real(kind=rp) abs
   integer max0,min0
    !  ..subroutine references..
    !    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
    !  ..
    !  set constants
    one = 0.1e+01_rp
    con1 = 0.1e0_rp
    con9 = 0.9e0_rp
    con4 = 0.4e-01_rp
    half = 0.5e0_rp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  part 1: determination of the number of knots and their position    !
    !  **************************************************************     !
    !  given a set of knots we compute the least-squares spline sinf(x),  !
    !  and the corresponding sum of squared residuals fp=f(p=inf).        !
    !  if iopt=-1 sinf(x) is the requested approximation.                 !
    !  if iopt=0 or iopt=1 we check whether we can accept the knots:      !
    !    if fp <=s we will continue with the current set of knots.        !
    !    if fp > s we will increase the number of knots and compute the   !
    !       corresponding least-squares spline until finally fp<=s.       !
    !    the initial choice of knots depends on the value of s and iopt.  !
    !    if s=0 we have spline interpolation; in that case the number of  !
    !    knots equals nmax = m+k+1.                                       !
    !    if s > 0 and                                                     !
    !      iopt=0 we first compute the least-squares polynomial of        !
    !      degree k; n = nmin = 2*k+2                                     !
    !      iopt=1 we start with the set of knots found at the last        !
    !      call of the routine, except for the case that s > fp0; then    !
    !      we compute directly the least-squares polynomial of degree k.  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  determine nmin, the number of knots for polynomial approximation.
    nmin = 2*k1
    if(iopt.lt.0) go to 60
    !  calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol*s
    !  determine nmax, the number of knots for spline interpolation.
    nmax = m+k1
    if(s.gt.0.) go to 45
    !  if s=0, s(x) is an interpolating spline.
    !  test whether the required storage space exceeds the available one.
    n = nmax
    if(nmax.gt.nest) go to 420
    !  find the position of the interior knots in case of interpolation.
10  mk1 = m-k1
    if(mk1.eq.0) go to 60
    k3 = k/2
    i = k2
    j = k3+2
    if(k3*2.eq.k) go to 30
    do l=1,mk1
       t(i) = x(j)
       i = i+1
       j = j+1
    enddo
    go to 60
30  do l=1,mk1
       t(i) = (x(j)+x(j-1))*half
       i = i+1
       j = j+1
    enddo
    go to 60
    !  if s>0 our initial choice of knots depends on the value of iopt.
    !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
    !  polynomial of degree k which is a spline without interior knots.
    !  if iopt=1 and fp0>s we start computing the least squares spline
    !  according to the set of knots found at the last call of the routine.
45  if(iopt.eq.0) go to 50
    if(n.eq.nmin) go to 50
    fp0 = fpint(n)
    fpold = fpint(n-1)
    nplus = nrdata(n)
    if(fp0.gt.s) go to 60
50  n = nmin
    fpold = 0.0_rp
    nplus = 0
    nrdata(1) = m-2
    !  main loop for the different sets of knots. m is a save upper bound
    !  for the number of trials.
60  do 200 iter = 1,m
       if(n.eq.nmin) ier = -2
       !  find nrint, tne number of knot intervals.
       nrint = n-nmin+1
       !  find the position of the additional knots which are needed for
       !  the b-spline representation of s(x).
       nk1 = n-k1
       i = n
       do j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
       enddo
       !  compute the b-spline coefficients of the least-squares spline
       !  sinf(x). the observation matrix a is built up row by row and
       !  reduced to upper triangular form by givens transformations.
       !  at the same time fp=f(p=inf) is computed.
       fp = 0.0_rp
       !  initialize the observation matrix a.
       do i=1,nk1
          z(i) = 0.0_rp
          do j=1,k1
             a(i,j) = 0.0_rp
          enddo
       enddo
       l = k1
       do it=1,m
          !  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
          !  search for knot interval t(l) <= xi < t(l+1).
85        if(xi.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
          !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
90        call fpbspl(t,n,k,xi,l,h)
          do i=1,k1
             q(it,i) = h(i)
             h(i) = h(i)*wi
          enddo
          !  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
             j = j+1
             piv = h(i)
             if(piv.eq.0.) go to 110
             !  calculate the parameters of the givens transformation.
             call fpgivs(piv,a(j,1),cos,sin)
             !  transformations to right hand side.
             call fprota(cos,sin,yi,z(j))
             if(i.eq.k1) go to 120
             i2 = 1
             i3 = i+1
             do i1 = i3,k1
                i2 = i2+1
                !  transformations to left hand side.
                call fprota(cos,sin,h(i1),a(j,i2))
             enddo
110       continue
          !  add contribution of this row to the sum of squares of residual
          !  right hand sides.
120       fp = fp+yi**2
       enddo
       if(ier.eq.(-2)) fp0 = fp
       fpint(n) = fp0
       fpint(n-1) = fpold
       nrdata(n) = nplus
       !  backward substitution to obtain the b-spline coefficients.
       call fpback(a,z,nk1,k1,c,nest)
       !  test whether the approximation sinf(x) is an acceptable solution.
       if(iopt.lt.0) go to 440
       fpms = fp-s
       if(abs(fpms).lt.acc) go to 440
       !  if f(p=inf) < s accept the choice of knots.
       if(fpms.lt.0.) go to 250
       !  if n = nmax, sinf(x) is an interpolating spline.
       if(n.eq.nmax) go to 430
       !  increase the number of knots.
       !  if n=nest we cannot increase the number of knots because of
       !  the storage capacity limitation.
       if(n.eq.nest) go to 420
       !  determine the number of knots nplus we are going to add.
       if(ier.eq.0) go to 140
       nplus = 1
       ier = 0
       go to 150
140    npl1 = nplus*2
       rn = nplus
       if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
       nplus = min0(nplus*2,max0(npl1,nplus/2,1))
150    fpold = fp
       !  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
       !  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
       fpart = 0.0_rp
       i = 1
       l = k2
       new = 0
       do 180 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
160       term = 0.0_rp
          l0 = l-k2
          do j=1,k1
             l0 = l0+1
             term = term+c(l0)*q(it,j)
          enddo
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
180    continue
       fpint(nrint) = fpart
       do l=1,nplus
          !  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
          !  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 10
          !  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
       enddo
       !  restart the computations with the new set of knots.
200 continue
    !  test whether the least-squares kth degree polynomial is a solution
    !  of our approximation problem.
250 if(ier.eq.(-2)) go to 440
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  part 2: determination of the smoothing spline sp(x).                !
    !  ***************************************************                 !
    !  we have determined the number of knots and their position.          !
    !  we now compute the b-spline coefficients of the smoothing spline    !
    !  sp(x). the observation matrix a is extended by the rows of matrix   !
    !  b expressing that the kth derivative discontinuities of sp(x) at    !
    !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     !
    !  ponding weights of these additional rows are set to 1/p.            !
    !  iteratively we then have to determine the value of p such that      !
    !  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    !
    !  the least-squares kth degree polynomial corresponds to p=0, and     !
    !  that the least-squares spline corresponds to p=infinity. the        !
    !  iteration process which is proposed here, makes use of rational     !
    !  interpolation. since f(p) is a convex and strictly decreasing       !
    !  function of p, it can be approximated by a rational function        !
    !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  !
    !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      !
    !  to calculate the new value of p such that r(p)=s. convergence is    !
    !  guaranteed by taking f1>0 and f3<0.                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  evaluate the discontinuity jump of the kth derivative of the
    !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
       call fpdisc(t,n,k2,b,nest)
    !  initial value for p.
    p1 = 0.0_rp
    f1 = fp0-s
    p3 = -one
    f3 = fpms
    p = 0.0_rp
    do i=1,nk1
       p = p+a(i,1)
    enddo
    rn = nk1
    p = rn/p
    ich1 = 0
    ich3 = 0
    n8 = n-nmin
    !  iteration process to find the root of f(p) = s.
    do 360 iter=1,maxit
       !  the rows of matrix b with weight 1/p are rotated into the
       !  triangularised observation matrix a which is stored in g.
       pinv = one/p
       do i=1,nk1
          c(i) = z(i)
          g(i,k2) = 0.0_rp
          do j=1,k1
             g(i,j) = a(i,j)
          enddo
       enddo
       do 300 it=1,n8
          !  the row of matrix b is rotated into triangle by givens transformation
          do i=1,k2
             h(i) = b(it,i)*pinv
          enddo
          yi = 0.0_rp
          do j=it,nk1
             piv = h(1)
             !  calculate the parameters of the givens transformation.
             call fpgivs(piv,g(j,1),cos,sin)
             !  transformations to right hand side.
             call fprota(cos,sin,yi,c(j))
             if(j.eq.nk1) go to 300
             i2 = k1
             if(j.gt.n8) i2 = nk1-j
             do i=1,i2
                !  transformations to left hand side.
                i1 = i+1
                call fprota(cos,sin,h(i1),g(j,i1))
                h(i) = h(i1)
             enddo
             h(i2+1) = 0.0_rp
          enddo
300    continue
       !  backward substitution to obtain the b-spline coefficients.
       call fpback(g,c,nk1,k2,c,nest)
       !  computation of f(p).
       fp = 0.0_rp
       l = k2
       do it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.0_rp
          do j=1,k1
             l0 = l0+1
             term = term+c(l0)*q(it,j)
          enddo
          fp = fp+(w(it)*(term-y(it)))**2
       enddo
       !  test whether the approximation sp(x) is an acceptable solution.
       fpms = fp-s
       if(abs(fpms).lt.acc) go to 440
       !  test whether the maximal number of iterations is reached.
       if(iter.eq.maxit) go to 400
       !  carry out one more step of the iteration process.
       p2 = p
       f2 = fpms
       if(ich3.ne.0) go to 340
       if((f2-f3).gt.acc) go to 335
       !  our initial choice of p is too large.
       p3 = p2
       f3 = f2
       p = p*con4
       if(p.le.p1) p=p1*con9 + p2*con1
       go to 360
335    if(f2.lt.0.) ich3=1
340    if(ich1.ne.0) go to 350
       if((f1-f2).gt.acc) go to 345
       !  our initial choice of p is too small
       p1 = p2
       f1 = f2
       p = p/con4
       if(p3.lt.0.0_rp) go to 360
       if(p.ge.p3) p = p2*con1 + p3*con9
       go to 360
345    if(f2.gt.0.0_rp) ich1=1
       !  test whether the iteration process proceeds as theoretically
       !  expected.
350    if(f2.ge.f1 .or. f2.le.f3) go to 410
       !  find the new value for p.
       p = fprati(p1,f1,p2,f2,p3,f3)
360 continue
    !  error codes and messages.
400 ier = 3
    go to 440
410 ier = 2
    go to 440
420 ier = 1
    go to 440
430 ier = -1
440 return
  end subroutine fpcurf
  !
  !
  !
  subroutine fpdisc(t,n,k2,b,nest)
    !  subroutine fpdisc calculates the discontinuity jumps of the kth
    !  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
    !  ..scalar arguments..
    integer n,k2,nest
    !  ..array arguments..
    real(kind=rp) t(n),b(nest,k2)
    !  ..local scalars..
    real(kind=rp) an,fac,prod
    integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
    !  ..local array..
    real(kind=rp) h(12)
    !  ..
    k1 = k2-1
    k = k1-1
    nk1 = n-k1
    nrint = nk1-k
    an = nrint
    fac = an/(t(nk1+1)-t(k1))
    do l=k2,nk1
       lmk = l-k1
       do j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
       enddo
       lp = lmk
       do j=1,k2
          jk = j
          prod = h(j)
          do i=1,k
             jk = jk+1
             prod = prod*h(jk)*fac
          enddo
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
       enddo
    enddo
    return
  end subroutine fpdisc
  !
  !
  !
  subroutine fpgivs(piv,ww,cos,sin)
    !  subroutine fpgivs calculates the parameters of a givens
    !  transformation .
    !  ..
    !  ..scalar arguments..
    real(kind=rp) piv,ww,cos,sin
    !  ..local scalars..
    real(kind=rp) dd,one,store
    !  ..function references..
    real(kind=rp) abs,sqrt
    !  ..
    one = 0.1e+01_rp
    store = abs(piv)
    if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
    if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
    cos = ww/dd
    sin = piv/dd
    ww = dd
    return
  end subroutine fpgivs
  !
  !
  !
  subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
    !  subroutine fpknot locates an additional knot for a spline of degree
    !  k and adjusts the corresponding parameters,i.e.
    !    t     : the position of the knots.
    !    n     : the number of knots.
    !    nrint : the number of knotintervals.
    !    fpint : the sum of squares of residual right hand sides
    !            for each knot interval.
    !    nrdata: the number of data points inside each knot interval.
    !  istart indicates that the smallest data point at which the new knot
    !  may be added is x(istart+1)
    !  ..
    !  ..scalar arguments..
    integer m,n,nrint,nest,istart
    !  ..array arguments..
    real(kind=rp) x(m),t(nest),fpint(nest)
    integer nrdata(nest)
    !  ..local scalars..
    real(kind=rp) an,am,fpmax
    integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,next,nrx,number
    !  ..
    k = (n-nrint-1)/2
    !  search for knot interval t(number+k) <= x <= t(number+k+1) where
    !  fpint(number) is maximal on the condition that nrdata(number)
    !  not equals zero.
    fpmax = zero
    jbegin = istart
    do j=1,nrint
       jpoint = nrdata(j)
       if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
       fpmax = fpint(j)
       number = j
       maxpt = jpoint
       maxbeg = jbegin
10     jbegin = jbegin+jpoint+1
    enddo
    !  let coincide the new knot t(number+k+1) with a data point x(nrx)
    !  inside the old knot interval t(number+k) <= x <= t(number+k+1).
    ihalf = maxpt/2+1
    nrx = maxbeg+ihalf
    next = number+1
    if(next.gt.nrint) go to 40
    !  adjust the different parameters.
    do j=next,nrint
       jj = next+nrint-j
       fpint(jj+1) = fpint(jj)
       nrdata(jj+1) = nrdata(jj)
       jk = jj+k
       t(jk+1) = t(jk)
    enddo
40  nrdata(number) = ihalf-1
    nrdata(next) = maxpt-ihalf
    am = maxpt
    an = nrdata(number)
    fpint(number) = fpmax*an/am
    an = nrdata(next)
    fpint(next) = fpmax*an/am
    jk = next+k
    t(jk) = x(nrx)
    n = n+1
    nrint = nrint+1
    return
  end subroutine fpknot
  !
  !
  !
  function fprati(p1,f1,p2,f2,p3,f3)
    !
    !  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
    !  gives the value of p such that the rational interpolating function
    !  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
    !  ..
    !  ..scalar arguments..
    real(kind=rp) p1,f1,p2,f2,p3,f3,fprati
    !  ..local scalars..
    real(kind=rp) h1,h2,h3,p
    !  ..
    if(p3.gt.0.0_rp) go to 10
    !  value of p in case p3 = infinity.
    p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
    go to 20
    !  value of p in case p3 ^= infinity.
10  h1 = f1*(f2-f3)
    h2 = f2*(f3-f1)
    h3 = f3*(f1-f2)
    p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
    !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
20  if(f2.lt.0.) go to 30
    p1 = p2
    f1 = f2
    go to 40
30  p3 = p2
    f3 = f2
40  fprati = p
    return
  end function fprati
  !
  !
  !
  subroutine fprota(cos,sin,a,b)
    !  subroutine fprota applies a givens rotation to a and b.
    !  ..
    !  ..scalar arguments..
    real(kind=rp) cos,sin,a,b
    ! ..local scalars..
    real(kind=rp) stor1,stor2
    !  ..
    stor1 = a
    stor2 = b
    b = cos*stor2+sin*stor1
    a = cos*stor1-sin*stor2
    !
    return
  end subroutine fprota
  !
  !
  !
  subroutine splev(t,n,c,k,x,y,wrk,ier)
    !
    ! the subroutine splev evaluates a spline s(x) of degree k,
    ! given in its b-spline representation.
    !
    ! tomas groenstedt modified/adapted 1998_08_28:
    !
    !  * the routine to extrapolate points using the assumption of 
    !  constant derivative 
    !
    !  * the routine only to work with scalar input 
    !
    !  calling sequence:
    !     call splev(t,n,c,k,x,y,ier)
    !
    !  input parameters:
    !    t    : array,length n, which contains the position of the knots.
    !    n    : integer, giving the total number of knots of s(x).
    !    c    : array,length n, which contains the b-spline coefficients.
    !    k    : integer, giving the degree of s(x).
    !    x    : scalar, , which contains the points where s(x) must
    !           be evaluated.
    !
    !  output parameter:
    !    y    : scalar giving the value of s(x) at the different
    !           points.
    !    ier  : error flag
    !    ier =   0 : normal return
    !    ier = 1,2 : extrapolation occurred
    !    ier =  10 : invalid input data (see restrictions)
    !
    !  restrictions:
    !    m >= 1
    !    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
    !
    !  other subroutines required: fpbspl.
    !
    !  references :
    !    de boor c  : on calculating with b-splines, j. approximation theory
    !                 6 (1972) 50-62.
    !    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
    !                 applics 10 (1972) 134-149.
    !    dierckx p. : curve and surface fitting with splines, monographs on
    !                 numerical analysis, oxford university press, 1993.
    !
    !  author :
    !    p.dierckx
    !    dept. computer science, k.u.leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  latest update : march 1987
    !
    !  ..scalar arguments..
       integer n,k,m,ier,interval_flag
    !  ..array arguments..
    real(kind=rp) t(n),c(n),x,y,wrk(n)
    !  ..local scalars..
    integer i,j,k1,l,ll,l1,nk1
    real(kind=rp) arg,arg1(1),der1(1),sp,tb,te
    !  ..local array..
    real(kind=rp) h(6)
    !  ..
    !  before starting computations a data check is made. if the input data
    !  are invalid control is immediately repassed to the calling program.
    ier = 10
    interval_flag = 0
    ier = 0
    !  fetch tb and te, the boundaries of the approximation interval.
    k1 = k+1
    nk1 = n-k1
    tb = t(k1)
    te = t(nk1+1)
    l = k1
    l1 = l+1
    !  main loop for the different points.
    !
    arg = x
    if(arg.lt.tb) then
       interval_flag = -1
       arg = tb
       arg1(1) = tb
       call splder(t,n,c,k,1,arg1,der1,1,wrk,ier)          
       ier = 1
    endif
    if(arg.gt.te) then
       interval_flag = +1
       arg = te
       arg1(1) = te
       call splder(t,n,c,k,1,arg1,der1,1,wrk,ier)          
       ier = 2
    endif
    !
    do  !  search for knot interval t(l) <= arg < t(l+1) 
       if(arg.lt.t(l1) .or. l.eq.nk1) exit
       l = l1
       l1 = l+1
    enddo
    !  evaluate the non-zero b-splines at arg.
    call fpbspl(t,n,k,arg,l,h)
    !  find the value of s(x) at x=arg.
    sp = zero
    ll = l-k1
    do j=1,k1
       ll = ll+1
       sp = sp+c(ll)*h(j)
    enddo
    if(interval_flag.eq.0) then
       y = sp
    elseif(interval_flag.eq.-1) then
       !
       y = sp - der1(1)*(tb-x)
    elseif(interval_flag.eq.1) then
       !          
       y = sp + der1(1)*(x-te)
    endif
100 return
  end subroutine splev
  !
  ! splev2: 
  !    same as splev. used to circumvent compiler bug.
  !
  subroutine splev2(t,n,c,k,x,y,wrk,ier)
    !
    ! the subroutine splev2 evaluates a spline s(x) of degree k,
    ! given in its b-spline representation.
    !
    ! tomas groenstedt modified/adapted 1998_08_28:
    !
    !  * the routine to extrapolate points using the assumption of 
    !  constant derivative 
    !
    !  * the routine only to work with scalar input 
    !
    !  calling sequence:
    !     call splev2(t,n,c,k,x,y,ier)
    !
    !  input parameters:
    !    t    : array,length n, which contains the position of the knots.
    !    n    : integer, giving the total number of knots of s(x).
    !    c    : array,length n, which contains the b-spline coefficients.
    !    k    : integer, giving the degree of s(x).
    !    x    : scalar, , which contains the points where s(x) must
    !           be evaluated.
    !
    !  output parameter:
    !    y    : scalar giving the value of s(x) at the different
    !           points.
    !    ier  : error flag
    !    ier =   0 : normal return
    !    ier = 1,2 : extrapolation occurred
    !    ier =  10 : invalid input data (see restrictions)
    !
    !  restrictions:
    !    m >= 1
    !    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
    !
    !  other subroutines required: fpbspl.
    !
    !  references :
    !    de boor c  : on calculating with b-splines, j. approximation theory
    !                 6 (1972) 50-62.
    !    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
    !                 applics 10 (1972) 134-149.
    !    dierckx p. : curve and surface fitting with splines, monographs on
    !                 numerical analysis, oxford university press, 1993.
    !
    !  author :
    !    p.dierckx
    !    dept. computer science, k.u.leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  latest update : march 1987
    !
    !  ..scalar arguments..
       integer n,k,m,ier,interval_flag
    !  ..array arguments..
    real(kind=rp) t(n),c(n),x,y,wrk(n)
    !  ..local scalars..
    integer i,j,k1,l,ll,l1,nk1
    real(kind=rp) arg,arg1(1),der1(1),sp,tb,te
    !  ..local array..
    real(kind=rp) h(6)
    !  ..
    !  before starting computations a data check is made. if the input data
    !  are invalid control is immediately repassed to the calling program.
    ier = 10
    interval_flag = 0
    ier = 0
    !  fetch tb and te, the boundaries of the approximation interval.
    k1 = k+1
    nk1 = n-k1
    tb = t(k1)
    te = t(nk1+1)
    l = k1
    l1 = l+1
    !  main loop for the different points.
    !
    arg = x
    if(arg.lt.tb) then
       interval_flag = -1
       arg = tb
       arg1(1) = tb
       call splder(t,n,c,k,1,arg1,der1,1,wrk,ier)          
       ier = 1
    endif
    if(arg.gt.te) then
       interval_flag = +1
       arg = te
       arg1(1) = te
       call splder(t,n,c,k,1,arg1,der1,1,wrk,ier)          
       ier = 2
    endif
    !
    do  !  search for knot interval t(l) <= arg < t(l+1) 
       if(arg.lt.t(l1) .or. l.eq.nk1) exit
       l = l1
       l1 = l+1
    enddo
    !  evaluate the non-zero b-splines at arg.
    call fpbspl(t,n,k,arg,l,h)
    !  find the value of s(x) at x=arg.
    sp = zero
    ll = l-k1
    do j=1,k1
       ll = ll+1
       sp = sp+c(ll)*h(j)
    enddo
    if(interval_flag.eq.0) then
       y = sp
    elseif(interval_flag.eq.-1) then
       !
       y = sp - der1(1)*(tb-x)
    elseif(interval_flag.eq.1) then
       !          
       y = sp + der1(1)*(x-te)
    endif
100 return
  end subroutine splev2
  !
  !
  !
  subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,y,z,wrk,lwrk,iwrk,kwrk,ier)
    !  subroutine bispev evaluates in a point (x,y) a bivariate
    !  spline s(x,y) of degrees kx and ky, given in the b-spline
    !  representation.
    !
    !  calling sequence:
    !     call bispev(tx,nx,ty,ny,c,kx,ky,x,y,z,wrk,lwrk,iwrk,kwrk,ier)
    !
    !  input parameters:
    !   tx    : real array, length nx, which contains the position of the
    !           knots in the x-direction.
    !   nx    : integer, giving the total number of knots in the x-direction
    !   ty    : real array, length ny, which contains the position of the
    !           knots in the y-direction.
    !   ny    : integer, giving the total number of knots in the y-direction
    !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
    !           b-spline coefficients.
    !   kx,ky : integer values, giving the degrees of the spline.
    !   x     : real scalar
    !           before entry x must be set to the x co-ordinate
    !   y     : real scalar
    !           before entry y must be set to the y co-ordinate
    !   wrk   : real array of dimension lwrk. used as workspace.
    !   lwrk  : integer, specifying the dimension of wrk.
    !           lwrk >= kx+ky+2
    !   iwrk  : integer array of dimension kwrk. used as workspace.
    !   kwrk  : integer, specifying the dimension of iwrk. kwrk >= 2.
    !
    !  output parameters:
    !   z     : real scalar
    !           on succesful exit z contains the value of s(x,y)
    !           at the point x,y
    !   ier   : integer error flag
    !   ier=0 : normal return
    !   ier=10: invalid input data (see restrictions)
    !
    !  restrictions:
    !   lwrk>=kx+ky+2, kwrk>=2
    !   tx(kx+1) <= x  <= tx(nx-kx)
    !   ty(ky+1) <= y  <= ty(ny-ky)
    !
    !  other subroutines required:
    !    fpbisp,fpbspl
    !
    !  references :
    !    de boor c : on calculating with b-splines, j. approximation theory
    !                6 (1972) 50-62.
    !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
    !                applics 10 (1972) 134-149.
    !    dierckx p. : curve and surface fitting with splines, monographs on
    !                 numerical analysis, oxford university press, 1993.
    !
    !  author :
    !    p.dierckx
    !    dept. computer science, k.u.leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  latest update : march 1987
    !  modified by tomas groenstedt 1998 
    !  (vector evaluation removed (mx = 1, my = 1))
    !
    !  ..scalar arguments..
    integer nx,ny,kx,ky,lwrk,kwrk,ier
    !  ..array arguments..
    integer iwrk(kwrk)
    real(kind=rp) tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x,y,z,wrk(lwrk)
    !  ..local scalars..
    integer i,iw,lwest
    !  ..
    !  before starting computations a data check is made. if the input data
    !  are invalid control is immediately repassed to the calling program.
    ier = 10
    lwest = kx+ky+2
    if(lwrk.lt.lwest) go to 100
    if(kwrk.lt.2) go to 100
60  ier = 0
    iw = kx+2
    call fpbisp(tx,nx,ty,ny,c,kx,ky,x,y,z,wrk(1),wrk(iw),iwrk(1),iwrk(1+1))
    if(x.lt.tx(kx+1).or.x.gt.tx(nx-kx).or.y.lt.ty(ky+1).or.y.gt.ty(ny-ky)) then
       z = extrap_spline(x,y)
    endif
100 return
  contains
    function extrap_spline(x,y)
      real(kind=rp), intent(in) :: x,y
      real(kind=rp) :: extrap_spline
      ! locals       
      real(kind=rp) :: xmax,ymax,xmin,ymin,xi,yi,zi,der,dx,dy
      !
      xmin = tx(kx+1)
      xmax = tx(nx-kx)
      ymin = ty(ky+1)
      ymax = ty(ny-ky)
      !
      dx = 0
      dy = 0
      !
      if(x.gt.xmax) then
         xi = tx(nx-kx)-(tx(nx-kx)-tx(nx-kx-1))/1000.0_rp
         call fpbisp(tx,nx,ty,ny,c,kx,ky,xi,y,zi,wrk(1)&
              & ,wrk(iw),iwrk(1),iwrk(2))
         der = (z-zi)/(xmax-xi)
         dx = (x-xmax)*der
         ier = 1
      endif
      if(x.lt.xmin) then
         xi = tx(kx+1)+(tx(kx+2)-tx(kx+1))/1000.0_rp
         call fpbisp(tx,nx,ty,ny,c,kx,ky,xi,y,zi,wrk(1)&
              & ,wrk(iw),iwrk(1),iwrk(2))
         der = (z-zi)/(xmin-xi)
         dx = (x-xmin)*der
         ier = 1
      endif
      if(y.gt.ymax) then
         yi = ty(ny-ky)-(ty(ny-ky)-ty(ny-ky-1))/1000.0_rp
         call fpbisp(tx,nx,ty,ny,c,kx,ky,x,yi,zi,wrk(1)&
              & ,wrk(iw),iwrk(1),iwrk(2))
         der = (z-zi)/(ymax-yi)
         dy = (y-ymax)*der
         ier = 1
      endif
      if(y.lt.ymin) then
         yi = ty(ky+1)+(ty(ky+2)-ty(ky+1))/1000.0_rp
         call fpbisp(tx,nx,ty,ny,c,kx,ky,x,yi,zi,wrk(1)&
              & ,wrk(iw),iwrk(1),iwrk(2))
         der = (z-zi)/(ymin-yi)
         dy = (y-ymin)*der
         ier = 1
      endif
      extrap_spline = z + dx + dy
      !
    end function extrap_spline
  end subroutine bispev
  !
  !
  !
  subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,y,z,wx,wy,lx,ly)
    !  ..scalar arguments..
    integer nx,ny,kx,ky
    !  ..array arguments..
    integer lx,ly
    real(kind=rp) tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x,y,z,wx(kx+1),wy(ky+1)
    !  ..local scalars..
    integer i1,j,j1,kx1,ky1,l,l1,l2,m,nkx1,nky1
    real(kind=rp) arg,sp,tb,te
    !  ..local arrays..
    real(kind=rp) h(6)
    !  ..subroutine references..
    !    fpbspl
    !  ..
    kx1 = kx+1
    nkx1 = nx-kx1
    tb = tx(kx1)
    te = tx(nkx1+1)
    l = kx1
    l1 = l+1
    arg = x
    if(arg.lt.tb) arg = tb
    if(arg.gt.te) arg = te
10  if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
    l = l1
    l1 = l+1
    go to 10
20  call fpbspl(tx,nx,kx,arg,l,h)
    lx = l-kx1
    do j=1,kx1
       wx(j) = h(j)
    enddo
    ky1 = ky+1
    nky1 = ny-ky1
    tb = ty(ky1)
    te = ty(nky1+1)
    l = ky1
    l1 = l+1
    arg = y
    if(arg.lt.tb) arg = tb
    if(arg.gt.te) arg = te
50  if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
    l = l1
    l1 = l+1
    go to 50
60  call fpbspl(ty,ny,ky,arg,l,h)
    ly = l-ky1
    do j=1,ky1
       wy(j) = h(j)
    enddo
80  continue
    l = lx*nky1
    do i1=1,kx1
       h(i1) = wx(i1)
    enddo
    l1 = l+ly
    sp = zero
    do i1=1,kx1
       l2 = l1
       do j1=1,ky1
          l2 = l2+1
          sp = sp+c(l2)*h(i1)*wy(j1)
       enddo
       l1 = l1+nky1
    enddo
    z = sp
    return
  end subroutine fpbisp
  !
  !
  !
  subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest&
       & ,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
    ! given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
    ! i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
    ! iate spline approximation s(x,y) of degrees kx and ky on the rect-
    ! angle xb <= x <= xe, yb <= y <= ye.
    ! if iopt = -1 regrid calculates the least-squares spline according
    ! to a given set of knots.
    ! if iopt >= 0 the total numbers nx and ny of these knots and their
    ! position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
    ! ally by the routine. the smoothness of s(x,y) is then achieved by
    ! minimalizing the discontinuity jumps in the derivatives of s(x,y)
    ! across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
    ! the amounth of smoothness is determined by the condition that f(p) =
    ! sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
    ! constant, called the smoothing factor.
    ! the fit is given in the b-spline representation (b-spline coefficients
    ! c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
    ! uated by means of subroutine bispev.
    !
    ! calling sequence:
    !     call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
    !    *  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
    !
    ! parameters:
    !  iopt  : integer flag. on entry iopt must specify whether a least-
    !          squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
    !          must be determined.
    !          if iopt=0 the routine will start with an initial set of knots
    !          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
    !          1,...,ky+1. if iopt=1 the routine will continue with the set
    !          of knots found at the last call of the routine.
    !          attention: a call with iopt=1 must always be immediately pre-
    !                     ceded by another call with iopt=1 or iopt=0 and
    !                     s.ne.0.
    !          unchanged on exit.
    !  mx    : integer. on entry mx must specify the number of grid points
    !          along the x-axis. mx > kx . unchanged on exit.
    !  x     : real array of dimension at least (mx). before entry, x(i)
    !          must be set to the x-co-ordinate of the i-th grid point
    !          along the x-axis, for i=1,2,...,mx. these values must be
    !          supplied in strictly ascending order. unchanged on exit.
    !  my    : integer. on entry my must specify the number of grid points
    !          along the y-axis. my > ky . unchanged on exit.
    !  y     : real array of dimension at least (my). before entry, y(j)
    !          must be set to the y-co-ordinate of the j-th grid point
    !          along the y-axis, for j=1,2,...,my. these values must be
    !          supplied in strictly ascending order. unchanged on exit.
    !  z     : real array of dimension at least (mx*my).
    !          before entry, z(my*(i-1)+j) must be set to the data value at
    !          the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
    !          unchanged on exit.
    !  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
    !  yb,ye   aries of the rectangular approximation domain.
    !          xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
    !          unchanged on exit.
    !  kx,ky : integer values. on entry kx and ky must specify the degrees
    !          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
    !          (kx=ky=3) splines. unchanged on exit.
    !  s     : real. on entry (in case iopt>=0) s must specify the smoothing
    !          factor. s >=0. unchanged on exit.
    !          for advice on the choice of s see further comments
    !  nxest : integer. unchanged on exit.
    !  nyest : integer. unchanged on exit.
    !          on entry, nxest and nyest must specify an upper bound for the
    !          number of knots required in the x- and y-directions respect.
    !          these numbers will also determine the storage space needed by
    !          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
    !          in most practical situation nxest = mx/2, nyest=my/2, will
    !          be sufficient. always large enough are nxest=mx+kx+1, nyest=
    !          my+ky+1, the number of knots needed for interpolation (s=0).
    !          see also further comments.
    !  nx    : integer.
    !          unless ier=10 (in case iopt >=0), nx will contain the total
    !          number of knots with respect to the x-variable, of the spline
    !          approximation returned. if the computation mode iopt=1 is
    !          used, the value of nx should be left unchanged between sub-
    !          sequent calls.
    !          in case iopt=-1, the value of nx should be specified on entry
    !  tx    : real array of dimension nmax.
    !          on succesful exit, this array will contain the knots of the
    !          spline with respect to the x-variable, i.e. the position of
    !          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
    !          position of the additional knots tx(1)=...=tx(kx+1)=xb and
    !          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
    !          if the computation mode iopt=1 is used, the values of tx(1),
    !          ...,tx(nx) should be left unchanged between subsequent calls.
    !          if the computation mode iopt=-1 is used, the values tx(kx+2),
    !          ...tx(nx-kx-1) must be supplied by the user, before entry.
    !          see also the restrictions (ier=10).
    !  ny    : integer.
    !          unless ier=10 (in case iopt >=0), ny will contain the total
    !          number of knots with respect to the y-variable, of the spline
    !          approximation returned. if the computation mode iopt=1 is
    !          used, the value of ny should be left unchanged between sub-
    !          sequent calls.
    !          in case iopt=-1, the value of ny should be specified on entry
    !  ty    : real array of dimension nmax.
    !          on succesful exit, this array will contain the knots of the
    !          spline with respect to the y-variable, i.e. the position of
    !          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
    !          position of the additional knots ty(1)=...=ty(ky+1)=yb and
    !          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
    !          if the computation mode iopt=1 is used, the values of ty(1),
    !          ...,ty(ny) should be left unchanged between subsequent calls.
    !          if the computation mode iopt=-1 is used, the values ty(ky+2),
    !          ...ty(ny-ky-1) must be supplied by the user, before entry.
    !          see also the restrictions (ier=10).
    !  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
    !          on succesful exit, c contains the coefficients of the spline
    !          approximation s(x,y)
    !  fp    : real. unless ier=10, fp contains the sum of squared
    !          residuals of the spline approximation returned.
    !  wrk   : real array of dimension (lwrk). used as workspace.
    !          if the computation mode iopt=1 is used the values of wrk(1),
    !          ...,wrk(4) should be left unchanged between subsequent calls.
    !  lwrk  : integer. on entry lwrk must specify the actual dimension of
    !          the array wrk as declared in the calling (sub)program.
    !          lwrk must not be too small.
    !           lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
    !            my*(ky+1) +u
    !           where u is the larger of my and nxest.
    !  iwrk  : integer array of dimension (kwrk). used as workspace.
    !          if the computation mode iopt=1 is used the values of iwrk(1),
    !          ...,iwrk(3) should be left unchanged between subsequent calls
    !  kwrk  : integer. on entry kwrk must specify the actual dimension of
    !          the array iwrk as declared in the calling (sub)program.
    !          kwrk >= 3+mx+my+nxest+nyest.
    !  ier   : integer. unless the routine detects an error, ier contains a
    !          non-positive value on exit, i.e.
    !   ier=0  : normal return. the spline returned has a residual sum of
    !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
    !            ive tolerance set to 0.001 by the program.
    !   ier=-1 : normal return. the spline returned is an interpolating
    !            spline (fp=0).
    !   ier=-2 : normal return. the spline returned is the least-squares
    !            polynomial of degrees kx and ky. in this extreme case fp
    !            gives the upper bound for the smoothing factor s.
    !   ier=1  : error. the required storage space exceeds the available
    !            storage space, as specified by the parameters nxest and
    !            nyest.
    !            probably causes : nxest or nyest too small. if these param-
    !            eters are already large, it may also indicate that s is
    !            too small
    !            the approximation returned is the least-squares spline
    !            according to the current set of knots. the parameter fp
    !            gives the corresponding sum of squared residuals (fp>s).
    !   ier=2  : error. a theoretically impossible result was found during
    !            the iteration proces for finding a smoothing spline with
    !            fp = s. probably causes : s too small.
    !            there is an approximation returned but the corresponding
    !            sum of squared residuals does not satisfy the condition
    !            abs(fp-s)/s < tol.
    !   ier=3  : error. the maximal number of iterations maxit (set to 20
    !            by the program) allowed for finding a smoothing spline
    !            with fp=s has been reached. probably causes : s too small
    !            there is an approximation returned but the corresponding
    !            sum of squared residuals does not satisfy the condition
    !            abs(fp-s)/s < tol.
    !   ier=10 : error. on entry, the input data are controlled on validity
    !            the following restrictions must be satisfied.
    !            -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
    !            nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
    !            lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
    !             my*(ky+1) +max(my,nxest),
    !            xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
    !            if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
    !                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
    !                        2*ky+2<=ny<=min(nyest,my+ky+1)
    !                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
    !                    the schoenberg-whitney conditions, i.e. there must
    !                    be subset of grid co-ordinates xx(p) and yy(q) such
    !                    that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
    !                           ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
    !            if iopt>=0: s>=0
    !                        if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
    !            if one of these conditions is found to be violated,control
    !            is immediately repassed to the calling program. in that
    !            case there is no approximation returned.
    !
    ! further comments:
    !   regrid does not allow individual weighting of the data-values.
    !   so, if these were determined to widely different accuracies, then
    !   perhaps the general data set routine surfit should rather be used
    !   in spite of efficiency.
    !   by means of the parameter s, the user can control the tradeoff
    !   between closeness of fit and smoothness of fit of the approximation.
    !   if s is too large, the spline will be too smooth and signal will be
    !   lost ; if s is too small the spline will pick up too much noise. in
    !   the extreme cases the program will return an interpolating spline if
    !   s=0 and the least-squares polynomial (degrees kx,ky) if s is
    !   very large. between these extremes, a properly chosen s will result
    !   in a good compromise between closeness of fit and smoothness of fit.
    !   to decide whether an approximation, corresponding to a certain s is
    !   satisfactory the user is highly recommended to inspect the fits
    !   graphically.
    !   recommended values for s depend on the accuracy of the data values.
    !   if the user has an idea of the statistical errors on the data, he
    !   can also find a proper estimate for s. for, by assuming that, if he
    !   specifies the right s, regrid will return a spline s(x,y) which
    !   exactly reproduces the function underlying the data he can evaluate
    !   the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
    !   for example, if he knows that the statistical errors on his z(i,j)-
    !   values is not greater than 0.1, he may expect that a good s should
    !   have a value not larger than mx*my*(0.1)**2.
    !   if nothing is known about the statistical error in z(i,j), s must
    !   be determined by trial and error, taking account of the comments
    !   above. the best is then to start with a very large value of s (to
    !   determine the least-squares polynomial and the corresponding upper
    !   bound fp0 for s) and then to progressively decrease the value of s
    !   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
    !   and more carefully as the approximation shows more detail) to
    !   obtain closer fits.
    !   to economize the search for a good s-value the program provides with
    !   different modes of computation. at the first call of the routine, or
    !   whenever he wants to restart with the initial set of knots the user
    !   must set iopt=0.
    !   if iopt=1 the program will continue with the set of knots found at
    !   the last call of the routine. this will save a lot of computation
    !   time if regrid is called repeatedly for different values of s.
    !   the number of knots of the spline returned and their location will
    !   depend on the value of s and on the complexity of the shape of the
    !   function underlying the data. if the computation mode iopt=1
    !   is used, the knots returned may also depend on the s-values at
    !   previous calls (if these were smaller). therefore, if after a number
    !   of trials with different s-values and iopt=1, the user can finally
    !   accept a fit as satisfactory, it may be worthwhile for him to call
    !   regrid once more with the selected value for s but now with iopt=0.
    !   indeed, regrid may then return an approximation of the same quality
    !   of fit but with fewer knots and therefore better if data reduction
    !   is also an important objective for the user.
    !   the number of knots may also depend on the upper bounds nxest and
    !   nyest. indeed, if at a certain stage in regrid the number of knots
    !   in one direction (say nx) has reached the value of its upper bound
    !   (nxest), then from that moment on all subsequent knots are added
    !   in the other (y) direction. this may indicate that the value of
    !   nxest is too small. on the other hand, it gives the user the option
    !   of limiting the number of knots the routine locates in any direction
    !   for example, by setting nxest=2*kx+2 (the lowest allowable value for
    !   nxest), the user can indicate that he wants an approximation which
    !   is a simple polynomial of degree kx in the variable x.
    !
    !  other subroutines required:
    !    fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
    !    fpknot
    !
    !  references:
    !   dierckx p. : a fast algorithm for smoothing data on a rectangular
    !                grid while using spline functions, siam j.numer.anal.
    !                19 (1982) 1286-1304.
    !   dierckx p. : a fast algorithm for smoothing data on a rectangular
    !                grid while using spline functions, report tw53, dept.
    !                computer science,k.u.leuven, 1980.
    !   dierckx p. : curve and surface fitting with splines, monographs on
    !                numerical analysis, oxford university press, 1993.
    !
    !  author:
    !    p.dierckx
    !    dept. computer science, k.u. leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  creation date : may 1979
    !  latest update : march 1989
    !
    !  ..
    !  ..scalar arguments..
    real(kind=rp) xb,xe,yb,ye,s,fp
    integer iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier
    !  ..array arguments..
    real(kind=rp) x(mx),y(my),z(mx*my),tx(nxest),ty(nyest),c((nxest-kx-1)&
         & *(nyest-ky-1)),wrk(lwrk)
    integer iwrk(kwrk)
    !  ..local scalars..
    real(kind=rp) tol
    integer i,j,jwrk,kndx,kndy,knrx,knry,kwest,kx1,kx2,ky1,ky2,lfpx&
         & ,lfpy,lwest,lww,maxit,nc,nminx,nminy,mz
    !  ..function references..
    integer max0
    !  ..subroutine references..
    !    fpregr,fpchec
    !  ..
    !  we set up the parameters tol and maxit.
    maxit = 20
    tol = 0.1e-04
    !  before starting computations a data check is made. if the
    ! input data are invalid, control is immediately repassed to the calling
    ! program.
    ier = 10
    if(kx.le.0 .or. kx.gt.5) go to 70
    kx1 = kx+1
    kx2 = kx1+1
    if(ky.le.0 .or. ky.gt.5) go to 70
    ky1 = ky+1
    ky2 = ky1+1
    if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
    nminx = 2*kx1
    if(mx.lt.kx1 .or. nxest.lt.nminx) go to 70
    nminy = 2*ky1
    if(my.lt.ky1 .or. nyest.lt.nminy) go to 70
    mz = mx*my
    nc = (nxest-kx1)*(nyest-ky1)
    lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+my*ky1&
         & +max0(nxest,my)
    kwest = 3+mx+my+nxest+nyest
    if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 70
    if(xb.gt.x(1) .or. xe.lt.x(mx)) go to 70
    do i=2,mx
       if(x(i-1).ge.x(i)) go to 70
    enddo
    if(yb.gt.y(1) .or. ye.lt.y(my)) go to 70
    do i=2,my
       if(y(i-1).ge.y(i)) go to 70
    enddo
    if(iopt.ge.0) go to 50
    if(nx.lt.nminx .or. nx.gt.nxest) go to 70
    j = nx
    do i=1,kx1
       tx(i) = xb
       tx(j) = xe
       j = j-1
    enddo
    call fpchec(x,mx,tx,nx,kx,ier)
    if(ier.ne.0) go to 70
    if(ny.lt.nminy .or. ny.gt.nyest) go to 70
    j = ny
    do i=1,ky1
       ty(i) = yb
       ty(j) = ye
       j = j-1
    enddo
    call fpchec(y,my,ty,ny,ky,ier)
    if(ier) 70,60,70
50  if(s.lt.0.) go to 70
    if(s.eq.0. .and. (nxest.lt.(mx+kx1) .or. nyest.lt.(my+ky1)) ) goto 70
    ier = 0
    !  we partition the working space and determine the spline approximation
60  lfpx = 5
    lfpy = lfpx+nxest
    lww = lfpy+nyest
    jwrk = lwrk-4-nxest-nyest
    knrx = 4
    knry = knrx+mx
    kndx = knry+my
    kndy = kndx+nxest
    call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest&
         & ,tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4)&
         & ,wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx)&
         & ,iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)
70  return
  end subroutine regrid
  !
  !
  !
  subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
       !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
    !  the derivative of order nu of a spline s(x) of degree k,given in
    !  its b-spline representation.
    !
    !  calling sequence:
    !     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
    !
    !  input parameters:
    !    t    : array,length n, which contains the position of the knots.
    !    n    : integer, giving the total number of knots of s(x).
    !    c    : array,length n, which contains the b-spline coefficients.
    !    k    : integer, giving the degree of s(x).
    !    nu   : integer, specifying the order of the derivative. 0<=nu<=k
    !    x    : array,length m, which contains the points where the deriv-
    !           ative of s(x) must be evaluated.
    !    m    : integer, giving the number of points where the derivative
    !           of s(x) must be evaluated
    !    wrk  : real array of dimension n. used as working space.
    !
    !  output parameters:
    !    y    : array,length m, giving the value of the derivative of s(x)
    !           at the different points.
    !    ier  : error flag
    !      ier = 0 : normal return
    !      ier =10 : invalid input data (see restrictions)
    !
    !  restrictions:
    !    0 <= nu <= k
    !    m >= 1
    !    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
    !
    !  other subroutines required: fpbspl
    !
    !  references :
    !    de boor c : on calculating with b-splines, j. approximation theory
    !                6 (1972) 50-62.
    !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
    !                applics 10 (1972) 134-149.
    !   dierckx p. : curve and surface fitting with splines, monographs on
    !                numerical analysis, oxford university press, 1993.
    !
    !  author :
    !    p.dierckx
    !    dept. computer science, k.u.leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  latest update : march 1987
    !
    !  ..scalar arguments..
    integer n,k,nu,m,ier
    !  ..array arguments..
    real(kind=rp) t(n),c(n),x(m),y(m),wrk(n)
    !  ..local scalars..
    integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
    real(kind=rp) ak,arg,fac,sp,tb,te
    !  ..local arrays ..
    real(kind=rp) h(6)
    !  before starting computations a data check is made. if the input data
    !  are invalid control is immediately repassed to the calling program.
    ier = 10
    if(nu.lt.0 .or. nu.gt.k) go to 200
    if(m-1) 200,30,10
10  do i=2,m
       if(x(i).lt.x(i-1)) go to 200
    enddo
30  ier = 0
    !  fetch tb and te, the boundaries of the approximation interval.
    k1 = k+1
    nk1 = n-k1
    tb = t(k1)
    te = t(nk1+1)
    !  the derivative of order nu of a spline of degree k is a spline of
    !  degree k-nu,the b-spline coefficients wrk(i) of which can be found
    !  using the recurrence scheme of de boor.
    l = 1
    kk = k
    nn = n
    do i=1,nk1
       wrk(i) = c(i)
    enddo
    if(nu.eq.0) go to 100
    nk2 = nk1
    do j=1,nu
       ak = kk
       nk2 = nk2-1
       l1 = l
       do 50 i=1,nk2
          l1 = l1+1
          l2 = l1+kk
          fac = t(l2)-t(l1)
          if(fac.le.zero) go to 50
          wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50   continue
       l = l+1
       kk = kk-1
    enddo
    if(kk.ne.0) go to 100
    !  if nu=k the derivative is a piecewise constant function
    j = 1
    do i=1,m
       arg = x(i)
70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
       l = l+1
       j = j+1
       go to 70
80     y(i) = wrk(j)
    enddo
    go to 200
100 l = k1
    l1 = l+1
    k2 = k1-nu
    !  main loop for the different points.
    do i=1,m
       !  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
        !  search for knot interval t(l) <= arg < t(l+1)
140     if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
        !  evaluate the non-zero b-splines of degree k-nu at arg.
150     call fpbspl(t,n,kk,arg,l,h)
        !  find the value of the derivative at x=arg.
        sp = 0.
        ll = l-k1
        do j=1,k2
           ll = ll+1
           sp = sp+wrk(ll)*h(j)
        enddo
        y(i) = sp
    enddo
200 return
  end subroutine splder
  !
  !
  !
  subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest&
       & ,tol,maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,reducy,fpintx&
       & ,fpinty,lastdi,nplusx,nplusy,nrx,nry,nrdatx,nrdaty,wrk,lwrk&
       & ,ier)
       !  ..
    !  ..scalar arguments..
    real(kind=rp) xb,xe,yb,ye,s,tol,fp,fp0,fpold,reducx,reducy
    integer iopt,mx,my,mz,kx,ky,nxest,nyest,maxit,nc,nx,ny,lastdi&
         & ,nplusx,nplusy,lwrk,ier
    !  ..array arguments..
    real(kind=rp) x(mx),y(my),z(mz),tx(nxest),ty(nyest),c(nc),fpintx(nxest)&
         & ,fpinty(nyest),wrk(lwrk)
    integer nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
    !  ..local scalars
    real(kind=rp) acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,half,con1,con9,con4
    integer i,ich1,ich3,ifbx,ifby,ifsx,ifsy,iter,j,kx1,kx2,ky1,ky2,k3&
         & ,l,lax,lay,lbx,lby,lq,lri,lsx,lsy,mk1,mm,mpm,mynx,ncof&
         & ,nk1x,nk1y,nmaxx,nmaxy,nminx,nminy,nplx,nply,npl1,nrintx&
         & ,nrinty,nxe,nxk,nye
    !  ..function references..
    real(kind=rp) abs
    integer max0,min0
    !  ..subroutine references..
    !    fpgrre,fpknot
    !  ..
    !   set constants
    one = 1
    half = 0.5E0_rp
    con1 = 0.1E0_rp
    con9 = 0.9E0_rp
    con4 = 0.4E-01_rp
    !  we partition the working space.
    kx1 = kx+1
    ky1 = ky+1
    kx2 = kx1+1
    ky2 = ky1+1
    lsx = 1
    lsy = lsx+mx*kx1
    lri = lsy+my*ky1
    mm = max0(nxest,my)
    lq = lri+mm
    mynx = nxest*my
    lax = lq+mynx
    nxk = nxest*kx2
    lbx = lax+nxk
    lay = lbx+nxk
    lby = lay+nyest*ky2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! part 1: determination of the number of knots and their position.     !
    ! ****************************************************************     !
    !  given a set of knots we compute the least-squares spline sinf(x,y), !
    !  and the corresponding sum of squared residuals fp=f(p=inf).         !
    !  if iopt=-1  sinf(x,y) is the requested approximation.               !
    !  if iopt=0 or iopt=1 we check whether we can accept the knots:       !
    !    if fp <=s we will continue with the current set of knots.         !
    !    if fp > s we will increase the number of knots and compute the    !
    !       corresponding least-squares spline until finally fp<=s.        !
    !    the initial choice of knots depends on the value of s and iopt.   !
    !    if s=0 we have spline interpolation; in that case the number of   !
    !    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               !
    !    if s>0 and                                                        !
    !     *iopt=0 we first compute the least-squares polynomial of degree  !
    !      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       !
    !     *iopt=1 we start with the knots found at the last call of the    !
    !      routine, except for the case that s > fp0; then we can compute  !
    !      the least-squares polynomial directly.                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  determine the number of knots for polynomial approximation.
    nminx = 2*kx1
    nminy = 2*ky1
    if(iopt.lt.0) go to 120
    !  acc denotes the absolute tolerance for the root of f(p)=s.
    acc = tol*s
    !  find nmaxx and nmaxy which denote the number of knots in x- and y-
    !  direction in case of spline interpolation.
    nmaxx = mx+kx1
    nmaxy = my+ky1
    !  find nxe and nye which denote the maximum number of knots
    !  allowed in each direction
    nxe = min0(nmaxx,nxest)
    nye = min0(nmaxy,nyest)
    if(s.gt.zero) go to 100
    !  if s = 0, s(x,y) is an interpolating spline.
    nx = nmaxx
    ny = nmaxy
    !  test whether the required storage space exceeds the available one.
    if(ny.gt.nyest .or. nx.gt.nxest) go to 420
    !  find the position of the interior knots in case of interpolation.
    !  the knots in the x-direction.
    mk1 = mx-kx1
    if(mk1.eq.0) go to 60
    k3 = kx/2
    i = kx1+1
    j = k3+2
    if(k3*2.eq.kx) go to 40
    do l=1,mk1
       tx(i) = x(j)
       i = i+1
       j = j+1
    enddo
    go to 60
40  do l=1,mk1
       tx(i) = (x(j)+x(j-1))*half
       i = i+1
       j = j+1
    enddo
    !  the knots in the y-direction.
60  mk1 = my-ky1
    if(mk1.eq.0) go to 120
    k3 = ky/2
    i = ky1+1
    j = k3+2
    if(k3*2.eq.ky) go to 80
    do l=1,mk1
       ty(i) = y(j)
       i = i+1
       j = j+1
    enddo
    go to 120
80  do l=1,mk1
       ty(i) = (y(j)+y(j-1))*half
       i = i+1
       j = j+1
    enddo
    go to 120
    !  if s > 0 our initial choice of knots depends on the value of iopt.
100 if(iopt.eq.0) go to 115
    if(fp0.le.s) go to 115
    !  if iopt=1 and fp0 > s we start computing the least- squares spline
    !  according to the set of knots found at the last call of the routine.
    !  we determine the number of grid coordinates x(i) inside each knot
    !  interval (tx(l),tx(l+1)).
    l = kx2
    j = 1
    nrdatx(1) = 0
    mpm = mx-1
    do 105 i=2,mpm
       nrdatx(j) = nrdatx(j)+1
       if(x(i).lt.tx(l)) go to 105
       nrdatx(j) = nrdatx(j)-1
       l = l+1
       j = j+1
       nrdatx(j) = 0
105 continue
    !  we determine the number of grid coordinates y(i) inside each knot
    !  interval (ty(l),ty(l+1)).
    l = ky2
    j = 1
    nrdaty(1) = 0
    mpm = my-1
    do 110 i=2,mpm
       nrdaty(j) = nrdaty(j)+1
       if(y(i).lt.ty(l)) go to 110
       nrdaty(j) = nrdaty(j)-1
       l = l+1
       j = j+1
       nrdaty(j) = 0
110 continue
    go to 120
    !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
    !  polynomial of degree kx in x and ky in y (which is a spline without
    !  interior knots).
115 nx = nminx
    ny = nminy
    nrdatx(1) = mx-2
    nrdaty(1) = my-2
    lastdi = 0
    nplusx = 0
    nplusy = 0
    fp0 = zero
    fpold = zero
    reducx = zero
    reducy = zero
120 mpm = mx+my
    ifsx = 0
    ifsy = 0
    ifbx = 0
    ifby = 0
    p = -one
    !  main loop for the different sets of knots.mpm=mx+my is a save upper
    !  bound for the number of trials.
    do 250 iter=1,mpm
       if(nx.eq.nminx .and. ny.eq.nminy) ier = -2
       !  find nrintx (nrinty) which is the number of knot intervals in the
       !  x-direction (y-direction).
       nrintx = nx-nminx+1
       nrinty = ny-nminy+1
       !  find ncof, the number of b-spline coefficients for the current set
       !  of knots.
       nk1x = nx-kx1
       nk1y = ny-ky1
       ncof = nk1x*nk1y
       !  find the position of the additional knots which are needed for the
       !  b-spline representation of s(x,y).
       i = nx
       do j=1,kx1
          tx(j) = xb
          tx(i) = xe
          i = i-1
       enddo
       i = ny
       do j=1,ky1
          ty(j) = yb
          ty(i) = ye
          i = i-1
       enddo
       !  find the least-squares spline sinf(x,y) and calculate for each knot
       !  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
       !  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
       !  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
       !  value belonging to that interval.
       !  fp gives the total sum of squared residuals.
       call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty&
            & ,ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2&
            & ,wrk(lsx),wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay)&
            & ,wrk(lbx),wrk(lby),nrx,nry)
       if(ier.eq.(-2)) fp0 = fp
       !  test whether the least-squares spline is an acceptable solution.
       if(iopt.lt.0) go to 440
       fpms = fp-s
       if(abs(fpms) .lt. acc) go to 440
       !  if f(p=inf) < s, we accept the choice of knots.
       if(fpms.lt.0.) go to 300
       !  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating
       ! spline.
       if(nx.eq.nmaxx .and. ny.eq.nmaxy) go to 430
       !  increase the number of knots.
       !  if nx=nxe and ny=nye we cannot further increase the number
       !  of knots because of the storage capacity limitation.
       if(nx.eq.nxe .and. ny.eq.nye) go to 420
       ier = 0
       !  adjust the parameter reducx or reducy according to the
       ! direction in which the last added knots were located.
       if(lastdi) 150,170,160
150    reducx = fpold-fp
       go to 170
160    reducy = fpold-fp
       !  store the sum of squared residuals for the current set of knots.
170    fpold = fp
       !  find nplx, the number of knots we should add in the x-direction.
       nplx = 1
       if(nx.eq.nminx) go to 180
       npl1 = nplusx*2
       rn = nplusx
       if(reducx.gt.acc) npl1 = rn*fpms/reducx
       nplx = min0(nplusx*2,max0(npl1,nplusx/2,1))
       !  find nply, the number of knots we should add in the y-direction.
180    nply = 1
       if(ny.eq.nminy) go to 190
       npl1 = nplusy*2
       rn = nplusy
       if(reducy.gt.acc) npl1 = rn*fpms/reducy
       nply = min0(nplusy*2,max0(npl1,nplusy/2,1))
190    if(nplx-nply) 210,200,230
200    if(lastdi.lt.0) go to 230
210    if(nx.eq.nxe) go to 230
       !  addition in the x-direction.
       lastdi = -1
       nplusx = nplx
       ifsx = 0
       do l=1,nplusx
          !  add a new knot in the x-direction
          call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1)
          !  test whether we cannot further increase the number of
          ! knots in the
          !  x-direction.
          if(nx.eq.nxe) go to 250
       enddo
       go to 250
230    if(ny.eq.nye) go to 210
       !  addition in the y-direction.
       lastdi = 1
       nplusy = nply
       ifsy = 0
       do l=1,nplusy
          !  add a new knot in the y-direction.
          call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1)
          !  test whether we cannot further increase the number
          ! of knots in the y-direction.
          if(ny.eq.nye) go to 250
       enddo
       !  restart the computations with the new set of knots.
250 continue
    !  test whether the least-squares polynomial is a solution of
    !  our approximation problem.
300 if(ier.eq.(-2)) go to 440
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! part 2: determination of the smoothing spline sp(x,y)                 !
    ! *****************************************************                 !
    !  we have determined the number of knots and their position. we now    !
    !  compute the b-spline coefficients of the smoothing spline sp(x,y).   !
    !  this smoothing spline varies with the parameter p in such a way that !
    !    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)              !
    !  is a continuous, strictly decreasing function of p. moreover the     !
    !  least-squares polynomial corresponds to p=0 and the least-squares    ! 
    !  spline to p=infinity. iteratively we then have to determine the      !
    !  positive value of p such that f(p)=s. the process which is proposed  !
    !  here makes use of rational interpolation. f(p) is approximated by a  !
    !  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)   !
    !  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) !
    !  are used to calculate the new value of p such that r(p)=s.           !
    !  convergence is guaranteed by taking f1 > 0 and f3 < 0.               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  initial value for p.
    p1 = zero
    f1 = fp0-s
    p3 = -one
    f3 = fpms
    p = one
    ich1 = 0
    ich3 = 0
    !  iteration process to find the root of f(p)=s.
    do 350 iter = 1,maxit
       !  find the smoothing spline sp(x,y) and the corresponding sum of
       !  squared residuals fp.
       call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty&
            & ,ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2&
            & ,wrk(lsx),wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay)&
            & ,wrk(lbx),wrk(lby),nrx,nry)
       !  test whether the approximation sp(x,y) is an acceptable solution.
       fpms = fp-s
       if(abs(fpms).lt.acc) go to 440
       !  test whether the maximum allowable number of iterations has been
       !  reached.
       if(iter.eq.maxit) go to 400
       !  carry out one more step of the iteration process.
       p2 = p
       f2 = fpms
       if(ich3.ne.0) go to 320
       if((f2-f3).gt.acc) go to 310
       !  our initial choice of p is too large.
       p3 = p2
       f3 = f2
       p = p*con4
       if(p.le.p1) p = p1*con9 + p2*con1
       go to 350
310    if(f2.lt.zero) ich3 = 1
320    if(ich1.ne.0) go to 340
       if((f1-f2).gt.acc) go to 330
       !  our initial choice of p is too small
       p1 = p2
       f1 = f2
       p = p/con4
       if(p3.lt.zero) go to 350
       if(p.ge.p3) p = p2*con1 + p3*con9
       go to 350
       !  test whether the iteration process proceeds as theoretically
       !  expected.
330    if(f2.gt.zero) ich1 = 1
340    if(f2.ge.f1 .or. f2.le.f3) go to 410
       !  find the new value of p.
       p = fprati(p1,f1,p2,f2,p3,f3)
350 continue
    !  error codes and messages.
400 ier = 3
    go to 440
410 ier = 2
    go to 440
420 ier = 1
    go to 440
430 ier = -1
    fp = zero
440 return
  end subroutine fpregr
  !
  !
  !
  subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty&
       & ,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right&
       & ,q,ax,ay,bx,by,nrx,nry)
       !  ..
    !  ..scalar arguments..
    real(kind=rp) p,fp
    integer ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,kx1&
         & ,kx2,ky1,ky2
    !  ..array arguments..
    real(kind=rp) x(mx),y(my),z(mz),tx(nx),ty(ny),c(nc),spx(mx,kx1),spy(my&
         & ,ky1),right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),ay(ny,ky2)&
         & ,by(ny,ky2),fpx(nx),fpy(ny)
    integer nrx(mx),nry(my)
    !  ..local scalars..
    real(kind=rp) arg,cos,fac,pinv,piv,sin,term,one,half
    integer i,ibandx,ibandy,ic,iq,irot,it,iz,i1,i2,i3,j,k,k1,k2,l,l1&
         & ,l2,ncof,nk1x,nk1y,nrold,nroldx,nroldy,number,numx,numx1&
         & ,numy,numy1,n1
    !  ..local arrays..
    real(kind=rp) h(7)
    !  ..subroutine references..
    !    fpback,fpbspl,fpgivs,fpdisc,fprota
    !  ..
    !  the b-spline coefficients of the smoothing spline are calculated as
    !  the least-squares solution of the over-determined linear system of
    !  equations  (ay) c (ax)' = q       where
    !
    !               |   (spx)    |            |   (spy)    |
    !        (ax) = | ---------- |     (ay) = | ---------- |
    !               | (1/p) (bx) |            | (1/p) (by) |
    !
    !                                | z  ' 0 |
    !                            q = | ------ |
    !                                | 0  ' 0 |
    !
    !  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
    !                b-spline coefficients.
    !       z      : the my x mx matrix which contains the function values.
    !       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
    !                matrices according to the least-squares problems in
    !                the x- and y-direction.
    !       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
    !                matrices which contain the discontinuity jumps of the
    !                derivatives of the b-splines in the x- and y-direction.
    one = 1
    half = 0.5_rp
    nk1x = nx-kx1
    nk1y = ny-ky1
    if(p.gt.zero) pinv = one/p
    !  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
    !  the value of p whether the matrices (spx),(spy),(bx) and (by) still
    !  must be determined.
    if(ifsx.ne.0) go to 50
    !  calculate the non-zero elements of the matrix (spx) which is the
    !  observation matrix according to the least-squares spline approximat-
    !  ion problem in the x-direction.
    l = kx1
    l1 = kx2
    number = 0
    do it=1,mx
       arg = x(it)
10     if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
       l = l1
       l1 = l+1
       number = number+1
       go to 10
20     call fpbspl(tx,nx,kx,arg,l,h)
       do i=1,kx1
          spx(it,i) = h(i)
       enddo
       nrx(it) = number
    enddo
    ifsx = 1
50  if(ifsy.ne.0) go to 100
    !  calculate the non-zero elements of the matrix (spy) which is the
    !  observation matrix according to the least-squares spline approximat-
    !  ion problem in the y-direction.
    l = ky1
    l1 = ky2
    number = 0
    do it=1,my
       arg = y(it)
60     if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
       l = l1
       l1 = l+1
       number = number+1
       go to 60
70     call fpbspl(ty,ny,ky,arg,l,h)
       do i=1,ky1
          spy(it,i) = h(i)
       enddo
       nry(it) = number
    enddo
    ifsy = 1
100 if(p.le.zero) go to 120
    !  calculate the non-zero elements of the matrix (bx).
    if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110
    call fpdisc(tx,nx,kx2,bx,nx)
    ifbx = 1
    !  calculate the non-zero elements of the matrix (by).
110 if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
    call fpdisc(ty,ny,ky2,by,ny)
    ifby = 1
    !  reduce the matrix (ax) to upper triangular form (rx) using givens
    !  rotations. apply the same transformations to the rows of matrix q
    !  to obtain the my x (nx-kx-1) matrix g.
    !  store matrix (rx) into (ax) and g into q.
120 l = my*nk1x
    !  initialization.
    do i=1,l
       q(i) = zero
    enddo
    do i=1,nk1x
       do j=1,kx2
          ax(i,j) = zero
       enddo
    enddo
    l = 0
    nrold = 0
    !  ibandx denotes the bandwidth of the matrices (ax) and (rx).
    ibandx = kx1
    do 270 it=1,mx
       number = nrx(it)
150    if(nrold.eq.number) go to 180
       if(p.le.zero) go to 260
       ibandx = kx2
       !  fetch a new row of matrix (bx).
       n1 = nrold+1
       do j=1,kx2
          h(j) = bx(n1,j)*pinv
       enddo
       !  find the appropriate column of q.
       do j=1,my
          right(j) = zero
       enddo
       irot = nrold
       go to 210
       !  fetch a new row of matrix (spx).
180    h(ibandx) = zero
       do j=1,kx1
          h(j) = spx(it,j)
       enddo
       !  find the appropriate column of q.
       do j=1,my
          l = l+1
          right(j) = z(l)
       enddo
       irot = number
       !  rotate the new row of matrix (ax) into triangle.
210    do 240 i=1,ibandx
          irot = irot+1
          piv = h(i)
          if(piv.eq.zero) go to 240
          !  calculate the parameters of the givens transformation.
          call fpgivs(piv,ax(irot,1),cos,sin)
          !  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do j=1,my
             iq = iq+1
             call fprota(cos,sin,right(j),q(iq))
          enddo
          !  apply that transformation to the columns of (ax).
          if(i.eq.ibandx) go to 250
          i2 = 1
          i3 = i+1
          do j=i3,ibandx
             i2 = i2+1
             call fprota(cos,sin,h(j),ax(irot,i2))
          enddo
 240   continue
 250   if(nrold.eq.number) go to 270
 260   nrold = nrold+1
       go to 150
270 continue
    !  reduce the matrix (ay) to upper triangular form (ry) using givens
    !  rotations. apply the same transformations to the columns of matrix g
    !  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
    !  store matrix (ry) into (ay) and h into c.
    ncof = nk1x*nk1y
    !  initialization.
    do i=1,ncof
       c(i) = zero
    enddo
    do i=1,nk1y
       do j=1,ky2
          ay(i,j) = zero
       enddo
    enddo
    nrold = 0
    !  ibandy denotes the bandwidth of the matrices (ay) and (ry).
    ibandy = ky1
    do 420 it=1,my
       number = nry(it)
300    if(nrold.eq.number) go to 330
       if(p.le.zero) go to 410
       ibandy = ky2
       !  fetch a new row of matrix (by).
       n1 = nrold+1
       do j=1,ky2
          h(j) = by(n1,j)*pinv
       enddo
       !  find the appropiate row of g.
       do j=1,nk1x
          right(j) = zero
       enddo
       irot = nrold
       go to 360
       !  fetch a new row of matrix (spy)
330    h(ibandy) = zero
       do j=1,ky1
          h(j) = spy(it,j)
       enddo
       !  find the appropiate row of g.
       l = it
       do j=1,nk1x
          right(j) = q(l)
          l = l+my
       enddo
       irot = number
       !  rotate the new row of matrix (ay) into triangle.
360    do 390 i=1,ibandy
          irot = irot+1
          piv = h(i)
          if(piv.eq.zero) go to 390
          !  calculate the parameters of the givens transformation.
          call fpgivs(piv,ay(irot,1),cos,sin)
          !  apply that transformation to the colums of matrix g.
          ic = irot
          do j=1,nk1x
             call fprota(cos,sin,right(j),c(ic))
             ic = ic+nk1y
          enddo
          !  apply that transformation to the columns of matrix (ay).
          if(i.eq.ibandy) go to 400
          i2 = 1
          i3 = i+1
          do j=i3,ibandy
             i2 = i2+1
             call fprota(cos,sin,h(j),ay(irot,i2))
          enddo
390    continue
400       if(nrold.eq.number) go to 420
410       nrold = nrold+1
        go to 300
420 continue
    !  backward substitution to obtain the b-spline coefficients as the
    !  solution of the linear system    (ry) c (rx)' = h.
    !  first step: solve the system  (ry) (c1) = h.
    k = 1
    do i=1,nk1x
       call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
       k = k+nk1y
    enddo
    !  second step: solve the system  c (rx)' = (c1).
    k = 0
    do j=1,nk1y
       k = k+1
       l = k
       do i=1,nk1x
          right(i) = c(l)
          l = l+nk1y
       enddo
       call fpback(ax,right,nk1x,ibandx,right,nx)
       l = k
       do i=1,nk1x
          c(l) = right(i)
          l = l+nk1y
       enddo
    enddo
    !  calculate the quantities
    !    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
    !    fp = sumi=1,mx(sumj=1,my(res(i,j)))
    !    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
    !                  tx(r+kx) <= x(i) <= tx(r+kx+1)
    !    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
    !                  ty(r+ky) <= y(j) <= ty(r+ky+1)
    fp = zero
    do i=1,nx
       fpx(i) = zero
    enddo
    do i=1,ny
       fpy(i) = zero
    enddo
    nk1y = ny-ky1
    iz = 0
    nroldx = 0
    !  main loop for the different grid points.
    do 550 i1=1,mx
       numx = nrx(i1)
       numx1 = numx+1
       nroldy = 0
       do 540 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
          iz = iz+1
          !  evaluate s(x,y) at the current grid point by making the sum of the
          !  cross products of the non-zero b-splines at (x,y), multiplied with
          !  the appropiate b-spline coefficients.
          term = zero
          k1 = numx*nk1y+numy
          do l1=1,kx1
             k2 = k1
             fac = spx(i1,l1)
             do l2=1,ky1
                k2 = k2+1
                term = term+fac*spy(i2,l2)*c(k2)
             enddo
             k1 = k1+nk1y
          enddo
          !  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
          !  adjust the different parameters.
          fp = fp+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*half
          if(numy.eq.nroldy) go to 530
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
530       nroldy = numy
          if(numx.eq.nroldx) go to 540
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
540    continue
          nroldx = numx
550 continue
    return
  end subroutine fpgrre
  !
  !
  !
  subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest&
       & ,nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk&
       & ,ier)
    ! given the set of data points (x(i),y(i),z(i)) and the set of positive
    ! numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
    ! iate spline approximation s(x,y) of degrees kx and ky on the rect-
    ! angle xb <= x <= xe, yb <= y <= ye.
    ! if iopt = -1 surfit calculates the weighted least-squares spline
    ! according to a given set of knots.
    ! if iopt >= 0 the total numbers nx and ny of these knots and their
    ! position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
    ! ally by the routine. the smoothness of s(x,y) is then achieved by
    ! minimalizing the discontinuity jumps in the derivatives of s(x,y)
    ! across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
    ! the amounth of smoothness is determined by the condition that f(p) =
    ! sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
    ! ative constant, called the smoothing factor.
    ! the fit is given in the b-spline representation (b-spline coefficients
    ! c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
    ! uated by means of subroutine bispev.
    !
    ! calling sequence:
    !     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
    !    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
    !
    ! parameters:
    !  iopt  : integer flag. on entry iopt must specify whether a weighted
    !          least-squares spline (iopt=-1) or a smoothing spline (iopt=0
    !          or 1) must be determined.
    !          if iopt=0 the routine will start with an initial set of knots
    !          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
    !          1,...,ky+1. if iopt=1 the routine will continue with the set
    !          of knots found at the last call of the routine.
    !          attention: a call with iopt=1 must always be immediately pre-
    !                     ceded by another call with iopt=1 or iopt=0.
    !          unchanged on exit.
    !  m     : integer. on entry m must specify the number of data points.
    !          m >= (kx+1)*(ky+1). unchanged on exit.
    !  x     : real array of dimension at least (m).
    !  y     : real array of dimension at least (m).
    !  z     : real array of dimension at least (m).
    !          before entry, x(i),y(i),z(i) must be set to the co-ordinates
    !          of the i-th data point, for i=1,...,m. the order of the data
    !          points is immaterial. unchanged on exit.
    !  w     : real array of dimension at least (m). before entry, w(i) must
    !          be set to the i-th value in the set of weights. the w(i) must
    !          be strictly positive. unchanged on exit.
    !  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
    !  yb,ye   aries of the rectangular approximation domain.
    !          xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
    !  kx,ky : integer values. on entry kx and ky must specify the degrees
    !          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
    !          (kx=ky=3) splines. unchanged on exit.
    !  s     : real. on entry (in case iopt>=0) s must specify the smoothing
    !          factor. s >=0. unchanged on exit.
    !          for advice on the choice of s see further comments
    !  nxest : integer. unchanged on exit.
    !  nyest : integer. unchanged on exit.
    !          on entry, nxest and nyest must specify an upper bound for the
    !          number of knots required in the x- and y-directions respect.
    !          these numbers will also determine the storage space needed by
    !          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
    !          in most practical situation nxest = kx+1+sqrt(m/2), nyest =
    !          ky+1+sqrt(m/2) will be sufficient. see also further comments.
    !  nmax  : integer. on entry nmax must specify the actual dimension of
    !          the arrays tx and ty. nmax >= nxest, nmax >=nyest.
    !          unchanged on exit.
    !  eps   : real.
    !          on entry, eps must specify a threshold for determining the
    !          effective rank of an over-determined linear system of equat-
    !          ions. 0 < eps < 1.  if the number of decimal digits in the
    !          computer representation of a real number is q, then 10**(-q)
    !          is a suitable value for eps in most practical applications.
    !          unchanged on exit.
    !  nx    : integer.
    !          unless ier=10 (in case iopt >=0), nx will contain the total
    !          number of knots with respect to the x-variable, of the spline
    !          approximation returned. if the computation mode iopt=1 is
    !          used, the value of nx should be left unchanged between sub-
    !          sequent calls.
    !          in case iopt=-1, the value of nx should be specified on entry
    !  tx    : real array of dimension nmax.
    !          on succesful exit, this array will contain the knots of the
    !          spline with respect to the x-variable, i.e. the position of
    !          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
    !          position of the additional knots tx(1)=...=tx(kx+1)=xb and
    !          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
    !          if the computation mode iopt=1 is used, the values of tx(1),
    !          ...,tx(nx) should be left unchanged between subsequent calls.
    !          if the computation mode iopt=-1 is used, the values tx(kx+2),
    !          ...tx(nx-kx-1) must be supplied by the user, before entry.
    !          see also the restrictions (ier=10).
    !  ny    : integer.
    !          unless ier=10 (in case iopt >=0), ny will contain the total
    !          number of knots with respect to the y-variable, of the spline
    !          approximation returned. if the computation mode iopt=1 is
    !          used, the value of ny should be left unchanged between sub-
    !          sequent calls.
    !          in case iopt=-1, the value of ny should be specified on entry
    !  ty    : real array of dimension nmax.
    !          on succesful exit, this array will contain the knots of the
    !          spline with respect to the y-variable, i.e. the position of
    !          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
    !          position of the additional knots ty(1)=...=ty(ky+1)=yb and
    !          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
    !          if the computation mode iopt=1 is used, the values of ty(1),
    !          ...,ty(ny) should be left unchanged between subsequent calls.
    !          if the computation mode iopt=-1 is used, the values ty(ky+2),
    !          ...ty(ny-ky-1) must be supplied by the user, before entry.
    !          see also the restrictions (ier=10).
    !  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
    !          on succesful exit, c contains the coefficients of the spline
    !          approximation s(x,y)
    !  fp    : real. unless ier=10, fp contains the weighted sum of
    !          squared residuals of the spline approximation returned.
    !  wrk1  : real array of dimension (lwrk1). used as workspace.
    !          if the computation mode iopt=1 is used the value of wrk1(1)
    !          should be left unchanged between subsequent calls.
    !          on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
    !          contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
    !          with d(i) the i-th diagonal element of the reduced triangular
    !          matrix for calculating the b-spline coefficients. it includes
    !          those elements whose square is less than eps,which are treat-
    !          ed as 0 in the case of presumed rank deficiency (ier<-2).
    !  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
    !          the array wrk1 as declared in the calling (sub)program.
    !          lwrk1 must not be too small. let
    !            u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
    !            ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
    !            if(bx.le.by) b1 = bx, b2 = b1+v-ky
    !            if(bx.gt.by) b1 = by, b2 = b1+u-kx  then
    !          lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
    !  wrk2  : real array of dimension (lwrk2). used as workspace, but
    !          only in the case a rank deficient system is encountered.
    !  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
    !          the array wrk2 as declared in the calling (sub)program.
    !          lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2
    !          where u,v and b2 are as above. if there are enough data
    !          points, scattered uniformly over the approximation domain
    !          and if the smoothing factor s is not too small, there is a
    !          good chance that this extra workspace is not needed. a lot
    !          of memory might therefore be saved by setting lwrk2=1.
    !          (see also ier > 10)
    !  iwrk  : integer array of dimension (kwrk). used as workspace.
    !  kwrk  : integer. on entry kwrk must specify the actual dimension of
    !          the array iwrk as declared in the calling (sub)program.
    !          kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
    !  ier   : integer. unless the routine detects an error, ier contains a
    !          non-positive value on exit, i.e.
    !   ier=0  : normal return. the spline returned has a residual sum of
    !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
    !            ive tolerance set to 0.001 by the program.
    !   ier=-1 : normal return. the spline returned is an interpolating
    !            spline (fp=0).
    !   ier=-2 : normal return. the spline returned is the weighted least-
    !            squares polynomial of degrees kx and ky. in this extreme
    !            case fp gives the upper bound for the smoothing factor s.
    !   ier<-2 : warning. the coefficients of the spline returned have been
    !            computed as the minimal norm least-squares solution of a
    !            (numerically) rank deficient system. (-ier) gives the rank.
    !            especially if the rank deficiency which can be computed as
    !            (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
    !            curate. they could also seriously depend on the value of
    !            eps.
    !   ier=1  : error. the required storage space exceeds the available
    !            storage space, as specified by the parameters nxest and
    !            nyest.
    !            probably causes : nxest or nyest too small. if these param-
    !            eters are already large, it may also indicate that s is
    !            too small
    !            the approximation returned is the weighted least-squares
    !            spline according to the current set of knots.
    !            the parameter fp gives the corresponding weighted sum of
    !            squared residuals (fp>s).
    !   ier=2  : error. a theoretically impossible result was found during
    !            the iteration proces for finding a smoothing spline with
    !            fp = s. probably causes : s too small or badly chosen eps.
    !            there is an approximation returned but the corresponding
    !            weighted sum of squared residuals does not satisfy the
    !            condition abs(fp-s)/s < tol.
    !   ier=3  : error. the maximal number of iterations maxit (set to 20
    !            by the program) allowed for finding a smoothing spline
    !            with fp=s has been reached. probably causes : s too small
    !            there is an approximation returned but the corresponding
    !            weighted sum of squared residuals does not satisfy the
    !            condition abs(fp-s)/s < tol.
    !   ier=4  : error. no more knots can be added because the number of
    !            b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
    !            the number of data points m.
    !            probably causes : either s or m too small.
    !            the approximation returned is the weighted least-squares
    !            spline according to the current set of knots.
    !            the parameter fp gives the corresponding weighted sum of
    !            squared residuals (fp>s).
    !   ier=5  : error. no more knots can be added because the additional
    !            knot would (quasi) coincide with an old one.
    !            probably causes : s too small or too large a weight to an
    !            inaccurate data point.
    !            the approximation returned is the weighted least-squares
    !            spline according to the current set of knots.
    !            the parameter fp gives the corresponding weighted sum of
    !            squared residuals (fp>s).
    !   ier=10 : error. on entry, the input data are controlled on validity
    !            the following restrictions must be satisfied.
    !            -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
    !            nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
    !            xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
    !            lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
    !            kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
    !            if iopt=-1: 2*kx+2<=nx<=nxest
    !                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
    !                        2*ky+2<=ny<=nyest
    !                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
    !            if iopt>=0: s>=0
    !            if one of these conditions is found to be violated,control
    !            is immediately repassed to the calling program. in that
    !            case there is no approximation returned.
    !   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
    !            space for computing the minimal least-squares solution of
    !            a rank deficient system of linear equations. ier gives the
    !            requested value for lwrk2. there is no approximation re-
    !            turned but, having saved the information contained in nx,
    !            ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
    !            the dimension of the array wrk2 accordingly, the user can
    !            continue at the point the program was left, by calling
    !            surfit with iopt=1.
    !
    ! further comments:
    !  by means of the parameter s, the user can control the tradeoff
    !   between closeness of fit and smoothness of fit of the approximation.
    !   if s is too large, the spline will be too smooth and signal will be
    !   lost ; if s is too small the spline will pick up too much noise. in
    !   the extreme cases the program will return an interpolating spline if
    !   s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
    !   very large. between these extremes, a properly chosen s will result
    !   in a good compromise between closeness of fit and smoothness of fit.
    !   to decide whether an approximation, corresponding to a certain s is
    !   satisfactory the user is highly recommended to inspect the fits
    !   graphically.
    !   recommended values for s depend on the weights w(i). if these are
    !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
    !   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
    !   sqrt(2*m)). if nothing is known about the statistical error in z(i)
    !   each w(i) can be set equal to one and s determined by trial and
    !   error, taking account of the comments above. the best is then to
    !   start with a very large value of s ( to determine the least-squares
    !   polynomial and the corresponding upper bound fp0 for s) and then to
    !   progressively decrease the value of s ( say by a factor 10 in the
    !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
    !   approximation shows more detail) to obtain closer fits.
    !   to choose s very small is strongly discouraged. this considerably
    !   increases computation time and memory requirements. it may also
    !   cause rank-deficiency (ier<-2) and endager numerical stability.
    !   to economize the search for a good s-value the program provides with
    !   different modes of computation. at the first call of the routine, or
    !   whenever he wants to restart with the initial set of knots the user
    !   must set iopt=0.
    !   if iopt=1 the program will continue with the set of knots found at
    !   the last call of the routine. this will save a lot of computation
    !   time if surfit is called repeatedly for different values of s.
    !   the number of knots of the spline returned and their location will
    !   depend on the value of s and on the complexity of the shape of the
    !   function underlying the data. if the computation mode iopt=1
    !   is used, the knots returned may also depend on the s-values at
    !   previous calls (if these were smaller). therefore, if after a number
    !   of trials with different s-values and iopt=1, the user can finally
    !   accept a fit as satisfactory, it may be worthwhile for him to call
    !   surfit once more with the selected value for s but now with iopt=0.
    !   indeed, surfit may then return an approximation of the same quality
    !   of fit but with fewer knots and therefore better if data reduction
    !   is also an important objective for the user.
    !   the number of knots may also depend on the upper bounds nxest and
    !   nyest. indeed, if at a certain stage in surfit the number of knots
    !   in one direction (say nx) has reached the value of its upper bound
    !   (nxest), then from that moment on all subsequent knots are added
    !   in the other (y) direction. this may indicate that the value of
    !   nxest is too small. on the other hand, it gives the user the option
    !   of limiting the number of knots the routine locates in any direction
    !   for example, by setting nxest=2*kx+2 (the lowest allowable value for
    !   nxest), the user can indicate that he wants an approximation which
    !   is a simple polynomial of degree kx in the variable x.
    !
    !  other subroutines required:
    !    fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
    !
    !  references:
    !   dierckx p. : an algorithm for surface fitting with spline functions
    !                ima j. numer. anal. 1 (1981) 267-283.
    !   dierckx p. : an algorithm for surface fitting with spline functions
    !                report tw50, dept. computer science,k.u.leuven, 1980.
    !   dierckx p. : curve and surface fitting with splines, monographs on
    !                numerical analysis, oxford university press, 1993.
    !
    !  author:
    !    p.dierckx
    !    dept. computer science, k.u. leuven
    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
    !
    !  creation date : may 1979
    !  latest update : march 1987
    !
    !  ..
    !  ..scalar arguments..
    real(kind=rp) xb,xe,yb,ye,s,eps,fp
    integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
    !  ..array arguments..
    real(kind=rp) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c((nxest-kx-1)*(nyest&
         & -ky-1)),wrk1(lwrk1),wrk2(lwrk2)
    integer iwrk(kwrk)
    !  ..local scalars..
    real(kind=rp) tol
    integer i,ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx,lby&
         & ,lco,lf,lff,lfp,lh,lq,lsx,lsy,lwest,maxit,ncest,nest,nek&
         & ,nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk
    !  ..function references..
    integer max0
    !  ..subroutine references..
    !    fpsurf
    !  ..
    !  we set up the parameters tol and maxit.
    maxit = 20
    tol = 0.1e-04
    !  before starting computations a data check is made. if the
    ! input data are invalid,control is immediately repassed to the
    ! calling program.
    ier = 10
    if(eps.le.0. .or. eps.ge.1.) go to 70
    if(kx.le.0 .or. kx.gt.5) go to 70
    kx1 = kx+1
    if(ky.le.0 .or. ky.gt.5) go to 70
    ky1 = ky+1
    kmax = max0(kx,ky)
    km1 = kmax+1
    km2 = km1+1
    if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
    if(m.lt.(kx1*ky1)) go to 70
    nminx = 2*kx1
    if(nxest.lt.nminx .or. nxest.gt.nmax) go to 70
    nminy = 2*ky1
    if(nyest.lt.nminy .or. nyest.gt.nmax) go to 70
    nest = max0(nxest,nyest)
    nxk = nxest-kx1
    nyk = nyest-ky1
    ncest = nxk*nyk
    nmx = nxest-nminx+1
    nmy = nyest-nminy+1
    nrint = nmx+nmy
    nreg = nmx*nmy
    ib1 = kx*nyk+ky1
    jb1 = ky*nxk+kx1
    ib3 = kx1*nyk+1
    if(ib1.le.jb1) go to 10
    ib1 = jb1
    ib3 = ky1*nxk+1
10  lwest = ncest*(2+ib1+ib3)+2*(nrint+nest*km2+m*km1)+ib3
    kwest = m+nreg
    if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 70
    if(xb.ge.xe .or. yb.ge.ye) go to 70
    do i=1,m
       if(w(i).le.0.) go to 70
       if(x(i).lt.xb .or. x(i).gt.xe) go to 70
       if(y(i).lt.yb .or. y(i).gt.ye) go to 70
    enddo
    if(iopt.ge.0) go to 50
    if(nx.lt.nminx .or. nx.gt.nxest) go to 70
    nxk = nx-kx1
    tx(kx1) = xb
    tx(nxk+1) = xe
    do i=kx1,nxk
       if(tx(i+1).le.tx(i)) go to 70
    enddo
    if(ny.lt.nminy .or. ny.gt.nyest) go to 70
    nyk = ny-ky1
    ty(ky1) = yb
    ty(nyk+1) = ye
    do i=ky1,nyk
       if(ty(i+1).le.ty(i)) go to 70
    enddo
    go to 60
50  if(s.lt.0.) go to 70
60  ier = 0
    !  we partition the working space and determine the spline
    ! approximation
    kn = 1
    ki = kn+m
    lq = 2
    la = lq+ncest*ib3
    lf = la+ncest*ib1
    lff = lf+ncest
    lfp = lff+ncest
    lco = lfp+nrint
    lh = lco+nrint
    lbx = lh+ib3
    nek = nest*km2
    lby = lbx+nek
    lsx = lby+nek
    lsy = lsx+m*km1
    call fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,eps&
         & ,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,nreg,nx,tx,ny&
         & ,ty,c,fp,wrk1(1),wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff)&
         & ,wrk1(la),wrk1(lq),wrk1(lbx),wrk1(lby),wrk1(lsx),wrk1(lsy)&
         & ,wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
70  return
  end subroutine surfit
  !
  !
  !
  subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy,s,nxest,nyest&
     & ,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest,nrest,nx0,tx,ny0&
     & ,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx,spy,h,index,nummer&
     & ,wrk,lwrk,ier)
    !  ..
    !  ..scalar arguments..
    real(kind=rp) xb,xe,yb,ye,s,eta,tol,fp,fp0
    integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3,nc&
         & ,intest,nrest,nx0,ny0,lwrk,ier
    !  ..array arguments..
    real(kind=rp) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest)&
         & ,coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax&
         & ,km2),by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
    integer index(nrest),nummer(m)
    !  ..local scalars..
    real(kind=rp) acc,arg,cos,dmax,fac1,fac2,fpmax,fpms,f1,f2,f3,hxi,p,pinv&
         & ,piv,p1,p2,p3,sigma,sin,sq,store,wi,x0,x1,y0,y1,zi,eps,rn&
         & ,one,con1,con9,con4,half,ten
    integer i,iband,iband1,iband3,iband4,ibb,ichang,ich1,ich3,ii,in&
         & ,irot,iter,i1,i2,i3,j,jrot,jxy,j1,kx,kx1,kx2,ky,ky1,ky2,l&
         & ,la,lf,lh,lwest,lx,ly,l1,l2,n,ncof,nk1x,nk1y,nminx,nminy&
         & ,nreg,nrint,num,num1,nx,nxe,nxx,ny,nye,nyy,n1,rank
    !  ..local arrays..
    real(kind=rp) hx(6),hy(6)
    !  ..function references..
    real(kind=rp) abs,sqrt
    integer min0
    !  ..subroutine references..
    !    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota
    !  ..
    !  set constants
    one = 0.1e+01_rp
    con1 = 0.1e0_rp
    con9 = 0.9e0_rp
    con4 = 0.4e-01_rp
    half = 0.5e0_rp
    ten = 0.1e+02_rp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! part 1: determination of the number of knots and their position.     !
    ! ****************************************************************     !
    ! given a set of knots we compute the least-squares spline sinf(x,y),  !
    ! and the corresponding weighted sum of squared residuals fp=f(p=inf). !
    ! if iopt=-1  sinf(x,y) is the requested approximation.                !
    ! if iopt=0 or iopt=1 we check whether we can accept the knots:        !
    !   if fp <=s we will continue with the current set of knots.          !
    !   if fp > s we will increase the number of knots and compute the     !
    !      corresponding least-squares spline until finally  fp<=s.        !
    ! the initial choice of knots depends on the value of s and iopt.      !
    !   if iopt=0 we first compute the least-squares polynomial of degree  !
    !     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        !
    !     fp0=f(0) denotes the corresponding weighted sum of squared       !
    !     residuals                                                        !
    !   if iopt=1 we start with the knots found at the last call of the    !
    !     routine, except for the case that s>=fp0; then we can compute    !
    !     the least-squares polynomial directly.                           !
    ! eventually the independent variables x and y (and the corresponding  !
    ! parameters) will be switched if this can reduce the bandwidth of the !
    ! system to be solved.                                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  ichang denotes whether(1) or not(-1) the directions have been inter-
    !  changed.
    ichang = -1
    x0 = xb
    x1 = xe
    y0 = yb
    y1 = ye
    kx = kxx
    ky = kyy
    kx1 = kx+1
    ky1 = ky+1
    nxe = nxest
    nye = nyest
    eps = sqrt(eta)
    if(iopt.lt.0) go to 20
    !  calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol*s
    if(iopt.eq.0) go to 10
    if(fp0.gt.s) go to 20
    !  initialization for the least-squares polynomial.
10  nminx = 2*kx1
    nminy = 2*ky1
    nx = nminx
    ny = nminy
    ier = -2
    go to 30
20  nx = nx0
    ny = ny0
    !  main loop for the different sets of knots. m is a save upper bound
    !  for the number of trials.
30  do 420 iter=1,m
       !  find the position of the additional knots which are needed for the
       !  b-spline representation of s(x,y).
       l = nx
       do i=1,kx1
          tx(i) = x0
          tx(l) = x1
          l = l-1
       enddo
       l = ny
       do i=1,ky1
          ty(i) = y0
          ty(l) = y1
          l = l-1
       enddo
       !  find nrint, the total number of knot intervals and nreg, the number
       !  of panels in which the approximation domain is subdivided by the
       !  intersection of knots.
       nxx = nx-2*kx1+1
       nyy = ny-2*ky1+1
       nrint = nxx+nyy
       nreg = nxx*nyy
       !  find the bandwidth of the observation matrix a.
       !  if necessary, interchange the variables x and y, in order to obtain
       !  a minimal bandwidth.
       iband1 = kx*(ny-ky1)+ky
       l = ky*(nx-kx1)+kx
       if(iband1.le.l) go to 130
       iband1 = l
       ichang = -ichang
       do i=1,m
          store = x(i)
          x(i) = y(i)
          y(i) = store
       enddo
       store = x0
       x0 = y0
       y0 = store
       store = x1
       x1 = y1
       y1 = store
       n = min0(nx,ny)
       do i=1,n
          store = tx(i)
          tx(i) = ty(i)
          ty(i) = store
       enddo
       n1 = n+1
       if(nx-ny) 80,120,100
80     do i=n1,ny
          tx(i) = ty(i)
       enddo
       go to 120
100    do i=n1,nx
          ty(i) = tx(i)
       enddo
120    l = nx
       nx = ny
       ny = l
       l = nxe
       nxe = nye
       nye = l
       l = nxx
       nxx = nyy
       nyy = l
       l = kx
       kx = ky
       ky = l
       kx1 = kx+1
       ky1 = ky+1
130    iband = iband1+1
       !  arrange the data points according to the panel they belong to.
       call fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
       !  find ncof, the number of b-spline coefficients.
       nk1x = nx-kx1
       nk1y = ny-ky1
       ncof = nk1x*nk1y
       !  initialize the observation matrix a.
       do i=1,ncof
          f(i) = 0.
          do j=1,iband
             a(i,j) = 0.
          enddo
       enddo
       !  initialize the sum of squared residuals.
       fp = 0.
       !  fetch the data points in the new order. main loop for the
       !  different panels.
       do 250 num=1,nreg
          !  fix certain constants for the current panel; jrot records the column
          !  number of the first non-zero element in a row of the observation
          !  matrix according to a data point of the panel.
          num1 = num-1
          lx = num1/nyy
          l1 = lx+kx1
          ly = num1-lx*nyy
          l2 = ly+ky1
          jrot = lx*nk1y+ly
          !  test whether there are still data points in the panel.
          in = index(num)
150       if(in.eq.0) go to 250
          !  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
          !  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
          call fpbspl(tx,nx,kx,x(in),l1,hx)
          !  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
          call fpbspl(ty,ny,ky,y(in),l2,hy)
          !  store the value of these b-splines in spx and spy respectively.
          do i=1,kx1
             spx(in,i) = hx(i)
          enddo
          do i=1,ky1
             spy(in,i) = hy(i)
          enddo
          !  initialize the new row of observation matrix.
          do i=1,iband
             h(i) = 0.
          enddo
          !  calculate the non-zero elements of the new row by making the cross
          !  products of the non-zero b-splines in x- and y-direction.
          i1 = 0
          do i=1,kx1
             hxi = hx(i)
             j1 = i1
             do j=1,ky1
                j1 = j1+1
                h(j1) = hxi*hy(j)*wi
             enddo
             i1 = i1+nk1y
          enddo
          !  rotate the row into triangle by givens transformations .
          irot = jrot
          do 220 i=1,iband
             irot = irot+1
             piv = h(i)
             if(piv.eq.0.) go to 220
             !  calculate the parameters of the givens transformation.
             call fpgivs(piv,a(irot,1),cos,sin)
             !  apply that transformation to the right hand side.
             call fprota(cos,sin,zi,f(irot))
             if(i.eq.iband) go to 230
             !  apply that transformation to the left hand side.
             i2 = 1
             i3 = i+1
             do j=i3,iband
                i2 = i2+1
                call fprota(cos,sin,h(j),a(irot,i2))
             enddo
220       continue
          !  add the contribution of the row to the sum of squares of residual
          !  right hand sides.
230       fp = fp+zi**2
          !  find the number of the next data point in the panel.
240       in = nummer(in)
          go to 150
250    continue
       !  find dmax, the maximum value for the diagonal elements in the reduced
       !  triangle.
       dmax = 0.
       do 260 i=1,ncof
          if(a(i,1).le.dmax) go to 260
          dmax = a(i,1)
 260    continue
       !  check whether the observation matrix is rank deficient.
       sigma = eps*dmax
       do i=1,ncof
          if(a(i,1).le.sigma) go to 280
       enddo
       !  backward substitution in case of full rank.
       call fpback(a,f,ncof,iband,c,nc)
       rank = ncof
       do i=1,ncof
          q(i,1) = a(i,1)/dmax
       enddo
       go to 300
       !  in case of rank deficiency, find the minimum norm solution.
       !  check whether there is sufficient working space
280    lwest = ncof*iband+ncof+iband
       if(lwrk.lt.lwest) go to 780
       do i=1,ncof
          ff(i) = f(i)
          do j=1,iband
             q(i,j) = a(i,j)
          enddo
       enddo
       lf =1
       lh = lf+ncof
       la = lh+iband
       call fprank(q,ff,ncof,iband,nc,sigma,c,sq,rank,wrk(la),wrk(lf),wrk(lh))
       do i=1,ncof
          q(i,1) = q(i,1)/dmax
       enddo
       !  add to the sum of squared residuals, the contribution of reducing
       !  the rank.
       fp = fp+sq
300    if(ier.eq.(-2)) fp0 = fp
       !  test whether the least-squares spline is an acceptable solution.
       if(iopt.lt.0) go to 820
       fpms = fp-s
       if(abs(fpms).le.acc) if(fp) 815,815,820
       !  test whether we can accept the choice of knots.
       if(fpms.lt.0.) go to 430
       !  test whether we cannot further increase the number of knots.
       if(ncof.gt.m) go to 790
       ier = 0
       !  search where to add a new knot.
       !  find for each interval the sum of squared residuals fpint for the
       !  data points having the coordinate belonging to that knot interval.
       !  calculate also coord which is the same sum, weighted by the position
       !  of the data points considered.
310    do i=1,nrint
          fpint(i) = 0.0_rp
          coord(i) = 0.0_rp
       enddo
       do 360 num=1,nreg
          num1 = num-1
          lx = num1/nyy
          l1 = lx+1
          ly = num1-lx*nyy
          l2 = ly+1+nxx
          jrot = lx*nk1y+ly
          in = index(num)
 330      if(in.eq.0) go to 360
          store = 0.0_rp
          i1 = jrot
          do i=1,kx1
             hxi = spx(in,i)
             j1 = i1
            do j=1,ky1
               j1 = j1+1
               store = store+hxi*spy(in,j)*c(j1)
            enddo
            i1 = i1+nk1y
         enddo
         store = (w(in)*(z(in)-store))**2
         fpint(l1) = fpint(l1)+store
         coord(l1) = coord(l1)+store*x(in)
         fpint(l2) = fpint(l2)+store
         coord(l2) = coord(l2)+store*y(in)
         in = nummer(in)
         go to 330
360    continue
       !  find the interval for which fpint is maximal on the condition that
       !  there still can be added a knot.
 370   l = 0
       fpmax = 0.0_rp
       l1 = 1
       l2 = nrint
       if(nx.eq.nxe) l1 = nxx+1
       if(ny.eq.nye) l2 = nxx
       if(l1.gt.l2) go to 810
       do 380 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 380
          l = i
          fpmax = fpint(i)
 380   continue
       !  test whether we cannot further increase the number of knots.
       if(l.eq.0) go to 785
       !  calculate the position of the new knot.
       arg = coord(l)/fpint(l)
       !  test in what direction the new knot is going to be added.
       if(l.gt.nxx) go to 400
       !  addition in the x-direction.
       jxy = l+kx1
       fpint(l) = 0.
       fac1 = tx(jxy)-arg
       fac2 = arg-tx(jxy-1)
       if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
       j = nx
       do i=jxy,nx
          tx(j+1) = tx(j)
          j = j-1
       enddo
       tx(jxy) = arg
       nx = nx+1
       go to 420
       !  addition in the y-direction.
400    jxy = l+ky1-nxx
       fpint(l) = 0.
       fac1 = ty(jxy)-arg
       fac2 = arg-ty(jxy-1)
       if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
       j = ny
       do i=jxy,ny
          ty(j+1) = ty(j)
          j = j-1
       enddo
       ty(jxy) = arg
       ny = ny+1
       !  restart the computations with the new set of knots.
420 continue
    !  test whether the least-squares polynomial is a solution of our
    !  approximation problem.
430 if(ier.eq.(-2)) go to 830
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! part 2: determination of the smoothing spline sp(x,y)                !
    ! *****************************************************                !
    ! we have determined the number of knots and their position. we now    !
    ! compute the b-spline coefficients of the smoothing spline sp(x,y).   !
    ! the observation matrix a is extended by the rows of a matrix,        !
    ! expressing that sp(x,y) must be a polynomial of degree kx in x and   !
    ! ky in y. the corresponding weights of these additional rows are set  !
    ! to 1./p.  iteratively we than have to determine the value of p       !
    ! such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           !
    ! we already know that the least-squares polynomial corresponds to     !
    ! p=0  and that the least-squares spline corresponds to p=infinity.    !
    ! the iteration process which is proposed here makes use of rational   !
    ! interpolation. since f(p) is a convex and strictly decreasing        !
    ! function of p, it can be approximated by a rational function r(p)=   !
    ! (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values !
    ! of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the !
    ! new value of p such that r(p)=s. convergence is guaranteed by taking !
    ! f1 > 0 and f3 < 0.                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kx2 = kx1+1
    !  test whether there are interior knots in the x-direction.
    if(nk1x.eq.kx1) go to 440
    !  evaluate the discotinuity jumps of the kx-th order derivative of
    !  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
    call fpdisc(tx,nx,kx2,bx,nmax)
440 ky2 = ky1 + 1
    !  test whether there are interior knots in the y-direction.
    if(nk1y.eq.ky1) go to 450
    !  evaluate the discontinuity jumps of the ky-th order derivative of
    !  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
    call fpdisc(ty,ny,ky2,by,nmax)
    !  initial value for p.
450 p1 = 0.0_rp
    f1 = fp0-s
    p3 = -one
    f3 = fpms
    p = 0.0_rp
    do i=1,ncof
       p = p+a(i,1)
    enddo
    rn = ncof
    p = rn/p
    !  find the bandwidth of the extended observation matrix.
    iband3 = kx1*nk1y
    iband4 = iband3 +1
    ich1 = 0
    ich3 = 0
    !  iteration process to find the root of f(p)=s.
    do 770 iter=1,maxit
       pinv = one/p
       !  store the triangularized observation matrix into q.
       do i=1,ncof
          ff(i) = f(i)
          do j=1,iband
             q(i,j) = a(i,j)
          enddo
          ibb = iband+1
          do j=ibb,iband4
             q(i,j) = 0.0_rp
          enddo
       enddo
       if(nk1y.eq.ky1) go to 560
       !  extend the observation matrix with the rows of a matrix, expressing
       !  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
       do 550 i=ky2,nk1y
          ii = i-ky1
          do 550 j=1,nk1x
             !  initialize the new row.
             do l=1,iband
                h(l) = 0.0_rp
             enddo
             !  fill in the non-zero elements of the row. jrot records the column
             !  number of the first non-zero element in the row.
             do l=1,ky2
                h(l) = by(ii,l)*pinv
             enddo
             zi = 0.0_rp
             jrot = (j-1)*nk1y+ii
             !  rotate the new row into triangle by givens transformations without
             !  square roots.
             do irot=jrot,ncof
                piv = h(1)
                i2 = min0(iband1,ncof-irot)
                if(piv.eq.0.) if(i2) 550,550,520
                !  calculate the parameters of the givens transformation.
                call fpgivs(piv,q(irot,1),cos,sin)
                !  apply that givens transformation to the right hand side.
                call fprota(cos,sin,zi,ff(irot))
                if(i2.eq.0) go to 550
                !  apply that givens transformation to the left hand side.
                do l=1,i2
                   l1 = l+1
                   call fprota(cos,sin,h(l1),q(irot,l1))
                enddo
520             do l=1,i2
                   h(l) = h(l+1)
                enddo
                h(i2+1) = 0.0_rp
             enddo
550    continue
560    if(nk1x.eq.kx1) go to 640
       !  extend the observation matrix with the rows of a matrix expressing
       !  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
       do 630 i=kx2,nk1x
          ii = i-kx1
          do 630 j=1,nk1y
             !  initialize the new row
             do l=1,iband4
                h(l) = 0.0_rp
             enddo
             !  fill in the non-zero elements of the row. jrot records the column
             !  number of the first non-zero element in the row.
             j1 = 1
             do l=1,kx2
                h(j1) = bx(ii,l)*pinv
                j1 = j1+nk1y
             enddo
             zi = 0.0_rp
             jrot = (i-kx2)*nk1y+j
             !  rotate the new row into triangle by givens transformations .
             do irot=jrot,ncof
                piv = h(1)
                i2 = min0(iband3,ncof-irot)
                if(piv.eq.0.) if(i2) 630,630,600
                !  calculate the parameters of the givens transformation.
                call fpgivs(piv,q(irot,1),cos,sin)
                !  apply that givens transformation to the right hand side.
                call fprota(cos,sin,zi,ff(irot))
                if(i2.eq.0) go to 630
                !  apply that givens transformation to the left hand side.
                do l=1,i2
                   l1 = l+1
                   call fprota(cos,sin,h(l1),q(irot,l1))
                enddo
600             do l=1,i2
                   h(l) = h(l+1)
                enddo
                h(i2+1) = 0.0_rp
             enddo
630    continue
       !  find dmax, the maximum value for the diagonal elements in the
       !  reduced triangle.
640    dmax = 0.0_rp
       do 650 i=1,ncof
          if(q(i,1).le.dmax) go to 650
          dmax = q(i,1)
650    continue
       !  check whether the matrix is rank deficient.
       sigma = eps*dmax
       do i=1,ncof
          if(q(i,1).le.sigma) go to 670
       enddo
       !  backward substitution in case of full rank.
       call fpback(q,ff,ncof,iband4,c,nc)
       rank = ncof
       go to 675
       !  in case of rank deficiency, find the minimum norm solution.
670    lwest = ncof*iband4+ncof+iband4
       if(lwrk.lt.lwest) go to 780
       lf = 1
       lh = lf+ncof
       la = lh+iband4
       call fprank(q,ff,ncof,iband4,nc,sigma,c,sq,rank,wrk(la),wrk(lf),wrk(lh))
675    do i=1,ncof
          q(i,1) = q(i,1)/dmax
       enddo
       !  compute f(p).
       fp = 0.0_rp
       do 720 num = 1,nreg
          num1 = num-1
          lx = num1/nyy
          ly = num1-lx*nyy
          jrot = lx*nk1y+ly
          in = index(num)
690       if(in.eq.0) go to 720
          store = 0.0_rp
          i1 = jrot
          do i=1,kx1
             hxi = spx(in,i)
             j1 = i1
             do j=1,ky1
                j1 = j1+1
                store = store+hxi*spy(in,j)*c(j1)
             enddo
             i1 = i1+nk1y
          enddo
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 690
720    continue
       !  test whether the approximation sp(x,y) is an acceptable solution.
       fpms = fp-s
       if(abs(fpms).le.acc) go to 820
       !  test whether the maximum allowable number of iterations has been
       !  reached.
       if(iter.eq.maxit) go to 795
       !  carry out one more step of the iteration process.
       p2 = p
       f2 = fpms
       if(ich3.ne.0) go to 740
       if((f2-f3).gt.acc) go to 730
       !  our initial choice of p is too large.
       p3 = p2
       f3 = f2
       p = p*con4
       if(p.le.p1) p = p1*con9 + p2*con1
       go to 770
730    if(f2.lt.0.) ich3 = 1
740    if(ich1.ne.0) go to 760
       if((f1-f2).gt.acc) go to 750
       !  our initial choice of p is too small
       p1 = p2
       f1 = f2
       p = p/con4
       if(p3.lt.0.) go to 770
       if(p.ge.p3) p = p2*con1 + p3*con9
       go to 770
750    if(f2.gt.0.) ich1 = 1
       !  test whether the iteration process proceeds as theoretically
       !  expected.
760    if(f2.ge.f1 .or. f2.le.f3) go to 800
       !  find the new value of p.
       p = fprati(p1,f1,p2,f2,p3,f3)
770 continue
    !  error codes and messages.
780 ier = lwest
    go to 830
785 ier = 5
    go to 830
790 ier = 4
    go to 830
795 ier = 3
    go to 830
800 ier = 2
    go to 830
810 ier = 1
    go to 830
815 ier = -1
    fp = 0.
820 if(ncof.ne.rank) ier = -rank
    !  test whether x and y are in the original order.
830 if(ichang.lt.0) go to 930
    !  if not, interchange x and y once more.
    l1 = 1
    do i=1,nk1x
       l2 = i
       do j=1,nk1y
          f(l2) = c(l1)
          l1 = l1+1
          l2 = l2+nk1x
       enddo
    enddo
    do i=1,ncof
       c(i) = f(i)
    enddo
    do i=1,m
       store = x(i)
       x(i) = y(i)
       y(i) = store
    enddo
    n = min0(nx,ny)
    do i=1,n
       store = tx(i)
       tx(i) = ty(i)
       ty(i) = store
    enddo
    n1 = n+1
    if(nx-ny) 880,920,900
880 do i=n1,ny
       tx(i) = ty(i)
    enddo
    go to 920
900 do i=n1,nx
       ty(i) = tx(i)
    enddo
920 l = nx
    nx = ny
    ny = l
930 if(iopt.lt.0) go to 940
    nx0 = nx
    ny0 = ny
940 return
  end subroutine fpsurf
  !
  !
  !
  subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
    !  subroutine fprank finds the minimum norm solution of a least-
    !  squares problem in case of rank deficiency.
    !
    !  input parameters:
    !    a : array, which contains the non-zero elements of the observation
    !        matrix after triangularization by givens transformations.
    !    f : array, which contains the transformed right hand side.
    !    n : integer,wich contains the dimension of a.
    !    m : integer, which denotes the bandwidth of a.
    !  tol : real value, giving a threshold to determine the rank of a.
    !
    !  output parameters:
    !    c : array, which contains the minimum norm solution.
    !   sq : real value, giving the contribution of reducing the rank
    !        to the sum of squared residuals.
    ! rank : integer, which contains the rank of matrix a.
    !
    !  ..scalar arguments..
    integer n,m,na,rank
    real(kind=rp) tol,sq
    !  ..array arguments..
    real(kind=rp) a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
    !  ..local scalars..
    integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
    real(kind=rp) cos,fac,piv,sin,yi
    double precision store,stor1,stor2,stor3
    !  ..function references..
    integer min0
    !  ..subroutine references..
    !    fpgivs,fprota
    !  ..
    m1 = m-1
    !  the rank deficiency nl is considered to be the number of sufficient
    !  small diagonal elements of a.
    nl = 0
    sq = 0.0_rp
    do 90 i=1,n
       if(a(i,1).gt.tol) go to 90
       !  if a sufficient small diagonal element is found, we put it to
       !  zero. the remainder of the row corresponding to that zero diagonal
       !  element is then rotated into triangle by givens rotations .
       !  the rank deficiency is increased by one.
       nl = nl+1
       if(i.eq.n) go to 90
       yi = f(i)
       do j=1,m1
          h(j) = a(i,j+1)
       enddo
       h(m) = 0.0_rp
       i1 = i+1
       do ii=i1,n
          i2 = min0(n-ii,m1)
          piv = h(1)
          if(piv.eq.0.) go to 30
          call fpgivs(piv,a(ii,1),cos,sin)
          call fprota(cos,sin,yi,f(ii))
          if(i2.eq.0) go to 70
          do j=1,i2
             j1 = j+1
             call fprota(cos,sin,h(j1),a(ii,j1))
             h(j) = h(j1)
          enddo
          go to 50
30        if(i2.eq.0) go to 70
          do j=1,i2
             h(j) = h(j+1)
          enddo
50        h(i2+1) = 0.0_rp
       enddo
       !  add to the sum of squared residuals the contribution of deleting
       !  the row with small diagonal element.
70     sq = sq+yi**2
90  continue
    !  rank denotes the rank of a.
    rank = n-nl
    !  let b denote the (rank*n) upper trapezoidal matrix which can be
    !  obtained from the (n*n) upper triangular matrix a by deleting
    !  the rows and interchanging the columns corresponding to a zero
    !  diagonal element. if this matrix is factorized using givens
    !  transformations as  b = (r) (u)  where
    !    r is a (rank*rank) upper triangular matrix,
    !    u is a (rank*n) orthonormal matrix
    !  then the minimal least-squares solution c is given by c = b' v,
    !  where v is the solution of the system  (r) (r)' v = g  and
    !  g denotes the vector obtained from the old right hand side f, by
    !  removing the elements corresponding to a zero diagonal element of a.
    !  initialization.
    do i=1,rank
        do j=1,m
          aa(i,j) = 0.0_rp
       enddo
    enddo
    !  form in aa the upper triangular matrix obtained from a by
    !  removing rows and columns with zero diagonal elements. form in ff
    !  the new right hand side by removing the elements of the old right
    !  hand side corresponding to a deleted row.
    ii = 0
    do 120 i=1,n
       if(a(i,1).le.tol) go to 120
       ii = ii+1
       ff(ii) = f(i)
       aa(ii,1) = a(i,1)
       jj = ii
       kk = 1
       j = i
       j1 = min0(j-1,m1)
       if(j1.eq.0) go to 120
       do 110 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 110
          kk = kk+1
          jj = jj-1
          aa(jj,kk) = a(j,k+1)
110    continue
120 continue
    !  form successively in h the columns of a with a zero diagonal element.
    ii = 0
    do 200 i=1,n
        ii = ii+1
        if(a(i,1).gt.tol) go to 200
        ii = ii-1
        if(ii.eq.0) go to 200
        jj = 1
        j = i
        j1 = min0(j-1,m1)
        do 130 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 130
          h(jj) = a(j,k+1)
          jj = jj+1
 130    continue
        do 140 kk=jj,m
          h(kk) = 0.0_rp
 140    continue
!  rotate this column into aa by givens transformations.
        jj = ii
        do i1=1,ii
          j1 = min0(jj-1,m1)
          piv = h(1)
          if(piv.ne.0.) go to 160
          if(j1.eq.0) go to 200
          do j2=1,j1
             j3 = j2+1
             h(j2) = h(j3)
          enddo
          go to 180
160       call fpgivs(piv,aa(jj,1),cos,sin)
          if(j1.eq.0) go to 200
          kk = jj
          do j2=1,j1
             j3 = j2+1
             kk = kk-1
             call fprota(cos,sin,h(j3),aa(kk,j3))
             h(j2) = h(j3)
          enddo
180       jj = jj-1
          h(j3) = 0.0_rp
       enddo
200    continue
       !  solve the system (aa) (f1) = ff
       ff(rank) = ff(rank)/aa(rank,1)
       i = rank-1
       if(i.eq.0) go to 230
       do j=2,rank
        store = ff(i)
        i1 = min0(j-1,m1)
        k = i
        do ii=1,i1
           k = k+1
           stor1 = ff(k)
           stor2 = aa(i,ii+1)
           store = store-stor1*stor2
        enddo
        stor1 = aa(i,1)
        ff(i) = store/stor1
        i = i-1
     enddo
     !  solve the system  (aa)' (f2) = f1
230  ff(1) = ff(1)/aa(1,1)
     if(rank.eq.1) go to 260
     do j=2,rank
        store = ff(j)
        i1 = min0(j-1,m1)
        k = j
        do ii=1,i1
           k = k-1
           stor1 = ff(k)
           stor2 = aa(k,ii+1)
           store = store-stor1*stor2
        enddo
        stor1 = aa(j,1)
        ff(j) = store/stor1
     enddo
     !  premultiply f2 by the transpoze of a.
260  k = 0
     do i=1,n
        store = 0.0_rp
        if(a(i,1).gt.tol) k = k+1
        j1 = min0(i,m)
        kk = k
        ij = i+1
        do 270 j=1,j1
           ij = ij-1
           if(a(ij,1).le.tol) go to 270
           stor1 = a(ij,j)
           stor2 = ff(kk)
           store = store+stor1*stor2
           kk = kk-1
270     continue
        c(i) = store
     enddo
     !  add to the sum of squared residuals the contribution of putting
     !  to zero the small diagonal elements of matrix (a).
     stor3 = 0.0_rp
     do 310 i=1,n
        if(a(i,1).gt.tol) go to 310
        store = f(i)
        i1 = min0(n-i,m1)
        if(i1.eq.0) go to 300
        do j=1,i1
           ij = i+j
           stor1 = c(ij)
           stor2 = a(i,j+1)
           store = store-stor1*stor2
        enddo
300     fac = a(i,1)*c(i)
        stor1 = a(i,1)
        stor2 = c(i)
        stor1 = stor1*stor2
        stor3 = stor3+stor1*(stor1-store-store)
310  continue
     fac = stor3
     sq = sq+fac
     return
  end subroutine fprank
  !
  !
  !
  subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
    !  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
    !  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
    !  to. for each panel a stack is constructed  containing the numbers
    !  of data points lying inside; index(j),j=1,2,...,nreg points to the
    !  first data point in the jth panel while nummer(i),i=1,2,...,m gives
    !  the number of the next data point in the panel.
    !  ..
    !  ..scalar arguments..
    integer m,kx,ky,nx,ny,nreg
    !  ..array arguments..
    real(kind=rp) x(m),y(m),tx(nx),ty(ny)
    integer nummer(m),index(nreg)
    !  ..local scalars..
    real(kind=rp) xi,yi
    integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
    !  ..
    kx1 = kx+1
    ky1 = ky+1
    nk1x = nx-kx1
    nk1y = ny-ky1
    nyy = nk1y-ky
    do i=1,nreg
       index(i) = 0
    enddo
    do im=1,m
       xi = x(im)
       yi = y(im)
       l = kx1
       l1 = l+1
20     if(xi.lt.tx(l1) .or. l.eq.nk1x) go to 30
       l = l1
       l1 = l+1
       go to 20
30     k = ky1
       k1 = k+1
40     if(yi.lt.ty(k1) .or. k.eq.nk1y) go to 50
       k = k1
       k1 = k+1
       go to 40
50     num = (l-kx1)*nyy+k-ky
       nummer(im) = index(num)
       index(num) = im
    enddo
    return
  end subroutine fporde
  !
  !
  !
!  subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk&
!       & ,iwrk,kwrk,ier)
!    !  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
!    !  ,my the partial derivative ( order nux,nuy) of a bivariate spline
!    !  s(x,y) of degrees kx and ky, given in the b-spline representation.
!    !
!    !  calling sequence:
!    !     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
!    !    * iwrk,kwrk,ier)
!    !
!    !  input parameters:
!    !   tx    : real array, length nx, which contains the position of the
!    !           knots in the x-direction.
!    !   nx    : integer, giving the total number of knots in the x-direction
!    !   ty    : real array, length ny, which contains the position of the
!    !           knots in the y-direction.
!    !   ny    : integer, giving the total number of knots in the y-direction
!    !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!    !           b-spline coefficients.
!    !   kx,ky : integer values, giving the degrees of the spline.
!    !   nux   : integer values, specifying the order of the partial
!    !   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
!    !   x     : real array of dimension (mx).
!    !           before entry x(i) must be set to the x co-ordinate of the
!    !           i-th grid point along the x-axis.
!    !           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
!    !   mx    : on entry mx must specify the number of grid points along
!    !           the x-axis. mx >=1.
!    !   y     : real array of dimension (my).
!    !           before entry y(j) must be set to the y co-ordinate of the
!    !           j-th grid point along the y-axis.
!    !           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
!    !   my    : on entry my must specify the number of grid points along
!    !           the y-axis. my >=1.
!    !   wrk   : real array of dimension lwrk. used as workspace.
!    !   lwrk  : integer, specifying the dimension of wrk.
!    !           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
!    !   iwrk  : integer array of dimension kwrk. used as workspace.
!    !   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
!    !
!    !  output parameters:
!    !   z     : real array of dimension (mx*my).
!    !           on succesful exit z(my*(i-1)+j) contains the value of the
!    !           specified partial derivative of s(x,y) at the point
!    !           (x(i),y(j)),i=1,...,mx;j=1,...,my.
!    !   ier   : integer error flag
!    !    ier=0 : normal return
!    !    ier=10: invalid input data (see restrictions)
!    !
!    !  restrictions:
!    !   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my
!    !   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
!    !   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!    !   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!    !
!    !  other subroutines required:
!    !    fpbisp,fpbspl
!    !
!    !  references :
!    !    de boor c : on calculating with b-splines, j. approximation theory
!    !                6 (1972) 50-62.
!    !   dierckx p. : curve and surface fitting with splines, monographs on
!    !                numerical analysis, oxford university press, 1993.
!    !
!    !  author :
!    !    p.dierckx
!    !    dept. computer science, k.u.leuven
!    !    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!    !
!    !  latest update : march 1989
!    !
!    !  ..scalar arguments..
!    integer nx,ny,kx,ky,nux,nuy,mx,my,lwrk,kwrk,ier
!    !  ..array arguments..
!    integer iwrk(kwrk)
!    real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),wrk(lwrk)
!    !  ..local scalars..
!    integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,m,m0,m1,nc&
!         & ,nkx1,nky1,nxx,nyy
!    real ak,fac
!    !  ..
!    !  before starting computations a data check is made. if the input data
!    !  are invalid control is immediately repassed to the calling program.
!    ier = 10
!    kx1 = kx+1
!    ky1 = ky+1
!    nkx1 = nx-kx1
!    nky1 = ny-ky1
!    nc = nkx1*nky1
!    if(nux.lt.0 .or. nux.ge.kx) go to 400
!    if(nuy.lt.0 .or. nuy.ge.ky) go to 400
!    lwest = nc +(kx1-nux)*mx+(ky1-nuy)*my
!    if(lwrk.lt.lwest) go to 400
!    if(kwrk.lt.(mx+my)) go to 400
!    if(mx-1) 400,30,10
!10  do i=2,mx
!       if(x(i).lt.x(i-1)) go to 400
!    enddo
!30  if(my-1) 400,60,40
!40  do i=2,my
!       if(y(i).lt.y(i-1)) go to 400
!    enddo
!60  ier = 0
!    nxx = nkx1
!    nyy = nky1
!    kkx = kx
!    kky = ky
!    !  the partial derivative of order (nux,nuy) of a bivariate spline of
!    !  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
!    !  we calculate the b-spline coefficients of this spline
!    do i=1,nc
!       wrk(i) = c(i)
!    enddo
!    if(nux.eq.0) go to 200
!    lx = 1
!    do j=1,nux
!       ak = kkx
!       nxx = nxx-1
!       l1 = lx
!       m0 = 1
!       do 90 i=1,nxx
!          l1 = l1+1
!          l2 = l1+kkx
!          fac = tx(l2)-tx(l1)
!          if(fac.le.0.) go to 90
!          do m=1,nyy
!             m1 = m0+nyy
!             wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
!             m0  = m0+1
!          enddo
!90     continue
!       lx = lx+1
!       kkx = kkx-1
!    enddo
!200 if(nuy.eq.0) go to 300
!    ly = 1
!    do j=1,nuy
!       ak = kky
!       nyy = nyy-1
!       l1 = ly
!       do 220 i=1,nyy
!          l1 = l1+1
!          l2 = l1+kky
!          fac = ty(l2)-ty(l1)
!          if(fac.le.0.) go to 220
!          m0 = i
!          do m=1,nxx
!             m1 = m0+1
!             wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
!             m0  = m0+nky1
!          enddo
!220    continue
!       ly = ly+1
!       kky = kky-1
!    enddo
!    !
!    m0 = nyy
!    m1 = nky1
!    do m=2,nxx
!       do i=1,nyy
!          m0 = m0+1
!          m1 = m1+1
!          wrk(m0) = wrk(m1)
!       enddo
!       m1 = m1+nuy
!    enddo
!    !  we partition the working space and evaluate the partial derivative
!300 iwx = 1+nxx*nyy
!    iwy = iwx+mx*(kx1-nux)
!    call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky,x&
!         & ,mx,y,my,z,wrk(iwx),wrk(iwy),iwrk(1),iwrk(mx+1))
!400 return
!  end subroutine parder


end module interpolation
