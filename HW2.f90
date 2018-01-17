!****************************************************************************
!  HW2.f90 
!****************************************************************************

	program HW2
  	  ! Variables
	  integer nmax
	  real*8 pmax,tol
	  complex*16 p0
	  open(unit=7,file='hw2_input.txt',status='old',action='read')
	  open(unit=11,file='hw2_output.txt',action='write',position='append')
	  read(7,*)nmax,pmax,tol,p0

	  call hw2_newton(nmax,pmax,p0,tol)
	end program HW2
!****************************************************************************
!  subroutine hw2_newton
!****************************************************************************
	subroutine hw2_newton(nmax,pmax,p0,tol)
	  integer nmax,o,k
	  real*8 pmax,tol
	  complex*16 p0,p,f,fp
	  o=11
	  write(o,"(50('-'))")
10	  write(6,*)'Input k to choose function'
	  write(6,*)'1 f(x)=x*x*x+4.d0*x*x-10.d0'
	  write(6,*)'2 f(x)=dcos(x)-x'
	  read(5,*)k
	  if(k == 1)then
		write(o,*)'  f(x)=x*x*x+4.d0*x*x-10.d0'
	  else if(k == 2)then
		write(o,*)'  f(x)=dcos(x)-x'
	  else 
		write(6,*)'goto'
		goto 10
	  end if

  	  write(o,"(50('-'))")
	  write(o,"(2X,'namx = ',I4)")nmax
	  write(o,"(2X,'pamx =',F8.3)")pmax
	  write(o,"(2X,'P0 = ','(',E15.9,1X,E15.9,')')")p0
	  write(o,"(2X,'TOL = ',F10.8)")tol
	  write(o,"(50('-'))")

	  ! Body of HW2
	  call newton(k, nmax,p0,pmax,tol,p)
	  write(o,"(50('-'))")
	  write(o,*)'Answer is',p
	end subroutine hw2_newton
!****************************************************************************
!  subroutine newton
!****************************************************************************
	subroutine newton(k, nmax,p0,pmax,tol,p)
	  integer nmax,o,k
	  real*8 pmax,tol
	  complex*16 p0,p,f,fp
      o=11
	  write(o,"(2X,'n',6X,'pn')")
	  write(o,"(2X,I2,2X,'(',E15.9,1X,E15.9,')')")0,p0
	  do i=1, nmax
	    p=p0-f(k,p0)/fp(k,p0)
		write(o,"(2X,I2,2X,'(',E15.9,1X,E15.9,')')")i,p
	    if(cdabs((p-p0)/p)<tol)then
	      return
	    elseif(cdabs(p)>pmax)then
		  write(o,*)'Divergent'
	      stop
	    else
	      p0=p
	    end if
      end do
      write(o,*)'Over Nmax'
	  return
	end subroutine newton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function f(k,x)
	  complex*16 x,f
	  if(k==1)then
		f = x*x*x+4.d0*x*x-10.d0
	  elseif(k==2)then
		f=cdcos(x)-x
	  else
		stop
	  end if
	end function f

	function fp(k,x)
	  complex*16 x,fp
	  if(k==1)then
		fp = 3.d0*x*x+8.d0*x
	  elseif(k==2)then
		fp=-cdsin(x)-1.d0
	  else
		stop
	  end if
	end function fp
