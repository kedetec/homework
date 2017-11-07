!  HW2.f90 
!****************************************************************************

	program HW2

	implicit none

	! Variables
	integer nmax,k
	real*8 pmax,p0,tol,p
	open(unit=7,file='hw2.txt',status='old',position='append')
11	write(6,"(50('_'))")
	write(7,"(50('_'))")
10	write(6,*)'Input k to choose function'
	write(6,*)'1 f(x)=x*x*x+4.d0*x*x-10.d0'
	write(6,*)'2 f(x)=dcos(x)-x'
	read(5,*)k
	if(k == 1)then
		write(6,*)'  f(x)=x*x*x+4.d0*x*x-10.d0'
		write(7,*)'f(x)=x*x*x+4.d0*x*x-10.d0'
	else if(k == 2)then
		write(6,*)'  f(x)=dcos(x)-x'
		write(7,*)'f(x)=dcos(x)-x'
	else 
		write(6,*)'goto'
		goto 10
	end if

    write(6,*)'Input nmax'
	read(5,*)nmax
	write(6,*)'Input pmax'
	read(5,*)pmax
	write(6,*)'Input p0 '
	read(5,*) p0 
	write(6,*)'Input tol'
	read(5,*)tol

	write(6,"(40('_'))")
	write(6,"(A10,I15)")'nmax =',nmax
	write(6,"(A10,F15.10)")'pmax =',pmax
	write(6,"(A10,F15.10)")'tol =',tol
	write(6,"(A10,F15.10)")'p0 =',p0
	write(6,"(40('_'))")

	write(7,"(40('_'))")
	write(7,"(A10,I15)")'nmax =',nmax
	write(7,"(A10,F15.10)")'pmax =',pmax
	write(7,"(A10,F15.10)")'tol =',tol
	write(7,"(A10,F15.10)")'p0 =',p0
	write(7,"(40('_'))")

	! Body of HW2
	call newton(k, nmax,p0,pmax,tol,p)
	write(6,"(40('_'))")
	write(6,*)'Answer is',p
	write(7,"(40('_'))")
	write(7,*)'Answer is',p

	write(6,*)'1-go on or other quit'
	read(5,*)k
	if(k == 1)then
		goto 11
	else 
	end if
	end program HW2


	function f1(k,x)
	real*8 x
	if(k==1)then
		f1 = x - (x*x*x+4.d0*x*x-10.d0)/(3.d0*x*x+8.d0*x)
	elseif(k==2)then
		f1=x-(dcos(x)-x)/(-dsin(x)-1)
	else
		stop
	end if
	end function f1

	subroutine newton(k, nmax,p0,pmax,tol,p)
	integer nmax
	real*8 p0,pmax,TOL,p

	write(6,"(A3,5x,A10)")'n','p'
	write(7,"(A3,5x,A10)")'n','p'
	do i=1, nmax
	  write(6,"(I3,5x,F15.10)")i-1,p0
	  write(7,"(I3,5x,F15.10)")i-1,p0
	  p=f1(k, p0)
	  if(dabs((p-p0)/p)<tol)then
		write(6,"(I3,5x,F15.10)")i,p
		write(7,"(I3,5x,F15.10)")i,p
	    return
	  elseif(dabs(p)>pmax)then
	    write(6,*)'Divergent'
		write(7,*)'Divergent'
	    stop
	  else
	    p0=p
	  end if
    end do
    write(6,*)'Over Nmax'
	write(7,*)'Over Nmax'
	stop
	end subroutine newton
