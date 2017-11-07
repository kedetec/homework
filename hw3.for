!  HW3.for 
!
!  FUNCTIONS:
!	HW3      - Entry point of console application.
	program HW3
	implicit none

	! Variables
	integer nmax,k
	real*8 tol
	complex*16 p
	open(unit=7,file='hw3.txt',position='append')
11	write(6,"(80('_'))")
	write(7,"(80('-'))")
	write(7,*)'f(x) = 16.d0*x*x*x*x-40.d0*x*x*x+5.d0*x*x+20.d0*x+6.d0'
      write(6,*)'Input nmax'
	read(5,*)nmax
	write(6,*)'Input tol'
	read(5,*)tol

	write(6,"(40('_'))")
	write(6,"(2x,'nmax = ',I4,'  tol = ',F10.6)")
     c nmax,tol
	write(6,"(40('_'))")

	write(7,"(80('-'))")
	write(7,"(2x,'nmax = ',I4,'  tol = ',F10.6)")
     c nmax,tol
	write(7,"(80('-'))")

	! Body of HW3
	call muller(nmax,tol,p)
	write(6,"(40('_'))")
	write(6,*)'Answer is',p
	write(7,"(80('-'))")
	write(7,*)'Answer is',p

	write(6,*)'1-go on or other quit'
	read(5,*)k
	if(k == 1)then
		goto 11
	else 
	end if
	end program HW3


	function f(x)
	complex*16 x,f
		f = 16.d0*x*x*x*x-40.d0*x*x*x+5.d0*x*x+20.d0*x+6.d0	    
	end function f

	subroutine muller(nmax,tol,p)
	integer nmax
	real*8 tol
	complex*16 p,x0,x1,x2,h1,h2,h,d,t1,t2,t,f,b
	write(6,*)'Input x0,x1,x2'
	read(5,*)x0,x1,x2

	write(7,"('  x0 =',F5.2,'  x1 =',F5.2,'  x2 =',F5.2)")
     c dreal(x0),dreal(x1),dreal(x2)
	write(7,"(2(' '),'i',15(' '),'x',20x,'p(x)')")
	write(7,"(80('-'))")
	
	h1 = x1 - x0
	h2 = x2 - x1
	t1 = (f(x1)-f(x0))/h1
	t2 = (f(x2)-f(x1))/h2
	d = (t2-t1)/(h1+h2)    

	do i=0, nmax
	    b=t2+h2*d
	    t = cdsqrt(b*b-4*f(x2)*d)
		if(cdabs(b-t)<cdabs(b+t))then
		  h=-2.d0*f(x2)/(b+t)
		else
		  h=-2.d0*f(x2)/(b-t)
		end if
		p=x2+h
	write(7,"(2x,I2,F12.6,' +',F12.6,'i',3x,F12.6,' +',F12.6,'i')")
     c i+2,dreal(p),aimag(p),dreal(f(p)),aimag(f(p))

		if(cdabs(h)<tol)then	    	
			return
		end if
		
		x0=x1
		x1=x2
		x2=p
		h1 = x1 - x0
		h2 = x2 - x1
		t1 = (f(x1)-f(x0))/h1
		t2 = (f(x2)-f(x1))/h2
		d = (t2-t1)/(h1+h2)
	end do
      write(6,*)'Over Nmax'
	stop
	end subroutine muller