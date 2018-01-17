!  HW3.for 
!
!  FUNCTIONS:
!	HW3      - Entry point of console application.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	program HW3

	  integer nmax,o
	  real*8 tol,pmax
	  complex*16 x0,x1,x2,p
	  open(unit=7,file='hw3_input.txt',status='old',action='read')
	  open(unit=11,file='hw3_output.txt',action='write',
     c      position='append')
	  o=6
	  read(7,*)nmax,pmax,tol,x0,x1,x2
	  call hw3_muller(nmax,tol,x0,x1,x2,p)
	
	  write(o,"(80('-'))")
	  write(o,*)'f(x) = 16.d0*x*x*x*x-40.d0*x*x*x+5.d0*x*x+20.d0*x+6.d0'
	  write(o,"(2x,'nmax = ',I4,'  tol = ',F10.6)")
     c nmax,tol
	  write(o,"(80('-'))")
	  write(o,*)'Answer is',p
	end program HW3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  function f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine hw3_muller(nmax,tol,x0,x1,x2,p)
	  integer nmax,o
	  real*8 tol
	  complex*16 p,x0,x1,x2,h1,h2,h,d,t1,t2,t,f,b
	  o=6
	  write(o,"('  x0 =',F5.2,'  x1 =',F5.2,'  x2 =',F5.2)")
     c       dreal(x0),dreal(x1),dreal(x2)
	  write(o,"(2(' '),'i',15(' '),'x',20x,'p(x)')")
	  write(o,"(80('-'))")
	
	  h1 = x1 - x0
	  h2 = x2 - x1
	  t1 = (f(x1)-f(x0))/h1
	  t2 = (f(x2)-f(x1))/h2
	  d = (t2-t1)/(h1+h2)    

	  do i=0, nmax
	    b=t2+h2*d
	    t = cdsqrt(b*b-4*f(x2)*d)
		if(cdabs(b-t) .le. cdabs(b+t))then
		  h=-2.d0*f(x2)/(b+t)
		else
		  h=-2.d0*f(x2)/(b-t)
		end if
		p=x2+h
c	write(o,*)p,f(p)
	    write(o,"(2x,I2,F12.6,' +',F12.6,'i',3x,F12.6,' +',F12.6,'i')")
     c       i+2,dreal(p),dimag(p),dreal(f(p)),dimag(f(p))

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
	end subroutine hw3_muller
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  function f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function f(x)
	  complex*16 x,f
	  f = 16.d0*x*x*x*x-40.d0*x*x*x+5.d0*x*x+20.d0*x+6.d0	    
	end function f