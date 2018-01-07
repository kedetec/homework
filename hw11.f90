!****************************************************************************
!
!  PROGRAM: hw11
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
	
	program hw11
	  integer n
	  complex*16 a,b
	  open(unit=8,file='hw11_input.txt',status='old',action='read')
	  open(unit=11,file='hw11_output.txt',action='write',position='append')
	  read(8,*)n,a,b
	  call runge_kuta(a,b,n)
	end program hw11

	subroutine runge_kuta(a,b,n)
	  integer n,i,o
	  complex*16 a,b,k1,k2,k3,k4,h,y1,y2,t,fp,f
	  o=6
	  y1=dcmplx(0.5,0.d0)
	  h=(b-a)/n
	  write(o,"(50('*'))")
	  write(o,"(10X,'Runge-Kutta')")
	  write(o,"(40X,'h =',F4.2)")real(h)
	  write(o,"('ti',10X,'Exact',6X,'OrderFour',5X,'Error')")	  
	  write(o,"(50('-'))")
	  do i=0,n
	    t=a+i*h
	    k1=fp(t,y1)
		k2=fp(t+h/2,y1+k1*h/2)
		k3=fp(t+h/2,y1+k2*h/2)
		k4=fp(t+h,y1+k3*h)
		y2=y1+h/6*(k1+2*k2+2*k3+k4)

!		k1=h*fp(t,y1)
!		k2=h*fp(t+h/2,y1+k1/2)
!		k3=h*fp(t+h/2,y1+k2/2)
!		k4=h*fp(t+h,y1+k3)
!		y2=y1+(k1+2*k2+2*k3+k4)/6
		write(o,100)real(t),real(f(t)),real(y1),cdabs(f(t)-y1)
		y1=y2
	  end do  
	  write(o,"(50('-'))")
	  !complex
	  200 FORMAT('(', E11.5,2X,E11.5')')  
	  !real
	  100 FORMAT(F3.1,3X,E13.7,2X,E13.7,2X,E13.7)
	end subroutine runge_kuta

	function fp(t,y)
	  complex*16 t,y,fp
	  fp=y-t*t+1.d0
	end function fp

	function f(t)
	  complex*16 t,f
	  f=-0.5*cdexp(t)+t*t+2.d0*t+1.d0
	end function f
