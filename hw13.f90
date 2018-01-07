!****************************************************************************
!
!  PROGRAM: hw13
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program hw13
	  integer n
	  complex*16 a,b
	  open(unit=8,file='hw13_input.txt',status='old',action='read')
	  open(unit=11,file='hw13_output.txt',action='write',position='append')
	  read(8,*)n,a,b
	  call runge_kuta(2,a,b,n)
	end program hw13
!****************************************************************************
! subroutine runge_kuta
!****************************************************************************
	subroutine runge_kuta(nf,a,b,n)
	  integer n,nf,i,o
	  complex*16 a,b,k(4,nf),h,y1(nf),y2(nf),t,fp
	  o=11
	  y1(1)=dcmplx(-0.4,0.d0)
	  y1(2)=dcmplx(-0.6,0.d0)
	  h=(b-a)/n
	  write(o,"(50('*'))")
	  write(o,"(10X,'Second-order initial value')")
	  write(o,"(40X,'h =',F4.2)")real(h)
	  write(o,"('ti',10X,'Y',12X,'Y prime')")	  
	  write(o,"(50('-'))")
	  do i=0,n
	    t=a+i*h
		do j=1,nf
	      k(1,j)=h*fp(j,t,y1(1),y1(2))
		end do
		do j=1,nf
		  k(2,j)=h*fp(j,t+h/2,y1(1)+k(1,1)/2,y1(2)+k(1,2)/2)
		end do
		do j=1,nf
		  k(3,j)=h*fp(j,t+h/2,y1(1)+k(2,1)/2,y1(2)+k(2,2)/2)
		end do
		do j=1,nf
		  k(4,j)=h*fp(j,t+h,y1(1)+k(3,1),y1(2)+k(3,2))
		end do
		do j=1,nf
		  y2(j)=y1(j)+(k(1,j)+2.d0*k(2,j)+2.d0*k(3,j)+k(4,j))/6
		end do
		write(o,100)real(t),real(y1(2)),real(y1(1))
	
		y1(1)=y2(1)
		y1(2)=y2(2)		
	  end do  
	  100 FORMAT(F3.1,3X,E14.7,2X,E14.7)
	end subroutine runge_kuta
!****************************************************************************
! function fp function prime
!****************************************************************************
	function fp(nf,t,y1,y2)
	  integer nf
	  complex*16 t,y1,y2,fp
	  if(nf .eq. 1)then
	    fp=y2
	  else if(nf .eq. 2)then
	    fp=cdexp(2*t)*cdsin(t)-2.d0*y1+2.d0*y2
	  else
	  end if	  
	end function fp




