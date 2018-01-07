!****************************************************************************
!
!  PROGRAM: hw12
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program hw12
	  integer n
	  complex*16 a,b
	  open(unit=8,file='hw12_input.txt',status='old',action='read')
	  open(unit=11,file='hw12_output.txt',action='write',position='append')
	  read(8,*)n,a,b
	  call runge_kuta(2,a,b,n)
	end program hw12
!****************************************************************************
! subroutine runge_kuta
!****************************************************************************
	subroutine runge_kuta(nf,a,b,n)
	  integer n,nf,i,o
	  complex*16 a,b,k(4,nf),h,y1(nf),y2(nf),t,fp,f
	  o=11
	  y1(1)=dcmplx(0.d0,0.d0)
	  y1(2)=dcmplx(0.d0,0.d0)
	  h=(b-a)/n
	  write(o,"(50('*'))")
	  write(o,"(2X,'f(1)=-3.375*cdexp(-2*t)+1.875*cdexp(-0.4*t)+1.5')")
  	  write(o,"(2X,'f(2)=-2.25*cdexp(-2*t)+2.25*cdexp(-0.4*t)')")
	  write(o,"(40X,'h =',F4.2)")real(h)
	  write(o,"('ti',5X,'Exact F(1)',5X,'OrderFour',5X,'Exact F(2)',5X,'OrderFour')")	  
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
		write(o,100)real(t),real(f(1,t)),real(y1(1)),real(f(2,t)),real(y1(2))
	
		y1(1)=y2(1)
		y1(2)=y2(2)

		100 FORMAT(F3.1,3X,E13.7,2X,E13.7,2X,E13.7,2X,E13.7)
	  end do  
	end subroutine runge_kuta
!****************************************************************************
! function fp function prime
!****************************************************************************
	function fp(nf,t,y1,y2)
	  integer nf
	  complex*16 t,y1,y2,fp
	  if(nf .eq. 1)then
	    fp=-4*y1+3*y2+6.d0
	  else if(nf .eq. 2)then
	    fp=-2.4*y1+1.6*y2+3.6
	  else
	  end if	  
	end function fp
!****************************************************************************
! function f
!****************************************************************************
	function f(nf,t)
	  integer nf
	  complex*16 t,f
	  if(nf .eq. 1)then
	    f=-3.375*cdexp(-2*t)+1.875*cdexp(-0.4*t)+1.5
	  else if(nf .eq. 2)then
	    f=-2.25*cdexp(-2*t)+2.25*cdexp(-0.4*t)
	  else
	  end if
	end function f