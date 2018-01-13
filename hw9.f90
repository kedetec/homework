!****************************************************************************
!
!  PROGRAM: hw9
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program hw9
	  integer n
	  complex*16 a,b
	  open(unit=8,file='hw9_input.txt',status='old',action='read')
  	  open(unit=12,file='hw9_output.data',action='write',position='append')
	  open(unit=11,file='hw9_output.txt',action='write',position='append')
	  read(8,*)n,a,b
	  call hw9_newtondd(n,a,b)	

	end program hw9
!****************************************************************************
! subroutine hw9_newtondd(n,a,b)
!****************************************************************************
	subroutine hw9_newtondd(n,a,b)
	  integer n,i,j,o
	  real*8 maxerr,err,errx
	  complex*16 a,b,x(0:n),f(0:n,0:n),h,y,p,t1,t2,fun,pn
	  call ndd(n,a,b,x,f)
	  ! calculate error maxmum
	  h=(b-a-dble(0.2))/n
	  maxerr=0.d0
	  y=a
	  do i=1,200
!	    y=a+0.1+i*h
        y=y+dble(0.1)
		p=pn(n,x,f,y)
		err=cdabs((p-fun(y))/fun(y))
		write(6,*)i,err
		if(maxerr .le. err)then
		  maxerr=err
		  errx=real(y)
		end if
	  end do
	  write(11,"('a=('F6.3,',',F6.3,')',2X,'b=('F6.3,',',F6.3,')',4X,'n=',I3)")a,b,n
	  write(11,"('y=',E20.13,2X,'maxerror=',E20.13)")errx,maxerr
	  !output data with divided in 999 parts
	  h=(b-a-dble(0.06))/999.d0
	  do i=0,999
	    t1=fun(a+i*h)
		t2=pn(n,x,f,(a+i*h))
		err=cdabs((t2-t1)/t1)
		write(12,"(I4,2X,5E20.12)")i,t1,t2,err
		if(maxerr .le. err)then
		  maxerr=err
		  write(6,"(F16.13)")maxerr
		end if
	  end do
	end subroutine hw9_newtondd

	subroutine ndd(n,a,b,x,f)
	  integer n,i,j,o
	  real*8 maxerr,err,errx
	  complex*16 a,b,x(0:n),f(0:n,0:n),h,y,p,t1,t2,fun,pn
	  h=(b-a)/n
	  do i=0,n
	    x(i)=a+i*h
	    f(i,0)=fun(x(i))
	  end do
	  do i=1,n
	    do j=1,i
		  f(i,j)=(f(i,j-1)-f(i-1,j-1))/(x(i)-x(i-j))
	    end do
	  end do
	end subroutine ndd
!****************************************************************************
! function polynomial
!****************************************************************************
	function pn(n,x,f,y)
	  integer n,i,j
	  complex*16 x(0:n),f(0:n,0:n),y,tt,pn
	  pn=f(0,0)
	  do i=1,n
	    tt=f(i,i)
	    do j=0,i-1
		  tt=tt*(y-x(j))
		end do
		pn=pn+tt
	  end do
	end function pn
!****************************************************************************
! function polynomial
!****************************************************************************
	function fun(x)
	  complex*16 x,fun
	  fun=cdcos(x)
!	  fun=1.d0/x
	end function fun
