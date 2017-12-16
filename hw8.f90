!****************************************************************************
!
!  PROGRAM: hw8
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
	program hw8
	integer n
	open(unit=8,file='hw8_input.txt',status='old',action='read')
	open(unit=11,file='hw8_output.txt',action='write',position='append')
	read(8,*)n
	call hw8_nonliner(n)
	close(8)
	close(11)
	end program hw8
!*******************************************************************
!  Subroutine hw8_nonliner
!*******************************************************************
	subroutine hw8_nonliner(n)
	integer n,i,j,nmax
	real*8 tol,pmax,xo(n),xx(n),a(n,n),b(n),xn(n),me,tt
	read(8,*)nmax,tol,pmax
	read(8,*)(xo(i),i=1,n)
	!output
	write(11,"(80('-'))")
	write(11,*)' 3*x1-dcos(x2*x3)-0.5 = 0'
	write(11,*)' x1*x1-81*(x2+0.1)*(x2+0.1)+dsin(x))+1.06 = 0'
	write(11,*)' dexp(-x1*x2)+20.d0*x3+(10*pi-3.d0)/3.d0 = 0'
	write(11,"(80('-'))")
	write(11,"(2X,'tol =',F9.6,3X,'pmax =',F9.3)")tol,pmax
	write(11,"(80('-'))")
	write(11,"(2X,'k',10X,'x1',10X,'X2',10X,'x3')")
	write(11,"(I2,3e15.6)")0,(xo(j),j=1,n)
	do i=1,nmax
	  me=0.d0	  
	  call cm(3,xo,a)
	  call cf(3,xo,b)
	  call gebs(3,a,b,xx)
	  do j=1,3
	    xn(j)=xo(j)+xx(j)
		tt=dabs(xx(j)/xn(j))
		if(me .le. tt)then
		  me=tt
		end if
		if(dabs(xx(j)) .gt. pmax)then
		  write(6,*)'Divergence!'
		  write(11,*)'Divergence!'
		  stop
		end if
	  end do	  
	  do j=1,3
	    xo(j)=xn(j)
	  end do
	  write(11,"(I2,3e15.6)")i,(xo(j),j=1,n)
	  if(me .le. tol)then
	    write(6,*)'Find solutions'
		write(11,*)'Find solutions'
		exit
	  end if
	end do
	if(i .gt. nmax)then
	  write(11,*)'Over iteration times!'
	end if
	end subroutine hw8_nonliner
!*******************************************************************
!  Subroutine calculate matrix
!*******************************************************************
	subroutine cm(n,xo,a)
	  integer n,i,j
	  real*8 a(n,n),xo(n)
	  do i=1,n
	    do j=1,n
	      a(i,j)=fp(i,j,xo)
		end do
	  end do
	  return
	end subroutine cm
!*******************************************************************
!  Subroutine calculate f
!*******************************************************************
	subroutine cf(n,xo,b)
	  integer n,i,j
	  real*8 b(n),xo(n)
	  do i=1,n
	    b(i)=-f(i,xo)
	  end do
	  return
	end subroutine cf
!*******************************************************************
!  Subroutine Gauss Elimination and backward substitution
!*******************************************************************
	subroutine gebs(n,a,b,xx)
	  integer n,i,j,k
	  real*8 a(n,n),b(n),xx(n),r,sum
  	  do i=1,n
	    do j=i+1,n
	      r=a(j,i)/a(i,i)
		  do k=1,n
		    a(j,k)=a(j,k)-r*a(i,k)
		  end do
		  b(j)=b(j)-r*b(i)
	    end do
	  end do
	  xx(n)=b(n)/a(n,n)
	  do i=n-1,1,-1
	    sum=0.d0
		do j=i+1,n
		  sum=sum+a(i,j)*xx(j)
		end do
		xx(i)=(b(i)-sum)/a(i,i)
	  end do	
	end subroutine gebs
!*******************************************************************
!  functions
!*******************************************************************
	function f(i,xo)
	  integer i
	  real*8 f,xo(3),pi,t
	  pi=3.141592654
	  if(i .eq. 1)then
	    f=3*xo(1)-dcos(xo(2)*xo(3))-0.5
	  else if(i .eq. 2)then
	    f=xo(1)*xo(1)-81*(xo(2)+0.1)*(xo(2)+0.1)+dsin(xo(3))+1.06
	  else if(i .eq. 3)then 
	    f=dexp(-xo(1)*xo(2))+20.d0*xo(3)+(10*pi-3.d0)/3.d0
	  else
	  end if
	end function f

!*******************************************************************
!  functions prime
!*******************************************************************
	function fp(i,j,xo)
	  integer i,j,k
	  real*8 fp1,xo(3)	  
	  if(i .eq. 1)then
		if(j .eq. 1)then
	      fp=3.d0
		else if(j .eq. 2)then
		  fp=xo(3)*dsin(xo(2)*xo(3))
		else if(j .eq. 3)then
		  fp=xo(2)*dsin(xo(2)*xo(3))
		else
		  ! parameter error
		end if
	  else if(i .eq. 2)then
	    if(j .eq. 1)then
	      fp=2.d0*xo(1)
		else if(j .eq. 2)then
		  fp=-81.d0*2.d0*(xo(2)+0.1)
		else if(j .eq. 3)then
		  fp=dcos(xo(3))
		else
		  ! parameter error
		end if
	  else if(i .eq. 3)then
	    if(j .eq. 1)then
	      fp=-xo(2)*dexp(-xo(1)*xo(2))
		else if(j .eq. 2)then
		  fp=-xo(1)*dexp(-xo(1)*xo(2))
		else if(j .eq. 3)then
		  fp=20.d0
		else
		  ! parameter error
		end if
	  else
	  end if
	end function fp
