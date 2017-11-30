!  hw5.f90 
!
!  FUNCTIONS:
!	hw5      - Entry point of console application.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	program hw5
	integer n
	! Variables
	open(unit=11,file='hw5output.txt',action='write',position='append')
	open(unit=7,file='input_findx.txt',status='old',action='read')
	read(7,*)n
	call hw5_findx(n)
	end program hw5

	subroutine hw5_findx(n)
	integer n,i,j
	complex*16 a(n,n),b(n),c(n)
	do i=1,n
	  read(7,*)(a(i,j),j=1,n),b(i)	  
	end do
	write(11,"(60('-'))")
	write(11,"('findx')")
	write(11,"(60('-'))")	
	call output_matrix(n,n,a)
	write(11,"('B')")
	call output_matrix(1,n,b)
	call matrix_findx(n,a,b,c)
	write(11,"('X')")
	call output_matrix(1,n,c)
	close(7)
	! check the result
	open(unit=7,file='input_findx.txt',status='old',action='read')
	read(7,*)n
	do i=1,n
	  read(7,*)(a(i,j),j=1,n),b(i)	  
	end do
	call check_axb(n,a,c,b)
	end subroutine hw5_findx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! matrix_findx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine matrix_findx(n,a,b,x)
	integer n,i,j,k
	complex*16 a(n,n),aa(n),b(n),x(n),amax,c0,sum,t
	c0=dcmplx(0.d0,0.d0)
	! generate aa(n)
	do i=1,n
	  aa(i)=a(i,1)
	  do j=1,n
	    if(cdabs(a(i,j)) .gt. cdabs(aa(i)))then
		  aa(i)=a(i,j)
		end if
	  end do
	  if(aa(i) .eq. c0)then
	    write(6,*)'amax = 0. No answer!'
		stop
	  end if
	end do
	! Gauss elimination and Backward substitution
	do i=1,n-1
	  amax=a(i,i)/aa(i)
	  k=i
	  do j=i+1,n
	    if(cdabs(a(j,i)/aa(j)) .gt. cdabs(amax))then
		  amax=a(j,i)/aa(j)
		  k=j
		end if
	  end do
	  if(k .gt. i)then
  	    call exchange_col(i,k,n,a)
		call exchange_ab(aa(i),aa(k))
		call exchange_ab(b(i),b(k))
	  end if
	  ! elimination
	  do j=i+1,n	 
	    r=a(j,i)/a(i,i)
		do k=1,n
		  a(j,k)=a(j,k)-r*a(i,k)
		end do
		b(j)=b(j)-r*b(i)
	  end do  
	end do
	! Backward substitution
	x(n)=b(n)/a(n,n)
	do i=n-1,1,-1
	  sum=c0
	  do j=i+1,n
	    sum=sum+a(i,j)*x(j)
	  end do
	  x(i)=(b(i)-sum)/a(i,i)
	end do
	end subroutine matrix_findx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exchange a and b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine exchange_ab(a,b)
	complex*16 a,b,t
    t=a
	a=b
	b=t
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exchange matrix's col i and col j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine exchange_col(i,j,n,a)
	integer i,j,n,k
	complex*16 a(n,n),t
    do k=1,n
	  call exchange_ab(a(i,k),a(j,k))
	end do
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output m*n matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine output_matrix(m,n,c)
	integer m,n,i,j
	complex*16 c(m,n)
    do i=1,m
	  write(11,"(1x,E20.6,$)")(c(i,j),j=1,n)
	  write(11,*)
	end do	
	return
	end subroutine output_matrix
!******************************************************************************
! check a*x=b
! return m--0 failure
!         --1 successful
!******************************************************************************
	subroutine check_axb(n,a,x,b)
	integer n,m,i,j
	complex*16 a(n,n),x(n),b(n),sum,c0,t
	c0=dcmplx(0.d0,0.d0)
	
	do i=1,n
	  sum=c0
	  do j=1,n
	    sum=sum+a(i,j)*x(j)
	  end do	  
	  if(sum .ne. b(i))then
		if(cdabs((sum-b(i))/sum) .gt. 0.0000001)then
	      write(6,*)'A * X = B is wrong!'
		  return
		end if
	  end if
	end do
	write(6,*)'A * X = B is correct!'
	return
	end subroutine check_axb
