!  hw6.f90 
!
!  FUNCTIONS:
!	hw6      - Entry point of console application.
!

!****************************************************************************
	program hw6

	implicit none

	! Variables
	integer n,ic

	open(unit=11,file='hw6output.txt',action='write',position='append')
	write(6,*)'Matrix calculation:'
	write(6,*)'  1-Matrix add'
	write(6,*)'  2-Matrix dotproduct'
	write(6,*)'  3-Matrix determinant'
	write(6,*)'  4-Matrix inverse'
	write(6,*)'  5-Matrix findx'
	write(6,*)'Input value between 1 and 5!'
10	read(5,*)ic
	if(ic .eq. 1)then
	  open(unit=7,file='input_add.txt',status='old',action='read')
	  read(7,*)n
	  call hw6_add(n)
	else if(ic .eq. 2)then
	  open(unit=7,file='input_dotproduct.txt',status='old',action='read')
	  read(7,*)n
	  call hw6_dotproduct(n)
	else if(ic .eq. 3)then
	  open(unit=7,file='input_determinant.txt',status='old',action='read')
	  read(7,*)n
	  call hw6_determinant(n)
	else if(ic .eq. 4)then
	  open(unit=7,file='input_inverse.txt',status='old',action='read')
	  read(7,*)n
	  call hw6_inverse(n)
	else if(ic .eq. 5)then
	  open(unit=7,file='input_findx.txt',status='old',action='read')
	  read(7,*)n
	  call hw6_findx(n)
!	else if(ic .eq. 6)then
	  
	else
	  write(6,*)'Input value between 1 and 5!'
	  goto 10
	end if
	close(unit=11)
	end program hw6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Homework subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine hw6_add(n)
	integer n,i,j
	complex*16 a(n,n),b(n,n),c(n,n)
	do i=1,n
	  read(7,*)(a(i,j),j=1,n)	
	end do
	do i=1,n
	  read(7,*)(b(i,j),j=1,n)		  
	end do
	write(11,"(60('-'))")
	write(11,"('A','+','B')")
	write(11,"(60('-'))")	
	call output_matrix(n,n,a)
	write(11,"(4x,'+')")
	call output_matrix(n,n,b)
	write(11,"(4x,'=')")
	call matrix_add(n,a,b,c)
	call output_matrix(n,n,c)
	end subroutine hw6_add

	subroutine hw6_dotproduct(n)
	integer n,i,j
	complex*16 a(n,n),b(n,n),c(n,n)
	do i=1,n
	  read(7,*)(a(i,j),j=1,n)		  
	end do
	do i=1,n
	  read(7,*)(b(i,j),j=1,n)		  
	end do
	write(11,"(60('-'))")
	write(11,"('A','*','B')")
	write(11,"(60('-'))")	
	call output_matrix(n,n,a)
	write(11,"(4x,'*')")
	call output_matrix(n,n,b)
	write(11,"(4x,'=')")
	call matrix_dotproduct(n,a,b,c)
	call output_matrix(n,n,c)
	end subroutine hw6_dotproduct
	
	subroutine hw6_determinant(n)
	integer n,i,j
	complex*16 a(n,n),c
	do i=1,n
	  read(7,*)(a(i,j),j=1,n)		  
	end do
	write(11,"(60('-'))")
	write(11,"('determinant')")
	write(11,"(60('-'))")	
	call output_matrix(n,n,a)
	write(11,"(2x,'determinant = ',$)")
	call matrix_determinant(n,a,c)
	call output_matrix(1,1,c)
	end subroutine hw6_determinant

	subroutine hw6_inverse(n)
	integer n,i,j
	complex*16 a(n,n),c(n,n),d(n,n)
	do i=1,n
	  read(7,*)(a(i,j),j=1,n)		  
	end do
	write(11,"(60('-'))")
	write(11,"('inverse')")
	write(11,"(60('-'))")	
	call output_matrix(n,n,a)
	call matrix_inverse(n,a,c)
	write(11,"(2x,'inverse = ')")
	call output_matrix(n,n,c)
	end subroutine hw6_inverse

	subroutine hw6_findx(n)
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
	end subroutine hw6_findx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine matrix_add(n,a,b,c)
	integer i,j
	complex*16 a(n,n),b(n,n),c(n,n)
	do i=1,n
	  do j=1,n
	    c(i,j) = a(i,j)+b(i,j)
	  end do
	end do
	return
	end subroutine matrix_add

	subroutine matrix_dotproduct(n,a,b,c)
	integer i,j
	complex*16 a(n,n),b(n,n),c(n,n),sum,c0
	c0=cmplx(0.d0,0.d0)
	do i=1,n
	  do j=1,n
	    sum = c0
		do k=1,n
	      sum = sum+(a(i,k)*b(k,j))
		end do
		c(i,j) = sum
	  end do
	end do
	return
	end subroutine matrix_dotproduct

	subroutine matrix_determinant(n,a,det)
	integer n,i,j,m
	complex*16 a(n,n),det,c1
	c1=cmplx(1.d0,0.d0)
    call matrix_lu1(n,a,m)
	det = c1*m
	do i=1,n
	  det=det*a(i,i)
	end do
	return
	end subroutine matrix_determinant

	subroutine matrix_inverse(n,a,c)
	integer i,j
	complex*16 a(n,n),c(n,n),b(n),x(n),c0,c1
	c0=cmplx(0.d0,0.d0)
	c1=cmplx(1.d0,0.d0)
	do i=1,n
	  b(i)=0
	end do
	call matrix_lu(n,a)
	do i=1,n
	  b(i)=c1
	  call matrix_fsbs(n,a,b,x)
	  do j=1,n
	    c(i,j)=x(j)	    
	  end do
	  b(i)=c0
	end do
	return
	end subroutine matrix_inverse

	subroutine matrix_findx(n,a,b,x)
	integer n,i,j
	complex*16 a(n,n),b(n),x(n)
	call matrix_lu(n,a)
	call matrix_fsbs(n,a,b,x)
	return
	end subroutine matrix_findx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output m*n matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine output_matrix(m,n,c)
	integer m,n,i,j
	complex*16 c(m,n)
    do i=1,m
	  write(11,"(1x,F10.6,$)")(c(i,j),j=1,n)
	  write(11,*)
	end do
	
	return
	end subroutine output_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n*n matrix LU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine matrix_lu(n,a)
	integer i,j,k,n
	complex*16 a(n,n),r,c0
	c0=cmplx(0.d0,0.d0)
	! Gauss Elimination
	do i=1,n-1
	  ! a(i,i)==0,find a(j,i) .ne. 0
	  ! change row i and row j
	  if(a(i,i) .eq. c0)then
	    do j=i+1,n
		  if(a(j,i) .ne. c0)then
		    do k=1,n
			  r=a(i,k)
			  a(i,k)=a(j,k)
			  a(j,k)=r			
			end do
			exit
		  end if
		end do
		if(j .gt. n)then
		  write(6,*)'No answer'
		  return
		end if
	  end if
	  ! elimination
	  do j=i+1,n	 
	    r=a(j,i)/a(i,i)
		do k=i,n
		  a(j,k)=a(j,k)-r*a(i,k)
		end do
		a(j,i)=r
	  end do  
	end do
	end subroutine matrix_lu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n*n matrix LU1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine matrix_lu1(n,a,m)
	integer i,j,k,n,m
	complex*16 a(n,n),r,c0
	c0=cmplx(0.d0,0.d0)
	m=1
	! Gauss Elimination
	do i=1,n-1
	  ! a(i,i)==0,find a(j,i) .ne. 0
	  ! change row i and row j
	  if(a(i,i) .eq. c0)then
	    do j=i+1,n
		  if(a(j,i) .ne. c0)then
		    do k=1,n
			  r=a(i,k)
			  a(i,k)=a(j,k)
			  a(j,k)=r			
			end do
			m=m*(-1)
			exit
		  end if
		end do
		if(j .gt. n)then
		  write(6,*)'No answer'
		  return
		end if
	  end if
	  ! elimination
	  do j=i+1,n	 
	    r=a(j,i)/a(i,i)
		do k=i,n
		  a(j,k)=a(j,k)-r*a(i,k)
		end do
!		a(j,i)=r
	  end do  
	end do
	end subroutine matrix_lu1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n*n matrix forward and backward substitution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine matrix_fsbs(n,a,b,x)
	integer n,i,j
	complex*16 a(n,n),b(n),x(n),y(n),sum,c0
	c0=cmplx(0.d0,0.d0)
	! Forward substitution
	y(1)=b(1)
	do i=2,n
	  sum=c0
	  do j=1,i-1
	    sum=sum+a(i,j)*y(j)
	  end do
	  y(i)=b(i)-sum
	end do
	! Backward substitution
	if(a(n,n) .eq. c0)then
	  write(6,*)'No answer'
	  return
	end if
	x(n)=y(n)/a(n,n)
	do i=n-1,1,-1
	  sum=c0
	  do j=i+1,n
	    sum=sum+a(i,j)*x(j)
	  end do
	  x(i)=(y(i)-sum)/a(i,i)
	end do
	end subroutine matrix_fsbs