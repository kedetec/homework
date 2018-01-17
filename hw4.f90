!  hw4.f90 
!
!  FUNCTIONS:
!	hw4      - Entry point of console application.
	program hw4
	  implicit none
	  integer n,k
100	  write(6,*)'input number to choose inputfile: 1 2 3'
	  read(5,*)k
	  if(k .eq. 1)then
	    open(unit=7,file='hw4input1.txt',status='old',action='read')
	  else if(k .eq. 2)then
	    open(unit=7,file='hw4input2.txt',status='old',action='read')
	  else if(k .eq. 3)then
	    open(unit=7,file='hw4input3.txt',status='old',action='read')
	  else
	    goto 100
	  end if
	  read(7,*)n
	  call hw4ge(n)
	  close(7)
	  write(6,*)'quit press 1, else press others'
	  read(5,*)k
	  if(k .ne. 1)then	 
	    goto 100
	  end if
	end program hw4

	subroutine hw4ge(n)
	  implicit none
	  integer n,i,j,k,o
	  complex*16 a(n,n),b(n),x(n),r,sum,c0
	  c0=dcmplx(0.d0,0.d0)
	  open(unit=11,file='hw4output2.txt',action='write',position='append')
	  o=11
	  do i=1,n
	    read(7,*)(a(i,j),j=1,n),b(i)	
	  end do

	  write(o,*)"coefficient"
	  write(o,"(50('-'))")
	  call outarr(11,n,a)	
	  write(o,"(50('-'))")
	  write(o,*)'b'
	  write(o,100)(b(j),j=1,n)
	  write(o,*)
	  ! Gauss Elimination
	  do i=1,n-1
	    ! a(i,i)==0,find a(j,i) .ne. 0
	    ! change row i and row j
	    if(a(i,i) .eq. c0)then
	      do j=i+1,n
		    if(a(j,i) .ne. c0)then
		      do k=i,n
			    r=a(i,k)
			    a(i,k)=a(j,k)
			    a(j,k)=r			
			  end do
			  r=b(i)
			  b(i)=b(j)
		  	  b(j)=r
			  exit
		    end if
		  end do
		  if(j .gt. n)then
		    write(o,*)'No answer'
		    return
		  end if
	    end if

	  ! elimination
	    do j=i+1,n
	      r=a(j,i)/a(i,i)
		  do k=i+1,n
		    a(j,k)=a(j,k)-r*a(i,k)
		  end do
		    b(j)=b(j)-r*b(i)
	    end do  
	  end do
	  ! Backward substitution
	  if(a(n,n) .eq. 0)then
	    write(o,*)'No answer'
	    return
	  end if
	  x(n)=b(n)/a(n,n)
	  do i=n-1,1,-1
	    sum=c0
	    do j=i+1,n
	      sum=sum+a(i,j)*x(j)
	    end do
	    x(i)=(b(i)-sum)/a(i,i)
	  end do
	! output solutions
      write(o,"(50('-'))")
	  write(o,*)'x'
      write(o,100)(x(i),i=1,n)
	  write(o,*)
	  close(11)
  	  return
	  100 FORMAT('(',E16.9,2X,E16.9,')')
	end subroutine hw4ge
	
	subroutine outarr(m,n,a)
	  implicit none
	  integer m,n,i,j
	  complex*16 a(n,n)
	  do i=1,n
	    write(m,100)(a(i,j),j=1,n)
	    write(m,*)
	  end do
	100 FORMAT('(',E16.9,2X,E16.9,')',$)
	end subroutine outarr