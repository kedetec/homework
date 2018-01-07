!****************************************************************************
!
!  PROGRAM: hw10
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program hw10
	  integer nf,n,m,o
	  complex*16 a,b,s
	  open(unit=8,file='hw10_input.txt',status='old',action='read')
  	  open(unit=11,file='hw10_output.data',action='write',position='append')	
	  read(8,*)m,nf,n,a,b
	  o=11
	  if(nf .eq. 1)then
	    write(o,*)"(f=cdcos(x)*cdexp(x))"
	  else if(nf .eq. 2)then
	    write(o,*)'f=(100.d0/(x*x))*cdsin(10/x)'
	  else if(nf .eq. 3)then
	    write(o,*)'f=cdexp(-x*x)'
	  else if(nf .eq. 4)then
	    write(o,*)'f=cdcos(x)'
	  else
	  end if
	  write(o,"(2X,'a = ','(',2F8.6,'),',2X,'b = ','(',2F8.6,'),',2X,'n = ',I3)")a,b,n
	  if(m .eq. 1)then
	    write(o,*)'Somposite Simpson''s Rule'
		call compsimpson(nf,n,a,b,s)
	  else if(m .eq. 2)then
		write(o,*)'Gauss Quadrature'
		call gauss(nf,n,a,b,s)
	  endif
	  write(o,"(2X,'(',2E20.12,')')")s
	  write(o,"(70('*'))")
	end program hw10
!******************************************************************************
! subroutine compsimpson
!******************************************************************************
	subroutine compsimpson(nf,n,a,b,s)
	  integer nf,n,i,oe
	  complex*16 a,b,h,x(0:n),s,f
      h=(b-a)/n
	  do i=0,n
	    x(i)=a+i*h
	  end do
	  s=(h/3)*(f(nf,x(0))+f(nf,x(n)))
	  oe=1
	  do i=1,n-1
	    if(oe .eq. 1)then
		  s=s+(h/3)*4*f(nf,x(i))
		  oe=0
		else
		  s=s+(h/3)*2*f(nf,x(i))
		  oe=1
		end if
	  end do
	end subroutine compsimpson
!******************************************************************************
! subroutine gauss
!******************************************************************************
	subroutine gauss(nf,n,a,b,s)
	  integer nf,n,i
	  complex*16 x(n),c(n),a,b,t(n),fx,f,s
	  if(n .eq. 2)then
	    open(unit=9,file='gauss_coefficient_2.txt',status='old',action='read')
	  else if(n .eq. 3)then
	    open(unit=9,file='gauss_coefficient_3.txt',status='old',action='read')
	  else if(n .eq. 4)then
	    open(unit=9,file='gauss_coefficient_4.txt',status='old',action='read')
	  else if(n .eq. 5)then
	    open(unit=9,file='gauss_coefficient_5.txt',status='old',action='read')
	  else
	  end if
	  read(9,*)(x(i),i=1,n)
	  read(9,*)(c(i),i=1,n)
		
	  do i=1,n
	    t(i)=fx(a,b,x(i))
	  end do
	  s=dcmplx(0,0)
	  do i=1,n
	    s=s+c(i)*f(nf,t(i))
	  end do
	  s=s*((b-a)/2)
	end subroutine gauss
!******************************************************************************
! function f
!******************************************************************************
	function f(nf,x)
	  integer nf
	  complex*16 f,x
	  if(nf .eq. 1)then
	    f=cdcos(x)*cdexp(x)
	  else if(nf .eq. 2)then
	    f=(100.d0/(x*x))*cdsin(10/x)
	  else if(nf .eq. 3)then
	    f=cdexp(-x*x)
	  else if(nf .eq. 4)then
	    f=cdcos(x)
	  else
	  end if	  
	end function f
!******************************************************************************
! function fx 
! fx passed to function f() as parameter x
!******************************************************************************
	function fx(a,b,t)	  
	  complex*16 a,b,t,fx
	  fx=((b-a)*t+(b+a))/2
	end function fx

