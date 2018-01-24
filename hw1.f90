!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HW1
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	program hw1
	  integer nmax,np,k
	  real*8 tol
	  complex*16 a,a1,b,b1,q,h,c(10,2),f
	  open(unit=8,file='hw1_input.txt',status='old',action='read')
	  open(unit=11,file='hw1_output.txt',action='write',position='append')
	  read(8,*)nmax,np,tol,a,b
	  call hw1_binsea(nmax,np,a,b,tol)
	end program hw1

	subroutine hw1_binsea(nmax,np,a,b,tol)
	  integer nmax,np,k,o
	  real*8 tol
	  complex*16 a,a1,b,b1,q,h,c(10,2),f
	  o=6
	  h=(b-a)/np
	  j=0
	  a1=a
	  do i=0,np
	    a1=a+i*h
	    call hassol(a1,a1+h,k)
	    if(k .eq. 1)then
	      j=j+1
	      c(j,1)=a1
	      c(j,2)=a1+h
	    else
	    end if
	  end do
	  do i=1,j
	    call binsea(nmax,c(i,1),c(i,2),tol,q)
	    write(o,"(2X,E15.8,2X,E15.8,2X,E15.8)")real(c(i,1)),real(c(i,2)),real(q)
	  end do
	end subroutine hw1_binsea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine hassol
! k--0 has no answer
!  --1 has answer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine hassol(a,b,k)
          integer k
	  complex*16 a,b,t,f
	  if((real(f(a))*real(f(b))) .gt. 0.d0)then
	    k=0
	    return
	  end if	  
	  t=f((a+b)/2.d0)
	  if(real(f(a))*real(t) .gt. 0.d0)then
	    if(cdabs(a) .gt. cdabs(t))then
	      k=1
	    else
	      k=0
	    end if
	  else
	    if(cdabs(f(b)) .gt. cdabs(t))then
	      k=1
	    else
	      k=0
	    end if
	  end if
	end subroutine hassol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine binsea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine binsea(nmax,a,b,tol,q)
	  integer nmax,i
	  real*8 tol
	  complex*16 a,b,q,t,f
	  do i=1,nmax
	    t=(a+b)/2.d0
	    if((cdabs((b-t)/t))<tol)then
	      q=t
	      return
	    end if

	    if((real(f(a))*real(f(t))) .GT. 0.d0)then
	      a=t
	    else
	      b=t
	    end if
	  end do
	  write(6,*)'iteration times over',nmax
	  q=t
	  return
	end subroutine binsea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function f(x)
	  complex*16 x,f
	  !f=x*x*x+4.d0*x*x-10.d0
	  f=cdcos(x)
	end function f