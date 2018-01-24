!  hw7.f90 
!
!  FUNCTIONS:
!	hw7      - Gauss-Seidel Method
!****************************************************************************

	program hw7
	  integer n
	  open(unit=7,file='hw7_input.txt',status='old',action='read')
	  open(unit=11,file='hw7_output.txt',action='write',position='append')
	  read(7,*)n
	  call hw7_gauss_seidel(n)
	end program hw7
!***************************************************************************
! subroutine hw5_gauss_seidel
!***************************************************************************
	subroutine hw7_gauss_seidel(n)
	  integer n,i,j,k,nmax
	  complex*16 a(n,n),b(n),xo(n),xn(n),c0,sumo,sumn
	  real*8 tol,pmax,xt,t0,w
	  c0=dcmplx(0.d0,0.d0)
	  !init
	  do i=1,n
	    do j=1,n
	      a(i,j)=c0
	    end do
	    b(i)=c0
	    xo(i)=c0
	    xn(i)=c0
	  end do
	  !read data from input_file
	  do i=1,n
	    read(7,*)(a(i,j),j=1,n),b(i)	  
	  end do	
	  read(7,*)
	  read(7,*)nmax,tol,pmax,w
          read(7,*) 
	  read(7,*)(xo(i),i=1,n)
	  !output parameter
	  write(11,"(120('-'))")
	  write(11,"('w = ', e15.6,'  tol =',e15.6)")w,tol
	  write(11,"(120('*'))")
	  write(11,"(2X,'k',18X,'x1',38X,'x2',35X,'x3')")
	  write(11,"(120('-'))")
	  do k=1,nmax
	    do i=1,n
	      sumo=c0
	      sumn=c0
	      do j=1,i-1
	        sumn=sumn+a(i,j)*xn(j)
	      end do
	      do j=i+1,n
	        sumo=sumo+a(i,j)*xo(j)
	      end do
	      xn(i)=(1.d0-w)*xo(i)+w*((b(i)-sumn-sumo)/a(i,i))	  
	    end do
  	    !control
	    sumn=c0 !maxmun
	    do i=1,n
	      t0=0.d0
	      xt=cdabs((xo(i)-xn(i))/xn(i))
	      if(xt .gt. pmax) then 
		!divergent
		write(6,*)'divegent'
		stop
	      end if
	      if(xt>t0) then
		t0=xt
	      end if
	    end do
	    ! xn(i)<-xo(i)
	    do i=1,n
	      xo(i)=xn(i)
	    end do
!	    write(11,"('(',2e17.10,')', 2X,'|', 2X,'(',2e17.10, ')',2X,'|', 2X,'(',2e17.10,')')")(xo(i),i=1,n)
	    write(11,"(I3,$)")k
	    do i=1,n
	      write(11,"(2X,'(',2e17.10,')',4X,$)")xo(i)
	    end do
	    write(11,*)
	    if(t0<tol)then
	      write(6,*)'find solutions'
	      exit
	    end if
	  end do
	  write(6,*)(xo(i),i=1,n)
	end subroutine hw7_gauss_seidel