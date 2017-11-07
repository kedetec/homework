! file: hw.f.
! Author: Yazhou.Lee  lee.san.cn@gmail.com


	!MAIN PROGRAM
	PROGRAM HOMEWORK0
	INTEGER NMAX, N, L
	real*8 a, b, q, TOL, x, sa, sb
	
	! Input some initial variables
	write(6,*)'Input a: '
	read(5,*)a
	write(6,*)'Input b: '
	read(5,*)b
	write(6,*)'Input NMAX: '
	read(5,*)NMAX
	write(6,*)'Input TOL: '
	read(5,*)TOL
	write(6,*)'Input N: '
	read(5,*)N

	! write the variables in hw0.txt
	OPEN(UNIT=7, FILE='hw0.txt')
	write(7,*)'Input a = ', a
	write(7,*)'Input b =  ', b
	write(7,*)'Input NMAX = ', NMAX
	write(7,*)'Input TOL = ', TOL
	write(7,*)'Input N = ',N
	
	x = (b - a)/dble(N)
	sa = a
	sb = sa + x
	DO i=1,N
	CALL HASSOLUTION(sa, sb, L)
	if(L == 1) then 
		write(6,*)'Has answer between: ',sa,'--',sb
		write(7,*)'Has answer between: ',sa,'--',sb
		! CALL SUBROUTINE
		CALL BISECT(NMAX, sa, sb, TOL, q)
		write(6,*)'Answer is: ',q
		write(7,*)'Answer is: ',q
	else
		! Do nothing
		write(6,*)'Has no answer between: ',sa,'--',sb
	ENDIF
	sa = sb
	sb = sb + x
	END DO

	CLOSE(7)
	STOP

	END PROGRAM HOMEWORK0 
