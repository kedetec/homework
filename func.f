! file: hw0-bisect.f. 
! Author: Yazhou.Lee  lee.san.cn@gmail.com


	! FUNCTION
	FUNCTION ftan(x)
	real*8 ftan, x
	ftan = dtan(x)
	return
	END FUNCTION

	! FUNCTION
	FUNCTION poly(x)
	real*8 poly, x
	poly = x*x*x + 4*x*x - 10.d0
	return
	END FUNCTION

	! the interl between inter_l and inter_r is small enough
	! answer: 0 -- has no solution
	!         1 -- has solution
	SUBROUTINE HASSOLUTION(a, b, l)
	INTEGER l
	real*8 a, b, t, c, d
	answer = 1

	c = ftan(a)
	d = ftan(b)
	IF(c*d>0)then
		l = 0
		return 
	ENDIF
	t = ftan((a+b)/2.d0)

	if(t*a > 0.d0)then
		if(dabs(t) .gt. dabs(c))then 
		l = 0
		else 
		l = 1
		end if
	else
		if(dabs(t) .gt. dabs(d))then 
		l = 0
		else 
		l = 1
		end if
	endif
	return
	END SUBROUTINE

	!SUBROUTINE, binary search
	! it is supposed that the equation has only one solution 
	! between inter_l and inter_r
	SUBROUTINE BISECT(NMAX, a, b, TOL, q)
	INTEGER NMAX
	real*8 a, b, q, TOL, t

	!check if has no answer
	if(ftan(a)*ftan(b)>0) then
		write(6,*)'No answer.'
		write(7,*)'No answer'
		return
	endif

	DO i=1, NMAX
		t = (a+b)/2.d0
		if(dabs(b-a) < TOL) then !
			l = t
			return
		ENDIF
		
		if(ftan(t)*ftan(b) > 0.d0) then
			b = t
		else 
			a = t
		end if

	ENDDO
	write(6,*)'Iteration times are over .',NMAX
	answer = t
	return
	END SUBROUTINE

	!SUNROUTINE
	SUBROUTINE CHOOSEFUNC(M, f)
		INTEGER M
		character(len=50) squit, stan, spoly, sindicator
		PROCEDURE f
		OPEN(UNIT=8, FILE='hw0-choosefun.txt')
		
		squit = "0 quit the program"
		stan = '1 tan(x).'
		spoly= '2 Polynomial x*x*x + 4*x*x - 10.d0'
		sindicator= 'Please choose your operation by keyboard.'
		
!		write(8,*)'0 quit the program'
!		write(8,*)'1 tan(x).'
!		write(8,*)'2 Polynomial:  x*x*x + 4*x*x - 10.d0'			
!		write(8,*)'Please choose your operation by keyboard.'

		DO
			write(6,*)'0 quit the program'
			write(6,*)'1 tan(x).'
			write(6,*)'2 Polynomial:  x*x*x + 4*x*x - 10.d0'			
			write(6,*)'Please choose your operation by keyboard.'

			read(5,*)M
			if(M==0)then
				write(6,*)'Byebye'
				STOP
			else if(M==1)then
				write(6,*)'tan is choosen'
				f = ftan
				EXIT
			else if(M==2)then
				write(6,*)'poly is choosen'
				f = poly
				EXIT
			else
			END IF
		END DO
		return
	END SUBROUTINE