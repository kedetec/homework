! file: hw0-bisect.f. 
! Author: Yazhou.Lee  lee.san.cn@gmail.com

	! FUNCTION
	FUNCTION ftan(x)
	real*8 ftan, x
	ftan = tan(x)
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
	SUBROUTINE HASSOLUTION(inter_l, inter_r, answer)
	INTEGER answer
	real*8 inter_l, inter_r, t, l, r
	answer = 1

	l = ftan(inter_l)
	r = ftan(inter_r)
	IF(l*r>0)then
		answer = 0
		return 
	ENDIF
	t = ftan((inter_l+inter_r)/2.d0)

	if(t*l > 0.d0)then
		if(dabs(t) .gt. dabs(l))then 
		answer = 0
		else 
		answer = 1
		end if
	else
		if(dabs(t) .gt. dabs(r))then 
		answer = 0
		else 
		answer = 1
		end if
	endif
	return
	END SUBROUTINE

	!SUBROUTINE, binary search
	! it is supposed that the equation has only one solution 
	! between inter_l and inter_r
	SUBROUTINE BISECT(NMAX, inter_l, inter_r, TOL, answer)
	INTEGER NMAX
	real*8 inter_l, inter_r, answer, TOL, t

	!check if has no answer
	if(ftan(inter_l)*ftan(inter_r)>0) then
		write(6,*)'No answer.'
		write(7,*)'No answer'
		return
	endif

	DO i=1, NMAX
		t = (inter_l+inter_r)/2.d0
		if(dabs(inter_r-inter_l) < TOL) then !
			answer = t
			return
		ENDIF
		
		if(ftan(t)*ftan(inter_r) > 0.d0) then
			inter_r = t
		else 
			inter_l = t
		end if

	ENDDO
	write(6,*)'Iteration times are over .',NMAX
	answer = t
	return
	END SUBROUTINE
