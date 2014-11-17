!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Davalos, Jadurani M.	!
!	201118835		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cs131mp2
implicit none

! declare variables
	integer, parameter ::  ikind=selected_real_kind(p=18)
	double precision, dimension(4)	::	y_old, y_new  ! size 4 for the 4 equations
	double precision	::	epsilon, h, period, t_max, t_inc, t_in
	double precision	::	F
	integer	::	m
	
	y_old(1) = 1.2d0
    y_old(2) = 0.0d0
    y_old(3) = 0.0d0
    y_old(4) = -1.049358d0
    
!	print *, y_old


	period = 6.19216933
	t_in = 0
    t_max = period*2
	t_inc = period/10
	epsilon = 0.0000001d0
    h = epsilon**(0.25)
	m = 0


	do while( t_in < t_max )
		call RFK45(y_old, y_new, t_in, t_in + t_inc, epsilon, h, m, 4)
		print *, "Y_NEW: ", y_new
        print *, "M: ", m
        m = 0
		print *, "ITERATION ", t_in
    end do

!	!print *, y
!    !print *, "eq 2 or w", F(2, y)

end program cs131mp2


! ++++++++++++++++++++++++++++++++ SBRTN RFK45 ++++++++++++++++++++++++++++++++ !
subroutine RFK45(y_old, y_new, t_in, t_out, epsilon, h, m, n)

   	integer, parameter ::  ikind=selected_real_kind(p=18)
	double precision, dimension(4)	::	y_old, y_new, y_dum, yhat_new, r, alpha_min_arr
	double precision, dimension(4,6)	::	k
	double precision	::	t_in, t_out, h, epsilon, F, alpha
    integer	::	m, p, n, ayos
	

	!print *, "I'm HERE!", epsilon		

	do while( t_in < t_out )

		!print *, "T_IN ", t_in, " T_OUT ",t_out 
		if(t_in + h > t_out) then
			h = t_out - t_in
        end if

! for kp1
		do p=1,n
        	k(p,1) = F(p, y_old)
        end do
		
!for kp2
		do p=1,n
			y_dum(p) = y_old(p) + (1.0d0/4) * h * k(p,1)
        end do

		do p=1,n
			k(p,2) = F(p, y_dum)
        end do

!for kp3
		do p=1,n
			y_dum(p) = y_old(p) + (h/32) * (3*k(p,1) + 9*k(p,2))
        end do

		do p=1,n
			k(p,3) = F(p, y_dum)
        end do

!for kp4
		do p=1,n
			y_dum(p) = y_old(p) + (h/2197) * (  1932*k(p,1) - 7200*k(p,2) + 7296*k(p,3)  )
        end do
        
		do p=1,n
			k(p,4) = F(p, y_dum)
        end do

!for kp5
		do p=1,n
			y_dum(p) = y_old(p) + ( (439.0d0/216) * h * k(p,1) ) - & 
            			& ( 8 * h * k(p,2) ) + ( (3680.0d0/513) * h * k(p,3) ) - ( (845.0d0/4104) * h * k(p,4) )
        end do		

		do p=1,n
			k(p,5) = F(p, y_dum)
        end do

!for kp6
		do p=1,n
			y_dum(p) = y_old(p) - ( (8.0d0/27) * h * k(p,1) ) + ( 2 * h * k(p,2)  ) - &
            			& ( (3544.0d0/2565) * h * k(p,3) ) + ( (1859.0d0/4104) * h * k(p,4) ) - ( (11.0d0/40) * h * k(p,5) )
        end do
        
		do p=1,n
			k(p,6) = F(p, y_dum)
        end do

! computing y_new and yhat_new
		do p=1,n
        	y_new(p) = y_old(p) +   h * ( (25.0d0/216)*k(p,1) + (1408.0d0/2565)*k(p,3) + &
            				& (2197.0d0/4104)*k(p,4) - (1.0d0/5)*k(p,5) )
            yhat_new(p) = y_old(p) +  h*( (16.0d0/135)*k(p,1) + (6656.0d0/12825)*k(p,3) + & 
            				& (28561.0d0/56430)*k(p,4) - (9.0d0/50)*k(p,5) + (2.0d0/55)*k(p,6) )
        end do

		
		!print *, "Y_OLD: ", y_old
		!print *, "Yhat_NEW: ", yhat_new
		!print *, "Y_NEW: ", y_new
        !print *, "I'M H: ", h
! SEE IF INTEGRATION STEP IS OK -> compute r and alpha
		
		do p=1,n
        	r(p) = abs( yhat_new(p) - y_new(p))/h
            !print *, "I'M FUCKING Rp!!", r(p)

			if (r(p) /= 0.0d0) then
            	alpha_min_arr(p)= 0.84*(epsilon/r(p))**(0.25)
            else
				alpha_min_arr(p) = 4
			end if
        end do
        
		alpha = minval(alpha_min_arr)
		!print *, "ARR : ", alpha_min_arr
		!print *, "ALPHA: ", alpha
        !print *, "EPSILON: ", epsilon
! Find new h for failed step or for next step
		
		ayos = 1        
		do p=1,n
			if (r(p) > epsilon) then
            	ayos = 0
            end if
		end do

		if (ayos == 1) then
        	t_in = t_in + h
            y_old = y_new
        end if


		!print *, "OLD H!! ", h
		if(alpha <= 0.1) then
			h = 0.1*h
		else if (alpha >= 4.0d0) then
			h = 4.0d0*h
		else
			h = alpha*h
		end if

		m = m+1
		!print *, "NEW H!! ", h
        if (mod(m,25)==0) then
			print *, "HERE I M!! ", m 
		end if
    end do
	
end subroutine RFK45


function F(p, v)
implicit none

!declare variables
	integer	::	p	!function number
   	integer, parameter ::  ikind=selected_real_kind(p=18)
	double precision, dimension(4)	::	v	!contains y_old.  ignore t na lang; it ain't gonna be used anyway
	double precision	::	F, a, b, p1, p2

	a = 1/82.45
    b = 1 - a
    p1 = sqrt( ( (v(1) + a)**2 + (v(2))**2  )**3 )
    p2 = sqrt( ( (v(1) - b)**2 + (v(2))**2  )**3 )

!	!print *, "hello I'm v", v
	select case (p)
		case (1)
        	F = v(3)
		case (2)
        	F = v(4)
		case (3)
        	F = 2 * v(4) + v(1) - b * ( v(1) + a)/p1 - a * (v(1) - b)/p2
		case (4)
			F = -2 * v(3) + v(2) - b * v(2)/p1 - a * v(2)/p2
    end select

	F=F
return
end function F
