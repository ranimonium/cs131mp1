!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Davalos, Jadurani M.	!
!	201118835		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cs131mp1
implicit none

	! declare variables
	real, allocatable, dimension(:,:)	::	A_matrix
	integer					::	t=0		! tth matrix/case; just a counter
	integer					::	n		! size of tth matrix
	integer					::	row, col	! counters in reading contents for A_matrix
	integer, allocatable, dimension(:)	::	p, q		! permutation vector.  p for rows, q for cols
	real					::	det
	real, dimension(3)			::	mat_norm_A, mat_norm_AInv	!contains norms for A and A inverse		
	integer					::	ierror
	character (len=15)			::	inputfile,	outputfile

	! declare functions
	real					::	norm1, normF, normInf

	
	! tell sir about my identity
	print *, 'Submitted by:	Davalos, Jadurani M.	201118835'

	! promt user for name of input file
	print *, 'Enter name of input file: '
	read *,inputfile

	! read input file	
	open(unit=10, file=inputfile,  iostat=ierror)

	! check for error in reading input file
	if ( ierror /= 0 ) then
		print *, "Failed to open file!"
	stop
	end if

	! promt user for name of output file
	print *, 'Enter name of output file: '
	read *,outputfile

	! read output file
	open(unit=99, file=outputfile)

	! start reading matrices
	readfile : do
		
		! read size of A_matrix to allocate
		read(unit=10, fmt=*,  iostat=ierror) n
	
		if (ierror < 0) then
!        	print *,"END OF FILE"
			exit readfile
		end if

		t = t + 1
        
		! allocate spaces for A_matrix, p, q
	        allocate(A_matrix(n,n))
    		allocate(p(n))
    		allocate(q(n))		

		! fill up A_matrix            
		do row = 1,n
			read(10,*) (A_matrix(row,col),col=1,n)
		end do
    
		

		! calculate for mat_normA
        	mat_norm_A(1) = norm1(A_matrix, n)
        	mat_norm_A(2) = normF(A_matrix, n)
        	mat_norm_A(3) = normInf(A_matrix, n)



        	write(99,*)'****** Matrix',t,' ******'
        
		! write contents of A_matrix
        	write(99, *) 'a) Input matrix'
		do row = 1,n
			write(99, '(10ES30.15) ') ( A_matrix(row,col), col=1, n)
		end do
        

        	! compute for det and transformed A
		call GJRM(A_matrix,n,0.0000000001, det, p, q)


		! write transformed A_matrix
	        write(99, *) 'b) Transformed matrix'
	        do row = 1,n
			write(99, '(10ES30.15) ') ( A_matrix(row,col), col=1, n)
		end do

		call unscramble(A_matrix,p,q,n)
        
        	! write vectors p and q
        	write(99, *) 'c) The permutation vectors p and q'
        	write(99, *) '		p = [', p, ']'
        	write(99, *) '		q = [', q, ']'
        
        	! write A_matrix inverse
        	write(99, *) 'd) Inverse matrix'
        	do row = 1,n
			write(99, '(10ES30.15) ') ( A_matrix(row,col), col=1, n)
		end do
        
        	! write determinant
        	write(99, *) 'e) Determinant of matrix'
		write(99, *) '		det = ', det
	        
		! compute norms for A inverse
		mat_norm_AInv(1) = norm1(A_matrix, n)
		mat_norm_AInv(2) = normF(A_matrix, n)
		mat_norm_AInv(3) = normInf(A_matrix, n)
            
		! condition numbers based on 1-, F-, inf-norms
		write(99, *) 'f) Condition numbers'
        	write(99, *) '1-norm based condition number:		', 
		write(99, '(10ES30.15)') mat_norm_A(1)*mat_norm_AInv(1)
		write(99, *) 'F-norm based condition number:		'
		write(99, '(10ES30.15)') mat_norm_A(2)*mat_norm_AInv(2)
	        write(99, *) 'Infinity-norm based condition number:		' 
		write(99, '(10ES30.15)') mat_norm_A(3)*mat_norm_AInv(3)
            
		! deallocate spaces
		deallocate(A_matrix)
		deallocate(p)
		deallocate(q)
        
		write (99,*)'********************************'
		write (99,*)''
        
    end do readfile !end of reading matrices


    close(unit=10) ! close input file
    close(unit=99) ! close outputfile

	print *, "Computation and write successful.  See ", outputfile
    
end program cs131mp1

! ++++++++++++++++++++++++++++++++ SBRTN GJRM ++++++++++++++++++++++++++++++++ !

subroutine GJRM(A_matrix, n, e, det, p, q)
implicit none

	integer				::		n
	real, dimension(n,n)		::		A_matrix
	real				::		a_max, swappie, a_rowpivot
	real				::		e
	real				::		det
	integer				::		row, col, pivot, temp_row, temp_col, count		! i=row, j=col, pivot=k, r=temp,
	integer, dimension(n)		::		p, q	! permutation vector.  p for row, q for col

! initializations
	det	=	1
	do count=1,n
		p(count)=count
	q(count)=count
	end do


! Gauss-Jordan reduction
	do pivot=1,n

		! search pivot element, total pivoting
	    	a_max = A_matrix(pivot, pivot)
			temp_row = pivot
	        temp_col = pivot
	        do row=pivot,n
			do col=pivot,n
	            	if (abs(A_matrix(row, col)) > abs(a_max)) then
	                	a_max = A_matrix(row, col)
				temp_row = row
				temp_col = col
	                end if
			end do
		end do
	
		! test for singularity
	    	if ( abs(a_max) < e)  then
			det = 0
			return
		else
	        	det = det * a_max
		end if
	
		! interchange rows, if necessary
	    	if ( temp_row /= pivot ) then
	
	        	do col=1,n
	            		swappie = A_matrix(pivot,col)
	                	A_matrix(pivot,col) =	A_matrix(temp_row, col)
	                	A_matrix(temp_row,col)	= 	swappie
			end do
	
			swappie = p(pivot)
			p(pivot) = p(temp_row)
			p(temp_row) = swappie
	
			det = (-1)*det  !because row interchange yields the negative of current det
		end if
	
		! interchange cols, if necessary
	    	if ( temp_col /= pivot ) then
	
	        	do row=1,n
	        	    	swappie = A_matrix(row,pivot)
	                	A_matrix(row,pivot) =	A_matrix(row,temp_col)
	                	A_matrix(row,temp_col)	= 	swappie
			end do
	
			swappie = q(pivot)
			q(pivot) = q(temp_col)
			q(temp_col) = swappie
	
			det = (-1)*det  !because cols interchange yields the negative of current det
		end if
	
		! normalize pivot row
	    	do col=1,n
	        	A_matrix(pivot,col) = A_matrix(pivot,col)/a_max
	        end do
	
	        A_matrix(pivot, pivot) = 1.0/a_max !i don't get why
	
	
		! reduce pivot column to unit vector e
	    	do row=1, n
	        	if( row /= pivot ) then
	        	    	a_rowpivot = A_matrix(row,pivot)
	                	do col=1,n
					A_matrix(row,col) = A_matrix(row,col) - a_rowpivot*A_matrix(pivot,col)
				end do
	            		A_matrix(row,pivot) = (-1)*(a_rowpivot/a_max)
	            	end if	
		end do
	end do ! Gauss-Jordan reduction

end subroutine GJRM

! ++++++++++++++++++++++++++++++++ SBRTN UNSCRAMBLE ++++++++++++++++++++++++++++++++ !	

subroutine unscramble(A_matrix, p, q, n)
implicit none
	integer				::		n
	integer				::		row, col
	real, dimension(n,n)		::		A_matrix
	integer, dimension(n)		::		p, q		! permutation vector.  p for row, q for col
	real, dimension(n)		::		t		! permutation vector.  t for swapping

	! unscramble
    	do row=1,n
        	do col=1,n
            		t(p(col)) = A_matrix(row,col)
            	end do
		do col=1,n
			A_matrix(row,col) = t(col)
		end do
        end do

	do col=1,n
        	do row=1,n
            		t(q(row)) = A_matrix(row,col)
		end do

		do row=1,n
			A_matrix(row,col) = t(row)
		end do
	end do

end subroutine unscramble

! ++++++++++++++++++++++++++++++++ FCN NORM1 ++++++++++++++++++++++++++++++++ !	

function norm1(M_matrix, n)
implicit none
	integer			::	n, row, col
	real, dimension(n,n)	::	M_matrix
	real			::	norm1, sum

	norm1=0.0 ! just to init some value
    
	do col=1,n
		sum = 0.0
		do row=1,n
        		sum = sum + abs(M_matrix(row, col))
        	end do
		if(col == 1) then
        		norm1=sum
        	else if(sum>norm1) then
			norm1=sum
		end if
	end do

	norm1=norm1
return
end function norm1


! ++++++++++++++++++++++++++++++++ FCN NORMF ++++++++++++++++++++++++++++++++ !	


function normF(M_matrix, n)
implicit none
	integer			::	n, row, col
	real, dimension(n,n)	::	M_matrix
	real			::	normF, sum

	normF=0.0 ! just to init some value
	sum = 0.0
    
	do row=1,n
		do col=1,n
        		sum = sum + (M_matrix(row, col))**2
		end do
	end do

    normF=sqrt(sum)
return
end function normF


! ++++++++++++++++++++++++++++++++ FCN NORMINF ++++++++++++++++++++++++++++++++ !	

function normInf(M_matrix, n)
implicit none
	integer			::	n, row, col
	real, dimension(n,n)	::	M_matrix
	real			::	normInf, sum

	normInf=0.0 ! just to init some value
    
	do row=1,n
    		sum = 0.0
		do col=1,n
        		sum = sum + abs(M_matrix(row, col))
        	end do
		if(row == 1) then
        		normInf=sum
        	else if(sum>normInf) then
			normInf=sum
		end if
	end do

	normInf=normInf
return
end function normInf

