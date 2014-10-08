program cs131mp1
implicit none

	! declare variables
	real, allocatable, dimension(:,:)	::	A_matrix
	integer	::	test_size	! number of test matrices
    integer	::	t			! tth matrix/case; just a counter
    integer ::	n			! size of tth matrix
	integer	::	row, col	! counters in reading contents for A_matrix



	! read input file	
	open(10, file='tin.txt')

	! read output file
    open(99, file='tout.txt')

	! read number of test matrices
	read(10, *) test_size

	! start reading matrices
	do t=1,test_size
		
		! read size of A_matrix to allocate
		read(10, *) n
		! allocate nxn spaces for A_matrix
        allocate(A_matrix(n,n))
    
		! fill up A_matrix            
		do row = 1,n
!	    	do col = 1,n
!			READ(11,*) (a(row,col),col=1,max_cols)
			read(10,*) (A_matrix(row,col),col=1,n)
 !       	read(10,*) A_matrix(row,col)
!			end do
		end do
    

		! LET THE ALGO BEGIN

!		call GJRM(A_matrix,n,0.0000000001)

		! END OF ALGO


		! write contents of A_matrix
		do row = 1,n
!	   		do col = 1,n
        		write(99, '(10ES30.10) ') ( A_matrix(row,col), col=1, n)
!			end do
		end do
	PRINT *, A_matrix(1,3)

	    deallocate(A_matrix)
	write (99,*)'END OF MATRIX'
        
    end do !end of reading matrices
end program cs131mp1

! ++++++++++++++++++++++++++++++++ !

subroutine GJRM(A_matrix, n, e)

	integer								::		n
	real, dimension(n,n)				::		A_matrix
    real								::		a_max, swappie
	real								::		e
	real								::		det
	integer				 				::		row, col, pivot, temp, count			! i=row, j=col, pivot=k, r=temp,
    integer, dimension(n)	::		p,	q		! permutation vectors.  p for places, q for swapping


! initializations
	det	=	1
	do count=1,n
		p(count)=count
    end do


! Gauss-Jordan reduction
	do pivot=1,n

	! search pivot element
    	a_max = A_matrix(pivot, pivot)
        temp = pivot
        do count=pivot, n
			if (abs(A_matrix(count, pivot)) > abs(a_max)) then
				a_max = A_matrix(count, pivot)
                temp = count
			end if
		end do

	! test for singularity
    	if ( abs(a_max) < e)  then
        	det = 0
            return
		else
        	det = det * a_max
		end if

	! interchange rows, if necessary

    	if ( temp /= pivot ) then
        	do col=1,n
                print *, "BEFORE SWAP: "
                print *, "A_matrix(pivot,col): ", A_matrix(pivot,col)
                print *, "A_matrix(temp,col): ", A_matrix(temp,col)
            	swappie = A_matrix(pivot,col)
                A_matrix(pivot,col) =	A_matrix(temp, col)
                A_matrix(temp,col)	= 	swappie
                print *, "AFTER SWAP: "
                print *, "A_matrix(pivot,col): ", A_matrix(pivot,col)
                print *, "A_matrix(temp,col): ", A_matrix(temp,col)
				print *, ""
			end do


			swappie = p(pivot)
            p(pivot) = p(temp)
            p(temp) = swappie

            det = (-1)*det  !because row interchange yields the negative of current det
		end if

	! normalize pivot row
    	do col=1,n
        	A_matrix(pivot,col) = A_matrix(pivot,col)/a_max
        end do
		print *, "Akk=", A_matrix(pivot, pivot)
        A_matrix(pivot, pivot) = 1.0/a_max !i don't get why
		print *, "Akk=", A_matrix(pivot, pivot)

        
	end do ! Gauss-Jordan reduction


!print *, n

!		do row = 1,n
!	   		do col = 1,n
!        		print *, A_matrix(row,col)
!			end do
!		end do

	


end subroutine GJRM