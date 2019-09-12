
!==================================================================================================
!	Lanczos algorithm to find some lowest eigenstates. Real version.    
!==================================================================================================

subroutine Lanczos(ntot, nev, EgnVal, EgnVec)
!   ntot is the dimension of Hilbert space
!	nev  is the number of eigenstates.

	use MyPrec
	implicit none
	integer ntot, nev, ncv
	real(rp) EgnVal(nev)       ! eigenvalues
    real(rp) EgnVec(ntot, nev) ! eigenvectors
	
	integer p, q, i, j, k, err
	integer TSKM			! the maximum size of Krylov space
	integer TSKR			! the convergent size of Krylov space
	integer EgnConv

    real(rp) bet
	real(rp), allocatable :: EV(:, :)		! unitary transformation matrix of tridiagonalized matrix
	real(rp), allocatable :: ER(:)			! relative error
	real(rp), allocatable :: HB(:), HD(:)   ! diagonal elements / eigenvalues
	real(rp), allocatable :: HC(:), HE(:)   ! sub-diagonal elements
    real(rp), allocatable :: GE(:)
	real(rp), allocatable :: LV(:, :)
	real(rp), allocatable :: HP(:)
	real(rp), allocatable :: FV(:)
    
    ncv = 999

	allocate(LV(ntot, 2)) !cost memory
	allocate(HP(ntot))    !cost memory
	allocate(FV(ntot))    !cost memory
	allocate(HB(ncv))
	allocate(HC(ncv))
	allocate(HD(ncv))
	allocate(HE(ncv))
    allocate(GE(ncv))
    allocate(ER(nev))

    do k = 1, nev
    
	    call Random_Number(FV)
        call RLancOrth(ntot, FV, EgnVec, k - 1)

        TSKM = ncv  ! maximum size of the Krylov space
        if (TSKM > ntot - k + 1) TSKM = ntot - k + 1

	    LV(:, 1) = FV    !vk
	    LV(:, 2) = 0._rp !vk+1
	    p = 1
	    q = 2
	    bet = 0._rp
	    EgnConv = 0

	    do i = 1, TSKM

		    call RMatVecProd(ntot, HP, LV(:, p))  !H|v_n>
		    HE(i) = bet
		    HP = HP - HE(i) * LV(:, q) !H|v_n>-beta*|v_{n-1}>

		    HD(i) = Dot_Product(LV(:, p), HP) !a_n=<v_n|H|v_n>

		    LV(:, q) = HP - HD(i) * LV(:, p) !|v_{n+1}>=H|v_n>-a_n|v_n>-beta*|v_{n-1}>
		    !LV(:, q) = LV(:, q) - Dot_Product(LV(:, p), LV(:, q)) * LV(:, p) !Reorthogonalize |v_{n+1}> to |v_n>, don't need to do that.
		    bet = sqrt(Dot_Product(LV(:, q), LV(:, q))) !b_n
		    LV(:, q) = LV(:, q) * (1._rp / (bet + 1.E-32_rp)) !normalize |v_{n+1}>
            call RLancOrth(ntot, LV(:, q), EgnVec, k - 1) !Reorthonormalize |v_{n+1}> to lowest eigenstates when we are targeting excited states.

		    TSKR = i

		    HB = HD
		    HC = HE
            allocate(EV(TSKR, TSKR))
	        EV = 0._rp
	        do j = 1, TSKR
		        EV(j, j) = 1._rp
            end do
            !finds the eigenvalues and eigenvectors of a SYMMETRIC TRIDIAGONAL matrix by the QL method
            call TQL2(TSKR, TSKR, HB, HC, EV, err)

		    if (err > 0) then
			    write(ErrMsg, "('TQL2 ERR > 0, ERR = ', I11)") err
			    call ErrOut(ErrMsg)
            end if

            GE(TSKR) = HB(1) !only take care of the lowest energy because we target the eigenstates one by one.

            !convergent criteria of the lowest eigenvalue.
            if(TSKR > 1) then
                ER(k) = abs(GE(TSKR) - GE(TSKR - 1)) / abs(GE(TSKR))
                if(ER(k) < 1.E-14_rp) EgnConv = 1
                open(21, file = trim(Flog), access = "append")
		        write(21, "('[', I3, ']', (I5), 3ES<rp_ow>.<rp_od>)") k, i, HB(1), ER(k), HB(2)
                close(21)
                if(abs(HB(2)-HB(1)) < rp_zero) then
                    write(ErrMsg, "('orthogonal lost in the [', I3, ']th eigenvalue')") k
                    call ErrOut(ErrMsg)
                end if
            end if
        
            deallocate(EV)

		    if (EgnConv == 1 .or. TSKR == ntot) exit

            !a simple trick to switch between |v_n> and |v_{n-1}>, and we don't need to store all the Krylov vectors.
            p = 3 - p
		    q = 3 - q

        end do

        !now we know the number of steps (TSKR) to get convergent eigenvalues, the next step is to calculate the corresponding eigenvector.
	    allocate(EV(TSKR, TSKR))
	    EV = 0._rp
	    do j = 1, TSKR
		    EV(j, j) = 1._rp
	    end do
	    HB = HD
	    HC = HE
	    call TQL2(TSKR, TSKR, HB, HC, EV, err)
        if (err > 0) then
		    write(ErrMsg, "('TQL2 ERR > 0, ERR = ', I11)") err
		    call ErrOut(ErrMsg)
        end if
        !EV is the unitary matrix which can transform the tridiagonal matrix to a diagonal matrix.

        !in the following, we want to get the corresponding eigenvector.
        !Since we haven't stored all the Krylov vectors, we need to redo the Lanczos iteration with the same |v_0>.
	    LV(:, 1) = FV
	    LV(:, 2) = 0._rp
	    p = 1
	    q = 2
	    bet = 0._rp

        EgnVal(k) = HB(1)
  
	    do i = 1, TSKR
  
            EgnVec(:, k) = EgnVec(:, k) + EV(i, 1) * LV(:, p) !the columns of V*U stores the eigenvectors
  
		    call RMatVecProd(ntot, HP, LV(:, p))
		    HE(i) = bet
		    HP = HP - HE(i) * LV(:, q)
  
		    HD(i) = Dot_Product(LV(:, p), HP)
  
		    LV(:, q) = HP - HD(i) * LV(:, p)
		    !LV(:, q) = LV(:, q) - Dot_Product(LV(:, p), LV(:, q)) * LV(:, p)
		    bet = sqrt(Dot_Product(LV(:, q), LV(:, q)))
		    LV(:, q) = LV(:, q) * (1._rp / (bet + 1.E-32_rp))
            call RLancOrth(ntot, LV(:, q), EgnVec, k - 1)
  
		    p = 3 - p
		    q = 3 - q
  
        end do
    
        EgnVec(:, k) = EgnVec(:, k) * (1._rp / sqrt(Dot_Product(EgnVec(:, k), EgnVec(:, k)))) !normalization of eigenvectors
  
	    deallocate(EV)
    
    end do

	deallocate(LV)
	deallocate(HP)
	deallocate(FV)
	deallocate(HB)
	deallocate(HC)
	deallocate(HD)
	deallocate(HE)
    deallocate(GE)
	deallocate(ER)

end subroutine

!==================================================================================================
    
subroutine RMatVecProd(ntot, y, x)
! y = H*x, sparse matrix-vector product

    use MyPrec
    use MyHmlt
	implicit none

	integer ntot, i, j
	real(rp) y(ntot), x(ntot)

    do i = 1, ntot
		y(i) = 0._rp
		do j = IHL(i), IHR(i)
            y(i) = y(i) + RHE(j) * x(IHI(j))
        end do
	end do

end subroutine

!===================================================================================================================

subroutine RLancOrth(ntot, FV, EV, NEO)
!Orthonormalization of Krylov vectors with lowest eigenstates of H which have been calculated. 
!if NEO = 0, it only do the normalization.

	use MyPrec
	implicit none
	integer NEO, i, ntot
	real(rp) FV(ntot), EV(ntot, NEO)

	do i = 1, NEO
		FV = FV - EV(:, i) * Dot_Product(EV(:, i), FV)
	end do
	FV = FV / sqrt(Dot_Product(FV, FV))

end subroutine

!==================================================================================================
