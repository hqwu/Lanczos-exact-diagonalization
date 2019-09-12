
!======================================================================

subroutine HmltElmt

    use MyPrec
    use MyHmlt
    implicit none
    integer i, j, swt, err

    integer IHTL, IHT1
    integer,  allocatable :: IHIL(:)  !store the column indices of nonzero matrix elements in one row
    real(rp), allocatable :: RHEL(:)  !store the nonzero matrix elements in one row

    !temporary variables, we don't know the exact number of nonzero matrix elements by advance.
    !smart but cumbersome
    integer,  allocatable :: IHI1(:), IHI2(:)
    real(rp), allocatable :: RHE1(:), RHE2(:)
    
    call PrintMsg("Start generating Hmlt elements")
    
    IHT = 0    !counting the number of nonzero matrix elements.
    IHT1 = IHM !you can guess a number which closes to the actual number of nonzero matrix elements.
    swt  = 1

    allocate(RHEL(IHM), RHE1(IHT1), &
             IHIL(IHM), IHI1(IHT1), stat = err)
    if(err .ne. 0) then
        write(ErrMsg, "('allocating error in HmltElmt')")
        call ErrOut(ErrMsg)
    end if

    do i = 1, DmHb
        
        IHTL = 0
        call HmltDiagElmt   (i, BsWf(i), IHTL, IHIL, RHEL) ! SzSz interaction
        call HmltOffDiagElmt(i, BsWf(i), IHTL, IHIL, RHEL) ! S+S- interaction
        call QuickSortIR(1, 1, IHTL, IHIL, RHEL) !ascending all the nonzero matrix elements according to their column indices.
        call MergeIndex(IHTL, IHIL, RHEL) !merge the matrix elements with same column indices or remove the matrix elements with zero value.
        
        IHL(i) = IHT + 1    !starting index of nonzero elements which belong to the same row i
        IHR(i) = IHT + IHTL !ending index of nonzero elements which belong to the same row i
        IHT    = IHT + IHTL !counting the total number of nonzero matrix elements
        
        if(IHT1 >= IHT) then
            if(swt == 1) then
                IHI1(IHL(i) : IHR(i)) = IHIL(1 : IHTL)
                RHE1(IHL(i) : IHR(i)) = RHEL(1 : IHTL)
            else
                IHI2(IHL(i) : IHR(i)) = IHIL(1 : IHTL)
                RHE2(IHL(i) : IHR(i)) = RHEL(1 : IHTL)
            end if
        else !if IHT1 < IHT, we need to allocate a larger memory to place the nonzero elements.
            do while(IHT1 < IHT)
                IHT1 = IHT1 * 2
            end do
            if(swt == 1) then
                allocate(IHI2(IHT1))
                allocate(RHE2(IHT1))
                swt = 2
                IHI2(1 : IHL(i) - 1) = IHI1(1 : IHL(i) - 1)
                RHE2(1 : IHL(i) - 1) = RHE1(1 : IHL(i) - 1)
                deallocate(IHI1)
                deallocate(RHE1)
                IHI2(IHL(i) : IHR(i)) = IHIL(1 : IHTL)
                RHE2(IHL(i) : IHR(i)) = RHEL(1 : IHTL)
            else if(swt == 2) then
                allocate(IHI1(IHT1))
                allocate(RHE1(IHT1))
                swt = 1
                IHI1(1 : IHL(i) - 1) = IHI2(1 : IHL(i) - 1)
                RHE1(1 : IHL(i) - 1) = RHE2(1 : IHL(i) - 1)
                deallocate(IHI2)
                deallocate(RHE2)
                IHI1(IHL(i) : IHR(i)) = IHIL(1 : IHTL)
                RHE1(IHL(i) : IHR(i)) = RHEL(1 : IHTL)
            end if
        end if
        
    end do
    
    deallocate(IHIL, RHEL, stat = err)
    allocate(IHI(IHT), RHE(IHT), stat = err)
    if(swt == 1) then
        IHI(1 : IHT) = IHI1(1 : IHT)
        RHE(1 : IHT) = RHE1(1 : IHT)
        deallocate(IHI1, RHE1, stat = err)
    else
        IHI(1 : IHT) = IHI2(1 : IHT)
        RHE(1 : IHT) = RHE2(1 : IHT)
        deallocate(IHI2, RHE2, stat = err)
    end if
    if(err .ne. 0) then
        write(ErrMsg, "('deallocating error in HmltElmt')")
        call ErrOut(ErrMsg)
    end if
    
    open(21, file = trim(Flog), access = "append")
    write(21, "(A28, I12)") "# Nonzero matrix elements :", IHT
    close(21)
    
    call PrintMsg("End   generating Hmlt matrix elements")
    
end subroutine

!======================================================================

subroutine HmltDiagElmt(Idx, Cfg, IHTL, IHIL, RHEL)

    use MyPrec
    use MyHmlt
    implicit none
    integer i, j, Idx, Cfg
    integer IHIL(IHM), IHTL
    real(rp) DgEm, RHEL(IHM)
    
    DgEm = 0._rp
    do i = 1, NSite
    do j = i + 1, NSite
        if(abs(JzzMat(i, j)) > rp_zero) then
            DgEm = DgEm + JzzMat(i, j) &
                        * (dble(ibits(Cfg, i-1, 1)) - 0.5_rp) &
                        * (dble(ibits(Cfg, j-1, 1)) - 0.5_rp)
        end if
    end do
    end do

    if(abs(DgEm) > rp_zero) then
        IHTL = IHTL + 1
        RHEL(1) = DgEm
        IHIL(1) = Idx
    end if
    
end subroutine

!======================================================================

subroutine HmltOffDiagElmt(Idx, Cfg, IHTL, IHIL, RHEL)

    use MyPrec
    use MyHmlt
    implicit none
    integer i, j, Idx, Lct, Cfg
    integer IHTL, IHIL(IHM)
    real(rp) HopElmt, RHEL(IHM)
        
    do i = 1, NSite
    do j = i + 1, NSite
        HopElmt = JpmMat(i, j) * 0.5_rp
        if(abs(HopElmt) > rp_zero .and. &
        ibits(Cfg, i-1, 1) + ibits(Cfg, j-1, 1) == 1) then
            !flip two spins and find the representative state
            call BasisIndex(ibchng(ibchng(Cfg, i-1), j-1), DmHb, BsWf, Lct)
            if(Lct > 0) then
                IHTL = IHTL + 1
                IHIL(IHTL) = Lct
                RHEL(IHTL) = HopElmt
            end if
        end if
    end do
    end do
    
end subroutine

!======================================================================

subroutine BasisIndex(Cfg, Dm, Bs, NIdx)
!bisection search

    use MyPrec
    implicit none
    integer l, r, m, Dm, NIdx
    integer Bs(Dm), Cfg

    l = 1
    r = Dm   
    do while(l <= r)
    m = (l + r) / 2
    if(Bs(m) < Cfg) then
        l = m + 1
    else if(Bs(m) > Cfg) then
        r = m - 1
    else
        exit
    end if
    end do

    if(l <= r) then
        NIdx = m
    else
        NIdx = -1
    end if

end subroutine

!======================================================================

subroutine IMySwap(a, b)

    implicit none
    integer a, b, c
    
    c = a
    a = b
    b = c
    
end subroutine

subroutine RMySwap(a, b)

    use MyPrec
    implicit none
    real(rp) a, b, c
    
    c = a
    a = b
    b = c
    
end subroutine
    
!==================================================================================================
    
subroutine MergeIndex(IHTL, IHIL, IHEL)

	use MyPrec
	implicit none
	integer IHTL, IHIL(IHTL), i, j
	real(rp) IHEL(IHTL)

	if (IHTL <= 0) return

    !merge the elements with the same column indices
	j = 1
	do i = 2, IHTL
	if (IHIL(j) == IHIL(i)) then
		IHEL(j) = IHEL(j) + IHEL(i)
	else
		j = j + 1
		IHEL(j) = IHEL(i)
		IHIL(j) = IHIL(i)
	end if
	end do

	IHTL = j

    !remove zeros
	j = 0
	do i = 1, IHTL
	if (abs(IHEL(i)) > rp_zero) then
		j = j + 1
		IHEL(j) = IHEL(i)
		IHIL(j) = IHIL(i)
	end if
	end do

	IHTL = j

end subroutine
    
!==================================================================================================
