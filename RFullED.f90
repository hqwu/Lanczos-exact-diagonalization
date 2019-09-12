
!======================================================================
!	Lanczos exact diagonalization
!======================================================================
    
subroutine RFullED

    use MyPrec
    use MyHmlt
    implicit none
    integer i, k, err
    integer(8) timst, timed
    real(rp), allocatable :: EgnVal(:)
    real(rp), allocatable :: mfv(:)
    character(64) CTmp

    call TimingBgn(timst)
    call PrintMsg("Start Full Exact Diagonalization")
    
    allocate(ALLENG(2**NSite))
    ALLENG = 0._rp
    NEGTH = 0
    
    !if k = 1, 1, you will get only the eigen information of Mz=0 block
    do k = 1, NHBlk
    
        NeUp = BlkNeUp(k)
        DmHb = BlkHDim(k)

        allocate(BsWf(DmHb))
        allocate(EgnVal(DmHb))
        allocate(RHMat(DmHb, DmHb))
        EgnVal = 0._rp
        RHMat  = 0._rp
    
        call BlockBasis
        call HmltElmt_Full

        call PrintMsg("Calculating Eigenvalues and Eigenvectors")
        allocate(mfv(DmHb))
        !use accumulating orthogonal similarity transformations 
        !to get tridiagonal matrix and then use QL algorithm to diagonalize it.
        call RS(DmHb, DmHb, RHMat, EgnVal, RHMat, mfv, err)
        deallocate(mfv)
        
        !Gather all the eigenvalues from all blocks.
        do i = 1, DmHb
            ALLENG(NEGTH + i) = EgnVal(i)
        end do
        NEGTH = NEGTH + DmHb

        !output the eigenvalues of the current block.
        write(CTmp, "('Block', I2.2, '.txt')") k
        open(20, file = "FullED_EE_" // trim(CTmp))
        do i = 1, DmHb
            write(20, "(2I5, ES<rp_ow>.<rp_od>)") k, i, EgnVal(i)
        end do
        close(20)
    
        deallocate(BsWf)
        deallocate(EgnVal)
        deallocate(RHMat)
    
    end do

    !Sort the eigenvalues from different blocks into 
    !ascending order and then calculate the specific heat
    call QuickSortR(-1, 1, 2 ** NSite, ALLENG)
    call SpecificHeat
    
    deallocate(ALLENG)
    
    call PrintMsg("End Full Exact Diagonalization")
    call TimingEnd(timst, timed)

end subroutine

!======================================================================

subroutine HmltElmt_Full

    use MyPrec
    use MyHmlt
    implicit none
    integer Cfg
    integer i, j, n, Lct
    
    call PrintMsg("Start generating Hamiltionian matrix")
    
    RHMat = 0._rp
    do n = 1, DmHb
    Cfg = BsWf(n)
    do i = 1, NSite
    do j = i + 1, NSite
        if(abs(JzzMat(i, j)) > rp_zero) then !SzSz interaction
            RHMat(n, n) = RHMat(n, n) + JzzMat(i, j) * &
            (dble(ibits(Cfg, i-1, 1)) - 0.5_rp) * (dble(ibits(Cfg, j-1, 1)) - 0.5_rp)
        end if
        if(abs(JpmMat(i, j)) > rp_zero .and. &
        ibits(Cfg, i-1, 1) + ibits(Cfg, j-1, 1) == 1) then !S+S- interaction
            !use the bisection method to find the nonzero matrix elements.
            call BasisIndex(ibchng(ibchng(Cfg, i-1), j-1), DmHb, BsWf, Lct)
            RHMat(Lct, n) = RHMat(Lct, n) + JpmMat(i, j) * 0.5_rp
        end if
    end do
    end do
    end do
    
    !open(20, file = "HmltMatrix.txt")
    !do i = 1, DmHb
    !    write(20, "(<DmHb>F12.8)") (RHMat(i, j), j = 1, DmHb)
    !end do
    !close(20)
    
    call PrintMsg("End  generating Hamiltionian matrix")
    
end subroutine

!======================================================================

subroutine SpecificHeat

    use MyPrec
    use MyHmlt
    implicit none
    integer i, k
    real(rp) Temp, TBet, Cv
    real(rp) prfc, Cv2
    
    open(60, file = "FullED_ET_SpeCv.txt")
    do k = 1, 10000 !temperature slices
        Temp = 0.01_rp + 0.001_rp * dble(k)
        TBet = 1._rp / Temp
        Cv   = 0._rp
        Cv2  = 0._rp
        prfc = 0._rp
        do i = 1, NEGTH
            if(exp(-TBet * (ALLENG(i) - ALLENG(1))) > 1.E-14_rp) then
                Cv  = Cv  + ALLENG(i) * exp(-TBet * (ALLENG(i) - ALLENG(1)))
                Cv2 = Cv2 + (ALLENG(i) ** 2) * exp(-TBet * (ALLENG(i) - ALLENG(1)))
                prfc = prfc + exp(-TBet * (ALLENG(i) - ALLENG(1))) !partition function
            else
                exit
            end if
        end do
        Cv = (Cv2 / prfc - (Cv / prfc)** 2) / (dble(NSite) * (Temp ** 2))
        write(60, "(2ES<rp_ow>.<rp_od>)") Temp, Cv
    end do
    close(60)    
    
end subroutine

!======================================================================
