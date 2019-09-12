
!======================================================================
!	Lanczos exact diagonalization
!======================================================================
    
subroutine RLancED

    use MyPrec
    use MyHmlt
    implicit none
    integer i, k, NEV
    integer(8) timst, timed
    real(rp), allocatable :: EgnVal(:)
    real(rp), allocatable :: EgnVec(:, :)
    character(64) CTmp

    call TimingBgn(timst)
    call PrintMsg("Start Lanczos Exact Diagonalization")

    GGSE = 1.E+12_rp
    NDeg = 0

    !if k = 1, 1, you can get only the eigen information of Mz=0 block
    do k = 1, NHBlk
    
        NeUp = BlkNeUp(k)
        DmHb = BlkHDim(k)

        allocate(BsWf(DmHb))

        call BlockBasis
        call MySMSFInit(DmHb)
        call HmltElmt
    
        NEV = 6 !the maximum number of states that you want to calculate
        if(NEV > DmHb) NEV = DmHb

        allocate(EgnVal(NEV))
        allocate(EgnVec(DmHb, NEV))
        EgnVal = 0._rp
        EgnVec = 0._rp
    
        call PrintMsg("Calculating Eigenvalues and Eigenvectors")
        call Lanczos(DmHb, NEV, EgnVal, EgnVec)
        call MySMSFFinl
    
        if(GGSE > EgnVal(1) .and. abs(GGSE - EgnVal(1)) > rp_zero) then
            !the true ground state may not be in the Mz=0 block.
            GGSE = EgnVal(1)
            NDeg = 0
            SpnET = 0._rp
        end if

        do i = 1, NEV
            if(abs(EgnVal(i) - GGSE) / abs(GGSE) < 1.E-10_rp) then
                NDeg = NDeg + 1
            end if
        end do

        write(CTmp, "('Block', I2.2, '.txt')") k
        open(20, file = "LancED_EE_" // trim(CTmp))
        do i = 1, NEV
            write(20, "(2I5, ES<rp_ow>.<rp_od>)") k, i, EgnVal(i)
        end do
        close(20)

        if(k == 1) then !please make sure the ground state is in the Mz=0 subspace.

        !calculate the static spin correlation
        call CorrFunc(DmHb, NDeg, EgnVal, EgnVec)

        open(30, file = "LancED_ET_SpnET_R.txt")
        do i = 1, NSite
            write(30, "(I5, 3ES<rp_ow>.<rp_od>)") i, SpnET(i, 1), SpnET(i, 2), SpnET(i, 1) + SpnET(i, 2)
        end do
        close(30)

        !fourier transform to get the static spin correlation in momentum space
        call FourierTransform1D

        end if
    
        deallocate(BsWf)
        deallocate(EgnVal)
        deallocate(EgnVec)
    
    end do
    
    call PrintMsg("End Lanczos Exact Diagonalization")
    call TimingEnd(timst, timed)

end subroutine

!======================================================================
    
subroutine TimingBgn(timst)

    use MyPrec
    implicit none
    integer(8) timst

    call system_clock(timst)

end subroutine
    
!======================================================================

subroutine TimingEnd(timst, timed)

    use MyPrec
    implicit none
    integer(8) timst, timed
    real(rp) AvgTime

    call system_clock(timed)
    
    AvgTime = TimeIntrvl(timst, timed)

    call PresentMoment
    open(21, file = trim(FLog), access = "append")
    write(21, 10) DatTim, "LancED total time:", &
                  AvgTime, AvgTime / 3600
10  format(A<DatTimLen>, 4x, A18, F10.4, "s = ", F5.2, "h")
    close(21)

end subroutine 

!======================================================================
    
subroutine BlockBasis

    use MyPrec
    use MyHmlt
    implicit none
    integer i, k

    k = 0
    !for example Ns=8, Nup=4, it ranges from 00001111 to 11110000.
    do i = ISHFT(1, NeUp) - 1, ISHFT(1, NSite - NeUp) * (ISHFT(1, NeUp) - 1)
        if(popcnt(i) == NeUp) then
            k = k + 1
            BsWf(k) = i
        end if
    end do
    
    if(k .ne. DmHb) then
        write(ErrMsg, "('BlockBasis k < DmHb, k = ', I11, ', DmHb = ', I11)") k, DmHb
        call ErrOut(ErrMsg)
    end if

end subroutine
    
!==================================================================================================

subroutine FourierTransform1D

    use MyPrec
    use MyHmlt
    implicit none
    integer i, m
    real(rp) rk, spnk

    open(140, file = "LancED_ET_StotET_K.txt")
    do m = 0, NSite
        rk = (2._rp * rp_pi / dble(NSite)) * dble(m)
        spnk = 0._rp
        do i = 1, NSite
            spnk = spnk + (SpnET(i, 1) + SpnET(i, 2)) * cos(dble(i - 1) * rk)
        end do
        spnk = spnk / dble(NSite)
        write(140, "(2ES<rp_ow>.<rp_od>)") rk, spnk
    end do
    close(140)

end subroutine

!======================================================================  
