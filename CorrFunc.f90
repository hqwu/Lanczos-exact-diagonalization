
!======================================================================
    
subroutine CorrFunc(ND, NEV, EE, EV)

    use MyPrec
    use MyHmlt
    implicit none
    
    integer i, k, ND, NEV, err
    real(rp) EE(NEV), EV(ND, NEV)
    real(rp), allocatable :: YV(:)
    
    call PrintMsg("Start calculating static spin correlations")
    
    allocate(YV(ND), stat = err)
    if(err .ne. 0) then
        write(ErrMsg, "('allocating error in CorrFunc')")
        call ErrOut(ErrMsg)
    end if

    YV = 0._rp
    
    do k = 1, NEV

        do i = 1, NSite
            call Operate_SizSjz(EV(:, k), YV, ND, 1, i)
            SpnET(i, 1) = SpnET(i, 1) + Dot_Product(EV(:, k), YV)
        end do

        do i = 1, NSite
            call Operate_SipSjm(EV(:, k), YV, ND, 1, i)
            SpnET(i, 2) = SpnET(i, 2) + Dot_Product(EV(:, k), YV)
        end do

    end do
    
    deallocate(YV, stat = err)
    if(err .ne. 0) then
        write(ErrMsg, "('deallocating error in CorrFunc')")
        call ErrOut(ErrMsg)
    end if
    
    call PrintMsg("End  calculating static spin correlations")
    
end subroutine

!======================================================================

subroutine Operate_SizSjz(xv, yv, nd, i, j)

    use MyPrec
    use MyHmlt
    implicit none
    integer Cfg
    integer m, nd, i, j
    real(rp) xv(nd), yv(nd)

    yv = 0._rp
    do m = 1, nd
        Cfg = BsWf(m)
        yv(m) = yv(m) + (dble(ibits(Cfg, i-1, 1)) - 0.5_rp) * &
                      (dble(ibits(Cfg, j-1, 1)) - 0.5_rp) * xv(m)
    end do
    
end subroutine 
    
!======================================================================

subroutine Operate_SipSjm(xv, yv, nd, i, j)

    use MyPrec
    use MyHmlt
    implicit none
    integer xL, xR, Lct
    integer i, j, m, nd
    integer Cfg
    real(rp) xv(nd), yv(nd)

    yv = 0._rp
    do m = 1, nd
        Cfg = BsWf(m)
        if(i == j) then
            yv(m) = yv(m) + xv(m) * 0.5_rp
        else 
        if(ibits(Cfg, i-1, 1) + ibits(Cfg, j-1, 1) == 1) then
            call BasisIndex(ibchng(ibchng(Cfg, i-1), j-1), DmHb, BsWf, Lct)
            yv(Lct) = yv(Lct) + xv(m) * 0.5_rp
        end if
        end if
    end do
    
end subroutine

!======================================================================
