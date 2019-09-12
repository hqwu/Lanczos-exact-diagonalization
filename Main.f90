!#####################################################
! PROGRAM: Lanczos Exact Diagonalization
! TYPE   : program
! PURPOSE: Lanczos exact diagonalization study spin
!          model. total Sz-conserved symmetry is used.
! I/O    :
! VERSION: 2019-09-09
! AUTHOR : Han-Qing Wu, Sun Yat-Sen Univeristy, Guangzhou
!          wuhanq3@mail.sysu.edu.cn
!#####################################################

!for linux, please include the following files.
!for visual studio with windows, do not include them when you create a project to compile.

!include "MyPrec.f90"
!include "QuickSort.f90"
!include "MyHmlt.f90"
!include "RLancED.f90"
!include "RFullED.f90"
!include "RS.f90"
!include "HmltElmt.f90"
!include "Lanczos.f90"
!include "CorrFunc.f90"

program Main

    use MyHmlt
    implicit none

    call MyPrecInit
    call MyRndmInit
    call FileIOInit
    call MyHmltInit
    call ModelInit
    call BlockList
    if(IfFED == 1) then
    call RFullED
    else
    call RLancED
    end if
    call MyHmltFinl

end program

!#####################################################

subroutine MyRndmInit

    implicit none
    integer(8) nclock
    integer i, nsd
    integer, allocatable :: oldseed(:), newseed(:)

    call random_seed()
    call random_seed(size=nsd)

    allocate(oldseed(nsd))
    allocate(newseed(nsd))

    call random_seed(get = oldseed(1:nsd))
    call system_clock(nclock)

    do i = 1, nsd
        newseed(i) = (oldseed(i) + ibits(nclock, 3*i, 3)) &
                     * ibits(nclock, 10 + 3*i, 3)
    end do
    
    call random_seed(put = newseed(1 : nsd))

    deallocate(oldseed)
    deallocate(newseed)

end subroutine

!=====================================================
    
subroutine FileIOInit

    use MyPrec
    implicit none
    character(64) CTmp

    write(CTmp, "('1D_Spin_Model.txt')")
    FLog = "Info_log_" // trim(CTmp)
    FWrn = "Info_wrn_" // trim(CTmp)
    FErr = "Info_err_" // trim(CTmp)

end subroutine

!=====================================================

subroutine ModelInit

    use MyPrec
    use MyHmlt
    implicit none
    integer i, j

    JzzMat = 0._rp
    JpmMat = 0._rp

    if(DmLt == 1 .and. CfLt == 1) then

    !nearest-neighbor exchange interaction
    do i = 1, NSite - 1
        JzzMat(i, i + 1) = HeisJ
        JpmMat(i, i + 1) = HeisJ * Anstr
    end do
    JzzMat(1, NSite) = HeisJ
    JpmMat(1, NSite) = HeisJ * Anstr
    
    !second nearest-neighbor exchange interaction
    do i = 1, NSite - 2
        JzzMat(i, i + 2) = HeisJp
        JpmMat(i, i + 2) = HeisJp * Anstr
    end do
    JzzMat(1, NSite - 1) = HeisJp
    JpmMat(1, NSite - 1) = HeisJp * Anstr
    JzzMat(2, NSite) = HeisJp
    JpmMat(2, NSite) = HeisJp * Anstr
    
    else
        write(*, *) "Fail to run the code !"
        stop
    end if
    
end subroutine

!=====================================================

subroutine BlockList

    use MyPrec
    use MyHmlt
    implicit none
    integer i, k
    integer, external :: Binomial
    
    k = 0
    do i = 0, NSite
        k = k + 1
        BlkHDim(k) = Binomial(NSite, i)
        BlkNeUp(k) = i
    end do
    
    call QuickSortII(-1, 1, k, BlkHDim, BlkNeUp)
    
end subroutine

!=====================================================
