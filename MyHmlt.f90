
!======================================================================
!	Define global variables and physical observables
!======================================================================

module MyHmlt

    use MyPrec
    implicit none

    !dimension and configuration of your lattices
    !finite lattice size, the number of spin-up spins
    !dimension of Hilbert space, ground state degeneracy
    !the number of Mz block, the number of total energy state been calculated
    integer DmLt, CfLt
    integer NSite, NeUp, DmHb, NDeg, NHBlk, NEGTH
    
    !tight binding model parameter
    !ground state energy
    real(rp) HeisJ, HeisJp, Anstr
    real(rp) GGSE

    !store the hopping and interaction info
    !static spin correlation function
    real(rp), allocatable :: JzzMat(:, :)
    real(rp), allocatable :: JpmMat(:, :)
    real(rp), allocatable :: SpnET(:, :)
    real(rp), allocatable :: ALLENG(:)
    
    !sparse matrix storage format for Lanzos ED
    integer, parameter :: IHM = 65536
    integer IHT
    integer,  allocatable :: IHI(:)
    integer,  allocatable :: IHL(:)
    integer,  allocatable :: IHR(:)
    real(rp), allocatable :: RHE(:)
    
    !binary representation of basis
    integer, allocatable :: BsWf(:)
    integer, allocatable :: BlkNeUp(:)
    integer, allocatable :: BlkHDim(:)
    
    !Hamiltonian matrix for full ED
    integer IfFED
    real(rp), allocatable :: RHMat(:, :)

end module

!======================================================================

subroutine MyHmltInit

    use MyPrec
    use MyHmlt
    implicit none
    integer err
    integer, external :: Binomial
    
    open(10, file = "param.in")
    read(10, *) DmLt
    read(10, *) CfLt
    read(10, *) IfFED
    read(10, *) NSite
    read(10, *) HeisJ
    read(10, *) HeisJp
    read(10, *) Anstr
    close(10)

    NHBlk = NSite + 1
    
    allocate(JzzMat(NSite, NSite), &
             JpmMat(NSite, NSite), &
             SpnET(NSite, 2), &
             BlkNeUp(NHBlk), &
             BlkHDim(NHBlk), stat = err)
    if(err .ne. 0) then
        write(*, *) "allocating error in MyHmltInit1!"
        stop
    end if

    JzzMat = 0._rp
    JpmMat = 0._rp
    SpnET  = 0._rp

end subroutine

!======================================================================

subroutine MyHmltFinl

    use MyHmlt
    implicit none
    integer err
    
    deallocate(JzzMat, JpmMat, SpnET, stat = err)
    if(err .ne. 0) then
        write(*, *) "deallocating error in MyHmltFinl !"
        stop
    end if

end subroutine

!======================================================================
    
subroutine MySMSFInit(ndim)

    use MyHmlt
    implicit none
    integer ndim, err
    
    allocate(IHL(ndim), &
             IHR(ndim), stat = err)
    if(err .ne. 0) then
        write(ErrMsg, "('allocating error in MySMSFInit')")
        call ErrOut(ErrMsg)
    end if
    
end subroutine
    
!======================================================================   
 
subroutine MySMSFFinl

    use MyHmlt
    implicit none

    if(allocated(IHL)) deallocate(IHL)
    if(allocated(IHR)) deallocate(IHR)
    if(allocated(IHI)) deallocate(IHI)
    if(allocated(RHE)) deallocate(RHE)
    
end subroutine
    
!======================================================================

function Binomial(n, m)
!combination C_n^m

    use MyPrec
    implicit none
    integer Binomial
    integer n, m, l, i
    real(16) a

    if (n < 0 .or. m < 0 .or. m > n) then
        write(WrnMsg, "('Binomial(n, m) n = ', I11, ', m = ', I11)") n, m
        call WrnOut(WrnMsg)
        return
    end if

    l = m
    if (l > n / 2) l = n - l

    a = 1._16
    do i = 1, l
        a = a * (n - l + i) / (l + 1 - i)
    end do
    Binomial = int(a + 0.5_16)

end function

!======================================================================
