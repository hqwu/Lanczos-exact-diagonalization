
!======================================================================
!	Real/Integer precision and some constant variables
!======================================================================

module MyPrec

    implicit none

    !define some constants of real datatype
    integer, parameter :: rp = selected_real_kind(15,  307)
    integer rp_ow, rp_od
    real(rp) rp_pi, rp_prec, rp_tiny
    real(rp) rp_huge, rp_gold, rp_zero
    complex(rp) ImagUnit
    
    !define some parameters for analyzing computing time
    integer(rp) TimRat, TimMax
    real(rp) TimThh
    integer, parameter :: DatTimLen = 19
    character(8)  ccyymmdd
    character(10) hhmmss
    character(19) DatTim
    
    !define some files for output
    character(64)  FLog, FErr, FWrn
    character(128) ErrMsg, WrnMsg

    contains

    function expi(phi)
        implicit none
        real(rp) phi
        complex(rp) expi

        expi = cmplx(cos(phi), sin(phi), rp)
    
    end function
    
    function TimeIntrvl(timst, timed)
    
        implicit none
        integer(8) timst, timed, tims
        real(8) TimeIntrvl
        
        tims = timed - timst
        if(tims < 0) tims = tims + TimMax + 1
        TimeIntrvl = tims * 1._8 / TimRat
        
    end function
    
    subroutine PresentMoment
    
        implicit none
        
        call date_and_time(ccyymmdd, hhmmss)
        write(DatTim, 11) ccyymmdd(1 : 4), ccyymmdd(5 : 6), &
        ccyymmdd(7 : 8), hhmmss(1 : 2), hhmmss(3 : 4), hhmmss(5 : 6)
11      format(A4, "-", A2, "-", A2, " ", A2, ":", A2, ":", A2)
                          
    end subroutine
    
    subroutine PrintMsg(AMsg)

        implicit none
        character(*) AMsg

        call PresentMoment
        open(21, file = trim(Flog), access = "append")
        write(21, *) DatTim, "  ",  AMsg
        close(21)
    
    end subroutine
    
    subroutine ErrOut(AMsg)

        implicit none
        character(*) AMsg

        open(49, file = trim(FErr), access = "append")
        write(49, "(A)") "ERR: " // trim(AMsg)
        close(49)
    
    end subroutine
    
    subroutine WrnOut(AMsg)

        implicit none
        character(*) AMsg

        open(49, file = trim(FWrn), access = "append")
        write(49, "(A)") "WRN: " // trim(AMsg)
        close(49)

    end subroutine

end module
    
!======================================================================
    
subroutine MyPrecInit

    use MyPrec
    implicit none
    integer(rp) tims

    rp_pi    = acos(-1._rp)
    rp_tiny  = tiny(rp_pi)
    rp_huge  = huge(rp_pi)
    rp_gold  = (sqrt(5._rp) - 1._rp) / 2
    ImagUnit = cmplx(0._rp, 1._rp, rp)
    
    call MachineRealPrec
    
    if(rp == 4) then
        rp_ow = 16
        rp_od = 8
        rp_zero = 1.E-4_rp
    else if(rp == 8) then
        rp_ow = 25
        rp_od = 16
        rp_zero = 1.E-12_rp
    else if(rp == 16) then
        rp_ow = 44
        rp_od = 34
        rp_zero = 1.E-30_rp
    else
        write(*, *) "RealPreInit rp = ", rp
        stop
    end if
    
    call system_clock(tims, TimRat, TimMax)
    TimThh = TimMax * 1._8 / TimRat

end subroutine
    
!======================================================================
    
subroutine MachineRealPrec

    use MyPrec
    implicit none
	real(rp) x, y, d
	integer i

    d = 2._rp
    x = 1._rp
    y = x + 1._rp

    do i = 1, 100
    do while (y > x .and. 1._rp + x > 1._rp)
        y = x
        x = x / d
    end do
    x = y
    y = x + 1._rp
    d = sqrt(d)
    end do

    rp_prec = x

end subroutine
    
!======================================================================
