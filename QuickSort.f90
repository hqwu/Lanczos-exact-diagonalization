    
!==================================================================================================
    
recursive subroutine QuickSortR(IncDec, BgnIdx, EndIdx, RKeyVec)

    use MyPrec
    implicit none
    integer BgnIdx, EndIdx, IncDec, l, r
    integer RKeyVal
    real(rp) RKeyVec(BgnIdx : EndIdx)

    if (BgnIdx >= EndIdx) return

    l = BgnIdx
    r = EndIdx + 1
    RKeyVal = RKeyVec(BgnIdx)

    do while(.true.)

        do while(IncDec == 1)
            l = l + 1
            if (l >= EndIdx .or. RKeyVec(l) > RKeyVal) exit
        end do
        
        do while(IncDec == -1)
            l = l + 1
            if (l >= EndIdx .or. RKeyVec(l) < RKeyVal) exit
        end do
        
        do while(IncDec == 1)
            r = r - 1
            if (r <= BgnIdx .or. RKeyVec(r) < RKeyVal) exit
        end do

        do while(IncDec == -1)
            r = r - 1
            if (r <= BgnIdx .or. RKeyVec(r) > RKeyVal) exit
        end do

        if (r <= l) exit

        if (l /= r) call RMySwap(RKeyVec(l),  RKeyVec(r) )

    end do

    if (BgnIdx /= r) call RMySwap(RKeyVec (BgnIdx), RKeyVec(r) )

    if (r - 1 > BgnIdx) call QuickSortR(IncDec, BgnIdx, r - 1,  RKeyVec(BgnIdx))
    if (r + 1 < EndIdx) call QuickSortR(IncDec, r + 1,  EndIdx, RKeyVec(r + 1) )

end subroutine
    
!==================================================================================================
    
recursive subroutine QuickSortII(IncDec, BgnIdx, EndIdx, IKeyVec, ISibVec1)

    use MyPrec
    implicit none
    integer BgnIdx, EndIdx, IncDec, l, r
    integer IKeyVal
    integer IKeyVec(BgnIdx : EndIdx)
    integer ISibVec1(BgnIdx : EndIdx)

    if (BgnIdx >= EndIdx) return

    l = BgnIdx
    r = EndIdx + 1
    IKeyVal = IKeyVec(BgnIdx)

    do while(.true.)

        do while(IncDec == 1)
            l = l + 1
            if (l >= EndIdx .or. IKeyVec(l) > IKeyVal) exit
        end do
        
        do while(IncDec == -1)
            l = l + 1
            if (l >= EndIdx .or. IKeyVec(l) < IKeyVal) exit
        end do
        
        do while(IncDec == 1)
            r = r - 1
            if (r <= BgnIdx .or. IKeyVec(r) < IKeyVal) exit
        end do

        do while(IncDec == -1)
            r = r - 1
            if (r <= BgnIdx .or. IKeyVec(r) > IKeyVal) exit
        end do

        if (r <= l) exit

        if (l /= r) call IMySwap(IKeyVec(l),  IKeyVec(r) )
        if (l /= r) call IMySwap(ISibVec1(l), ISibVec1(r))

    end do

    if (BgnIdx /= r) call IMySwap(IKeyVec (BgnIdx), IKeyVec(r) )
    if (BgnIdx /= r) call IMySwap(ISibVec1(BgnIdx), ISibVec1(r))

    if (r - 1 > BgnIdx) call QuickSortII(IncDec, BgnIdx, r - 1,  IKeyVec(BgnIdx), ISibVec1(BgnIdx))
    if (r + 1 < EndIdx) call QuickSortII(IncDec, r + 1,  EndIdx, IKeyVec(r + 1),  ISibVec1(r + 1) )

end subroutine
    
!==================================================================================================

recursive subroutine QuickSortIR(IncDec, BgnIdx, EndIdx, IKeyVec, RSibVec)

	use MyPrec
	implicit none
	integer BgnIdx, EndIdx, IncDec, l, r
	integer IKeyVal
	integer IKeyVec(BgnIdx : EndIdx)
	real(rp) RSibVec(BgnIdx : EndIdx)

	if (BgnIdx >= EndIdx) return

	l = BgnIdx
	r = EndIdx + 1

	IKeyVal = IKeyVec(BgnIdx)

	do while(.true.)

        do while(IncDec == 1)
            l = l + 1
            if (l >= EndIdx .or. IKeyVec(l) > IKeyVal) exit
        end do
        
        do while(IncDec == -1)
            l = l + 1
            if (l >= EndIdx .or. IKeyVec(l) < IKeyVal) exit
        end do
        
        do while(IncDec == 1)
            r = r - 1
            if (r <= BgnIdx .or. IKeyVec(r) < IKeyVal) exit
        end do

        do while(IncDec == -1)
            r = r - 1
            if (r <= BgnIdx .or. IKeyVec(r) > IKeyVal) exit
        end do

		if (r <= l) exit

		if (l /= r) call IMySwap(IKeyVec(l), IKeyVec(r))
		if (l /= r) call RMySwap(RSibVec(l), RSibVec(r))

	end do

	if (BgnIdx /= r) call IMySwap(IKeyVec(BgnIdx), IKeyVec(r))
	if (BgnIdx /= r) call RMySwap(RSibVec(BgnIdx), RSibVec(r))

	if (r - 1 > BgnIdx) call QuickSortIR(IncDec, BgnIdx, r - 1,  IKeyVec(BgnIdx), RSibVec(BgnIdx))
	if (r + 1 < EndIdx) call QuickSortIR(IncDec, r + 1,  EndIdx, IKeyVec(r + 1),  RSibVec(r + 1) )

end subroutine

!==================================================================================================