!======================================
module SearchSort
use VarPrecision
!======================================
contains
!======================================
  function SimpleSearch(val, list) result(outIndx)
    implicit none
    integer, intent(in) :: val
    integer, intent(in) :: list(:)
    integer :: outIndx
    integer :: upper, lower, curIndx

    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)

    curIndx = lower
    do while( list(curIndx) /= val )
      curIndx = curIndx + 1
      if(curIndx >= upper) then
        exit
      endif
    enddo
 
    outIndx = curIndx

  end function
!======================================
  function BinarySearch(val, list) result(outIndx)
    implicit none
    integer, intent(in) :: val
    integer, intent(in) :: list(:)
    integer :: outIndx
    integer :: upper, lower, curIndx

    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)

    curIndx = 1
!    write(*,*) curIndx, list(curIndx), lower, upper
    do while( list(curIndx) /= val )
      curIndx = nint(0.5E0_dp*(lower + upper))
      if(list(curIndx) < val) then
        lower = curIndx
      elseif(list(curIndx) > val) then
        upper = curIndx
      endif

!      write(*,*) "target=", val
!      write(*,*) curIndx, list(curIndx), lower, upper
      if( lower >= upper ) then
        exit
      endif
    enddo
 
    outIndx = curIndx

  end function
!======================================
  subroutine FloatOrder(startIndx, list)
    implicit none
    integer, intent(in) :: startIndx
    real(dp), intent(inout) :: list(:)


  end subroutine
!======================================
end module
!======================================
