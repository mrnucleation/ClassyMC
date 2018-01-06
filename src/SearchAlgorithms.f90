!======================================
module Search
use VarPrecision
!======================================
contains
!======================================
  function BinarySearch(val, list) result(outIndx)
    implicit none
    integer, intent(in) :: val
    real(dp), intent(in) :: list(:)
    integer :: outIndx
    integer :: upper, lower, curIndx

    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)

    curIndx = 1
    do while( list(curIndx) /= val )
      curIndx = nint(0.5E0_dp*(lower + upper))
      if(list(curIndx) < val) then
        lower = curIndx
      elseif(list(curIndx) > val) then
        upper = curIndx
      endif

      if( lower >= upper ) then
        exit
      endif
    enddo
 
    outIndx = curIndx

  end function
!======================================
end module
!======================================
