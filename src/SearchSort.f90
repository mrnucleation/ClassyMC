!======================================
module SearchSort
use VarPrecision
!======================================
contains
!======================================
  subroutine Swap(a, b)
    implicit none
    integer, intent(inout) :: a,b
    integer :: c

    c = a
    a = b
    b = c

  end subroutine
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
        curIndx = 0
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
        curIndx = 0
        exit
      endif
    enddo
 
    outIndx = curIndx

  end function
!======================================
  recursive subroutine QSort(list)
    implicit none
    integer, intent(inout) :: list(:)
    integer :: piv

    if( size(list) > 1 ) then
      piv = QSort_Partition(list)
      call QSort( list(:piv-1) )
      call QSort( list(piv+1:) )
    endif


  end subroutine
!======================================
  function QSort_Partition(list) result(indx)
    implicit none
    integer, intent(inout) :: list(:)
    integer :: lo, hi
    integer :: val, indx

    integer :: i, j

    lo = LBound(list, 1)
    hi = UBound(list, 1)
    val = list(1)
  
    i = 0
    j = size(list) + 1
    do
      j = j - 1
      do 
        if(list(j) <= val) then
          exit
        endif
        j = j - 1
      enddo
      i = i + 1
      do
        if (list(i) >= val) then
          exit
        endif
        i = i + 1
      enddo
      if( i < j ) then
        call swap( list(i), list(j) )
      elseif(i == j) then
        indx = i + 1
        return
      else
        indx = i
        return
      endif
    enddo

  end function
!======================================
end module
!======================================
