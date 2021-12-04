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
  function IsSorted(list) result(sorted)
    !Checks if a list is already sorted from lowest to highest.
    implicit none
    integer, intent(in) :: list(:)
    logical :: sorted
    integer :: i
    integer :: upper, lower

    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)
    sorted = .true.
    if(lower == upper) then
      return
    endif
    do i = lower+1, upper
      if(list(i-1) > list(i)) then
        sorted = .false.
        return
      endif
    enddo
  end function
!======================================
  function SimpleSearch(val, list) result(outIndx)
    implicit none
    integer, intent(in) :: val
    integer, intent(in) :: list(:)
    integer :: i, outIndx
    integer :: upper, lower, listSize

    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)


    listSize = size(list)
    if(listSize < 1) then
      outIndx = 0
      return
    endif

!    write(*,*) curIndx, list(curIndx)
    outindx = -1
    do i = lower, upper
      if( list(i) == val ) then
        outIndx = i
        exit
      endif
    enddo
 

  end function
!======================================
  recursive function BinarySearch(val, list) result(outIndx)
    implicit none
    integer, intent(in) :: val
    integer, intent(in) :: list(:)
    integer :: outIndx
    integer :: upper, lower, curIndx
    integer :: listSize, loop


    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)


    listSize = size(list)
    if(listSize < 1) then
      outIndx = 0
      return
    endif


    curIndx = nint(0.5E0_dp*(lower + upper))
    if(list(curIndx) == val) then
      outIndx = curIndx
    elseif(list(curIndx) < val) then
      outIndx = BinarySearch(val, list(curIndx+1:upper))
      if(outIndx /= 0) then
        outIndx = outIndx + curIndx 
      endif
    else
      outIndx = BinarySearch(val, list(lower:curIndx-1))
      if(outIndx /= 0) then
        outIndx = outIndx + lower - 1
      endif
    endif

  end function
!======================================
  function BinarySearch_Old(val, list) result(outIndx)
    implicit none
    integer, intent(in) :: val
    integer, intent(in) :: list(:)
    integer :: outIndx
    integer :: upper, lower, curIndx
    integer :: listSize, loop


    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)

    curIndx = 1

    listSize = size(list)
    if(listSize < 1) then
      write(0,*) "Critical Error! A list size of 0 has been passed to the sort function!"
      error stop 
    endif


    loop = 0    
    do while( list(curIndx) /= val )
!      write(*,*) "-------------"
!      write(*,*) list
!      write(*,*)
!      write(2,*) "sort",val, curIndx, list(curIndx), lower, upper
      if(loop > listSize) then
        curIndx = 0
        exit
      endif


      if(upper-lower > 1) then
          curIndx = nint(0.5E0_dp*(lower + upper))
          if(list(curIndx) < val) then
            lower = curIndx
          elseif(list(curIndx) > val) then
            upper = curIndx
          endif
      else
          if(list(lower) == val) then
            curIndx = lower
          elseif(list(upper) == val) then
            curIndx = upper
          else
            curIndx = 0
          endif
          exit
      endif
!      write(*,*) curIndx, list(curIndx), lower, upper


!      write(*,*) "target=", val
      if( lower >= upper ) then
        curIndx = 0
        exit
      endif
      loop = loop + 1
    enddo

!    write(*,*) "Finish"
    if(curIndx /= 0) then
      if(list(curIndx) /= val) then
        curIndx = 0
      endif
    endif
    outIndx = curIndx

  end function
!======================================
recursive subroutine QSort(list)
  implicit none
  integer, intent(inout) ::  list(:)
  real(dp) :: x, t
  integer :: lower, upper
  integer :: i, j


  lower = LBound(list, dim=1)
  upper = UBound(list, dim=1)

  if(size(list(lower:upper)) < 2) then
    return
  endif
  x = list( (lower+upper) / 2 )
  i = lower
  j = upper
  do
    do while (list(i) < x)
      i=i+1
    end do
    do while (x < list(j))
      j=j-1
    end do
    if (i >= j) exit
!    t = list(i);  list(i) = list(j);  list(j) = t
    call Swap(list(i), list(j))
    i=i+1
    j=j-1
  end do
  if (lower < i-1) call QSort(list(lower:i-1))
  if (j+1 < upper)  call QSort(list(j+1:upper))
end subroutine QSort
!======================================
  recursive subroutine QSortOld(list)
    implicit none
    integer, intent(inout) :: list(:)
    integer :: piv

    if( size(list) > 1 ) then
      piv = QSort_Partition(list)
      call QSortOld( list(:piv-1) )
      call QSortOld( list(piv+1:) )
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
  subroutine BubbleSort(list)
    implicit none
    integer, intent(inout) ::  list(:)
    integer :: lower, upper
    integer :: i, j, n
    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)
    n = size(list) 
     
    do i =lower, upper-1 
      do j= lower, n-i-1 
        if(list(j) > list(j+1)) then
          call Swap(list(j), list(j+1))
        endif
      enddo
    enddo

  end subroutine

!======================================
  recursive logical function csr(a, left, right,n) result(swapped)
    implicit none
    integer, intent(in) :: left, right,n
    integer, intent(inout) :: a(n)
    integer :: lo, hi, mid
    integer :: temp
    logical :: lefthalf,righthalf
!
    swapped = .FALSE.
    if (right <= left) return
    lo = left   !Store the upper and lower bounds of list for
    hi = right  !Recursion later
!
    do while (lo < hi)
!   Swap the pair of elements if hi < lo
       if (a(hi) < a(lo)) then
          swapped = .TRUE.
          temp = a(lo)
          a(lo) = a(hi)
          a(hi) = temp
       endif
       lo = lo + 1
       hi = hi - 1
    end do
!   Special case if array is an odd size (not even)
    if (lo == hi)then
       if(a(hi+1) .lt. a(lo))then
           swapped = .TRUE.
           temp = a(hi+1)
           a(hi+1) = a(lo)
           a(lo) = temp
       endif
    endif
    mid = (left + right) / 2 ! Bisection point
    lefthalf = csr(a, left, mid,n)
    righthalf = csr(a, mid + 1, right,n)
    swapped = swapped .or. lefthalf .or. righthalf
  end function csr
! 
  subroutine circle_sort(list)
    use iso_c_binding, only: c_ptr, c_loc
    implicit none
    integer, target,intent(inout) :: list(:)
    integer :: n, lower, upper
    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)
    n = size(list(lower:upper)) 
 
    do while ( csr(list, 1, n,n))
! This is the canonical algorithm. However, if you want to
! speed it up, count the iterations and when you have approached
! 0.5*ln(n) iterations, perform a binary insertion sort then exit the loop.
    end do
  end subroutine circle_sort
 
!======================================
end module
!======================================
