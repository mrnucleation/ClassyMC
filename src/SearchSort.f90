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
    integer :: upper, lower, curIndx, listSize

    lower = LBound(list, dim=1)
    upper = UBound(list, dim=1)

    curIndx = lower

    listSize = size(list)
    if(listSize < 1) then
      outIndx = 0
      return
!      write(0,*) "Critical Error! A list size of 0 has been passed to the sort function!"
!      error stop 
    endif

!    write(*,*) curIndx, list(curIndx)
    do 
      if( list(curIndx) == val ) then
        exit
      endif

      curIndx = curIndx + 1   
      if(curIndx > upper) then
        curIndx = 0
        exit
      endif
    enddo
 
    outIndx = curIndx

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
  if(size(list) < 1) then
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
    t = list(i);  list(i) = list(j);  list(j) = t
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
end module
!======================================
