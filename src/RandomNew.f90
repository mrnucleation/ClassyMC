!=======================================================
      module RandomGen
      use VarPrecision
      integer :: initseed
!=======================================================
      contains
!=======================================================
      real(dp) function grnd()
        implicit none
        real(dp) :: r

        call RANDOM_NUMBER(r)
        grnd = r
       
      end function
!=======================================================
      subroutine sgrnd(seed)
      implicit none
      integer, intent(in) :: seed
      integer :: i,n
      integer, allocatable :: tempSeed(:)
      
      call RANDOM_SEED(size=n)      
      allocate(tempSeed(1:n))
      tempSeed = seed + 37 * (/ (i - 1, i = 1, n) /)

      call RANDOM_SEED(put=tempSeed)
     
      deallocate(tempSeed)
     
      end subroutine      
!=======================================================
      real(dp) function Gaussian() result(num)
      use Constants
      use VarPrecision
      implicit none
      real(dp) :: y1, w, x1, x2

      x1 = grnd()
      x2 = grnd()
      y1 = sqrt( -2.0d0 * log(x1) ) * cos(two_pi*x2)
      num = y1

      end function
!=========================================================
      end module
