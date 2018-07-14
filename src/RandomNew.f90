!=======================================================
  module RandomGen
      use VarPrecision
      integer :: initseed = -1
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
!========================================================            
      subroutine Generate_UnitSphere(x,y,z)
      use VarPrecision
      implicit none
      real(dp), intent(out) :: x,y,z
      real(dp) :: u_12_sq, u1, u2
      
      u_12_sq = 2E0
      do while(u_12_sq .ge. 1)
       u1 = 2E0 * grnd() - 1E0
       u2 = 2E0 * grnd() - 1E0
       u_12_sq = u1 * u1 + u2 * u2
      enddo
 
      x = 2E0 * u1 * sqrt(1E0 - u_12_sq)
      y = 2E0 * u2 * sqrt(1E0 - u_12_sq)
      z = (1E0 - 2E0 * u_12_sq)
      
      end subroutine
!=======================================================
  function ListRNG(list, norm) result(bin)
    use Constants
    use VarPrecision
    implicit none
    real(dp), intent(in) :: list(:)
    real(dp), intent(in), optional :: norm
    integer :: bin, nSel
    real(dp) :: ran_num, intSum

    if( present(norm) ) then
      ran_num = grnd()*norm
    else
      ran_num = grnd()
    endif

    nSel = 1 
    intSum = list(1)
    do while(intSum < ran_num)
      nSel = nSel + 1 
      intSum = intSum + list(nSel)
    enddo

    bin = nSel
  end function
!=========================================================
end module
