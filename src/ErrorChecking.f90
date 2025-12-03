!======================================
! Module with functions designed to check for common
! programming errors such as Inf, NaN, etc.
!======================================
module ErrorChecking
use VarPrecision
!======================================
contains
!======================================
  function IsNan(a) result(nan)
    implicit none
    real(dp), intent(in) :: a
    logical :: nan


    if(a /= a) then
      nan = .true.
    else
      nan = .false.
    endif

  end function
!======================================
  function IsInf(a) result(inf)
    implicit none
    real(dp), intent(in) :: a
    logical :: inf
    real(dp), parameter :: posinf = HUGE(dp)
    real(dp), parameter :: neginf = -HUGE(dp)


    if(a > posinf) then
      inf = .true.
    else if(a < neginf) then
      inf = .true.
    else
      inf = .false.
    endif

  end function
!======================================
end module
!======================================
