!==================================================
module Box_Utility
use VarPrecision
use Template_SimBox, only: SimBox

!==================================================
contains
!==================================================
  subroutine FindAtom(box, globindx, atmIndx)
    implicit none
    class(SimBox), intent(in) :: box
    integer, intent(in) :: globIndx
    integer, intent(out) :: atmIndx


  end subroutine
!==================================================
end module
!==================================================

