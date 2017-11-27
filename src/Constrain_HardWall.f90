!=============================================================
! The purpose of this constraint is to provide the user a way to constrain all or parts of the system by a flat wall.
!
module Constrain_HardWall
  use ConstraintTemplate, only: constraint
  use VarPrecision

  type, public, extends(constraint) :: hardwall
    integer, allocatable :: atmTypes(:)
    logical :: wallAxis(1:3) = [.false., .false., .false.]
    real(dp) :: xhi, xlo
    real(dp) :: yhi, ylo
    real(dp) :: zhi, zlo
    contains
      procedure, pass :: CheckInitialConstraint => HardWall_CheckInitialConstraint
      procedure, pass :: ShiftCheck => HardWall_ShiftCheck
      procedure, pass :: SwapInCheck => HardWall_SwapInCheck
      procedure, pass :: SwapOutCheck => HardWall_SwapOutCheck
  end type

!=============================================================
  contains
!=============================================================
  subroutine HardWall_CheckInitialConstraint(self)
    implicit none
    class(hardwall), intent(in) :: self
  end subroutine
!=============================================================
  subroutine HardWall_ShiftCheck(self, trialBox, disp, accept)
    use CoordinateTypes, only: Displacement
    use Template_SimBox, only: SimBox
    implicit none
    class(hardwall), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    integer :: iDisp

    accept = .true.
    do iDisp = 1, size(disp)
      if(any(self%atmTypes == disp(iDisp)%atmindx) ) then
         !Check to see if the x-axis wall condition has been violated
        if( self%wallAxis(1) ) then
          if( (disp(iDisp)%x_new > self%xhi) .or.  (disp(iDisp)%x_new > self%xlo)) then
            accept = .false.            
            return
          endif 
        endif
 
         !Check to see if the y-axis wall condition has been violated
        if( self%wallAxis(2) ) then
          if( (disp(iDisp)%y_new > self%yhi) .or.  (disp(iDisp)%y_new > self%ylo)) then
            accept = .false.            
            return
          endif 
        endif

         !Check to see if the z-axis wall condition has been violated
        if( self%wallAxis(3) ) then
          if( (disp(iDisp)%z_new > self%zhi) .or.  (disp(iDisp)%z_new > self%zlo)) then
            accept = .false.            
            return
          endif 
        endif
      endif       
    enddo




  end subroutine
!=============================================================
  subroutine HardWall_SwapInCheck(self)
    implicit none
    class(hardwall), intent(in) :: self
  end subroutine
!=============================================================
  subroutine HardWall_SwapOutCheck(self)
    implicit none
    class(hardwall), intent(in) :: self
  end subroutine
!=============================================================
end module
!=============================================================
