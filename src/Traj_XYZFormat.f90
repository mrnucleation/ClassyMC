!====================================================================
module Traj_XYZ
  use VarPrecision
  use TrajectoryTemplate, only: trajectory

  type, public, extends(trajectory) :: trajXYZ
    contains
       procedure, pass :: WriteFrame => TrajXYZ_WriteFrame
  end type
!====================================================================
  contains
!====================================================================
  subroutine TrajXYZ_WriteFrame(self) 
    use BoxData, only: BoxArray
    implicit none
    class(trajXYZ), intent(in) :: self
    integer :: iAtom, jDim, boxNum
    integer :: nDim

    boxNum = self%boxNum
    nDim = BoxArray(boxNum)%box%nDimension

    write(self%fileUnit, *) BoxArray(boxNum)%box%nAtoms
    write(self%fileUnit, *) 
    do iAtom = 1, BoxArray(boxNum)%box%nAtoms
      write(self%fileUnit, *) "Ar", (BoxArray(boxNum)%box%atoms(jDim, iAtom), jDim=1,nDim)
    enddo


  end subroutine
!====================================================================
end module
!====================================================================
