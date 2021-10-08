!===================================================================
module Template_IntraMiscIntra
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes

  type, public, extends(Intra_FF) :: MiscIntra_FF
    integer, private :: moltype = 0
    contains
      procedure, pass :: DetailedECalc
      procedure, pass :: SetMolType
      procedure, pass :: GetMolType
  end type

  contains
!===================================================================
  subroutine DetailedECalc(self, curbox, atompos, E_T, accept)
    implicit none
    class(MiscIntra_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    accept = .true.
    E_T = 0E0_dp

  end subroutine
!===================================================================
  subroutine SetMolType(self, moltype)
    implicit none
    class(MiscIntra_FF), intent(inout) :: self
    integer, intent(in) :: molType

    self%moltype = molType

  end subroutine
!===================================================================
  function GetMolType(self) result(molType)
    implicit none
    class(MiscIntra_FF), intent(in) :: self
    integer :: molType

    molType = self%moltype 

  end function
!===================================================================
end module
!===================================================================
