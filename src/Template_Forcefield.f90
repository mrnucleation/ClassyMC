module Template_ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public :: forcefield
    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor 
      procedure, pass :: DetailedECalc 
      procedure, pass :: DiffECalc
      procedure, pass :: ShiftECalc_Single
      procedure, pass :: ShiftECalc_Multi
      procedure, pass :: NewECalc
      procedure, pass :: OldECalc
      procedure, pass :: ExchangeECalc
      procedure, pass :: ProcessIO
      procedure, pass :: GetCutOff
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(forcefield), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, E_T)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T

  end subroutine
!============================================================================
  subroutine DiffECalc(self, curbox, disp, tempList, tempNNei, E_Diff)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    real(dp) :: E_Half

    if(disp(1)%newAtom .and. disp(1)%oldAtom) then
      if(disp(1)%oldAtmIndx == disp(1)%atmIndx) then
        call self % ShiftECalc_Single(curbox, disp, E_Diff)
      else
        E_Diff = 0E0_dp
        call self % NewECalc(curbox, disp, tempList, tempNNei, E_Half)
        E_Diff = E_Diff + E_Half
        call self % OldECalc(curbox, disp, E_Half)
        E_Diff = E_Diff + E_Half
      endif

      return
    endif

    if(disp(1)%newAtom) then
      call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff)
      return
    endif

    if(disp(1)%oldAtom) then
      call self % OldECalc(curbox, disp, E_Diff)
      return
    endif


  end subroutine
!=============================================================================+
  subroutine ShiftECalc_Single(self, curbox, disp, E_Diff)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff

  end subroutine
!=============================================================================+
  subroutine ShiftECalc_Multi(self, curbox, disp, E_Diff)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inout) :: E_Diff

  end subroutine
!=============================================================================+
  subroutine NewECalc(self, curbox, disp, tempList, tempNNei, E_Diff)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff

  end subroutine
!=============================================================================+
  subroutine OldECalc(self, curbox, disp, E_Diff)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
  end subroutine
!=============================================================================+
  subroutine ExchangeECalc(self, curbox, atmIndx, newType, E_Diff)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_Diff
    integer, intent(in) :: atmIndx, newType
  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    implicit none
    class(forcefield), intent(inout) :: self
    character(len=*), intent(in) :: line

  end subroutine
!=============================================================================+
  function GetCutOff(self) result(rCut)
    implicit none
    class(forcefield), intent(inout) :: self
    real(dp) :: rCut

!    write(*,*) self%rCut
    rCut = self%rCut
  end function
!=============================================================================+
end module
!=============================================================================+
