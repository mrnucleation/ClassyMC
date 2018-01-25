!=============================================================================+
module FF_ThermoIntegration
  use Template_ForceField, only: ForceField

  use VarPrecision
  use ForceFieldData, only: EnergyCalculator
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(forcefield) :: thermointegration
!    real(dp) :: rCut, rCutSq
    real(dp) :: lambda = 0E0_dp
    integer :: ECalc1 = -1
    integer :: ECalc2 = -1
    real(dp) :: E1 = 0E0_dp
    real(dp) :: E2 = 0E0_dp
    contains
      procedure, pass :: Constructor => ThermoInt_Constructor
      procedure, pass :: DetailedECalc => ThermoInt_DetailedECalc
!      procedure, pass :: DiffECalc
      procedure, pass :: ShiftECalc_Single => ThermoInt_ShiftECalc_Single
!      procedure, pass :: ShiftECalc_Multi
      procedure, pass :: NewECalc => ThermoInt_NewECalc
      procedure, pass :: OldECalc => ThermoInt_OldECalc
      procedure, pass :: VolECalc => ThermoInt_VolECalc
      procedure, pass :: ProcessIO => ThermoInt_ProcessIO
      procedure, pass :: GetCutOff => ThermoInt_GetCutOff
  end type

  contains
!=============================================================================+
  subroutine ThermoInt_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(thermointegration), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine ThermoInt_DetailedECalc(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    implicit none
    class(thermointegration), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    real(dp) :: ESub

    accept = .true.
    E_T = 0E0_dp

    call EnergyCalculator(self%ECalc1) % Method % DetailedECalc(curbox, ESub, accept)
    E_T = E_T + (1E0_dp - self%lambda) * ESub

    call EnergyCalculator(self%ECalc2) % Method % DetailedECalc(curbox, ESub, accept)
    E_T = E_T + self%lambda * ESub

    write(nout, *) "Thermo Integration Energy:", E_T
  end subroutine
!=============================================================================+
  subroutine ThermoInt_ShiftECalc_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(thermointegration), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    real(dp) :: ESub

    accept = .true.
    E_Diff = 0E0_dp

    call EnergyCalculator(self%ECalc1) % Method %ShiftECalc_Single(curbox, disp, ESub, accept)
    E_Diff = E_Diff + (1E0_dp - self%lambda) * ESub
    if(.not. accept) then
      return
    endif

    call EnergyCalculator(self%ECalc2) % Method %ShiftECalc_Single(curbox, disp, ESub, accept)
    if(.not. accept) then
      return
    endif
    E_Diff = E_Diff + self%lambda * ESub

  end subroutine
!=============================================================================+
  subroutine ThermoInt_NewECalc(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(thermointegration), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    real(dp) :: ESub

    accept = .true.
    E_Diff = 0E0_dp

    call EnergyCalculator(self%ECalc1) % Method % NewECalc(curbox, disp, tempList, tempNNei, ESub, accept)

    E_Diff = E_Diff + (1E0_dp - self%lambda) * ESub
    if(.not. accept) then
      return
    endif

    call EnergyCalculator(self%ECalc2) % Method % NewECalc(curbox, disp, tempList, tempNNei, ESub, accept)
    if(.not. accept) then
      return
    endif
    E_Diff = E_Diff + self%lambda * ESub

  end subroutine
!=============================================================================+
  subroutine ThermoInt_OldECalc(self, curbox, disp, E_Diff)
    implicit none
    class(thermointegration), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff

    real(dp) :: ESub

    E_Diff = 0E0_dp

    call EnergyCalculator(self%ECalc1) % Method %OldECalc(curbox, disp, ESub)
    E_Diff = E_Diff + (1E0_dp - self%lambda) * ESub

    call EnergyCalculator(self%ECalc2) % Method %OldECalc(curbox, disp, ESub)
    E_Diff = E_Diff + self%lambda * ESub


  end subroutine
!=============================================================================+
  subroutine ThermoInt_VolECalc(self, curbox, scalars, E_Diff)
    implicit none
    class(thermointegration), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(in) :: scalars(:)
    real(dp), intent(inOut) :: E_Diff

    real(dp) :: ESub
    call EnergyCalculator(self%ECalc1) % Method %VolECalc(curbox, scalars, ESub)
    E_Diff = E_Diff + (1E0_dp - self%lambda) * ESub


    call EnergyCalculator(self%ECalc2) % Method %VolECalc(curbox, scalars, ESub)
    E_Diff = E_Diff + self%lambda * ESub


  end subroutine
!=============================================================================+
  subroutine ThermoInt_ProcessIO(self, line)
    implicit none
    class(thermointegration), intent(inout) :: self
    character(len=*), intent(in) :: line

  end subroutine
!=============================================================================+
  function ThermoInt_GetCutOff(self) result(rCut)
    implicit none
    class(thermointegration), intent(inout) :: self
    real(dp) :: rCut

!    write(*,*) self%rCut
    rCut = self%rCut
  end function
!=============================================================================+
end module
!=============================================================================+
