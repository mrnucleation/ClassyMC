!=============================================================================+
module FF_Hybrid
  use Template_ForceField, only: ForceField
  use VarPrecision
  use ForceFieldData, only: EnergyCalculator
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(forcefield) :: Pair_Hybrid
!    real(dp) :: rCut, rCutSq
    integer :: NFFields = 1
    integer, allocatable :: ECalcIndx(:)
    real(dp), allocatable :: EDiff(:)
    contains
!      procedure, pass :: Constructor => Hybrid_Constructor
      procedure, pass :: DetailedECalc => Hybrid_DetailedECalc
!      procedure, pass :: DiffECalc
!      procedure, pass :: ShiftECalc_Single => Hybrid_ShiftECalc_Single
!      procedure, pass :: ShiftECalc_Multi
!      procedure, pass :: NewECalc => Hybrid_NewECalc
!      procedure, pass :: OldECalc => Hybrid_OldECalc
!      procedure, pass :: VolECalc => Hybrid_VolECalc
      procedure, pass :: ProcessIO => Hybrid_ProcessIO
      procedure, pass :: GetCutOff => Hybrid_GetCutOff
      procedure, pass :: Update => Hybrid_Update
  end type

  contains
!=============================================================================+
!  subroutine Hybrid_Constructor(self)
!    use Common_MolInfo, only: nMolTypes
!    implicit none
!    class(Pair_Hybrid), intent(inout) :: self
!
!
!  end subroutine
!=============================================================================+
  subroutine Hybrid_DetailedECalc(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    integer :: iField

    real(dp) :: ESub

    accept = .true.
    E_T = 0E0_dp

    do iField = 1, self%NFFields 
      call EnergyCalculator(self%ECalcIndx(iField)) % Method % DetailedECalc(curbox, ESub, accept)
      E_T = E_T + ESub
    enddo


    write(nout, *) "Hybrid Forcefield Energy:", E_T
  end subroutine
!============================================================================
  subroutine Hybrid_DiffECalc(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(inout) :: disp(:)
    integer, intent(in), pointer :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    integer :: iField
    real(dp) :: E_Half
    real(dp) :: ESub

    accept = .true.
    E_Diff = 0E0_dp
    self%EDiff = 0E0_dp

    do iField = 1, self%NFFields 
        call EnergyCalculator(self%ECalcIndx(iField)) % Method %DiffECalc(curbox, disp, tempList, tempNNei, ESub, accept)
      if(.not. accept) then
        self%EDiff = 0E0_dp
        return
      endif
      self%EDiff(iField) = ESub
      E_Diff = E_Diff + ESub
    enddo

  end subroutine
!=============================================================================+
  subroutine Hybrid_ProcessIO(self, line)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

    integer :: intVal
    integer :: lineStat, iField
    character(len=30) :: command

    call GetXCommand(line, command, 1, lineStat)
    select case(trim(adjustl(command)))
      case("nfields")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) intVal
        self%NFFields = intVal
        if(allocated(self%ECalcIndx)) then
          deallocate(self%ECalcIndx)
          deallocate(self%EDiff)
        endif
        allocate(self%ECalcIndx(1:intVal))
        allocate(self%EDiff(1:intVal))
        self%ECalcIndx = -1
        self%EDiff = 0E0_dp

      case("forcefields")
        do iField = 1, self%NFFields 
          call GetXCommand(line, command, 1+iField, lineStat)
          read(command, *) intVal
          self%ECalcIndx(iField) = intVal
        enddo

      case default
        lineStat = -1
    end select


  end subroutine
!=============================================================================+
  function Hybrid_GetCutOff(self) result(rCut)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    integer :: iField
    real(dp) :: rCut, rCutCur

    rCut = 0E0_dp

    do iField = 1, self%NFFields
      rCutCur = EnergyCalculator(self%ECalcIndx(iField)) % Method % GetCutOff()
      rCut = max(rCut, rCutCur)
    enddo

    self%rCut = rCut
    self%rCutSq = rCut*rCut

  end function
!=============================================================================+
  subroutine Hybrid_Update(self)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self

    self%EDiff = 0E0_dp
  end subroutine
!=============================================================================+
end module
!=============================================================================+
