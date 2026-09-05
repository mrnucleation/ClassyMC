!=============================================================================+
module FF_Hybrid
  use Template_ForceField, only: ForceField
  use VarPrecision
  use ForceFieldData, only: EnergyCalculator
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(forcefield) :: Pair_Hybrid
!    real(dp) :: rCut, rCutSq
    integer :: NFFields = 0
    integer, allocatable :: ECalcIndx(:)
    real(dp), allocatable :: EDiff(:)
    contains
      procedure, pass :: Constructor => Hybrid_Constructor
      procedure, pass :: DetailedECalc => Hybrid_DetailedECalc
      procedure, pass :: DiffECalc => Hybrid_DiffECalc
      procedure, pass :: HasComputeForces => Hybrid_HasComputeForces
      procedure, pass :: ComputeForces => Hybrid_ComputeForces
      procedure, pass :: ProcessIO => Hybrid_ProcessIO
!      procedure, pass :: ShiftECalc_Single => Hybrid_ShiftECalc_Single
!      procedure, pass :: NewECalc => Hybrid_NewECalc
!      procedure, pass :: OldECalc => Hybrid_OldECalc
!      procedure, pass :: OrthoVolECalc => Hybrid_OrthoVolECalc
      procedure, pass :: GetCutOff => Hybrid_GetCutOff
      procedure, pass :: Update => Hybrid_Update
  end type

  contains
!=============================================================================+
  subroutine Hybrid_Constructor(self)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self

    self%NFFields = 0
    self%rCut = 0E0_dp
    self%rCutSq = 0E0_dp

  end subroutine
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
      write(nout,*) "SubEnergy", iField,":", ESub
      E_T = E_T + ESub
    enddo


    write(nout, *) "Hybrid Forcefield Energy:", E_T
  end subroutine
!============================================================================
  subroutine Hybrid_DiffECalc(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(inout), target :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    integer :: iField
    real(dp) :: ESub

    accept = .true.
    E_Diff = 0E0_dp
    if(allocated(self%EDiff)) then
      self%EDiff = 0E0_dp
    endif

    ! Check if this is a volume change move
    do iField = 1, self%NFFields 
      call EnergyCalculator(self%ECalcIndx(iField)) % Method % DiffECalc(curbox, disp, tempList, tempNNei, ESub, accept)
      if(.not. accept) then
        if(allocated(self%EDiff)) self%EDiff = 0E0_dp
        return
      endif
      if(allocated(self%EDiff)) self%EDiff(iField) = ESub
      E_Diff = E_Diff + ESub
    enddo

    ! For volume changes, sub-forcefields return E_new (not E_new - E_old)
    ! The hybrid handles the E_Inter subtraction here
    select type(disp)
      class is(OrthoVolChange)
        E_Diff = E_Diff - curbox%E_Inter
        curbox%dETable = curbox%dETable - curbox%ETable
      class is(NewState_IsoMol)
        E_Diff = E_Diff - curbox%E_Inter
        curbox%dETable = curbox%dETable - curbox%ETable
    end select



  end subroutine
!=============================================================================+
  function Hybrid_HasComputeForces(self) result(hasForce)
    implicit none
    class(Pair_Hybrid), intent(in) :: self
    logical :: hasForce

    hasForce = .true.
  end function
!=============================================================================+
  subroutine Hybrid_ComputeForces(self, curbox, forces, coords)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(out) :: forces(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer :: iField
    real(dp) :: fsub(size(forces,1), size(forces,2))

    forces = 0E0_dp
    do iField = 1, self%NFFields
      if(EnergyCalculator(self%ECalcIndx(iField)) % Method % HasComputeForces()) then
        call EnergyCalculator(self%ECalcIndx(iField)) % Method % ComputeForces(curbox, fsub, coords)
      else
        call EnergyCalculator(self%ECalcIndx(iField)) % Method % ComputeForcesFiniteDiff(curbox, fsub, coords)
      endif
      forces = forces + fsub
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
          ! Mark sub-forcefields so they don't subtract E_Inter themselves
          EnergyCalculator(intVal) % Method % isSubPair = .true.
        enddo

      case default
        lineStat = -1
    end select


  end subroutine
!=============================================================================+
! The following subroutines are commented out because they are not part of the
! base forcefield template interface. They can be enabled once the base template
! is updated to include these methods.
!=============================================================================+
!  subroutine Hybrid_ShiftECalc_Single(self, curbox, disp, E_Diff, accept, tempList, tempNNei)
!    implicit none
!    class(Pair_Hybrid), intent(inout) :: self
!    class(SimBox), intent(inout) :: curbox
!    type(Displacement), intent(inout) :: disp(:)
!    integer, intent(in), target, optional :: tempList(:,:), tempNNei(:)
!    real(dp), intent(inOut) :: E_Diff
!    logical, intent(out) :: accept
!    integer :: iField
!    real(dp) :: ESub
!
!    accept = .true.
!    E_Diff = 0E0_dp
!    self%EDiff = 0E0_dp
!
!    do iField = 1, self%NFFields
!      call EnergyCalculator(self%ECalcIndx(iField)) % Method % ShiftECalc_Single(curbox, disp, ESub, accept, tempList, tempNNei)
!      if(.not. accept) then
!        self%EDiff = 0E0_dp
!        return
!      endif
!      self%EDiff(iField) = ESub
!      E_Diff = E_Diff + ESub
!    enddo
!
!  end subroutine
!=============================================================================+
!  subroutine Hybrid_NewECalc(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
!    implicit none
!    class(Pair_Hybrid), intent(inout) :: self
!    class(SimBox), intent(inout) :: curbox
!    type(Addition), intent(inout) :: disp(:)
!    integer, intent(in) :: tempList(:,:), tempNNei(:)
!    real(dp), intent(inOut) :: E_Diff
!    logical, intent(out) :: accept
!    integer :: iField
!    real(dp) :: ESub
!
!    accept = .true.
!    E_Diff = 0E0_dp
!    self%EDiff = 0E0_dp
!
!    do iField = 1, self%NFFields
!      call EnergyCalculator(self%ECalcIndx(iField)) % Method % NewECalc(curbox, disp, tempList, tempNNei, ESub, accept)
!      if(.not. accept) then
!        self%EDiff = 0E0_dp
!        return
!      endif
!      self%EDiff(iField) = ESub
!      E_Diff = E_Diff + ESub
!    enddo
!
!  end subroutine
!=============================================================================+
!  subroutine Hybrid_OldECalc(self, curbox, disp, E_Diff)
!    implicit none
!    class(Pair_Hybrid), intent(inout) :: self
!    class(SimBox), intent(inout) :: curbox
!    type(Deletion), intent(inout) :: disp(:)
!    real(dp), intent(inOut) :: E_Diff
!    integer :: iField
!    real(dp) :: ESub
!
!    E_Diff = 0E0_dp
!    self%EDiff = 0E0_dp
!
!    do iField = 1, self%NFFields
!      call EnergyCalculator(self%ECalcIndx(iField)) % Method % OldECalc(curbox, disp, ESub)
!      self%EDiff(iField) = ESub
!      E_Diff = E_Diff + ESub
!    enddo
!
!  end subroutine
!=============================================================================+
!  subroutine Hybrid_OrthoVolECalc(self, curbox, disp, E_Diff, accept)
!    implicit none
!    class(Pair_Hybrid), intent(inout) :: self
!    class(SimBox), intent(inout) :: curbox
!    type(OrthoVolChange), intent(inout) :: disp(:)
!    real(dp), intent(inOut) :: E_Diff
!    logical, intent(out) :: accept
!    integer :: iField
!    real(dp) :: ESub
!
!    accept = .true.
!    E_Diff = 0E0_dp
!    self%EDiff = 0E0_dp
!
!    do iField = 1, self%NFFields
!      call EnergyCalculator(self%ECalcIndx(iField)) % Method % OrthoVolECalc(curbox, disp, ESub, accept)
!      if(.not. accept) then
!        self%EDiff = 0E0_dp
!        return
!      endif
!      self%EDiff(iField) = ESub
!      E_Diff = E_Diff + ESub
!    enddo
!
!  end subroutine
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
  subroutine Hybrid_Update(self, accept)
    implicit none
    class(Pair_Hybrid), intent(inout) :: self
    logical, intent(in) :: accept

    integer :: iField

    ! Call the update function for each sub-forcefield, which will allow them to reset any internal information after a successful move.
    do iField = 1, self%NFFields
      call EnergyCalculator(self%ECalcIndx(iField)) % Method % Update(accept)
    enddo

    if(allocated(self%EDiff)) self%EDiff = 0E0_dp


  end subroutine
!=============================================================================+
end module
!=============================================================================+
