!====================================================================
!This module contains the EnergyFloor constraint that prevents the system from
!going below a given energy value
!====================================================================
module Constrain_EnergyFloor
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use SimpleSimBox, only: SimpleBox
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout
  use Units, only: outEngUnit

  type, public, extends(constraint) :: EnergyFloor
    integer :: boxID = -1
    integer :: eStyle = 1
    real(dp) :: E_Min = -1E50_dp
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => EnergyFloor_Constructor
      procedure, pass :: ProcessIO => EnergyFloor_ProcessIO
      procedure, pass :: PostEnergy => EnergyFloor_PostEnergy
      procedure, pass :: Prologue => EnergyFloor_Prologue
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine EnergyFloor_Constructor(self, boxID)
    use BoxData, only: BoxArray
    implicit none
    class(EnergyFloor), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: AllocateStat
    integer :: nMolMax


    self%boxID = boxID
    self%parent => BoxArray(boxID) % box 

    IF (AllocateStat /= 0) STOP "Allocation Error in EnergyFloor Constraint"
  end subroutine
!=====================================================================
  subroutine EnergyFloor_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(EnergyFloor), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine EnergyFloor_PostEnergy(self, trialBox, disp, E_Diff, accept)
    implicit none
    class(EnergyFloor), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(in) :: E_Diff
    logical, intent(out) :: accept
    integer :: iDisp, nNew
    real(dp) :: E_New


    select type(trialBox)
      class is(SimpleBox)
        E_New = trialBox % GetNewEnergy(E_Diff)
        select case(self%eStyle)
!          case(0) ! Total
          case(1) ! Per Mol
              nNew = trialBox % GetNewMolCount(disp)
              E_New = E_New/real(nNew, dp)
!          case(2) ! Per Atom
        end select
    end select

    if(E_New < self%E_Min) then
      accept = .false.
    else
      accept = .true.
    endif


  end subroutine
!=============================================================
  subroutine EnergyFloor_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ParallelVar, only: nout
    implicit none
    class(EnergyFloor), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 2, lineStat)
    select case(trim(adjustl(command)))
      case("total") ! Total
          self%eStyle = 0
      case("mol") ! Per Mol
          self%eStyle = 1
      case("atom") ! Per Atom
          self%eStyle = 2
    end select


    call GetXCommand(line, command, 3, lineStat)
    read(command, *) realVal
    self%E_Min = realVal*outEngUnit

  end subroutine
!====================================================================
!  subroutine EnergyFloor_Maintenance(self)
!    implicit none
!    class(EnergyFloor), intent(inout) :: self
!
!  end subroutine
!====================================================================
  subroutine EnergyFloor_Prologue(self)
    use ParallelVar, only: nout
    use Units, only: outEngUnit
    implicit none
    class(EnergyFloor), intent(inout) :: self
    logical :: accept


    write(nout, *) "Applying an Energy Floor:", self%E_Min/outEngUnit


  end subroutine
!=====================================================================
end module
!=====================================================================
