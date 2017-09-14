module MoleculeTranslation
use CoordinateTypes
use MoveClassDef
use SimBoxDef, only: SimBox
use VarPrecision

  type, public, extends(MCMove) :: MolTranslate
    real(dp), :: max_dist
    type(Displacement), allocatable :: disp
    contains
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
  end type

 contains

  subroutine GeneratePosition(self, disp)
    implicit none
    class(MCMove), intent(in) :: self
    type(Displacement), intent(inout) :: disp
    real(dp) :: dx, dy, dz
      dx = self % max_dist(nType) * (2E0_dp*grnd() - 1E0_dp)
      dy = self % max_dist(nType) * (2E0_dp*grnd() - 1E0_dp)
      dz = self % max_dist(nType) * (2E0_dp*grnd() - 1E0_dp)
  end subroutine

  subroutine FullMove(self, trialBox)
    implicit none
    class(MCMove), intent(in) :: self
    class(SimBox), intent(inout) :: trialBox
    logical :: accept


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
       !Choose Molecule
    dx = self % max_dist(nType) * (2E0_dp*grnd() - 1E0_dp)
    dy = self % max_dist(nType) * (2E0_dp*grnd() - 1E0_dp)
    dz = self % max_dist(nType) * (2E0_dp*grnd() - 1E0_dp)
 
    !Check Constraint
    if(.not. accept) then
      return
    endif


    !Energy Calculation


    !Accept/Reject
    if(accept) then
      self % accpt = self % atmps + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(disp)

    endif

  end subroutine


end module
