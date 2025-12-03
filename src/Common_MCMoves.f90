!======================================================
module MCMoveData
use VarPrecision
use MoveClassDef, only: MCMove

!  integer :: nMovesTypes = 0
  type MoveArray 
    class(MCMove), allocatable :: Move
  end type

  type(MoveArray), allocatable, target :: Moves(:)
  real(dp), allocatable :: MoveProb(:)

end module
!======================================================
