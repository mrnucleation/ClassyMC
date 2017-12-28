!==============================================================
module TrajData
  use TrajectoryTemplate, only: trajectory

  type TrjArray
    class(trajectory), allocatable:: traj
  end type

  type(TrjArray), allocatable, target  :: TrajArray(:)

end module
!==============================================================
