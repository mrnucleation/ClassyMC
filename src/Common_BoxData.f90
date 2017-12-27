!==============================================================
module BoxData
  use SimpleSimBox, only: SimpleBox
!  use ConstraintTemplate

  type BxArray 
    class(SimpleBox), allocatable:: box
  end type

  type(BxArray), allocatable, target  :: BoxArray(:)
!  type(constrainArray), allocatable :: Constrain(:)

end module
!==============================================================
