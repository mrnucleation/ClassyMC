!==============================================================
module BoxData
  use SimBoxDef
  use ConstraintTemplate

  type BxArray 
    class(SimBox), allocatable:: box
  end type

  type(BxArray), allocatable, target  :: BoxArray(:)
  type(constrainArray), allocatable :: Constrain(:)

end module
!==============================================================
