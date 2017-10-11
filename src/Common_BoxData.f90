!==============================================================
module BoxData
  use SimBoxDef
  use ConstraintTemplate

  type(SimBox), allocatable :: BoxArray(:)
  type(constrainArray), allocatable :: Constrain(:)

end module
!==============================================================
