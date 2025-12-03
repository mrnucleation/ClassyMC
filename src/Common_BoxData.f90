!==============================================================
module BoxData
  use SimpleSimBox, only: SimpleBox

  type BxArray 
    class(SimpleBox), allocatable:: box
  end type

  type(BxArray), allocatable, target  :: BoxArray(:)

end module
!==============================================================
