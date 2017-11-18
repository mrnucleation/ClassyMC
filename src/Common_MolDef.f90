!================================================================
module StructureTypes
  use VarPrecision 

  type AtomDef 
    character(len=5) :: Symb
    real(dp) :: mass
  end type

  type BondDef 
    integer :: mem1, mem2
  end type

  type AngleDef 
    integer :: mem1, mem2, mem3
  end type

end module
!================================================================
module Common_MolDef
  use VarPrecision 
  use StructureTypes

  integer :: nMolTypes
  integer :: nAtomTypes

  type(AtomDef), allocatable :: AtomData(:)
  type(BondDef), allocatable :: bondData(:)


end module
!================================================================
