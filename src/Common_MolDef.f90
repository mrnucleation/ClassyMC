!================================================================
module StructureTypes
  use VarPrecision 

  type AtomDef 
    character(len=5) :: Symb
    real(dp) :: mass
  end type

end module
!================================================================
module Common_MolDef
  use VarPrecision 

  type AtomDef 
    character(len=5) :: Symb
    real(dp) :: mass
  end type

  integer :: nMolTypes
  integer :: nAtomTypes

  type(AtomDef), allocatable :: AtomData(:)
  type(BondDef), allocatable :: bondData(:)


end module
!================================================================
