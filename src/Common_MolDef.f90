!================================================================
module StructureTypes
  use VarPrecision 


  type AtomDef 
    character(len=5) :: Symb
    real(dp) :: mass
  end type


  type BondMem
    integer :: bondType
    integer :: mem1, mem2
  end type


  type AngleMem 
    integer :: angType
    integer :: mem1, mem2, mem3
  end type


  type TorsMem
    integer :: torsType
    integer :: mem1, mem2, mem3, mem4
  end type


  type MolDef 
    integer :: nAtoms = 1
    integer, allocatable :: atomType(:)

    integer :: nBonds = 0
    type(BondMem), allocatable :: bond(:)

    integer :: nAngles = 0
    type(AngleMem), allocatable :: angle(:)

    integer :: nTors = 0
    type(TorsMem), allocatable :: torsion(:)
  end type



end module
!================================================================
module Common_MolInfo
  use VarPrecision 
  use StructureTypes

  integer :: nMolTypes = 1
  integer :: nAtomTypes = -1
  integer :: nBondTypes = -1


  type(MolDef), allocatable :: MolData(:)
  type(AtomDef), allocatable :: AtomData(:)
!  type(BondDef), allocatable :: BondData(:)



end module
!================================================================
