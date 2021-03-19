!================================================================
module StructureTypes
  use VarPrecision 
  use Template_MolConstructor, only: MolConstructor
  use Template_IntraBond, only: Bond_FF
  use Template_IntraAngle, only: Angle_FF
  use Template_IntraTorsion, only: Torsion_FF


  type AtomDef 
    character(len=5) :: Symb
    real(dp) :: mass
  end type

  type BondDef 
    class(Bond_FF), allocatable :: bondFF
  end type

  type AngleDef 
    class(Angle_FF), allocatable :: angleFF
  end type

  type TorsionDef 
    class(Torsion_FF), allocatable :: torsionFF
  end type

  type MiscDef 
    class(*), allocatable :: miscFF
  end type

  type BondMem
    integer :: bondType
    integer :: mem1, mem2
  end type


  type AngleMem 
    integer :: angleType
    integer :: mem1, mem2, mem3
  end type

  type TorsMem
    integer :: torsType
    integer :: mem1, mem2, mem3, mem4
  end type

  type MolDef 
    logical :: ridgid = .true.
    integer :: nAtoms = 1

    class(MolConstructor), allocatable :: molConstruct
    integer, allocatable :: atomType(:)

    integer :: nBonds = 0
    type(BondMem), allocatable :: bond(:)

    integer :: nAngles = 0
    type(AngleMem), allocatable :: angle(:)

    integer :: nTors = 0
    type(TorsMem), allocatable :: torsion(:)

    integer :: nMisc = 0
    integer, allocatable :: misc(:)
  end type


end module
!================================================================
module Common_MolInfo
  use VarPrecision 
  use StructureTypes

  integer :: nMolTypes = 1
  integer :: nAtomTypes = -1
  integer :: nBondTypes = -1
  integer :: nAngleTypes = -1
  integer :: nTorsionTypes = -1
  integer :: nMiscTypes = -1

  integer :: mostAtoms = -1

  type(MolDef), allocatable :: MolData(:)
  type(AtomDef), allocatable :: AtomData(:)
  type(BondDef), allocatable :: BondData(:)
  type(AngleDef), allocatable :: AngleData(:)
  type(TorsionDef), allocatable :: TorsionData(:)
  type(MiscDef), allocatable :: MiscData(:)

end module
!================================================================
