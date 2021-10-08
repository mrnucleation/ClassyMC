!================================================================
module StructureTypes
  use VarPrecision 
  use Template_MolConstructor, only: MolConstructor
  use Template_IntraBond, only: Bond_FF
  use Template_IntraAngle, only: Angle_FF
  use Template_IntraTorsion, only: Torsion_FF
  use Template_IntraMiscIntra, only: MiscIntra_FF


  !These containers are used to define all the various atom, bond, angle,etc.
  !types used in the systems.  Multiple bonds in the same molecule can use
  !the same type so this is more a species type definition rather than the definition
  !of indidual bond definition.
   
  !Defines atom types used for Intermolecular purposes.
  type AtomDef 
    character(len=5) :: Symb
    real(dp) :: mass
  end type

  !Defines traditional bond potentials
  type BondDef 
    class(Bond_FF), allocatable :: bondFF
  end type

  !Defines traditional bond angle potentials.
  type AngleDef 
    class(Angle_FF), allocatable :: angleFF
  end type

  !Defines both proper and improper torsional potentials.
  type TorsionDef 
    class(Torsion_FF), allocatable :: torsionFF
  end type

  !Defines an Intra Molecular Interaction that does not fall under
  !the other categories.  Used for things such as 1-5 non-bonded,
  !1-4 non-bonded, etc.
  type MiscDef 
    class(MiscIntra_FF), allocatable :: miscFF
  end type


  !The following containers contain information about specific bonds
  !within a molecule. :
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


  !This type holds information about all the various information
  !required to construct and simulate a molecule
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
!    integer, allocatable :: misc(:)
    type(MiscDef), allocatable :: miscdata(:)
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
!  integer :: nMiscTypes = -1
  integer :: mostAtoms = -1

  type(MolDef), allocatable :: MolData(:)
  type(AtomDef), allocatable :: AtomData(:)
  type(BondDef), allocatable :: BondData(:)
  type(AngleDef), allocatable :: AngleData(:)
  type(TorsionDef), allocatable :: TorsionData(:)
!  type(MiscDef), allocatable :: MiscData(:)

end module
!================================================================
