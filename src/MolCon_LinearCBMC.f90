!==========================================================================================
! Simple Regrowth Object
!==========================================================================================
module MolCon_LinearCBMC
  use CoordinateTypes, only: Perturbation, Addition
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: LinearCBMC
    real(dp), allocatable :: tempcoords(:, :)
    integer :: centralAtom = 1
    integer, allocatable :: sideAtoms(:)
    contains
!      procedure, public, pass :: Constructor => LinearCBMC_Constructor
      procedure, public, pass :: Prologue => LinearCBMC_Prologue
      procedure, public, pass :: GenerateConfig => LinearCBMC_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine LinearCBMC_Prologue(self)
    use Common_MolInfo, only: MolData, BondData, nMolTypes
    use MolSearch, only: FindBond
    use ParallelVar, only: nout
    implicit none
    class(LinearCBMC), intent(inout) :: self
!    integer, intent(in) :: molType
    integer :: iType, iBond, iAtom, atm, curMax, nAtoms
    integer, allocatable :: freq(:)


  end subroutine
!==========================================================================================
  subroutine LinearCBMC_GenerateConfig(self, trialBox, disp, probconstruct, insPoint)
    use Common_MolInfo, only: MolData, BondData, AngleData, nMolTypes
    use MolSearch, only: FindBond, FindAngle
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone
    implicit none
    class(LinearCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)

    integer :: bondType, angleType, molType
    integer :: atm1, atm2,atm3, iDisp, iAtom
    real(dp), intent(out) :: probconstruct
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r, theta
    real(dp) :: r1, r2
    real(dp) :: prob
    real(dp) :: ang1, ang2

    probconstruct = 1E0_dp
    select type(disp)
      class is(Addition)
        molType = disp(1)%molType
      class default
        stop "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
    end select

t


  end subroutine
!==========================================================================================
end module
!==========================================================================================
