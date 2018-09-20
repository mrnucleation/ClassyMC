!==========================================================================================
! The Simple Regrowth Object
!==========================================================================================
module MolCon_SimpleRegrowth
  use CoordinateTypes, only: Displacement, Perturbation
  use MasterTemplate, only: classyClass
  use VarPrecision


  type, public, extends(MolConstructor) :: SimpleRegrowth
    integer :: molType
    real(dp), allocatable :: tempcoords(:, :)
    contains
      procedure, public, pass :: Constructor => SimpleRegrwoth_Constructor
      procedure, public, pass :: GenerateConfig => SimpleRegrwoth_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine SimpleRegrowth_Constructor(self, molType)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(MolConstructor), intent(inout) :: self
    integer, intent(in) :: molType
    integer :: iType, maxAtoms

    self%molType = molType

    allocate( tempcoords(1:3, 1:MolData(molType)%nAtoms) )
  end subroutine
!==========================================================================================
  subroutine SimpleRegrowth_GenerateConfig(self, trialBox, disp, probconstruct)
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimpleBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct = 1E0_dp
    real(dp), dimension(1:3) :: v1, v2
    real(dp) :: dx, dy, dz, r

    select case(nAtoms(nType))
      case(1)
        continue
      case(2)
        Atm1 = pathArray(nType)%path(1, 1)
        Atm2 = pathArray(nType)%path(1, 2)        
        bondType = bondArray(nType,1)%bondType
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        tempcoords(1,Atm2) = r * dx
        tempcoords(2,Atm2) = r * dy
        tempcoords(3,Atm2) = r * dz
    case(3)
      Atm1 = pathArray(nType)%path(1, 1)
      Atm2 = pathArray(nType)%path(1, 2)
      Atm3 = pathArray(nType)%path(1, 3)
      call FindBond(nType, Atm1, Atm2, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      call GenerateBondLength(r, k_bond, r_eq, Prob)
      call Generate_UnitSphere(dx, dy, dz)
      v1(1) = -r*dx
      v1(2) = -r*dy
      v1(3) = -r*dz
      tempcoords(1,atm2) = r * dx
      tempcoords(2,atm2) = r * dy
      tempcoords(3,atm2) = r * dz
      call FindBond(nType, Atm2, Atm3, bondType)
      k_bond = bondData(bondType)%k_eq
      r_eq = bondData(bondType)%r_eq
      call GenerateBondLength(r, k_bond, r_eq, Prob)
      bendType = bendArray(nType,1)%bendType
      call GenerateBendAngle(ang, bendType, Prob)
      call Generate_UnitCone(v1, r, ang, v2)
      tempcoords(1, atm3) = tempcoords(1,atm2) + v2(1) 
      tempcoords(2, atm3) = tempcoords(2,atm2) + v2(2)
      tempcoords(3, atm3) = tempcoords(3,atm2) + v2(3)
      case default

    end select


  end subroutine
!==========================================================================================
end module
!==========================================================================================
