!==========================================================================================
! The Simple Regrowth Object
!==========================================================================================
module MolCon_SimpleRegrowth
  use CoordinateTypes, only: Perturbation, Addition
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: SimpleRegrowth
    real(dp), allocatable :: tempcoords(:, :)
    contains
!      procedure, public, pass :: Constructor => SimpleRegrowth_Constructor
      procedure, public, pass :: Prologue => SimpleRegrowth_Prologue
      procedure, public, pass :: GenerateConfig => SimpleRegrowth_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine SimpleRegrowth_Prologue(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(SimpleRegrowth), intent(inout) :: self
!    integer, intent(in) :: molType
    integer :: iType, maxAtoms


    allocate( self%tempcoords(1:3, 1:MolData(self%molType)%nAtoms) )
  end subroutine
!==========================================================================================
  subroutine SimpleRegrowth_GenerateConfig(self, trialBox, disp, probconstruct,insPoint )
    use Common_MolInfo, only: MolData, BondData, nMolTypes
    use MolSearch, only: FindBond
    use RandomGen, only: Generate_UnitSphere
    implicit none
    class(SimpleRegrowth), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
!    class(SimpleBox), intent(inout) :: trialBox
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)

    integer :: bondType, molType
    integer :: atm1, atm2, iDisp
    real(dp), intent(out) :: probconstruct
    real(dp), dimension(1:3) :: v1, v2
    real(dp) :: dx, dy, dz, r
    real(dp) :: prob

    probconstruct = 1E0_dp
    select type(disp)
!      class is(DisplacementNew)
!        molType = disp(1)%molType
      class is(Addition)
        molType = disp(1)%molType
      class default
        stop "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
    end select

    select case( MolData(molType)%nAtoms )
      case(1)
        atm1 = 1
        self%tempcoords(1, Atm1) = 0E0_dp
        self%tempcoords(2, Atm1) = 0E0_dp
        self%tempcoords(3, Atm1) = 0E0_dp

      case(2)
        atm1 = 1
        atm2 = 2
        self%tempcoords(1, Atm1) = 0E0_dp
        self%tempcoords(2, Atm1) = 0E0_dp
        self%tempcoords(3, Atm1) = 0E0_dp
        call FindBond(molType, atm1, atm2, bondType)
        call BondData(bondType) % bondFF % GenerateDist(r, prob)
        call Generate_UnitSphere(dx, dy, dz)
        self%tempcoords(1, Atm2) = r * dx
        self%tempcoords(2, Atm2) = r * dy
        self%tempcoords(3, Atm2) = r * dz
!    case(3)
!      Atm1 = pathArray(nType)%path(1, 1)
!      Atm2 = pathArray(nType)%path(1, 2)
!      Atm3 = pathArray(nType)%path(1, 3)
!      call FindBond(nType, Atm1, Atm2, bondType)
!      k_bond = bondData(bondType)%k_eq
!      r_eq = bondData(bondType)%r_eq
!      call GenerateBondLength(r, k_bond, r_eq, Prob)
!      call Generate_UnitSphere(dx, dy, dz)
!      v1(1) = -r*dx
!      v1(2) = -r*dy
!      v1(3) = -r*dz
!      tempcoords(1,atm2) = r * dx
!      tempcoords(2,atm2) = r * dy
!      tempcoords(3,atm2) = r * dz
!      call FindBond(nType, Atm2, Atm3, bondType)
!      k_bond = bondData(bondType)%k_eq
!      r_eq = bondData(bondType)%r_eq
!      call GenerateBondLength(r, k_bond, r_eq, Prob)
!      bendType = bendArray(nType,1)%bendType
!      call GenerateBendAngle(ang, bendType, Prob)
!      call Generate_UnitCone(v1, r, ang, v2)
!      tempcoords(1, atm3) = tempcoords(1,atm2) + v2(1) 
!      tempcoords(2, atm3) = tempcoords(2,atm2) + v2(2)
!      tempcoords(3, atm3) = tempcoords(3,atm2) + v2(3)
      case default

    end select

    select type(disp)
!      class is(DisplacementNew)
!        molType = disp(1)%molType
      class is(Addition)
        do iDisp = 1, MolData(molType)%nAtoms
          if( present(insPoint) ) then
            disp(iDisp)%x_new = self%tempcoords(1, iDisp) + insPoint(1)
            disp(iDisp)%y_new = self%tempcoords(2, iDisp) + insPoint(2)
            disp(iDisp)%z_new = self%tempcoords(3, iDisp) + insPoint(3)
          else
            disp(iDisp)%x_new = self%tempcoords(1, iDisp)
            disp(iDisp)%y_new = self%tempcoords(2, iDisp)
            disp(iDisp)%z_new = self%tempcoords(3, iDisp)
          endif

        enddo
    end select


  end subroutine
!==========================================================================================
end module
!==========================================================================================
