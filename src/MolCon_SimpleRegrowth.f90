!==========================================================================================
! Simple Regrowth Object
!==========================================================================================
module MolCon_SimpleRegrowth
  use CoordinateTypes, only: Perturbation, Addition
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: SimpleRegrowth
    real(dp), allocatable :: tempcoords(:, :)
    integer :: centralAtom = 1
    integer, allocatable :: sideAtoms(:)
    contains
!      procedure, public, pass :: Constructor => SimpleRegrowth_Constructor
      procedure, public, pass :: Prologue => SimpleRegrowth_Prologue
      procedure, public, pass :: GenerateConfig => SimpleRegrowth_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine SimpleRegrowth_Prologue(self)
    use Common_MolInfo, only: MolData, BondData, nMolTypes
    use MolSearch, only: FindBond
    use ParallelVar, only: nout
    implicit none
    class(SimpleRegrowth), intent(inout) :: self
!    integer, intent(in) :: molType
    integer :: iType, iBond, iAtom, atm, curMax, nAtoms
    integer, allocatable :: freq(:)

    nAtoms = MolData(self%molType)%nAtoms
    allocate( self%tempcoords(1:3, 1:nAtoms) )

    !Locate the central atom by determining the atom that appears the most frequently
    !in the bond topology.
    if(MolData(self%molType)%nAtoms > 1) then
      allocate( freq(1:nAtoms) )  
      freq = 0
      do iBond = 1, size(MolData(self%molType)%bond)
        atm = MolData(self%molType)%bond(iBond)%mem1
        freq(atm) = freq(atm) + 1
        atm = MolData(self%molType)%bond(iBond)%mem2
        freq(atm) = freq(atm) + 1
      enddo

      atm = 1
      curMax = 0
      do iAtom = 1, nAtoms
        if(curMax < freq(iAtom)) then
          curMax = freq(iAtom)
          atm = iAtom
        endif
      enddo
      self%centralAtom = atm
      allocate( self%sideAtoms(1:nAtoms) )

      atm = 0
      do iAtom = 1, nAtoms
        if(iAtom /= self%centralAtom) then
          atm = atm + 1
          self%sideAtoms(atm) = iAtom
        endif
      enddo 
      deallocate(freq)
    endif
    write(nout, *) "Central Atom Identified:", self%centralAtom

  end subroutine
!==========================================================================================
  subroutine SimpleRegrowth_GenerateConfig(self, trialBox, disp, probconstruct,insPoint )
    use Common_MolInfo, only: MolData, BondData, AngleData, nMolTypes
    use MolSearch, only: FindBond, FindAngle
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone
    implicit none
    class(SimpleRegrowth), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
!    class(SimpleBox), intent(inout) :: trialBox
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
      case(3)
        Atm1 = self%sideAtoms(1)
        Atm2 = self%centralAtom
        Atm3 = self%sideAtoms(2)
        self%tempcoords(1, atm2) = 0E0_dp
        self%tempcoords(2, atm2) = 0E0_dp
        self%tempcoords(3, atm2) = 0E0_dp
        call FindBond(molType, atm1, atm2, bondType)
        call BondData(bondType) % bondFF % GenerateDist(r, prob)
        call Generate_UnitSphere(dx, dy, dz)
        v1(1) = r*dx
        v1(2) = r*dy
        v1(3) = r*dz
        self%tempcoords(1, atm1) = r * dx
        self%tempcoords(2, atm1) = r * dy
        self%tempcoords(3, atm1) = r * dz
        call FindBond(molType, atm2, atm3, bondType)
        call BondData(bondType) % bondFF % GenerateDist(r, prob)
        call FindAngle(molType, atm1, atm2, atm3, angleType)
        call AngleData(angleType) % angleFF % GenerateDist(theta, prob)
        call Generate_UnitCone(v1, r, theta, v2)
        self%tempcoords(1, atm3) = v2(1) 
        self%tempcoords(2, atm3) = v2(2)
        self%tempcoords(3, atm3) = v2(3)
!      case(4)
!        Atm1 = self%sideAtoms(1)
!        Atm2 = self%centralAtom
!        Atm3 = self%sideAtoms(2)
!        Atm4 = self%sideAtoms(3)
!        self%tempcoords(1, atm2) = 0E0_dp
!        self%tempcoords(2, atm2) = 0E0_dp
!        self%tempcoords(3, atm2) = 0E0_dp
!
!        call FindBond(molType, atm1, atm2, bondType)
!        call BondData(bondType) % bondFF % GenerateDist(r1, prob)
!        call Generate_UnitSphere(dx, dy, dz)
!        v1(1) = -r1*dx
!        v1(2) = -r1*dy
!        v1(3) = -r1*dz
!        self%tempcoords(1, atm2) = r1 * dx
!        self%tempcoords(2, atm2) = r1 * dy
!        self%tempcoords(3, atm2) = r1 * dz
!        call FindBond(nType, Atm2, Atm3, bondType)
!        call BondData(bondType) % bondFF % GenerateDist(r2, prob)
!
!        call FindBond(nType, Atm2, Atm4, bondType)
!        call BondData(bondType) % bondFF % GenerateDist(r3, prob)
!        call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
!        call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
!        call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
!          wBending(iRosen) = 0d0
!          call GenerateTwoBranches(ang1, ang2, dihed, bendType, bendType2, bendType3, wBending(iRosen))
!          call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
!          rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
!          rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
!          rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
!          rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
!          rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
!          rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
      case default
        stop "Simple regrowth is not valid for this many atoms"
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
