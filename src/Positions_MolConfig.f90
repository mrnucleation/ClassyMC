!==============================================================
module MolConfig
  use CoordinateTypes, only: Displacement

  type SimpleMolCoords
    real(dp) :: x, y, z
  end type

!==============================================================
contains
!==============================================================
  subroutine SimpleMol_ConfigGen(disp)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    type(Displacement), intent(inout) :: disp(:)



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
        rosenTrial(iRosen)%x(Atm2) = r * dx
        rosenTrial(iRosen)%y(Atm2) = r * dy
        rosenTrial(iRosen)%z(Atm2) = r * dz
      case(3)
        Atm1 = pathArray(nType)%path(1, 1)
        Atm2 = pathArray(nType)%path(1, 2)
        Atm3 = pathArray(nType)%path(1, 3)
        call FindBond(nType, Atm1, Atm2, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        v1%x = -r*dx
        v1%y = -r*dy
        v1%z = -r*dz
        rosenTrial(iRosen)%x(atm2) = r * dx
        rosenTrial(iRosen)%y(atm2) = r * dy
        rosenTrial(iRosen)%z(atm2) = r * dz
        call FindBond(nType, Atm2, Atm3, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r, k_bond, r_eq, Prob)
        bendType = bendArray(nType,1)%bendType
  !          k_bend = bendData(bendType)%k_eq
  !          ang_eq = bendData(bendType)%ang_eq
  !          call GenerateBendAngle(ang, k_bend, ang_eq, Prob)
        call GenerateBendAngle(ang, bendType, Prob)
        call Generate_UnitCone(v1, r, ang, v2)
        rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
        rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
        rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
      case(4)
        Atm1 = pathArray(nType)%path(1, 1)
        Atm2 = pathArray(nType)%path(1, 2)
        Atm3 = pathArray(nType)%path(2, 1)
        Atm4 = pathArray(nType)%path(3, 1)
        call FindBond(nType, Atm1, Atm2, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r1, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        v1%x = -r1*dx
        v1%y = -r1*dy
        v1%z = -r1*dz
        rosenTrial(iRosen)%x(atm2) = r1 * dx
        rosenTrial(iRosen)%y(atm2) = r1 * dy
        rosenTrial(iRosen)%z(atm2) = r1 * dz
        call FindBond(nType, Atm2, Atm3, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r2, k_bond, r_eq, Prob)
        call FindBond(nType, Atm2, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r3, k_bond, r_eq, Prob)
        call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
        call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
        call FindAngle(nType, Atm3, Atm2, Atm4, bendType3)
        wBending(iRosen) = 0d0
        call GenerateTwoBranches(ang1, ang2, dihed, bendType, bendType2, bendType3, wBending(iRosen))
        call Generate_UnitPyramid(v1, r2, r3, ang1, ang2, dihed, v2, v3)
        rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
        rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
        rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
        rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
        rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
        rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
      case(5)
        Atm1 = pathArray(nType)%path(1, 1)
        Atm2 = pathArray(nType)%path(1, 2)
        Atm3 = pathArray(nType)%path(2, 1)
        Atm4 = pathArray(nType)%path(3, 1)
        Atm5 = pathArray(nType)%path(4, 1)
        call FindBond(nType, Atm1, Atm2, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r1, k_bond, r_eq, Prob)
        call Generate_UnitSphere(dx, dy, dz)
        v1%x = -r1*dx
        v1%y = -r1*dy
        v1%z = -r1*dz
        rosenTrial(iRosen)%x(atm2) = r1 * dx
        rosenTrial(iRosen)%y(atm2) = r1 * dy
        rosenTrial(iRosen)%z(atm2) = r1 * dz
        call FindBond(nType, Atm2, Atm3, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r2, k_bond, r_eq, Prob)
        call FindBond(nType, Atm2, Atm4, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r3, k_bond, r_eq, Prob)
        call FindBond(nType, Atm2, Atm5, bondType)
        k_bond = bondData(bondType)%k_eq
        r_eq = bondData(bondType)%r_eq
        call GenerateBondLength(r4, k_bond, r_eq, Prob)
        call FindAngle(nType, Atm1, Atm2, Atm3, bendType)
        call FindAngle(nType, Atm1, Atm2, Atm4, bendType2)
        call FindAngle(nType, Atm1, Atm2, Atm5, bendType3)
        call FindAngle(nType, Atm3, Atm2, Atm4, bendType4)
        call FindAngle(nType, Atm4, Atm2, Atm5, bendType5)
        call FindAngle(nType, Atm5, Atm2, Atm3, bendType6)
        wBending(iRosen) = 0d0
        call GenerateThreeBranches(ang1, ang2, ang3, dihed, dihed2, &
             bendType, bendType2, bendType3, bendType4, bendType5, bendType6, wBending(iRosen))
        call Generate_UnitTetrahedral(v1, r2, r3, r4, ang1, ang2, ang3, dihed, dihed2, v2, v3, v4)
        rosenTrial(iRosen)%x(atm3) = rosenTrial(iRosen)%x(atm2) + v2%x 
        rosenTrial(iRosen)%y(atm3) = rosenTrial(iRosen)%y(atm2) + v2%y
        rosenTrial(iRosen)%z(atm3) = rosenTrial(iRosen)%z(atm2) + v2%z
        rosenTrial(iRosen)%x(atm4) = rosenTrial(iRosen)%x(atm2) + v3%x 
        rosenTrial(iRosen)%y(atm4) = rosenTrial(iRosen)%y(atm2) + v3%y
        rosenTrial(iRosen)%z(atm4) = rosenTrial(iRosen)%z(atm2) + v3%z
        rosenTrial(iRosen)%x(atm5) = rosenTrial(iRosen)%x(atm2) + v4%x 
        rosenTrial(iRosen)%y(atm5) = rosenTrial(iRosen)%y(atm2) + v4%y
        rosenTrial(iRosen)%z(atm5) = rosenTrial(iRosen)%z(atm2) + v4%z
      case default
        write(*,*) "Error! Molecule has too many atoms for a simple regrowth"
        stop
      end select

  end subroutine
!==============================================================
end module
!==============================================================
