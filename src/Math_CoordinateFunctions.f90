!=========================================================================
!  Module containing common mathematical opperations
!=========================================================================
module Math
  use VarPrecision

  !----------------------------------------------------------------------
  contains
  !----------------------------------------------------------------------
  function AngleFromVectors(v1, v2) result(theta)
     !Standard dot product angle calculation for computing the angle between
     !vectors.
    implicit none
    real(dp), intent(in) :: v1(:), v2(:)
    real(dp) :: theta
    real(dp) :: dotprod, r1, r2

    r1 = norm2(v1)
    if(r1 <= 0.0) error stop "Vector of zero length has been passed to Vector function!"
    r2 = norm2(v2)
    if(r2 <= 0.0) error stop "Vector of zero length has been passed to Vector function!"
    dotprod = dot_product(v1, v2)
    theta = dotprod/(r1*r2)
    theta = acos(theta)

  end function
  !----------------------------------------------------------------------
  function DihedralAngle(ang1, ang2, ang3) result(dihed)
     !Computes the dihedral angle made by three different angles.
    implicit none
    real(dp), intent(in) :: ang1, ang2, ang3
    real(dp) :: dihed

    dihed = acos((cos(ang3) - cos(ang1)*cos(ang2))/(sin(ang1)*sin(ang2)))
  end function
  !----------------------------------------------------------------------
  function DihedralAngle_FromVecs(v1, v2, v3) result(dihed)
     !Computes the dihedral angle made by three different vectors
    implicit none
    real(dp), intent(in) :: v1(1:3), v2(1:3), v3(1:3)
    real(dp) :: dihed
    real(dp) :: theta12, theta13, theta23

    theta12 = AngleFromVectors(v1, v2)
    theta13 = AngleFromVectors(v1, v3)
    theta23 = AngleFromVectors(v2, v3)
    dihed = DihedralAngle(theta12, theta13, theta23)

  end function
  !----------------------------------------------------------------------
end module
!=========================================================================
