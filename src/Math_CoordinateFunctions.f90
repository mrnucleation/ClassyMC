!=========================================================================
!  Module containing common mathematical opperations
!=========================================================================
module Math_Coordinates
  use VarPrecision

  !----------------------------------------------------------------------
  contains

  !----------------------------------------------------------------------
  subroutine GetPerpendicularVector(v1, v2, phi)
     !Produces a vector, v2, that is perpendicular to the input v1 vector.
     ! v1 => The Input Reference Vector which the user needs a perpendicular vector towards
     ! v2 => The Output Vector that will be perpendicular to
     ! phi => Optional Rotational angle. Since
    use ClassyConstants, only: pi, two_pi
    use RandomGen, only: grnd
    implicit none
    real(dp), dimension(1:3), intent(in) :: v1
    real(dp), dimension(1:3), intent(out) :: v2
    real(dp), intent(in), optional :: phi
    real(dp) :: tors_angle
    real(dp), parameter :: bond_ang = 0.5E0_dp*pi
    real(dp) :: r1, r_proj, coeff2, coeff3
    real(dp) :: s_term, c_term

    r1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
    r1 = sqrt(r1)

    if( present(phi) ) then
      tors_angle = phi
    else
      tors_angle = two_pi*grnd()
    endif
        
    s_term = sin(tors_angle)
    c_term = cos(tors_angle)      
    r_proj = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
        
    !Since the bond_ang is fixed at pi/2 the first coefficient drops to 0 and the
    !second's sine term is fixed at 1.
!    coeff1 = (1E0_dp/r1)*cos(bond_ang)
!    coeff2 = (1E0_dp/r_proj)*sin(bond_ang)
    coeff2 = 1E0_dp/r_proj
    coeff3 = coeff2/r1
    v2(1) = -coeff2*c_term*v1(2) - coeff3*s_term*v1(1)*v1(3)
    v2(2) = +coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
    v2(3) =                      + coeff3*s_term*(r_proj*r_proj)

  end subroutine
  !----------------------------------------------------------------------
  !     The purpose of this function is to generate a vector for an atom
  !     given a fixed bond dist, angle, and torsion with a given vector (v1).
  !     The coordinate is created using a relative orthonormal framework given by these vectors
  !     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
  !     Where w1 is a normalized v1, w2 is perpendicular to v1, and w3 is generated 
  !     by the cross product w3 = (w1 x w2).
  !     Using these vectors, the new vector(v2) is calculated using a rotational matrix
   subroutine Generate_RelativeVector(v1,r2,bond_ang,phi,v2)
      use ClassyConstants, only: two_pi, pi
      use CoordinateTypes
      use RandomGen, only: grnd
      implicit none
      real(dp), dimension(1:3), intent(in) :: v1
      real(dp), dimension(1:3), intent(out) :: v2
      real(dp), intent(in) :: bond_ang, r2, phi
      real(dp) :: r1
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj

      r1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
      r1 = sqrt(r1)

!      phi = two_pi*grnd()
        
      s_term = sin(phi)
      c_term = cos(phi)      
      r_proj = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
        
      coeff1 = (r2/r1)*cos(bond_ang)
      coeff2 = (r2/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r1

      v2(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v2(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v2(3) = coeff1*v1(3)                       + coeff3*s_term*(r_proj*r_proj)
         
  end subroutine
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
     !Computes the dihedral angle made by three different bond angles.
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
