!=======================================================
  module RandomGen
      use VarPrecision
      integer :: initseed = -1
!=======================================================
      contains
!=======================================================
      real(dp) function grnd()
        implicit none
        real(dp) :: r

        call RANDOM_NUMBER(r)
        grnd = r
       
      end function
!=======================================================
      subroutine sgrnd(seed)
      implicit none
      integer, intent(in) :: seed
      integer :: i,n
      integer, allocatable :: tempSeed(:)
      
      call RANDOM_SEED(size=n)      
      allocate(tempSeed(1:n))
      tempSeed = seed + 37 * (/ (i - 1, i = 1, n) /)

      call RANDOM_SEED(put=tempSeed)
     
      deallocate(tempSeed)
     
      end subroutine      
!=======================================================
      real(dp) function Gaussian() result(num)
      use ClassyConstants
      use VarPrecision
      implicit none
      real(dp) :: y1, w, x1, x2

      x1 = grnd()
      x2 = grnd()
      y1 = sqrt( -2.0E0_dp * log(x1) ) * cos(two_pi*x2)
      num = y1
      end function
!========================================================            
      subroutine Generate_UnitSphere(x,y,z)
      use VarPrecision
      implicit none
      real(dp), intent(out) :: x,y,z
      real(dp) :: u_12_sq, u1, u2
      
      u_12_sq = 2E0_dp
      do while(u_12_sq .ge. 1E0_dp)
       u1 = 2E0_dp * grnd() - 1E0_dp
       u2 = 2E0_dp * grnd() - 1E0_dp
       u_12_sq = u1 * u1 + u2 * u2
      enddo
 
      x = 2E0_dp * u1 * sqrt(1E0_dp - u_12_sq)
      y = 2E0_dp * u2 * sqrt(1E0_dp - u_12_sq)
      z = (1E0_dp - 2E0_dp * u_12_sq)
      
      end subroutine
!==========================================================================     
!     The purpose of this function is to generate a random position for an atom
!     given a fixed bond angle (bond_ang) with a given vector (v1) and a given distance (r2).
!     The torsional rotational angle is chosen randomly from 0 to 2pi.
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vector(v2) is calculated using a rotational matrix
 
      subroutine Generate_UnitCone(v1,r2,bond_ang,v2)
      use ClassyConstants      
      use CoordinateTypes
      implicit none
      real(dp), dimension(1:3), intent(in) :: v1
      real(dp), dimension(1:3), intent(out) :: v2
      real(dp), intent(in) :: bond_ang, r2
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj

      r1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()
        
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)      
      r_proj = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
        
      coeff1 = (r2/r1)*cos(bond_ang)
      coeff2 = (r2/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r1

      v2(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v2(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v2(3) = coeff1*v1(3)                       + coeff3*s_term*(r_proj*r_proj)

         
      end subroutine
  
!==========================================================================     
!     The purpose of this function is to generate two random position for two atoms
!     given two fixed bond angles (bond_ang1 and bond_ang2), one dihedral angle (dihed) which 
!     is the angle between two planes made by two angles with a given vector (v1) and two 
!     given distances (r2 and r3).
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vectors(v2 and v3) is calculated using a rotational matrix

      subroutine Generate_UnitPyramid(v1, r2, r3, bond_ang1, bond_ang2, dihed, v2, v3)
      use ClassyConstants      
      use CoordinateTypes
      implicit none
      real(dp), dimension(1:3), intent(in) :: v1
      real(dp), dimension(1:3), intent(out) :: v2, v3
      real(dp), intent(in) :: bond_ang1, bond_ang2, dihed, r2, r3
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj

      r1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()
        
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)      
      r_proj = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
        
      coeff1 = (r2/r1)*cos(bond_ang1)
      coeff2 = (r2/r_proj)*sin(bond_ang1)
      coeff3 = coeff2/r1

      v2(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v2(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v2(3) = coeff1*v1(3)                      + coeff3*s_term*(r_proj*r_proj)
  
      tors_angle = tors_angle + dihed
      s_term = sin(tors_angle)
      c_term = cos(tors_angle) 

      coeff1 = (r3/r1)*cos(bond_ang2)
      coeff2 = (r3/r_proj)*sin(bond_ang2)
      coeff3 = coeff2/r1
  
      v3(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v3(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v3(3) = coeff1*v1(3)                      + coeff3*s_term*(r_proj*r_proj)
         
      end subroutine
  
!==========================================================================     
!     The purpose of this function is similar to the UnitCone function however
!     in this case the torsional angle is fixed in addition to the distance and
!     bond angle which means there is exactly a single point in space which
!     corresponds to these constraints.
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x2,y2,z2)   w2=(-y2,x2,0)  w3=(-x2*z2, -y2*z2, x2^2 + y2^2)
      subroutine Generate_UnitTorsion(v1,v2,r3,bond_ang,tors_angle,v3)
      use ClassyConstants      
      use CoordinateTypes
      implicit none
      real(dp), dimension(1:3), intent(in) :: v1,v2
      real(dp), dimension(1:3), intent(out) :: v3
      real(dp), intent(in) :: bond_ang, tors_angle,r3
      real(dp) :: x1_u, y1_u, z1_u
      real(dp) :: x1_s, y1_s, z1_s
      real(dp) :: r2, rot_angle
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj

      r2 = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3)
      r2 = sqrt(r2)
      r_proj = sqrt(v2(1)*v2(1) + v2(2)*v2(2))

!            
      x1_u =  v1(1) - v2(1)
      y1_u =  v1(2) - v2(2)
      z1_u =  v1(3) - v2(3)

!     Calculate the v1 vector's w2 and w3 components in the new orthonormal framework.
!      x1_s =  ( v2(1)*x1_u + v2(2)*y1_u + v2(3) * z1_u)/(r2)
      y1_s =  (-v2(2)*x1_u + v2(1)*y1_u) / r_proj
      z1_s =  (-v2(1)*v2(3)*x1_u - v2(2)*v2(3)*y1_u + r_proj*r_proj*z1_u) / (r_proj*r2)


        
!     Calculate the torsional rotation angle for the new v3 vector from the v1 components
      rot_angle = atan2(z1_s, y1_s)
!      write(2,*) rot_angle*180E0/pi, tors_angle*180E0/pi, (rot_angle + tors_angle)*180E0/pi
      rot_angle = rot_angle + tors_angle

!     Rescale the angle      
!      do while(rot_angle .lt. 0E0) 
!        rot_angle = rot_angle + two_pi
!      enddo
!      do while(rot_angle .gt. two_pi) 
!        rot_angle = rot_angle - two_pi
!      enddo
        
      s_term = sin(rot_angle)
      c_term = cos(rot_angle)      

        
      coeff1 = (r3/r2)*cos(bond_ang)
      coeff2 = (r3/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r2

      v3(1) = coeff1*v2(1) - coeff2*c_Term*v2(2) - coeff3*s_Term*v2(1)*v2(3)
      v3(2) = coeff1*v2(2) + coeff2*c_Term*v2(1) - coeff3*s_Term*v2(2)*v2(3)
      v3(3) = coeff1*v2(3)                      + coeff3*s_Term*(r_proj*r_proj)

         
      end subroutine
!========================================================================== 
!     The purpose of this function is to generate three random position for three atoms
!     given three fixed bond angles (bond_ang1 and bond_ang2 and bond_ang3), 
!     two dihedral angle (dihed1 and dihed2) where dihed1 is the angle between two planes made by 
!     first and second angles and dihed2 is the angle between two planes made by second and third angles
!     with a given vector (v1) and three given distances (r2 and r3 and r4).
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vectors(v2 and v3) is calculated using a rotational matrix

      subroutine Generate_UnitTetrahedral(v1, r2, r3, r4, bond_ang1, bond_ang2, bond_ang3, dihed1, dihed2, v2, v3, v4)
      use ClassyConstants      
      use CoordinateTypes
      implicit none
      real(dp), dimension(1:3), intent(in) :: v1
      real(dp), dimension(1:3), intent(out) :: v2, v3, v4
      real(dp), intent(in) :: bond_ang1, bond_ang2, bond_ang3, dihed1, dihed2, r2, r3, r4
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term,c_term
      real(dp) :: coeff1,coeff2,coeff3  
      real(dp) :: r_proj
  
      r1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()
        
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)      
      r_proj = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
        
      coeff1 = (r2/r1)*cos(bond_ang1)
      coeff2 = (r2/r_proj)*sin(bond_ang1)
      coeff3 = coeff2/r1

      v2(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v2(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v2(3) = coeff1*v1(3)                      + coeff3*s_term*(r_proj*r_proj)
  
      tors_angle = tors_angle + dihed1
      s_term = sin(tors_angle)
      c_term = cos(tors_angle) 

      coeff1 = (r3/r1)*cos(bond_ang2)
      coeff2 = (r3/r_proj)*sin(bond_ang2)
      coeff3 = coeff2/r1
  
      v3(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v3(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v3(3) = coeff1*v1(3)                      + coeff3*s_term*(r_proj*r_proj)
  
      tors_angle = tors_angle + dihed2
      s_term = sin(tors_angle)
      c_term = cos(tors_angle) 

      coeff1 = (r4/r1)*cos(bond_ang3)
      coeff2 = (r4/r_proj)*sin(bond_ang3)
      coeff3 = coeff2/r1
  
      v4(1) = coeff1*v1(1) - coeff2*c_Term*v1(2) - coeff3*s_Term*v1(1)*v1(3)
      v4(2) = coeff1*v1(2) + coeff2*c_term*v1(1) - coeff3*s_term*v1(2)*v1(3)
      v4(3) = coeff1*v1(3)                      + coeff3*s_term*(r_proj*r_proj)
  
      end subroutine

!=======================================================
  function ListRNG(list, norm) result(bin)
    use ClassyConstants
    use VarPrecision
    implicit none
    real(dp), intent(in) :: list(:)
    real(dp), intent(in), optional :: norm
    integer :: bin, nSel
    real(dp) :: ran_num, intSum

    if( present(norm) ) then
      ran_num = grnd()*norm
    else
      ran_num = grnd()
    endif

    nSel = 1 
    intSum = list(1)
    do while(intSum < ran_num)
      nSel = nSel + 1 
      intSum = intSum + list(nSel)
    enddo

    bin = nSel
  end function
!=========================================================
end module
