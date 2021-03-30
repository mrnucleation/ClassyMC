!=============================================================================+
module Template_IntraTorsion
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes

  type, public, extends(Intra_FF) :: Torsion_FF
    contains
      procedure, pass :: Constructor 
      procedure, pass :: ComputeTors
      procedure, pass :: ComputeProb
      procedure, pass :: EFunc
      procedure, pass :: DetailedECalc 
      procedure, pass :: GenerateDist
      procedure, pass :: GenerateReverseDist
      procedure, pass :: ProcessIO
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
!    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Torsion_FF), intent(inout) :: self


  end subroutine
!=============================================================================+
  function ComputeTors(self, curbox, atompos) result(angle)
    implicit none
    class(Torsion_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp) :: angle

    integer :: iAtom
    real(dp) :: p1(1:3)
    real(dp) :: x12, y12, z12
    real(dp) :: x23, y23, z23
    real(dp) :: x34, y34, z34
    real(dp) :: vx1, vy1, vz1
    real(dp) :: vx2, vy2, vz2
    real(dp) :: vx3, vy3, vz3
    real(dp) :: r1, r3
    real(dp) :: dot1, dot2


    !Our goal here is to construct a sequence of vectors that are orthogonal
    !to each other with one of the vectors lined up along the central bond, the 2-3 bond,
    !and one of the orthoganal vectors lined up on top of the 1-2 bond. 
    !Once we have that we can simply take the dot product between this new
    !vector system and the remaining 3-4 vector to get the torsional angle. 



    x12 = atompos(1, 2) - atompos(1, 1)
    y12 = atompos(2, 2) - atompos(2, 1)
    z12 = atompos(3, 2) - atompos(3, 1)
    call curbox % Boundary(x12, y12, z12)
 
    x23 = atompos(1, 3) - atompos(1, 2) 
    y23 = atompos(2, 3) - atompos(2, 2) 
    z23 = atompos(3, 3) - atompos(3, 2) 
    call curbox % Boundary(x23, y23, z23)

    x34 = atompos(1, 4) - atompos(1, 3)
    y34 = atompos(2, 4) - atompos(2, 3)
    z34 = atompos(3, 4) - atompos(3, 3)
    call curbox % Boundary(x34, y34, z34)

!   Calculate the vector, v1 = <r12 x r23>
    vx1 =   y12*z23 - z12*y23
    vy1 = -(x12*z23 - z12*x23)
    vz1 =   x12*y23 - y12*x23
   
!   Calculate the vector, v2 = <r3 x r2>
    vx2 =   y23*z34 - z23*y34
    vy2 = -(x23*z34 - z23*x34)
    vz2 =   x23*y34 - y23*x34   

!   Calculate the vector, v3 = <v1 x r2> to create an orthonormal framework
    vx3 =   vy1*z23 - vz1*y23
    vy3 = -(vx1*z23 - vz1*x23)
    vz3 =   vx1*y23 - vy1*x23  

    !Normalize the v vectors.
    r1 = vx1**2 + vy1**2 + vz1**2
    r3 = vx3**2 + vy3**2 + vz3**2
    r1 = sqrt(r1)
    r3 = sqrt(r3)
    
    !Take the dot with respect to the remaining vectors
    !and compute the angle.  Keep in mind we want the angle
    !to span (0, 2pi)
    dot1 = vx1*vx2 + vy1*vy2 + vz1*vz2
    dot2 = vx2*vx3 + vy2*vy3 + vz2*vz3            
    dot1 = dot1/(r1)
    dot2 = dot2/(r3)    
    angle = atan2(dot2,dot1)


  end function
!=============================================================================+
  function EFunc(self, angle) result(E_Torsion)
    implicit none
    class(Torsion_FF), intent(in) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Torsion

    E_Torsion = 0E0_dp
  end function
!==========================================================================
  function ComputeProb(self, beta, val) result(probgen)
    implicit none
    class(Torsion_FF), intent(in) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: val
    real(dp) :: probgen
    real(dp) :: E_Val

    E_Val = self%EFunc(val)
    probgen = exp(-beta*E_Val)
  end function
!=============================================================================+
  subroutine DetailedECalc(self, curbox, atompos, E_T, accept)
    implicit none
    class(Torsion_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    real(dp), intent(in) :: atompos(:, :)
    logical, intent(out) :: accept
    real(dp) :: angle

    accept = .true.
    angle = self%ComputeTors(curbox, atompos)
    E_T = self%EFunc(angle)
  end subroutine
!==========================================================================
  subroutine GenerateDist(self, beta, val, probgen, E_T)
    implicit none
    class(Torsion_FF), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T

    val = 0E0_dp
    E_T = 0E0_dp
    probgen = 1E0_dp

  end subroutine
!==========================================================================
  subroutine GenerateReverseDist(self, curbox, atompos, probgen)
    use RandomGen, only: grnd, Gaussian
    use ClassyConstants, only: two_pi
    implicit none
    class(Torsion_FF), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(out) :: probgen

    logical :: accept
    integer :: iAtom
    real(dp) :: E_Torsion
    real(dp) :: beta
    real(dp) :: reducepos(1:3,1:4)

    do iAtom = 1,4
      reducepos(1:3, iAtom) = atompos(1:3, iAtom) - atompos(1:3, 1)
      call curbox%Boundary(reducepos(1, iAtom), reducepos(2, iAtom), reducepos(3, iAtom))
    enddo



    beta = curbox%beta
    call self%DetailedECalc(curbox, reducepos(1:3,1:4), E_Torsion, accept)
    probgen = exp(-beta*E_Torsion)
  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Torsion_FF), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
end module
!=============================================================================+
