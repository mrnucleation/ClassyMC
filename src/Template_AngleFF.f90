!=============================================================================+
module Template_IntraAngle
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes

  type, public, extends(Intra_FF) :: Angle_FF
    real(dp) :: theta0
    contains
      procedure, pass :: Constructor 
      procedure, pass :: EFunc
      procedure, pass :: ComputeAngle
      procedure, pass :: DetailedECalc 
!      procedure, pass :: DiffECalc
      procedure, pass :: ProcessIO
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
!    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Angle_FF), intent(inout) :: self


  end subroutine
!=============================================================================+
  function EFunc(self, angle) result(E_Angle)
    implicit none
    class(Angle_FF), intent(inout) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Angle

    E_Angle = 0E0_dp
  end function
!=============================================================================+
  function ComputeAngle(self, curbox, atompos) result(theta)
    implicit none
    class(Angle_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp) :: theta
    
    real(dp) :: r1, rx1, ry1, rz1
    real(dp) :: r2, rx2, ry2, rz2
    rx1 = atompos(1, 1) - atompos(1, 2)
    ry1 = atompos(2, 1) - atompos(2, 2)
    rz1 = atompos(3, 1) - atompos(3, 2)
    call curbox % Boundary(rx1, ry1, rz1)
    r1 = sqrt(rx1*rx1 + ry1*ry1 + rz1*rz1)

    rx2 = atompos(1, 3) - atompos(1, 2)
    ry2 = atompos(2, 3) - atompos(2, 2)
    rz2 = atompos(3, 3) - atompos(3, 2)
    call curbox % Boundary(rx2, ry2, rz2)
    r2 = sqrt(rx2*rx2 + ry2*ry2 + rz2*rz2)
     
    theta = rx1*rx2 + ry1*ry2 + rz1*rz2
    theta = acos(theta/(r1*r2))

  end function
!=============================================================================+
  subroutine DetailedECalc(self, curbox, atompos, E_T, accept)
    implicit none
    class(Angle_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    real(dp), intent(in) :: atompos(:, :)
    logical, intent(out) :: accept

    real(dp) :: theta

    accept = .true.
    theta = self%ComputeAngle(curbox, atompos)
    E_T = self%EFunc(theta)

!    rx1 = atompos(1, 1) - atompos(1, 2)
!    ry1 = atompos(2, 1) - atompos(2, 2)
!    rz1 = atompos(3, 1) - atompos(3, 2)
!    call curbox % Boundary(rx1, ry1, rz1)
!    r1 = sqrt(rx1*rx1 + ry1*ry1 + rz1*rz1)

!    rx2 = atompos(1, 3) - atompos(1, 2)
!    ry2 = atompos(2, 3) - atompos(2, 2)
!    rz2 = atompos(3, 3) - atompos(3, 2)
!    call curbox % Boundary(rx2, ry2, rz2)
!    r2 = sqrt(rx2*rx2 + ry2*ry2 + rz2*rz2)
     
!    theta = rx1*rx2 + ry1*ry2 + rz1*rz2
!    theta = acos(theta/(r1*r2))

  end subroutine
!==========================================================================
!  subroutine GenerateDist(self, val, probgen)
!    implicit none
!    class(Angle_FF), intent(inout) :: self
!    real(dp), intent(out) :: val
!    real(dp), intent(out) :: probgen
!
!    val = 0E0_dp
!    probgen = 1E0_dp
!
!  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Angle_FF), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
end module
!=============================================================================+
