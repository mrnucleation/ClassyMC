!=============================================================================+
module Template_IntraBond
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use Template_Intra_FF, only: Intra_FF
  use CoordinateTypes

  type, public, extends(Intra_FF) :: Bond_FF
    real(dp) :: r0
    contains
      procedure, pass :: Constructor 
      procedure, pass :: ComputeProb
      procedure, pass :: ComputeBond
      procedure, pass :: EFunc
      procedure, pass :: DetailedECalc 
!      procedure, pass :: GenerateDist
      procedure, pass :: GenerateReverseDist
      procedure, pass :: ProcessIO
  end type

  contains
!=============================================================================+
  function ComputeBond(self, curbox, atompos) result(r)
!    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Bond_FF), intent(in) :: self
    class(simBox), intent(in) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp) :: rx, ry, rz, r

    rx = atompos(1, 1) - atompos(1, 2)
    ry = atompos(2, 1) - atompos(2, 2)
    rz = atompos(3, 1) - atompos(3, 2)
    call curbox % Boundary(rx, ry, rz)
    r = sqrt(rx*rx + ry*ry + rz*rz)
  end function
!==========================================================================
  function ComputeProb(self, beta, val) result(probgen)
    implicit none
    class(Bond_FF), intent(in) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: val
    real(dp) :: probgen
    real(dp) :: E_Val

    E_Val = self%EFunc(val)
    probgen = exp(-beta*E_Val)
  end function
!=============================================================================+
  function EFunc(self, dist) result(E_Bond)
!    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Bond_FF), intent(in) :: self
    real(dp), intent(in) :: dist
    real(dp) :: E_Bond

    E_Bond = 0E0_dp
  end function
!=============================================================================+
  subroutine Constructor(self)
!    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Bond_FF), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, atompos, E_T, accept)
    implicit none
    class(Bond_FF), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    real(dp) :: r

    accept = .true.
    r = self%ComputeBond(curbox, atompos)
    E_T = self%EFunc(r)

  end subroutine
!==========================================================================
!  subroutine GenerateDist(self, val, probgen)
!    implicit none
!    class(Bond_FF), intent(inout) :: self
!    real(dp), intent(out) :: val
!    real(dp), intent(out) :: probgen
!
!    val = 0E0_dp
!    probgen = 1E0_dp
!
!  end subroutine
!==========================================================================
  subroutine GenerateReverseDist(self, curbox, atompos, probgen)
    use RandomGen, only: grnd, Gaussian
    use ClassyConstants, only: two_pi
    implicit none
    class(Bond_FF), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(out) :: probgen

    logical :: accept
    integer :: iAtom
    real(dp) :: beta
    real(dp) :: E_Bond, reducepos(1:3,1:2)

    do iAtom = 1,2
      reducepos(1:3, iAtom) = atompos(1:3, iAtom) - atompos(1:3, 1)
      call curbox%Boundary(reducepos(1, iAtom), reducepos(2, iAtom), reducepos(3, iAtom))
    enddo

    beta = curbox%beta
    call self%DetailedECalc(curbox, reducepos(1:3,1:2), E_Bond, accept)
    probgen = exp(-beta*E_Bond)
  end subroutine
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(Bond_FF), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
end module
!=============================================================================+
