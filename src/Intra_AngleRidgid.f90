!=============================================================================+
module IntraAngle_Ridgid
  use Template_IntraAngle, only: Angle_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  type, public, extends(Angle_FF) :: RidgidAngle
!    real(dp) :: r0
!    real(dp) :: k0
    contains
      procedure, pass :: Constructor => RidgidAngle_Constructor
!      procedure, pass :: DetailedECalc => RidgidAngle_DetailedECalc
      procedure, pass :: GenerateDist => RidgidAngle_GenerateDist
      procedure, pass :: ProcessIO => RidgidAngle_ProcessIO
  end type

  contains
!=============================================================================+
  function EFunc(self, angle) result(E_Angle)
    implicit none
    class(RidgidAngle), intent(inout) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Angle

    E_Angle = 0E0_dp
  end function
!=============================================================================+
  subroutine RidgidAngle_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RidgidAngle), intent(inout) :: self


  end subroutine
!==========================================================================
  subroutine RidgidAngle_GenerateTrial(self, beta, val, bounds )
    use RandomGen, only: grnd, Gaussian
    use ClassyConstants, only: two_pi
    implicit none
    class(RidgidAngle), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(in), optional :: bounds(1:2)
    real(dp) :: EGen
    real(dp) :: lb, ub

    if(.not. present(bounds)) then
      lb = 0E0_dp
      ub = two_pi
    else
      lb = bounds(1)
      ub = bounds(2)
    endif
    val = self%theta0
  end subroutine
!==========================================================================
  subroutine RidgidAngle_GenerateDist(self, beta, val, probgen, E_T)
    implicit none
    class(RidgidAngle), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T

!    val = self%theta0
    call self%GenerateTrial(beta, val)
    probgen = 1E0_dp
    if(present(E_T)) then     
      E_T = 0E0_dp
    endif

  end subroutine
!=============================================================================+
  subroutine RidgidAngle_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(RidgidAngle), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%theta0
    self%theta0 = self%theta0*inAngUnit

    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the ridgid Angle style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
