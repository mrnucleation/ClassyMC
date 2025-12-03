!=============================================================================+
module IntraBond_Ridgid
  use Template_IntraBond, only: Bond_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(Bond_FF) :: RidgidBond
!    real(dp) :: r0
    contains
      procedure, pass :: Constructor => RidgidBond_Constructor
      procedure, pass :: DetailedECalc => RidgidBond_DetailedECalc
      procedure, pass :: GenerateDist => RidgidBond_GenerateDist
      procedure, pass :: ProcessIO => RidgidBond_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine RidgidBond_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RidgidBond), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine RidgidBond_DetailedECalc(self, curbox, atompos, E_T, accept)
    implicit none
    class(RidgidBond), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    real(dp), intent(in) :: atompos(:, :)
    logical, intent(out) :: accept

    E_T = 0E0_dp
    accept = .true.
  end subroutine
!==========================================================================
  subroutine RidgidBond_GenerateTrial(self, beta, val, bounds )
    use RandomGen, only: grnd, Gaussian
    implicit none
    class(RidgidBond), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(in), optional :: bounds(1:2)
    real(dp) :: E_Gen, rMax
    real(dp) :: lb, ub

    if(.not. present(bounds)) then
      lb = 0E0_dp
      ub = self%r0
    else
      lb = bounds(1)
      ub = bounds(2)
    endif

    val = self%r0

  end subroutine
!==========================================================================
  subroutine RidgidBond_GenerateDist(self, beta, val, probgen, E_T)
    implicit none
    class(RidgidBond), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T

    val = self%r0
    probgen = 1E0_dp
    if(present(E_T)) then
      E_T = 0E0_dp
    endif

  end subroutine
!=============================================================================+
  subroutine RidgidBond_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand
    use Units, only: inLenUnit
    implicit none
    class(RidgidBond), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%r0
    self%r0 = self%r0* inLenUnit

    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the ridgid bond style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
