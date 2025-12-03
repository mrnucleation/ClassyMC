!=============================================================================+
!
!=============================================================================+
module IntraTorsion_CHARMM
  use Template_IntraTorsion, only: Torsion_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  type, public, extends(Torsion_FF) :: CHARMMTorsion
    integer :: nParameters
    real(dp) :: const = 0E0_dp
    real(dp), allocatable :: k(:)
    real(dp), allocatable :: delta(:)
    contains
      procedure, pass :: Constructor => CHARMMTorsion_Constructor
      procedure, pass :: DetailedECalc => CHARMMTorsion_DetailedECalc
      procedure, pass :: DiffECalc => CHARMMTorsion_DiffECalc
      procedure, pass :: EnergyFunction => CHARMMTorsion_EnergyFunction
      procedure, pass :: GenerateDist => CHARMMTorsion_GenerateDist
      procedure, pass :: ProcessIO => CHARMMTorsion_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine CHARMMTorsion_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(CHARMMTorsion), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine CHARMMTorsion_DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(CHARMMTorsion), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    E_T = 0E0_dp
    accept = .true.
  end subroutine
!============================================================================
  subroutine CHARMMTorsion_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(CHARMMTorsion), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine CHARMMTorsion_GenerateDist(self, val, probgen)
    implicit none
    class(CHARMMTorsion), intent(inout) :: self
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen

    val = self%ang0
    probgen = 1E0_dp

  end subroutine
!=============================================================================+
  function CHARMMTorsion_EnergyFunction(self, angle) result(E_Tors)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(CHARMMTorsion), intent(inout) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Tors

    E_Tors = self%const



  end function
!=============================================================================+
  subroutine CHARMMTorsion_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(CHARMMTorsion), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat
    integer :: iPar, nParmeter


    call CountCommands(line, nParameter)


    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%const
    allocate(self%k(1:nParameter-1))
    allocate(self%delta(1:nParameter-1))

    do iPar = 3, nParameters+1
    enddo
    self%ang0 = self%ang0*inAngUnit
    if(lineStat /= 0) then
      write(*,*) "Missing input rquired for the ridgid Torsion style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
