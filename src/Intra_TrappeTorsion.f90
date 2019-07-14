!=============================================================================+
!
!=============================================================================+
module IntraTorsion_TRAPPE
  use Template_IntraTorsion, only: Torsion_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  integer, parameter 
  type, public, extends(Torsion_FF) :: TRAPPETorsion
    real(dp) :: const = 0E0_dp
    real(dp), allocatable :: k(1:3)
    contains
      procedure, pass :: Constructor => TRAPPETorsion_Constructor
      procedure, pass :: DetailedECalc => TRAPPETorsion_DetailedECalc
      procedure, pass :: DiffECalc => TRAPPETorsion_DiffECalc
      procedure, pass :: EnergyFunction => TRAPPETorsion_EnergyFunction
      procedure, pass :: GenerateDist => TRAPPETorsion_GenerateDist
      procedure, pass :: ProcessIO => TRAPPETorsion_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine TRAPPETorsion_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(TRAPPETorsion), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine TRAPPETorsion_DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    E_T = 0E0_dp
    accept = .true.
  end subroutine
!============================================================================
  subroutine TRAPPETorsion_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine TRAPPETorsion_GenerateDist(self, beta, val, probgen)
    use ClassyConstants, only: two_pi
    use RandomGen, only: grnd
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp) :: E_Tors

    do while(.true.) 
      val = two_pi * grnd()
      E_Tors = self%EnergyFunction(val)

      
    enddo

  end subroutine
!=============================================================================+
  function TRAPPETorsion_EnergyFunction(self, angle) result(E_Tors)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Tors
    integer :: iTerm

    E_Tors = self%const
    E_Tors = E_Tors + self%k(1)*(1E0_dp+cos(angle))
    E_Tors = E_Tors + self%k(2)*(1E0_dp-cos(2E0_dp*angle))
    E_Tors = E_Tors + self%k(3)*(1E0_dp+cos(3E0_dp*angle))



  end function
!=============================================================================+
  subroutine TRAPPETorsion_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
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
