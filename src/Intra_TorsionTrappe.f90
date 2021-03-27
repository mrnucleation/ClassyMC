!=============================================================================+
! Torsional Functional form for OPLS, Trappe, and other equivalent Torsional
! styles.
!=============================================================================+
module IntraTorsion_TRAPPE
  use Template_IntraTorsion, only: Torsion_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  type, public, extends(Torsion_FF) :: TRAPPETorsion
    integer :: nParameters
!    real(dp) :: const = 0E0_dp
    real(dp), allocatable :: k(:)
    contains
      procedure, pass :: Constructor => TRAPPETorsion_Constructor
      procedure, pass :: EFunc => TRAPPETorsion_EFunc
!      procedure, pass :: DetailedECalc => TRAPPETorsion_DetailedECalc
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
  subroutine TRAPPETorsion_GenerateDist(self, beta, val, probgen, E_T)
    use ClassyConstants, only: two_pi
    use RandomGen, only: grnd
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T
    real(dp) :: E_Tors

    do
      val = two_pi * grnd()
      E_Tors = self%EFunc(val)
      probgen = exp(-beta*E_Tors)
      if(probgen < grnd()) then
        exit
      endif
    enddo

    if(present(E_T)) then     
      E_T = E_Tors
    endif
  end subroutine
!=============================================================================+
  function TRAPPETorsion_EFunc(self, angle) result(E_Tors)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Tors
    integer :: iTerm
    real(dp) :: signterm

    E_Tors = self%k(0)
    do iTerm = 1, self%nParameters-1
      if( mod(iTerm, 2) == 0) then
        signterm = -1E0_dp
      else
        signterm = 1E0_dp
      endif
      E_Tors = E_Tors + self%k(iTerm)*(1.0E0_dp+signterm*cos(real(iTerm,dp)*angle))
    enddo

  end function
!=============================================================================+
  subroutine TRAPPETorsion_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand, CountCommands
    use Units, only: inLenUnit, inEngUnit
    implicit none
    class(TRAPPETorsion), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat
    integer :: iPar, nParameter


    call CountCommands(line, nParameter)
    !Subtract 1 because the first element of the line is the "Trappe" command
    self%nParameters = nParameter-1

    allocate(self%k(0:self%nParameters-1))
!    allocate(self%delta(1:nParameter-1))

    do iPar = 1, self%nParameters
        call GetXCommand(line, command, iPar+1, lineStat)
        if(lineStat /= 0) then
          write(0,*) "Missing input rquired for the Trappe Torsion style"
          error stop
        endif
        read(command, *) self%k(iPar-1)
        self%k(iPar-1) = self%k(iPar-1)*inEngUnit
    enddo
!    write(*,*) self%k(0: self%nParameters-1)
  end subroutine
!=============================================================================+
end module
!=============================================================================+
