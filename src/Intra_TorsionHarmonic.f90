!=============================================================================+
! Torsional Functional form for OPLS, Trappe, and other equivalent Torsional
! styles.
!=============================================================================+
module IntraTorsion_Harmonic
  use Template_IntraTorsion, only: Torsion_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Units, only: inAngUnit

  type, public, extends(Torsion_FF) :: HarmonicTorsion
    integer :: nParameters
    real(dp) :: k0, ang0
    contains
      procedure, pass :: Constructor => HarmonicTorsion_Constructor
      procedure, pass :: EFunc => HarmonicTorsion_EFunc
!      procedure, pass :: DetailedECalc => HarmonicTorsion_DetailedECalc
      procedure, pass :: GenerateDist => HarmonicTorsion_GenerateDist
      procedure, pass :: ProcessIO => HarmonicTorsion_ProcessIO
  end type

  contains
!=============================================================================+
  subroutine HarmonicTorsion_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(HarmonicTorsion), intent(inout) :: self


  end subroutine
!============================================================================
  subroutine HarmonicTorsion_DiffECalc(self, curbox, disp, E_Diff, accept)
    implicit none
    class(HarmonicTorsion), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    accept = .true.
    E_Diff = 0E0_dp

  end subroutine
!==========================================================================
  subroutine HarmonicTorsion_GenerateDist(self, beta, val, probgen, E_T)
    use ClassyConstants, only: two_pi
    use RandomGen, only: grnd
    implicit none
    class(HarmonicTorsion), intent(inout) :: self
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
  function HarmonicTorsion_EFunc(self, angle) result(E_Tors)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(HarmonicTorsion), intent(in) :: self
    real(dp), intent(in) :: angle
    real(dp) :: E_Tors


    E_Tors = 0.5E0_dp * self%k0 * (angle - self%ang0)**2

  end function
!=============================================================================+
  subroutine HarmonicTorsion_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    use Units, only: inLenUnit, inEngUnit
    implicit none
    class(HarmonicTorsion), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%ang0
    self%ang0 = self%ang0*inAngUnit

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) self%k0
    self%k0 = self%k0*inEngUnit

  end subroutine
!=============================================================================+
end module
!=============================================================================+
