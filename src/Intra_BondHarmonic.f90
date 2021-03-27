!=============================================================================+
module IntraBond_Harmonic
  use Template_IntraBond, only: Bond_FF
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(Bond_FF) :: HarmonicBond
!    real(dp) :: r0
    real(dp) :: k0
    contains
      procedure, pass :: EFunc => HarmonicBond_EFunc
      procedure, pass :: Constructor => HarmonicBond_Constructor
!      procedure, pass :: DetailedECalc => HarmonicBond_DetailedECalc
      procedure, pass :: GenerateDist => HarmonicBond_GenerateDist
      procedure, pass :: GenerateReverseDist => HarmonicBond_GenerateReverseDist
      procedure, pass :: ProcessIO => HarmonicBond_ProcessIO
  end type

  contains
!=============================================================================+
  function HarmonicBond_EFunc(self, dist) result(E_Bond)
    implicit none
    class(HarmonicBond), intent(inout) :: self
    real(dp), intent(in) :: dist
    real(dp) :: E_Bond

    E_Bond = 0.5E0_dp * self%k0 * (dist - self%r0)**2
  end function
!=============================================================================+
  subroutine HarmonicBond_Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(HarmonicBond), intent(inout) :: self


  end subroutine
!==========================================================================
  subroutine HarmonicBond_GenerateDist(self, beta, val, probgen, E_T)
    use RandomGen, only: grnd, Gaussian
    implicit none
    class(HarmonicBond), intent(inout) :: self
    real(dp), intent(in) :: beta
    real(dp), intent(out) :: val
    real(dp), intent(out) :: probgen
    real(dp), intent(out), optional :: E_T
    real(dp) :: E_Gen, rMax

    val = 0E0_dp
    probgen = 1E0_dp
    rMax = self%r0+5.0/sqrt(self%k0*beta)
    do
      val = Gaussian()
      val = self%r0 + val/sqrt(self%k0*beta)
      if( (val < 0.0E0_dp) .or. (val > self%r0+5.0/sqrt(self%k0*beta))  ) then
        cycle
      endif
      probgen = (val/rMax)**2
      if(probgen < grnd()) then
        exit
      endif
    enddo
    E_Gen = self%EFunc(val)
    probgen = exp(-beta*E_Gen)

    if(present(E_T)) then     
      E_T = E_Gen
    endif

  end subroutine
!==========================================================================
  subroutine HarmonicBond_GenerateReverseDist(self, curbox, atompos, probgen)
    use RandomGen, only: grnd, Gaussian
    use ClassyConstants, only: two_pi
    implicit none
    class(HarmonicBond), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp) :: beta
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(out) :: probgen

    logical :: accept
    real(dp) :: E_Bond

    beta = curbox%beta
    call self%DetailedECalc(curbox, atompos, E_Bond, accept)
    probgen = exp(-beta*E_Bond)
  end subroutine
!=============================================================================+
  subroutine HarmonicBond_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    use Units, only: inLenUnit, inEngUnit
    implicit none
    class(HarmonicBond), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

    !Input style
    !harmonic r0 k0

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) self%r0
    self%r0 = self%r0*inLenUnit

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) self%k0
    self%k0 = self%k0*inEngUnit

    if(lineStat /= 0) then
      write(0,*) "Missing input rquired for the harmonic bond style"
    endif
  end subroutine
!=============================================================================+
end module
!=============================================================================+
