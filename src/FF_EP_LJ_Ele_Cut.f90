!================================================================================
! Easy Pair Template for general pair forcefields. The purpose of this class is
! to act as an inheritable object that can allow the user to quickly set up
! any standard pair distance based forcefields by simply creating a child object
! from this class. By overwriting the PairFunction class and parameter settings
! the procedures for all cut-off based forcefields are automatically defined.
!================================================================================
module FF_EasyEP_LJ_Ele_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use ClassyConstants, only: coulombConst

  type, extends(EasyPair_Cut) :: EP_LJ_Ele_Cut
!    real(dp), allocatable :: rMin(:)
!    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq

    real(dp), allocatable :: eps(:)
!    real(dp), parameter :: coul = 1.671009770E5_dp
    real(dp), allocatable :: sig(:)
    real(dp), allocatable :: q(:)
!    real(dp), allocatable :: rMin(:)

    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: qTable(:,:)
!    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    real(dp) :: rLJCut, rLJCutSq
    real(dp) :: rQCut, rQCutSq
    contains
      procedure, pass :: Constructor => Constructor_EP_LJ_Ele_Cut
      procedure, pass :: PairFunction => PairFunction_EP_LJ_Ele_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_LJ_Ele_Cut
      procedure, pass :: Prologue => Prologue_EP_LJ_Ele_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EP_LJ_Ele_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_LJ_Ele_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%eps(1:nAtomTypes), stat = AllocateStat)
    allocate(self%sig(1:nAtomTypes), stat = AllocateStat)
    allocate(self%q(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)

    allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%sigTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%qTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%eps = 4E0_dp
    self%sig = 1E0_dp
    self%q = 0E0_dp
    self%rMin = 0.5E0_dp

    self%epsTable = 4E0_dp
    self%sigTable = 1E0_dp
    self%qTable = 0E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2
    self%rLJCut = 5E0_dp
    self%rLJCutSq = 5E0_dp**2
    self%rQCut = 5E0_dp
    self%rQCutSq = 5E0_dp**2
    IF (AllocateStat /= 0) error STOP "Allocation in the LJ/Cut Pair Style"

  end subroutine
  !=============================================================================+
  function PairFunction_EP_LJ_Ele_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_LJ_Ele_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: ep, sig_sq, q, LJ, Ele, r

    E_Pair = 0E0_dp
    if(rsq < self%rLJCutSq) then
      ep = self % epsTable(atmType1, atmType2)
      sig_sq = self % sigTable(atmType1, atmType2)          
      LJ = (sig_sq/rsq)**3
      LJ = ep * LJ * (LJ-1E0_dp)
      E_Pair = E_Pair + LJ
    endif

    if(rsq < self%rQCutSq) then
      q = self % qTable(atmType1, atmType2)
      if(abs(q) > 1E-15_dp) then
        r = sqrt(rsq)
        Ele = q/r
        E_Pair = E_Pair + Ele
      endif
    endif


  end function
  !=====================================================================
  subroutine ProcessIO_EP_LJ_Ele_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_LJ_Ele_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, q, rCut, rMin
  

    call GetXCommand(line, command, 1, lineStat)

    call CountCommands(line, nPar)
    select case(trim(adjustl(command)))
      case("ljrcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rLJCut = rCut
        self % rLJCutSq = rCut*rCut
        self % rCut = max(self%rLJCut, self%rQCut)
        self % rCutSq = self % rCut * self % rCut

      case("qrcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rQCut = rCut
        self % rQCutSq = rCut*rCut
        self % rCut = max(self%rLJCut, self%rQCut)
        self % rCutSq = self % rCut * self % rCut

      case("charge")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1

        call GetXCommand(line, command, 3, lineStat)
        read(command, *) q
        self%q(type1) = q
        do jType = 1, nAtomTypes
          self%qTable(type1, jType) = q * self%q(jType) * coulombConst
          self%qTable(jType, type1) = q * self%q(jType) * coulombConst
        enddo


      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      select case(nPar)
        case(5)
          read(line, *) type1, ep, sig, q, rMin
          ep = ep * inEngUnit
          sig = sig * inLenUnit
          rMin = rMin * inLenUnit


          self%eps(type1) = ep
          self%sig(type1) = sig
          self%q(type1) = q
          self%rMin(type1) = rMin

          do jType = 1, nAtomTypes
            self%epsTable(type1, jType) = 4E0_dp * sqrt(ep * self%eps(jType))
            self%epsTable(jType, type1) = 4E0_dp * sqrt(ep * self%eps(jType))


            self%sigTable(type1, jType) = (0.5E0_dp * (sig + self%sig(jType)))**2
            self%sigTable(jType, type1) = (0.5E0_dp * (sig + self%sig(jType)))**2

            self%qTable(type1, jType) = q * self%q(jType) * coulombConst
            self%qTable(jType, type1) = q * self%q(jType) * coulombConst

            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case default
          lineStat = -1
          write(0,*) line
          stop "Unknown Command Found"
      end select
    endif

  end subroutine
  !=====================================================================
  subroutine Prologue_EP_LJ_Ele_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(EP_LJ_Ele_Cut), intent(inout) :: self
    integer :: i, j
    
    
    write(nout,*) "Cutoffs (LJ Q Global):", self%rLJCut, self%rQCut, self%rCut
    write(nout,*) "Cutoffs Sq (LJ Q Global):", self%rLJCutSq, self%rQCutSq, self%rCutSq
    write(nout,*)

!    do i = 1, nAtomTypes
!      write(nout, *) (self%epsTable(i,j), j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%sigTable(i,j), j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%rMinTable(i,j), j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%qTable(i,j)/coulombConst, j=1,nAtomTypes)
!    enddo

  end subroutine

  !=====================================================================
end module
!=====================================================================
