!===================================================================
! Easy Pair Template for LJ pair forcefields. 
!===================================================================
module FF_EasyEP_Pedone_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_Pedone_Cut
!    real(dp), allocatable :: rMin(:)
!    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq


    logical :: implicitSolvent = .false.
    real(dp), allocatable :: repul_tab(:,:), rEqTable(:,:)
    real(dp), allocatable :: qTable(:,:), alpha_Tab(:,:), D_Tab(:,:)
    real(dp), allocatable :: bornRad(:)
    real(dp), allocatable :: q(:)
    contains
      procedure, pass :: Constructor => Constructor_EP_Pedone_Cut
      procedure, pass :: SolventFunction => EP_Pedone_SolventFunction
      procedure, pass :: PairFunction => PairFunction_EP_Pedone_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_Pedone_Cut
!      procedure, pass :: TailCorrection => TailCorrection_EP_Pedone_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EP_Pedone_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Pedone_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    allocate(self%repul_tab(1:nAtomTypes,1:natomtypes), stat = AllocateStat)
    allocate(self%rEqTable(1:nAtomTypes,1:natomtypes), stat = AllocateStat)
    allocate(self%qTable(1:natomtypes,1:natomtypes), stat = AllocateStat)
    allocate(self%alpha_tab(1:nAtomTypes, 1:natomtypes), stat = AllocateStat)
    allocate(self%D_tab(1:nAtomTypes, 1:natomtypes), stat = AllocateStat)

    allocate(self%q(1:nAtomTypes), stat = AllocateStat)
    self%repul_tab = 0E0_dp
    self%rEqTable = 0E0_dp
    self%qTable = 0E0_dp
    self%alpha_tab = 0E0_dp
    self%D_tab = 0E0_dp

    self%q(1) = 1.2E0_dp

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

  end subroutine
!=======================================================================
  pure function EP_Pedone_SolventFunction(self,r, q_ij, born1, born2) result(f)
    class(EP_Pedone_Cut), intent(in) :: self
    real(dp), intent(in) :: r, q_ij, born1, born2
    real(dp) :: f

    f = sqrt(r*r + born1*born2*exp(-r*r/(4E0_dp*born1*born2) ) ) 
    f = -0.5E0_dp*(1E0_dp-1E0_dp/dieletric)*q_ij / f
  end function
  !=============================================================================+
  function PairFunction_EP_Pedone_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_Pedone_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: q_ij, alpha, delta, repul_C, r_eq
    real(dp) :: r
    real(dp) :: E_Ele, E_LJ, E_Morse, E_Solvent
    real(dp) :: born1, born2

    r_eq = self%rEqTable(atmType2, atmType1)
    q_ij = self%qTable(atmType2, atmType1)
    alpha = self%alpha_Tab(atmType2, atmType1)
    delta = self%D_Tab(atmType2, atmType1)
    repul_C = self%repul_tab(atmType2, atmType1)

    E_LJ = (1E0_dp/rsq)**6
    E_LJ = repul_C * E_LJ
 
    r = sqrt(rsq)
    E_Ele = q_ij/r
    E_Solvent = 0E0_dp
    if(self%implicitSolvent) then
      born1 = self%bornRad(atmType1)
      born2 = self%bornRad(atmType2)
      E_Solvent = self%solventFunction(r, q_ij, born1, born2)
    endif

    E_Morse = 1E0_dp - exp(-alpha*(r-r_eq))
    E_Morse = delta*(E_Morse*E_Morse - 1E0_dp)

    E_Pair = E_LJ + E_Ele + E_Morse + E_Solvent


  end function
  !=====================================================================
  subroutine ProcessIO_EP_Pedone_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use ClassyConstants, only: coulombConst
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_Pedone_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    logical :: logicVal
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: q, repC, delta, alpha, rCut, rMin, rEq
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut

      case("implicitsolvent")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) self%implicitSolvent

      case("born")
        if(.not. allocated(self%bornRad)) then
          allocate(self%bornRad(1:nAtomTypes))
        endif
        call CountCommands(line, nPar)
        
        do jType = 1, nAtomTypes
          call GetXCommand(line, command, jType+1, lineStat)
          read(command, *) self%bornRad(jType)
        enddo


      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      call CountCommands(line, nPar)
      select case(nPar)
        case(7)
          read(line, *) type1, repC, rEq, alpha, delta, q, rMin
          self%q(type1) = q
          repC = repC * inEngUnit/inLenUnit**12
          rEq = rEq * inLenUnit
          alpha = alpha / inLenUnit**2
          delta = delta * inEngUnit
          rMin = rMin * inLenUnit

          self%repul_tab(1, type1) = repC
          self%repul_tab(type1, 1) = repC

          self%rEqTable(1, type1) = rEq
          self%rEqTable(type1, 1) = rEq

          self%alpha_tab(1, type1) = alpha
          self%alpha_tab(type1, 1) = alpha

          self%D_tab(1, type1) = delta
          self%D_tab(type1, 1) = delta


          do jType = 1, nAtomTypes
            self%qTable(type1, jType) = q * self%q(jType) * coulombConst
            self%qTable(jType, type1) = q * self%q(jType) * coulombConst

            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case default
          lineStat = -1
      end select
    endif


  end subroutine
  !=====================================================================
end module
!=====================================================================
