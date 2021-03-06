!========================================================
module MCMove_UB_Simple
use CoordinateTypes, only: DisplacementNew, Addition, Deletion
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: UB_Simple
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: inatmps = 1E-30_dp
    real(dp) :: inaccpt = 0E0_dp
    real(dp) :: outatmps = 1E-30_dp
    real(dp) :: outaccpt = 0E0_dp
    real(dp) :: avbmcRad = 1.5E0_dp
    real(dp) :: avbmcVol = 0E0_dp
    type(Addition) :: newPart(1:1)
    type(Deletion) :: oldPart(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => UB_Simp_Constructor
!      procedure, pass :: GeneratePosition => UB_Simp_GeneratePosition
      procedure, pass :: FullMove => UB_Simp_FullMove
      procedure, pass :: SwapIn => UB_Simp_SwapIn
      procedure, pass :: SwapOut => UB_Simp_SwapOut
!      procedure, pass :: Maintenance => UB_Simp_Maintenance
      procedure, pass :: Prologue => UB_Simp_Prologue
      procedure, pass :: Epilogue => UB_Simp_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine UB_Simp_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(UB_Simple), intent(inout) :: self
!    integer :: iType, maxAtoms



!    allocate( self%tempNNei(1) )
!    allocate( self%tempList(200, 1) )
  end subroutine
!========================================================
!  subroutine UB_Simp_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(UB_Simple), intent(in) :: self
!    type(DisplacementNew), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine UB_Simp_FullMove(self, trialBox, accept) 
    use RandomGen, only: grnd
    implicit none
    class(UB_Simple), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    if(grnd() > 0.5E0_dp) then
      call self % SwapIn(trialBox, accept)
    else
      call self % SwapOut(trialBox, accept)
    endif

  end subroutine
!===============================================
  subroutine UB_Simp_SwapIn(self, trialBox, accept) 
    use Box_Utility, only: FindAtom
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(UB_Simple), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, radius
    real(dp) :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    self % inatmps = self % inatmps + 1E0_dp
    accept = .true.

!    integer(kind=atomIntType) :: molType, atmIndx, molIndx
!    real(dp) :: x_new, y_new, z_new
!    integer :: listIndex = -1
    if(trialBox%NMol(1) + 1 > trialBox%NMolMax(1)) then
      accept = .false.
      return
    endif

    nMove = trialBox%NMol(1) + 1
    nMove = trialbox%MolStartIndx(nMove)
    !Choose an atom to serve as the target for the new molecule.
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp)
    call FindAtom(trialbox, rawIndx, nTarget)

    !Choose the position relative to the target atom 
    call Generate_UnitSphere(dx, dy, dz)
    radius = self % avbmcRad * grnd()**(1.0E0_dp/3.0E0_dp)
    dx = radius * dx
    dy = radius * dy
    dz = radius * dz

    self%newPart(1)%molType = trialBox%MolType(nMove)
    self%newPart(1)%molIndx = trialBox%MolIndx(nMove)
    self%newPart(1)%atmIndx = nMove

    self%newPart(1)%x_new = trialBox%atoms(1, nTarget) + dx
    self%newPart(1)%y_new = trialBox%atoms(2, nTarget) + dy
    self%newPart(1)%z_new = trialBox%atoms(3, nTarget) + dz

    call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%newPart(1), &
                                              nCount, self%avbmcRad)

    self%newPart(1)%listIndex = 1

    !Check Constraint
    accept = trialBox % CheckConstraint( self%newPart(1:1) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%newPart(1:1), self%tempList, &
                                              self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif

    Prob = real(trialBox%nAtoms, dp) * self%avbmcVol
    Prob = Prob/(real(nCount, dp) * real(trialBox%nAtoms+1, dp))
!    write(*,*) "Prob In", Prob, E_Diff, trialBox%nAtoms, self%avbmcVol, nCount, trialBox%nAtoms+1

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%newPart(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % inaccpt = self % inaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%newPart(1:1), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine UB_Simp_SwapOut(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Box_Utility, only: FindAtom
    implicit none
    class(UB_Simple), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, nNei, nCount
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE
    real(dp) :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    self % outatmps = self % outatmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp)
    call FindAtom(trialbox, rawIndx, nMove)
    if(trialBox%NMol(1) - 1 < trialBox%NMolMin(1)) then
!      write(*,*) "Bounds Rejection"
      accept = .false.
      return
    endif



    self%oldPart(1)%molType = trialBox%MolType(nMove)
    self%oldPart(1)%molIndx = trialBox%MolIndx(nMove)
!    self%oldPart(1)%atmIndx = nMove

    !Check Constraint
    accept = trialBox % CheckConstraint( self%oldPart(1:1) )
    if(.not. accept) then
!      write(*,*) "Constraint Rejection"
!      write(*,*) "============================================"
      return
    endif

    !Energy Calculation
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
!      write(*,*) "Energy Rejection"
      return
    endif

    nNei = trialBox % NeighList(1) % GetNeighCount (nMove, self%avbmcRad)

    Prob = real(nNei, dp) * real(trialBox%nAtoms, dp)
    Prob = Prob/(real(trialBox%nAtoms-1, dp) * self%avbmcVol)
!    write(*,*) "Prob Out:", Prob, trialBox%nAtoms, self%avbmcVol, nNei, trialBox%nAtoms-1

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%oldPart(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % outaccpt = self % outaccpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(self%oldPart(1)%molIndx)
    endif


  end subroutine
!=========================================================================
  subroutine UB_Simp_Maintenance(self)
    implicit none
    class(UB_Simple), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine UB_Simp_Prologue(self)
    use ParallelVar, only: nout
    use ClassyConstants, only: pi
    implicit none
    class(UB_Simple), intent(inout) :: self

    self%avbmcVol = (4E0_dp/3E0_dp)*pi*self%avbmcRad**3
!    write(*,*) self%avbmcVol

    allocate( self%tempNNei(1) )
    allocate( self%tempList(200, 1) )
  end subroutine
!=========================================================================
  subroutine UB_Simp_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(UB_Simple), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "UB Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "UB Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "UB Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "UB Out Moves Accepted: ", nint(self%outaccpt)
    write(nout,"(1x,A,I15)") "UB Out Moves Attempted: ", nint(self%outatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "UB Out Acceptance Rate: ", accptRate

    write(nout,"(1x,A,I15)") "UB In Moves Accepted: ", nint(self%inaccpt)
    write(nout,"(1x,A,I15)") "UB In Moves Attempted: ", nint(self%inatmps)
    accptRate = 1E2_dp * self%outaccpt/self%outatmps
    write(nout,"(1x,A,F15.8)") "UB In Acceptance Rate: ", accptRate
 

  end subroutine
!========================================================
end module
!========================================================
