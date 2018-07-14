!========================================================
module MCMove_AVBMC_Simple
use CoordinateTypes, only: DisplacementNew, Addition
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: AVBMC_Simple
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: avbmcRad = 1.5E0_dp
    real(dp) :: avbmcVol = 0E0_dp
    type(Addition) :: newPart(1:1)
    type(Deletion) :: oldPart(1:1)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => AVBMC_Simp_Constructor
!      procedure, pass :: GeneratePosition => AVBMC_Simp_GeneratePosition
      procedure, pass :: FullMove => AVBMC_Simp_FullMove
      procedure, pass :: SwapIn => AVBMC_SwapIn
      procedure, pass :: SwapOut => AVBMC_SwapOut
!      procedure, pass :: Maintenance => AVBMC_Simp_Maintenance
      procedure, pass :: Prologue => AVBMC_Simp_Prologue
      procedure, pass :: Epilogue => AVBMC_Simp_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine AVBMC_Simp_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(AVBMC_Simple), intent(inout) :: self
!    integer :: iType, maxAtoms



    allocate( self%tempNNei(1) )
    allocate( self%tempList(200, 1) )
  end subroutine
!========================================================
!  subroutine AVBMC_Simp_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(AVBMC_Simple), intent(in) :: self
!    type(DisplacementNew), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine AVBMC_Simp_FullMove(self, trialBox, accept) 
    implicit none
    class(AVBMC_Simple), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    if(grnd() > 0.5E0_dp) then
      call self % SwapIn(trialBox, accept)
    else
      call self % SwapOut(trialBox, accept)
    endif

  end subroutine
!===============================================
  subroutine AVBMC_Simp_SwapIn(self, trialBox, accept) 
    use Box_Utility, only: FindAtom
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd, Generate_UnitSphere
    implicit none
    class(AVBMC_Simple), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nTarget, nType, rawIndx, iConstrain
    integer :: CalcIndex, nMove
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE, radius
    real(dp) :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
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

    !If the particle moved a large distance get a temporary neighborlist
    call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%newPart(1) &
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
    Prob = Prob/(real(trialBox%nCount, dp) * real(trialBox%nAtoms+1, dp))

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%newPart(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%newPart(1:1), self%tempList, self%tempNNei)
    endif

  end subroutine
!===============================================
  subroutine AVBMC_Simp_SwapOut(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Box_Utility, only: FindAtom
    implicit none
    class(AVBMC_Simple), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, nNei
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE
    real(dp), parameter :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp)
    call FindAtom(trialbox, rawIndx, nMove)

    self%oldPart(1)%molType = trialBox%MolType(nMove)
    self%oldPart(1)%molIndx = trialBox%MolIndx(nMove)
    self%oldPart(1)%atmIndx = nMove

    !Check Constraint
    accept = trialBox % CheckConstraint( self%oldPart(1:1) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
!    call trialbox% EFunc % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%oldPart(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif

    nNei = trialBox % NeighList(1) % GetNeighCount (nMove, self%avbmcRad)


    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%oldPart(1:1))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % DeleteMol(self%oldPart(1:1)%molIndx)
    endif


    endif

  end subroutine

!=========================================================================
  subroutine AVBMC_Simp_Maintenance(self)
    implicit none
    class(AVBMC_Simple), intent(inout) :: self
 

  end subroutine
!=========================================================================
  subroutine AVBMC_Simp_Prologue(self)
    use ParallelVar, only: nout
    use Constants, only: pi
    implicit none
    class(AVBMC_Simple), intent(inout) :: self

    self%avbmcVol = (4E0_dp/3E0_dp)*pi*self%avbmcRad**3

  end subroutine
!=========================================================================
  subroutine AVBMC_Simp_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(AVBMC_Simple), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "AVBMC Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "AVBMC Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "AVBMC Acceptance Rate: ", accptRate
!    write(nout,"(1x,A,F15.8)") "Final Maximum Displacement: ", self%max_dist
 

  end subroutine
!========================================================
end module
!========================================================
