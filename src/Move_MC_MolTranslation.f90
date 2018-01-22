!========================================================
module MCMove_MolTranslation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: MolTranslate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp

    real(dp), allocatable :: molAccpt(:)
    real(dp), allocatable :: molAtmps(:)
    real(dp), allocatable :: max_dist(:)
    type(Displacement), allocatable :: disp(:)

!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => MolTrans_Constructor
      procedure, pass :: GeneratePosition => MolTrans_GeneratePosition
      procedure, pass :: FullMove => MolTrans_FullMove
      procedure, pass :: Maintenance => MolTrans_Maintenance
!      procedure, pass :: Prologue => MolTrans_Prologue
      procedure, pass :: Epilogue => MolTrans_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine MolTrans_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(MolTranslate), intent(inout) :: self
    integer :: iType, maxAtoms

    maxAtoms = -1
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo

    allocate( self%disp(1:maxAtoms) )
    allocate( self%max_dist(1:nMolTypes) )

    allocate( self%molAtmps(1:nMolTypes) )
    allocate( self%molAccpt(1:nMolTypes) )

    allocate( self%tempNNei(1:maxAtoms) )
    allocate( self%tempList(1000, 1:maxAtoms) )
  end subroutine
!========================================================
  subroutine MolTrans_GeneratePosition(self, disp)
    use RandomGen, only: grnd
    implicit none
    class(MolTranslate), intent(in) :: self
    type(Displacement), intent(inout) :: disp
    real(dp) :: dx, dy, dz
!      dx = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_dist * (2E0_dp*grnd() - 1E0_dp)
  end subroutine
!===============================================
  subroutine MolTrans_FullMove(self, trialBox, accept) 
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use Common_MolInfo, only: nMolTypes, MolData
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Box_Utility, only: FindMolecule
    implicit none
    class(MolTranslate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: nMove, nType, nAtoms, nFirstAtom, rawIndx, iConstrain
    integer :: iAtom, curAtom
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE


    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nMove)
!    write(*,*) nMove, rawIndx
    nType = trialBox%molType(nMove)
    nAtoms = MolData(nType)%nAtoms
    nFirstAtom = trialBox % MolStartIndx(nMove)
    self % molAtmps(nType) = self % molAtmps(nType) + 1E0_dp

    dx = self % max_dist(nType) * (2E0_dp * grnd() - 1E0_dp)
    dy = self % max_dist(nType) * (2E0_dp * grnd() - 1E0_dp)
    dz = self % max_dist(nType) * (2E0_dp * grnd() - 1E0_dp)

    do iAtom = 1, nAtoms
      curAtom = iAtom + nFirstAtom - 1
      self%disp(iAtom)%newatom = .true.
      self%disp(iAtom)%molType = nType
      self%disp(iAtom)%molIndx = nMove
      self%disp(iAtom)%atmIndx = curAtom

      self%disp(iAtom)%oldatom = .true.
      self%disp(iAtom)%oldMolType = nType  
      self%disp(iAtom)%oldMolIndx = nMove
      self%disp(iAtom)%oldAtmIndx = curAtom
      
      self%disp(iAtom)%x_new = trialBox%atoms(1, curAtom) + dx
      self%disp(iAtom)%y_new = trialBox%atoms(2, curAtom) + dy
      self%disp(iAtom)%z_new = trialBox%atoms(3, curAtom) + dz
    enddo

    !If the particle moved a large distance get a temporary neighborlist
    if(any([dx,dy,dz] > neighSkin)) then
      do iAtom = 1, nAtoms
        call trialBox % NeighList(1) % GetNewList(1, self%tempList, self%tempNNei, self%disp(iAtom))
        self%disp(iAtom)%newlist = .true.
      enddo
    else
      do iAtom = 1, nAtoms
        self%disp(iAtom)%newlist = .false.
        self%disp(iAtom)%listIndex = nMove
      enddo
    endif

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialbox % EFunc % Method % DiffECalc(trialBox, self%disp(1:nAtoms), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, 1E0_dp, self%disp(1:nAtoms))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      self % molAccpt(nType) = self % molAccpt(nType) + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine MolTrans_Maintenance(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(MolTranslate), intent(inout) :: self
    real(dp), parameter :: limit = 3.0E0_dp
    
    integer :: iType
      
    if(self%atmps .lt. 0.5E0_dp) then
      return
    endif

    do iType = 1, nMolTypes
      if(self%molAccpt(iType)/self%molAtmps(iType) .gt. 0.5E0_dp) then
        if(self%max_dist(iType)*1.01E0_dp .lt. limit) then
          self%max_dist(iType) = self%max_dist(iType)  * 1.01E0_dp
        else 
          self%max_dist(iType) = limit       
        endif
      else
        self%max_dist(iType) = self%max_dist(iType) * 0.99E0_dp
      endif
    enddo

  end subroutine
!=========================================================================
  subroutine MolTrans_Epilogue(self)
    use Common_MolInfo, only: nMolTypes
    use ParallelVar, only: nout
    implicit none
    class(MolTranslate), intent(inout) :: self
    integer :: iType
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Molecule Translation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Molecule Translation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Molecule Translation Acceptance Rate: ", accptRate
 
    write(nout,"(1x,A)") "By Type: "
    do iType = 1, nMolTypes
      write(nout, *) "Type ", iType,":", 1E2_dp * self%molAccpt(iType)/self%molAtmps(iType)
    enddo


  end subroutine
!========================================================
end module
!========================================================
