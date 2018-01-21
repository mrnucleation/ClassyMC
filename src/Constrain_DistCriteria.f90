!====================================================================
!This module contains the Stilinger distance criteria that is used to 
!enforce clustering. 
!====================================================================
module Constrain_DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Displacement
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout

  type, public, extends(constraint) :: DistCriteria
    integer :: neighList = 1
    integer :: molType = 1
    real(dp) :: rCut, rCutSq
    integer :: boxID = 1

    logical, allocatable :: flipped(:)
    logical, allocatable :: clustMemb(:)
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => DistCrit_Constructor
      procedure, pass :: CheckInitialConstraint => DistCrit_CheckInitialConstraint
!      procedure, pass :: DiffCheck
      procedure, pass :: ShiftCheck => DistCrit_ShiftCheck
      procedure, pass :: NewCheck => DistCrit_NewCheck
      procedure, pass :: OldCheck => DistCrit_OldCheck
      procedure, pass :: ProcessIO => DistCrit_ProcessIO
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine DistCrit_Constructor(self, boxID)
    use BoxData, only: BoxArray
    implicit none
    class(DistCriteria), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: AllocateStat
    integer :: nMolMax

    self%boxID = boxID
    self%parent => BoxArray(boxID) % box 

    nMolMax = self % parent % NMolMax(self%molType)

    allocate(self%flipped(1:nMolMax), stat = AllocateStat)
    allocate(self%clustMemb(1:nMolMax), stat = AllocateStat)

    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"
  end subroutine
!=====================================================================
  subroutine DistCrit_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    integer :: totalMol, nNew, nClust, neiIndx
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    integer :: molIndx, molType
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    self%flipped = .false.
    self%clustMemb = .false.

    totalMol = trialBox%NMol(self%molType)
    if(totalMol < 2) then
      return
    endif


     !Seed the initial cluter check by adding the first particle in the array
     !to the cluster
    iMol = trialBox % MolGlobalIndx(self%molType, 1)
    molIndx = trialBox % MolGlobalIndx(self%molType, iMol)
    self%clustMemb(molIndx) = .true.

    nClust = 1
    do iLimit = 1, totalMol
      nNew = 0
      do iMol = 1, totalMol
        if(.not. self%clustMemb(iMol)) then
          cycle
        endif
        molIndx = trialBox % MolGlobalIndx(self%molType, iMol)
        iAtom = trialBox % MolStartIndx(molIndx)

         !If the member flag is true, but the flipped flag is false
         !the neighbors of this molecule have not been checked.
        if( self%clustMemb(iMol) .neqv. self%flipped(iMol)) then
          do jMol = 1, totalMol
            if(.not. self%clustMemb(jMol) )then
              molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
              jAtom = trialBox % MolStartIndx(molIndx)
              rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
              ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
              rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
              rsq = rx*rx + ry*ry + rz*rz
              if(rsq < self%rCutSq) then
                self%clustMemb(jMol) = .true.
                nNew = nNew + 1
                nClust = nClust + 1
              endif
            endif
          enddo
          self % flipped(iMol) = .true.
        endif
      enddo

       ! If no new molecules were added, the algorithm has hit a dead end
       ! and the cluster is broken.  
       if(nNew <= 0) then
        exit
      endif

       !If every molecule has been added then no further calculations are
       !needed.
      if(nClust >= totalMol) then
        exit
      endif
    enddo

     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
    if( (nNew <= 0) .or. (nClust < totalMol) ) then
      accept = .false.
      return
    endif

    accept = .true.

  end subroutine
!=====================================================================
  subroutine DistCrit_ShiftCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    integer :: i, startIndx, molIndx, jMolIndx
    integer :: totalMol, nNew, nClust, neiIndx
    integer :: iMol,jMol,jNei, iAtom, jAtom, iLimit
    integer :: molType, dispIndx
    integer :: neiList(1:600), nOldNei, nNewNei

    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    totalMol = trialBox%NMol(self%molType)
    if(totalMol < 2) then
      return
    endif

    do i = 1, size(disp)
      if(disp(i)%molType == self%molType) then
        if(disp(i)%atmIndx == trialBox%MolStartIndx(disp(i)%molIndx)) then
          accept = .false.
          dispIndx = i
          startIndx = disp(i)%molIndx
          exit
        endif
      endif
    enddo
 
    if(accept) then
      return
    endif

    self%flipped = .false.
    self%clustMemb = .false.
    iAtom = trialBox % MolStartIndx(startIndx)

    nOldNei = 0 
    do jNei = 1, trialBox%NeighList(self%neighlist)%nNeigh(iAtom)
      jAtom = trialBox%NeighList(self%neighlist)%list(jNei, iAtom)
      rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
      ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
      rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
      call trialBox%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < self%rCutSq) then
        jMol = trialBox%MolIndx(jAtom)
        nOldNei = nOldNei + 1
        neiList(nOldNei) = jMol
      endif
    enddo
  
    if(nOldNei <= 0) then
      write(nout, *) "WARNING! Catestrophic error in distance criteria detected!"
      write(nout, *) "No neighbors were found for the old position for a given"
      write(nout, *) "shift move. Cluster criteria was not properly maintained!"
    endif


     !Seed the initial cluter check by adding the first particle in the array
     !to the cluster
    self%clustMemb(startIndx) = .true.
    self%flipped(startIndx) = .true.
    nClust = 1
    nNew = 0
    nNewNei = 0

    do jNei = 1, trialBox%NeighList(self%neighlist)%nNeigh(iAtom)
      jAtom = trialBox%NeighList(self%neighlist)%list(jNei, iAtom)
      jMol = trialBox%MolSubIndx(jAtom)
      rx = disp(dispIndx)%x_new - trialBox%atoms(1, jAtom)
      ry = disp(dispIndx)%y_new - trialBox%atoms(2, jAtom)
      rz = disp(dispIndx)%z_new - trialBox%atoms(3, jAtom)
      call trialBox%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < self%rCutSq) then
        self%clustMemb(jMol) = .true.
        nNew = nNew + 1
        nClust = nClust + 1
        if( any(jMol == neiList(1:nOldNei)) ) then
          nNewNei = nNewNei + 1
        endif
      endif
    enddo

    if( (nNew <= 0) ) then
      accept = .false.
      return
    elseif(nNewNei == nOldNei) then
      accept = .true.
      return
    endif


    do iLimit = 1, totalMol-1
      nNew = 0
      do iMol = 1, totalMol
        if(.not. self%clustMemb(iMol)) then
          cycle
        endif
        molIndx = trialBox % MolGlobalIndx(self%molType, iMol)
        iAtom = trialBox % MolStartIndx(molIndx)

         !If the member flag is true, but the flipped flag is false
         !the neighbors of this molecule have not been checked.
        if( self%clustMemb(iMol) .neqv. self%flipped(iMol)) then
          do jNei = 1, trialBox%NeighList(self%neighlist)%nNeigh(iAtom)
            jAtom = trialBox%NeighList(self%neighlist)%list(jNei, iAtom)
            jMol = trialBox%MolSubIndx(jAtom)
            if(.not. self%clustMemb(jMol) )then
              molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
              jAtom = trialBox % MolStartIndx(molIndx)
              rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
              ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
              rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
              if(rsq < self%rCutSq) then
                self%clustMemb(jMol) = .true.
                nNew = nNew + 1
                nClust = nClust + 1
                if( any(jMol == neiList(1:nOldNei)) ) then
                  nNewNei = nNewNei + 1
                endif
              endif
            endif
          enddo
          self % flipped(iMol) = .true.
        endif
      enddo

       ! If no new molecules were added, the algorithm has hit a dead end
       ! and the cluster is broken.  
      if(nNew <= 0) then
        exit
      endif

      if(nNewNei == nOldNei) then
        exit
      endif

       !If every molecule has been added then no further calculations are
       !needed.
      if(nClust >= totalMol) then
        exit
      endif
    enddo

     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
    if( (nNew <= 0) .or. (nClust < totalMol) ) then
!      write(*,*) "Fail 2", nNew, nClust
      accept = .false.
      return
    endif

    accept = .true.


  end subroutine
!=====================================================================
  subroutine DistCrit_NewCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=====================================================================
  subroutine DistCrit_OldCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine DistCrit_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ParallelVar, only: nout
    implicit none
    class(DistCriteria), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%molType = intVal
    call GetXCommand(line, command, 3, lineStat)
    read(command, *) realVal
    self%rCut = realVal
    self%rCutSq = realVal*realVal
    write(nout, *) "Distance Criteria:",self%rCut
    write(nout, *) "Distance Criteria (SQ):",self%rCutSq
    call GetXCommand(line, command, 4, lineStat)
    read(command, *) intVal
    self%neighList = intVal

  end subroutine
!=====================================================================
end module
!=====================================================================
