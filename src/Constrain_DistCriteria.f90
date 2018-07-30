!====================================================================
!This module contains the Stilinger distance criteria that is used to 
!enforce clustering. 
!====================================================================
module Constrain_DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Displacement, Perturbation
  use CoordinateTypes, only: DisplacementNew, Deletion, Addition
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout

  type, public, extends(constraint) :: DistCriteria
    integer :: neighList = 1
    integer :: molType = 1
    integer :: atomNum = 1
    real(dp) :: rCut, rCutSq
    integer :: boxID = 1

    logical, allocatable :: flipped(:)
    logical, allocatable :: clustMemb(:)
    logical, allocatable :: topoList(:, :)
    logical, allocatable :: newTopoList(:, :)
    integer, allocatable :: newList(:)
    integer, allocatable :: newList2(:)
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => DistCrit_Constructor
      procedure, pass :: CheckInitialConstraint => DistCrit_CheckInitialConstraint
      procedure, pass :: DiffCheck => DistCrit_DiffCheck
      procedure, pass :: ShiftCheck => DistCrit_ShiftCheck
      procedure, pass :: NewCheck => DistCrit_NewCheck
      procedure, pass :: OldCheck => DistCrit_OldCheck
!      procedure, pass :: CheckCluster => DistCrit_CheckCluster
      procedure, pass :: ProcessIO => DistCrit_ProcessIO
      procedure, pass :: Maintenance => DistCrit_Maintenance
      procedure, pass :: Update => DistCrit_Update
      procedure, pass :: Epilogue => DistCrit_Epilogue
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
    allocate(self%topoList(1:nMolMax, 1:nMolMax), stat = AllocateStat)
    allocate(self%newTopoList(1:nMolMax, 1:nMolMax), stat = AllocateStat)
    allocate(self%newlist(1:nMolMax), stat = AllocateStat)
    allocate(self%newlist2(1:nMolMax), stat = AllocateStat)

    IF (AllocateStat /= 0) STOP "Allocation Error in Distance Constraint"
  end subroutine
!=====================================================================
  subroutine DistCrit_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    integer :: totalMol, nNew, nClust, neiIndx
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    integer :: molIndx, molType, nNext, nNext2, iNext
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    self%flipped = .false.
    self%clustMemb = .false.
    self%topoList = .false.

    totalMol = trialBox%NMol(self%molType)
    if(totalMol < 2) then
      return
    endif

     !Build the topology list that will be used through out the simulation
    do iMol = 1, totalMol-1
      molIndx = trialBox % MolGlobalIndx(self%molType, iMol)
      iAtom = trialBox % MolStartIndx(molIndx) + self%atomNum - 1
!      write(*,*) iMol, molIndx

      do jMol = iMol+1, totalMol
        molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
        jAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
        rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
        ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
        rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
        call trialBox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) iMol, jMol, rx, ry, rz, rsq
        if(rsq < self%rCutSq ) then
          self%topoList(iMol, jMol) = .true.
          self%topoList(jMol, iMol) = .true.
        endif
      enddo
    enddo
 
!    do iMol = 1, totalMol
!      write(*,*) (self%topoList(iMol, jMol), jMol=1,totalMol)
!    enddo

    if(all(self%topoList .eqv. .false.) ) then
      accept = .false.
      write(nout,*) "Detailed Cluster Criteria Check Failed!"
      return
    endif

     !Seed the initial cluter check by adding the first particle in the array
     !to the cluster
!    molIndx = trialBox % MolGlobalIndx(self%molType, 1)
    self%clustMemb(1) = .true.
    self%newlist(1) = 1
    nNew = 1
    nClust = 1
    
    do iLimit = 1, totalMol
      nNext = nNew
      nNew = 0
      do iNext = 1, nNext
        iMol = self%newlist(iNext)
        do jMol = 1, totalMol
          if(self%topoList(jMol, iMol) )then
            if(.not. self%clustMemb(jMol) )then
              self%clustMemb(jMol) = .true.
              nNew = nNew + 1
              nClust = nClust + 1
              self%newlist2(nNew) = jMol
            endif
          endif
        enddo
      enddo

       !If every molecule has been added then no further calculations are
       !needed.
      if(nClust >= totalMol) then
        exit
      endif 
       ! If no new molecules were added, the algorithm has hit a dead end
       ! and the cluster is broken.  
      if(nNew <= 0) then
        exit
      endif

      self%newlist(1:nNew) = self%newlist2(1:nNew)
    enddo

     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
    if( (nNew <= 0) .or. (nClust < totalMol) ) then
      accept = .false.
      write(nout,*) "Detailed Cluster Criteria Check Failed!"
      return
    endif

    write(nout,*) "Detailed Cluster Criteria Check Succeeded!"
    accept = .true.

  end subroutine
!=============================================================
  subroutine DistCrit_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
!    type(Displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    integer :: iDisp
    integer :: totalMol, nNew, nClust, neiIndx, startMol
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    integer :: molIndx, molType, nNext, nNext2, iNext
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    totalMol = trialBox%NMol(self%molType)

    !This section creates the topology list of the new state using information
    !based on the what kind of perturbation was performed.
    select type(disp)
      class is(Displacement)
        stop "Diff Done"
        !Called when a particle is either moved or replaced with another particle
        if(disp(1)%newAtom .and. disp(1)%oldAtom) then
          if(disp(1)%molType == disp(1)%oldMolType) then
            call self % ShiftCheck(trialBox, disp, accept)
          elseif(self%molType == disp(1)%molType) then
            call self % NewCheck(trialBox, disp, accept)
          elseif(self%molType == disp(1)%oldMolType) then
            call self % OldCheck(trialBox, disp, accept)
          endif
          return
        endif

        !Called when a particle is added to a system
        if(disp(1)%newAtom) then
          call self % NewCheck(trialBox, disp, accept)
          return
        endif

        !Called when a particle is removed from a system.
        if(disp(1)%oldAtom) then
          call self % OldCheck(trialBox, disp, accept)
          return
        endif

      class is(DisplacementNew)
        self%newTopoList = self%topoList 
        accept = .true.
        do iDisp = 1, size(disp)
          if( disp(iDisp)%molType == self%molType ) then
            molIndx = disp(iDisp)%molIndx
            iAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
            if( disp(iDisp)%atmIndx == iAtom) then
              accept = .false.
              iMol = disp(iDisp)%molIndx
              do jMol = 1, totalMol
                if(iMol /= jMol) then
                  self%newTopoList(jMol, iMol) = .false.
                  self%newTopoList(iMol, jMol) = .false.
                endif
              enddo
              do jMol = 1, totalMol
                if(iMol /= jMol) then
                  molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
                  jAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
                  rx = disp(iDisp)%x_new - trialBox%atoms(1, jAtom)
                  ry = disp(iDisp)%y_new - trialBox%atoms(2, jAtom)
                  rz = disp(iDisp)%z_new - trialBox%atoms(3, jAtom)
                  call trialBox%Boundary(rx, ry, rz)
                  rsq = rx*rx + ry*ry + rz*rz
!                  write(*,*) jMol, rx, ry, rz, rsq
                  if(rsq < self%rCutSq ) then
                    self%newTopoList(jMol, iMol) = .true.
                    self%newTopoList(iMol, jMol) = .true.
                  endif
                endif
              enddo
            endif
          endif
        enddo

        if(accept) then
          return
        endif

        self%clustMemb = .false.
        self%clustMemb(1) = .true.
        self%newlist(1) = 1
        nNew = 1
        nClust = 1

      class is(Addition)
        self%newTopoList = self%topoList 
        accept = .true.
        do iDisp = 1, size(disp)
          if( disp(iDisp)%molType == self%molType ) then
            molIndx = disp(iDisp)%molIndx
            iAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
            if( disp(iDisp)%atmIndx == iAtom) then
              accept = .false.
              iMol = disp(iDisp)%molIndx
              do jMol = 1, totalMol
                if(iMol /= jMol) then
                  self%newTopoList(jMol, iMol) = .false.
                  self%newTopoList(iMol, jMol) = .false.
                endif
              enddo
              do jMol = 1, totalMol
                if(iMol /= jMol) then
                  molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
                  jAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
                  rx = disp(iDisp)%x_new - trialBox%atoms(1, jAtom)
                  ry = disp(iDisp)%y_new - trialBox%atoms(2, jAtom)
                  rz = disp(iDisp)%z_new - trialBox%atoms(3, jAtom)
                  call trialBox%Boundary(rx, ry, rz)
                  rsq = rx*rx + ry*ry + rz*rz
!                  write(*,*) jMol, rx, ry, rz, rsq
                  if(rsq < self%rCutSq ) then
                    self%newTopoList(jMol, iMol) = .true.
                    self%newTopoList(iMol, jMol) = .true.
                  endif
                endif
              enddo
            endif
          endif
        enddo

        if(accept) then
          return
        endif


      class is(Deletion)
        !molType, atmIndx, molIndx
        self%newTopoList = self%topoList 
        accept = .true.
        do iDisp = 1, size(disp)
          if( disp(iDisp)%molType == self%molType ) then
            if( disp(iDisp)%atmIndx == self%atomNum ) then
              accept = .false.
              iMol = disp(iDisp)%molIndx
              do jMol = 1, totalMol
                if(self%newTopoList(jMol, iMol)) then
                  self%newTopoList(jMol, iMol) = .false.
                  self%newTopoList(iMol, jMol) = .false.
                endif
              enddo
            endif
          endif
        enddo


        if(accept) then
          return
        endif

        startMol = 1
        do iMol = 1, totalMol
          do jMol = iMol+1, totalMol
            if(self%newTopoList(jMol, iMol)) then
              startMol = iMol
              exit
            endif

          enddo
          if(self%newTopoList(jMol, iMol)) then
            exit
          endif
        enddo
        self%clustMemb = .false.
        self%clustMemb(startMol) = .true.
        self%newlist(1) = startMol
        nNew = 1
        nClust = 1

      class default
        stop "Distance criteria is not compatiable with this perturbation type."
    end select
    
!    write(*,*)
!    do iMol = 1, totalMol
!      write(*,*) (self%newTopoList(jMol, iMol), jMol=1,totalMol)
!    enddo

    !Using the newly constructed topology list, check to see if the new cluster satisfies
    !the cluster criteria. 
    do iLimit = 1, totalMol
      nNext = nNew
      nNew = 0
      do iNext = 1, nNext
        iMol = self%newlist(iNext)
        do jMol = 1, totalMol
          if(self%newTopoList(jMol, iMol) )then
            if(.not. self%clustMemb(jMol) )then
              self%clustMemb(jMol) = .true.
              nNew = nNew + 1
              nClust = nClust + 1
              self%newlist2(nNew) = jMol
            endif
          endif
        enddo
      enddo

       !If every molecule has been added then no further calculations are
       !needed.
      if(nClust >= totalMol) then
        exit
      endif 
       ! If no new molecules were added, the algorithm has hit a dead end
       ! and the cluster is broken.  
      if(nNew <= 0) then
        exit
      endif

      self%newlist(1:nNew) = self%newlist2(1:nNew)
    enddo

     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
    if( (nNew <= 0) .or. (nClust < totalMol) ) then
      accept = .false.
!      write(*,*) .false.
      return
    endif
    accept = .true.
!    write(*,*) .true.

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
    iMol = trialBox % MolSubIndx(startIndx)

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
      write(nout, *) iMol, iAtom
      stop
    endif


     !Seed the initial cluter check by adding the first particle in the array
     !to the cluster
    self%clustMemb(iMol) = .true.
    self%flipped(iMol) = .true.
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
    integer :: i, jMol, jAtom, dispIndx, molIndx, totalMol
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    do i = 1, size(disp)
      if(disp(i)%MolType == self%molType) then
        if(disp(i)%AtmIndx == trialBox%MolStartIndx(disp(i)%molIndx)) then
          accept = .false.
          dispIndx = i
          exit
        endif
      endif
    enddo
 
    if(accept) then
      return
    endif

    totalMol = trialBox%NMol(self%molType)

    do jMol = 1, totalMol  
      molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
      jAtom = trialBox % MolStartIndx(molIndx)
      rx = disp(dispIndx)%x_new - trialBox%atoms(1, jAtom)
      ry = disp(dispIndx)%y_new - trialBox%atoms(2, jAtom)
      rz = disp(dispIndx)%z_new - trialBox%atoms(3, jAtom)
      call trialBox%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < self%rCutSq) then
        accept = .true.
        return
      endif
    enddo



  end subroutine
!=====================================================================
  subroutine DistCrit_OldCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    integer :: i, startIndx, molIndx, jMolIndx
    integer :: totalMol, nNew, nClust, neiIndx
    integer :: iMol,jMol,jNei, iAtom, jAtom, iLimit
    integer :: molType, dispIndx, nMol
    integer :: neiList(1:600), nOldNei, nNewNei

    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    totalMol = trialBox%NMol(self%molType)
    if(totalMol-1 < 2) then
      return
    endif

    do i = 1, size(disp)
      if(disp(i)%oldMolType == self%molType) then
        if(disp(i)%oldAtmIndx == trialBox%MolStartIndx(disp(i)%molIndx)) then
          accept = .false.
          dispIndx = i
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
    nMol = trialBox % MolSubIndx(startIndx)

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
    nClust = 1
    nNew = 0
    nNewNei = 1
    self%clustMemb( neiList(1) ) = .true.

    do iLimit = 1, totalMol-1
      nNew = 0
      do iMol = 1, totalMol
        if( (.not. self%clustMemb(iMol))  .and. (iMol /= nMol ) ) then
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
            if( (.not. self%clustMemb(jMol)) .and. (jMol /= nMol)  )then
              rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
              ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
              rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
              call trialBox%Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
!              write(*,*) rx, ry, rz, rsq
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

!      write(*,*) iLimit, nClust, nNew, nNewNei
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
      if(nClust >= totalMol-1) then
        exit
      endif
    enddo

!    write(*,*) "Blah3"
     ! If no new particles were added or the limit has been hit without finding all the molecules
     ! then a disconnect in the cluster network was created and the criteria has not been satisfied. 
    if( (nNew <= 0) .or. (nClust < totalMol-1) ) then
!      write(*,*) "Fail 2", nNew, nClust
      accept = .false.
      return
    endif

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
!====================================================================
  subroutine DistCrit_Maintenance(self)
    implicit none
    class(DistCriteria), intent(inout) :: self

  end subroutine
!====================================================================
  subroutine DistCrit_Epilogue(self)
    implicit none
    class(DistCriteria), intent(inout) :: self
    logical :: accept

    call self % CheckInitialConstraint(self%parent, accept)

  end subroutine
!=============================================================
  subroutine DistCrit_Update(self)
    implicit none
    class(DistCriteria), intent(inout) :: self

!    write(*,*) "Update"
    self%topoList = self%newTopoList

  end subroutine
!=====================================================================
end module
!=====================================================================
