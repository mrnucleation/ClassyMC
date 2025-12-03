!====================================================================
!This module contains the Stilinger distance criteria that is used to 
!enforce clustering. 
!====================================================================
module Constrain_DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use CoordinateTypes, only: Displacement, Deletion, Addition
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
    class(SimBox), intent(inout) :: trialBox
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

      do jMol = iMol+1, totalMol
        molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
        jAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
        rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
        ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
        rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
        call trialBox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq ) then
          self%topoList(iMol, jMol) = .true.
          self%topoList(jMol, iMol) = .true.
        endif
      enddo
    enddo


!    write(*,*) "REAL!"
!    do iMol = 1, totalMol
!      write(*,*) (self%topoList(jMol, iMol), jMol=1,totalMol)
!    enddo
!    write(*,*)

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
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    logical :: leave
    integer :: iDisp
    integer :: totalMol, nNew, nClust, neiIndx, startMol, topMol
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    integer :: molIndx, molType, nNext, nNext2, iNext
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    totalMol = trialBox%NMol(self%molType)

    !This section creates the topology list of the new state using information
    !based on the what kind of perturbation was performed.
    select type(disp)
       !----------------------------------------------------------------------------
      class is(Displacement)
        self%newTopoList = self%topoList 
        accept = .true.
!        write(*,*) "Disp"
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
!          write(*,*) "No move"
          return
        endif

        self%clustMemb = .false.
        self%clustMemb(1) = .true.
        self%newlist(1) = 1
        nNew = 1
        nClust = 1

       !----------------------------------------------------------------------------
      class is(Addition)
        self%newTopoList = self%topoList 
        accept = .true.
!        write(*,*) "Add"
        do iDisp = 1, size(disp)
          if( disp(iDisp)%molType == self%molType ) then
            molIndx = disp(iDisp)%molIndx
            iAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
            if( disp(iDisp)%atmIndx == iAtom) then
              accept = .false.
              iMol = disp(iDisp)%molIndx
!              do jMol = 1, totalMol
!                if(iMol /= jMol) then
!                  self%newTopoList(jMol, iMol) = .false.
!                  self%newTopoList(iMol, jMol) = .false.
!                endif
!              enddo
              do jMol = 1, totalMol
                if(iMol /= jMol) then
!                  write(*,*) iMol, jMol
                  molIndx = trialBox % MolGlobalIndx(self%molType, jMol)
                  jAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
                  rx = disp(iDisp)%x_new - trialBox%atoms(1, jAtom)
                  ry = disp(iDisp)%y_new - trialBox%atoms(2, jAtom)
                  rz = disp(iDisp)%z_new - trialBox%atoms(3, jAtom)
                  call trialBox%Boundary(rx, ry, rz)
                  rsq = rx*rx + ry*ry + rz*rz
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
        totalMol = totalMol + 1
       !----------------------------------------------------------------------------
      class is(Deletion)
        !molType, atmIndx, molIndx
!        write(*,*) "---------------------------------------"
!        write(*,*) "Del"
        self%newTopoList = self%topoList 
        accept = .true.
        do iDisp = 1, size(disp)
          if( disp(iDisp)%molType == self%molType ) then
            molIndx = disp(iDisp)%molIndx
            iAtom = trialBox % MolStartIndx(molIndx) + self%atomNum  - 1
            accept = .false.
            iMol = disp(iDisp)%molIndx
!              write(*,*) "Del", iMol
            do jMol = 1, totalMol
              if(self%newTopoList(jMol, iMol)) then
                self%newTopoList(jMol, iMol) = .false.
                self%newTopoList(iMol, jMol) = .false.
              endif
            enddo
          endif
        enddo


        if(accept) then
          return
        endif

        startMol = 1
        leave = .false.
        do iMol = 1, totalMol
          do jMol = iMol+1, totalMol
            if(self%newTopoList(jMol, iMol)) then
              startMol = iMol
              leave = .true.
              exit
            endif

          enddo
          if(leave) then
            exit
          endif
        enddo
        self%clustMemb = .false.
        self%clustMemb(startMol) = .true.
        self%newlist(1) = startMol
        nNew = 1
        nClust = 1

       !----------------------------------------------------------------------------
      class default
        stop "Distance criteria is not compatiable with this perturbation type."
       !----------------------------------------------------------------------------
    end select

!    write(*,*) "=============================="
!    do iMol = 1, totalMol
!      write(*,*) (self%newTopoList(jMol, iMol), jMol=1,totalMol)
!    enddo
!    write(*,*)

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
!    write(*,*) "Result:", nClust, totalMol, nNew
    select type(disp)
      class is(Deletion)
!          if( (nNew <= 0) .or. (nClust < totalMol-1) ) then
          if( nClust < totalMol-1 ) then
            accept = .false.
            return
          endif
!      write(*,*) "=============================="
!      do iMol = 1, totalMol
!        write(*,*) (self%topoList(jMol, iMol), jMol=1,totalMol)
!      enddo
!      write(*,*)
!      do iMol = 1, totalMol
!        write(*,*) (self%newTopoList(jMol, iMol), jMol=1,totalMol)
!      enddo
!      write(*,*)

      topMol = trialBox % NMol(self%molType)
      iMol = disp(1) % molIndx
     ! write(*,*) topMol, iMol

      if(iMol == topMol) then
        self%newTopoList(topMol,:) = .false.
        self%newTopoList(:,topMol) = .false.       
      else
     
      do jMol = 1, totalMol
        if(self%newTopoList(jMol,topMol)) then
          if(jMol /= iMol) then
            self%newTopoList(jMol,iMol) = .true.
            self%newTopoList(iMol,jMol) = .true.            
          endif
        else
          self%newTopoList(jMol,iMol) = .false.
          self%newTopoList(iMol,jMol) = .false.
        endif        
      enddo
      
      self%newTopoList(iMol,iMol) = .false.

      self%newTopoList(topMol,:) = .false.
      self%newTopoList(:,topMol) = .false.      

      endif
!      do iMol = 1, totalMol
!        write(*,*) (self%newTopoList(jMol, iMol), jMol=1,totalMol)
!      enddo
!      write(*,*)



      class default
!          if( (nNew <= 0) .or. (nClust < totalMol) ) then
          if( nClust < totalMol ) then
            accept = .false.
            return
          endif
    end select
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

    call GetXCommand(line, command, 5, lineStat)
    read(command, *) intVal
    self%atomnum = intVal

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

    self%topoList = self%newTopoList

  end subroutine
!=====================================================================
end module
!=====================================================================
