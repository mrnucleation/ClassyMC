!===================================================================================
! This module contains a simple neighborlist
!===================================================================================
module RSqListDef
use VarPrecision
use CoordinateTypes
use Template_SimBox, only: SimBox
use SimpleSimBox, only: SimpleBox
use Template_NeighList, only: NeighListDef

  type, public, extends(NeighListDef) :: RSqList
!      logical :: Sorted = .false.
!      logical :: Strict = .false.
!      integer, allocatable :: list(:,:)
!      integer, allocatable :: nNeigh(:)
!      integer :: maxNei
!      real(dp) :: rCut, rCutSq
!      logical :: restrictType = .false.
!      integer, allocatable :: allowed(:)
      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => RSqList_Constructor 
      procedure, pass :: BuildList => RSqList_BuildList 
      procedure, pass :: GetNewList => RSqList_GetNewList
      procedure, pass :: AddMol => RSqList_AddMol
      procedure, pass :: GetNeighCount => RSqList_GetNeighCount
!      procedure, pass :: TransferList
      procedure, pass :: DeleteMol => RSqList_DeleteMol
      procedure, pass :: Prologue => RSqList_Prologue
      procedure, pass :: Update => RSqList_Update
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine RSqList_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 1.2E0_dp  !Used to estimate an approximate volume of 
    integer :: AllocateStatus

    self%parent => BoxArray(parentID)%box
    if(.not. allocated(self%parent%atoms) ) then
      stop
    endif

!     If no rCut value is given by the subroutine call attempt to pull
!     the rSq value from the parent box's energy function. The assumption being
!     that the neighborlist is used for the energy calculation routines.
    if( present(rCut) ) then
      self % rCut = rCut
      self % rCutSq = rCut * rCut
      self % maxNei = ceiling(rCut**3/atomRadius**3)
    else
      if(self%rCut > 0E0_dp) then
        self % rCutSq = (self%rCut)**2
        self % maxNei = ceiling(self%rCut**3/atomRadius**3)

      else
        self % rCut = self % parent % EFunc % Method % GetCutOff() + neighSkin
        self % rCutSq = (self%rCut)**2
        self % maxNei = ceiling(self%rCut**3/atomRadius**3)
      endif
    endif
 
    write(nout,*) "Neighbor List CutOff:", self%rCut
    if(self%maxNei > self%parent%nMaxAtoms) then
      self%maxNei = self%parent%nMaxAtoms
    endif
    write(nout,*) "Neighbor List Maximum Neighbors:", self%maxNei

    allocate( self%list(1:self%maxNei, 1:self%parent%nMaxAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:self%parent%nMaxAtoms), stat=AllocateStatus )

    if(.not. allocated(self%allowed) ) then
      allocate(self%allowed(1:nAtomTypes), stat=AllocateStatus )
      self%allowed = .true.
    endif

    self%list = 0
    self%nNeigh = 0 
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    self%restrictType = .false.
  end subroutine
!===================================================================================
  subroutine RSqList_Prologue(self)
    implicit none
    class(RSqList), intent(inout) :: self


    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine RSqList_Update(self)
    implicit none
    class(RSqList), intent(inout) :: self


    call self%DumpList(2)
  end subroutine
!===================================================================================
  subroutine RSqList_BuildList(self)
    implicit none
    class(RSqList), intent(inout) :: self


    call Builder_RSq(self%parent)
  end subroutine
!===================================================================================
  function RSqList_GetNeighCount(self, nAtom, rCount) result(nCount)
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: nAtom
    real(dp), intent(in) ,optional :: rCount
    integer :: iNei
    real(dp) :: rCut, rCutSq
    real(dp) :: rx, ry, rz, rsq
    integer :: nCount, jAtom

    if(present(rCount)) then
      rCut = rCount
      rCutSq = rCount*rCount
    else
      rCut = self%rCut
      rCutSq = self%rCutSq
    endif

    nCount = 0
    do iNei = 1, self % nNeigh(nAtom)
      jAtom = self%list(iNei, nAtom) 
      rx = self%parent%atoms(1, nAtom) - self%parent%atoms(1, jAtom)
      ry = self%parent%atoms(2, nAtom) - self%parent%atoms(2, jAtom)
      rz = self%parent%atoms(3, nAtom) - self%parent%atoms(3, jAtom)
      call self%parent%Boundary(rx,ry,rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < rCutSq) then
        nCount = nCount + 1
      endif
    enddo

  end function
!===================================================================================
  subroutine RSqList_AddMol(self, disp, tempList, tempNNei)
    implicit none
    class(RSqList), intent(inout) :: self
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)

    select type(disp)

      class is(Addition)
!        write(*,*) "Add"
        call UpdateList_AddMol_RSq(self%parent, disp, tempList, tempNNei)
    end select

  end subroutine
!===================================================================================
  subroutine RSqList_DeleteMol(self, molIndx, topIndx)
    use Common_MolInfo, only: nMolTypes, MolData
    use SearchSort, only: BinarySearch, SimpleSearch
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: molIndx, topIndx
    integer :: iAtom, iNei, jNei, nType
    integer :: nStart, topStart
    integer :: atmIndx, topAtom
    integer :: curNei, curIndx, nNei

    nStart = self % parent % MolStartIndx(molIndx)
    topStart = self % parent % MolStartIndx(topIndx)
    nType = self % parent % MolType(nStart)

!    write(*,*) "Delete"
    do iAtom = 1, MolData(nType)%nAtoms
      atmIndx = nStart + iAtom - 1
      topAtom = topStart + iAtom - 1
      !Remove the deleted from the list of it's neighbors
      do iNei = 1, self % nNeigh(atmIndx)
        curNei = self % list(iNei, atmIndx)
        nNei = self%nNeigh(curNei)

        if(nNei <= 1) then
          self%nNeigh(curNei) = 0
          self%list(:, curNei) = 0
          cycle
        endif

        if(self%sorted) then
          curIndx = BinarySearch( atmIndx, self%list(1:nNei, curNei) )
        else
          curIndx = SimpleSearch( atmIndx, self%list(1:nNei, curNei) )
        endif
        if(curIndx /= 0) then
          if(nNei > 2) then
            self%list(1:nNei, curNei ) = [self%list(1:curIndx-1, curNei), &
                                          self%list(curIndx+1:nNei, curNei) ]
          else
            if(curIndx == 1) then
              self%list(1, curNei) = self%list(2,curNei)
            endif
          endif
          self%nNeigh(curNei) = self%nNeigh(curNei) - 1 
        endif
      enddo
      
      !Move the top atom into the spot formerly taken up by the deleted atom.
      do iNei = 1, self % nNeigh(topAtom)
        self%list(iNei, atmIndx) = self%list(iNei, topAtom)
      enddo
      self%nNeigh(atmIndx) = self%nNeigh(topAtom)

      !Re-index the neighbor's of the top atom to it's new array location
      do iNei = 1, self % nNeigh(topAtom)
        curNei = self % list(iNei, topAtom)
        nNei = self%nNeigh(curNei)
        if(nNei == 0 ) then
          cycle
        endif
        if(self%sorted) then
          curIndx = BinarySearch( topAtom, self%list(1:nNei, curNei) )
        else
          curIndx = SimpleSearch( topAtom, self%list(1:nNei, curNei) )
        endif
        if(curIndx /= 0) then
          self % list(curIndx, curNei) = atmIndx
        endif
      enddo
      self%nNeigh(topAtom) = 0
    enddo

    self % sorted = .false.

  end subroutine
!===================================================================================
  subroutine RSqList_GetNewList(self, iDisp, tempList, tempNNei, disp, nCount, rCount)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: iDisp
!    type(Displacement), intent(inout) :: disp
    class(Perturbation), intent(inout) :: disp
    integer, intent(inout) :: tempList(:,:), tempNNei(:)
    integer, optional :: nCount
    real(dp), optional :: rCount
    integer :: jType, jAtom, j, iAtom
    integer :: jUp, jLow, molStart, molIndx, jMol
    real(dp) :: xn, yn, zn
    real(dp) :: rx, ry, rz, rsq

    if(present(nCount)) then
      nCount = 0
    endif
    select type(disp)
      class is (Addition)
!        disp % newlist = .true.
        disp % listIndex = iDisp
        iAtom = disp%atmIndx
        molIndx = self%parent%MolIndx(iAtom)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
      class is (DisplacementNew)
        disp % newlist = .true.
        disp % listIndex = iDisp
        iAtom = disp%atmIndx
        molIndx = self%parent%MolIndx(iAtom)
        xn = disp%x_new
        yn = disp%y_new
        zn = disp%z_new
    end select

    templist(:, iDisp) = 0
    tempNNei(iDisp) = 0

    molStart = 1
    do jType = 1, nMolTypes
      do jAtom = 1, self%parent%nMaxAtoms
        if( self%parent%MolSubIndx(jAtom) == molIndx ) then
          cycle
        endif
        if( self%parent%MolSubIndx(jAtom) > self%parent%NMol(self%parent%MolType(jAtom)) ) then
          cycle
        endif
        rx = xn - self%parent%atoms(1, jAtom)
        ry = yn - self%parent%atoms(2, jAtom)
        rz = zn - self%parent%atoms(3, jAtom)
        call self%parent%Boundary(rx,ry,rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          tempNNei(iDisp) = tempNNei(iDisp) + 1
          templist(tempNNei(iDisp), iDisp) = jAtom
        endif
        if(present(rCount)) then
          if(rsq < rCount*rCount) then
            nCount = nCount + 1
          endif
        endif
      enddo
      molStart = molStart + self%parent%NMolMax(jType)
    enddo
    write(*,*) "NewList:", templist(1:tempNNei(iDisp), iDisp)
!    if(present(nCount)) then
!        write(*,*) nCount
!    endif
  end subroutine
!====================================================================
  subroutine RSqList_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetAllCommands, GetXCommand,maxLineLen
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: i, intVal, nPar
    real(dp) :: realVal

    character(len=30) :: command 
    character(len=30), allocatable :: parlist(:)


    lineStat = 0
    call GetXCommand(line, command, 6, lineStat)
    select case( trim(adjustl(command)) )
      case("rcut")
        call GetXCommand(line, command, 7, lineStat)
        read(command,*) realVal
        self%rCut = realVal
        self%rCutSq = realVal * realVal

      case("restricttype")
        call GetAllCommands(line, parlist, nPar, lineStat)
        self%restrictType = .true.
        if(.not. allocated(self%allowed) ) then
          allocate(self%allowed(1:nAtomTypes) )
        endif
        self%allowed = .false.
        do i = 7, size(parList)
          read(parList(i), *) intVal
          self%allowed(intVal) = .true.
        enddo

      case default
        lineStat = -1
    end select

  end subroutine
!===================================================================================
! End Type Bound
!===================================================================================
  subroutine Builder_RSq(trialBox)
    use Common_MolInfo, only: nMolTypes, MolData
    use ParallelVar, only: nout
    implicit none
    class(SimpleBox), intent(inout) :: trialBox
    integer :: iList
    integer :: iType, jType, iAtom, jAtom, j
    integer :: iUp, iLow, jUp, jLow, molStart, jMolStart, jMolEnd, atmType
    real(dp) :: rx, ry, rz, rsq


    do iList = 1, size(trialBox%NeighList)
      trialBox%NeighList(iList)%nNeigh = 0
      trialBox%NeighList(iList)%list = 0
    enddo

    do iAtom = 1, trialBox%nMaxAtoms-1
      if( trialBox%MolSubIndx(iAtom) > trialBox%NMol(trialBox%MolType(iAtom)) ) then
        cycle
      endif
      do jAtom = iAtom+1, trialBox%nMaxAtoms
        if( trialBox%MolSubIndx(jAtom) > trialBox%NMol(trialBox%MolType(jAtom)) ) then
          cycle
        endif
        rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
        ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
        rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
        call trialBox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        do iList = 1, size(trialBox%NeighList)
          if( trialBox % NeighList(iList) % restrictType ) then
            atmType = trialBox % atomType(iAtom)
            if( trialBox%NeighList(iList)%allowed(atmType)  ) then
              cycle
            endif

            atmType = trialBox % atomType(jAtom)
            if( trialBox%NeighList(iList)%allowed(atmType)  ) then
              cycle
            endif
          endif
          if( rsq <= trialBox%NeighList(iList)%rCutSq ) then 
            trialBox%NeighList(iList)%nNeigh(iAtom) = trialBox%NeighList(iList)%nNeigh(iAtom) + 1
            if(trialBox%NeighList(iList)%nNeigh(iAtom) > trialBox%NeighList(iList)%maxNei) then
              write(nout, *) "Neighborlist overflow!"
            endif
            trialBox%NeighList(iList)%list( trialBox%NeighList(iList)%nNeigh(iAtom), iAtom ) = jAtom

            trialBox%NeighList(iList)%nNeigh(jAtom) = trialBox%NeighList(iList)%nNeigh(jAtom) + 1
            if(trialBox%NeighList(iList)%nNeigh(jAtom) > trialBox%NeighList(iList)%maxNei) then
              write(nout, *) "Neighborlist overflow!"
            endif
            trialBox%NeighList(iList)%list( trialBox%NeighList(iList)%nNeigh(jAtom), jAtom ) = iAtom
          endif
        enddo
      enddo  
    enddo


    do iList = 1, size(trialBox%NeighList)
      trialBox%NeighList(iList)%sorted = .true.
    enddo
  end subroutine
!===================================================================================
  subroutine UpdateList_AddMol_RSq(trialBox, disp, tempList, tempNNei)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: trialBox
    class(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    integer :: iList, iDisp, iAtom, iNei, nNei, neiIndx
    real(dp) :: rx, ry, rz, rsq

!    write(*,*) tempNNei(:)
!    write(*,*) tempList(:,1)

    do iList = 1, size(trialBox%NeighList)
      if(iList == 1) then
        do iDisp = 1, size(disp)
          iAtom = disp(iDisp)%atmIndx
          trialBox % NeighList(iList) % nNeigh(iAtom) = tempNNei(iDisp)
          do iNei = 1, tempNNei(iDisp)
            neiIndx = tempList(iNei, iDisp)
!            write(*,*) iAtom, neiIndx
            trialBox % NeighList(iList) % list(iNei, iAtom) =  neiIndx
            trialBox % NeighList(iList) % list( trialBox%NeighList(iList)%nNeigh(neiIndx)+1, neiIndx ) = iAtom
            trialBox%NeighList(iList)%nNeigh(neiIndx)= trialBox%NeighList(iList)%nNeigh(neiIndx) + 1
          enddo
        enddo
      endif
!      write(*,*) "N", trialBox%NeighList(iList)%nNeigh(:)     
    enddo
  end subroutine
!===================================================================================
end module
!===================================================================================
