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
      class(SimpleBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => RSqList_Constructor 
      procedure, pass :: BuildList => RSqList_BuildList 
      procedure, pass :: GetNewList => RSqList_GetNewList
      procedure, pass :: DeleteMol => RSqList_DeleteMol
  end type

!===================================================================================
  contains
!===================================================================================
  subroutine RSqList_Constructor(self, parentID, rCut)
    use BoxData, only: BoxArray
    use Common_NeighData, only: neighSkin
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: parentID
    real(dp), intent(in), optional :: rCut
    real(dp), parameter :: atomRadius = 0.5E0_dp  !Used to estimate an approximate volume of 
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

    allocate( self%list(1:self%maxNei, 1:self%parent%nMaxAtoms), stat=AllocateStatus )
    allocate( self%nNeigh(1:self%parent%nMaxAtoms), stat=AllocateStatus )

    self%list = 0
    self%nNeigh = 0 
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  end subroutine
!===================================================================================
  subroutine RSqList_BuildList(self)
    implicit none
    class(RSqList), intent(inout) :: self

    call Builder_RSq(self%parent)
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

    do iAtom = 1, MolData(nType)%nAtoms
      atmIndx = nStart + iAtom - 1
      topAtom = topStart + iAtom - 1
      !Remove the deleted from the list of it's neighbors
      do iNei = 1, self % nNeigh(atmIndx)
        curNei = self % list(iNei, atmIndx)
        nNei = self%nNeigh(curNei)
        if(self%sorted) then
          curIndx = BinarySearch( atmIndx, self%list(1:nNei, curNei) )
        else
          curIndx = SimpleSearch( atmIndx, self%list(1:nNei, curNei) )
        endif
!        write(*,*) curIndx, self%list(curIndx, curNei)
        self%list(1:nNei, curNei ) = [self%list(1:curIndx-1, curNei), &
                                      self%list(curIndx+1:nNei, curNei) ]
        self%nNeigh(curNei) = self%nNeigh(curNei) - 1 
      enddo

      !Re-index the top atom to it's new location
      do iNei = 1, self % nNeigh(topAtom)
        curNei = self % list(iNei, topAtom)
        nNei = self%nNeigh(curNei)
        if(self%sorted) then
          curIndx = BinarySearch( topAtom, self%list(1:nNei, curNei) )
        else
          curIndx = SimpleSearch( topAtom, self%list(1:nNei, curNei) )
        endif
        self % list(curIndx, curNei) = atmIndx

      enddo
    enddo

  end subroutine

!===================================================================================
  subroutine RSqList_GetNewList(self, iDisp, tempList, disp)
    use Common_NeighData, only: tempList
    implicit none
    class(RSqList), intent(inout) :: self
    integer, intent(in) :: iDisp
    type(Displacement), intent(inout) :: disp
    integer :: iList
    integer :: iType, jType, iAtom, jAtom, j
    integer :: iUp, iLow, jUp, jLow, molStart, jMolStart, jMolEnd
    real(dp) :: rx, ry, rz, rsq



    disp % newlist = .true.
    disp % listIndex = iDisp

    templist(:, iDisp) = 0
    tempNNei(iDisp) = 0

  end subroutine
!===================================================================================
! End Type Bound
!===================================================================================
  subroutine Builder_RSq(trialBox)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(SimpleBox), intent(inout) :: trialBox
    integer :: iList
    integer :: iType, jType, iAtom, jAtom, j
    integer :: iUp, iLow, jUp, jLow, molStart, jMolStart, jMolEnd
    real(dp) :: rx, ry, rz, rsq


    do iList = 1, size(trialBox%NeighList)
      trialBox%NeighList(iList)%nNeigh = 0
      trialBox%NeighList(iList)%list = 0
    enddo

    molStart = 1
    do iType = 1, nMolTypes
!      write(*,*) "molstart", molStart
      iLow = trialBox%MolStartIndx(molStart)
      iUp = trialBox%MolEndIndx(molStart + trialBox%NMol(iType) - 1) 
!          write(*,*) "i", iType, iLow, iUp
      do iAtom = iLow, iUp
        jMolStart = molStart
        do jType = iType, nMolTypes
          jMolEnd = 1
          do j = 1, jType-1
            jMolEnd = jMolEnd + trialBox%NMolMax(j)
          enddo 
          if(iType == jType) then
            if( trialBox%MolIndx(iAtom)+1 <= jMolEnd + trialBox%NMol(iType) -1  ) then
              jLow = trialBox%MolStartIndx( trialBox%MolIndx(iAtom)+1 )
            else
              cycle
            endif
          else
            jLow = trialBox%TypeFirst(jType)
          endif
          jUp = trialBox%MolEndIndx( jMolEnd + trialBox%NMol(jType) - 1 )
 
!          write(*,*) "i", iType, iLow, iUp
!          write(*,*) "j", jType, jLow, jUp
          do jAtom = jLow, jUp
!            write(*,*) iAtom, jAtom
            rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
            ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
            rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
            call trialBox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            do iList = 1, size(trialBox%NeighList)
              if( rsq <= trialBox%NeighList(iList)%rCutSq ) then 
!                write(*,*) "NEIGH"
!                if(trialBox%NeighList(iList)%nNeigh(iAtom) > trialBox%NeighList(iList)%maxNei) then
!                  write(*,*) "ERROR! NeighList Overflow!"
!                endif
!                if(trialBox%NeighList(iList)%nNeigh(jAtom) > trialBox%NeighList(iList)%maxNei) then
!                  write(*,*) "ERROR! NeighList Overflow!"
!                endif
 
                trialBox%NeighList(iList)%nNeigh(iAtom) = trialBox%NeighList(iList)%nNeigh(iAtom) + 1
                trialBox%NeighList(iList)%list( trialBox%NeighList(iList)%nNeigh(iAtom), iAtom ) = jAtom

                trialBox%NeighList(iList)%nNeigh(jAtom) = trialBox%NeighList(iList)%nNeigh(jAtom) + 1
                trialBox%NeighList(iList)%list( trialBox%NeighList(iList)%nNeigh(jAtom), jAtom ) = iAtom
              endif
            enddo        
          enddo
        enddo
      enddo  

      molStart = molStart + trialBox%NMolMax(iType) 
    enddo


  end subroutine
!===================================================================================
end module
!===================================================================================
