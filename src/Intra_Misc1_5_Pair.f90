!=============================================================================+
! Intra-molecular potential designed to compute 1-5 pairs within a molecule
!=============================================================================+
module Misc_IntraPair_1_5
  use VarPrecision
  use Template_SimBox, only: SimBox
  use Template_IntraMiscIntra, only: MiscIntra_FF
  use CoordinateTypes

  type, public, extends(MiscIntra_FF) :: Pair_1_5
!    integer, private :: molType
    integer, private, allocatable :: neilist(:,:)
    integer, private, allocatable :: nNei(:)
    contains
      procedure, pass :: Prologue => Pair15_Prologue
      procedure, pass :: DetailedECalc => Pair15_DetailedECalc
  end type

  contains
!=============================================================================+
  subroutine Pair15_Prologue(self)
    use Common_MolInfo, only: MolData, BondData
    use SearchSort, only: QSort
    implicit none
    class(Pair_1_5), intent(inout) :: self
    logical, allocatable :: flipped(:)
    logical, allocatable :: clustMemb(:)
    logical, allocatable :: topoList(:, :)
    integer, allocatable :: newList(:)
    integer, allocatable :: newList2(:)

    integer :: iLimit, iBond, iNext
    integer :: molType, iAtom, jAtom, kAtom
    integer :: nAtoms
    integer :: mem1, mem2
    integer :: nNext, nNew, nClust

    molType = self%GetMolType()
    nAtoms = MolData(moltype)%natoms

    allocate( self%neilist(1:nAtoms, 1:nAtoms) )
    allocate( self%nNei(1:nAtoms) )
    self%neilist = 0
    self%nNei = 0

    allocate(flipped(1:nAtoms))
    allocate(clustMemb(1:nAtoms))
    allocate(topolist(1:nAtoms,1:nAtoms))
    allocate(newList(1:nAtoms))
    allocate(newList2(1:nAtoms))

    topolist = .false.
    do iBond = 1, MolData(molType)%nBonds
      mem1 = MolData(molType)%bond(iBond)%mem1
      mem2 = MolData(molType)%bond(iBond)%mem2
      topolist(mem1, mem2) = .true.
      topolist(mem2, mem1) = .true.
    enddo


    !Borrowed the distance criteria code to map out
    !the 1-5 interactions by performing a simulatnious graph expansion
    !There's some inefficiencies with the loop, but eh...who cares
    !this gets ran once at the start of the program and never again
    !It would require a molecule of insane size to cause problems.
    do iAtom = 1, nAtoms
      nNew = 1
      newlist = 0
      newlist2 = 0
      newlist(1) = iAtom
      nClust = 1
      clustMemb = .false.
      clustMemb(iAtom) = .true.
      do iLimit = 1, nAtoms
        nNext = nNew
        nNew = 0
        do iNext = 1, nNext
          jAtom = newlist(iNext)
          do kAtom = 1, nAtoms
            if(topoList(jAtom, kAtom) )then
              if(.not. clustMemb(kAtom) )then
                clustMemb(kAtom) = .true.
                nNew = nNew + 1
                nClust = nClust + 1
                newlist2(nNew) = kAtom
                if(iLimit >= 4) then
                  self%nNei(iAtom) = self%nNei(iAtom) + 1
                  self%neilist(self%nNei(iAtom), iAtom) = kAtom
                endif
              endif
            endif
          enddo
        enddo
!        write(*,*) newlist(1:nNext)
!        write(*,*) newlist2(1:nNew)
!        write(*,*)

         !If every atom has been added then no further calculations are
         !needed.
        if(nClust >= nAtoms) then
          exit
        endif 
         ! If no new atoms were added, the algorithm has hit a dead end
         ! and the molecule is broken.  
        if(nNew <= 0) then
          write(0,*) "Topology Error! An atom was defined, but not connected to"
          write(0,*) "molecule!"
          write(0,*) "MolType:", molType
          stop 
        endif

        newlist(1:nNew) = newlist2(1:nNew)
      enddo
      call QSort( self%neilist(1:self%nNei(iAtom),iAtom) )
    enddo

    deallocate(flipped)
    deallocate(clustMemb)
    deallocate(topoList)
    deallocate(newList)
    deallocate(newList2)


  end subroutine
!===================================================================
  subroutine Pair15_DetailedECalc(self, curbox, atompos, E_T, accept)
    use Common_MolInfo, only: MolData
    use ForceFieldData, only: ECalcArray
    use SimpleSimBox, only: SimpleBox
    implicit none
    class(Pair_1_5), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(in) :: atompos(:, :)
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    integer :: molType, nAtoms, iAtom, jNei, jAtom
    integer :: atmType1, atmType2 
    real(dp) :: rx, ry, rz, rsq

    class(ECalcArray), pointer :: epointer => null()

      
    select type(curbox)
      class is(SimpleBox)
        call curbox%GetEFunc(epointer)
    end select

    molType = self%GetMolType()
    nAtoms = MolData(molType)%natoms


    accept = .true.
    E_T = 0E0_dp
    do iAtom = 1, nAtoms
      atmType1 = MolData(molType)%AtomType(iAtom)
      do jNei = 1, self%nNei(iAtom)
        jAtom = self%neilist(jNei, iAtom)
        if(iAtom <= jAtom) cycle
        atmType2 = MolData(molType) % AtomType(jAtom)
        rx = atompos(1, iAtom) - atompos(1, jAtom)
        ry = atompos(2, iAtom) - atompos(2, jAtom)
        rz = atompos(3, iAtom) - atompos(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        E_T = E_T + epointer % method % SinglePair(rsq, atmtype2, atmtype1)
      enddo
    enddo

  end subroutine
!=============================================================================+
  subroutine Pair15_ProcessIO(self, line)
    use Input_Format, only: maxLineLen, GetXCommand 
    implicit none
    class(Pair_1_5), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat

!    call GetXCommand(line, command, 2, lineStat)
!    read(command, *) self%theta0
!    self%theta0 = self%theta0*inAngUnit

!    if(lineStat /= 0) then
!      write(*,*) "Missing input rquired for the ridgid Angle style"
!    endif
  end subroutine

!=============================================================================+
end module
!=============================================================================+
