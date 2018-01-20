!====================================================================
!This module contains the Stilinger distance criteria that is used to 
!enforce clustering. 
!====================================================================
module Constrain_DistanceCriteria
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Displacement
  use Template_SimBox, only: SimBox

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

    totalMol = trialBox%NMol(molType)
    if(totalMol < 2) then
      return
    endif

    neiIndx = self%neighList

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
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    integer :: molType

    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    totalMol = trialBox%NMol(molType)
    if(totalMol < 2) then
      return
    endif

    do i = 1, size(disp)
      if(disp(i)%molType == self%molType) then
        if(disp(i)%atmIndx == 1) then
          accept = .false.
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

     !Seed the initial cluter check by adding the first particle in the array
     !to the cluster
    self%clustMemb(molIndx) = .true.
    self%flipped(molIndx) = .true.
    nClust = nClust + 1
    iAtom = trialBox % MolStartIndx(molIndx)
    nNew = 0
    do jMol = 1, totalMol
      jMolIndx = trialBox % MolGlobalIndx(self%molType, jMol)
      jAtom = trialBox % MolStartIndx(jMolIndx)
      rx = trialBox%atoms(1, iAtom) - trialBox%atoms(1, jAtom)
      ry = trialBox%atoms(2, iAtom) - trialBox%atoms(2, jAtom)
      rz = trialBox%atoms(3, iAtom) - trialBox%atoms(3, jAtom)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < self%rCutSq) then
        self%clustMemb(jMol) = .true.
        nNew = nNew + 1
        nClust = nClust + 1
      endif
    enddo

    if( (nNew <= 0) .or. (nClust < totalMol) ) then
      accept = .false.
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
    implicit none
    class(DistCriteria), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % neighList = intVal

      case("moleculetype")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % moltype = intVal

      case("rcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % rCut = realVal
        self % rCutSq = realVal**2

      case default
        lineStat = -1
    end select

  end subroutine
!=====================================================================
end module
!=====================================================================
