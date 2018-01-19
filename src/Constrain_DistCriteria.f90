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
    integer :: neighList = -1
    integer :: molType = -1
    real(dp) :: rCut, rCutSq
    integer :: boxID = -1

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
    class(DistCriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    integer :: totalMol, nNew, nClust, neiIndx
    integer :: iMol,jMol, iAtom, jAtom, iLimit
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    self%flipped = .false.
    self%clustMemb = .false.

    totalMol = trialBox%NMol(molType)
    neiIndx = self%neighList


    nClust = 0
    do iLimit = 1, totalMol
      nNew = 0
      do iMol = 1, totalMol
        iAtom = 
        if( self%clustMemb(iMol) /= self%flipped) then
          do jMol = 1, totalMol
            if(.not. self%clustMemb(jMol) )then
              
            endif
          enddo
        endif
      enddo

      if(nNew <= 0) then
        exit
      endif

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
    class(DistCriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept
    
    real(dp) :: rx, ry, rz, rsq
    

    accept = .true.

  end subroutine
!=====================================================================
  subroutine DistCrit_NewCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=====================================================================
  subroutine DistCrit_OldCheck(self, trialBox, disp, accept)
    implicit none
    class(DistCriteria), intent(in) :: self
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
