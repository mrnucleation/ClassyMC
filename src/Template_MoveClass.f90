!=========================================================================
module MoveClassDef
  use MasterTemplate, only: classyClass
  use SimpleSimBox, only: SimpleBox
  use VarPrecision

  type, public, extends(classyClass) :: MCMove
    real(dp) :: atmps = 1E-30_dp
    real(dp) :: accpt = 0E0_dp
    real(dp), allocatable :: boxProb(:)

    !Temporary Neighborlist Variables
    integer, allocatable :: tempNnei(:)
    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor
      procedure, pass :: GeneratePosition 
      procedure, pass :: FullMove
      procedure, pass :: GetAcceptRate
      procedure, pass :: GetBoxProb
      procedure, pass :: ProcessIO
      procedure, pass :: LoadBoxInfo
      procedure, pass :: UniformMoleculeSelect
      procedure, pass :: UniformAtomSelect
!      procedure, pass :: Maintenance 
  end type

 contains
!=========================================================================
  subroutine Constructor(self)
    implicit none
    class(MCMove), intent(inout) :: self
  end subroutine
!=========================================================================
  subroutine GeneratePosition(self, disp)
    use CoordinateTypes, only: Perturbation
    implicit none
    class(MCMove), intent(in) :: self
    class(Perturbation), intent(inout) :: disp
  end subroutine
!=========================================================================
  subroutine FullMove(self, trialBox, accept)
    class(MCMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept

    write(0,*) "FullMove for this move type has not been defined!"
    error stop
    accept = .true.
  end subroutine
!=========================================================================
  subroutine LoadBoxInfo(self, trialBox, disp)
    ! Places the information such as boxID into
    use CoordinateTypes, only: Perturbation
    class(MCMove), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    class(Perturbation), intent(inout) :: disp(:)
    integer :: iDisp

    do iDisp = 1, size(disp)
      disp(iDisp)%boxID = trialBox%GetBoxID()
    enddo

  end subroutine
!=========================================================================
  function GetAcceptRate(self) result(rate)
    class(MCMove), intent(in) :: self
    real(dp) :: rate

    if(self%atmps > 0E0_dp) then
      rate = 1E2_dp*self%accpt/self%atmps
    else
      rate = 0E0_dp
    endif

    return
  end function
!=========================================================================
  subroutine GetBoxProb(self, boxProb)
    implicit none
    class(MCMove), intent(inout) :: self
    integer :: iBox
    real(dp), intent(inout) :: boxProb(:)

    do iBox = 1, size(boxProb)
      boxProb(iBox) = self%boxProb(iBox)
    enddo



  end subroutine
!=========================================================================
! Pulls a random atom's index from the box. This is used when all atoms can be
! targeted by this move.
!
  function UniformAtomSelect(self, trialBox) result(nAtom)
    use Box_Utility, only: FindAtom
    use RandomGen, only: grnd
    implicit none
    class(MCMove), intent(inout) :: self
    class(SimpleBox), intent(in) :: trialBox
    integer :: nAtom
    integer :: rawIndx

    rawIndx = floor( trialBox%nAtoms * grnd() + 1E0_dp )
    call FindAtom(trialbox, rawIndx, nAtom)

  end function
!=========================================================================
! Pulls a random molecule's index from the box. This is used when all molecules can be
! targeted by this move.
!
  function UniformMoleculeSelect(self, trialBox, restrict) result(nMol)
    use Box_Utility, only: FindMolecule
    use RandomGen, only: grnd, ListRNG
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(MCMove), intent(in) :: self
    class(SimpleBox), intent(in) :: trialBox
    logical, intent(in), optional :: restrict(1:nMolTypes)
    integer :: nMol, nType, iType, nActive
    integer :: rawIndx
    real(dp) :: weights(1:nMolTypes), norm

    if(.not. present(restrict)) then
      rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp )
      call FindMolecule(trialbox, rawIndx, nMol)
    else
      nActive = 0
      do iType = 1, nMolTypes
        if( restrict(iType) ) then
          nActive = nActive + 1
          weights(iType) = real(trialBox%NMol(iType), dp)
          norm = norm + weights(iType)
        else
          weights(iType) = 0E0_dp
        endif
        nType = ListRng(weights(1:nMolTypes), norm)
        nMol = floor( trialBox%NMol(nType) * grnd() + 1E0_dp )
        nMol = trialBox%MolGlobalIndx(nType, nMol)
      enddo
    endif

  end function

!=========================================================================
  subroutine ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen
    implicit none
    class(MCMove), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat

    lineStat = 0
  end subroutine
!=========================================================================
!  subroutine Maintenance(self)
!    implicit none
!    class(MCMove), intent(inout) :: self
!  end subroutine
!=========================================================================
end module
!=========================================================================
