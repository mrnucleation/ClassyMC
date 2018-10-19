!====================================================================
!This module contains the FreezeType constraint that prevents an entire molecule type from
!moving.
!====================================================================
module Constrain_FreezeType
  use VarPrecision
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use CoordinateTypes, only: DisplacementNew, Deletion, Addition
  use Template_SimBox, only: SimBox
  use ParallelVar, only: nout

  type, public, extends(constraint) :: FreezeType
    integer :: boxID = -1
    integer :: molType = -1
    class(SimBox), pointer :: parent => null()
    contains
      procedure, pass :: Constructor => FreezeType_Constructor
!      procedure, pass :: CheckInitialConstraint => FreezeType_CheckInitialConstraint
      procedure, pass :: DiffCheck => FreezeType_DiffCheck
!      procedure, pass :: CheckCluster => FreezeType_CheckCluster
      procedure, pass :: ProcessIO => FreezeType_ProcessIO
!      procedure, pass :: Maintenance => FreezeType_Maintenance
!      procedure, pass :: Update => FreezeType_Update
!      procedure, pass :: Epilogue => FreezeType_Epilogue
  end type
!=====================================================================
  contains
!=====================================================================
  subroutine FreezeType_Constructor(self, boxID)
    use BoxData, only: BoxArray
    implicit none
    class(FreezeType), intent(inout) :: self
    integer, intent(in) :: boxID
    integer :: AllocateStat
    integer :: nMolMax


    self%boxID = boxID
    self%parent => BoxArray(boxID) % box 

    IF (AllocateStat /= 0) STOP "Allocation Error in FreezeType Constraint"
  end subroutine
!=====================================================================
  subroutine FreezeType_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(FreezeType), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!=============================================================
  subroutine FreezeType_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(FreezeType), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    integer :: iDisp




    !This section creates the topology list of the new state using information
    !based on the what kind of perturbation was performed.
    select type(disp)
       !----------------------------------------------------------------------------
      class is(DisplacementNew)
        do iDisp = 1, size(disp)
          if(disp(iDisp)%molType == self%molType) then
            accept = .false.
            return
          endif
        enddo
       !----------------------------------------------------------------------------
      class is(Addition)
        do iDisp = 1, size(disp)
          if(disp(iDisp)%molType == self%molType) then
            accept = .false.
            return
          endif
        enddo

       !----------------------------------------------------------------------------
      class is(Deletion)
        do iDisp = 1, size(disp)
          if(disp(iDisp)%molType == self%molType) then
            accept = .false.
            return
          endif
        enddo

       !----------------------------------------------------------------------------
      class default
        stop "Distance criteria is not compatiable with this perturbation type."
       !----------------------------------------------------------------------------
    end select
    accept = .true.

  end subroutine
!=============================================================
  subroutine FreezeType_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ParallelVar, only: nout
    implicit none
    class(FreezeType), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%molType = intVal

  end subroutine
!====================================================================
!  subroutine FreezeType_Maintenance(self)
!    implicit none
!    class(FreezeType), intent(inout) :: self
!
!  end subroutine
!====================================================================
  subroutine FreezeType_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(FreezeType), intent(inout) :: self
    logical :: accept


    write(nout, *) "Freezing Molecule Type", self%molType


  end subroutine
!====================================================================
!  subroutine FreezeType_Epilogue(self)
!    implicit none
!    class(FreezeType), intent(inout) :: self
!    logical :: accept
!
!  end subroutine
!=====================================================================
end module
!=====================================================================
