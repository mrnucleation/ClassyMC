!=============================================================
! This module constraints the total number of molecules allowed in a box
! It is primarily used in a multi component system where you wish
! to allow
!=============================================================
module Constrain_MolTotal
  use ConstraintTemplate, only: constraint
  use CoordinateTypes, only: Perturbation
  use CoordinateTypes, only: Deletion, Addition
  use Template_SimBox, only: SimBox
  use VarPrecision

  type, public, extends(constraint) :: MolTotal
    integer, private :: nLimit

    contains
      procedure, pass :: CheckInitialConstraint => MolTotal_CheckInitialConstraint
      procedure, pass :: DiffCheck => MolTotal_DiffCheck
      procedure, pass :: ProcessIO => MolTotal_ProcessIO
          
    end type

!=============================================================
  contains
!=============================================================
  subroutine MolTotal_CheckInitialConstraint(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(MolTotal), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    logical, intent(out) :: accept


    if(trialBox%nMolTotal > self%nLimit) then
      accept = .false.
    endif

    

  end subroutine
!=============================================================
  subroutine MolTotal_DiffCheck(self, trialBox, disp, accept)
    implicit none
    class(MolTotal), intent(inout) :: self
    class(SimBox), intent(inout) :: trialBox
    class(Perturbation), intent(in) :: disp(:)
    logical, intent(out) :: accept
    integer :: nNewTotal


    !This section creates the topology list of the new state using information
    !based on the what kind of perturbation was performed.




    select type(disp)
      class is(Addition)
        nNewTotal = trialBox%nMolTotal + 1
      class is(Deletion)
        nNewTotal = trialBox%nMolTotal - 1
      class default
        nNewTotal = trialBox%nMolTotal
    end select

    if(nNewTotal > self%nLimit) then
      accept = .false.
    endif
    accept = .true.

  end subroutine

!=============================================================
  subroutine MolTotal_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, LowerCaseLine, GetXCommand
    use ParallelVar, only: nout
    implicit none
    class(MolTotal), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    real(dp) :: realVal
    integer :: newLimit

    lineStat = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) newLimit
    if(newLimit < 1) then
      write(0,*) "ERROR! THe user specified an invalid MolTotal Limit!"
      stop
    endif
    self%nLimit = newLimit
    write(nout, *) "Total Molecule Limit:", newLimit

   
  end subroutine
!=============================================================
end module
!=============================================================
