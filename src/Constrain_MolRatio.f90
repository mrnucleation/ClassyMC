!=============================================================
! The purpose of this constraint is to provide the user a way to control the
! ratio between two different molecule types.  This is useful for forcefields
! where molecular units (IE SiO2) are represented non-bonded potentials
!=============================================================
module Constrain_MolRatio
  use ConstraintTemplate, only: constraint
  use Template_SimBox, only: SimBox
  use VarPrecision

  type, public, extends(constraint) :: MolRatio
    integer :: molType1, molType2

    contains
      procedure, pass :: CheckInitialConstraint => MolRatio_CheckInitialConstraint
!      procedure, pass :: ShiftCheck => MolRatio_ShiftCheck
      procedure, pass :: ProcessIO => MolRatio_ProcessIO
          
    end type

!=============================================================
  contains
!=============================================================
  subroutine MolRatio_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(MolRatio), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    integer :: iAtom
    
    accept = .true.

  end subroutine
!=============================================================
  subroutine MolRatio_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, LowerCaseLine, GetXCommand
    use ParallelVar, only: nout
    implicit none
    class(MolRatio), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    integer :: i, intVal
    real(dp) :: realVal
    character(len=30), allocatable :: parlist(:)
    integer :: buffer, nPar

    lineStat = 0
    call GetAllCommands(line, parlist, nPar, lineStat)
    
    buffer = 0
    i = 2
    do while(i <= size(parlist))

      select case(trim(adjustl(parlist(buffer))))
        case("x")
          self%wallAxis(1) = .true.
          read(parlist(i+1), *) self%xlo
          read(parlist(i+2), *) self%xhi
          buffer = 2

        case("y")
          self%wallAxis(2) = .true.
          read(parlist(i+1), *) self%ylo
          read(parlist(i+2), *) self%yhi
          buffer = 2

        case("z")
          self%wallAxis(3) = .true.
          read(parlist(i+1), *) self%zlo
          read(parlist(i+2), *) self%zhi
          buffer = 2

        case default
          lineStat = -1
          return

      end select

      if(buffer > 0) then
        i = i + buffer
      endif

    enddo
  end subroutine

!=============================================================
end module
!=============================================================
