!=============================================================
! The purpose of this constraint is to provide the user a way to constrain all or parts of the system by a flat wall.
!
module Constrain_HardWall
  use ConstraintTemplate, only: constraint
  use Template_SimBox, only: SimBox
  use VarPrecision

  type, public, extends(constraint) :: hardwall
    integer, allocatable :: atmTypes(:)
    logical :: wallAxis(1:3) = [.false., .false., .false.]
    real(dp) :: xhi, xlo
    real(dp) :: yhi, ylo
    real(dp) :: zhi, zlo
    contains
      procedure, pass :: CheckInitialConstraint => HardWall_CheckInitialConstraint
!      procedure, pass :: ShiftCheck => HardWall_ShiftCheck
      procedure, pass :: ProcessIO => HardWall_ProcessIO
          
    end type

!=============================================================
  contains
!=============================================================
  subroutine HardWall_CheckInitialConstraint(self, trialBox, accept)
    implicit none
    class(hardwall), intent(inout) :: self
    class(SimBox), intent(in) :: trialBox
    logical, intent(out) :: accept

    integer :: iAtom
    
    accept = .true.

!    do iAtom = 1, trialBox%
!    enddo

  end subroutine
!=============================================================
  subroutine HardWall_ProcessIO(self, line, lineStat)
    use Input_Format, only: maxLineLen, LowerCaseLine, GetAllCommands
    use ParallelVar, only: nout
    implicit none
    class(HardWall), intent(inout) :: self
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
