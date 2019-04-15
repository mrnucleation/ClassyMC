!===========================================
module Exceptions
!===========================================
  contains
!===========================================
  subroutine IOError(errCode, errStr)
  implicit none
  integer, intent(in) :: errCode
  character(len=*), intent(in) :: errStr


  select case(errCode)
    case(-1) ! Unknown command
    case(-2) ! Duplicate command
      write(*,"(A,1x,A,1x,A)") "ERROR! The", errStr, "command has already been used and can not be called twice"
  end select

  stop
  end subroutine
!===========================================
  subroutine HardError(errCode, errStr)
  use ParallelVar, only: nout
  implicit none
  integer, intent(in) :: errCode
  character(len=*), intent(in) :: errStr


  select case(errCode)
    case(-1) !Important Array not Allocated
      write(nout, "(A,1x,A,1x,A)") "ERROR!", errStr, " has not been allocated! Please use the appropriate command"
      write(*,"(A,1x,A,1x,A)") "ERROR!", errStr, " has not been allocated! Please use the appropriate command"
  end select

  stop "A fatal error has been encountered.  Please check your input and try again"
  end subroutine
!===========================================
  subroutine ThrowWarning()

  end subroutine
!===========================================
  subroutine Argument(errCode, errStr)
  implicit none
  integer, intent(in) :: errCode
  character(len=*), intent(in) :: errStr


  select case(errCode)
    case(-1) ! Unknown command
    case(-2) ! Duplicate command
      write(*,"(A,1x,A,1x,A)") "ERROR! The", errStr, "command has already been used and can not be called twice"
  end select

  stop
  end subroutine
!===========================================

!===========================================
end module
!===========================================
