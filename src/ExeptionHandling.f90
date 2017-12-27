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
  subroutine HardError()

  end subroutine
!===========================================
  subroutine ThrowWarning()

  end subroutine
!===========================================
end module
!===========================================
