!===========================================================================
module C_F_Routines
  use iso_c_binding, only: c_ptr, c_char, c_associated, c_f_pointer, c_int
  public :: C_F_String
!===========================================================================
contains
!===========================================================================
  function C_F_String(c_str) result(f_str)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char
    type(c_ptr), intent(in) :: c_str
    character(:,kind=c_char), pointer :: f_str
    character(kind=c_char), pointer :: arr(:)
    interface
      ! steal std c library function rather than writing our own.
      function strlen(s) bind(c, name='strlen')
        use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
        implicit none
        !----
        type(c_ptr), intent(in), value :: s
        integer(c_size_t) :: strlen
      end function strlen
    end interface
    !****
    call c_f_pointer(c_str, arr, [strlen(c_str)])
    call get_scalar_pointer(size(arr), arr, f_str)
  end function c_f_string
!===========================================================================
  subroutine get_scalar_pointer(scalar_len, scalar, ptr)
    use, intrinsic :: iso_c_binding, only: c_char
    integer, intent(in) :: scalar_len
    character(kind=c_char,len=scalar_len), intent(in), target :: scalar(1)
    character(:,kind=c_char), intent(out), pointer :: ptr
    !***
    ptr => scalar(1)
  end subroutine get_scalar_pointer
  !===========================================================================
end module
