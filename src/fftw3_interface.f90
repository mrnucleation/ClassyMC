!================================================================================
! Minimal FFTW3 Fortran Interface for P3M
! This provides the basic FFTW3 functions needed for P3M electrostatics
!================================================================================
module fftw3_interface
  use, intrinsic :: iso_c_binding
  implicit none

  ! FFTW planning rigor flags
  integer(C_INT), parameter :: FFTW_ESTIMATE = 64
  integer(C_INT), parameter :: FFTW_MEASURE = 0
  integer(C_INT), parameter :: FFTW_PATIENT = 32
  integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8

  interface
    ! Create 3D real-to-complex plan
    function fftw_plan_dft_r2c_3d(n0, n1, n2, in, out, flags) bind(C, name='fftw_plan_dft_r2c_3d')
      import :: C_INT, C_PTR
      integer(C_INT), value :: n0, n1, n2
      type(C_PTR), value :: in, out
      integer(C_INT), value :: flags
      type(C_PTR) :: fftw_plan_dft_r2c_3d
    end function

    ! Execute real-to-complex transform
    subroutine fftw_execute_dft_r2c(plan, in, out) bind(C, name='fftw_execute_dft_r2c')
      import :: C_PTR
      type(C_PTR), value :: plan, in, out
    end subroutine

    ! Destroy plan
    subroutine fftw_destroy_plan(plan) bind(C, name='fftw_destroy_plan')
      import :: C_PTR
      type(C_PTR), value :: plan
    end subroutine
  end interface

end module fftw3_interface

