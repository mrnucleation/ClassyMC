!======================================================
module CoordinateTypes
  use VarPrecision

  type Displacement
    logical :: newAtom = .false.
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
    real(dp) :: x_new, y_new, z_new
  end type

end module
!======================================================
module ParallelVar
  integer :: myid, p_size, ierror, tag, nout, seed
end module  
!======================================================
module Constants
  use VarPrecision
        
  real(dp), parameter :: pi = 4E0_dp * atan(1E0_dp)
  real(dp), parameter :: two_pi = 8E0_dp * atan(1E0_dp)
end module  

!======================================================

