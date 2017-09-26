!======================================================
      module CoordinateTypes
      use VarPrecision

      type Displacement
        integer(kind=atomIntType) :: molType, atmIndx, molIndx
        real(dp) :: x_new, y_new, z_new
      end type


      end module
!======================================================
      module ParallelVar
        integer :: myid, p_size, ierror, tag, nout, seed
      end module  

!======================================================

