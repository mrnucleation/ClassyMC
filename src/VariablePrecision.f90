!======================================================
!     This module 
      module VarPrecision
        integer, parameter :: dp = kind(0.0d0)
        integer, parameter :: qp = selected_real_kind(33, 4931)
!        integer, parameter :: dp = 8
        integer, parameter :: loopInt = 8
        integer, parameter :: atomIntType = kind(0)
      end module
!======================================================
