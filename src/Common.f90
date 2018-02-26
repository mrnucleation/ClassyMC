!======================================================
module CoordinateTypes
  use VarPrecision

  type :: Perturbation
    
  end type

  type, extends(Perturbation) :: Displacement
    logical :: newAtom = .false.
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
    real(dp) :: x_new, y_new, z_new

    logical :: oldAtom = .false.
    integer(kind=atomIntType) :: OldmolType, OldatmIndx, OldmolIndx

    logical :: newList = .false.
    integer :: listIndex = -1
  end type

  type, extends(Perturbation) :: DisplacementNew
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
    real(dp) :: x_new, y_new, z_new

    integer(kind=atomIntType) :: OldmolType, OldatmIndx, OldmolIndx

    logical :: newList = .false.
    integer :: listIndex = -1
  end type

  type, extends(Perturbation) :: Exchange
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
    real(dp) :: x_new, y_new, z_new
    integer(kind=atomIntType) :: OldmolType, OldatmIndx, OldmolIndx

    logical :: newList = .false.
    integer :: listIndex = -1
  end type

  type, extends(Perturbation) :: Deletion
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
  end type

  type, extends(Perturbation) :: Addition
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
    real(dp) :: x_new, y_new, z_new
    integer :: listIndex = -1
  end type

  type, extends(Perturbation) :: VolChange
    real(dp) :: volNew, volOld
  end type


end module
!======================================================
module ParallelVar
  integer :: myid, p_size, ierror, tag, nout, seed
  integer :: paraMode = 1
end module  
!======================================================
module SimControl
  integer :: simType = 1
  integer(kind=8) :: nCycles = 0
  integer(kind=8) :: nMoves = 0

  integer :: screenFreq = 1000
  logical :: printBox = .true.
  logical :: printAcc = .true.

end module  
!======================================================

