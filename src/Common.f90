!======================================================
module CoordinateTypes
  use VarPrecision
 
  !Base Perturbation Class
  type :: Perturbation
    
  end type

  !Move type where a single particle's position is changed.
  type, extends(Perturbation) :: Displacement
    integer(kind=atomIntType) :: molType, atmIndx, molIndx, atmSubIndx
    real(dp) :: x_new, y_new, z_new
    logical :: newList = .false.
    integer :: listIndex = -1
  end type

  !Move type where one particle is removed and another is inserted
!  type, extends(Perturbation) :: Exchange
!    integer(kind=atomIntType) :: molType, atmIndx, molIndx
!    real(dp) :: x_new, y_new, z_new
!    integer(kind=atomIntType) :: OldmolType, OldatmIndx, OldmolIndx
!
!    logical :: newList = .false.
!    integer :: listIndex = -1
!  end type

  !Move type where one particle is removed 
  type, extends(Perturbation) :: Deletion
    integer(kind=atomIntType) :: molType, molIndx
  end type
  
  !Move type where one particle is created
  type, extends(Perturbation) :: Addition
    integer(kind=atomIntType) :: molType, atmIndx, molIndx
    real(dp) :: x_new, y_new, z_new
    integer :: listIndex = -1
  end type

  !Move type where an atom of one type is exchanged for an atom of another type.
  type, extends(Perturbation) :: Exchange
    type(Addition) :: inAtom
    type(Deletion) :: outAtom
  end type

  !Move type where the volume of the entire box is changed.
  type, extends(Perturbation) :: VolChange
    real(dp) :: volNew, volOld
  end type

  !Move type where the vol
  type, extends(VolChange) :: OrthoVolChange
    real(dp) :: xScale, yScale, zScale
  end type

  !Move type where the vol
  type, extends(VolChange) :: TriVolChange
    real(dp) :: xScale, yScale, zScale
  end type



end module
!======================================================
module ParallelVar
  integer :: p_size, ierror, tag, seed
  integer :: myid = 0
  integer :: stderr = 0
  integer :: nout =6
  integer :: paraMode = 1
end module  
!======================================================
module SimControl
  use VarPrecision
  integer :: simType = 1
  integer(kind=8) :: nCycles = 0
  integer(kind=8) :: nMoves = 0

  integer :: screenFreq = 1000
  logical :: printBox = .true.
  logical :: printAcc = .true.

  real(dp) :: TimeStart = 0E0_dp
  real(dp) :: TimeEnd = 0E0_dp
end module  
!======================================================

