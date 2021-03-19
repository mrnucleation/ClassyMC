!======================================================
module CoordinateTypes
  use VarPrecision
 
  !Base Perturbation Class
  type :: Perturbation
  end type

  !Move type where a single particle's position is changed.
  type, extends(Perturbation) :: Displacement
    integer(kind=atomIntType) :: molType, atomsubIndx
    integer(kind=atomIntType) :: molIndx, atmIndx   
    real(dp) :: x_new, y_new, z_new
    logical :: newList = .false.
    integer :: listIndex = -1
  end type


  !Move type where one particle is removed 
  type, extends(Perturbation) :: Deletion
    integer(kind=atomIntType) :: molType
    integer(kind=atomIntType) :: molIndx, atmIndx   
  end type
  
  !Move type where one particle is created
  type, extends(Perturbation) :: Addition
    integer(kind=atomIntType) :: molType
    integer(kind=atomIntType) :: molIndx, atmIndx   
    real(dp) :: x_new, y_new, z_new
    integer :: listIndex = -1
  end type

  !Move type where an atom of one type is exchanged for an atom of another type.
  type, extends(Perturbation) :: AtomExchange
    integer(kind=atomIntType) :: newAtmIndx, newType
    integer(kind=atomIntType) :: oldAtmIndx, oldType
  end type

  !Move type where the volume of the entire box is changed.
  type, extends(Perturbation) :: VolChange
    real(dp) :: volNew, volOld
  end type

  !Move type where the volume of the entire box by scaling each side by different scale factors
  type, extends(VolChange) :: OrthoVolChange
    real(dp) :: xScale, yScale, zScale
  end type

  !Not Implimented Yet
  type, extends(VolChange) :: TriVolChange
    real(dp) :: xScale, yScale, zScale
  end type



end module
!======================================================
module ParallelVar
  integer :: p_size = 1
  integer :: ierror, tag, seed
  integer :: myid = 0
  integer :: stderr = 0
  integer :: nout = 6
  integer :: paraMode = 1
end module  
!======================================================
module SimControl
  use VarPrecision
  integer :: simType = 1
  integer(kind=8) :: nCycles = 0
  integer(kind=8) :: nMoves = 0

  integer :: energyCheck = -1
  integer :: screenFreq = 1000
  integer :: configFreq = 100
  logical :: printBox = .true.
  logical :: printAcc = .true.

  real(dp) :: TimeStart = 0E0_dp
  real(dp) :: TimeEnd = 0E0_dp

  !Minimization Parameters
  real(dp) :: ETol = 1E-5_dp
  real(dp) :: ForceTol = 1E-5_dp
  real(dp) :: lrate =  1E-5_dp
end module  
!======================================================

