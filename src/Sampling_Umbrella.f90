!====================================================================
module UmbrellaRule
  use VarPrecision
  use CoordinateTypes, only: Displacement
  use AcceptRuleTemplate, only: acceptrule
 
  type, public, extends(acceptrule) :: Umbrella

    integer :: nBiasVar = 0
    integer, allocatable :: AnalysisIndex(:)

    integer :: umbrellaLimit = 0
    real(dp), allocatable :: UBias(:)
    real(dp), allocatable :: UHist(:)
    real(dp), allocatable :: UBinSize(:)
    real(dp), allocatable :: indexCoeff(:)
    real(dp), allocatable :: binMin(:), binMax(:)
    real(dp), allocatable :: binIndx(:)

    character(len=50) :: fileName = "umbrella.dat"

    contains
       procedure, pass :: MakeDecision => Umbrella_MakeDecision
       procedure, pass :: GetBiasIndex => Umbrella_GetBiasIndex
       procedure, pass :: Prologue => Umbrella_Prologue
!       procedure, pass :: Maintenance => Umbrella_Maintenance
!       procedure, pass :: ProcessIO => Umbrella_ProcessIO
  end type
!====================================================================
  contains
!====================================================================
  subroutine Umbrella_Constructor(self)
    implicit none
    class(Umbrella), intent(inout) :: self

    allocate( self%indexCoeff(1:self%nBiasVar) ) 
    allocate( self%binMax(1:self%nBiasVar) ) 
    allocate( self%binMin(1:self%nBiasVar) ) 
    allocate( self%binIndx(1:self%nBiasVar) ) 

  end subroutine
!====================================================================
  subroutine Umbrella_Prologue(self)
    implicit none
    class(Umbrella), intent(inout) :: self
    integer :: i,j

     ! Since the number of biasing variables is only known at run time, the bias matrix
     ! must be stored in a 1D array.  The index coefficient variables are used to emulate a N dimension matrix
     ! using the linear mapping equation of this form:
     ! U =  a1*x1 + a2*x2 + a3*x3 + .....
     ! Which maps a N dimension matrix onto a 1D array. 

    self%indexCoeff(1) = 1
    do i = 2, nBiasVariables 
      self%indexCoeff(i) = 1
      do j = 1, i-1
        self%indexCoeff(i) = self%indexCoeff(i) + self%indexCoeff(j) * (self%binMax(j) - self%binMin(j))
      enddo
    enddo      


  end subroutine
!====================================================================
  function Umbrella_MakeDecision(self, trialBox, E_Diff, inProb, disp) result(accept)
    use Template_SimBox, only: SimBox
    use RandomGen, only: grnd
    implicit none
    class(Umbrella), intent(in) :: self
    class(SimBox), intent(in) :: trialBox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(in) :: inProb
    real(dp), intent(in) :: E_Diff

    logical :: accept
    real(dp) :: biasE

    accept = .false.
    biasE = -trialBox%beta * E_Diff + log(inProb) 
    if(biasE > 0.0E0_dp) then
      accept = .true.
    elseif(biasE > log(grnd())) then
      accept = .true.
    endif

  end function
!==========================================================================
  function Umbrella_GetBiasIndex(self) result(biasIndx)
    implicit none
    class(Umbrella), intent(in) :: self
    integer :: biasIndx
    integer :: iBias
    
   ! Get the variables that are used for the biasing and figure out which histogram bin they
   ! fall into. 
    do iBias = 1, nBiasVariables

    enddo


   ! Using the bin values from each biasing variable, determine which
    biasIndx = 1
    do iBias = 1, nBiasVariables
!       write(*,*) biasIndx,self%indexCoeff(iBias), binIndx(iBias), binMin(iBias) 
      biasIndx = biasIndx + self%indexCoeff(iBias) * ( binIndx(iBias) - binMin(iBias) )
    enddo

!     write(*,*) biasIndx

  end function

!====================================================================
end module
!====================================================================
