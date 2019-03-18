#define __StdErr__ 0
!=========================================================================
module Anaylsis_DensityOfStates
  use AnaylsisClassDef, only: Analysis
  use VarPrecision

  type, public, extends(Analysis):: DensityOfStates
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

#ifdef PARALLEL
    logical :: parallel = .true.
#else
    logical, parameter :: parallel = .false.
#endif

    integer :: fileunit = 300
    integer :: boxNum = -1
    integer :: thermNum = -1
    integer :: writeNum = 0
    character(len=80) :: fileName = ""

    integer :: nBins
    integer :: misses = 0
    real(dp) :: EMin, EMax, dE
    real(dp), allocatable :: hist(:)
    real(dp), allocatable :: temphist(:)
    contains
!      procedure, pass :: Initialize
      procedure, pass :: Compute => DensOfStates_Compute
      procedure, pass :: Maintenance => DensOfStates_Maintenance
      procedure, pass :: ProcessIO => DensOfStates_ProcessIO
!      procedure, pass :: WriteInfo => DensOfStates_WriteInfo
!      procedure, pass :: GetResult => DensOfStates_GetResult
!      procedure, pass :: Finalize => DensOfStates_Finalize
      procedure, pass :: CastCommonType => DensOfStates_CastCommonType
      procedure, pass :: Epilogue => DensOfStates_Epilogue
  end type

 contains
!=========================================================================
  subroutine DensOfStates_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(DensityOfStates), intent(inout) :: self
    logical, intent(in) :: accept
    integer :: EBin
    real(dp) :: thermVal

    thermVal = BoxArray(self%boxNum) % box % GetThermo(self%thermNum)

    if(thermVal <= self%EMin) then
      self%misses = self%misses + 1
      return
    else if(thermVal >= self%EMax) then
      self%misses = self%misses + 1
      return
    endif
    EBin = floor( (thermVal-self%EMin)*self%dE ) + 1
    self%temphist(EBin) = self%temphist(EBin) + 1E0_dp

  end subroutine
!=========================================================================
  subroutine DensOfStates_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText
    use ParallelVar, only: myid, nout
    use Units, only: outEngUnit
    implicit none
    class(DensityOfStates), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: lineStat = 0
    integer :: intVal, iCharacter
    real(dp) :: realVal
    character(len=30) :: idString
    character(len=30), parameter :: energy = "energy"


    !Input format
    ! DensityOfStates (BoxNum) (Energy Min) (Energy Max) (nBins) (FileName)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal
    self%thermNum = BoxArray(self%boxNum) % box % ThermoLookUp(energy)

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) realVal
    self%EMin = realVal*outEngUnit

 
    call GetXCommand(line, command, 4, lineStat)
    read(command, *) realVal
    self%EMax = realVal*outEngUnit

    if(self%EMax <= self%EMin) then
      write(__StdErr__, *) "WARNING! Density of state's energy bounds are improperly set!"
      write(__StdErr__, *) "E_Max is less than E_Min!", self%EMax, self%EMin
      stop
    endif

    call GetXCommand(line, command, 5, lineStat)
    read(command, *) intVal
    self%nBins = intVal

    allocate(self%temphist(1:self%nBins) )
    self%temphist = 0E0_dp
#ifdef PARALLEL
    if(myid == 0) then
      allocate(self%hist(1:self%nBins) )
      self%hist = 0E0_dp
    endif
#endif


    !dE is stored as a reciprocal
    self%dE = real(self%nBins, dp)/(self%EMax-self%EMin)

    call GetXCommand(line, command, 6, lineStat)
    read(command, *) self%fileName
    do iCharacter = 1, len(self%filename)
      if(self%filename(iCharacter:iCharacter) == "&") then
        write(idString, *) myid
        self%filename = ReplaceText(self%filename, "&", trim(adjustl(idString)))
#ifdef PARALLEL
        self%parallel = .false.
#endif
        exit
      endif
    enddo

    do iCharacter = 1, len(self%filename)
      if(self%filename(iCharacter:iCharacter) == '"') then
        self%fileName(iCharacter:iCharacter) = " "
      endif
    enddo


    open(newunit=self%fileunit, file=trim(adjustl(self%fileName)) )

  end subroutine
!=========================================================================
  subroutine DensOfStates_Maintenance(self)
#ifdef PARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout
    use Units, only: outEngUnit
    implicit none
    integer :: ierror
    integer :: iBin
    class(DensityOfStates), intent(inout) :: self
    real(dp) :: eVal


    rewind(self%fileunit)
#ifdef PARALLEL
    
    if(self%parallel) then
        write(nout, *) "Stopping for density of states averaging"
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)       

        call MPI_REDUCE(self%hist, self%temphist, self%nBins, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

      if(myid == 0) then
        do iBin = 1, self%nBins
          eVal = ((iBin-1)/self%dE) + self%EMin
          write(self%fileunit, *) eVal/outEngUnit, self%hist(iBin)
        enddo
      endif

    endif
#endif

    do iBin = 1, self%nBins
      eVal = ((iBin-1)/self%dE) + self%EMin
      write(self%fileunit, *) eVal/outEngUnit, self%temphist(iBin)
    enddo
  end subroutine
!=========================================================================
!  subroutine DensOfStates_WriteInfo(self)
!    use ParallelVar, only: nout
!    implicit none
!    class(DensityOfStates), intent(inout) :: self
!
!    write(nout, *) trim(adjustl(self%varName)), " Average: ", self%GetResult(), "+\-", self%GetSTDev()
!  end subroutine
!=========================================================================
  subroutine DensOfStates_CastCommonType(self, anaVar)
    implicit none
    class(DensityOfStates), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
!      write(*,*) "Allocated as Real"
    endif

  end subroutine
!=========================================================================
  subroutine DensOfStates_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(DensityOfStates), intent(inout) :: self

    write(nout, *) "Density of States Bin Misses:", self%misses

  end subroutine
!=========================================================================
end module
!=========================================================================
