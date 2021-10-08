!=========================================================================
! Computes the 
!=========================================================================
module Anaylsis_MolFractionHist
use AnaylsisClassDef, only: Analysis
use VarPrecision

  type, public, extends(Analysis):: MolFractionHist
!    logical :: perMove = .false.
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    logical, private :: parallel = .false.
    integer, private :: nBins = 100
    integer, private :: boxNum = -1
    integer, private :: molType = 1

    real(dp), allocatable, private :: hist(:)
    real(dp), allocatable, private :: nSample(:)
    real(dp), allocatable, private :: tempHist(:)
    real(dp), allocatable, private :: tempnSample(:)

    integer, private :: fileUnit
    character(len=80), private :: fileName = ""

    contains
!      procedure, pass :: Initialize
      procedure, pass :: Prologue => MolFractionHist_Prologue
      procedure, pass :: Compute => MolFractionHist_Compute
!      procedure, pass :: CalcNewState => MolFractionHist_CalcNewState
      procedure, pass :: CastCommonType => MolFractionHist_CastCommonType
      procedure, pass :: Maintenance =>  MolFractionHist_Maintenance
      procedure, pass :: ProcessIO => MolFractionHist_ProcessIO
      procedure, pass :: GetResult => MolFractionHist_GetResult
  end type

 contains
!=========================================================================
  subroutine MolFractionHist_Prologue(self)
    use BoxData, only: BoxArray
    implicit none
    class(MolFractionHist), intent(inout) :: self
    integer :: nMax

    self%UpdateFreq = 1
    self%perMove = .true.

    nMax = BoxArray(self%boxNum) % box % maxMol
    self%nBins = nMax
    if(.not. allocated(self%hist) ) then
      allocate(self%hist(1:nMax))
      allocate(self%nSample(1:nMax))
      allocate(self%temphist(1:nMax))
      allocate(self%tempnSample(1:nMax))
      self%hist = 0E0_dp
      self%temphist = 0E0_dp
      self%nSample = 0E0_dp
      self%tempnSample = 0E0_dp
    endif


  end subroutine
!=========================================================================
  subroutine MolFractionHist_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(MolFractionHist), intent(inout) :: self
    logical, intent(in) :: accept

    integer :: nMol
    integer :: nTotal
    real(dp) :: molFrac


    nMol = BoxArray(self%boxNum) % box % NMol(self%molType)
    nTotal = BoxArray(self%boxNum) % box % nMolTotal
    molFrac = real(nMol, dp)/real(nTotal, dp)
    self%hist(nTotal) = self%hist(nTotal) + molFrac
    self%nSample(nTotal) = self%nSample(nTotal) + 1E0_dp
  end subroutine
!=========================================================================
  function MolFractionHist_GetResult(self) result(var)
    implicit none
    class(MolFractionHist), intent(in) :: self
    real(dp) :: var

    var = 0E0_dp
    stop "MolFractionHist does not support variable look up. "
  end function
!=========================================================================
  subroutine MolFractionHist_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText
    use ParallelVar, only: myid, nout, p_size
    implicit none
    class(MolFractionHist), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: idString, command
    integer :: lineStat = 0
    integer :: intVal, iCharacter

    !Format = (MolFractionHist) (Box Number) (Mol Type) (Number of Bins) (Filename)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    self%boxNum = intVal

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) intVal
    self%molType = intVal

    call GetXCommand(line, command, 4, lineStat)
    read(command, *) intVal
    self%MaintFreq = intVal

    call GetXCommand(line, command, 5, lineStat)
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
  subroutine MolFractionHist_Maintenance(self)
    use ClassyConstants, only: pi
#ifdef PARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout, p_size
    implicit none
    class(MolFractionHist), intent(inout) :: self
    integer :: ierror
    integer :: iBin, nProc, nAtoms, nTotal
    integer :: nSamples
    real(dp) :: frac, norm, volume

#ifdef PARALLEL
    
    if(self%parallel) then
        write(nout, *) "Stopping for block averaging"
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
        call MPI_REDUCE(self%temphist, self%hist, self%nBins, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 
        call MPI_REDUCE(self%tempnSample, self%nSample, self%nBins, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

    endif
#endif

    if(.not. self%parallel) then
      self%temphist = self%hist
      self%tempnSample = self%nSample
    endif

    if(self%parallel) then
      if(myid /= 0 ) then
        return
      endif
    endif


    rewind(self%fileunit)
    do iBin = 1, self%nBins
      if(self%tempnSample(iBin) > 1e-10) then
        write(self%fileunit, *) iBin, self%temphist(iBin)/(real(self%tempnSample(iBin), dp) )
!      else
!        write(self%fileunit, *) iBin, 0.0E0_dp
      endif
    enddo
    flush(self%fileunit)

  end subroutine
!=========================================================================
  subroutine MolFractionHist_CastCommonType(self, anaVar)
    implicit none
    class(MolFractionHist), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def

    def = 0
    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
