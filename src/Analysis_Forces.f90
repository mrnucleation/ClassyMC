!=========================================================================
#define __StdErr__ 0
!=========================================================================
module Analysis_ComputeForces

use AnaylsisClassDef, only: Analysis
use SimpleSimBox, only: SimpleBox
use ClassyConstants, only: pi
use VarPrecision

  !=================
  !  Class Variables
  !  logical continous => Determines if the forces are updated every step or only periodically.
  !  int boxNum => BoxArray index for the target box number.
  !  SimpleBox simbox boxNum => Pointer to the box object 
  !================
  type, public, extends(Analysis) :: rdf
!    integer :: maintFreq = 100
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    logical :: continous = .false.
    integer :: boxNum = 1
    class(SimpleBox), pointer :: simbox => null()

    real(dp), allocatable :: forcecomponent(:,:)
    contains
      procedure, pass :: Prologue => ComputeForces_Prologue
      procedure, pass :: Compute => ComputeForces_Compute
      procedure, pass :: ProcessIO => ComputeForces_ProcessIO
      procedure, pass :: ModifyIO => ComputeForces_ModifyIO
      procedure, pass :: Maintenance => ComputeForces_Maintenance
!      procedure, pass :: WriteInfo => ComputeForces_WriteInfo
      procedure, pass :: CastCommonType => ComputeForces_CastCommonType
  end type


 contains
!=========================================================================
  subroutine ComputeForces_Prologue(self)
    implicit none
    class(rdf), intent(inout) :: self
    character(len=30), parameter :: volume = "volume"

    if(allocated(self%hist)) then
      deallocate(self%hist)
    endif
    allocate(self%hist(1:self%bins))

     if(allocated(self%temphist)) then
      deallocate(self%temphist)
    endif
    allocate(self%temphist(1:self%bins)) 

    self%hist = 0E0_dp
    self%temphist = 0E0_dp

    self%perMove = .false.
    self%maintFreq = 10*self%UpdateFreq

    self%thermNum = self%simbox%ThermoLookUp(volume)
  end subroutine
!=========================================================================
  subroutine ComputeForces_Compute(self, accept)
    use BoxData, only: BoxArray
    implicit none
    class(rdf), intent(inout) :: self
    logical, intent(in) :: accept
 
    real(dp), pointer :: atoms(:,:) => null()
    integer, pointer :: neighlist(:,:) => null()
    integer, pointer :: nNeigh(:) => null()

    integer :: iAtom, jAtom, jNei
    integer :: nMaxAtoms
    integer :: bin
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r

    call self%simbox%GetCoordinates(atoms)
    nMaxAtoms = self%simbox%GetMaxAtoms()

    do iAtom = 1, nMaxAtoms
      if( .not. self%simbox%IsActive(iAtom) ) then
        cycle
      endif

      atmType1 = self%simbox % AtomType(iAtom)
    enddo


  end subroutine
!=========================================================================
  subroutine ComputeForces_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText, CountCommands
    use ParallelVar, only: myid
    implicit none
    class(rdf), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer :: lineStat = 0
    integer :: intVal, iCharacter, nPar
    real(dp) :: realVal
    character(len=30) :: idString, command
    character(len=80) :: tempStr, intStr

    call CountCommands(line, nPar)
    if(nPar /= 10) then
      write(tempStr, "(A)") "ERROR! The ComputeForces module was expecting 9 arguments, but received %s."
      write(intStr, *) nPar-1
      tempStr = ReplaceText(tempStr, "%s", trim(adjustl(intStr)))
      write(__StdErr__, "(A)") tempStr
      write(__StdErr__, "(A)") trim(adjustl(line))
      write(__StdErr__, "(A)") "Format: ComputeForces (BoxNumber) (UpdateFreq) (Write Freq) (Type 1) (Type 2) (rMin) (rMax) (nBins) (FileName)"
      stop
    endif
    !Format =  BoxNum (UpdateFreq) (Write Freq) (Type 1) (Type 2) (rMin) (rMax) (nBins) (FileName)
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal

    self%boxNum = intVal
    self%simbox => BoxArray(self%boxnum)%box

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) intVal
    self%UpdateFreq = intVal

    call GetXCommand(line, command, 4, lineStat)
    read(command, *) intVal
    self%maintFreq = intVal

    call GetXCommand(line, command, 5, lineStat)
    read(command, *) intVal
    self%Type1 = intVal

    call GetXCommand(line, command, 6, lineStat)
    read(command, *) intVal
    self%Type2 = intVal

    call GetXCommand(line, command, 7, lineStat)
    read(command, *) realVal
    self%rMin = realVal

    call GetXCommand(line, command, 8, lineStat)
    read(command, *) realVal
    self%rMax = realVal
    self%rMaxSq = realVal**2

    if(self%rMax < self%rMin) then
      write(__StdErr__, *) "ERROR! The ComputeForces module was given an r-max smaller than the corresponding r-min!"
      write(__StdErr__, *) "r-min:", self%rMin
      write(__StdErr__, *) "r-max:", self%rMax
      stop
    endif

    self%bins = 0
    call GetXCommand(line, command, 9, lineStat)
    read(command, *) intVal
    self%bins = intVal


    self%dr = real(self%bins,dp)/(self%rMax - self%rMin)

    call GetXCommand(line, command, 10, lineStat)
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
  subroutine ComputeForces_Maintenance(self)
    use ClassyConstants, only: pi
#ifdef PARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout, p_size
    implicit none
    class(rdf), intent(inout) :: self
    integer :: ierror
    integer :: iBin, nProc, nAtoms, nTotal
    integer :: nSamples
    real(dp) :: r, binVol, norm, volume

#ifdef PARALLEL
    
    if(self%parallel) then
        write(nout, *) "Stopping for block averaging"
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
        call MPI_REDUCE(self%temphist, self%hist, self%bins, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 

        call MPI_REDUCE(nSamples, self%nSamples, 1, &
                    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror) 


    endif
#endif

    if(.not. self%parallel) then

      self%temphist = self%hist
      nSamples = self%nSamples
    endif

    if(self%parallel) then
      if(myid /= 0 ) then
        return
      endif
    endif

    if(self%Type1 == self%Type2) then
      nTotal = self%simbox%CountAtoms(self%Type1)
    else
      nAtoms = self%simbox%CountAtoms(self%Type1)
      nTotal = nAtoms
      nAtoms = self%simbox%CountAtoms(self%Type2)
      nTotal = nTotal + nAtoms
    endif

    rewind(self%fileunit)
    volume = self % simbox % GetThermo(self%thermNum)
!    norm = sum(self%temphist(1:self%bins))
!    write(*,*) nTotal, volume, nTotal/volume
    do iBin = 1, self%bins
      r = real(iBin-1, dp)/self%dR
      binVol = (4E0_dp/3E0_dp)*pi*((r+1E0_dp/self%dR)**3 - r**3)
      r = r + 0.5E0_dp/self%dR
      write(self%fileunit, *) r, self%temphist(iBin)/(nSamples*binVol*nTotal**2/volume)
    enddo

  end subroutine
!=========================================================================
  subroutine ComputeForces_CastCommonType(self, anaVar)
    implicit none
    class(ComputeForces), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
!      write(*,*) "Allocated as Real"
    endif

  end subroutine
!====================================================================
  subroutine ComputeForces_ModifyIO(self, line, lineStat)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText, CountCommands
    use ParallelVar, only: myid
    implicit none
    class(ComputeForces), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat
    integer :: intVal
    character(len=30) :: command

    !Format = Modify Analysis Num (Modification Type) (Values)
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%listIndx = intVal
      case default

        lineStat = -1
    end select


  end subroutine

!=========================================================================
end module
!=========================================================================
