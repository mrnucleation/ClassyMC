!=========================================================================
#define __StdErr__ 0
!=========================================================================
module Analysis_AngleDistribution

use AnaylsisClassDef, only: Analysis
use SimpleSimBox, only: SimpleBox
use ClassyConstants, only: pi
use VarPrecision

  type, public, extends(Analysis) :: angledistribution
!    integer :: maintFreq = 100
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1
    logical :: parallel = .false.
    integer :: boxNum = 1


    !molType => Molecule Type ID
    !angleNum => The angle index as it appears in the MolData angle array
    !mem1 => Atom Relative index for the first member of the angle
    !mem2 => Atom Relative index for the second member of the angle
    !bins => Number of histogram bins
    !nSamples => Total number of 
    !rMax => Upper bound of the histogram, values above this will be ignored
    !dr => Bin size
    integer, private :: molType, angleNum
    integer, private :: angleType, mem1, mem2, mem3
    integer, private :: bins = 10
    integer, private:: nSamples = 0
    real(dp), private :: dtheta = 0.01E0_dp
    real(dp), private :: temppos(1:3, 1:3)

    !fileUnit => Fortran File ID Unit
    !fileName => Name of the file to write the output to.
    integer, private :: fileUnit
    character(len=80), private :: fileName = ""

    !hist => Local Histogram for Angle Distance
    !temphist => Temporary Placeholder Histogram for Angle Distance 
    !            used when gathering data across MPI instances of the simulation.
        !simbox => Pointer to the Simulation Box that data is being collected from
    real(dp), allocatable, private :: hist(:)
    real(dp), allocatable, private :: tempHist(:)
    class(SimpleBox), pointer :: simbox => null()
    contains
      procedure, pass :: Prologue => AngleDistribution_Prologue
      procedure, pass :: Compute => AngleDistribution_Compute
      procedure, pass :: ProcessIO => AngleDistribution_ProcessIO
      procedure, pass :: Maintenance => AngleDistribution_Maintenance
!      procedure, pass :: WriteInfo => AngleDistribution_WriteInfo
      procedure, pass :: CastCommonType => AngleDistribution_CastCommonType
  end type


 contains
!=========================================================================
  subroutine AngleDistribution_Prologue(self)
    implicit none
    class(AngleDistribution), intent(inout) :: self

    if(allocated(self%hist)) then
      deallocate(self%hist)
    endif
    allocate(self%hist(0:self%bins))

     if(allocated(self%temphist)) then
      deallocate(self%temphist)
    endif
    allocate(self%temphist(0:self%bins)) 

    self%hist = 0E0_dp
    self%temphist = 0E0_dp

    self%perMove = .false.
!    self%maintFreq = 10*self%UpdateFreq

  end subroutine
!===============================================================
  subroutine AngleDistribution_Compute(self, accept)
    use Common_MolInfo, only: MolData, AngleData
    use ClassyConstants, only: pi
    implicit none
    class(AngleDistribution), intent(inout) :: self
    logical, intent(in) :: accept
 
    real(dp), pointer :: atoms(:,:) => null()

    integer :: iMol, nMolecules
    integer :: bin, molIndx
    integer :: molStart, molEnd
    integer :: slice(1:2)
    real(dp) :: theta

    nMolecules = self%simbox%NMol(self%molType)
    do iMol = 1, nMolecules
      molIndx = self%simbox%MolGlobalIndx(self%molType, iMol)
      call self%simbox%GetMolData(molIndx, molStart=molStart, molEnd=molEnd)
      slice(1) = molStart
      slice(2) = molEnd
      call self%simbox%GetCoordinates(atoms, slice=slice)

      self%temppos(1:3, 1) = atoms(1:3, self%mem1)
      self%temppos(1:3, 2) = atoms(1:3, self%mem2)
      self%temppos(1:3, 3) = atoms(1:3, self%mem3)
      theta = AngleData(self%angletype) % angleFF % ComputeAngle(self%simbox, self%temppos(1:3,1:3))
      bin = floor(theta/self%dtheta) 
      self%hist(bin) = self%hist(bin) + 1.0E0_dp
      self%nSamples = self%nSamples + 1
    enddo


  end subroutine
!=========================================================================
  subroutine AngleDistribution_ProcessIO(self, line)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText, CountCommands
    use ParallelVar, only: myid
    use ClassyConstants, only: pi
    implicit none
    class(AngleDistribution), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer :: lineStat = 0
    integer :: intVal, iCharacter, nPar
    real(dp) :: realVal
    character(len=30) :: idString, command
    character(len=80) :: tempStr, intStr

    call CountCommands(line, nPar)
    if(nPar /= 8) then
      write(tempStr, "(A)") "ERROR! The AngleDistribution module was expecting 8 arguments, but received %s."
      write(intStr, *) nPar-1
      tempStr = ReplaceText(tempStr, "%s", trim(adjustl(intStr)))
      write(__StdErr__, "(A)") tempStr
      write(__StdErr__, "(A)") trim(adjustl(line))
      write(__StdErr__, "(A)") "Format: AngleDistribution (UpdateFreq) (Write Freq) (MolType) (Angle Number) (dtheta) (FileName)"
      error stop
    endif
    !Format =  BoxNum (UpdateFreq) (Write Freq) (MolType) (Angle Number) (dtheta) (FileName)
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
    self%molType = intVal

    call GetXCommand(line, command, 6, lineStat)
    read(command, *) intVal
    self%angleNum = intVal
    self%angleType = MolData(self%molType)%angle(self%angleNum)%angleType
    self%mem1 = MolData(self%molType)%angle(self%angleNum)%mem1
    self%mem2 = MolData(self%molType)%angle(self%angleNum)%mem2
    self%mem3 = MolData(self%molType)%angle(self%angleNum)%mem3

    call GetXCommand(line, command, 7, lineStat)
    read(command, *) realVal
    self%dtheta = realVal
    self%bins = ceiling(pi/self%dtheta)

    call GetXCommand(line, command, 8, lineStat)
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
  subroutine AngleDistribution_Maintenance(self)
    use ClassyConstants, only: pi
#ifdef PARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout, p_size
    implicit none
    class(AngleDistribution), intent(inout) :: self
    integer :: ierror
    integer :: iBin, nProc, nAtoms, nTotal
    integer :: nSamples
    real(dp) :: theta, binVol, norm, volume

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


    rewind(self%fileunit)
    do iBin = 0, self%bins-1
      theta = real(iBin, dp)*self%dtheta
      write(self%fileunit, *) theta, self%temphist(iBin)/(real(nSamples, dp) * self%dtheta)
    enddo

  end subroutine
!=========================================================================
  subroutine AngleDistribution_CastCommonType(self, anaVar)
    implicit none
    class(AngleDistribution), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
!      write(*,*) "Allocated as Real"
    endif

  end subroutine
!=========================================================================
end module
!=========================================================================
