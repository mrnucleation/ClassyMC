!=========================================================================
#define __StdErr__ 0
!=========================================================================
module Analysis_RDF

use AnaylsisClassDef, only: Analysis
use SimpleSimBox, only: SimpleBox
use ClassyConstants, only: pi
use VarPrecision

  type, public, extends(Analysis) :: rdf
!    integer :: maintFreq = 100
!    integer :: IOUnit = -1
!    integer :: UpdateFreq = -1

    logical :: parallel = .false.
    integer :: boxNum = 1
    integer :: listIndx = 1

    integer :: thermNum = -1

    integer :: type1, type2
    integer :: bins = 10
    integer :: nSamples = 0
    real(dp) :: rMin, rMax, rMaxSq
    real(dp) :: dr = 0.01E0_dp


    integer :: fileUnit
    character(len=80) :: fileName = ""

    real(dp), allocatable :: hist(:)
    real(dp), allocatable :: tempHist(:)
    class(SimpleBox), pointer :: simbox => null()
    contains
      procedure, pass :: Prologue => RDF_Prologue
      procedure, pass :: Compute => RDF_Compute
      procedure, pass :: ProcessIO => RDF_ProcessIO
      procedure, pass :: ModifyIO => RDF_ModifyIO
      procedure, pass :: Maintenance => RDF_Maintenance
!      procedure, pass :: WriteInfo => RDF_WriteInfo
      procedure, pass :: CastCommonType => RDF_CastCommonType
  end type


 contains
!=========================================================================
  subroutine RDF_Prologue(self)
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
  subroutine RDF_Compute(self, accept)
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
    call self%simbox%GetNeighborList(self%listIndx, neighlist, nNeigh)

    do iAtom = 1, nMaxAtoms
      if( .not. self%simbox%IsActive(iAtom) ) then
        cycle
      endif

      atmType1 = self%simbox % AtomType(iAtom)
      if( (self%Type1 /= atmType1) .and. (self%Type2 /= atmType1)) then
        cycle
      endif
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)
        if(jAtom <= iAtom) then
          cycle
        endif
        atmType2 = self%simbox % AtomType(jAtom)
        if( (self%Type1 /= atmType2) .and. (self%Type2 /= atmType2)) then
          cycle
        endif

        rx = atoms(1, iAtom)  -  atoms(1, jAtom)
        ry = atoms(2, iAtom)  -  atoms(2, jAtom)
        rz = atoms(3, iAtom)  -  atoms(3, jAtom)
        call self%simbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rMaxSq) then
          r = sqrt(rsq)
          bin = floor( (r-self%rMin)*self%dR) + 1
          self%hist(bin) = self%hist(bin) + 2.0E0_dp

        endif
      enddo
    enddo

    self%nSamples = self%nSamples + 1

  end subroutine
!=========================================================================
  subroutine RDF_ProcessIO(self, line)
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
      write(tempStr, "(A)") "ERROR! The RDF module was expecting 9 arguments, but received %s."
      write(intStr, *) nPar-1
      tempStr = ReplaceText(tempStr, "%s", trim(adjustl(intStr)))
      write(__StdErr__, "(A)") tempStr
      write(__StdErr__, "(A)") trim(adjustl(line))
      write(__StdErr__, "(A)") "Format: RDF (BoxNumber) (UpdateFreq) (Write Freq) (Type 1) (Type 2) (rMin) (rMax) (nBins) (FileName)"
      error stop
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
      write(__StdErr__, *) "ERROR! The RDF module was given an r-max smaller than the corresponding r-min!"
      write(__StdErr__, *) "r-min:", self%rMin
      write(__StdErr__, *) "r-max:", self%rMax
      error stop
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
#ifdef MPIPARALLEL
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
  subroutine RDF_Maintenance(self)
    use ClassyConstants, only: pi
#ifdef MPIPARALLEL
    use MPI
#endif
    use ParallelVar, only: myid, nout, p_size
    implicit none
    class(rdf), intent(inout) :: self
    integer :: ierror
    integer :: iBin, nProc, nAtoms, nTotal
    integer :: nSamples
    real(dp) :: r, binVol, norm, volume

#ifdef MPIPARALLEL
    
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
  subroutine RDF_CastCommonType(self, anaVar)
    implicit none
    class(RDF), intent(inout) :: self
    class(*), allocatable, intent(inout) :: anaVar
    real(dp) :: def


    if(.not. allocated(anaVar) ) then
      allocate(anaVar, source=def)
!      write(*,*) "Allocated as Real"
    endif

  end subroutine
!====================================================================
  subroutine RDF_ModifyIO(self, line, lineStat)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText, CountCommands
    use ParallelVar, only: myid
    implicit none
    class(RDF), intent(inout) :: self
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
