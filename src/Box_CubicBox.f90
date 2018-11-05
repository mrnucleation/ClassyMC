!==========================================================================================
module CubicBoxDef
  use SimpleSimBox, only: SimpleBox
  use VarPrecision
  use CoordinateTypes

  !Sim Box Definition
  type, public, extends(SimpleBox) :: CubeBox
    real(dp) :: boxL, boxL2
    contains
!      procedure, pass :: Constructor => Cube_Constructor
      procedure, pass :: LoadDimension => Cube_LoadDimension
      procedure, pass :: GetDimensions => Cube_GetDimensions
      procedure, pass :: Boundary => Cube_Boundary
      procedure, pass :: BoundaryNew => Cube_BoundaryNew
      procedure, pass :: ProcessIO => Cube_ProcessIO
      procedure, pass :: DumpData => Cube_DumpData
      procedure, pass :: Prologue => Cube_Prologue
      procedure, pass :: UpdateVolume => Cube_UpdateVolume
  end type

!==========================================================================================
  contains
!==========================================================================================
!  subroutine Cube_Constructor(self)
!    use Common_MolInfo
!    implicit none
!    class(CubeBox), intent(inout) :: self
!    integer :: AllocateStatus
!    integer :: iType, iMol, iAtom, atmIndx, molIndx, maxMol
! 
!    if( .not. allocated(self%NMolMin) ) then
!      write(*,*) "ERROR! The maximum and minimum molecules allowed in the box must be defined"
!      write(*,*) "prior to box initialization!"
!      stop 
!    endif
!
!    !First begin by computing the maximium number of atoms that the box can potentially contain
!    self%nMaxAtoms = 0
!    maxMol = 0
!   do iType = 1, nMolTypes
!      self%nMaxAtoms = self%nMaxAtoms + self%NMolMax(iType)*MolData(iType)%nAtoms
!      maxMol = maxMol + self%NMolMax(iType)
!    enddo
!
!    !Allocate the position and energy related arrays. 
!    allocate(self%atoms(1:3, 1:self%nMaxAtoms), stat=AllocateStatus)
!    allocate(self%ETable(1:self%nMaxAtoms), stat=AllocateStatus)
!    allocate(self%dETable(1:self%nMaxAtoms), stat=AllocateStatus)
!
!    !Allocate the arrays which contain the atom type and quick look up information.
!    allocate(self%AtomType(1:self%nMaxAtoms), stat=AllocateStatus)
!    allocate(self%MolIndx(1:self%nMaxAtoms), stat=AllocateStatus)
!    allocate(self%SubIndx(1:self%nMaxAtoms), stat=AllocateStatus)
!    allocate(self%MolStartIndx(1:self%nMaxAtoms), stat=AllocateStatus)
!
!    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"!
!
!    self%AtomType = 0
!    self%MolIndx = 0
!    self%SubIndx = 0
!    self%MolStartIndx = 0
!
!    atmIndx = 0
!    molIndx = 0
!    do iType = 1, nMolTypes
!      do iMol = 1, self%NMolMax(iType)
!        molIndx = molIndx + 1
!        self%MolStartIndx(molIndx) = atmIndx + 1
!        do iAtom = 1, MolData(iType)%nAtoms
!          atmIndx = atmIndx + 1
!          self%AtomType(atmIndx) = MolData(iType)%atomType(iAtom)
!          self%MolIndx(atmIndx)  = molIndx
!          self%SubIndx(atmIndx)  = iAtom
!        enddo
!      enddo 
!    enddo
!
!  end subroutine
!==========================================================================================
  subroutine Cube_LoadDimension(self, line, lineStat)
    use Input_Format, only: GetXCommand
    implicit none
    class(CubeBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: dummy
    real(dp) :: boxL

    lineStat = 0
    read(line, *) dummy, boxL
    self%boxL = boxL
    self%boxL2 = boxL/2.0E0_dp

  end subroutine
!==========================================================================================
  subroutine Cube_GetDimensions(self, list)
    use Input_Format, only: GetXCommand
    implicit none
    class(CubeBox), intent(inout) :: self
    real(dp), intent(out) :: list(:, :)

    integer :: iDimen

    do iDimen = 1, self%nDimension
      list(1, iDimen) = -self%boxL2
      list(2, iDimen) = self%boxL2
    enddo


  end subroutine
!==========================================================================================
  subroutine Cube_Boundary(self, rx, ry, rz)
  implicit none
  class(CubeBox), intent(in) :: self
  real(dp), intent(inout) :: rx, ry, rz 
!  real(dp) :: rx_new, ry_new, rz_new

  if(abs(rx) > self%boxL2) then
    rx = rx - sign(self%boxL, rx)
  endif

  if(abs(ry) > self%boxL2) then
    ry = ry - sign(self%boxL, ry)
  endif

  if(abs(rz) > self%boxL2) then
    rz = rz - sign(self%boxL, rz)
  endif

  end subroutine
!==========================================================================================
  subroutine Cube_BoundaryNew(self, rx, ry, rz, disp)
  implicit none
  class(CubeBox), intent(in) :: self
  class(Perturbation), intent(in) :: disp(:)
  real(dp), intent(inout) :: rx, ry, rz 
  real(dp) :: scaleFactor

  select type(disp)
    class is(OrthoVolChange)
      scaleFactor = disp(1)%xScale
  end select
  
  if(abs(rx) > self%boxL2*scaleFactor) then
    rx = rx - sign(self%boxL*scaleFactor, rx)
  endif

  if(abs(ry) > self%boxL2*scaleFactor) then
    ry = ry - sign(self%boxL*scaleFactor, ry)
  endif

  if(abs(rz) > self%boxL2*scaleFactor) then
    rz = rz - sign(self%boxL*scaleFactor, rz)
  endif

  end subroutine

!==========================================================================================
!  subroutine Cube_UpdatePosition(self, disp)
!    use CoordinateTypes
!    implicit none
!    class(CubeBox), intent(inout) :: self
!    integer :: iDisp, dispLen, dispIndx
!
!    dispLen = size(disp)
!    do iDisp = 1, dispLen
!      dispIndx = disp(iDisp) % atmIndx
!      call self%Boundary( disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new )
!      self % atoms(1, dispIndx) = disp(iDisp)%x_new
!      self % atoms(2, dispIndx) = disp(iDisp)%y_new
!      self % atoms(3, dispIndx) = disp(iDisp)%z_new
!    enddo
!
!  end subroutine
!==========================================================================================
  subroutine Cube_ProcessIO(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    use ParallelVar, only: nout
    implicit none

    class(CubeBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: intVal
    real(dp) :: realVal
    character(len=30) :: command, val

!    write(*,*) "Cube"
    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("boxlength")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % boxL = realVal
        self % boxL2 = realVal/2.0E0_dp

      case("buildfreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % maintFreq = intVal
        write(nout,*) "Neigh Update Frequency:", self % maintFreq



      case("chempotential")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self % chempot(intVal) = realVal

      case("energycalc")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % EFunc => EnergyCalculator(intVal)

      case("neighcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self%NeighList(intVal)%rCut = realVal

      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call self%NeighList(intVal)%ProcessIO(line, lineStat)

      case("pressure")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % pressure = realVal

      case("temperature")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % temperature = realVal
        self % beta = 1E0_dp/realVal

      case default
        lineStat = -1
    end select

  end subroutine
!==========================================================================================
  subroutine Cube_DumpData(self, filename)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(CubeBox), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: iType, iMol, iAtom, jType, subIndx, arrayIndx

    open(unit=50, file=trim(adjustl(filename)))

    write(50,*) "boxtype cube"
    write(50,*) "dimension", self%boxL
    write(50,*) "molmin", (self%NMolMin(iType), iType=1,nMolTypes)
    write(50,*) "molmax", (self%NMolMax(iType), iType=1,nMolTypes)
    write(50,*) "mol", (self%NMol(iType), iType=1,nMolTypes)

    do iType = 1, nMolTypes
      do iMol = 1, self%NMol(iType)
        do iAtom = 1, MolData(iType)%nAtoms
          subIndx = 0
          do jType = 1, iType-1
            subIndx = self%NMolMax(jType)
          enddo
          subIndx = subIndx + iMol
          arrayIndx = self%MolStartIndx(subIndx)
          arrayIndx = arrayIndx + iAtom - 1

          write(50,*) iType, iMol, iAtom, self%atoms(1,arrayIndx), &
                       self%atoms(2,arrayIndx), self%atoms(3,arrayIndx)                     
        enddo
      enddo
    enddo


    close(50)

  end subroutine
!==========================================================================================
  subroutine Cube_Prologue(self)
    use Common_MolInfo, only: nMolTypes
    use ParallelVar, only: nout
    implicit none
    class(CubeBox), intent(inout) :: self
    logical :: accept
    integer :: iAtom, iDims, iConstrain, iType
    integer :: iMol, molStart

    call self % ComputeEnergy
    call self % NeighList(1) % BuildList

    do iAtom = 1, self%nMaxAtoms
      if( self%MolSubIndx(iAtom) > self%NMol(self%MolType(iAtom)) ) then
        cycle
      endif
      do iDims = 1, self%nDimension
        if(abs(self%atoms(iDims, iAtom)) > self%boxL2) then
          write(nout, *) "Warning! Particle out of bounds!"
          write(nout, *) "Particle Number:", iAtom
          write(nout, *) "Box Length:", self%boxL2
          write(nout, *) self%atoms(:, iAtom)
          stop
        endif
      enddo
    enddo
    
    if( size(self%Constrain) > 0 ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % Prologue
        call self%Constrain(iConstrain) % method % CheckInitialConstraint(self, accept)
      enddo
      if(.not. accept) then
        write(nout,*) "Initial Constraints are not statisfied!"
        stop
      endif
    endif

    self%nMolTotal = 0
    do iType = 1, nMolTypes    
      self%nMolTotal = self%nMolTotal + self % NMol(iType)
    enddo
    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout,*) "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal



    self%volume = self%boxL**3
    write(nout,*) "Box ", self%boxID, " Pressure: ", self%pressure
    write(nout,*) "Box ", self%boxID, " Volume: ", self%volume
    write(nout,*) "Box ", self%boxID, " Number Density: ", self%nMolTotal/self%boxL**3


    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Initial Energy: ", self % ETotal

    do iMol = 1, self%maxMol
      call self%GetMolData(iMol, molStart=molStart)
      if( self%MolSubIndx(molStart) <= self%NMol(self%MolType(molStart)) ) then
        call self%ComputeCM(iMol)
      endif
    enddo

  end subroutine
!==========================================================================================
  subroutine CubeBox_Epilogue(self)
    use ParallelVar, only: nout
    use Units, only: outEngUnit
    implicit none
    class(CubeBox), intent(inout) :: self
    integer :: iConstrain
    real(dp) :: E_Culm

    E_Culm = self%ETotal

    write(nout,*) "--------Box", self%boxID , "Energy---------"
    call self % ComputeEnergy
    call self % NeighList(1) % BuildList

    write(nout, *) "Final Energy:", self % ETotal/outEngUnit
    if(self%ETotal /= 0) then
      if( abs((E_Culm-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit
      endif
    else
      if( abs(E_Culm-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit
      endif
    endif

    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout,*) "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal
    self%volume = self%boxL**3
    write(nout,*) "Box ", self%boxID, " Pressure: ", self%pressure
    write(nout,*) "Box ", self%boxID, " Volume: ", self%volume
    write(nout,*) "Box ", self%boxID, " Number Density: ", self%nMolTotal/self%boxL**3
    if( size(self%Constrain) > 0 ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % Epilogue
      enddo
    endif


  end subroutine
!==========================================================================================
  subroutine Cube_UpdateVolume(self, disp)
    implicit none
    class(CubeBox), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)

    select type(disp)
      class is(OrthoVolChange)
        self%volume = disp(1)%volNew
        self%boxL = self%boxL*(disp(1)%volNew/disp(1)%volOld)**(1E0_dp/3E0_dp)
        self%boxL2 = self%boxL/2E0_dp
    end select
  end subroutine
!==========================================================================================

end module
!==========================================================================================
