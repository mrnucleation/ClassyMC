!==========================================================================================
module OrthoBoxDef
  use SimpleSimBox, only: SimpleBox
  use VarPrecision
  use CoordinateTypes


  !Sim Box Definition
  type, public, extends(SimpleBox) :: OrthoBox
    real(dp) :: boxLx, boxLx2
    real(dp) :: boxLy, boxLy2
    real(dp) :: boxLz, boxLz2
    contains
      procedure, pass :: Prologue => Ortho_Prologue
      procedure, pass :: Epilogue => Ortho_Epilogue
      procedure, pass :: LoadDimension => Ortho_LoadDimension
      procedure, pass :: GetDimensions => Ortho_GetDimensions
      procedure, pass :: Boundary => Ortho_Boundary
      procedure, pass :: BoundaryNew => Ortho_BoundaryNew
      procedure, pass :: ProcessIO => Ortho_ProcessIO
      procedure, pass :: DumpData => Ortho_DumpData
      procedure, pass :: UpdateVolume => Ortho_UpdateVolume

      procedure, pass :: GetRealCoords => Ortho_GetRealCoords
      procedure, pass :: GetReducedCoords => Ortho_GetReducedCoords
  end type

!==========================================================================================
  contains
!==========================================================================================
  subroutine Ortho_Prologue(self)
    use Common_MolInfo, only: nMolTypes, MolData
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(OrthoBox), intent(inout) :: self
    logical :: accept
    integer :: iType, iMol, iAtom, jType, subIndx, arrayIndx
    integer :: iConstrain, molStart, molEnd

    call self % ComputeEnergy
    call self % NeighList(1) % BuildList


    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Initial Energy: ", self % ETotal/outEngUnit
    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol


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

          call self%Boundary(self%atoms(1,arrayIndx), self%atoms(2,arrayIndx), self%atoms(3,arrayIndx))
        enddo
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
    write(nout, "(1x,A,I2,A,8I10)") "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout, "(1x,A,I2,A,I10)") "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal



    self%volume = self%boxLx * self%boxLy * self%boxLz
    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Pressure: ", self%pressure
    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Dimensions: ", self%volume
    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Volume: ", self%volume
    write(nout, "(1x,A,I2,A,E15.8)") "Box ", self%boxID, " Number Density: ", self%nMolTotal/self%volume


    write(nout, "(1x,A,I2,A,E15.8,1x,A)") "Box ", self%boxID, " Initial Energy: ", self % ETotal/outEngUnit, engStr
    write(nout, "(1x,A,I2,A,E15.8,1x,A)") "Box ", self%boxID, " Initial Energy (Per Mol): ", self % ETotal/(outEngUnit*self%nMolTotal), engStr
    do iMol = 1, self%maxMol
      call self%GetMolData(iMol, molStart=molStart)
      if( self%MolSubIndx(molStart) <= self%NMol(self%MolType(molStart)) ) then
        call self%ComputeCM(iMol)
      endif
    enddo


  end subroutine
!==========================================================================================
  subroutine Ortho_Epilogue(self)
    use ParallelVar, only: nout
    use Units, only: outEngUnit, engStr
    implicit none
    class(OrthoBox), intent(inout) :: self
    integer :: iConstrain
    real(dp) :: E_Culm

    E_Culm = self%ETotal

    write(nout,*) "--------Box", self%boxID , "Energy---------"
    call self % ComputeEnergy
    call self % NeighList(1) % BuildList

    write(nout, *) "Final Energy:", self % ETotal/outEngUnit, engStr
    write(nout, *) "Final Energy (Per Mol):", self % ETotal/(self%nMolTotal*outEngUnit), engStr
    if(self%ETotal /= 0) then
      if( abs((E_Culm-self%ETotal)/self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit, engStr
      endif
    else
      if( abs(E_Culm-self%ETotal) > 1E-7_dp ) then
        write(nout, *) "ERROR! Large energy drift detected!"
        write(nout, *) "Box: ", self%boxID
        write(nout, *) "Culmative Energy: ", E_Culm/outEngUnit, engStr
        write(nout, *) "Final Energy: ", self%ETotal/outEngUnit, engStr
        write(nout, *) "Difference: ", (self%ETotal-E_Culm)/outEngUnit, engStr
      endif
    endif

    write(nout,*) "Box ", self%boxID, " Molecule Count: ", self % NMol
    write(nout,*) "Box ", self%boxID, " Total Molecule Count: ", self % nMolTotal
    self%volume = self%boxLx * self%boxLy * self%boxLz
    write(nout,*) "Box ", self%boxID, " Lengths: ", self%boxLx,  self%boxLy, self%boxLz
    write(nout,*) "Box ", self%boxID, " Pressure: ", self%pressure
    write(nout,*) "Box ", self%boxID, " Volume: ", self%volume
    write(nout,*) "Box ", self%boxID, " Number Density: ", self%nMolTotal/self%volume
    if( size(self%Constrain) > 0 ) then
      do iConstrain = 1, size(self%Constrain)
        call self%Constrain(iConstrain) % method % Epilogue
      enddo
    endif


  end subroutine

!==========================================================================================
  subroutine Ortho_LoadDimension(self, line, lineStat)
    use Input_Format, only: GetXCommand
    use ParallelVar, only: nout
    implicit none
    class(OrthoBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: dummy
    real(dp) :: boxLx, boxLy, boxLz

    lineStat = 0

    read(line, *) dummy, boxLx, boxLy, boxLz
    write(nout,*) dummy, boxLx, boxLy, boxLz
    self%boxLx = boxLx
    self%boxLy = boxLy
    self%boxLz = boxLz

    self%boxLx2 = boxLx * 0.5E0_dp
    self%boxLy2 = boxLy * 0.5E0_dp
    self%boxLz2 = boxLz * 0.5E0_dp


  end subroutine
!==========================================================================================
  subroutine Ortho_GetDimensions(self, list)
    use Input_Format, only: GetXCommand
    implicit none
    class(OrthoBox), intent(inout) :: self
    real(dp), intent(inout) :: list(:, :)

    integer :: iDimen

    list = 0E0_dp
    list(1, 1) = -self%boxLx2
    list(2, 1) = self%boxLx2

    list(1, 2) = -self%boxLy2
    list(2, 2) = self%boxLy2

    list(1, 3) = -self%boxLz2
    list(2, 3) = self%boxLz2
  end subroutine
!==========================================================================================
  subroutine Ortho_Boundary(self, rx, ry, rz)
    implicit none
    class(OrthoBox), intent(in) :: self
    real(dp), intent(inout) :: rx, ry, rz 

    if(abs(rx) > self%boxLx2) then
      rx = rx - sign(self%boxLx, rx)
    endif
 
    if(abs(ry) > self%boxLy2) then
      ry = ry - sign(self%boxLy, ry)
    endif

    if(abs(rz) > self%boxLz2) then
      rz = rz - sign(self%boxLz, rz)
    endif

  end subroutine
!==========================================================================================
  subroutine Ortho_BoundaryNew(self, rx, ry, rz, disp)
  implicit none
  class(OrthoBox), intent(in) :: self
  class(Perturbation), intent(in) :: disp(:)
  real(dp), intent(inout) :: rx, ry, rz 
  real(dp) :: scaleX, scaleY, scaleZ

  select type(disp)
    class is(OrthoVolChange)
      scaleX = disp(1)%xScale
      scaleY = disp(1)%yScale
      scaleZ = disp(1)%zScale
  end select
  
  if(abs(rx) > self%boxLx2*scaleX) then
    rx = rx - sign(self%boxLx*scaleX, rx)
  endif

  if(abs(ry) > self%boxLy2*scaleY) then
    ry = ry - sign(self%boxLy*scaleY, ry)
  endif

  if(abs(rz) > self%boxLz2*scaleZ) then
    rz = rz - sign(self%boxLz*scaleZ, rz)
  endif

  end subroutine
!==========================================================================================
!  subroutine Ortho_UpdatePosition(self, disp)
!    use CoordinateTypes
!    implicit none
!    class(OrthoBox), intent(inout) :: self
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
  subroutine Ortho_ProcessIO(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    use Units, only: inPressUnit
    implicit none

    class(OrthoBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: intVal
    real(dp) :: realVal
    character(len=30) :: command, val

!    write(*,*) "Ortho"
    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("buildfreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self % maintFreq = intVal
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

      case("temperature")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % temperature = realVal
        self % beta = 1E0_dp/realVal

      case("chempot")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self % chempot(intVal) = realVal

      case("pressure")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self % pressure = realVal*inPressUnit

      case("neighcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self % NeighList(intVal) % rCut = realVal

      case("neighlist")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call self%NeighList(intVal)%ProcessIO(line, lineStat)

      case default
        lineStat = -1
    end select

  end subroutine
!==========================================================================================
  subroutine Ortho_DumpData(self, filename)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(OrthoBox), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: j, iType, iMol, iAtom, jType, subIndx, arrayIndx

    open(unit=50, file=trim(adjustl(filename)))

    write(50,*) "boxtype ortho"
    write(50,*) "dimension", self%boxLx, self%boxLy, self%boxLz
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

          write(50,*) iType, iMol, iAtom, (self%atoms(j,arrayIndx), j=1,3)                     
        enddo
      enddo
    enddo


    close(50)

  end subroutine
!==========================================================================================
  subroutine Ortho_UpdateVolume(self, disp)
    implicit none
    class(OrthoBox), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)

    select type(disp)
      class is(OrthoVolChange)
        self%volume = disp(1)%volNew
        self%boxLx = self%boxLx*disp(1)%xScale
        self%boxLx2 = self%boxLx2*disp(1)%xScale

        self%boxLy = self%boxLy*disp(1)%yScale
        self%boxLy2 = self%boxLy2*disp(1)%yScale

        self%boxLz = self%boxLz*disp(1)%zScale
        self%boxLz2 = self%boxLz2*disp(1)%zScale
    end select
  end subroutine
!==========================================================================================
  subroutine Ortho_GetReducedCoords(self,realCoords,reducedCoords )
    implicit none
    class(OrthoBox), intent(inout) :: self
    real(dp), intent(in) :: realCoords(:)
    real(dp), intent(out) :: reducedCoords(1:3)
    integer :: i

    reducedCoords(1) = (realCoords(1)+self%boxLx2)/self%boxLx
    reducedCoords(2) = (realCoords(2)+self%boxLy2)/self%boxLy
    reducedCoords(3) = (realCoords(3)+self%boxLz2)/self%boxLz


  end subroutine
!==========================================================================================
  subroutine Ortho_GetRealCoords(self, reducedCoords, realCoords)
    implicit none
    class(OrthoBox), intent(inout) :: self
    real(dp), intent(in) :: reducedCoords(:)
    real(dp), intent(out) :: realCoords(1:3)
    integer :: i

    realCoords(1) = self%boxLx*reducedCoords(1) - self%boxLx2
    realCoords(2) = self%boxLy*reducedCoords(2) - self%boxLy2
    realCoords(3) = self%boxLz*reducedCoords(3) - self%boxLz2


  end subroutine
!==========================================================================================
end module
!==========================================================================================
