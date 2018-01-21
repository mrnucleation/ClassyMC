!==========================================================================================
module OrthoBoxDef
  use SimpleSimBox, only: SimpleBox
!  use NeighListDef, only: NeighList
  use VarPrecision
!  use ForcefieldData, only: ECalcArray
  use CoordinateTypes, only: Displacement


  !Sim Box Definition
  type, public, extends(SimpleBox) :: OrthoBox
    real(dp) :: boxLx, boxLx2
    real(dp) :: boxLy, boxLy2
    real(dp) :: boxLz, boxLz2
    contains
!      procedure, pass :: Constructor => Ortho_Constructor
      procedure, pass :: LoadDimension => Ortho_LoadDimension
      procedure, pass :: UpdateEnergy => Ortho_UpdateEnergy
      procedure, pass :: Boundary => Ortho_Boundary
!      procedure, pass :: UpdatePosition => Ortho_UpdatePosition
      procedure, pass :: IOProcess => Ortho_IOProcess
      procedure, pass :: DumpData => Ortho_DumpData
  end type

!==========================================================================================
  contains
!==========================================================================================
  subroutine Ortho_LoadDimension(self, line, lineStat)
    use Input_Format, only: GetXCommand
    implicit none
    class(OrthoBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: dummy
    real(dp) :: boxLx, boxLy, boxLz

    lineStat = 0

    read(line, *) dummy, boxLx, boxLy, boxLz
    self%boxLx = boxLx
    self%boxLy = boxLy
    self%boxLz = boxLz
    self%boxLx2 = boxLx/2.0E0_dp
    self%boxLy2 = boxLy/2.0E0_dp
    self%boxLz2 = boxLz/2.0E0_dp

  end subroutine
!==========================================================================================
  subroutine Ortho_Boundary(self, rx, ry, rz)
  implicit none
  class(OrthoBox), intent(in) :: self
  real(dp), intent(inout) :: rx, ry, rz 
  real(dp) :: rx_new, ry_new, rz_new

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
  subroutine Ortho_UpdateEnergy(self, E_Diff)
    implicit none
    class(OrthoBox), intent(inout) :: self
    real(dp), intent(in) :: E_Diff

    self % ETotal = self % ETotal + E_Diff
    self % ETable = self % ETable + self % dETable

  end subroutine
!==========================================================================================
!  subroutine Ortho_UpdatePosition(self, disp)
!    use CoordinateTypes
!    implicit none
!    class(OrthoBox), intent(inout) :: self
!    type(Displacement), intent(inout) :: disp(:)
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
  subroutine Ortho_IOProcess(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    implicit none

    class(OrthoBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line   

    integer :: intVal
    real(dp) :: realVal
    character(len=30) :: command, val

!    write(*,*) "Cube"
    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
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

      case("neighcut")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        call GetXCommand(line, command, 6, lineStat)
        read(command, *) realVal
        self % NeighList(intVal) % rCut = realVal

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
    integer :: iType, iMol, iAtom, jType, subIndx, arrayIndx

    open(unit=50, file=trim(adjustl(filename)))

    write(50,*) "boxtype cube"
    write(50,*) "dimension", self%boxLx,  self%boxLy,  self%boxLz 
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
end module
!==========================================================================================
