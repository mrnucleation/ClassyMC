!==========================================================================================
module APeriodicOrthoDef
  use OrthoBoxDef, only: OrthoBox
  use SimpleSimBox, only: SimpleBox_ProcessIO
  use VarPrecision
  use CoordinateTypes

  ! Extension of OrthoBox that allows mixed periodic/non-periodic boundaries.
  ! Each dimension (X, Y, Z) can independently be set as periodic or non-periodic.
  ! For example, X and Y periodic but Z non-periodic for a slab geometry.
  !
  ! Input format in .clssy files:
  !   boxtype aperiodic_ortho
  !   dimension  Lx  Ly  Lz
  !   boundary   p   p   f
  !
  ! where 'p' = periodic, 'f' = fixed (non-periodic)
  ! Default: all dimensions are periodic (same as OrthoBox)

  type, public, extends(OrthoBox) :: APeriodicOrthoBox
    logical :: periodicX = .true.
    logical :: periodicY = .true.
    logical :: periodicZ = .true.
    contains
      procedure, pass :: Boundary => APeriodic_Boundary
      procedure, pass :: BoundaryNew => APeriodic_BoundaryNew
      procedure, pass :: DumpData => APeriodic_DumpData
      procedure, pass :: ProcessIO => APeriodic_ProcessIO
      procedure, pass :: LoadBoundaryFlags => APeriodic_LoadBoundaryFlags
  end type

!==========================================================================================
  contains
!==========================================================================================
  subroutine APeriodic_Boundary(self, rx, ry, rz)
    implicit none
    class(APeriodicOrthoBox), intent(in) :: self
    real(dp), intent(inout), optional :: rx, ry, rz

    ! Only apply periodic wrapping on dimensions marked as periodic
    if(present(rx)) then
      if(self%periodicX) then
        if(abs(rx) > self%boxLx2) then
          rx = rx - sign(self%boxLx, rx)
        endif
      endif
    endif

    if(present(ry)) then
      if(self%periodicY) then
        if(abs(ry) > self%boxLy2) then
          ry = ry - sign(self%boxLy, ry)
        endif
      endif
    endif

    if(present(rz)) then
      if(self%periodicZ) then
        if(abs(rz) > self%boxLz2) then
          rz = rz - sign(self%boxLz, rz)
        endif
      endif
    endif

  end subroutine
!==========================================================================================
  subroutine APeriodic_BoundaryNew(self, rx, ry, rz, disp)
    implicit none
    class(APeriodicOrthoBox), intent(in) :: self
    class(Perturbation), intent(in) :: disp(:)
    real(dp), intent(inout), optional :: rx, ry, rz
    real(dp) :: scaleX, scaleY, scaleZ

    select type(disp)
      class is(OrthoVolChange)
        scaleX = disp(1)%xScale
        scaleY = disp(1)%yScale
        scaleZ = disp(1)%zScale
    end select

    ! Only apply periodic wrapping on dimensions marked as periodic
    if(present(rx)) then
      if(self%periodicX) then
        if(abs(rx) > self%boxLx2*scaleX) then
          rx = rx - sign(self%boxLx*scaleX, rx)
        endif
      endif
    endif

    if(present(ry)) then
      if(self%periodicY) then
        if(abs(ry) > self%boxLy2*scaleY) then
          ry = ry - sign(self%boxLy*scaleY, ry)
        endif
      endif
    endif

    if(present(rz)) then
      if(self%periodicZ) then
        if(abs(rz) > self%boxLz2*scaleZ) then
          rz = rz - sign(self%boxLz*scaleZ, rz)
        endif
      endif
    endif

  end subroutine
!==========================================================================================
  subroutine APeriodic_LoadBoundaryFlags(self, line, lineStat)
    use Input_Format, only: GetXCommand
    use ParallelVar, only: nout
    implicit none
    class(APeriodicOrthoBox), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, flagX, flagY, flagZ

    lineStat = 0
    read(line, *) dummy, flagX, flagY, flagZ

    ! Parse X boundary
    select case(trim(adjustl(flagX)))
      case("p", "P", "periodic", "true", "t", "T", "1")
        self%periodicX = .true.
      case("f", "F", "fixed", "false", "n", "N", "0")
        self%periodicX = .false.
      case default
        write(0,*) "ERROR! Unrecognized boundary flag for X: ", trim(adjustl(flagX))
        write(0,*) "Use 'p' for periodic or 'f' for fixed (non-periodic)"
        lineStat = -1
        return
    end select

    ! Parse Y boundary
    select case(trim(adjustl(flagY)))
      case("p", "P", "periodic", "true", "t", "T", "1")
        self%periodicY = .true.
      case("f", "F", "fixed", "false", "n", "N", "0")
        self%periodicY = .false.
      case default
        write(0,*) "ERROR! Unrecognized boundary flag for Y: ", trim(adjustl(flagY))
        write(0,*) "Use 'p' for periodic or 'f' for fixed (non-periodic)"
        lineStat = -1
        return
    end select

    ! Parse Z boundary
    select case(trim(adjustl(flagZ)))
      case("p", "P", "periodic", "true", "t", "T", "1")
        self%periodicZ = .true.
      case("f", "F", "fixed", "false", "n", "N", "0")
        self%periodicZ = .false.
      case default
        write(0,*) "ERROR! Unrecognized boundary flag for Z: ", trim(adjustl(flagZ))
        write(0,*) "Use 'p' for periodic or 'f' for fixed (non-periodic)"
        lineStat = -1
        return
    end select

    write(nout,*) "Boundary conditions set:"
    write(nout,*) "  X: ", merge("periodic    ", "non-periodic", self%periodicX)
    write(nout,*) "  Y: ", merge("periodic    ", "non-periodic", self%periodicY)
    write(nout,*) "  Z: ", merge("periodic    ", "non-periodic", self%periodicZ)

  end subroutine
!==========================================================================================
  subroutine APeriodic_DumpData(self, filename)
    use Common_MolInfo, only: nMolTypes, MolData
    implicit none
    class(APeriodicOrthoBox), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: j, iType, iMol, iAtom, jType, subIndx, arrayIndx
    character(len=1) :: flagX, flagY, flagZ

    open(unit=50, file=trim(adjustl(filename)))

    write(50,*) "boxtype aperiodic_ortho"
    write(50,*) "dimension", self%boxLx, self%boxLy, self%boxLz

    flagX = merge("p", "f", self%periodicX)
    flagY = merge("p", "f", self%periodicY)
    flagZ = merge("p", "f", self%periodicZ)
    write(50,*) "boundary", " ", flagX, " ", flagY, " ", flagZ

    write(50,*) "molmin", (self%NMolMin(iType), iType=1,nMolTypes)
    write(50,*) "molmax", (self%NMolMax(iType), iType=1,nMolTypes)
    write(50,*) "mol", (self%NMol(iType), iType=1,nMolTypes)

    do iType = 1, nMolTypes
      do iMol = 1, self%NMol(iType)
        do iAtom = 1, MolData(iType)%nAtoms
          subIndx = 0
          do jType = 1, iType-1
            subIndx = subIndx + self%NMolMax(jType)
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
  subroutine APeriodic_ProcessIO(self, line, lineStat)
    use CoordinateTypes
    use Input_Format, only: maxLineLen, GetXCommand, LowerCaseLine
    use ForcefieldData, only: EnergyCalculator
    use Units, only: inPressUnit
    implicit none

    class(APeriodicOrthoBox), intent(inout) :: self
    integer, intent(out) :: lineStat
    character(len=maxLineLen), intent(in) :: line

    integer :: intVal
    real(dp) :: realVal
    character(len=30) :: command, val

    lineStat = 0
    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case default
        call SimpleBox_ProcessIO(self, line, lineStat)
    end select

  end subroutine
!==========================================================================================
end module
!==========================================================================================
