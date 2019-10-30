!====================================================================
module Traj_POSCAR
  use VarPrecision
  use TrajectoryTemplate, only: trajectory

  type, public, extends(trajectory) :: TrajPOSCAR
    logical :: recenter = .false.
!    integer :: fileUnit = -1
!    integer :: boxNum = -1
!    integer :: outFreq = 5000
!    character(len=50) :: fileName
    integer :: xLen
    real(dp), allocatable :: boxDim(:,:)
    contains
       procedure, pass :: WriteFrame => TrajPOSCAR_WriteFrame
!       procedure, pass :: SetUnit
!       procedure, pass :: SetBox
!       procedure, pass :: SetFileName
!       procedure, pass :: SetFreq
!       procedure, pass :: OpenFile
!       procedure, pass :: CloseFile
       procedure, pass :: Prologue => TrajPOSCAR_Prologue
       procedure, pass :: Epilogue => TrajPOSCAR_Epilogue

  end type
!====================================================================
  contains
!====================================================================
  subroutine TrajPOSCAR_WriteFrame(self, iCycle) 
    use BoxData, only: BoxArray
    use CubicBoxDef, only: CubeBox
    use SimpleSimBox, only: SimpleBox
    use OrthoBoxDef, only: OrthoBox
    use Common_MolInfo, only: AtomData, nAtomTypes
    implicit none
    class(TrajPOSCAR), intent(inout) :: self
    integer(kind=8), intent(in) :: iCycle
    integer :: iAtom, jDim, boxNum, i, j
    integer :: nDim, atomType, molType
    real(dp) :: cm(1:3)
    real(dp) :: Lx, Ly, Lz
    character(len=5) :: atmSymb
    class(SimpleBox), pointer :: box



    boxNum = self%boxNum
    box => BoxArray(boxNum)%box
    nDim = box%nDimension

    write(self%fileUnit, "(A)") "Test"
    write(self%fileUnit, *) "1.0"
    select type(box)
      class is(CubeBox)
        call box%GetDimensions(self%boxDim(1:2, 1:3))
        cm(1) = self%boxDim(1,1)
        cm(2) = self%boxDim(1,2)
        cm(3) = self%boxDim(1,3)
        Lx = self%boxDim(2,1) - self%boxDim(1,1)
        Ly = self%boxDim(2,2) - self%boxDim(1,2)
        Lz = self%boxDim(2,3) - self%boxDim(1,3)

        write(self%fileUnit, *) Lx, 0.0, 0.0
        write(self%fileUnit, *) 0.0, Ly, 0.0
        write(self%fileUnit, *) 0.0, 0.0, Lz
        


      class is(OrthoBox)
        call box%GetDimensions(self%boxDim(1:2, 1:3))
        cm(1) = self%boxDim(1,1)
        cm(2) = self%boxDim(1,2)
        cm(3) = self%boxDim(1,3)
        Lx = self%boxDim(2,1) - self%boxDim(1,1)
        Ly = self%boxDim(2,2) - self%boxDim(1,2)
        Lz = self%boxDim(2,3) - self%boxDim(1,3)

        write(self%fileUnit, *) Lx, 0.0, 0.0
        write(self%fileUnit, *) 0.0, Ly, 0.0
        write(self%fileUnit, *) 0.0, 0.0, Lz

      class default
        stop "Traj_POSCAR does not know how to handle this box type"

    end select

    write(self%fileUnit, "(9999(A,1x))") (AtomData(j)%Symb, j=1,nAtomTypes)
    write(self%fileUnit, *) box%nAtoms
    write(self%fileUnit, "(A)") "Cartesian"

    do iAtom = 1, box%nMaxAtoms
      molType = box%MolType(iAtom)
      atomType = box%AtomType(iAtom)
      atmSymb = AtomData(atomType)%Symb

!      write(*, *) iAtom, atomType, (box%atoms(jDim, iAtom), jDim=1,nDim)
      if(box%NMol(molType) >= box%MolSubIndx(iAtom) ) then
        write(self%fileUnit, *) (box%atoms(jDim, iAtom)-cm(jDim), jDim=1,nDim), atmSymb
      endif
    enddo

    flush(self%fileUnit)

  end subroutine
!====================================================================
  subroutine TrajPOSCAR_Prologue(self) 
    use BoxData, only: BoxArray
    implicit none
    class(TrajPOSCAR), intent(inout) :: self
    integer(kind=8) :: iCycle = 0

    select type(box => BoxArray(self%boxnum)%box)
      class default
        self%xLen = 2
        if(.not. allocated(self%boxDim)) then
          allocate( self%boxDim(1:2, 1:box%nDimension) )
        endif

    end select

    call self % WriteFrame(iCycle)
  end subroutine
!====================================================================
  subroutine TrajPOSCAR_Epilogue(self) 
    implicit none
    class(TrajPOSCAR), intent(inout) :: self
    integer(kind=8) :: iCycle = -1

    call self % WriteFrame(iCycle)
  end subroutine
!====================================================================
end module
!====================================================================
