!====================================================================
module Traj_AENet
  use VarPrecision
  use TrajectoryTemplate, only: trajectory

  type, public, extends(trajectory) :: AENetDump
    logical :: recenter = .false.
!    integer :: fileUnit = -1
!    integer :: boxNum = -1
!    integer :: outFreq = 5000
!    character(len=50) :: fileName
    integer :: xLen
    real(dp), allocatable :: boxDim(:,:)
    character(len=50) :: boxStr
    contains
       procedure, pass :: WriteFrame => AENetDump_WriteFrame
!       procedure, pass :: SetUnit
!       procedure, pass :: SetBox
!       procedure, pass :: SetFileName
!       procedure, pass :: SetFreq
!       procedure, pass :: OpenFile
!       procedure, pass :: CloseFile
       procedure, pass :: Prologue => AENetDump_Prologue
       procedure, pass :: Epilogue => AENetDump_Epilogue

  end type
!====================================================================
  contains
!====================================================================
  subroutine AENetDump_WriteFrame(self, iCycle) 
    use BoxData, only: BoxArray
    use Common_MolInfo, only: AtomData
    implicit none
    class(AENetDump), intent(inout) :: self
    integer(kind=8), intent(in) :: iCycle
    integer :: iAtom, jDim, boxNum, i
    integer :: nDim, atomType, molType


    boxNum = self%boxNum
    nDim = BoxArray(boxNum)%box%nDimension

    write(self%fileUnit, *) 
    write(self%fileUnit, *) 
    write(self%fileUnit, *) "Atoms" 

    write(self%fileUnit, "(A)") "ITEM: NUMBER OF ATOMS "
    write(self%fileUnit, *) BoxArray(boxNum)%box%nAtoms

!        write(self%fileUnit, "(A)") self%boxstr
    write(self%fileUnit, "(A)") "ITEM: BOX BOUNDS pp pp pp"
    call BoxArray(boxNum)%box%GetDimensions(self%boxdim)
    do jDim = 1, nDim
      write(self%fileUnit, *) (self%boxdim(i, jDim), i=1,self%xLen)
    enddo

    write(self%fileUnit, "(A)") "ITEM: ATOMS id type x y z  "

    do iAtom = 1, BoxArray(boxNum)%box%nMaxAtoms
      molType = BoxArray(boxNum)%box%MolType(iAtom)
      atomType = BoxArray(boxNum)%box%AtomType(iAtom)

!      write(*, *) iAtom, atomType, (BoxArray(boxNum)%box%atoms(jDim, iAtom), jDim=1,nDim)
      if(BoxArray(boxNum)%box%NMol(molType) >= BoxArray(boxNum)%box%MolSubIndx(iAtom) ) then
        write(self%fileUnit, *) iAtom, atomType, (BoxArray(boxNum)%box%atoms(jDim, iAtom), jDim=1,nDim)
      endif
    enddo

    flush(self%fileUnit)

  end subroutine
!====================================================================
  subroutine AENetDump_Prologue(self) 
    use BoxData, only: BoxArray
    implicit none
    class(AENetDump), intent(inout) :: self

    select type(box => BoxArray(self%boxnum)%box)
      class default
        self%xLen = 2
!        self%boxstr = "ITEM: BOX BOUNDS xx yy zz"
        self%boxstr = "ITEM: BOX BOUNDS pp pp pp"
        allocate( self%boxDim(1:2, 1:box%nDimension) )

    end select

    call self % WriteFrame
  end subroutine
!====================================================================
  subroutine AENetDump_Epilogue(self) 
    implicit none
    class(AENetDump), intent(inout) :: self

    call self % WriteFrame
  end subroutine
!====================================================================
end module
!====================================================================
