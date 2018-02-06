!====================================================================
module Traj_XYZ
  use VarPrecision
  use TrajectoryTemplate, only: trajectory

  type, public, extends(trajectory) :: trajXYZ
    logical :: padding = .false.
    logical :: recenter = .false.
!    integer :: fileUnit = -1
!    integer :: boxNum = -1
!    integer :: outFreq = 5000
!    character(len=50) :: fileName
    contains
       procedure, pass :: WriteFrame => TrajXYZ_WriteFrame
!       procedure, pass :: SetUnit
!       procedure, pass :: SetBox
!       procedure, pass :: SetFileName
!       procedure, pass :: SetFreq
!       procedure, pass :: OpenFile
!       procedure, pass :: CloseFile
       procedure, pass :: Prologue => TrajXYZ_Prologue
       procedure, pass :: Epilogue => TrajXYZ_Epilogue

  end type
!====================================================================
  contains
!====================================================================
  subroutine TrajXYZ_WriteFrame(self) 
    use BoxData, only: BoxArray
    use Common_MolInfo, only: AtomData
    implicit none
    class(trajXYZ), intent(in) :: self
    integer :: iAtom, jDim, boxNum
    integer :: nDim, atomType, molType
    real(dp) :: xcm, ycm, zcm


    boxNum = self%boxNum
    nDim = BoxArray(boxNum)%box%nDimension
    if(self%padding) then
      write(self%fileUnit, *) BoxArray(boxNum)%box%nMaxAtoms
    else
      write(self%fileUnit, *) BoxArray(boxNum)%box%nAtoms
    endif
    write(self%fileUnit, *) 

    do iAtom = 1, BoxArray(boxNum)%box%nMaxAtoms
      molType = BoxArray(boxNum)%box%MolType(iAtom)
      atomType = BoxArray(boxNum)%box%AtomType(iAtom)
      if(BoxArray(boxNum)%box%NMol(molType) < BoxArray(boxNum)%box%MolSubIndx(iAtom) ) then
        if(.not. self%padding) then
          cycle
        else
          write(self%fileUnit, *) AtomData(atomType)%symb, (1E15_dp, jDim=1,nDim)
        endif
      else
        write(self%fileUnit, *) AtomData(atomType)%symb, (BoxArray(boxNum)%box%atoms(jDim, iAtom), jDim=1,nDim)
      endif
      atomType = BoxArray(boxNum)%box%AtomType(iAtom)
    enddo


  end subroutine
!====================================================================
  subroutine TrajXYZ_Prologue(self) 
    implicit none
    class(trajXYZ), intent(inout) :: self

    call self % WriteFrame
  end subroutine
!====================================================================
  subroutine TrajXYZ_Epilogue(self) 
    implicit none
    class(trajXYZ), intent(inout) :: self

    call self % WriteFrame
  end subroutine
!====================================================================
end module
!====================================================================
!Notes: If the padding flag is true, the function will write dummy coordinates as a placeholder
!       for atoms which are not currently in the system, but have the potential of being added
!       by a swap move.  This is done because some visualization software prefer to have a fixed
!       number of atoms for each simulation frame. 
!====================================================================
