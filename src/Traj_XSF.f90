!====================================================================
module Traj_XSF
  use VarPrecision
  use TrajectoryTemplate, only: trajectory

  type, public, extends(trajectory) :: trajXSF
    logical :: padding = .false.
    logical :: recenter = .false.
!    integer :: fileUnit = -1
!    integer :: boxNum = -1
!    integer :: outFreq = 5000
!    character(len=50) :: fileName
    contains
       procedure, pass :: WriteFrame => TrajXSF_WriteFrame
!       procedure, pass :: SetUnit
!       procedure, pass :: SetBox
!       procedure, pass :: SetFileName
!       procedure, pass :: SetFreq
!       procedure, pass :: OpenFile
!       procedure, pass :: CloseFile
       procedure, pass :: Prologue => TrajXSF_Prologue
       procedure, pass :: Epilogue => TrajXSF_Epilogue

  end type
!====================================================================
  contains
!====================================================================
  subroutine TrajXSF_WriteFrame(self, iCycle) 
    use BoxData, only: BoxArray
    use SimpleSimBox, only: SimpleBox
    use CubicBoxDef, only: CubeBox 
    use OrthoBoxDef, only: OrthoBox 
    use Common_MolInfo, only: AtomData
    implicit none
    class(trajXSF), intent(inout) :: self
    integer(kind=8), intent(in) :: iCycle
    integer :: iAtom, jDim, boxNum
    integer :: nDim, atomType, molType
    real(dp) :: cm(1:3)
    real(dp), parameter :: eVConv = 8.6173303E-5_dp
    class(SimpleBox), pointer :: box
    real(dp) :: Lx, Ly, Lz
    real(dp) :: dimensions(1:2, 1:3)


    boxNum = self%boxNum
    box => BoxArray(boxNum)%box
    nDim = box%nDimension
    write(self%fileUnit, "(A, F25.10, A)") "# total energy =", box%ETotal * eVConv , " ev"
    write(self%fileUnit, *) 

    select type(box)
      type is(SimpleBox)
        write(self%fileUnit, *) "ATOMS"
        cm = 1E80_dp
        do iAtom = 1, box%nMaxAtoms
          molType = box%MolType(iAtom)
          atomType = box%AtomType(iAtom)
          if(box%NMol(molType) >= box%MolSubIndx(iAtom) ) then
            do jDim =1, 3
              if(box%atoms(jDim, iAtom) < cm(jDim)) then
                cm(jDim) = box%atoms(jDim, iAtom) - 1e-5
              endif
            enddo
          endif
        enddo

      class is(CubeBox)
        write(self%fileUnit, *) "CRYSTAL"
        write(self%fileUnit, *) "PRIMVEC"
        call box%GetDimensions(dimensions(1:2, 1:3))
        cm(1) = dimensions(1,1)
        cm(2) = dimensions(1,2)
        cm(3) = dimensions(1,3)
        Lx = dimensions(2,1) - dimensions(1,1)
        Ly = dimensions(2,2) - dimensions(1,2)
        Lz = dimensions(2,3) - dimensions(1,3)

        write(self%fileUnit, *) Lx, 0.0, 0.0
        write(self%fileUnit, *) 0.0, Ly, 0.0
        write(self%fileUnit, *) 0.0, 0.0, Lz
        

        write(self%fileUnit, *) "PRIMCOORD"
        write(self%fileUnit, *) box%nAtoms, 1

      class is(OrthoBox)
        write(self%fileUnit, *) "CRYSTAL"
        write(self%fileUnit, *) "PRIMVEC"
        call box%GetDimensions(dimensions(1:2, 1:3))
        cm(1) = dimensions(1,1)
        cm(2) = dimensions(1,2)
        cm(3) = dimensions(1,3)
        Lx = dimensions(2,1) - dimensions(1,1)
        Ly = dimensions(2,2) - dimensions(1,2)
        Lz = dimensions(2,3) - dimensions(1,3)

        write(self%fileUnit, *) Lx, 0.0, 0.0
        write(self%fileUnit, *) 0.0, Ly, 0.0
        write(self%fileUnit, *) 0.0, 0.0, Lz
        

        write(self%fileUnit, *) "PRIMCOORD"
        write(self%fileUnit, *) box%nAtoms, 1


      class default
        stop "Traj_XSF does not know how to handle this box type"

    end select



    do iAtom = 1, box%nMaxAtoms
      molType = box%MolType(iAtom)
      atomType = box%AtomType(iAtom)
      if(box%NMol(molType) < box%MolSubIndx(iAtom) ) then
        cycle
      else
        write(self%fileUnit, *) AtomData(atomType)%symb, (box%atoms(jDim, iAtom)-cm(jDim), jDim=1,nDim), (0.0, jDim=1,3)
      endif
    enddo


  end subroutine
!====================================================================
  subroutine TrajXSF_Prologue(self) 
    implicit none
    class(trajXSF), intent(inout) :: self
    integer(kind=8) :: iCycle = 0

!    call self % WriteFrame(iCycle)
  end subroutine
!====================================================================
  subroutine TrajXSF_Epilogue(self) 
    implicit none
    class(trajXSF), intent(inout) :: self
    integer(kind=8) :: iCycle = -1

    call self % WriteFrame(iCycle)
  end subroutine
!====================================================================
end module
!====================================================================
!Notes: If the padding flag is true, the function will write dummy coordinates as a placeholder
!       for atoms which are not currently in the system, but have the potential of being added
!       by a swap move.  This is done because some visualization software prefer to have a fixed
!       number of atoms for each simulation frame. 
!====================================================================
