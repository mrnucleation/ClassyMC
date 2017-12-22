!==================================================
module Box_Utility
use VarPrecision
use Template_SimBox, only: SimBox

!==================================================
contains
!==================================================
  subroutine FindAtom(box, globindx, atmIndx)
    implicit none
    class(SimBox), intent(in) :: box
    integer, intent(in) :: globIndx
    integer, intent(out) :: atmIndx
    integer :: iType, nSum, nType

    atmIndx = globindx
    do iType = 1, nMolTypes
      nSum = nSum + box%NMol(iType)*MolData(iType)%nAtoms
      if(globIndx < nSum) then
        nType = iType
        exit
      endif
    enddo

    atmIndx = globIndx
    do iType = 1, nType-1
      atmIndx = atmIndx + (box%NMolMax(iType) - box%NMol(iType))*MolData(iType)%nAtoms
    enddo

    
  end subroutine
!==================================================
  subroutine FindMolecule(box, globindx, molIndx)
    implicit none
    class(SimBox), intent(in) :: box
    integer, intent(in) :: globIndx
    integer, intent(out) :: molIndx
    integer :: iType, nSum, nType

    atmIndx = globindx
    do iType = 1, nMolTypes
      nSum = nSum + box%NMol(iType)
      if(globIndx < nSum) then
        nType = iType
        exit
      endif
    enddo

    atmIndx = globIndx
    do iType = 1, nType-1
      atmIndx = atmIndx + box%NMolMax(iType) - box%NMol(iType)
    enddo

    
  end subroutine
!==================================================
end module
!==================================================

