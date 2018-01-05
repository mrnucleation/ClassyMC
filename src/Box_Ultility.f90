!==================================================
module Box_Utility
use VarPrecision
use SimpleSimBox, only: SimpleBox

!==================================================
contains
!==================================================
  subroutine FindAtom(box, globindx, atmIndx)
    use Common_MolInfo, only: MolData, nMoltypes
    implicit none
    class(SimpleBox), intent(in) :: box
    integer, intent(in) :: globIndx
    integer, intent(out) :: atmIndx
    integer :: iType, nSum, nType

    nSum = 0
    do iType = 1, nMolTypes
      nSum = nSum + box%NMol(iType)*MolData(iType)%nAtoms
      if(globIndx <= nSum) then
        nType = iType
        exit
      endif
    enddo

    atmIndx = globIndx
    if( nType == 1 ) then
      return
    endif

  
    do iType = 1, nType-1
      atmIndx = atmIndx + (box%NMolMax(iType) - box%NMol(iType))*MolData(iType)%nAtoms
    enddo

    
  end subroutine
!==================================================
  subroutine FindMolecule(box, globindx, molIndx)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(SimpleBox), intent(in) :: box
    integer, intent(in) :: globIndx
    integer, intent(out) :: molIndx
    integer :: iType, nSum, nType

    nSum = 0
    do iType = 1, nMolTypes
      nSum = nSum + box%NMol(iType)
      if(globIndx <= nSum) then
        nType = iType
        exit
      endif
    enddo

    molIndx = globIndx
    do iType = 1, nType-1
      molIndx = molIndx + box%NMolMax(iType) - box%NMol(iType)
    enddo

    
  end subroutine
!==================================================
  subroutine FindFirstEmptyMol(box, moltype, atmIndx)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(SimpleBox), intent(in) :: box
    integer, intent(in) :: molType
    integer, intent(out) :: atmIndx
    integer :: iType, molIndx


    molIndx = 0
    do iType = 1, molType - 1
      molIndx = molIndx + box%NMolMax(iType) 
    enddo

    molIndx = molIndx + box % NMol(molType) + 1
    
  end subroutine
!==================================================
end module
!==================================================

