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
    nType = -1
    molIndx = 0
    do iType = 1, nMolTypes
      nSum = nSum + box%NMol(iType)

!      write(*,*)
      if(globIndx <= nSum) then
        nType = iType
        exit
      endif
    enddo

    if(nType < 1) then
      stop "Catestrophic Error! Global Index passed to FindMolecule is invaild"
    endif

    molIndx = globIndx
    do iType = 1, nType-1
      molIndx = molIndx + box%NMolMax(iType) - box%NMol(iType)
    enddo

    
  end subroutine
!==================================================
  function FindFirstEmptyMol(box, moltype) result(atmIndx)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(SimpleBox), intent(in) :: box
    integer, intent(in) :: molType
    integer  :: atmIndx
    integer :: iType, molIndx


    molIndx = 0
    do iType = 1, molType - 1
      molIndx = molIndx + box%NMolMax(iType) 
    enddo

    molIndx = molIndx + box % NMol(molType) + 1

    atmIndx = molIndx
  end function
!==================================================
end module
!==================================================

