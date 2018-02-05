!===============================
module BoxPresets
use VarPrecision

!===============================
contains
!===============================
  subroutine Crystal(iLine, lineStore, boxNum, lineStat)
    use Input_Format, only: GetXCommand
    implicit none

    integer, intent(in) :: iLine, boxNum
    character(len=maxLineLen), intent(in) :: lineStore(:)   
    integer, intent(out) :: lineStat

    character(len=30) :: command
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(lineStore(iLine), command, 3, lineStat)

    select case(trim(adjustl(command)))
      case("fcc")
        call GetXCommand(lineStore(iLine), command, 4, lineStat)
        read(command, *) realVal
        call GetXCommand(lineStore(iLine), command, 5, lineStat)
        read(command, *) intVal
        call FCC(BoxNum, realVal, intVal)

      case default
        lineStat = -1

    end select


  end subroutine
!===============================
  subroutine FCC(BoxNum, latConst, replicate)
    use BoxData, only: BoxArray
    use OrthoBoxDef, only: OrthoBox

    implicit none
    integer, intent(in) :: BoxNum, replicate
    real(dp), intent(in) :: latConst

    integer :: nAtoms
    real(dp) :: boxL
    real(dp) :: baseCell(1:4, 1:3)

    allocate( OrthoBox::BoxArray(BoxNum)%box )

    baseCell(1, 1) = 0.0
    baseCell(1, 2) = 0.0
    baseCell(1, 3) = 0.0

    baseCell(2, 1) = 0.0
    baseCell(2, 2) = 0.5*latConst
    baseCell(2, 2) = 0.5*latConst

    baseCell(3, 1) = 0.5*latConst
    baseCell(3, 2) = 0.5*latConst
    baseCell(3, 3) = 0.0

    baseCell(4, 1) = 0.5*latConst
    baseCell(4, 2) = 0.0
    baseCell(4, 3) = 0.5*latConst





    boxL = replicate * latConst


  end subroutine
!===============================
end module
!===============================
