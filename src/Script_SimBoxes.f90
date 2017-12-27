!================================================================================
module Input_SimBoxes
  use VarPrecision
  use Input_Format
  use BoxData, only: BoxArray
  use SimpleSimBox, only: Simplebox
  use CubicBoxDef, only: CubeBox
  contains
!================================================================================
  subroutine Script_BoxType(line, boxNum, lineStat)
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: boxNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command 
    character(len=30) :: fileName      
    logical :: logicValue
    integer :: j
    integer :: intValue
    integer :: FFNum
    real(dp) :: realValue

    lineStat  = 0
    read(line, *) command
    select case(trim(adjustl(command)))
      case("nobox")
        allocate(SimpleBox::BoxArray(boxNum)%box)

      case("cube")
        allocate(CubeBox::BoxArray(boxNum)%box) 

      case default
        lineStat = -1
    end select

  end subroutine
!================================================================================
end module
!================================================================================
