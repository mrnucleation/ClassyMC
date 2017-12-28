
!================================================================================
module Input_NeighType
use Input_Format, only: GetXCommand
use Template_NeighList, only: NeighListDef
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_NeighType(line, lineStat)
    use BoxData, only: BoxArray
    use RSqListDef, only: RSqList
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: command
    integer :: i, intVal, boxNum

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) boxNum

    if( allocated(BoxArray(boxNum)%box%NeighList) ) then
      deallocate( BoxArray(boxNum)%box%NeighList )
    endif

    call GetXCommand(line, command, 3, lineStat)
    lineStat  = 0
    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(command)))
      case("rsqlist")
        call GetXCommand(line, command, 4, lineStat)
        read(command,*) intVal
        allocate( RSqList::BoxArray(boxNum)%box%NeighList(1:intVal) )

      case default
        lineStat = -1
        return
    end select


  end subroutine
!================================================================================
end module
!================================================================================
