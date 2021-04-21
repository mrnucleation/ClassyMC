!================================================================================
module Input_MiscType
use Input_Format, only: LowerCaseLine, GetXCommand
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_MiscType(line, MiscNum, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Common_MolInfo, only: MiscData
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: MiscNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, Misc_Type
!    character(len=30) :: fileName    
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(line, command, 1, lineStat)
    read(command, *) Misc_Type

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(Misc_Type)))
      case("15_pair")
!        allocate( :: MiscData(MiscNum)%miscFF)

      case default
        write(*,*) Misc_Type
        lineStat = -1
        return
    end select
    call MiscData(MiscNum)%miscFF%ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
