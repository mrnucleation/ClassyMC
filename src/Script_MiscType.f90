!================================================================================
module Input_MiscType
use Input_Format, only: LowerCaseLine, GetXCommand
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_MiscType(line, molType, MiscNum, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Common_MolInfo, only: MolData
    use Misc_IntraPair_1_5, only: Pair_1_5
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: MiscNum, molType
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
      case("1_5_pair")
        allocate(Pair_1_5 :: MolData(molType)% miscdata(MiscNum) % miscFF)

      case default
        write(*,*) Misc_Type
        lineStat = -1
        return
    end select
    call MolData(molType)%MiscData(MiscNum)%miscFF%SetMolType(molType)
    call MolData(molType)%MiscData(MiscNum)%miscFF%ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
