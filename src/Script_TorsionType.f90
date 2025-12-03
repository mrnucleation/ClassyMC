!================================================================================
module Input_TorsionType
use Input_Format, only: LowerCaseLine, GetXCommand
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_TorsionType(line, TorsionNum, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Common_MolInfo, only: TorsionData
    use IntraTorsion_Ridgid, only: RidgidTorsion 
    use IntraTorsion_TRAPPE, only: TRAPPETorsion
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: TorsionNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, Torsion_Type
!    character(len=30) :: fileName    
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(line, command, 1, lineStat)
    read(command, *) Torsion_Type

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(Torsion_Type)))
      case("ridgid")
        allocate(RidgidTorsion :: TorsionData(TorsionNum)%torsionFF)

      case("trappe")
        allocate(TRAPPETorsion :: TorsionData(TorsionNum)%torsionFF)

      case default
        write(0,*) "Unknown Torsional Angle Type Given in the Forcefield File!"
        write(0,*) Torsion_Type
        lineStat = -1
        error stop

    end select
    call TorsionData(TorsionNum)%torsionFF%ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
