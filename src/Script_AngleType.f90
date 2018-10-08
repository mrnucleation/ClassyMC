!================================================================================
module Input_AngleType
use Input_Format, only: LowerCaseLine, GetXCommand
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_AngleType(line, AngleNum, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Common_MolInfo, only: AngleData
    use IntraAngle_Ridgid, only: RidgidAngle 
!    use IntraAngle_Harmonic, only: HarmonicAngle
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: AngleNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, Angle_Type
!    character(len=30) :: fileName    
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(line, command, 1, lineStat)
    read(command, *) Angle_Type

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(Angle_Type)))
      case("ridgid")
        allocate(RidgidAngle :: AngleData(AngleNum)%angleFF)

!      case("harmonic")
!        allocate(HarmonicAngle :: AngleData(AngleNum)%angleFF)
!        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut style"
      case default
        write(*,*) Angle_Type
        lineStat = -1
        return
    end select
    call AngleData(AngleNum)%angleFF%ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
