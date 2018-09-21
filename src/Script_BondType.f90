!================================================================================
module Input_BondType
use Input_Format, only: LowerCaseLine
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_BondType(line, BondNum, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Common_MolInfo, only: BondData

    use IntraBond_Harmonic, only: HarmonicBond
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: BondNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, Bond_Type
!    character(len=30) :: fileName    
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(line, command, 2, lineStat)
    read(command, *) Bond_FF

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(Bond_FF)))
      case("ridgid")

      case("harmonic")
        allocate(HarmonicBond :: BondData(BondNum)%bondFF)
!        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut style"
      case default
!        write(*,*) "Here"
        lineStat = -1
        return
    end select
    call BondData(BondNum)%bondFF%ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
