
!================================================================================
module Input_FieldType
use Input_Format, only: LowerCaseLine
use VarPrecision
use ForcefieldData
!================================================================================
contains
!================================================================================
  subroutine Script_FieldType(line, FFNum, lineStat)
    use ForcefieldData, only: nForceFields
    use FF_Pair_LJ_Cut, only: Pair_LJ_Cut
    use FF_Pair_LJ_Cut_NoNei, only: Pair_LJ_Cut_NoNei
    use FF_Pair_Tersoff, only: Pair_Tersoff
    use ParallelVar, only: nout
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: FFNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, FF_Type
    character(len=30) :: fileName      
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    read(line, *) FF_Type

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(FF_Type)))
      case("lj_cut")
        allocate(Pair_LJ_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as LJ_Cut style"

      case("lj_cut_nonei")
        allocate(Pair_LJ_Cut_NoNei::EnergyCalculator(FFNum) % Method)
        write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as LJ_Cut (No Neighbor List) style"

      case("tersoff")
        allocate(Pair_Tersoff::EnergyCalculator(FFNum) % Method)
        write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as Tersoff style"

      case default
        write(*,*) "Here"
        lineStat = -1

    end select

  end subroutine
!================================================================================
end module
!================================================================================