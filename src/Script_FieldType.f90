
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
    use FF_Einstein, only: Pair_Einstein
    use FF_HardSphere, only: Pair_HardSphere
    use FF_Hybrid, only: Pair_Hybrid
    use FF_Pair_LJ_Cut, only: Pair_LJ_Cut
    use FF_Pair_LJ_Shift, only: Pair_LJ_Shift
!    use FF_Pair_LJ_Cut_NoNei, only: Pair_LJ_Cut_NoNei
    use FF_Pair_LJ_Q_Cut, only: Pair_LJ_Q_Cut
    use FF_Pair_Pedone_Cut, only: Pair_Pedone_Cut
    use FF_Pair_Tersoff, only: Pair_Tersoff
    use FF_ThermoIntegration, only: Pair_ThermoIntegration

#ifdef AENET
    use FF_AENet, only: Pair_AENet
#endif

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

    !Forcefield objects go here!
    select case(trim(adjustl(FF_Type)))
#ifdef AENET
      case("aenet")
        allocate(Pair_AENet::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as AENet style"
#endif

      case("einstein")
        allocate(Pair_Einstein::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Einstein Crystal"

      case("hardsphere")
        allocate(Pair_HardSphere::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Hard Sphere"

      case("hybrid")
        allocate(Pair_Hybrid::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as a Hybrid Forcefield"

      case("lj_cut")
        allocate(Pair_LJ_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut style"

      case("lj_shift")
        allocate(Pair_LJ_Shift::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut/Shift style"

      case("lj_q_cut")
        allocate(Pair_LJ_Q_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as a 12-6 LJ w/ Eletrostatic Cut style"

      case("pedone")
        allocate(Pair_Pedone_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Pedone Cut style"

      case("tersoff")
        allocate(Pair_Tersoff::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Tersoff style"

      case("thermointegration")
        allocate(Pair_ThermoIntegration::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Thermo Integration Style"

      case default
!        write(*,*) "Here"
        lineStat = -1

    end select

  end subroutine
!================================================================================
end module
!================================================================================
