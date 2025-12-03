
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
    use FF_Ext_LAMMPS, only: Pair_LAMMPS
    use FF_EasyEP_LJ_Cut, only: EP_LJ_Cut
    use FF_EasyEP_LJ_CutShift, only: EP_LJ_CutShift
    use FF_EasyEP_LJ_Ele_Cut, only: EP_LJ_Ele_Cut
    use FF_EasyEP_Pedone_Cut, only: EP_Pedone_Cut
!    use FF_Pair_LJ_Wall, only: Pair_LJ_Wall
!    use FF_Pair_LJ_Shift, only: Pair_LJ_Shift
!    use FF_Pair_LJ_Cut_NoNei, only: Pair_LJ_Cut_NoNei
!    use FF_Pair_LJ_Q_Cut, only: Pair_LJ_Q_Cut
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

        !Phasing out the lj_cut in favor of the EP style
!      case("lj_cut")
!        allocate(Pair_LJ_Cut::EnergyCalculator(FFNum) % Method)
!        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut style"

      case("lj_cut")
        allocate(EP_LJ_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut style (EasyPair)"

      case("lj_ele_cut")
        allocate(EP_LJ_Ele_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6-1 LJ/Q Cut style (EasyPair)"


      case("lj_shift")
!        allocate(Pair_LJ_Shift::EnergyCalculator(FFNum) % Method)
        allocate(EP_LJ_CutShift::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 12-6 LJ Cut/Shift style (EasyPair)"

!      case("lj_q_cut")
!        allocate(Pair_LJ_Q_Cut::EnergyCalculator(FFNum) % Method)
!        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as a 12-6 LJ w/ Eletrostatic Cut style"

!      case("lj_wall")
!        allocate(Pair_LJ_Cut::EnergyCalculator(FFNum) % Method)
!        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as 9-3 LJ Wall style"

      case("pedone")
        allocate(EP_Pedone_Cut::EnergyCalculator(FFNum) % Method)
!        allocate(Pair_Pedone_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Pedone Cut style"

        !Depreciating the old Pedone code in favor of the EP Style
!      case("pedone")
!        allocate(Pair_Pedone_Cut::EnergyCalculator(FFNum) % Method)
!        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as Pedone Cut style"


#ifdef LAMMPS
      case("lammps")
        allocate(Pair_LAMMPS::EnergyCalculator(FFNum) % Method)
        write(nout,"(1x,A,I2,A)") "Forcefield", FFNum, " allocated as LAMMPS Library Interface"
#endif

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
