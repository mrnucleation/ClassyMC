
!================================================================================
module Input_TrajType
use Input_Format, only: GetXCommand
use TrajData, only: TrajArray
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_TrajType(line, TrajNum, lineStat)
    use ParallelVar, only: nout
    use Traj_Lammps, only: LAMMPSDump
    use Traj_POSCAR, only: TrajPOSCAR
    use Traj_XYZ, only: trajXYZ
    use Traj_XSF, only: trajXSF
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: TrajNum
    integer, intent(out) :: lineStat

    character(len=30) :: command
    integer :: i, intVal

    lineStat  = 0
    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("dump")
        allocate(LAMMPSDump::TrajArray(TrajNum) % traj)

      case("poscar")
        allocate(TrajPOSCAR::TrajArray(TrajNum) % traj)

      case("xyz")
        allocate(trajXYZ::TrajArray(TrajNum) % traj)

      case("xsf")
        allocate(trajXSF::TrajArray(TrajNum) % traj)

      case default
        lineStat = -1
        return
    end select

!    write(*,*) line
!    call TrajArray(TrajNum) % traj % SetUnit(1000+TrajNum)

    call GetXCommand(line, command, 2, lineStat)
    read(command, *) intVal
    call TrajArray(TrajNum) % traj % SetBox(intVal)

    call GetXCommand(line, command, 3, lineStat)
    read(command, *) intVal
    call TrajArray(TrajNum) % traj % SetFreq(intVal)

    call GetXCommand(line, command, 4, lineStat)
    do i = 1, len(command)
      if(command(i:i) == '"') then
        command(i:i) = " "
      endif
    enddo
!    write(*,*) command
    call TrajArray(TrajNum) % traj % SetFileName(command)
    call TrajArray(TrajNum) % traj % OpenFile



  end subroutine
!================================================================================
end module
!================================================================================
