!================================================================================
module Input_Forcefield
  use VarPrecision
  use Input_Format
  use ForcefieldData
  contains
!================================================================================
  subroutine Script_Forcefield(line, lineStat)
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat


      character(len=30) :: dummy, command 
      character(len=30) :: fileName      
      logical :: logicValue
      integer :: j
      integer :: intValue
      integer :: FFNum
      real(dp) :: realValue


      lineStat  = 0
      read(line, *) dummy, FFNum, command
      call LowerCaseLine(command)
      select case(trim(adjustl(command)))
!        case("fieldtype")
!          call Script_FieldType(line, lineStat) 

        case("readfile")
          read(line,*) dummy, FFNum, command, fileName
          call EnergyCalculator(FFNum) % Method % ReadParFile(fileName)

        case default
          lineStat = -1
      end select


  end subroutine

!================================================================================
  subroutine Script_FieldType(line, FFNum, lineStat)
    use ForcefieldData, only: nForceFields
    use FF_Pair_LJ_Cut, only: Pair_LJ_Cut
    use FF_Pair_LJ_Cut_NoNei, only: Pair_LJ_Cut_NoNei
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
      call LowerCaseLine(FF_Type)
      select case(trim(adjustl(FF_Type)))
        case("lj_cut")
          allocate(Pair_LJ_Cut::EnergyCalculator(FFNum) % Method)
          write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as LJ_Cut style"
        case("lj_cut_nonei")
          allocate(Pair_LJ_Cut_NoNei::EnergyCalculator(FFNum) % Method)
          write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as LJ_Cut (No Neighbor List) style"
        case default
          write(*,*) "Here"
          lineStat = -1

      end select

  end subroutine
!================================================================================
end module
!================================================================================
