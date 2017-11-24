!================================================================================
module Input_Moves
  use VarPrecision
  use Input_Format
  use MCMoveData, only: Moves, MoveProb
  use MCMove_AtomTranslation, only: AtomTranslate
  contains
!================================================================================
  subroutine Script_MCMoves(line, moveNum, lineStat)
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: moveNum
    integer, intent(out) :: lineStat


      character(len=30) :: dummy, command 
      character(len=30) :: fileName      
      logical :: logicValue
      integer :: j
      integer :: FFNum
      real(dp) :: realValue


      lineStat  = 0
      read(line, *) command, realValue
      call LowerCaseLine(command)
      select case(trim(adjustl(command)))
        case("atomtranslation")
          allocate(AtomTranslate::Moves(moveNum)%move)
          MoveProb(moveNum) = realValue
        case default
          lineStat = -1
      end select


  end subroutine

!================================================================================
end module
!================================================================================
