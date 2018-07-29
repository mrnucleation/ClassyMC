!================================================================================
module Input_Moves
  use VarPrecision
  use Input_Format
  use MCMoveData, only: Moves, MoveProb

  use Move_AtomExchange, only: AtomExchange
  use MCMove_AtomTranslation, only: AtomTranslate
  use MCMove_Delete, only: MoveDelete
  use MCMove_UB_Simple, only: UB_Simple
  use Move_ThermoLambda, only: ThermoLambda

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
    MoveProb(moveNum) = realValue

    select case(trim(adjustl(command)))
      case("atomtranslation")
        allocate(AtomTranslate::Moves(moveNum)%move)

      case("atomexchange")
        allocate(AtomExchange::Moves(moveNum)%move)

      case("debugdelete")
        allocate(MoveDelete::Moves(moveNum)%move)

      case("thermolambda")
        allocate(ThermoLambda::Moves(moveNum)%move)

      case("ubswap")
        allocate(UB_Simple::Moves(moveNum)%move)

      case default
        lineStat = -1
    end select

  end subroutine
!================================================================================
end module
!================================================================================
