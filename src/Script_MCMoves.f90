!================================================================================
module Input_Moves
  use VarPrecision
  use Input_Format
  use MCMoveData, only: Moves, MoveProb

  use Move_AtomExchange, only: MC_AtomExchange
  use MCMove_AVBMC, only: AVBMC
  use MCMove_AnisoVol, only: AnisoVol
  use MCMove_AtomTranslation, only: AtomTranslate
  use MCMove_Basic_Swap, only: Basic_Swap
  use MCMove_Delete, only: MoveDelete
  use MCMove_MolTranslation, only: MolTranslate
  use MCMove_ParticleExchange, only: ParticleExchange
  use MCMove_PlaneRotation, only: PlaneRotate
  use MCMove_PlaneTranslation, only: PlaneTranslate
  use MCMove_PlaneAtomTranslation, only: PlaneAtomTranslate
  use MCMove_IsoVol, only: IsoVol
  use MCMove_UB_Swap, only: UB_Swap
  use MCMove_VolExchange, only: VolExchange
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
      case("avbmc")
        allocate(AVBMC::Moves(moveNum)%move)

      case("anisovol")
        allocate(AnisoVol::Moves(moveNum)%move)

      case("atomtranslation")
        allocate(AtomTranslate::Moves(moveNum)%move)

      case("atomexchange")
        allocate(MC_AtomExchange::Moves(moveNum)%move)

      case("basicswap")
        allocate(Basic_Swap::Moves(moveNum)%move)

      case("debugdelete")
        allocate(MoveDelete::Moves(moveNum)%move)

      case("moltranslation")
        allocate(MolTranslate::Moves(moveNum)%move)

      case("thermolambda")
        allocate(ThermoLambda::Moves(moveNum)%move)

      case("isovol")
        allocate(IsoVol::Moves(moveNum)%move)

      case("particleexchange")
        allocate(ParticleExchange::Moves(moveNum)%move)

      case("planerotate")
        allocate(PlaneRotate::Moves(moveNum)%move)

      case("planeatomtranslation")
        allocate(PlaneAtomTranslate::Moves(moveNum)%move)

      case("planetranslation")
        allocate(PlaneTranslate::Moves(moveNum)%move)

      case("volexchange")
        allocate(VolExchange::Moves(moveNum)%move)

      case("ubswap")
        allocate(UB_Swap::Moves(moveNum)%move)

      case default
        lineStat = -1
    end select

  end subroutine
!================================================================================
end module
!================================================================================
