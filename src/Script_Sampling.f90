!================================================================================
module Input_Sampling
  use VarPrecision
  use Input_Format
  use CommonSampling, only: sampling
!================================================================================
  contains
!================================================================================
  subroutine Script_SamplingType(iLine, lineStore, lineStat)
    use MetropolisRule, only: Metropolis
    use MinMetroRule, only: MinMetro
    use UmbrellaRule, only: Umbrella
    use UmbrellaWHAMRule, only: UmbrellaWHAM
    implicit none

    character(len=maxLineLen), allocatable :: lineStore(:)
    integer, intent(in) :: iLine
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command 
    character(len=30) :: fileName      
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(lineStore(iLine), command, 2, lineStat)
    if(lineStat < 0) then
      return
    endif
    call LowerCaseLine(command)
    select case(trim(adjustl(command)))
      case("metropolis")
        allocate(metropolis::sampling)

      case("min")
        allocate(MinMetro::sampling)

      case("umbrella")
        allocate(Umbrella::sampling)

      case("umbrellawham")
        allocate(UmbrellaWHAM::sampling)

      case default
        lineStat = -1
    end select

  end subroutine
!================================================================================
end module
!================================================================================
