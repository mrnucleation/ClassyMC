!================================================================================
module Input_Sampling
  use VarPrecision
  use Input_Format
  use CommonSampling, only: sampling
!================================================================================
  contains
!================================================================================
  subroutine Script_SamplingType(line, lineStat)
    use MetropolisRule, only: metropolis
    use MinMetroRule, only: MinMetro
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command 
    character(len=30) :: fileName      
    real(dp) :: realValue

    lineStat  = 0
    call GetXCommand(line, command, 2, lineStat)
    if(lineStat < 0) then
      return
    endif
    call LowerCaseLine(command)
    select case(trim(adjustl(command)))
      case("metropolis")
        allocate(metropolis::sampling)
      case("min")
        allocate(MinMetro::sampling)

      case default
        lineStat = -1
    end select

  end subroutine
!================================================================================
end module
!================================================================================
