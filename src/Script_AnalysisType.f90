
!================================================================================
module Input_AnalysisType
use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_AnalysisType(line, AnaNum, lineStat)
    use Input_Format, only: GetXCommand
    use AnalysisData, only: AnalysisArray
    use ParallelVar, only: nout
    use Analysis_RDF, only: rdf
    use Anaylsis_ThermAverage, only: ThermAverage
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: AnaNum
    integer, intent(out) :: lineStat

    character(len=30) :: command
    integer :: i, intVal

    lineStat  = 0
    call GetXCommand(line, command, 1, lineStat)

    !Safety check to ensure that the index number is within proper bounds
    select case(trim(adjustl(command)))
      case("rdf")
        allocate(rdf::AnalysisArray(AnaNum) % func)

      case("thermoaverage")
        allocate(thermAverage::AnalysisArray(AnaNum) % func)

      case default
        lineStat = -1
        return
    end select

    if(lineStat == -1) then
      return
    endif

    call AnalysisArray(AnaNum) % func % ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
