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

    use Analysis_AngleDistribution, only: AngleDistribution
    use Analysis_BondDistribution, only: BondDistribution
    use Analysis_TorsionDistribution, only: TorsionDistribution
    use Anaylsis_TotalSize, only: TotalSize
    use Anaylsis_BlockAverage, only: BlockAverage
    use Anaylsis_ClusterSize, only: ClusterSize
    use Anaylsis_DensityOfStates, only: DensityOfStates
    use Anaylsis_DistPair, only: DistPair
    use Analysis_RDF, only: rdf
    use Anaylsis_MolFractionHist, only: MolFractionHist
    use Anaylsis_ThermoAverage, only: ThermoAverage
    use Anaylsis_ThermoIntegration, only: ThermoIntegration
#ifdef EMBPYTHON
    use Anaylsis_PythonFunc, only: PythonFunc
#endif

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

      case("angledistribution")
        allocate(AngleDistribution::AnalysisArray(AnaNum) % func)

      case("bonddistribution")
        allocate(BondDistribution::AnalysisArray(AnaNum) % func)

      case("blockaverage")
        allocate(blockAverage::AnalysisArray(AnaNum) % func)

      case("clustersize")
        allocate(ClusterSize::AnalysisArray(AnaNum) % func)

      case("densityofstates")
        allocate(DensityOfStates::AnalysisArray(AnaNum) % func)

      case("distpair")
        allocate(DistPair::AnalysisArray(AnaNum) % func)

      case("molfractionhist")
        allocate(MolFractionHist::AnalysisArray(AnaNum) % func)

      case("thermoaverage")
        allocate(thermoAverage::AnalysisArray(AnaNum) % func)

      case("thermointegration")
        allocate(ThermoIntegration::AnalysisArray(AnaNum) % func)

      case("torsiondistribution")
        allocate(TorsionDistribution::AnalysisArray(AnaNum) % func)

      case("totalsize")
        allocate(TotalSize::AnalysisArray(AnaNum) % func)

      case("rdf")
        allocate(rdf::AnalysisArray(AnaNum) % func)

#ifdef EMBPYTHON
      case("python")
        allocate(PythonFunc::AnalysisArray(AnaNum) % func)
#endif

     case default
        lineStat = -1
        write(0,*) "ERROR! Unknown Analysis Type Given in Input File!!!"
        error stop

    end select

    if(lineStat == -1) then
      return
    endif

    AnalysisArray(AnaNum) % func % analyID = AnaNum
    call AnalysisArray(AnaNum) % func % ProcessIO(line)

  end subroutine
!================================================================================
end module
!================================================================================
