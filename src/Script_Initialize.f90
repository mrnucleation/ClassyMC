!=================================================
module Input_Initialize
use VarPrecision
contains
!=================================================
  subroutine Script_Initialize
    use ParallelVar, only: p_size, myid, nout
    use BoxData, only: BoxArray
    use TrajData, only: TrajArray
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: initSeed, sgrnd
    implicit none

    if(initSeed < 0) then
      call system_clock(initSeed)
      initSeed = mod(initSeed, 10000)
    endif
    initSeed = p_size*initSeed + myid
    write(nout,*) "Random Generator Seed:", initSeed
    call sgrnd(initSeed)



  end subroutine
!=================================================
end module
!=================================================
