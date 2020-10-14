!===========================================================================
#define __StdErr__ 0
!===========================================================================
module SimMinimize
  use ParallelVar, only: myid, ierror, nout
use VarPrecision
!===========================================================================
contains
!===========================================================================
  subroutine RunMinimize(boxnum)
#ifdef PARALLEL
    use MPI
#endif

    use BoxData, only: BoxArray
    use SimpleSimBox, only: SimpleBox
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use Output_DumpCoords, only: Output_DumpData
    use SimControl, only: nMoves, nCycles, screenfreq, configfreq, energyCheck, lRate, ETol, ForceTol
    use SimMonteCarlo, only: ScreenOut, Prologue, Epilogue, SafetyCheck, Trajectory
    use Units, only: outEngUnit

    implicit none
 
    integer, intent(in) :: boxNum

    logical :: accept
    integer(kind=8) :: iCycle, iMove
    class(SimpleBox), pointer :: trialBox => null()

    character(len=30), parameter :: estring = "energy"
    integer :: nBoxes, maxatoms, thermonum, iAtom
    real(dp) :: E_Curr, E_New
    real(dp) :: lrate_eff
    real(dp), allocatable :: temppos(:,:) 
    real(dp), pointer :: atoms(:,:)
    real(dp), pointer :: forces(:,:)

    trialBox => BoxArray(boxnum)%box
    
    iCycle = 0
    iMove = 0

    thermonum = trialBox%ThermoLookUp(estring)
    call Prologue
    call SafetyCheck
    nBoxes = size(BoxArray)

    call trialBox%ComputeForces

    maxatoms = trialBox%GetMaxAtoms()
    call trialBox%GetCoordinates(atoms)
    call trialBox%GetForceArray(forces)

    allocate( temppos(1:3,1:maxatoms) )

    !-------Main Monte Carlo Simulation Loop-------
    write(nout, *) "============================================"
    write(nout, *) "      Starting Miniminzation...  "
    write(nout, *) "============================================"

    flush(nout)

    call ScreenOut(iCycle, iMove)
    iCycle = 0
    lrate_eff = lrate
    do
      iCycle = iCycle + 1
      if(mod(iCycle, int(screenfreq,8)) == 0) then
        call ScreenOut(iCycle, iMove)
        flush(nout)
      endif

      temppos(1:3, 1:maxatoms) = atoms(1:3, 1:maxatoms)
      call trialBox%ComputeEnergy
      E_Curr = trialBox%GetThermo(thermonum)
      call trialBox%ComputeForces
      atoms(1:3, 1:maxatoms) = temppos(1:3,1:maxatoms) - lrate*forces(1:3,1:maxatoms)
      do iAtom = 1, maxatoms
        if(.not. trialbox%IsActive(iAtom)) cycle
        call trialbox%Boundary(atoms(1,iAtom), atoms(2,iAtom), atoms(3,iAtom))

      enddo
      call trialBox%ComputeEnergy

      E_New = trialBox%GetThermo(thermonum)
      if(E_New < E_Curr) then
        trialBox%forceoutofdate = .true.
      else
        write(nout,*) "New position worse than previous one, try lowering learning rate"
        atoms(1:3,1:maxatoms) = temppos(1:3,1:maxatoms) 
        exit
      endif

      call Trajectory(iCycle, iMove)
      call Output_DumpData

      if(abs(E_New-E_Curr) < ETol) exit

    enddo
    !-------End of Main Minimize Simulation Loop-------
 
    call ScreenOut(iCycle, iMove)
    write(nout,*) "======================================="
    write(nout,*) "     Simulation End"
    write(nout,*) "======================================="

    call Epilogue

#ifdef PARALLEL
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
#endif

    call Output_DumpData
      
  end subroutine
!===========================================================================
end module
!===========================================================================
