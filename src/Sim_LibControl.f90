!===========================================================================
! This file contains the code for the C interface to ClassyMC.
!===========================================================================
#define __StdErr__ 0
!===========================================================================
module Library
  use ParallelVar, only: myid, ierror, nout
  use VarPrecision

!===========================================================================
contains
  !===========================================================================
  subroutine Library_ReadScript(strlen, cfilename) bind(C,name='Classy_ReadScript')
    use C_F_Routines, only: C_F_String
    use ScriptInput, only: Script_ReadParameters
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_null_char, c_loc
    implicit none
    integer(c_int), intent(in), value :: strlen
    type(c_ptr), value, intent(in) :: cfilename
    character(len=strlen):: infile
    integer :: i, slen

    infile = c_f_string(cfilename)
    call Script_ReadParameters(infile)
 
  end subroutine
  !===========================================================================
  subroutine Library_FullSimulation() bind(C,name='Classy_FullSim')
    use SimMonteCarlo, only: RunMonteCarlo
    implicit none

    call RunMonteCarlo
 
  end subroutine
!===========================================================================
  subroutine Library_RunMove() bind(C,name='Classy_RunMove')
    use BoxData, only: BoxArray
    use CommonSampling, only: Sampling
    use MCMoveData, only: Moves, MoveProb
    use MoveClassDef, only: MCMove
    use SimMonteCarlo, only: Analyze, Update
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use RandomGen, only: sgrnd, grnd, ListRNG
    implicit none
    logical :: accept
    integer(kind=8) :: iCycle, iMove
    integer :: moveNum, boxNum, nBoxes
    real(dp), allocatable :: boxProb(:)


    nBoxes = size(BoxArray)

    moveNum = ListRNG(MoveProb) !Randomly select a move to perform
    select type( curMove => Moves(moveNum) % Move )
      class is (MCMultiBoxMove) ! Mutli Box Move
        call curMove % MultiBox (accept)
      class is (MCMove) ! Single Box Move
        if(nBoxes > 1) then
          call curMove % GetBoxProb(boxProb)
          boxNum = ListRNG(boxProb)
        else
          boxNum = 1
        endif
        call curMove % FullMove(BoxArray(boxNum)%box, accept)
    end select
    call Sampling%UpdateStatistics(accept)
    if(accept) then
      call Update(accept)
    endif

    call Analyze(iCycle, iMove, accept, .true.) !Per Move Analysis


  end subroutine
!===========================================================================
!  subroutine Library_ForceMove(movenum) bind(C, name='Classy_ForceMove')
!    implicit none
!  end subroutine
!===========================================================================
  subroutine Library_ScreenOut() bind(C,name='Classy_ScreenOut')
    use SimMonteCarlo, only: ScreenOut
    implicit none

    integer(kind=8) :: iMove, iCycle

    iCycle = 0
    iMove = 0
    call ScreenOut(iMove, iCycle)

 
  end subroutine
!===========================================================================
end module
!===========================================================================