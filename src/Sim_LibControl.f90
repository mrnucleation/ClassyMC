!===========================================================================
! This file contains the code for the C interface to ClassyMC.
!===========================================================================
#define __StdErr__ 0
!===========================================================================
module Library
  use ISO_C_BINDING, only: c_int, c_double, c_ptr, c_f_pointer, c_int, c_null_char, c_loc
  use ParallelVar, only: myid, ierror, nout
  use VarPrecision

  !forcefield push/pull
  logical :: ff_push = .false.
  logical :: ff_pull = .false.

  !sampling push/pull
  logical :: samp_push = .false.
  logical :: samp_pull = .false.

  real(dp), allocatable, private :: boxProb(:) 

!===========================================================================
!  xx_push and xx_pull ->  Logical flags used to communicate with an external computation library such as a Python function
!                          when push is set to true this implies data from the external library has been pushed to Classy
!                          when pull is set to true this imples data is ready from Classy that needs to be pulled by the external
!                          library
!                                                        
!  boxProb -> The probability weight array used to random select a simulation box
!===========================================================================
contains
  !===========================================================================
  subroutine Library_ReadScript(strlen, cfilename) bind(C,name='Classy_ReadScript')
    use C_F_Routines, only: C_F_String
    use ScriptInput, only: Script_ReadParameters
    implicit none
    integer(c_int), intent(in), value :: strlen
    type(c_ptr), value, intent(in) :: cfilename
    character(len=strlen):: infile
    integer :: i, slen

    infile = c_f_string(cfilename)
    call Script_ReadParameters(infile)
 
  end subroutine
  !===========================================================================
  subroutine Library_SetLogFile(strlen, cfilename) bind(C,name='Classy_SetLogFile')
    use C_F_Routines, only: C_F_String
    use ParallelVar, only: nout
    use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_null_char, c_loc
    implicit none
    integer(c_int), intent(in), value :: strlen
    type(c_ptr), value, intent(in) :: cfilename
    character(len=strlen):: outfile
    logical :: isopen
    integer :: i, slen

    inquire(unit=nout, opened=isopen)

    if(isopen) then
      close(nout)
    endif

    outfile = c_f_string(cfilename)

    open(newunit=nout, file=outfile)
 
  end subroutine
  !===========================================================================
  subroutine Library_FullSimulation() bind(C,name='Classy_FullSim')
    use Input_Initialize, only: Script_Initialize
    use SimControl, only: TimeStart, TimeEnd
    use SimMonteCarlo, only: RunMonteCarlo
    implicit none

    call CPU_TIME(TimeStart)
    call Script_Initialize
    call RunMonteCarlo
    call CPU_TIME(TimeEnd)
 
  end subroutine
  !===========================================================================
  subroutine Library_RunPrologue() bind(C,name='Classy_RunPrologue')
    use Input_Initialize, only: Script_Initialize
    use SimMonteCarlo, only: Prologue, SafetyCheck, Analyze
    implicit none
    integer(kind=8) :: iCycle, iMove
    logical :: accept
    iCycle = 0
    iMove = 0

    call Script_Initialize
    call Prologue
    call SafetyCheck
    call Analyze(iCycle, iMove, accept, .true.)
    call Analyze(iCycle, iMove, accept, .false.)
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
    if(.not. allocated(boxProb) ) then
      allocate( boxProb(1:nBoxes) )
    endif
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
  subroutine Library_ForceMove(movenum) bind(C, name='Classy_ForceMove')
  ! ---------------
  ! This subroutine is used to skip over the random selection process and perform a MC move that is specified by the user.
  ! ---------------
    use ISO_C_BINDING, only: c_int
    use BoxData, only: BoxArray
    use CommonSampling, only: Sampling
    use MCMoveData, only: Moves, MoveProb
    use MoveClassDef, only: MCMove
    use SimMonteCarlo, only: Analyze, Update
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use RandomGen, only: sgrnd, grnd, ListRNG
    implicit none
    integer(kind=c_int), intent(in), value :: movenum

    logical :: accept
    integer(kind=8) :: iCycle, iMove
    integer :: boxNum, nBoxes
    real(dp), allocatable :: boxProb(:)


    nBoxes = size(BoxArray)
    if(.not. allocated(boxProb) ) then
      allocate( boxProb(1:nBoxes) )
    endif
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
  subroutine Library_ScreenOut(cyclenum) bind(C,name='Classy_ScreenOut')
    ! ---------------
    ! This subroutine is used to direct Classy to write the standard mid simulation summary to screen. 
    ! ---------------
    use ISO_C_BINDING, only: c_long
    use SimMonteCarlo, only: ScreenOut
    implicit none
    integer(kind=c_long), intent(in), value :: cyclenum
    integer(kind=8) :: iMove, iCycle

    iCycle = cyclenum
    iMove = 0
    call ScreenOut(iCycle, iMove)

 
  end subroutine
!===========================================================================
  function Library_GetAtomCount(boxnum) result(nAtoms) bind(C,name='Classy_GetAtomCount')
    use BoxData, only: BoxArray
    use SimMonteCarlo, only: ScreenOut
    use SimpleSimBox, only: SimpleBox
    implicit none
    integer(kind=c_int), intent(in), value :: boxnum
    integer(kind=c_int) :: nAtoms
    class(SimpleBox), pointer :: simbox => null()
    write(*,*) boxnum

    simbox => BoxArray(boxnum)%box
    nAtoms = simbox%nAtoms
  end function
!===========================================================================
  subroutine Library_GetAtomPos(boxnum, nAtoms, atompos) bind(C,name='Classy_GetAtomPos')
    use BoxData, only: BoxArray
    use SimMonteCarlo, only: ScreenOut
    use SimpleSimBox, only: SimpleBox
    implicit none
    integer(c_int), intent(in), value :: boxnum, nAtoms
    real(c_double), intent(inout) :: atompos(1:3, nAtoms)

    integer :: iAtom
    class(SimpleBox), pointer :: simbox => null()
    real(dp), pointer :: atoms(:,:) => null()
    simbox => BoxArray(boxnum)%box
    call simbox%GetCoordinates(atoms)

    do iAtom = 1, simbox%nMaxAtoms
      if(simbox%IsActive(iAtom)) then
        atompos(1:3,iAtom) = atoms(1:3, iAtom)
      endif
    enddo

 
  end subroutine
!===========================================================================
  subroutine Library_GetAtomTypes(boxnum, nAtoms, outatomtypes) bind(C,name='Classy_GetAtomTypes')
    use BoxData, only: BoxArray
    use SimMonteCarlo, only: ScreenOut
    use SimpleSimBox, only: SimpleBox
    implicit none
    integer(c_int), intent(in), value :: boxnum, nAtoms
    integer(c_int), intent(inout) :: outatomtypes(nAtoms)

    integer :: iAtom
    class(SimpleBox), pointer :: simbox => null()
    integer, pointer :: atomtypes(:) => null()
    simbox => BoxArray(boxnum)%box
    call simbox%GetAtomTypes(atomtypes)

    do iAtom = 1, simbox%nMaxAtoms
      if(simbox%IsActive(iAtom)) then
        outatomtypes(iAtom) = atomtypes(iAtom)
      endif
    enddo

 
  end subroutine

!===========================================================================
end module
!===========================================================================
