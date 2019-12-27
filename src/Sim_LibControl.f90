
!===========================================================================
#define __StdErr__ 0
!===========================================================================
module Library
  !----------------------------------------------------------------
  !
  !.. module:: Library
  !   :platform: Unix, Windows, Linux
  !   :synopsis: Interface layer used for
  !
  !.. moduleauthor:: Troy D Loeffler <tloeffler@anl.gov>
  !   
  !   This file contains the code for the C interface to ClassyMC.
  !----------------------------------------------------------------
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
    !------------------------------------------------------------------------------------------
    ! Loads a specified script for Classy to read.  This works the exact same
    ! as loading a script via the command line.
    !
    ! C Function Signature
    ! void Classy_ReadScript(int strlen, char *str);
    !------------------------------------------------------------------------------------------
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
    !------------------------------------------------------------------------------------------
    ! Specifies a filename for Classy to write it's screen log to.
    !
    ! C Function Signature
    ! void Classy_SetLogfile(int strlen, char *str);
    !------------------------------------------------------------------------------------------
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
    !------------------------------------------------------------------------------------------
    ! Same as calling the "Run" command from a Classy Script. Performs a full MC simulation
    ! according to the input specified by the input script. 
    !
    ! C Function Signature
    ! void Classy_FullSimulation();
    !------------------------------------------------------------------------------------------
    use Input_Initialize, only: Script_Initialize
    use SimControl, only: TimeStart, TimeEnd
    use SimMonteCarlo, only: RunMonteCarlo
    implicit none

    call CPU_TIME(TimeStart)
!    call Script_Initialize
    call Library_RunPrologue
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
    ! ---------------
    ! This subroutine performed one MC move chosen randomly through Classy's normal internal move selection process.
    !
    ! C Function Signature
    ! void Classy_RunMove();
    ! ---------------
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
    ! This is useful if the user wishes to handle the move selection process outside of Classy.
    !
    ! Input
    !    movenum => The index of a MC move that will be performed.  
    !
    ! C Function Signature
    ! void Classy_ForceMove(int movenum);
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

    if(.not. allocated(Moves)) then

    endif


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
    !
    ! Input
    !    cyclenum => The cycle number that will be printed to the screen along with the box's statistical information.
    !                This is useful if another code is controlling the main simulation loop instead of classy and one
    !                uses to pass the loop index to be printed with the screen information.
    !
    ! C Function Signature
    ! void Classy_ScreenOut(long cyclenum);
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
!    write(*,*) boxnum

    simbox => BoxArray(boxnum)%box
    nAtoms = simbox%nAtoms
  end function
!===========================================================================
  subroutine Library_GetAtomPos(boxnum, nAtoms, atompos) bind(C,name='Classy_GetAtomPos')
    !---------------------------------------------------------
    ! 
    ! Input:
    !    boxnum =>  An integer
    !    nAtoms =>  The number of atoms used to define the size of the C array. Can be determined with the Library_GetAtomCount
    !               function. 
    ! Output:
    !    atompos => A Nx3 C-array containing the (x,y,z) positions of each atom.
    !
    !
    ! C Function Signature
    ! void Classy_GetAtomPos(int boxnum, int nAtoms, double atompos[][3]);
    !---------------------------------------------------------
    use BoxData, only: BoxArray
    use SimMonteCarlo, only: ScreenOut
    use SimpleSimBox, only: SimpleBox
    implicit none
    integer(c_int), intent(in), value :: boxnum, nAtoms
    real(c_double), intent(inout) :: atompos(1:3, nAtoms)

    integer :: iAtom, cnt
    class(SimpleBox), pointer :: simbox => null()
    real(dp), pointer :: atoms(:,:) => null()
    simbox => BoxArray(boxnum)%box
    call simbox%GetCoordinates(atoms)

    cnt = 0
    do iAtom = 1, simbox%nMaxAtoms
      if(simbox%IsActive(iAtom)) then
        cnt = cnt + 1
        atompos(1:3, cnt) = atoms(1:3, iAtom)
      endif
    enddo

 
  end subroutine
!===========================================================================
  subroutine Library_GetAtomTypes(boxnum, nAtoms, outatomtypes, stat) bind(C,name='Classy_GetAtomTypes')
    !---------------------------------------------------------
    ! Returns the atomtype type ids from a specfied simulation box
    ! 
    ! Input:
    !    boxnum =>  An integer that specifies the boxnum that has been requested by the external program
    !    nAtoms =>  The number of atoms used to define the size of the C array. Can be determined with the Library_GetAtomCount
    !               function. 
    ! Output:
    !    outatomtypes => An N-sized integer C-array containing the atomtypes of each atom.
    !
    !
    ! Variables: integer, iAtom => Integer used to loop over the atom array.
    !            pointer, simbox => Pointer used to point toward the box by
    !            pointer, atomtypes => Pointer used to collect the atomtype array from the simulation box
    !
    ! Error Code: stat==-1  >  Invalid Box Number was passed
    !             stat==-2  >  Box Array has not been defined
    !
    ! C Function Signature
    ! void Classy_GetAtomTypes(int boxnum, int nAtoms, int atompos[]);
    !---------------------------------------------------------
    use BoxData, only: BoxArray
    use SimMonteCarlo, only: ScreenOut
    use SimpleSimBox, only: SimpleBox
    implicit none
    integer(c_int), intent(in), value :: boxnum, nAtoms
    integer(c_int), intent(inout) :: stat
    integer(c_int), intent(inout) :: outatomtypes(nAtoms)

    integer :: iAtom, cnt
    class(SimpleBox), pointer :: simbox => null()
    integer, pointer :: atomtypes(:) => null()

 
    !Error Checks to ensure the input is actually valid
    if(.not. allocated(BoxArray)) then
      stat = -2
      return
    endif
    if( (boxnum < 0) .or. (boxnum > size(BoxArray))) then
      stat = -1
      return
    endif

    simbox => BoxArray(boxnum)%box
    call simbox%GetAtomTypes(atomtypes)

    cnt = 0
    do iAtom = 1, simbox%nMaxAtoms
      if(simbox%IsActive(iAtom)) then
        cnt = cnt + 1
        outatomtypes(cnt) = atomtypes(iAtom)
      endif
    enddo
    stat = 0

  end subroutine
!===========================================================================
  subroutine Library_GetNeighList(boxnum, listNum, nAtoms, outlist, stat) bind(C,name='Classy_GetNeighList')
    !---------------------------------------------------------
    ! Returns the atomtype type ids from a specfied simulation box
    ! 
    ! Input:
    !    boxnum =>  An integer that specifies the boxnum that has been requested by the external program
    !    nAtoms =>  The number of atoms used to define the size of the C array. Can be determined with the Library_GetAtomCount
    !               function. 
    ! Output:
    !    outatomtypes => An NxM-sized integer C-array containing the neighborlist for each atom.
    !
    !
    ! Variables: integer, iAtom => Integer used to loop over the atom array.
    !            pointer, simbox => Pointer used to point toward the box by
    !            pointer, atomtypes => Pointer used to collect the atomtype array from the simulation box
    !
    ! Error Code: stat==-1  >  Invalid Box Number was passed
    !             stat==-2  >  Box Array has not been defined
    !
    ! C Function Signature
    ! void Classy_GetNeighList(int boxnum, int nAtoms, int atompos[]);
    !---------------------------------------------------------
    use BoxData, only: BoxArray
    use SimMonteCarlo, only: ScreenOut
    use SimpleSimBox, only: SimpleBox
    implicit none
    integer(c_int), intent(in), value :: boxnum, listnum, nAtoms
    integer(c_int), intent(inout) :: stat
    integer(c_int), intent(inout) :: outlist(1:1000, nAtoms)

    integer :: iAtom, cnt
    class(SimpleBox), pointer :: simbox => null()
    integer, pointer :: nNei(:) => null()
    integer, pointer :: NeighList(:,:) => null()

 
    !Error Checks to ensure the input is actually valid
    if(.not. allocated(BoxArray)) then
      stat = -2
      return
    endif
    if( (boxnum < 0) .or. (boxnum > size(BoxArray))) then
      stat = -1
      return
    endif

    simbox => BoxArray(boxnum)%box
    call simbox%GetNeighborList(listNum, neighlist, nNei)

    cnt = 0
    do iAtom = 1, simbox%nMaxAtoms
      if(simbox%IsActive(iAtom)) then
        cnt = cnt + 1
        outlist(1:nNei(iAtom), cnt) = NeighList(1:nNei(iAtom), iAtom)
      endif
    enddo
    stat = 0

  end subroutine

!===========================================================================
end module
!===========================================================================
