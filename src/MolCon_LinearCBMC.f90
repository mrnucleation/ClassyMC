!==========================================================================================
! Simple Regrowth Object
!==========================================================================================
module MolCon_LinearCBMC
  use CoordinateTypes, only: Perturbation, Addition
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: LinearCBMC
    integer :: firstAtom = 1
    integer :: rosenbluth = 1
    real(dp), allocatable :: tempcoords(:, :)
    integer, allocatable :: patharray(:)
    integer, allocatable :: pathposition(:)
    integer, allocatable :: grown(:)
    integer, allocatable :: schedule(:)
    contains
!      procedure, public, pass :: Constructor => LinearCBMC_Constructor
      procedure, public, pass :: Prologue => LinearCBMC_Prologue
      procedure, public, pass :: GenerateConfig => LinearCBMC_GenerateConfig
      procedure, public, pass :: ReverseConfig => LinearCBMC_ReverseConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine LinearCBMC_Prologue(self)
    use Common_MolInfo, only: MolData, BondData, nMolTypes
    use MolSearch, only: FindBond
    use ParallelVar, only: nout
    implicit none
    class(LinearCBMC), intent(inout) :: self
!    integer, intent(in) :: molType
    integer :: iType, iBond, iAtom, curMax, nAtoms
    integer :: atm1, atm2
    integer, allocatable :: freq(:)
    integer :: iError = 0
    integer :: nextAtm, prevAtm, curAtm, iPath

    !Count the number of bonds each atom has.  This will be used to classify the atom
    !as either a Terminal, Linker, or Branch atom. 
    ! nBonds = 0 >> Isolated Atom
    ! nBonds = 1 >> Terminal Atom located on the end of a chain
    ! nBonds = 2 >> Linker Atom.  Located in the middle of a linear chain
    ! nBonds >= 3 >> Branch Atom.  Has two or more potential paths one can wander down.
    nAtoms = MolData(self%molType)%nAtoms
    if(MolData(self%molType)%nAtoms > 1) then
      allocate( freq(1:nAtoms) )  
      freq = 0
      do iBond = 1, size(MolData(self%molType)%bond)
        atm1 = MolData(self%molType)%bond(iBond)%mem1
        freq(atm1) = freq(atm1) + 1
        atm2 = MolData(self%molType)%bond(iBond)%mem2
        freq(atm2) = freq(atm2) + 1
      enddo
    else 
      iError = -1
    endif


    !If any atom has more than 3 bonds it's a branched molecule and not suited for this algorithm.
    if( any(freq > 2) ) then
      iError = -1
    endif

    !If there's no terminal atoms the molecule is likely cyclic
    if( all(freq /= 1) ) then
      iError = -1
    endif


    !Since this module is for linear molecules, if a branch is detected or
    !the atom is simply an isolated atom throw an error since it is not compatible with.
    !this module
    if(iError /= 0) then
      write(0,*) "Linear CBMC has been invoked on a molecule that is not linear. Stopping Classy."
      stop
    endif


    !Create the path array which is used to determine regrowth order.
    allocate( self%patharray(1:nAtoms) )  
    allocate( self%pathposition(1:nAtoms) )  
    self%patharray = 0
    self%pathposition = 0
    !Pick an terminal atom to start building the patharray
    prevAtm = 0
    do iAtom = 1,nAtoms
      if(freq(iAtom) == 1) then
        curatm = iAtom
        exit
      endif
    enddo

    self%patharray(1) = curatm
    self%pathposition(curatm) = 1
    iPath = 1
    ! Begin the building process by traversing from bond to bond.
    do iAtom = 1, nAtoms-1
      nextatm = 0
      do iBond = 1, size(MolData(self%molType)%bond)
        atm1 = MolData(self%molType)%bond(iBond)%mem1
        atm2 = MolData(self%molType)%bond(iBond)%mem2
        if(atm1 == curatm) then
          if(atm2 /= prevatm) then
            nextatm = atm2
            exit
          endif
        endif
        if(atm2 == curatm) then
          if(atm1 /= prevatm) then
            nextatm = atm1
            exit
          endif
        endif
      enddo
!      write(*,*) prevatm, curatm, nextatm
      if(nextatm /= 0) then
         iPath = iPath + 1
         self%patharray(iPath) = nextatm
         self%pathposition(nextatm) = iPath
         prevatm = curatm
         curatm = nextatm
      else
        write(0,*) "Linear CBMC: DEAD END ERROR! There is a problem in the molecule definition!"
        write(0,*) "Ensure your molecule's bonds are properly connected."
        stop
      endif
    enddo


    if(any(self%patharray == 0)) then
      write(0,*) "LINEAR CBMC: ERROR! There are atoms which are unaccounted for by the path building algorithm."
      write(0,*) "Ensure your molecule's bonds are properly connected."
      stop
    endif


!    write(*,*) self%patharray
!    write(*,*) self%pathposition




  end subroutine
!==========================================================================================
  subroutine LinearCBMC_GenerateConfig(self, trialBox, disp, probconstruct, insPoint)
    use Common_MolInfo, only: MolData, BondData, AngleData, nMolTypes
    use MolSearch, only: FindBond, FindAngle
    use RandomGen, only: Generate_UnitSphere, Generate_UnitCone
    implicit none
    class(LinearCBMC), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)

    integer :: bondType, angleType, molType
    integer :: atm1, atm2,atm3, iDisp, iAtom
    real(dp), intent(out) :: probconstruct
    real(dp), dimension(1:3) :: v1, v2, v3
    real(dp) :: dx, dy, dz, r, theta
    real(dp) :: r1, r2
    real(dp) :: prob
    real(dp) :: ang1, ang2

    probconstruct = 1E0_dp
    select type(disp)
      class is(Addition)
        molType = disp(1)%molType
      class default
        stop "Critical Errror! An invalid perturbation type has been passed into the regrowth function"
    end select



  end subroutine
!======================================================================================
  subroutine LinearCBMC_ReverseConfig(self, trialBox, probconstruct, accept)
    implicit none
    class(LinearCBMC), intent(inout) :: self
!    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct 
    logical, intent(out) :: accept

    accept = .true.
    probconstruct = 1E0_dp
  end subroutine
!==========================================================================================
  subroutine LinearCBMC_GasConfig(self,  probGas)
    implicit none
    class(MolConstructor), intent(inout) :: self
    real(dp), intent(out) :: probGas


    probGas = 1E0_dp
  end subroutine
!==========================================================================================
end module
!==========================================================================================
