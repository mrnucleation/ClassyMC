!=============================================================================+
! Python forcefield interface function. 
!=============================================================================+
module FF_Python
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use Input_Format, only: maxLineLen
  use VarPrecision
  use CoordinateTypes

  !Python Library functions
#ifdef AENET
!  use predict_lib, only: initialize_lib, get_energy_lib
  use aenet, only: aenet_atomic_energy, aenet_init, aenet_load_potential, &
                   aenet_Rc_min, aenet_Rc_max
  use geometry, only: geo_recip_lattice
  use input, only: InputData
#endif

  real(dp), parameter :: boltz = 8.6173303E-5_dp

  type, public, extends(forcefield) :: Pair_Python
!    real(dp) :: rCut, rCutSq
    contains
    !If -DAENET is not passed to the compiler this defaults back to the
    !template's functions which effiectively do nothing.
      procedure, pass :: Constructor => Constructor_Python
      procedure, pass :: DetailedECalc => DetailedECalc_Python
      procedure, pass :: DiffECalc => DiffECalc_Python
      procedure, pass :: VolECalc => VolECalc_Python_V2
      procedure, pass :: ProcessIO => ProcessIO_Python
      procedure, pass :: Prologue => Prologue_Python
!      procedure, pass :: GetCutOff
  end type

  contains
#ifdef AENET
!=============================================================================+
  subroutine Constructor_Python(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes, AtomData
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(Pair_Python), intent(inout) :: self

    type(InputData) :: inp
    character(len=100) :: str1, str2
    integer :: iBox, atomLimit, iType
    integer :: AllocateStat
    character(len=5), allocatable :: symbols(:)

    write(nout, *) "Initializing Python Energy Control"
    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    self%rMin = 0E0_dp
    self%rMinTable = 0E0_dp

  end subroutine
!==========================================================================+
  subroutine DetailedECalc_Python(self, curbox, E_T, accept)
!    use boxData, only: self%boxArray
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(Pair_Python), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    logical :: pbc = .false.
    integer :: i,j
    integer :: iAtom, jAtom
    integer :: atmType1, atmType2
    integer :: nTotalMol, stat
    integer :: nCurAtoms = 0
    real(dp) :: E_Atom
    real(dp) :: rx, ry, rz, rsq, rmin_ij
    real(dp), pointer :: atoms(:,:) => null()

    select type(curbox)
      class is(SimpleBox)
        call curbox%GetCoordinates(atoms)
    end select
    E_T = 0E0_dp
    accept = .true.

  end subroutine
!============================================================================
  subroutine DiffECalc_Python(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    use ClassyConstants, only: pi
    implicit none
    class(Pair_Python), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()


    accept = .true.
    E_Diff = 0E0_dp



    nTotalMol = curbox%nMolTotal
    select type(disp)
      class is(OrthoVolChange)
        call self%VolECalc(curbox, disp, E_Diff, accept)
        return
    end select

    curbox%dETable = 0E0_dp
    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)
    recalcList = 0
    nRecalc = 0

    !Check the rMin criteria first to ensure there is no overlap prior to
    !passing the configuration to Python. In addition create a list of atoms
    !whose interactions will need to be recomputed. 
    select type(disp)
      !-----------------------------------------------------
      class is(Displacement)
        do iDisp = 1, size(disp)
          iAtom = disp(iDisp)%atmIndx
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, nNeigh(iAtom)
            jAtom = neighlist(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)
            rmin_ij = self % rMinTable(atmType2, atmType1)          

            rx = disp(iDisp)%x_new  -  atoms(1, jAtom)
            ry = disp(iDisp)%y_new  -  atoms(2, jAtom)
            rz = disp(iDisp)%z_new  -  atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < rmin_ij) then
                accept = .false.
                return
            endif 
            if(rsq < self%rCutSq) then
              nRecalc = nRecalc + 1
              recalclist(nRecalc) = jAtom
              cycle
            endif

            !If the neighboring atom is not near the new position, check the old
            !one as well.  
            rx = atoms(1, iAtom)  -  atoms(1, jAtom)
            ry = atoms(2, iAtom)  -  atoms(2, jAtom)
            rz = atoms(3, iAtom)  -  atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq) then
              nRecalc = nRecalc + 1
              recalclist(nRecalc) = jAtom
            endif
          enddo
        enddo
      !-----------------------------------------------------
      class is(Addition)
        do iDisp = 1, size(disp)
          iAtom = disp(iDisp)%atmIndx
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, tempNNei(iDisp)
            jAtom = templist(jNei, iDisp)
            atmType2 = curbox % AtomType(jAtom)
            rmin_ij = self % rMinTable(atmType2, atmType1)          

            rx = disp(iDisp)%x_new  -  atoms(1, jAtom)
            ry = disp(iDisp)%y_new  -  atoms(2, jAtom)
            rz = disp(iDisp)%z_new  -  atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < rmin_ij) then
                accept = .false.
                return
            endif 
            if(rsq < self%rCutSq) then
              nRecalc = nRecalc + 1
              recalclist(nRecalc) = jAtom
            endif
          enddo
        enddo
      !-----------------------------------------------------
      class is(Deletion)
        call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)
        do iAtom = molStart, molEnd
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, nNeigh(iAtom)
            jAtom = neighlist(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)
            rx = atoms(1, iAtom)  -  atoms(1, jAtom)
            ry = atoms(2, iAtom)  -  atoms(2, jAtom)
            rz = atoms(3, iAtom)  -  atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq) then
              nRecalc = nRecalc + 1
              recalclist(nRecalc) = jAtom
            endif
          enddo
        enddo
    end select


  end subroutine
!=============================================================================+
  subroutine ProcessIO_Python(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use ParallelVar, only: nout
    implicit none
    class(Pair_Python), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=100) :: command
    character(len=maxLineLen) :: command2
    logical :: param = .false.
    integer :: jType, lineStat, stat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, rCut, rMin
  

    if(.not. allocated(self%inputfiles)) then
      allocate(self%inputfiles(1:nAtomTypes))
    endif

    call GetXCommand(line, command, 1, lineStat)
    select case(trim(adjustl(command)))
      case("network")
        call GetXCommand(line, command, 2, lineStat)
        read(command,*) type1
        call GetXCommand(line, command, 3, lineStat)
!        write(*,*) command
        self%inputfiles = ""
        read(command, *) self%inputfiles(type1)

        write(nout, *) "Loading potential from file: ", trim(adjustl(self%inputfiles(type1)))
        call aenet_load_potential(type1, self%inputfiles(type1), stat)
        self%rCut = aenet_Rc_max
        self%rCutSq = aenet_Rc_max * aenet_Rc_max
        return

      case default
        param = .true.
    end select

    if(param) then
      call CountCommands(line, nPar)
      select case(nPar)
        case(2)
          read(line, *) type1, rMin
          self%rMin(type1) = rMin
          do jType = 1, nAtomTypes
            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case(3)
          write(*,*) line
          read(line, *) type1, type2, rMin
          self%rMinTable(type1, type2) = rMin**2
          self%rMinTable(type2, type1) = rMin**2
        case default
          lineStat = -1
      end select
    endif


  end subroutine
!End Python Safety Block
!=============================================================================+
  subroutine Prologue_Python(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_Python), intent(inout) :: self
    integer :: atomLimit, iBox
    integer :: iType, AllocateStat
      atomLimit = 0
     do iBox = 1, size(BoxArray)
       if( atomLimit < boxArray(iBox) % box % GetMaxAtoms()) then
         atomLimit = boxArray(iBox) % box % GetMaxAtoms()
       endif
     enddo
     allocate(self%atomTypes(1:atomLimit))
     allocate(self%tempcoords(3, 1:atomLimit))

     write(nout, *) "Python Cutoff Distance:", self%rCut

!    if(.not. self%initialized) then
!      do iType = 1, nAtomTypes
!      enddo
!      self%initialized = .true.
!    endif

!    self%rCut = aenet_Rc_max
!    self%rCutSq = aenet_Rc_max * aenet_Rc_max
  end subroutine
!=============================================================================+
end module
!=============================================================================+
