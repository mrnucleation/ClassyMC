!================================================================================
! Easy Pair Template for general pair forcefields. The purpose of this class is
! to act as an inheritable object that can allow the user to quickly set up
! any standard pair distance based forcefields by simply creating a child object
! from this class. By overwriting the PairFunction class and parameter settings
! the procedures for all cut-off based forcefields are automatically defined.
!================================================================================
module FF_EasyEP_LJ_Cut
  use FF_EasyPair_Cut, only: EasyPair_Cut
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(EasyPair_Cut) :: EP_LJ_Cut
!    real(dp), allocatable :: rMin(:)
!    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    real(dp), allocatable :: eps(:)
    real(dp), allocatable :: sig(:)

    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)

                            !Array for book keeping the number of each atom type
                            !is contained within a molecule
    integer, allocatable :: PerMolTypeCount(:, :)
    integer, allocatable :: TypeCount(:)

    contains
      procedure, pass :: Constructor => Constructor_EP_LJ_Cut
      procedure, pass :: PairFunction => PairFunction_EP_LJ_Cut
      procedure, pass :: ProcessIO => ProcessIO_EP_LJ_Cut
      procedure, pass :: TailCorrection => TailCorrection_EP_LJ_Cut
      procedure, pass :: ComputeTypeCount => ComputeTypeCount_EP_LJ_Cut
      procedure, pass :: CountType => CountType_EP_LJ_Cut
      procedure, pass :: SumTail => SumTail_EP_LJ_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EP_LJ_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%eps(1:nAtomTypes), stat = AllocateStat)
    allocate(self%sig(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)

    allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%sigTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%eps = 4E0_dp
    self%sig = 1E0_dp
    self%rMin = 0.5E0_dp

    self%epsTable = 4E0_dp
    self%sigTable = 1E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the LJ/Cut Pair Style"

  end subroutine

  !=============================================================================+
  function PairFunction_EP_LJ_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    real(dp) :: ep, sig_sq, LJ

    ep = self % epsTable(atmType2, atmType1)
    sig_sq = self % sigTable(atmType2, atmType1)
    LJ = (sig_sq/rsq)
    LJ = LJ * LJ * LJ
    E_Pair = ep * LJ * (LJ-1E0_dp)

  end function
  !=============================================================================+
   !Sums up the tail corrections over all atoms in the system.
   ! Uses the following formula for tail corrections. 
   ! u_tail = 8/3 * pi * density * eps * sig^3 * ((1/3) * (sig/rcut)^9 - (sig/rcut)^3 )
   ! This dependent on the output of the CountType function as such that must be ran first.
  function SumTail_EP_LJ_Cut(self) result(E_Corr)
    use Common_MolInfo, only: nAtomTypes
    use ClassyConstants, only: pi
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    real(dp) :: E_Corr
    integer :: iType, jType
    real(dp) :: E_Type, E_TypeTotal
    real(dp) :: ep, sig, LJCube

    E_Corr = 0E0_dp
    do iType = 1, nAtomTypes
       E_Type = 0E0_dp
       E_TypeTotal = 0E0_dp
       do jType = 1, nAtomTypes
         ep = self%epsTable(jType, iType)*0.25E0_dp
         sig = sqrt(self%sigTable(jType, iType))
         LJCube = (sig/self%rCut)**3
         E_Type = (8E0_dp/3E0_dp) * pi * real(self%TypeCount(jType),dp) * ep * sig**3
         E_Type = E_Type*( (1E0_dp/3E0_dp)*LJCube**3 - LJCube)
         E_TypeTotal = E_TypeTotal + E_Type
       enddo
       ! We have the tail energy for a single atom of type i, but now we multiply by the number
       ! of atoms in the system of type i. 
       E_Corr = E_Corr + self%TypeCount(iType)*E_TypeTotal
    enddo
  end function
  !=============================================================================+
  ! Counts the number of atoms of a given type based on the number of molecules in the
  ! simulation box
  subroutine CountType_EP_LJ_Cut(self, NMol) 
    use Common_MolInfo, only: nAtomTypes, nMolTypes
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    integer, intent(in) :: NMol(:)
    integer :: iType, iMolType

    self%TypeCount = 0
    do iType = 1, nAtomTypes
      do iMolType = 1, nMolTypes
         self%TypeCount(iType) = self%TypeCount(iType) + NMol(iMolType)*self%PerMolTypeCount(iType, iMolType)
      enddo
    enddo
  end subroutine
  !======================================================================+
   subroutine TailCorrection_EP_LJ_Cut(self, curbox, E_Corr, disp) 
     use CoordinateTypes, only: Displacement, Addition, Deletion, OrthoVolChange
     use Common_MolInfo, only: nAtomTypes, nMolTypes
     implicit none
     class(EP_LJ_Cut), intent(inout) :: self
     class(Perturbation), intent(in), optional :: disp(:)
     class(SimBox), intent(inout) :: curbox
     real(dp), intent(out) :: E_Corr

     integer :: iType, iMolType, jType, molType
     integer :: NMol(1:nAtomTypes)
     real(dp) :: vol, volnew, E_Type, E_TypeTotal
     character(len=30), parameter :: volume="volume"

     if(.not. allocated(self%PerMolTypeCount)) then
       call self%ComputeTypeCount
     endif

     E_Corr = 0E0_dp
     if(.not. present(disp)) then
       !Tail Calculation for the total box energy. Used at the start
       !and end of simulation or when a total recompute is called.
       vol = curbox%GetThermo_String(volume)
       NMol(1:nMolTypes) = curbox%NMol(1:nMolTypes)
       call self%CountType(NMol)
       E_Corr = self%SumTail()
       E_Corr = E_Corr/vol
     else
       !This section computes the delta of the tail correction
       select type(disp)
         class is(Displacement)
            !For bulk tail corrections, the density of the system is unchanged
            !for a simple displacement move. As such the delta is simply 0
            return

         class is(Addition)
            vol = curbox%GetThermo_String(volume)
            molType = disp(1)%molType

            NMol(1:nMolTypes) = curbox%NMol(1:nMolTypes)
            call self%CountType(NMol)
            E_Corr = self%SumTail()/vol

            NMol(molType) = NMol(molType) + 1
            call self%CountType(NMol)
            E_Corr = self%SumTail()/vol - E_Corr

         class is(Deletion)
            vol = curbox%GetThermo_String(volume)
            molType = disp(1)%molType

            NMol(1:nMolTypes) = curbox%NMol(1:nMolTypes)
            call self%CountType(NMol)
            E_Corr = self%SumTail()/vol

            NMol(molType) = NMol(molType) - 1
            call self%CountType(NMol)
            E_Corr = self%SumTail()/vol - E_Corr

         class is(OrthoVolChange)
            vol = curbox%GetThermo_String(volume)
            volnew = disp(1)%volnew
            NMol(1:nMolTypes) = curbox%NMol(1:nMolTypes)
            call self%CountType(NMol)
            E_Corr = self%SumTail()
            E_Corr = E_Corr*(1E0_dp/volnew)

         class default
           write(0,*) "LJ Tail Corrections Currently not supported for this move type."
           error stop

       end select
     endif

  end subroutine
  !=============================================================================+
  ! Gathers information about the number of each atom type per molecule.  Primarily used in
  ! tail correction calculations. 
  subroutine ComputeTypeCount_EP_LJ_Cut(self)
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    integer :: AllocateStat
    integer :: iMol, iAtom, atmType 

    if(allocated(self%PerMolTypeCount) .or. allocated(self%TypeCount)) then
      error stop "TypeCount arrays are already allocated! This function can't be called more than once!"
    endif

    allocate(self%PerMolTypeCount(1:nAtomTypes, 1:nMolTypes), stat=AllocateStat )
    allocate(self%TypeCount(1:nAtomTypes), stat=AllocateStat )

    if(AllocateStat /= 0) error stop "Memory Allocation Error Encountered!"
    self%PerMolTypeCount = 0
    do iMol = 1, nMolTypes
      do iAtom = 1, MolData(iMol)%nAtoms
        atmType = MolData(iMol)%atomType(iAtom)
        self%PerMolTypeCount(atmType, iMol) = self%PerMolTypeCount(atmType, iMol) + 1
      enddo
    enddo

  end subroutine
  !=====================================================================
  subroutine ProcessIO_EP_LJ_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(EP_LJ_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: logicval
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, rCut, rMin
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("tailcorrection")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) logicval
        self%usetailcorrection = logicval

      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut

      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      call CountCommands(line, nPar)
      select case(nPar)
        case(4)
          read(line, *) type1, ep, sig, rMin
          ep = ep * inEngUnit
          sig = sig * inLenUnit
          rMin = rMin * inLenUnit

          self%eps(type1) = ep 
          self%sig(type1) = sig 
          self%rMin(type1) = rMin 


          do jType = 1, nAtomTypes
            self%epsTable(type1, jType) = 4E0_dp * sqrt(ep * self%eps(jType))
            self%epsTable(jType, type1) = 4E0_dp * sqrt(ep * self%eps(jType))

            self%sigTable(type1, jType) = (0.5E0_dp * (sig + self%sig(jType)))**2
            self%sigTable(jType, type1) = (0.5E0_dp * (sig + self%sig(jType)))**2

            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case(5)
          read(line, *) type1, type2, ep, sig, rMin
          self%epsTable(type1, type2) = ep
          self%epsTable(type2, type1) = ep

          self%sigTable(type1, type2) = sig
          self%sigTable(type2, type1) = sig

          self%rMinTable(type1, type2) = rMin**2
          self%rMinTable(type2, type1) = rMin**2


        case default
          lineStat = -1
          write(0,*) line
          stop "Unknown Command Found"
      end select
    endif

  end subroutine
  !=====================================================================
end module
!=====================================================================
