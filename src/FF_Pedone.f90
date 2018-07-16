!================================================================================
module FF_Pair_Pedone
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_Pedone
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    type AtomDefPedone
      character(len=20) :: atmName
      character(len=5) :: Symb
      real(dp) :: repul, rEq, q, alpha, delta, mass
    end type


    type(AtomDefPedone), allocatable :: pedoneData(:)
    real(dp), allocatable :: repul_tab(:,:), rEq_tab(:,:)
    real(dp), allocatable :: q_tab(:,:), alpha_Tab(:,:), D_Tab(:,:)

    contains
      procedure, pass :: Constructor => Constructor_Pedone
      procedure, pass :: DetailedECalc => Detailed_Pedone
      procedure, pass :: ShiftECalc_Single => Shift_Pedone_Single
      procedure, pass :: ShiftECalc_Multi => Shift_Pedone_Multi
      procedure, pass :: NewECalc => New_Pedone
      procedure, pass :: OldECalc => Old_Pedone
      procedure, pass :: ProcessIO => ProcessIO_Pedone
      procedure, pass :: Prologue => Prologue_Pedone
      procedure, pass :: GetCutOff => GetCutOff_Pedone
  end type

  contains
  !=============================================================================+
  subroutine Constructor_Pedone(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the Pedone Pair Style"

  end subroutine
  !===================================================================================
  subroutine Detailed_Pedone(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: E_LJ
    real(dp) :: rmin_ij      

    E_LJ = 0E0
    curbox%ETable = 0E0
    accept = .true.
    do iAtom = 1, curbox%nMaxAtoms-1
      atmType1 = curbox % AtomType(iAtom)
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
        cycle
      endif
      do jAtom = iAtom+1, curbox%nMaxAtoms
        if( curbox%MolSubIndx(jAtom) > curbox%NMol(curbox%MolType(jAtom)) ) then
          cycle
        endif
        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx**2 + ry**2 + rz**2
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          ep = self % epsTable(atmType1, atmType2)
          sig_sq = self % sigTable(atmType1, atmType2)          
          rmin_ij = self % rMinTable(atmType1, atmType2)          
          if(rsq < rmin_ij) then
            write(*,*) sqrt(rsq)
            write(*,*) iAtom, jAtom
            write(*,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
            write(*,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
            write(*,*) "ERROR! Overlaping atoms found in the current configuration!"
          endif 
          LJ = (sig_sq/rsq)**3
          LJ = ep * LJ * (LJ-1E0)              
          E_LJ = E_LJ + LJ
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + LJ
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + LJ 
        endif
      enddo
    enddo
  
!    write(nout,*) "Lennard-Jones Energy:", E_LJ
      
    E_T = E_LJ    
  end subroutine
  !=====================================================================
  subroutine Shift_Pedone_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(DisplacementNew), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        atmType2 = curbox % AtomType(jAtom)
        if(rsq < self%rCutSq) then
          rmin_ij = self % rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif 
          ep = self % epsTable(atmType2, atmType1)
          sig_sq = self % sigTable(atmType2, atmType1)  

          LJ = (sig_sq/rsq)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0_dp)              
          E_Diff = E_Diff + LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ
        endif

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          ep = self % epsTable(atmType2, atmType1)
          sig_sq = self % sigTable(atmType2, atmType1)  
          LJ = (sig_sq/rsq)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0_dp)              
          E_Diff = E_Diff - LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
        endif

      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine Shift_Pedone_Multi(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_Pedone), intent(inout) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
   
  end subroutine
  !=====================================================================
  subroutine New_Pedone(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.

    do iDisp = 1, dispLen
      if(.not. disp(iDisp)%newAtom) then
        cycle
      endif
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      listIndx = disp(iDisp)%listIndex
      if(disp(iDisp)%newlist) then
        maxNei = tempNNei(listIndx)
      else
        maxNei = curbox%NeighList(1)%nNeigh(listIndx)
      endif

      do jNei = 1, maxNei
        if(disp(iDisp)%newlist) then
          jAtom = tempList(jNei, listIndx)
        else
          jAtom = curbox%NeighList(1)%list(jNei, listIndx)
        endif
        if( any(jAtom == disp(:)%atmIndx) ) then
          cycle
        endif

        rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self%rMinTable(atmType2, atmType1)          
        
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
          ep = self%epsTable(atmType2, atmType1)
          sig_sq = self%sigTable(atmType2, atmType1)          

          LJ = (sig_sq/rsq)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0_dp)
          E_Diff = E_Diff + LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_Pedone(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp
    do iDisp = 1, size(disp)
      if(.not. disp(iDisp)%oldAtom) then
        cycle
      endif
      iAtom = disp(iDisp)%oldAtmIndx

      atmType1 = curbox % AtomType(iAtom) 
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)
        ep = self % epsTable(atmType1, atmType2)
        sig_sq = self % sigTable(atmType1, atmType2)          
!        rmin_ij = self % rMinTable(atmType1, atmType2)          

        rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          LJ = (sig_sq/rsq)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0_dp)
          E_Diff = E_Diff - LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine ProcessIO_Pedone(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, rCut, rMin
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
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

          self%eps(type1) = ep
          self%sig(type1) = sig
          self%rMin(type1) = rMin

          do jType = 1, nAtomTypes
            self%epsTable(type1, jType) = 4E0_dp * sqrt(ep * self%eps(jType))
            self%epsTable(jType, type1) = 4E0_dp * sqrt(ep * self%eps(jType))

            self%sigTable(type1, jType) = 0.5E0_dp * (sig + self%sig(jType) )
            self%sigTable(jType, type1) = 0.5E0_dp * (sig + self%sig(jType) )

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
      end select
    endif

  end subroutine
  !=============================================================================+
  function GetCutOff_Pedone(self) result(rCut)
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function

  !=====================================================================
  subroutine Prologue_Pedone(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_Pedone), intent(inout) :: self
    integer :: i, j

    do i = 1, nAtomTypes
      write(nout, *) (self%epsTable(i,j), j=1,nAtomTypes)
    enddo

    do i = 1, nAtomTypes
      write(nout, *) (self%sigTable(i,j), j=1,nAtomTypes)
    enddo

    do i = 1, nAtomTypes
      write(nout, *) (self%rMinTable(i,j), j=1,nAtomTypes)
    enddo



  end subroutine
  !=====================================================================
end module
!=====================================================================