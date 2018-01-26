!================================================================================
module FF_Pair_LJ_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_LJ_Cut
    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor => Constructor_LJ_Cut
      procedure, pass :: DetailedECalc => Detailed_LJ_Cut
      procedure, pass :: ShiftECalc_Single => Shift_LJ_Cut_Single
      procedure, pass :: ShiftECalc_Multi => Shift_LJ_Cut_Multi
      procedure, pass :: NewECalc => New_LJ_Cut
      procedure, pass :: OldECalc => Old_LJ_Cut
      procedure, pass :: ProcessIO => ProcessIO_LJ_Cut
      procedure, pass :: GetCutOff => GetCutOff_LJ_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_LJ_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%sigTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%epsTable = 4E0_dp
    self%sigTable = 1E0_dp
    self%rMinTable = 0.05E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_LJ_Cut(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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
        atmType2 = curbox % AtomType(jAtom)
        ep = self % epsTable(atmType1,atmType2)
        sig_sq = self % sigTable(atmType1,atmType2)          
        rmin_ij = self % rMinTable(atmType1,atmType2)          

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx**2 + ry**2 + rz**2
        if(rsq < self%rCutSq) then
          if(rsq < rmin_ij) then
            write(*,*) sqrt(rsq)
            write(*,*) iAtom, jAtom
            write(*,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
            write(*,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
            stop "ERROR! Overlaping atoms found in the current configuration!"
          endif 
          LJ = (sig_sq/rsq)**3
          LJ = ep * LJ * (LJ-1E0)              
          E_LJ = E_LJ + LJ
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + LJ
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + LJ 
        endif
      enddo
    enddo
      
    write(nout,*) "Lennard-Jones Energy:", E_LJ
      
    E_T = E_LJ    
  end subroutine
  !=====================================================================
  subroutine Shift_LJ_Cut_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
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
    curbox%dETable = 0E0_dp
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)
        ep = self % epsTable(atmType1, atmType2)
        sig_sq = self % sigTable(atmType1, atmType2)          
        rmin_ij = self % rMinTable(atmType1, atmType2)          

        rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif 
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
  subroutine Shift_LJ_Cut_Multi(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_LJ_Cut), intent(inout) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
   
  end subroutine
  !=====================================================================
  subroutine New_LJ_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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

        atmType2 = curbox % AtomType(jAtom)
        ep = self%epsTable(atmType1,atmType2)
        sig_sq = self%sigTable(atmType1,atmType2)          
        rmin_ij = self%rMinTable(atmType1,atmType2)          

        rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
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
  subroutine Old_LJ_Cut(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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
  subroutine ProcessIO_LJ_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetAllCommands, GetXCommand
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    character(len=*), intent(in) :: line
    character(len=30), allocatable :: parlist(:)
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2
    real(dp) :: ep, sig, rCut
  

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
      call GetAllCommands(line, parlist, lineStat)
      select case(size(parlist))
        case(3)
          read(line, *) type1, ep, sig
          do jType = 1, nAtomTypes
            if(jType == type1) then
              self%epsTable(type1, jType) = 4E0_dp * ep
              self%sigTable(type1, jType) = sig
            else
              self%epsTable(type1, jType) = 4E0_dp * sqrt(ep * self%epsTable(jType, jType))
              self%epsTable(jType, type1) = 4E0_dp * sqrt(ep * self%epsTable(jType, jType))

              self%sigTable(type1, jType) = 0.5E0_dp * (sig + self%sigTable(jType, jType) )
              self%sigTable(jType, type1) = 0.5E0_dp * (sig + self%sigTable(jType, jType) )
            endif
          enddo
        case(4)
          read(line, *) type1, type2, ep, sig
          self%epsTable(type1, type2) = 4E0_dp * ep
          self%epsTable(type2, type1) = 4E0_dp * ep

          self%sigTable(type1, type2) = sig
          self%sigTable(type2, type1) = sig

        case default
          lineStat = -1
      end select
      if( allocated(parlist) ) then 
        deallocate(parlist)
      endif
    endif


!    deallocate(parlist)
  end subroutine
  !=============================================================================+
  function GetCutOff_LJ_Cut(self) result(rCut)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
  !=====================================================================
end module
!=====================================================================
