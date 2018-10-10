!================================================================================
! Standard 12-6-1 Lennard-Jones with Eletrostatic pairwise forcefield.
!================================================================================
module FF_Pair_LJ_Q_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use Constants, only: coulombConst

  type, extends(forcefield) :: Pair_LJ_Q_Cut
    real(dp), allocatable :: eps(:)
!    real(dp), parameter :: coul = 1.671009770E5_dp
    real(dp), allocatable :: sig(:)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: rMin(:)

    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: qTable(:,:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    real(dp) :: rLJCut, rLJCutSq
    real(dp) :: rQCut, rQCutSq
    contains
      procedure, pass :: Constructor => Constructor_LJ_Q_Cut
      procedure, pass :: DetailedECalc => Detailed_LJ_Q_Cut
      procedure, pass :: DiffECalc => DiffECalc_LJ_Q_Cut
      procedure, pass :: ShiftECalc_Single => Shift_LJ_Q_Cut_Single
      procedure, pass :: ShiftECalc_Multi => Shift_LJ_Q_Cut_Multi
      procedure, pass :: NewECalc => New_LJ_Q_Cut
      procedure, pass :: OldECalc => Old_LJ_Q_Cut
      procedure, pass :: ProcessIO => ProcessIO_LJ_Q_Cut
      procedure, pass :: Prologue => Prologue_LJ_Q_Cut
      procedure, pass :: GetCutOff => GetCutOff_LJ_Q_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_LJ_Q_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%eps(1:nAtomTypes), stat = AllocateStat)
    allocate(self%sig(1:nAtomTypes), stat = AllocateStat)
    allocate(self%q(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)

    allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%sigTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%qTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%eps = 4E0_dp
    self%sig = 1E0_dp
    self%q = 0E0_dp
    self%rMin = 0.5E0_dp

    self%epsTable = 4E0_dp
    self%sigTable = 1E0_dp
    self%qTable = 0E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2
    self%rLJCut = 5E0_dp
    self%rLJCutSq = 5E0_dp**2
    self%rQCut = 5E0_dp
    self%rQCutSq = 5E0_dp**2
    IF (AllocateStat /= 0) STOP "Allocation in the LJ/Cut Pair Style"

  end subroutine
!============================================================================
  subroutine DiffECalc_LJ_Q_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

    select type(disp)
      class is(DisplacementNew)
         call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)

      class is(Addition)
!         write(*,*) size(tempList)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

!      class is(Displacement)
!        stop
      class default
        write(*,*) "Unknown Perturbation Type Encountered by the LJ_Q_Cut Pair Style."
    end select


  end subroutine
  !===================================================================================
  subroutine Detailed_LJ_Q_Cut(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: ep, sig_sq, q
    real(dp) :: LJ, Ele
    real(dp) :: E_LJ, E_Ele
    real(dp) :: rmin_ij      

    E_LJ = 0E0_dp
    E_Ele = 0E0_dp
    curbox%ETable = 0E0_dp
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
        if( curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)  ) then
          cycle
        endif

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx**2 + ry**2 + rz**2
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self % rMinTable(atmType1, atmType2)          
          if(rsq < rmin_ij) then
            write(*,*) sqrt(rsq)
            write(*,*) iAtom, jAtom
            write(*,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
            write(*,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
            write(*,*) "ERROR! Overlaping atoms found in the current configuration!"
          endif 

          if(rsq < self%rLJCutSq) then
            ep = self % epsTable(atmType1, atmType2)
            sig_sq = self % sigTable(atmType1, atmType2)          
            write(*,*) ep, sig_sq, rsq
            LJ = (sig_sq/rsq)**3
            LJ = ep * LJ * (LJ-1E0_dp)
            write(*,*) LJ
            E_LJ = E_LJ + LJ
            curbox%ETable(iAtom) = curbox%ETable(iAtom) + LJ
            curbox%ETable(jAtom) = curbox%ETable(jAtom) + LJ 
          endif

          if(rsq < self%rQCutSq) then
            q = self % qTable(atmType1, atmType2)
            r = sqrt(rsq)
            Ele = q/r
            E_Ele = E_Ele + Ele
            curbox%ETable(iAtom) = curbox%ETable(iAtom) + Ele
            curbox%ETable(jAtom) = curbox%ETable(jAtom) + Ele 
          endif
      
        endif
      enddo
    enddo
  
    write(nout,*) "Lennard-Jones Energy:", E_LJ
    write(nout,*) "Eletrostatic Energy:", E_Ele
      
    E_T = E_LJ + E_Ele
  end subroutine
  !=====================================================================
  subroutine Shift_LJ_Q_Cut_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(DisplacementNew), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: ep, sig_sq, q
    real(dp) :: LJ, Ele
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)

!      write(*,*) iAtom, curbox%NeighList(1)%nNeigh(iAtom)
!      write(*,*) iAtom, curbox%NeighList(1)%list(:, iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        atmType2 = curbox % AtomType(jAtom)
        ep = self % epsTable(atmType2, atmType1)
        sig_sq = self % sigTable(atmType2, atmType1)  
        q = self % qTable(atmType2, atmType1)

        rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) sqrt(rsq)
        if(rsq < self%rCutSq) then
          rmin_ij = self % rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif 

          if(rsq < self%rLJCutSq) then
            if(ep /= 0E0_dp) then
              LJ = (sig_sq/rsq)
              LJ = LJ * LJ * LJ
              LJ = ep * LJ * (LJ-1E0_dp)
              E_Diff = E_Diff + LJ
            endif
          endif

          if(rsq < self%rQCutSq) then
            r = sqrt(rsq)
            Ele = q/r
            E_Diff = E_Diff + Ele
            curbox % dETable(iAtom) = curbox % dETable(iAtom) + Ele
            curbox % dETable(jAtom) = curbox % dETable(jAtom) + Ele
          endif
        endif

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) sqrt(rsq)
        if(rsq < self%rCutSq) then
!          if(iAtom == jAtom) then
!            write(*,*) "NeighborList Error!", iAtom, jAtom
!          endif
          if(rsq < self%rLJCutSq) then
            if(ep /= 0E0_dp) then
              sig_sq = self % sigTable(atmType2, atmType1)  
              LJ = (sig_sq/rsq)
              LJ = LJ * LJ * LJ
              LJ = ep * LJ * (LJ-1E0_dp)              
              E_Diff = E_Diff - LJ
              curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
              curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
            endif
          endif

          if(rsq < self%rQCutSq) then
            if(q /= 0E0_dp) then
              r = sqrt(rsq)
              Ele = q/r
              E_Diff = E_Diff - Ele
              curbox % dETable(iAtom) = curbox % dETable(iAtom) - Ele
              curbox % dETable(jAtom) = curbox % dETable(jAtom) - Ele
            endif
          endif
        endif

      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine Shift_LJ_Q_Cut_Multi(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_LJ_Q_Cut), intent(inout) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
   
  end subroutine
  !=====================================================================
  subroutine New_LJ_Q_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: ep, sig_sq, q 
    real(dp) :: LJ, Ele
    real(dp) :: rmin_ij

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.

!    write(*,*) "Length:", size(tempNNei)
!    write(*,*) "Length:", size(tempList)
!    write(*,*)
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
!      write(*,*) iAtom
      atmType1 = curbox % AtomType(iAtom)

      listIndx = disp(iDisp)%listIndex
      maxNei = tempNNei(listIndx)
      do jNei = 1, maxNei
        jAtom = tempList(jNei, listIndx)
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

          if(rsq < self%rLJCutSq) then
            ep = self%epsTable(atmType2, atmType1)
            if(ep /= 0E0_dp) then
              sig_sq = self%sigTable(atmType2, atmType1)          
              LJ = (sig_sq/rsq)
              LJ = LJ * LJ * LJ
              LJ = ep * LJ * (LJ-1E0_dp)
              E_Diff = E_Diff + LJ
              curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
              curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ
            endif
          endif

          if(rsq < self%rQCutSq) then
            q = self%qTable(atmType2, atmType1)
            if(q /= 0E0_dp) then
              r = sqrt(rsq)
              Ele = q/r
              E_Diff = E_Diff + Ele
              curbox % dETable(iAtom) = curbox % dETable(iAtom) + Ele
              curbox % dETable(jAtom) = curbox % dETable(jAtom) + Ele
            endif
          endif
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_LJ_Q_Cut(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2, globIndx
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: ep, sig_sq, q
    real(dp) :: LJ, Ele
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

    globIndx = curBox % MolGlobalIndx(disp(1)%molType, disp(1)%molIndx)
    call curBox % GetMolData(globIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
      atmType1 = curbox % AtomType(iAtom) 
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)
!        rmin_ij = self % rMinTable(atmType1, atmType2)          

        rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then

          if(rsq < self%rLJCutSq) then
            ep = self % epsTable(atmType2, atmType1)
            if(ep /= 0E0_dp) then
              sig_sq = self % sigTable(atmType2, atmType1)
              LJ = (sig_sq/rsq)
              LJ = LJ * LJ * LJ
              LJ = ep * LJ * (LJ-1E0_dp)
              E_Diff = E_Diff - LJ
              curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
              curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
            endif
          endif

          if(rsq < self%rQCutSq) then
            q = self%qTable(atmType2, atmType1)
            if(q /= 0E0_dp) then
              r = sqrt(rsq)
              Ele = q/r
              E_Diff = E_Diff - Ele
              curbox % dETable(iAtom) = curbox % dETable(iAtom) - Ele
              curbox % dETable(jAtom) = curbox % dETable(jAtom) - Ele
            endif
          endif
        endif
      enddo
    enddo
!    write(*,*) E_Diff
  end subroutine
  !=====================================================================
  subroutine ProcessIO_LJ_Q_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, q, rCut, rMin
  

    call GetXCommand(line, command, 1, lineStat)

    call CountCommands(line, nPar)
    select case(trim(adjustl(command)))
      case("ljrcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rLJCut = rCut
        self % rLJCutSq = rCut*rCut
        self % rCut = max(self%rLJCut, self%rQCut)
        self % rCutSq = self % rCut * self % rCut

      case("qrcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rQCut = rCut
        self % rQCutSq = rCut*rCut
        self % rCut = max(self%rLJCut, self%rQCut)
        self % rCutSq = self % rCut * self % rCut

!      case("lj")
!        select case(nPar)
!          case(4)
!            self%eps(type1) = ep
!            self%sig(type1) = sig
!            self%rMin(type1) = rMin
!          case default
!            lineStat = -1
!            return
!        end select

      case("charge")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1

        call GetXCommand(line, command, 3, lineStat)
        read(command, *) q
        self%q(type1) = q
        do jType = 1, nAtomTypes
          self%qTable(type1, jType) = q * self%q(jType) * coulombConst
          self%qTable(jType, type1) = q * self%q(jType) * coulombConst
        enddo


      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      select case(nPar)
        case(5)
          read(line, *) type1, ep, sig, q, rMin

          self%eps(type1) = ep
          self%sig(type1) = sig
          self%q(type1) = q
          self%rMin(type1) = rMin

          do jType = 1, nAtomTypes
            self%epsTable(type1, jType) = 4E0_dp * sqrt(ep * self%eps(jType))
            self%epsTable(jType, type1) = 4E0_dp * sqrt(ep * self%eps(jType))

            self%sigTable(type1, jType) = (0.5E0_dp * (sig + self%sig(jType)))**2
            self%sigTable(jType, type1) = (0.5E0_dp * (sig + self%sig(jType)))**2

            self%qTable(type1, jType) = q * self%q(jType) * coulombConst
            self%qTable(jType, type1) = q * self%q(jType) * coulombConst

            self%rMinTable(type1, jType) = max(rMin, self%rMin(jType))**2
            self%rMinTable(jType, type1) = max(rMin, self%rMin(jType))**2
          enddo
        case default
          lineStat = -1
      end select
    endif

  end subroutine
  !=============================================================================+
  function GetCutOff_LJ_Q_Cut(self) result(rCut)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
  !=====================================================================
  subroutine Prologue_LJ_Q_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    integer :: i, j
    
    
!    write(nout,*) self%rLJCut, self%rQCut
!    write(nout,*) self%rLJCutSq, self%rQCutSq
!    write(nout,*)

!    do i = 1, nAtomTypes
!      write(nout, *) (self%epsTable(i,j), j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%sigTable(i,j), j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%rMinTable(i,j), j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%qTable(i,j)/coulombConst, j=1,nAtomTypes)
!    enddo

!    write(nout,*)
!    do i = 1, nAtomTypes
!      write(nout, *) (self%qTable(i,j), j=1,nAtomTypes)
!    enddo


  end subroutine
  !=====================================================================
end module
!=====================================================================
