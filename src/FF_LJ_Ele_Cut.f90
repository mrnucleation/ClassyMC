!================================================================================
! Standard 12-6-1 Lennard-Jones with Eletrostatic pairwise forcefield.
!================================================================================
module FF_Pair_LJ_Q_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes
  use ClassyConstants, only: coulombConst

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
!      procedure, pass :: ShiftECalc_Multi => Shift_LJ_Q_Cut_Multi
      procedure, pass :: NewECalc => New_LJ_Q_Cut
      procedure, pass :: OldECalc => Old_LJ_Q_Cut
      procedure, pass :: ManyBody => ManyBody_LJ_Q_Cut
      procedure, pass :: AtomExchange => AtomExchange_LJ_Q_Cut
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
    IF (AllocateStat /= 0) error STOP "Allocation in the LJ/Cut Pair Style"

  end subroutine
!============================================================================
  subroutine DiffECalc_LJ_Q_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

    select type(disp)
      class is(Displacement)
         call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)

      class is(Addition)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

      class is(AtomExchange)
         call self % AtomExchange( curbox, disp, E_Diff, accept)

      class default
        write(0,*) "Unknown Perturbation Type Encountered by the LJ_Q_Cut Pair Style."
        error stop
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
            write(0,*) sqrt(rsq)
            write(0,*) iAtom, jAtom
            write(0,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
            write(0,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
            write(0,*) "ERROR! Overlaping atoms found in the current configuration!"
          endif 

          if(rsq < self%rLJCutSq) then
            ep = self % epsTable(atmType1, atmType2)
            sig_sq = self % sigTable(atmType1, atmType2)          
            LJ = (sig_sq/rsq)**3
            LJ = ep * LJ * (LJ-1E0_dp)
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
    type(Displacement), intent(in) :: disp(:)
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
              curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
              curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ
            endif
          endif

          if(rsq < self%rQCutSq) then
            if(q /= 0E0_dp) then
              r = sqrt(rsq)
              Ele = q/r
              E_Diff = E_Diff + Ele
              curbox % dETable(iAtom) = curbox % dETable(iAtom) + Ele
              curbox % dETable(jAtom) = curbox % dETable(jAtom) + Ele
            endif
          endif
        endif

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then

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
  subroutine New_LJ_Q_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
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

    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
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
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: ep, sig_sq, q
    real(dp) :: LJ, Ele
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

!    globIndx = curBox % MolGlobalIndx(disp(1)%molType, disp(1)%molIndx)
!    call curBox % GetMolData(globIndx, molEnd=molEnd, molStart=molStart)
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

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
  end subroutine
  !=====================================================================
  subroutine AtomExchange_LJ_Q_Cut(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(AtomExchange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtomNew, iAtomOld, jAtom, remLen, jNei
    integer :: newType1, oldType1
    integer :: atmType2, globIndx
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: LJNew, LJOld
    real(dp) :: EleNew, EleOld
    real(dp) :: ep_old, ep_new
    real(dp) :: sig_sq_old, sig_sq_new
    real(dp) :: qij_old, qij_new
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp


    iAtomNew = disp(1) % newAtmIndx
    iAtomOld = disp(1) % oldAtmIndx
    newType1 = curbox % AtomType(iAtomNew)
    oldType1 = curbox % AtomType(iAtomOld) 
    do jNei = 1, curbox%NeighList(1)%nNeigh(iAtomOld)
      jAtom = curbox%NeighList(1)%list(jNei, iAtomOld)
      atmType2 = curbox % AtomType(jAtom)

      rx = curbox % atoms(1, iAtomOld) - curbox % atoms(1, jAtom)
      ry = curbox % atoms(2, iAtomOld) - curbox % atoms(2, jAtom)
      rz = curbox % atoms(3, iAtomOld) - curbox % atoms(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < self%rCutSq) then
          if(rsq < self%rLJCutSq) then
              ep_new = self % epsTable(atmType2, newType1)
              sig_sq_new = self % sigTable(atmType2, newType1)
              LJNew = (sig_sq_new/rsq)
              LJNew = LJNew * LJNew * LJNew
              LJNew = ep_new * LJNew * (LJNew-1E0_dp)

              ep_old = self % epsTable(atmType2, oldType1)
              sig_sq_old = self % sigTable(atmType2, oldType1)
              LJOld = (sig_sq_old/rsq)
              LJOld = LJOld * LJOld * LJOld
              LJOld = ep_old * LJOld * (LJOld-1E0_dp)
              E_Diff = E_Diff + LJNew - LJOld
              curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + LJNew - LJOld
              curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJNew - LJOld
          endif
          if(rsq < self%rQCutSq) then
              qij_new = self%qTable(atmType2, newtype1)
              qij_old = self%qTable(atmType2, oldtype1)
              EleNew = 0E0_dp
              EleOld = 0E0_dp
              if( (qij_new /= 0E0_dp) .or. (qij_old /= 0E0_dp) ) then
                  r = sqrt(rsq)
                  EleNew = qij_new/r
                  EleOld = qij_Old/r
                  E_Diff = E_Diff + EleNew - EleOld
                  curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + EleNew - EleOld
                  curbox % dETable(jAtom) = curbox % dETable(jAtom) + EleNew - EleOld
              endif
          endif
        endif
    enddo
  end subroutine
!=============================================================================+
  subroutine ManyBody_LJ_Q_Cut(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept) 
    implicit none
    class(Pair_LJ_Q_Cut), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    integer, intent(in) :: atmtype1
    integer, intent(in) :: atmtypes(:)
    real(dp), intent(in) :: pos1(:)
    real(dp), intent(in) :: posN(:,:)
    logical, intent(out) :: accept
    real(dp), intent(out) :: E_Many


    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType2
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: sig_sq, ep, q
    real(dp) :: E_Pair, Ele, LJ
    real(dp) :: rmin_ij      

    if( size(posN,2) /= size(atmtypes) ) then
      write(0,*) "ERROR! Size of posN does not match the atmtypes array!"
      error stop
    endif

    if(any(atmtypes < 1))then
      write(0,*) "Bad Array Error! Zero elements detected in atom type array!"
      write(0,*) atmtypes(1:size(atmtypes))
      write(0,*) size(atmtypes)
      error stop
    endif


    accept = .true.
    E_Many = 0E0_dp
    do jAtom = 1, size(posN, 2)
!      write(*,*) jAtom
      atmType2 = atmtypes(jAtom)
      rx = pos1(1) - posN(1, jAtom)
      ry = pos1(2) - posN(2, jAtom)
      rz = pos1(3) - posN(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rmin_ij = self % rMinTable(atmType2, atmType1)          
      rsq = rx*rx + ry*ry + rz*rz
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
            E_Many = E_Many + LJ
         endif
       endif

       if(rsq < self%rQCutSq) then
         q = self%qTable(atmType2, atmType1)
         if(q /= 0E0_dp) then
           r = sqrt(rsq)
           Ele = q/r
           E_Many = E_Many + Ele
         endif
       endif
    enddo

  end subroutine
  !=====================================================================
  subroutine ProcessIO_LJ_Q_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
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
          ep = ep * inEngUnit
          sig = sig * inLenUnit
          rMin = rMin * inLenUnit


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
