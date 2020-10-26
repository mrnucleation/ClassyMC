!================================================================================
module FF_Pair_LJ_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_LJ_Cut
    real(dp), allocatable :: eps(:)
    real(dp), allocatable :: sig(:)
    real(dp), allocatable :: rMin(:)

    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor => Constructor_LJ_Cut
      procedure, pass :: DetailedECalc => Detailed_LJ_Cut
      procedure, pass :: DiffECalc => DiffECalc_LJ_Cut
      procedure, pass :: ShiftECalc_Single => Shift_LJ_Cut_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_LJ_Cut_Multi
      procedure, pass :: NewECalc => New_LJ_Cut
      procedure, pass :: OldECalc => Old_LJ_Cut
      procedure, pass :: OrthoVolECalc => OrthoVol_LJ_Cut
      procedure, pass :: AtomExchange => AtomExchange_LJ_Cut
      procedure, pass :: ProcessIO => ProcessIO_LJ_Cut
      procedure, pass :: Prologue => Prologue_LJ_Cut
      procedure, pass :: GetCutOff => GetCutOff_LJ_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_LJ_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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
          LJ = ep * LJ * (LJ-1E0_dp)
          E_LJ = E_LJ + LJ
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + LJ
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + LJ 
        endif
      enddo
    enddo
  
!    write(nout,*) "Lennard-Jones Energy:", E_LJ
      
    E_T = E_LJ    
  end subroutine
!============================================================================
  subroutine DiffECalc_LJ_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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
!         write(*,*) size(tempList)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

      class is(OrthoVolChange)
         call self % OrthoVolECalc( curbox, disp, E_Diff, accept)

      class is(AtomExchange)
         call self % AtomExchange( curbox, disp, E_Diff, accept)

      class default
        write(*,*) "Unknown Perturbation Type Encountered by the LJ_Cut Pair Style."
        stop

    end select


  end subroutine

  !=====================================================================
  subroutine Shift_LJ_Cut_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Displacement), intent(in) :: disp(:)
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

!      write(*,*) iAtom, curbox%NeighList(1)%nNeigh(iAtom)
!      write(*,*) iAtom, curbox%NeighList(1)%list(:, iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) sqrt(rsq)

        atmType2 = curbox % AtomType(jAtom)
        ep = self % epsTable(atmType2, atmType1)
        sig_sq = self % sigTable(atmType2, atmType1)  
        if(rsq < self%rCutSq) then
          rmin_ij = self % rMinTable(atmType2, atmType1)          
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
!        write(*,*) sqrt(rsq)
        if(rsq < self%rCutSq) then
          if(iAtom == jAtom) then
            write(*,*) "NeighborList Error!", iAtom, jAtom
          endif
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
  subroutine New_LJ_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Addition), intent(in) :: disp(:)
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
          ep = self%epsTable(atmType2, atmType1)
          sig_sq = self%sigTable(atmType2, atmType1)          

          LJ = (sig_sq/rsq)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0_dp)
          E_Diff = E_Diff + LJ
!          write(*,*) iAtom, jAtom, rsq, LJ
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
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ, LJ2
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

!    globIndx = curBox % MolGlobalIndx(disp(1)%molType, disp(1)%molIndx)
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
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
!          write(*,*) iAtom, jAtom, rsq, LJ
          E_Diff = E_Diff - LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine OrthoVol_LJ_Cut(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(OrthoVolChange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iAtom, jAtom, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    integer :: molIndx1, molIndx2
    real(dp) :: dxi, dyi, dzi
    real(dp) :: dxj, dyj, dzj
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: rx2, ry2, rz2, rsq2
    real(dp) :: ep, sig_sq, ratio
    real(dp) :: LJ, LJ2
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp
    do iAtom = 1, curbox%nMaxAtoms-1
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
        cycle
      endif
      atmType1 = curbox % AtomType(iAtom)
      molIndx1 = curbox % MolIndx(iAtom)
      dxi = curbox % centerMass(1, molIndx1) * (disp(1)%xScale-1E0_dp)
      dyi = curbox % centerMass(2, molIndx1) * (disp(1)%yScale-1E0_dp)
      dzi = curbox % centerMass(3, molIndx1) * (disp(1)%zScale-1E0_dp)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        if(jAtom <= iAtom) then
          cycle
        endif
        molIndx2 = curbox % MolIndx(jAtom)
        dxj = curbox % centerMass(1,molIndx2) * (disp(1)%xScale-1E0_dp)
        dyj = curbox % centerMass(2,molIndx2) * (disp(1)%yScale-1E0_dp)
        dzj = curbox % centerMass(3,molIndx2) * (disp(1)%zScale-1E0_dp)

        rx =  curbox % atoms(1, iAtom) + dxi - curbox % atoms(1, jAtom) - dxj
        ry =  curbox % atoms(2, iAtom) + dyi - curbox % atoms(2, jAtom) - dyj
        rz =  curbox % atoms(3, iAtom) + dzi - curbox % atoms(3, jAtom) - dzj

        call curbox%BoundaryNew(rx, ry, rz, disp)
        rsq = rx*rx + ry*ry + rz*rz

        atmType2 = curbox % AtomType(jAtom)
        if(rsq < self%rCutSq) then
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

    E_Diff = E_Diff - curbox%ETotal
    curbox % dETable = curbox%dETable - curbox % ETable

  end subroutine
  !=====================================================================
  subroutine AtomExchange_LJ_Cut(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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
    real(dp) :: ep_old, ep_new
    real(dp) :: sig_sq_old, sig_sq_new
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

!    write(*,*) disp(1)%molType, disp(1)%molIndx
!    globIndx = curBox % MolGlobalIndx(disp(1)%molType, )
!    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

    iAtomNew = disp(1) % newAtmIndx
    iAtomOld = disp(1) % oldAtmIndx
    newType1 = curbox % AtomType(iAtomNew)
    oldType1 = curbox % AtomType(iAtomOld) 
!    write(*,*) iAtomNew, iAtomOld, newType1, oldType1
    do jNei = 1, curbox%NeighList(1)%nNeigh(iAtomOld)
        jAtom = curbox%NeighList(1)%list(jNei, iAtomOld)
        atmType2 = curbox % AtomType(jAtom)

!        rmin_ij = self % rMinTable(atmType1, atmType2)          

        rx = curbox % atoms(1, iAtomOld) - curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtomOld) - curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtomOld) - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          rmin_ij = self%rMinTable(atmType2, newType1)
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif

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
    enddo
  end subroutine
  !=====================================================================
  subroutine ProcessIO_LJ_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
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
      end select
    endif

  end subroutine
  !=============================================================================+
  function GetCutOff_LJ_Cut(self) result(rCut)
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
  !=====================================================================
  subroutine Prologue_LJ_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    integer :: i, j

!    do i = 1, nAtomTypes
!      write(nout, *) (self%epsTable(i,j), j=1,nAtomTypes)
!    enddo

!    do i = 1, nAtomTypes
!      write(nout, *) (self%sigTable(i,j), j=1,nAtomTypes)
!    enddo

!    do i = 1, nAtomTypes
!      write(nout, *) (self%rMinTable(i,j), j=1,nAtomTypes)
!    enddo


  end subroutine
  !=====================================================================
end module
!=====================================================================
