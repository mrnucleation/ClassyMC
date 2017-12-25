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
    self%rMinTable = 0.5E0_dp
    self%rCut = 3E0_dp
    self%rCutSq = 3E0_dp**2

!    write(*,*) 
    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_LJ_Cut(self, curbox, E_T)
    use ParallelVar, only: nout
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
    class(SimBox), intent(inout) :: curbox
      real(dp), intent(inOut) :: E_T
      integer :: iAtom,jAtom
      integer :: atmType1, atmType2
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: E_LJ
      real(dp) :: rmin_ij      

      E_LJ = 0E0
      curbox%ETable = 0E0
      do iAtom = 1, curbox%nAtoms-1
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = iAtom+1, curbox%nAtoms
          atmType2 = curbox % AtomType(jAtom)
          ep = self % epsTable(atmType1,atmType2)
          sig_sq = self % sigTable(atmType1,atmType2)          
          rmin_ij = self % rMinTable(atmType1,atmType2)          

!          write(*,*) ep, sig_sq, rmin_ij
          rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
          call curbox%Boundary(rx, ry, rz)
          rsq = rx**2 + ry**2 + rz**2
          if(rsq < self%rCutSq) then
            if(rsq < rmin_ij) then
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
  subroutine Shift_LJ_Cut_Single(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
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
!          write(*,*) rx, ry, rz
        call curbox%Boundary(rx, ry, rz)
!          write(*,*) rx, ry, rz
        rsq = rx*rx + ry*ry + rz*rz
!          write(*,*) rsq
        if(rsq < self%rCutSq) then
!            if(rsq < rmin_ij) then
!            endif 
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
!        write(*,*) rsq
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
      class(Pair_LJ_Cut), intent(in) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
   
  end subroutine
  !=====================================================================
  subroutine New_LJ_Cut(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_LJ_Cut), intent(in) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff
      integer :: iDisp, iAtom, jAtom, dispLen
      integer :: atmType1, atmType2
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      dispLen = size(disp)
      E_Diff = 0E0
      curbox%dETable = 0E0


      do iDisp = 1, dispLen
        iAtom = disp(iDisp)%atmIndx
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, curbox % nAtoms        
          atmType2 = curbox % AtomType(jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = self%rMinTable(atmType1,atmType2)          

          rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
          ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
          rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < self%rCutSq) then
            LJ = (sig_sq/rsq)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)              
            E_Diff = E_Diff - LJ
            curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
            curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
          endif
        enddo
      enddo
  end subroutine
  !=====================================================================
  subroutine Old_LJ_Cut(self, curbox, atmIndx, E_Diff)
    implicit none
      class(Pair_LJ_Cut), intent(in) :: self
      class(SimBox), intent(inout) :: curbox
      real(dp), intent(inOut) :: E_Diff
      integer, intent(in) :: atmIndx(:)
      integer :: iIndx, iAtom, jAtom, remLen
      integer :: atmType1, atmType2
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      remLen = size(atmIndx)
      E_Diff = 0E0
      curbox%dETable = 0E0


      do iIndx = 1, remLen
        iAtom = atmIndx(iIndx)
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, curbox % nAtoms        
          atmType1 = curbox % AtomType(jAtom)
          ep = self % epsTable(atmType1,atmType2)
          sig_sq = self % sigTable(atmType1,atmType2)          
          rmin_ij = self % rMinTable(atmType1,atmType2)          

          rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ
        enddo
      enddo
  end subroutine
  !=====================================================================
  subroutine ProcessIO_LJ_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetAllCommands
    implicit none
    class(Pair_LJ_Cut), intent(inout) :: self
    character(len=*), intent(in) :: line
    character(len=30), allocatable :: parlist(:)
    integer :: jType, lineStat
    integer :: type1, type2
    real(dp) :: ep, sig
  
    call GetAllCommands(line, parlist, lineStat)
    select case(size(parlist))
      case(3)
        read(line, *) type1, ep, sig
        do jType = 1, nAtomTypes
          self%epsTable(type1, jType) = ep
          self%epsTable(jType, type1) = ep

          self%sigTable(type1, jType) = sig
          self%sigTable(jType, type1) = sig
        enddo

      case(4)
        read(line, *) type1, type2, ep, sig
        self%epsTable(type1, type2) = ep
        self%epsTable(type2, type1) = ep

        self%sigTable(type1, type2) = sig
        self%sigTable(type2, type1) = sig
    end select


    deallocate(parlist)
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
