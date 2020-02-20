!================================================================================
! Pedone Pair Class for the modeling of Silicates
! For details please see the original paper
! J. Phys. Chem. B, 2006, 110 (24), pp 11780–11795  DOI: 10.1021/jp0611018
! This module also includes the implicit solvent model shown in
! Computational Materials Science  2018 Volume 149 202-207, DOI:10.1016/j.commatsci.2018.03.034.
!
!================================================================================
module FF_Pair_Pedone_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_Pedone_Cut
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)

 
    logical :: implicitSolvent = .false.
    real(dp), allocatable :: repul_tab(:,:), rEqTable(:,:)
    real(dp), allocatable :: qTable(:,:), alpha_Tab(:,:), D_Tab(:,:)
    real(dp), allocatable :: bornRad(:)
    real(dp), allocatable :: q(:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor => Constructor_Pedone_Cut
      procedure, pass :: SolventFunction => Pedone_SolventFunction
      procedure, pass :: DetailedECalc => Detailed_Pedone_Cut
      procedure, pass :: DiffECalc => DiffECalc_Pedone_Cut
      procedure, pass :: ShiftECalc_Single => Shift_Pedone_Cut_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_Pedone_Cut_Multi
      procedure, pass :: NewECalc => New_Pedone_Cut
      procedure, pass :: OldECalc => Old_Pedone_Cut
      procedure, pass :: AtomExchange => AtomExchange_Pedone_Cut
      procedure, pass :: ProcessIO => ProcessIO_Pedone_Cut
      procedure, pass :: Prologue => Prologue_Pedone_Cut
      procedure, pass :: GetCutOff => GetCutOff_Pedone_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_Pedone_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    allocate(self%repul_tab(1:nAtomTypes,1:natomtypes), stat = AllocateStat)
    allocate(self%rEqTable(1:nAtomTypes,1:natomtypes), stat = AllocateStat)
    allocate(self%qTable(1:natomtypes,1:natomtypes), stat = AllocateStat)
    allocate(self%alpha_tab(1:nAtomTypes, 1:natomtypes), stat = AllocateStat)
    allocate(self%D_tab(1:nAtomTypes, 1:natomtypes), stat = AllocateStat)

    allocate(self%q(1:nAtomTypes), stat = AllocateStat)
    self%repul_tab = 0E0_dp
    self%rEqTable = 0E0_dp
    self%qTable = 0E0_dp
    self%alpha_tab = 0E0_dp
    self%D_tab = 0E0_dp

    self%q(1) = 1.2E0_dp

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation Error in the Pedone Pair Style"

  end subroutine
!======================================================================================      
  pure function Pedone_SolventFunction(self,r, q_ij, born1, born2) result(f)
    class(Pair_Pedone_Cut), intent(in) :: self
    real(dp), intent(in) :: r, q_ij, born1, born2
    real(dp) :: f

    f = sqrt(r*r + born1*born2*exp(-r*r/(4E0_dp*born1*born2) ) ) 
    f = -0.5E0_dp*(1E0_dp-1E0_dp/dieletric)*q_ij / f
  end function
!============================================================================
  subroutine DiffECalc_Pedone_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
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
        write(*,*) "Unknown Perturbation Type Encountered by the Pedone_Cut Pair Style."
    end select


  end subroutine
  !===================================================================================
  subroutine Detailed_Pedone_Cut(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: born1, born2
    real(dp) :: q_ij, alpha, delta, repul_C, r_eq
    real(dp) :: LJ, Ele, Morse, Solvent
    real(dp) :: E_LJ, E_Ele, E_Morse, E_Solvent
    real(dp) :: rmin_ij      

    E_LJ = 0E0_dp
    E_Ele = 0E0_dp
    E_Morse = 0E0_dp
    E_Solvent = 0E0_dp
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
            accept = .false.
          endif 
          r_eq = self%rEqTable(atmType2, atmType1)
          q_ij = self%qTable(atmType2, atmType1)
          alpha = self%alpha_Tab(atmType2, atmType1)
          delta = self%D_Tab(atmType2, atmType1)
          repul_C = self%repul_tab(atmType2, atmType1)

!          write(*,*) r_eq, q_ij, alpha, delta, repul_C
          LJ = 0E0_dp

          LJ = (1E0_dp/rsq)**6
          LJ = repul_C * LJ
          E_LJ = E_LJ + LJ
 
          r = sqrt(rsq)
          Ele = q_ij/r
          Solvent = 0E0_dp
          if(self%implicitSolvent) then
             born1 = self%bornRad(atmType1)
             born2 = self%bornRad(atmType2)
             Solvent = self%solventFunction(r, q_ij, born1, born2)
             E_Solvent = E_Solvent + Solvent
          endif
          E_Ele = E_Ele + Ele

          Morse = 1E0_dp - exp(-alpha*(r-r_eq))
          Morse = delta*(Morse*Morse - 1E0_dp)
          E_Morse = E_Morse + Morse

          curbox%ETable(iAtom) = curbox%ETable(iAtom) + LJ + Ele + Morse + Solvent
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + LJ + Ele + Morse + Solvent
        endif
      enddo
    enddo
  
    write(nout,*) "Lennard-Jones Energy:", E_LJ
    write(nout,*) "Electrostatic Energy:", E_Ele
    write(nout,*) "Moorse Energy:", E_Morse
    if(self%implicitSolvent) then
      write(nout,*) "Solvation Energy:", E_Solvent
    endif
      
    E_T = E_LJ + E_Solvent + E_Ele + E_Morse
  end subroutine
  !=====================================================================
  subroutine Shift_Pedone_Cut_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: born1, born2
    real(dp) :: q_ij, alpha, delta, repul_C, r_eq
    real(dp) :: LJ, Ele, Morse, Solvent
    real(dp) :: rmin_ij      

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)
        repul_C = self % repul_tab(atmType2, atmType1)
        q_ij = self % qTable(atmType2, atmType1)
        delta = self % D_Tab(atmType2, atmType1)
        alpha = self%alpha_Tab(atmType2, atmType1)
        r_eq = self%rEqTable(atmType2, atmType1)
        rx = disp(iDisp)%x_new  -  atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) sqrt(rsq)
        if(rsq < self%rCutSq) then
          rmin_ij = self % rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif 

          LJ = (1E0_dp/rsq)**6
          LJ = repul_C * LJ
          E_Diff = E_Diff + LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ

          r = sqrt(rsq)
          Ele = q_ij/r
          E_Diff = E_Diff + Ele
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + Ele
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + Ele


          if(self%implicitSolvent) then
             born1 = self%bornRad(atmType1)
             born2 = self%bornRad(atmType2)
             Solvent = self%solventFunction(r, q_ij, born1, born2)
             E_Diff = E_Diff + Solvent
             curbox % dETable(iAtom) = curbox % dETable(iAtom) + Solvent
             curbox % dETable(jAtom) = curbox % dETable(jAtom) + Solvent
          endif

          Morse = 1E0_dp - exp(-alpha*(r-r_eq))
          Morse = delta*(Morse*Morse - 1E0_dp)
          E_Diff = E_Diff + Morse
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + Morse
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + Morse
        endif

        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then

          LJ = (1E0_dp/rsq)**6
          LJ = repul_C * LJ
          E_Diff = E_Diff - LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ

          r = sqrt(rsq)
          Ele = q_ij/r
          E_Diff = E_Diff - Ele
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - Ele
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - Ele


          if(self%implicitSolvent) then
             born1 = self%bornRad(atmType1)
             born2 = self%bornRad(atmType2)
             Solvent = self%solventFunction(r, q_ij, born1, born2)
             E_Diff = E_Diff - Solvent
             curbox % dETable(iAtom) = curbox % dETable(iAtom) - Solvent
             curbox % dETable(jAtom) = curbox % dETable(jAtom) - Solvent
          endif

          Morse = 1E0_dp - exp(-alpha*(r-r_eq))
          Morse = delta*(Morse*Morse - 1E0_dp)
          E_Diff = E_Diff - Morse
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - Morse
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - Morse
        endif
      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine New_Pedone_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: born1, born2
    real(dp) :: q_ij, alpha, delta, repul_C, r_eq
    real(dp) :: LJ, Ele, Morse, Solvent
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.

!    write(*,*) "Length:", size(tempNNei)
!    write(*,*) "Length:", size(tempList)
!    write(*,*)
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
          repul_C = self%repul_tab(atmType2, atmType1)
          q_ij = self%qTable(atmType2, atmType1)
          delta = self%D_Tab(atmType2, atmType1)
          alpha = self % alpha_Tab(atmType2, atmType1)
          r_eq = self % rEqTable(atmType2, atmType1)


          LJ = (1E0_dp/rsq)**6
          LJ = repul_C * LJ
          E_Diff = E_Diff + LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ

          r = sqrt(rsq)
          Ele = q_ij/r
          E_Diff = E_Diff + Ele
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + Ele
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + Ele


          if(self%implicitSolvent) then
             born1 = self%bornRad(atmType1)
             born2 = self%bornRad(atmType2)
             Solvent = self%solventFunction(r, q_ij, born1, born2)
             E_Diff = E_Diff + Solvent
             curbox % dETable(iAtom) = curbox % dETable(iAtom) + Solvent
             curbox % dETable(jAtom) = curbox % dETable(jAtom) + Solvent
          endif

          Morse = 1E0_dp - exp(-alpha*(r-r_eq))
          Morse = delta*(Morse*Morse - 1E0_dp)
          E_Diff = E_Diff + Morse
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + Morse
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + Morse
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_Pedone_Cut(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2, globIndx
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: born1, born2
    real(dp) :: q_ij, alpha, delta, repul_C, r_eq
    real(dp) :: LJ, Ele, Morse, Solvent
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

!    write(*,*) disp(1)%molType, disp(1)%molIndx
!    globIndx = curBox % MolGlobalIndx(disp(1)%molType, )
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
      atmType1 = curbox % AtomType(iAtom) 
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)

        rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          repul_C = self%repul_tab(atmType2, atmType1)
          LJ = (1E0_dp/rsq)**6
          LJ = repul_C * LJ
          E_Diff = E_Diff - LJ
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ

          q_ij = self%qTable(atmType2, atmType1)
          r = sqrt(rsq)
          Ele = q_ij/r
          E_Diff = E_Diff - Ele
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - Ele
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - Ele


          if(self%implicitSolvent) then
             born1 = self%bornRad(atmType1)
             born2 = self%bornRad(atmType2)
             Solvent = self%solventFunction(r, q_ij, born1, born2)
             E_Diff = E_Diff - Solvent
             curbox % dETable(iAtom) = curbox % dETable(iAtom) - Solvent
             curbox % dETable(jAtom) = curbox % dETable(jAtom) - Solvent
          endif

          delta = self%D_Tab(atmType2, atmType1)
          alpha = self % alpha_Tab(atmType2, atmType1)
          r_eq = self % rEqTable(atmType2, atmType1)
          Morse = 1E0_dp - exp(-alpha*(r-r_eq))
          Morse = delta*(Morse*Morse - 1E0_dp)
          E_Diff = E_Diff - Morse
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - Morse
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - Morse
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine AtomExchange_Pedone_Cut(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(AtomExchange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtomNew, iAtomOld, jAtom, remLen, jNei
    integer :: newType1, oldType1
    integer :: atmType2, globIndx
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: born1, born2
    real(dp) ::  alpha, delta, repul_C, r_eq
    real(dp) :: LJ
    real(dp) :: repul_C_new, repul_C_old
    real(dp) :: q_ij_new, LJNew, EleNew, MorseNew, SolventNew
    real(dp) :: q_ij_old, LJOld, EleOld, MorseOld, SolventOld
    real(dp) :: rmin_ij      

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)

    E_Diff = 0E0_dp


    iAtomNew = disp(1) % newAtmIndx
    iAtomOld = disp(1) % oldAtmIndx
    newType1 = curbox % AtomType(iAtomNew)
    oldType1 = curbox % AtomType(iAtomOld) 
    do jNei = 1, curbox%NeighList(1)%nNeigh(iAtomOld)
        jAtom = curbox%NeighList(1)%list(jNei, iAtomOld)
        atmType2 = curbox % AtomType(jAtom)

        rx = atoms(1, iAtomOld) - atoms(1, jAtom)
        ry = atoms(2, iAtomOld) - atoms(2, jAtom)
        rz = atoms(3, iAtomOld) - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          repul_C_new = self%repul_tab(atmType2, newType1)
          repul_C_old = self%repul_tab(atmType2, oldType1)
          LJ = (1E0_dp/rsq)**6
          LJNew = repul_C_new * LJ
          LJOld = repul_C_old * LJ
          E_Diff = E_Diff + LJNew - LJOld
          curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + LJNew - LJOld
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJNew - LJOld

          r = sqrt(rsq)
          q_ij_new = self%qTable(atmType2, newType1)
          q_ij_old = self%qTable(atmType2, oldType1)
          EleNew = q_ij_new/r
          EleOld = q_ij_old/r
          E_Diff = E_Diff + EleNew - EleOld
          curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + EleNew - EleOld
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + EleNew - EleOld


          if(self%implicitSolvent) then

             born1 = self%bornRad(newType1)
             born2 = self%bornRad(atmType2)
             SolventNew = self%solventFunction(r, q_ij_new, born1, born2)

             born1 = self%bornRad(oldType1)
             SolventOld = self%solventFunction(r, q_ij_old, born1, born2)

             E_Diff = E_Diff + SolventNew - SolventOld
             curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + SolventNew - SolventOld
             curbox % dETable(jAtom) = curbox % dETable(jAtom) + SolventNew - SolventOld
          endif

          delta = self%D_Tab(atmType2, oldType1)
          alpha = self % alpha_Tab(atmType2, oldType1)
          r_eq = self % rEqTable(atmType2, oldType1)
          MorseOld = 1E0_dp - exp(-alpha*(r-r_eq))
          MorseOld = delta*(MorseOld*MorseOld - 1E0_dp)

          delta = self%D_Tab(atmType2, newType1)
          alpha = self % alpha_Tab(atmType2, newType1)
          r_eq = self % rEqTable(atmType2, newType1)
          MorseNew = 1E0_dp - exp(-alpha*(r-r_eq))
          MorseNew = delta*(MorseNew*MorseNew - 1E0_dp)
          E_Diff = E_Diff + MorseNew - MorseOld
          curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + MorseNew - MorseOld
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + MorseNew - MorseOld

        endif
      enddo
  end subroutine

  !=====================================================================
  subroutine ProcessIO_Pedone_Cut(self, line)
    use Common_MolInfo, only: nAtomTypes
    use ClassyConstants, only: coulombConst
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    logical :: logicVal
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: q, repC, delta, alpha, rCut, rMin, rEq
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("rcut")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) rCut
        self % rCut = rCut
        self % rCutSq = rCut * rCut

      case("implicitsolvent")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) self%implicitSolvent

      case("born")
        if(.not. allocated(self%bornRad)) then
          allocate(self%bornRad(1:nAtomTypes))
        endif
        call CountCommands(line, nPar)
        
        do jType = 1, nAtomTypes
          call GetXCommand(line, command, jType+1, lineStat)
          read(command, *) self%bornRad(jType)
        enddo


      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      call CountCommands(line, nPar)
      select case(nPar)
        case(7)
          read(line, *) type1, repC, rEq, alpha, delta, q, rMin
          self%q(type1) = q
          repC = repC * inEngUnit/inLenUnit**12
          rEq = rEq * inLenUnit
          alpha = alpha / inLenUnit**2
          delta = delta * inEngUnit
          rMin = rMin * inLenUnit

          self%repul_tab(1, type1) = repC
          self%repul_tab(type1, 1) = repC

          self%rEqTable(1, type1) = rEq
          self%rEqTable(type1, 1) = rEq

          self%alpha_tab(1, type1) = alpha
          self%alpha_tab(type1, 1) = alpha

          self%D_tab(1, type1) = delta
          self%D_tab(type1, 1) = delta


          do jType = 1, nAtomTypes
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
  function GetCutOff_Pedone_Cut(self) result(rCut)
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
  !=====================================================================
  subroutine Prologue_Pedone_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_Pedone_Cut), intent(inout) :: self
    integer :: i, j

!    do i = 1, nAtomTypes
!      write(nout, *) (self%rMinTable(i,j), j=1,nAtomTypes)
!    enddo

  end subroutine
  !=====================================================================
end module
!=====================================================================
