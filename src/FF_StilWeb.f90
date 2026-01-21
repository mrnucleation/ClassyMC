!================================================================================
! Stillinger-Weber Potential for Monte Carlo Simulations
!================================================================================
module FF_Pair_StillWeb
  use CoordinateTypes
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use VarPrecision

  implicit none
  private

  type :: SW2Body
    real(dp) :: epsilon = 2.1683_dp
    real(dp) :: sigma = 2.0951_dp
    real(dp) :: parA = 7.049556277_dp
    real(dp) :: parB = 0.6022245584_dp
    real(dp) :: parp = 4.0_dp
    real(dp) :: parq = 0.0_dp
    real(dp) :: paraCut = 1.80_dp
    real(dp) :: rCut = 3.77_dp
    real(dp) :: rCutSq = 14.2_dp
  end type

  type :: SW3Body
    real(dp) :: epsilon = 2.1683_dp
    real(dp) :: sigma = 2.0951_dp
    real(dp) :: lambda = 21.0_dp
    real(dp) :: parGamma = 1.20_dp
    real(dp) :: paraCut = 1.80_dp
  end type

  type, public, extends(forcefield) :: Pair_StillWeb
    type(SW2Body), allocatable :: sw2(:,:)
    type(SW3Body), allocatable :: sw3(:,:,:)
    real(dp), allocatable :: rMinTable(:,:)
    real(dp) :: eVtokb = 11604.522_dp
    real(dp) :: energyConversion = 11604.522_dp
    
  contains
    procedure, pass :: Constructor => Constructor_StillWeb
    procedure, pass :: DetailedECalc => Detailed_StillWeb
    procedure, pass :: DiffECalc => DiffECalc_StillWeb
    procedure, pass :: ShiftECalc_Single => Shift_StillWeb_Single
    procedure, pass :: NewECalc => New_StillWeb
    procedure, pass :: OldECalc => Old_StillWeb
    procedure, pass :: ProcessIO => ProcessIO_StillWeb
    procedure, pass :: GetCutOff => GetCutOff_StillWeb
    procedure, pass :: Prologue => Prologue_StillWeb
    procedure, pass :: Phi2 => Phi2_StillWeb
    procedure, pass :: Phi3 => Phi3_StillWeb
  end type

contains
!==================================================================================
pure function Phi2_StillWeb(self, r, type1, type2) result(E2)
  implicit none
  class(Pair_StillWeb), intent(in) :: self
  real(dp), intent(in) :: r
  integer, intent(in) :: type1, type2
  real(dp) :: E2
  real(dp) :: epsv
  real(dp) :: sigv
  real(dp) :: Av
  real(dp) :: Bv
  real(dp) :: pv
  real(dp) :: qv
  real(dp) :: aval
  real(dp) :: sr
  real(dp) :: rcutv
  real(dp) :: argv
  real(dp) :: expv
  
  epsv = self%sw2(type1, type2)%epsilon
  sigv = self%sw2(type1, type2)%sigma
  Av = self%sw2(type1, type2)%parA
  Bv = self%sw2(type1, type2)%parB
  pv = self%sw2(type1, type2)%parp
  qv = self%sw2(type1, type2)%parq
  aval = self%sw2(type1, type2)%paraCut
  rcutv = aval * sigv
  
  if (r >= rcutv) then
    E2 = 0.0_dp
    return
  endif
  
  sr = sigv / r
  argv = sigv / (r - rcutv)
  if (argv > 50.0_dp) argv = 50.0_dp
  expv = exp(argv)
  
  if (qv == 0.0_dp) then
    E2 = epsv * Av * (Bv * sr**pv - 1.0_dp) * expv
  else
    E2 = epsv * Av * (Bv * sr**pv - sr**qv) * expv
  endif
end function
!==================================================================================
pure function Phi3_StillWeb(self, rij, rik, cosTheta, type1, type2, type3) result(E3)
  implicit none
  class(Pair_StillWeb), intent(in) :: self
  real(dp), intent(in) :: rij, rik, cosTheta
  integer, intent(in) :: type1, type2, type3
  real(dp) :: E3
  real(dp) :: epsv
  real(dp) :: sigv
  real(dp) :: lamv
  real(dp) :: gamv
  real(dp) :: aval
  real(dp) :: rcutv
  real(dp) :: argij
  real(dp) :: argik
  real(dp) :: expv
  real(dp) :: angterm
  
  epsv = self%sw3(type1, type2, type3)%epsilon
  sigv = self%sw3(type1, type2, type3)%sigma
  lamv = self%sw3(type1, type2, type3)%lambda
  gamv = self%sw3(type1, type2, type3)%parGamma
  aval = self%sw3(type1, type2, type3)%paraCut
  rcutv = aval * sigv
  
  if (rij >= rcutv .or. rik >= rcutv) then
    E3 = 0.0_dp
    return
  endif
  
  argij = gamv * sigv / (rij - rcutv)
  argik = gamv * sigv / (rik - rcutv)
  if (argij > 50.0_dp) argij = 50.0_dp
  if (argik > 50.0_dp) argik = 50.0_dp
  expv = exp(argij + argik)
  angterm = (cosTheta + 1.0_dp/3.0_dp)**2
  E3 = epsv * lamv * expv * angterm
end function
!=============================================================================
subroutine Constructor_StillWeb(self)
  use Common_MolInfo, only: nAtomTypes
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  integer :: i, j

  allocate(self%rMinTable(nAtomTypes, nAtomTypes))
  allocate(self%sw2(nAtomTypes, nAtomTypes))
  allocate(self%sw3(nAtomTypes, nAtomTypes, nAtomTypes))

  self%rMinTable = 0.25_dp
  
  do i = 1, nAtomTypes
    do j = 1, nAtomTypes
      self%sw2(i,j)%rCut = self%sw2(i,j)%paraCut * self%sw2(i,j)%sigma
      self%sw2(i,j)%rCutSq = self%sw2(i,j)%rCut**2
    enddo
  enddo

  self%rCut = 1.80_dp * 2.0951_dp
  self%rCutSq = self%rCut**2
end subroutine
!===================================================================================
subroutine Detailed_StillWeb(self, curbox, E_T, accept)
  use ParallelVar, only: nout
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  real(dp), intent(inout) :: E_T
  logical, intent(out) :: accept
  
  integer :: iAtom, jAtom, kAtom, atmType1, atmType2, atmType3
  real(dp) :: rxij, ryij, rzij, rijSq, rij
  real(dp) :: rxik, ryik, rzik, rikSq, rik
  real(dp) :: cosTheta, E2total, E3total
  real(dp), pointer :: atoms(:,:) => null()

  call curbox%GetCoordinates(atoms)
  accept = .true.
  E2total = 0.0_dp
  E3total = 0.0_dp
  curbox%ETable = 0.0_dp

  do iAtom = 1, curbox%nMaxAtoms - 1
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType1 = curbox%AtomType(iAtom)
    do jAtom = iAtom + 1, curbox%nMaxAtoms
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      atmType2 = curbox%AtomType(jAtom)
      
      rxij = atoms(1, jAtom) - atoms(1, iAtom)
      ryij = atoms(2, jAtom) - atoms(2, iAtom)
      rzij = atoms(3, jAtom) - atoms(3, iAtom)
      call curbox%Boundary(rxij, ryij, rzij)
      rijSq = rxij**2 + ryij**2 + rzij**2
      
      if (rijSq < self%sw2(atmType1, atmType2)%rCutSq) then
        if (rijSq < self%rMinTable(atmType1, atmType2)) then
          accept = .false.
          return
        endif
        rij = sqrt(rijSq)
        E2total = E2total + self%Phi2(rij, atmType1, atmType2)
      endif
    enddo
  enddo

  do iAtom = 1, curbox%nMaxAtoms
    if (.not. curbox%IsActive(iAtom)) cycle
    atmType1 = curbox%AtomType(iAtom)
    do jAtom = 1, curbox%nMaxAtoms
      if (jAtom == iAtom) cycle
      if (.not. curbox%IsActive(jAtom)) cycle
      if (curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)) cycle
      atmType2 = curbox%AtomType(jAtom)
      
      rxij = atoms(1, jAtom) - atoms(1, iAtom)
      ryij = atoms(2, jAtom) - atoms(2, iAtom)
      rzij = atoms(3, jAtom) - atoms(3, iAtom)
      call curbox%Boundary(rxij, ryij, rzij)
      rijSq = rxij**2 + ryij**2 + rzij**2
      if (rijSq >= self%sw2(atmType1, atmType2)%rCutSq) cycle
      rij = sqrt(rijSq)
      
      do kAtom = jAtom + 1, curbox%nMaxAtoms
        if (kAtom == iAtom) cycle
        if (.not. curbox%IsActive(kAtom)) cycle
        if (curbox%MolIndx(kAtom) == curbox%MolIndx(iAtom)) cycle
        atmType3 = curbox%AtomType(kAtom)
        
        rxik = atoms(1, kAtom) - atoms(1, iAtom)
        ryik = atoms(2, kAtom) - atoms(2, iAtom)
        rzik = atoms(3, kAtom) - atoms(3, iAtom)
        call curbox%Boundary(rxik, ryik, rzik)
        rikSq = rxik**2 + ryik**2 + rzik**2
        if (rikSq >= self%sw2(atmType1, atmType3)%rCutSq) cycle
        rik = sqrt(rikSq)
        
        cosTheta = (rxij*rxik + ryij*ryik + rzij*rzik) / (rij * rik)
        E3total = E3total + self%Phi3(rij, rik, cosTheta, atmType1, atmType2, atmType3)
      enddo
    enddo
  enddo

  E2total = E2total * self%energyConversion
  E3total = E3total * self%energyConversion

  write(nout,*) "Stillinger-Weber Energy:"
  write(nout,*) "  Two-body:", E2total
  write(nout,*) "  Three-body:", E3total
  write(nout,*) "  Total:", E2total + E3total

  E_T = E2total + E3total
end subroutine
!============================================================================
subroutine DiffECalc_StillWeb(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  class(simBox), intent(inout) :: curbox
  class(Perturbation), intent(inout), target :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  accept = .true.
  curbox%dETable = 0.0_dp
  E_Diff = 0.0_dp

  select type(disp)
    class is(Displacement)
      call self%ShiftECalc_Single(curbox, disp, E_Diff, accept)
    class is(Addition)
      call self%NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)
    class is(Deletion)
      call self%OldECalc(curbox, disp, E_Diff)
    class default
      error stop "Unknown perturbation type in SW"
  end select
end subroutine
!=====================================================================
! Differential energy calculation for single atom displacement
! For SW, moving atom i affects:
!   1. Two-body: all pairs (i, j)
!   2. Three-body with i as central atom: (i, j, k)
!   3. Three-body with neighbor m as central: (m, i, k) and (m, j, i)
!=====================================================================
subroutine Shift_StillWeb_Single(self, curbox, disp, E_Diff, accept)
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Displacement), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept
  
  integer :: iDisp, iAtom, jNei, jAtom, kNei, kAtom, mNei, mAtom
  integer :: atmType1, atmType2, atmType3, atmTypeM
  real(dp) :: rxij, ryij, rzij, rijSq, rij
  real(dp) :: rxik, ryik, rzik, rikSq, rik
  real(dp) :: rxijnew, ryijnew, rzijnew, rijSqnew, rijnew
  real(dp) :: rxiknew, ryiknew, rziknew, rikSqnew, riknew
  real(dp) :: rxmj, rymj, rzmj, rmjSq, rmj
  real(dp) :: rxmk, rymk, rzmk, rmkSq, rmk
  real(dp) :: rxmi, rymi, rzmi, rmiSq, rmi
  real(dp) :: rxminew, ryminew, rzminew, rmiSqnew, rminew
  real(dp) :: cosT, cosTnew
  real(dp) :: E2diff, E3diff, rCutSq, rCutSqik
  real(dp) :: E3old, E3new
  real(dp), pointer :: atoms(:,:) => null()
  integer, pointer :: nNeigh(:) => null()
  integer, pointer :: neighlist(:,:) => null()
  
  ! Track affected neighbors for recalculation
  integer :: nRecalc
  integer :: recalcList(500)

  call curbox%GetCoordinates(atoms)
  call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)

  E_Diff = 0.0_dp
  E2diff = 0.0_dp
  E3diff = 0.0_dp
  accept = .true.
  nRecalc = 0

  do iDisp = 1, size(disp)
    iAtom = disp(iDisp)%atmIndx
    atmType1 = curbox%AtomType(iAtom)

    ! ============================================================
    ! PART 1: Pairs and triplets with moving atom i as central
    ! ============================================================
    do jNei = 1, nNeigh(iAtom)
      jAtom = neighlist(jNei, iAtom)
      atmType2 = curbox%AtomType(jAtom)
      rCutSq = self%sw2(atmType1, atmType2)%rCutSq

      ! Old i-j distance
      rxij = atoms(1, jAtom) - atoms(1, iAtom)
      ryij = atoms(2, jAtom) - atoms(2, iAtom)
      rzij = atoms(3, jAtom) - atoms(3, iAtom)
      call curbox%Boundary(rxij, ryij, rzij)
      rijSq = rxij**2 + ryij**2 + rzij**2

      ! New i-j distance
      rxijnew = atoms(1, jAtom) - disp(iDisp)%x_new
      ryijnew = atoms(2, jAtom) - disp(iDisp)%y_new
      rzijnew = atoms(3, jAtom) - disp(iDisp)%z_new
      call curbox%Boundary(rxijnew, ryijnew, rzijnew)
      rijSqnew = rxijnew**2 + ryijnew**2 + rzijnew**2

      ! Check overlap
      if (rijSqnew < self%rMinTable(atmType1, atmType2)) then
        accept = .false.
        return
      endif

      ! Track affected neighbor for part 2
      if (rijSq < rCutSq .or. rijSqnew < rCutSq) then
        if (all(recalcList(1:nRecalc) /= jAtom)) then
          nRecalc = nRecalc + 1
          recalcList(nRecalc) = jAtom
        endif
      endif

      ! Two-body energy difference
      if (rijSq < rCutSq) then
        rij = sqrt(rijSq)
        E2diff = E2diff - self%Phi2(rij, atmType1, atmType2)
      endif
      if (rijSqnew < rCutSq) then
        rijnew = sqrt(rijSqnew)
        E2diff = E2diff + self%Phi2(rijnew, atmType1, atmType2)
      endif

      ! Three-body with i as central atom: triplets (i, j, k)
      do kNei = jNei + 1, nNeigh(iAtom)
        kAtom = neighlist(kNei, iAtom)
        atmType3 = curbox%AtomType(kAtom)
        rCutSqik = self%sw2(atmType1, atmType3)%rCutSq

        ! Old i-k distance
        rxik = atoms(1, kAtom) - atoms(1, iAtom)
        ryik = atoms(2, kAtom) - atoms(2, iAtom)
        rzik = atoms(3, kAtom) - atoms(3, iAtom)
        call curbox%Boundary(rxik, ryik, rzik)
        rikSq = rxik**2 + ryik**2 + rzik**2

        ! New i-k distance
        rxiknew = atoms(1, kAtom) - disp(iDisp)%x_new
        ryiknew = atoms(2, kAtom) - disp(iDisp)%y_new
        rziknew = atoms(3, kAtom) - disp(iDisp)%z_new
        call curbox%Boundary(rxiknew, ryiknew, rziknew)
        rikSqnew = rxiknew**2 + ryiknew**2 + rziknew**2

        ! Old triplet (i, j, k)
        if (rijSq < rCutSq .and. rikSq < rCutSqik) then
          rij = sqrt(rijSq)
          rik = sqrt(rikSq)
          cosT = (rxij*rxik + ryij*ryik + rzij*rzik) / (rij * rik)
          E3diff = E3diff - self%Phi3(rij, rik, cosT, atmType1, atmType2, atmType3)
        endif

        ! New triplet (i, j, k)
        if (rijSqnew < rCutSq .and. rikSqnew < rCutSqik) then
          rijnew = sqrt(rijSqnew)
          riknew = sqrt(rikSqnew)
          cosTnew = (rxijnew*rxiknew + ryijnew*ryiknew + rzijnew*rziknew) / (rijnew * riknew)
          E3diff = E3diff + self%Phi3(rijnew, riknew, cosTnew, atmType1, atmType2, atmType3)
        endif
      enddo
    enddo
  enddo

  ! ============================================================
  ! PART 2: Triplets where moving atom is NOT central
  ! For each neighbor m of moving atom i, we need to recalculate
  ! triplets (m, j, k) where the moving atom is j or k
  ! ============================================================
  do mNei = 1, nRecalc
    mAtom = recalcList(mNei)
    atmTypeM = curbox%AtomType(mAtom)
    
    ! Loop over all pairs of neighbors of m
    do jNei = 1, nNeigh(mAtom)
      jAtom = neighlist(jNei, mAtom)
      if (jAtom == mAtom) cycle
      
      ! Skip if j is not the moving atom (we only care about pairs involving moving atom)
      if (jAtom /= disp(1)%atmIndx) cycle
      
      atmType2 = curbox%AtomType(jAtom)
      rCutSq = self%sw2(atmTypeM, atmType2)%rCutSq
      
      ! Old m-i distance (j = moving atom = i)
      rxmi = atoms(1, jAtom) - atoms(1, mAtom)
      rymi = atoms(2, jAtom) - atoms(2, mAtom)
      rzmi = atoms(3, jAtom) - atoms(3, mAtom)
      call curbox%Boundary(rxmi, rymi, rzmi)
      rmiSq = rxmi**2 + rymi**2 + rzmi**2
      
      ! New m-i distance
      rxminew = disp(1)%x_new - atoms(1, mAtom)
      ryminew = disp(1)%y_new - atoms(2, mAtom)
      rzminew = disp(1)%z_new - atoms(3, mAtom)
      call curbox%Boundary(rxminew, ryminew, rzminew)
      rmiSqnew = rxminew**2 + ryminew**2 + rzminew**2
      
      ! Loop over third atoms k (neighbors of m, not moving atom, not m)
      do kNei = 1, nNeigh(mAtom)
        kAtom = neighlist(kNei, mAtom)
        if (kAtom == mAtom .or. kAtom == jAtom) cycle
        
        atmType3 = curbox%AtomType(kAtom)
        rCutSqik = self%sw2(atmTypeM, atmType3)%rCutSq
        
        ! m-k distance (unchanged)
        rxmk = atoms(1, kAtom) - atoms(1, mAtom)
        rymk = atoms(2, kAtom) - atoms(2, mAtom)
        rzmk = atoms(3, kAtom) - atoms(3, mAtom)
        call curbox%Boundary(rxmk, rymk, rzmk)
        rmkSq = rxmk**2 + rymk**2 + rzmk**2
        
        if (rmkSq >= rCutSqik) cycle
        rmk = sqrt(rmkSq)
        
        ! Old triplet (m, i_old, k)
        E3old = 0.0_dp
        if (rmiSq < rCutSq) then
          rmi = sqrt(rmiSq)
          cosT = (rxmi*rxmk + rymi*rymk + rzmi*rzmk) / (rmi * rmk)
          E3old = self%Phi3(rmi, rmk, cosT, atmTypeM, atmType2, atmType3)
        endif
        
        ! New triplet (m, i_new, k)
        E3new = 0.0_dp
        if (rmiSqnew < rCutSq) then
          rminew = sqrt(rmiSqnew)
          cosTnew = (rxminew*rxmk + ryminew*rymk + rzminew*rzmk) / (rminew * rmk)
          E3new = self%Phi3(rminew, rmk, cosTnew, atmTypeM, atmType2, atmType3)
        endif
        
        E3diff = E3diff + (E3new - E3old)
      enddo
    enddo
  enddo

  E_Diff = (E2diff + E3diff) * self%energyConversion
end subroutine
!======================================================================================
subroutine New_StillWeb(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Addition), intent(in) :: disp(:)
  integer, intent(in) :: tempList(:,:), tempNNei(:)
  real(dp), intent(inout) :: E_Diff
  logical, intent(out) :: accept

  E_Diff = 0.0_dp
  accept = .true.
end subroutine
!======================================================================================
subroutine Old_StillWeb(self, curbox, disp, E_Diff)
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  class(SimBox), intent(inout) :: curbox
  type(Deletion), intent(in) :: disp(:)
  real(dp), intent(inout) :: E_Diff

  E_Diff = 0.0_dp
end subroutine
!=============================================================================
subroutine ProcessIO_StillWeb(self, line)
  use Input_Format, only: CountCommands, GetXCommand, maxLineLen
  use Units, only: inLenUnit
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  character(len=maxLineLen), intent(in) :: line
  character(len=30) :: command
  integer :: lineStat, type1
  real(dp) :: realVal

  call GetXCommand(line, command, 1, lineStat)

  select case(trim(adjustl(command)))
    case("epsilon")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) realVal
      self%sw2(type1, type1)%epsilon = realVal
      self%sw3(type1, type1, type1)%epsilon = realVal

    case("sigma")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) realVal
      self%sw2(type1, type1)%sigma = realVal * inLenUnit
      self%sw3(type1, type1, type1)%sigma = realVal * inLenUnit
      self%sw2(type1, type1)%rCut = self%sw2(type1,type1)%paraCut * self%sw2(type1,type1)%sigma
      self%sw2(type1, type1)%rCutSq = self%sw2(type1, type1)%rCut**2

    case("lambda")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) type1
      call GetXCommand(line, command, 3, lineStat)
      read(command, *) realVal
      self%sw3(type1, type1, type1)%lambda = realVal

    case("units")
      call GetXCommand(line, command, 2, lineStat)
      select case(trim(command))
        case("ev", "eV")
          self%energyConversion = self%eVtokb
        case("kb")
          self%energyConversion = 1.0_dp
      end select

    case("rcut")
      call GetXCommand(line, command, 2, lineStat)
      read(command, *) realVal
      self%rCut = realVal * inLenUnit
      self%rCutSq = self%rCut**2

    case default
      lineStat = -1
  end select
end subroutine
!=============================================================================
subroutine Prologue_StillWeb(self)
  use ParallelVar, only: nout
  implicit none
  class(Pair_StillWeb), intent(inout) :: self

  write(nout,*) "Stillinger-Weber Potential"
  write(nout,*) "  Cutoff:", self%rCut
  write(nout,*) "  Energy conversion:", self%energyConversion
end subroutine
!=============================================================================
function GetCutOff_StillWeb(self) result(rCut)
  implicit none
  class(Pair_StillWeb), intent(inout) :: self
  real(dp) :: rCut
  rCut = self%rCut
end function
!=============================================================================
end module FF_Pair_StillWeb
