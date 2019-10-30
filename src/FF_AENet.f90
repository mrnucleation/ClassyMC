!=============================================================================+
! AENet forcefield interface function. The AENet files are not currently
! included with the Classy code base and must be compiled separately.
! To add AENet to the compilation apply the -DAENET flag in the Makefile
!=============================================================================+
module FF_AENet
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use Input_Format, only: maxLineLen
  use VarPrecision
  use CoordinateTypes

  !AENet Library functions
#ifdef AENET
!  use predict_lib, only: initialize_lib, get_energy_lib
  use aenet, only: aenet_atomic_energy, aenet_init, aenet_load_potential, &
                   aenet_Rc_min, aenet_Rc_max
  use geometry, only: geo_recip_lattice
  use input, only: InputData
#endif

  real(dp), parameter :: boltz = 8.6173303E-5_dp

  type, public, extends(forcefield) :: Pair_AENet
!    real(dp) :: rCut, rCutSq
    logical :: initialized = .false.
    integer :: neilistindx = 1
    integer, allocatable :: atomTypes(:)
    character(len=maxLineLen), allocatable :: inputfiles(:)
    real(dp), allocatable :: tempcoords(:,:)
    real(dp)  :: box(3,3)
    real(dp)  :: boxrecp(3,3)

    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    contains
    !If -DAENET is not passed to the compiler this defaults back to the
    !template's functions which effiectively do nothing.
#ifdef AENET
      procedure, pass :: Constructor => Constructor_AENet
      procedure, pass :: DetailedECalc => DetailedECalc_AENet
      procedure, pass :: Predict => DetailedECalc_AENet_Predict
      procedure, pass :: DiffECalc => DiffECalc_AENet
      procedure, pass :: VolECalc => VolECalc_AENet_V2
      procedure, pass :: ProcessIO => ProcessIO_AENet
#endif
      procedure, pass :: Prologue => Prologue_AENet
!      procedure, pass :: GetCutOff
  end type

  contains
#ifdef AENET
!=============================================================================+
  subroutine Constructor_AENet(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes, AtomData
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(Pair_AENet), intent(inout) :: self

    type(InputData) :: inp
    character(len=100) :: str1, str2
    integer :: iBox, atomLimit, iType
    integer :: AllocateStat
    character(len=5), allocatable :: symbols(:)

    write(nout, *) "Initializing AENet"
    if(.not. self%initialized) then
      allocate(symbols(1:nAtomTypes))
      do iType = 1, nAtomTypes
        symbols(iType) = adjustl(trim(AtomData(iType)%Symb))
      enddo
      call aenet_init(symbols, AllocateStat)
      deallocate(symbols)
    endif
!    call initialize_lib(str1, str2, inp)



    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    self%rMin = 0E0_dp
    self%rMinTable = 0E0_dp

  end subroutine
!==========================================================================+
  subroutine DetailedECalc_AENet(self, curbox, E_T, accept)
!    use boxData, only: self%boxArray
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(Pair_AENet), intent(inout) :: self
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
    self%atomTypes = 0

    do iAtom = 1, curbox%nMaxAtoms
      if( .not. curbox%IsActive(iAtom) ) then
        cycle
      endif
      nCurAtoms = 0
      atmType1 = curbox % AtomType(iAtom)
      do jAtom = 1, curbox%nMaxAtoms
        if( .not. curbox%IsActive(jAtom) ) cycle
        if(iAtom == jAtom) cycle
!        if( curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)  ) cycle
        rx = atoms(1, jAtom) - atoms(1, iAtom) 
        ry = atoms(2, jAtom) - atoms(2, iAtom) 
        rz = atoms(3, jAtom) - atoms(3, iAtom) 
        call curbox % Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self % rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            write(*,*) sqrt(rsq)
            write(*,*) iAtom, jAtom
            write(*,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
            write(*,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
            write(*,*) "ERROR! Overlaping atoms found in the current configuration!"
          endif 
          nCurAtoms = nCurAtoms + 1
          self%atomTypes(nCurAtoms) = atmType2

          !Reshift and store atomic coordinates such that the 
          !nearest image of the atom is used.
          self%tempcoords(1, nCurAtoms) = rx + atoms(1, iAtom) 
          self%tempcoords(2, nCurAtoms) = ry + atoms(2, iAtom) 
          self%tempcoords(3, nCurAtoms) = rz + atoms(3, iAtom)
        endif
      enddo
!      write(*,*) iAtom,atmType1, nCurAtoms, self%rCutSq
      call aenet_atomic_energy(atoms(1:3, iAtom), atmType1, nCurAtoms, &
                               self%tempcoords(1:3, 1:nCurAtoms), &
                               self%atomtypes(1:nCurAtoms), E_Atom, stat) 
      E_T = E_T + E_Atom
      curbox%ETable(iAtom) = E_Atom/boltz
    enddo
    write(nout, *) "Raw Neuro Net Energy:", E_T
    E_T = E_T / boltz
    write(nout, *) "Total Neuro Net Energy:", E_T

!    call self%Predict(curbox, E_T, accept)
!    stop

  end subroutine
!==========================================================================+
! Old Version of the DetailedCalc.  This is becoming obsolete and will be
! replaced.
  subroutine DetailedECalc_AENet_Predict(self, curbox, E_T, accept)
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(Pair_AENet), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    logical :: pbc = .false.
    integer :: i,j
    integer :: iAtom
    integer :: nTotalMol
    integer :: nCurAtoms = 0
    real(dp) :: Ecoh
    real(dp) :: xoffset, yoffset, zoffset
    real(dp) :: xmax, ymax, zmax
    real(dp) :: dx, dy, dz
    real(dp) :: xscale, yscale, zscale
    real(dp) :: tempdim(1:3,1:3)


    accept = .true.
    xoffset = 0E0_dp
    yoffset = 0E0_dp
    zoffset = 0E0_dp
    nCurAtoms = 0
    self%box = 0E0_dp
    self%tempcoords = 0E0_dp
    self%atomTypes = 1

    xmax = 0E0_dp
    ymax = 0E0_dp
    zmax = 0E0_dp
    do iAtom = 1, curbox%nMaxAtoms
      if( .not. curbox%IsActive(iAtom) ) then
        cycle
      endif
      nCurAtoms = nCurAtoms + 1
      self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
      self%tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom) 
      self%tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom) 
      self%tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom) 

      select type(curbox)
        type is(SimpleBox)
          if(xoffset > self%tempcoords(1, nCurAtoms)) then
            xoffset = self%tempcoords(1, nCurAtoms)
          endif
          if(yoffset > self%tempcoords(2, nCurAtoms)) then
            yoffset = self%tempcoords(2, nCurAtoms)
          endif
          if(zoffset > self%tempcoords(3, nCurAtoms)) then
            zoffset = self%tempcoords(3, nCurAtoms)
          endif
       end select

     enddo
     !Collect the self%box dimensions and determine the offset. 
     select type(curbox)
       class is(CubeBox)
         pbc = .true.
         call curbox%GetDimensions(tempdim)
         xoffset = tempdim(1, 1)
         yoffset = tempdim(1, 2)
         zoffset = tempdim(1, 3)
         self%box(1,1) = tempdim(2, 1) - tempdim(1, 1)
         self%box(2,2) = tempdim(2, 2) - tempdim(1, 2)
         self%box(3,3) = tempdim(2, 3) - tempdim(1, 3)

       class is(OrthoBox)
         pbc = .true.
         call curbox%GetDimensions(tempdim)
         xoffset = tempdim(1, 1)
         yoffset = tempdim(1, 2)
         zoffset = tempdim(1, 3)
         self%box(1,1) = tempdim(2, 1) - tempdim(1, 1)
         self%box(2,2) = tempdim(2, 2) - tempdim(1, 2)
         self%box(3,3) = tempdim(2, 3) - tempdim(1, 3)

     end select

     do iAtom = 1, nCurAtoms
       self%tempcoords(1, iAtom) = self%tempcoords(1, iAtom) - xoffset
       self%tempcoords(2, iAtom) = self%tempcoords(2, iAtom) - yoffset
       self%tempcoords(3, iAtom) = self%tempcoords(3, iAtom) - zoffset
       select type(curbox)
         class is(SimpleBox)
           if(xmax < self%tempcoords(1, iAtom)) then
             xmax = self%tempcoords(1, iAtom)
           endif
           if(ymax < self%tempcoords(2, iAtom)) then
             ymax = self%tempcoords(2, iAtom)
           endif
           if(zmax < self%tempcoords(3, iAtom)) then
             zmax = self%tempcoords(3, iAtom)
           endif
       end select
     enddo

     select type(curbox)
       type is(SimpleBox)
         self%box(1,1) = 3*(xmax + 1e-5)
         self%box(2,2) = 3*(ymax + 1e-5)
         self%box(3,3) = 3*(zmax + 1e-5)
     end select


     self%boxrecp(:,:) = geo_recip_lattice(self%box)
     self%tempcoords(1:3,1:nCurAtoms) = matmul(self%boxrecp(1:3,1:3), self%tempcoords(1:3, 1:nCurAtoms)) / (2E0_dp * pi)

!     do iAtom = 1, ncurAtoms
!          call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
!                                   E_i, stat)
!     enddo
!     call get_energy_lib(self%box(1:3, 1:3), nCurAtoms, &
!                         self%tempcoords(1:3, 1:nCurAtoms), &
!                         self%atomTypes(1:nCurAtoms), pbc, &
!                         Ecoh, E_T)

     write(nout, *) "Raw Neuro Net Energy:", E_T, Ecoh
!     E_T = E_T * curbox%nMolTotal / boltz
     E_T = E_T / boltz
     write(nout, *) "Total Neuro Net Energy:", E_T



  end subroutine
!============================================================================
  subroutine DiffECalc_AENet(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    use ClassyConstants, only: pi
    implicit none
    class(Pair_AENet), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    logical :: pbc=.false.
    integer :: nCurAtoms = 0
    integer :: iAtom, jAtom, jNei, j, stat
    integer :: atmType1, atmType2, molIndx1, molIndx2, nTotalMol
    integer :: iDisp, iRecalc, molStart, molEnd

    integer :: nRecalc
    integer :: recalcList(1:300)
    real(dp) :: E_T, ecoh
    real(dp) :: E_New, E_Old, E_Atom
    real(dp) :: rmin_ij
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: xoffset, yoffset, zoffset
    real(dp) :: xi, yi, zi
    real(dp) :: xj, yj, zj
    real(dp) :: xmax, ymax, zmax
    real(dp) :: tempdim(1:3,1:3)
    real(dp) :: tempatom(1:3)

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
    !passing the configuration to AENet. In addition create a list of atoms
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


     !Compute the old position's energy for the neighboring atoms contained within the recalc list.
     E_Old = 0E0_dp
     do iRecalc = 1, nRecalc
        nCurAtoms = 0
        iAtom = recalcList(iRecalc)
        atmType1 = curbox % AtomType(iAtom)
        do jNei = 1, nNeigh(iAtom)
          jAtom = neighlist(jNei, iAtom)
          rx = atoms(1, jAtom) - atoms(1, iAtom) 
          ry = atoms(2, jAtom) - atoms(2, iAtom) 
          rz = atoms(3, jAtom) - atoms(3, iAtom) 
          call curbox % Boundary(rx, ry, rz)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < self%rCutSq) then
            atmType2 = curbox % AtomType(jAtom)
            nCurAtoms = nCurAtoms + 1
            self%atomTypes(nCurAtoms) = atmType2
            self%tempcoords(1, nCurAtoms) = rx + atoms(1, iAtom) 
            self%tempcoords(2, nCurAtoms) = ry + atoms(2, iAtom) 
            self%tempcoords(3, nCurAtoms) = rz + atoms(3, iAtom)
          endif
        enddo
        call aenet_atomic_energy(atoms(1:3, iAtom), atmType1, nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomtypes(1:nCurAtoms),&
                                E_Atom, stat) 
        E_Old = E_Old + E_Atom
        curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_Atom/boltz

     enddo


     !Compute the new position's energy for the neighboring atoms contained within the recalc list.
     E_New = 0E0_dp
     do iRecalc = 1, nRecalc
        nCurAtoms = 0
        iAtom = recalcList(iReCalc)
        atmType1 = curbox % AtomType(iAtom)
        do jNei = 1, nNeigh(iAtom)
          jAtom = neighlist(jNei, iAtom)
          xj = atoms(1, jAtom)
          yj = atoms(2, jAtom)
          zj = atoms(3, jAtom)
          select type(disp)
            class is(Displacement)
              do iDisp = 1, size(disp)
                if(jAtom == disp(iDisp)%atmIndx) then
                  xj = disp(iDisp)%x_new
                  yj = disp(iDisp)%y_new
                  zj = disp(iDisp)%z_new
                  exit
                endif
              enddo

            class is(Deletion)
              molIndx2 = curbox%molIndx(jAtom)
              if(disp(1)%molIndx == molIndx2) then
                cycle
              endif

          end select
          rx = xj - atoms(1, iAtom)
          ry = yj - atoms(2, iAtom)
          rz = zj - atoms(3, iAtom)
          call curbox % Boundary(rx, ry, rz)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < self%rCutSq) then
            atmType2 = curbox % AtomType(jAtom)
            nCurAtoms = nCurAtoms + 1
            self%atomTypes(nCurAtoms) = atmType2
            !Storing atom positions as the nearest image to atom_i.
            self%tempcoords(1, nCurAtoms) = rx + atoms(1, iAtom)
            self%tempcoords(2, nCurAtoms) = ry + atoms(2, iAtom)
            self%tempcoords(3, nCurAtoms) = rz + atoms(3, iAtom)
          endif
        enddo
        select type(disp)
          class is(Addition)
            do iDisp = 1, size(disp)
              xj = disp(iDisp)%x_new
              yj = disp(iDisp)%y_new
              zj = disp(iDisp)%z_new
              rx = xj - atoms(1, iAtom)
              ry = yj - atoms(2, iAtom)
              rz = zj - atoms(3, iAtom)
              call curbox % Boundary(rx, ry, rz)
              rsq = rx*rx + ry*ry + rz*rz
              if(rsq < self%rCutSq) then
                jAtom = disp(iDisp)%atmindx
                atmType2 = curbox % AtomType(jAtom)
                nCurAtoms = nCurAtoms + 1
                self%atomTypes(nCurAtoms) = atmType2
                self%tempcoords(1, nCurAtoms) = rx + atoms(1, iAtom)
                self%tempcoords(2, nCurAtoms) = ry + atoms(2, iAtom)
                self%tempcoords(3, nCurAtoms) = rz + atoms(3, iAtom)
              endif
            enddo

        end select
        call aenet_atomic_energy(atoms(1:3, iAtom), atmType1, nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomtypes(1:nCurAtoms),&
                                E_Atom, stat) 
        E_New = E_New + E_Atom
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_Atom/boltz
     enddo

     !Now calculate the contribution of the atoms that were moved during this move.
     select type(disp)
      !-----------------------------------------------------
      class is(Displacement)
        do iDisp = 1, size(disp)
          nCurAtoms = 0
          iAtom = disp(iDisp)%atmIndx
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, nNeigh(iAtom)
            jAtom = neighlist(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)

            rx = atoms(1, jAtom) - disp(iDisp)%x_new
            ry = atoms(2, jAtom) - disp(iDisp)%y_new
            rz = atoms(3, jAtom) - disp(iDisp)%z_new
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq) then
              atmType2 = curbox % AtomType(jAtom)
              nCurAtoms = nCurAtoms + 1
              self%atomTypes(nCurAtoms) = atmType2
              self%tempcoords(1, nCurAtoms) = rx + disp(iDisp)%x_new
              self%tempcoords(2, nCurAtoms) = ry + disp(iDisp)%y_new
              self%tempcoords(3, nCurAtoms) = rz + disp(iDisp)%z_new
            endif

          enddo
          tempatom(1) = disp(iDisp)%x_new
          tempatom(2) = disp(iDisp)%y_new
          tempatom(3) = disp(iDisp)%z_new
          call aenet_atomic_energy(tempatom(1:3), atmType1, nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomtypes(1:nCurAtoms),&
                                E_Atom, stat) 
          E_New = E_New + E_Atom
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_Atom/boltz


          nCurAtoms = 0
          do jNei = 1, nNeigh(iAtom)
            jAtom = neighlist(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)
            rx = atoms(1, jAtom)  -  atoms(1, iAtom)
            ry = atoms(2, jAtom)  -  atoms(2, iAtom)
            rz = atoms(3, jAtom)  -  atoms(3, iAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq) then
              atmType2 = curbox % AtomType(jAtom)
              nCurAtoms = nCurAtoms + 1
              self%atomTypes(nCurAtoms) = atmType2
              self%tempcoords(1, nCurAtoms) = rx + atoms(1, iAtom)
              self%tempcoords(2, nCurAtoms) = ry + atoms(2, iAtom)
              self%tempcoords(3, nCurAtoms) = rz + atoms(3, iAtom)
            endif
          enddo
          call aenet_atomic_energy(atoms(1:3, iAtom), atmType1, nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomtypes(1:nCurAtoms),&
                                E_Atom, stat) 
          E_Old = E_Old + E_Atom
          curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_Atom/boltz
        enddo
      !-----------------------------------------------------
      class is(Addition)
        do iDisp = 1, size(disp)
          nCurAtoms = 0
          iAtom = disp(iDisp)%atmIndx
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, tempNNei(iDisp)
            jAtom = templist(jNei, iDisp)
            atmType2 = curbox % AtomType(jAtom)

            rx = atoms(1, jAtom) - disp(iDisp)%x_new
            ry = atoms(2, jAtom) - disp(iDisp)%y_new
            rz = atoms(3, jAtom) - disp(iDisp)%z_new
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq) then
              atmType2 = curbox % AtomType(jAtom)
              nCurAtoms = nCurAtoms + 1
              self%atomTypes(nCurAtoms) = atmType2
              self%tempcoords(1, nCurAtoms) = rx + disp(iDisp)%x_new
              self%tempcoords(2, nCurAtoms) = ry + disp(iDisp)%y_new
              self%tempcoords(3, nCurAtoms) = rz + disp(iDisp)%z_new
            endif

          enddo
          tempatom(1) = disp(iDisp)%x_new
          tempatom(2) = disp(iDisp)%y_new
          tempatom(3) = disp(iDisp)%z_new
          call aenet_atomic_energy(tempatom(1:3), atmType1, nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomtypes(1:nCurAtoms),&
                                E_Atom, stat) 
          E_New = E_New + E_Atom
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + E_Atom/boltz
        enddo
      !-----------------------------------------------------
      class is(Deletion)
        call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)
        do iAtom = molStart, molEnd
          atmType1 = curbox % AtomType(iAtom)
          nCurAtoms = 0
          do jNei = 1, nNeigh(iAtom)
            jAtom = neighlist(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)
            rx = atoms(1, jAtom)  -  atoms(1, iAtom)
            ry = atoms(2, jAtom)  -  atoms(2, iAtom)
            rz = atoms(3, jAtom)  -  atoms(3, iAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < self%rCutSq) then
              atmType2 = curbox % AtomType(jAtom)
              nCurAtoms = nCurAtoms + 1
              self%atomTypes(nCurAtoms) = atmType2
              self%tempcoords(1, nCurAtoms) = rx + atoms(1, iAtom)
              self%tempcoords(2, nCurAtoms) = ry + atoms(2, iAtom)
              self%tempcoords(3, nCurAtoms) = rz + atoms(3, iAtom)
            endif
          enddo
          call aenet_atomic_energy(atoms(1:3, iAtom), atmType1, nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomtypes(1:nCurAtoms),&
                                E_Atom, stat) 
          E_Old = E_Old + E_Atom
          curbox%dETable(iAtom) = curbox%dETable(iAtom) - E_Atom/boltz

        enddo

    end select
    E_Diff = (E_New - E_Old)/boltz


  end subroutine
!============================================================================
  subroutine VolECalc_AENet_V2(self, curbox, disp, E_Diff, accept)
    use ClassyConstants, only: pi
    implicit none
    class(Pair_AENet), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    class(OrthoVolChange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    logical :: pbc=.false.
    integer :: nCurAtoms = 0
    integer :: iAtom, jAtom, jNei, j
    integer :: atmType1, atmType2, molIndx1, molIndx2
    integer :: iDisp, nTotalMol, stat
    real(dp) :: E_T, E_Atom
    real(dp) :: rmin_ij
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: newcoords(1:3)
    real(dp) :: dx, dy, dz
    real(dp) :: dxj, dyj, dzj
    real(dp) :: xscale, yscale, zscale
    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    accept = .true.
    E_Diff = 0E0_dp
    curbox%dETable = 0E0_dp


    select type(curbox)
      class is(SimpleBox)
        call curbox%GetCoordinates(atoms)
        call curbox%GetNeighborList(self%neilistindx, neighlist, nNeigh)
    end select


    nTotalMol = curbox%nMolTotal
    select type(disp)
      class is(OrthoVolChange)
        xscale = disp(1)%xScale
        yscale = disp(1)%yScale
        zscale = disp(1)%zScale
    end select

    do iAtom = 1, curbox%nMaxAtoms
      if( .not. curbox%IsActive(iAtom) ) then
        cycle
      endif
      nCurAtoms = 0
      atmType1 = curbox % AtomType(iAtom)
      molIndx1 = curbox % MolIndx(iAtom)
      dx = curbox % centerMass(1, molIndx1) * (xScale-1E0_dp)
      dy = curbox % centerMass(2, molIndx1) * (yScale-1E0_dp)
      dz = curbox % centerMass(3, molIndx1) * (zScale-1E0_dp)

      newcoords(1) = atoms(1, iAtom) + dx
      newcoords(2) = atoms(2, iAtom) + dy
      newcoords(3) = atoms(3, iAtom) + dz
      do jNei = 1, nNeigh(iAtom)
        jAtom = NeighList(jNei, iAtom)
        molIndx2 = curbox % MolIndx(jAtom)
        atmType2 = curbox % AtomType(jAtom)
        dxj = curbox % centerMass(1, molIndx2) * (xScale-1E0_dp)
        dyj = curbox % centerMass(2, molIndx2) * (yScale-1E0_dp)
        dzj = curbox % centerMass(3, molIndx2) * (zScale-1E0_dp)
        rx = atoms(1, jAtom) + dxj - newcoords(1)
        ry = atoms(2, jAtom) + dyj - newcoords(2)
        rz = atoms(3, jAtom) + dzj - newcoords(3) 
!        rx = -rx
!        ry = -ry
!        rz = -rz
        call curbox%BoundaryNew(rx, ry, rz, disp)
        rsq = rx*rx + ry*ry + rz*rz
        rmin_ij = self % rMinTable(atmType2, atmType1)      
        if(rsq < rmin_ij) then
          accept = .false.
          return
        endif 
        if(rsq < self%rCutSq) then
          nCurAtoms = nCurAtoms + 1
          self%atomTypes(nCurAtoms) = atmType2

          !Reshift and store atomic coordinates such that the 
          !nearest image of the atom is used.
          self%tempcoords(1, nCurAtoms) = rx + newcoords(1)
          self%tempcoords(2, nCurAtoms) = ry + newcoords(2)
          self%tempcoords(3, nCurAtoms) = rz + newcoords(3)
        endif
      enddo

      call aenet_atomic_energy(newcoords(1:3), atmType1, nCurAtoms, &
                               self%tempcoords(1:3, 1:nCurAtoms), &
                               self%atomtypes(1:nCurAtoms), E_Atom, stat) 
      E_T = E_T + E_Atom
      curbox%dETable(iAtom) = E_Atom/boltz
    enddo
!    write(*,*) E_T, curbox%ETotal
    E_T = E_T / boltz
!    write(*,*) E_T, curbox%ETotal
    E_Diff = E_T - curbox%ETotal
!    E_Diff = E_Diff - curbox%ETotal
    curbox % dETable = curbox%dETable - curbox % ETable
  end subroutine
!============================================================================
  subroutine VolECalc_AENet(self, curbox, disp, E_Diff, accept)
    use ClassyConstants, only: pi
    implicit none
    class(Pair_AENet), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(OrthoVolChange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    logical :: pbc=.false.
    integer :: nCurAtoms = 0
    integer :: iAtom, jAtom, jNei, j
    integer :: atmType1, atmType2, molIndx1, molIndx2
    integer :: iDisp, nTotalMol
    real(dp) :: E_T, ecoh
    real(dp) :: rmin_ij
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: xoffset, yoffset, zoffset
    real(dp) :: dx, dy, dz
    real(dp) :: dxj, dyj, dzj
    real(dp) :: xscale, yscale, zscale
    real(dp) :: xmax, ymax, zmax
    real(dp) :: tempdim(1:3,1:3)
    real(dp), pointer :: atoms(:,:) => null()

    accept = .true.
    E_Diff = 0E0_dp

    call curbox%GetCoordinates(atoms)

    nTotalMol = curbox%nMolTotal
    select type(disp)
      class is(OrthoVolChange)
        xscale = disp(1)%xScale
        yscale = disp(1)%yScale
        zscale = disp(1)%zScale


    end select

    !Check the rMin criteria first to ensure there is no overlap prior to
    !passing the configuration to AENet. This only needs to be checked
    !if the box is being compressed.  If it is expanded and it previously
    !statisfied rmin constraints, it won't overlap on expansion.
    if(any([xscale, yscale, zscale] < 1E0_dp)) then
        select type(disp)
          !-----------------------------------------------------
          class is(OrthoVolChange)
            do iAtom = 1, curBox%nMaxAtoms
              if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
                cycle
              endif
              atmType1 = curbox % AtomType(iAtom)
              molIndx1 = curbox % MolIndx(iAtom)
              dx = curbox % centerMass(1, molIndx1) * (xScale-1E0_dp)
              dy = curbox % centerMass(2, molIndx1) * (yScale-1E0_dp)
              dz = curbox % centerMass(3, molIndx1) * (zScale-1E0_dp)
              do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
                jAtom = curbox%NeighList(1)%list(jNei, iAtom)
                molIndx2 = curbox % MolIndx(jAtom)
                atmType2 = curbox % AtomType(jAtom)
                dxj = curbox % centerMass(1, molIndx2) * (xScale-1E0_dp)
                dyj = curbox % centerMass(2, molIndx2) * (yScale-1E0_dp)
                dzj = curbox % centerMass(3, molIndx2) * (zScale-1E0_dp)
                rx = atoms(1, iAtom) + dx  -  atoms(1, jAtom) - dxj
                ry = atoms(2, iAtom) + dy  -  atoms(2, jAtom) - dyj
                rz = atoms(3, iAtom) + dz  -  atoms(3, jAtom) - dzj
                rsq = rx*rx + ry*ry + rz*rz
                call curbox%BoundaryNew(rx, ry, rz, disp)
                rmin_ij = self % rMinTable(atmType2, atmType1)      
                if(rsq < rmin_ij) then
                  accept = .false.
                  return
                endif 
              enddo
            enddo
        end select
    endif

    ! Convert the Classy Array into an Array that AENet can read.
    nCurAtoms = 0
    self%tempcoords = 0E0_dp
    self%atomTypes = 1
    select type(disp)
      class is(OrthoVolChange)
        do iAtom = 1, curbox%nMaxAtoms
          if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
            cycle
          endif
          molIndx1 = curbox % MolIndx(iAtom)
          dx = curbox % centerMass(1, molIndx1) * (xScale-1E0_dp)
          dy = curbox % centerMass(2, molIndx1) * (yScale-1E0_dp)
          dz = curbox % centerMass(3, molIndx1) * (zScale-1E0_dp)
          nCurAtoms = nCurAtoms + 1
          self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
          self%tempcoords(1, nCurAtoms) = atoms(1, iAtom) + dx
          self%tempcoords(2, nCurAtoms) = atoms(2, iAtom) + dy
          self%tempcoords(3, nCurAtoms) = atoms(3, iAtom) + dz
        enddo
    end select

    accept = .true.
    xoffset = 0E0_dp
    yoffset = 0E0_dp
    zoffset = 0E0_dp
    select type(curbox)
      type is(SimpleBox)
        do iAtom = 1, nCurAtoms
          if(xoffset > self%tempcoords(1, iAtom)) then
            xoffset = self%tempcoords(1, iAtom)
          endif
          if(yoffset > self%tempcoords(2, iAtom)) then
            yoffset = self%tempcoords(2, iAtom)
          endif
          if(zoffset > self%tempcoords(3, iAtom)) then
            zoffset = self%tempcoords(3, iAtom)
          endif
      enddo 
    end select
    self%box = 0E0_dp

    xmax = 0E0_dp
    ymax = 0E0_dp
    zmax = 0E0_dp

     !Collect the self%box dimensions and determine the offset. 
     select type(curbox)
       class is(CubeBox)
         pbc = .true.
         call curbox%GetDimensions(tempdim)
         ! self%box dimensions given back as xlow, xhigh, ylow, yhigh, etc.  
         ! Need to be converted to AENet format
         xoffset = tempdim(1, 1)*xScale
         yoffset = tempdim(1, 2)*yScale
         zoffset = tempdim(1, 3)*zScale
         self%box(1,1) = (tempdim(2, 1) - tempdim(1, 1))*xScale
         self%box(2,2) = (tempdim(2, 2) - tempdim(1, 2))*yScale
         self%box(3,3) = (tempdim(2, 3) - tempdim(1, 3))*zScale

       class is(OrthoBox)
         pbc = .true.
         call curbox%GetDimensions(tempdim)
         ! self%box dimensions given back as xlow, xhigh, ylow, yhigh, etc.  
         ! Need to be converted to AENet format
         xoffset = tempdim(1, 1)*xScale
         yoffset = tempdim(1, 2)*yScale
         zoffset = tempdim(1, 3)*zScale
         self%box(1,1) = (tempdim(2, 1) - tempdim(1, 1))*xScale
         self%box(2,2) = (tempdim(2, 2) - tempdim(1, 2))*yScale
         self%box(3,3) = (tempdim(2, 3) - tempdim(1, 3))*zScale

     end select

     do iAtom = 1, nCurAtoms
       self%tempcoords(1, iAtom) = self%tempcoords(1, iAtom) - xoffset
       self%tempcoords(2, iAtom) = self%tempcoords(2, iAtom) - yoffset
       self%tempcoords(3, iAtom) = self%tempcoords(3, iAtom) - zoffset
       select type(curbox)
         type is(SimpleBox)
           if(xmax < self%tempcoords(1, iAtom)) then
             xmax = self%tempcoords(1, iAtom)
           endif
           if(ymax < self%tempcoords(2, iAtom)) then
             ymax = self%tempcoords(2, iAtom)
           endif
           if(zmax < self%tempcoords(3, iAtom)) then
             zmax = self%tempcoords(3, iAtom)
           endif
       end select
     enddo

     select type(curbox)
       type is(SimpleBox)
         self%box(1,1) = xmax + 1e-2
         self%box(2,2) = ymax + 1e-2
         self%box(3,3) = zmax + 1e-2
!       class default
!         do iAtom = 1, nCurAtoms
!           self%tempcoords(1, iAtom) = self%tempcoords(1, iAtom)/self%box(1,1)
!           self%tempcoords(2, iAtom) = self%tempcoords(2, iAtom)/self%box(2,2)
!           self%tempcoords(3, iAtom) = self%tempcoords(3, iAtom)/self%box(3,3)
!         enddo
     end select


     self%boxrecp = 0E0_dp
     self%boxrecp(:,:) = geo_recip_lattice(self%box)
     self%tempcoords(1:3,1:nCurAtoms) = matmul(self%boxrecp(1:3,1:3), self%tempcoords(1:3,1:nCurAtoms))/ (2E0_dp * pi)
!     call get_energy_lib(self%box(1:3,1:3), nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomTypes(1:nCurAtoms), pbc, Ecoh, E_T)
     E_T = E_T / boltz
     E_Diff = E_T - curbox%ETotal

  end subroutine
!=============================================================================+
!End -DAENET safety block
!=============================================================================+
  subroutine ProcessIO_AENet(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use ParallelVar, only: nout
    implicit none
    class(Pair_AENet), intent(inout) :: self
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

        write(nout, *) "Loading potential from file: ", self%inputfiles(type1)
        call aenet_load_potential(type1, self%inputfiles(type1), stat)
        self%rCut = aenet_Rc_max
        self%rCutSq = aenet_Rc_max * aenet_Rc_max

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

!=============================================================================+
#endif 
!End AENet Safety Block
!=============================================================================+
  subroutine Prologue_AENet(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_AENet), intent(inout) :: self
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

     write(nout, *) "AENet Cutoff Distance:", self%rCut

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
