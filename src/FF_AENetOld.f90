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
  use VarPrecision
  use CoordinateTypes

  !AENet Library functions
#ifdef AENET
  use predict_lib, only: initialize_lib, get_energy_lib
!  use aenet, only: aenet_init, aenet_atomic_energy, &
!                   aenet_Rc_min, aenet_Rc_max,
  use geometry, only: geo_recip_lattice
  use input, only: InputData
#endif

  real(dp), parameter :: boltz = 8.6173303E-5_dp

  type, public, extends(forcefield) :: AENet
!    real(dp) :: rCut, rCutSq
!    logical :: initialized = .false.
    integer, allocatable :: atomTypes(:)
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
      procedure, pass :: DiffECalc => DiffECalc_AENet
#endif
      procedure, pass :: ProcessIO => ProcessIO_AENet
      procedure, pass :: Prologue => Prologue_AENet
!      procedure, pass :: GetCutOff
  end type

  contains
#ifdef AENET
!=============================================================================+
  subroutine Constructor_AENet(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(AENet), intent(inout) :: self

    type(InputData) :: inp
    character(len=100) :: str1, str2
    integer :: iBox, atomLimit, iType
    integer :: AllocateStat
     write(nout, *) "Initializing AENet"
!     inp = read_InpPredict(inFile)
     call initialize_lib(str1, str2, inp)
!     call initialize(str1, str2, indata)
!     call aenet_init(inp%typeName, stat)
!     do iType = 1, inp%nTypes
!       call aenet_load_potential(itype, inp%netFile(itype), stat)
!     enddo

!     self%rCut = aenet_Rc_max
!     self%rCutSq = aenet_Rc_max * aenet_Rc_max
     allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
     allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
     self%rMin = 0E0_dp
     self%rMinTable = 0E0_dp

     self%rCut = 4.0E0_dp
     self%rCutSq = 4.0E0_dp**2


  end subroutine
!==========================================================================+
  subroutine DetailedECalc_AENet(self, curbox, E_T, accept)
!    use boxData, only: self%boxArray
    use ClassyConstants, only: pi
    use ParallelVar, only: nout
    implicit none
    class(AENet), intent(inout) :: self
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
!    self%atomTypes = 1

    xmax = 0E0_dp
    ymax = 0E0_dp
    zmax = 0E0_dp
    do iAtom = 1, curbox%nMaxAtoms
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
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
     call get_energy_lib(self%box(1:3, 1:3), nCurAtoms, &
                         self%tempcoords(1:3, 1:nCurAtoms), &
                         self%atomTypes(1:nCurAtoms), pbc, &
                         Ecoh, E_T)

     write(nout, *) "Raw Neuro Net Energy:", E_T, Ecoh
!     E_T = E_T * curbox%nMolTotal / boltz
     E_T = E_T / boltz
     write(nout, *) "Total Neuro Net Energy:", E_T

  end subroutine
!============================================================================
  subroutine DiffECalc_AENet(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    use ClassyConstants, only: pi
    implicit none
    class(AENet), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
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

    accept = .true.
    E_Diff = 0E0_dp

    nTotalMol = curbox%nMolTotal
    select type(disp)
      class is(OrthoVolChange)
        xscale = disp(1)%xScale
        yscale = disp(1)%yScale
        zscale = disp(1)%zScale

      class default
        xscale = 1E0_dp
        yscale = 1E0_dp
        zscale = 1E0_dp

    end select

    !Check the rMin criteria first to ensure there is no overlap prior to
    !passing the configuration to AENet
    select type(disp)
      !-----------------------------------------------------
      class is(Displacement)
        do iDisp = 1, size(disp)
          iAtom = disp(iDisp)%atmIndx
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            jAtom = curbox%NeighList(1)%list(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)
            rmin_ij = self % rMinTable(atmType2, atmType1)          

            rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
            ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
            rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < rmin_ij) then
                accept = .false.
                return
            endif 
          enddo
        enddo
      !-----------------------------------------------------
      class is(Addition)
        do iDisp = 1, size(disp)
          iAtom = disp(iDisp)%atmIndx
          atmType1 = curbox % AtomType(iAtom)
          do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            jAtom = curbox%NeighList(1)%list(jNei, iAtom)
            atmType2 = curbox % AtomType(jAtom)
            rmin_ij = self % rMinTable(atmType2, atmType1)          

            rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
            ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
            rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
            call curbox%Boundary(rx, ry, rz)
            rsq = rx*rx + ry*ry + rz*rz
            if(rsq < rmin_ij) then
                accept = .false.
                return
            endif 
          enddo
        enddo
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
            rx = curbox %atoms(1, iAtom) + dx  -  curbox % atoms(1, jAtom) - dxj
            ry = curbox %atoms(2, iAtom) + dy  -  curbox % atoms(2, jAtom) - dyj
            rz = curbox %atoms(3, iAtom) + dz  -  curbox % atoms(3, jAtom) - dzj
            rsq = rx*rx + ry*ry + rz*rz
            rmin_ij = self % rMinTable(atmType2, atmType1)      
            if(rsq < rmin_ij) then
              accept = .false.
              return
            endif 
          enddo
        enddo


    end select

    ! Convert the Classy Array into an Array that AENet can read.
    nCurAtoms = 0
    self%tempcoords = 0E0_dp
    self%atomTypes = 1
    select type(disp)
!       -----------------------------------------------------
      class is(Displacement)
        do iAtom = 1, curbox%nMaxAtoms
          if(curbox%MolIndx(iAtom) == disp(1)%molIndx) then
            do iDisp = 1, size(disp)
              if(disp(iDisp)%atmIndx == iAtom) then
                nCurAtoms = nCurAtoms + 1
                self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
                self%tempcoords(1, nCurAtoms) = disp(iDisp)%x_new
                self%tempcoords(2, nCurAtoms) = disp(iDisp)%y_new
                self%tempcoords(3, nCurAtoms) = disp(iDisp)%z_new
                exit
              endif
            enddo

          elseif( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
            cycle
          else
            nCurAtoms = nCurAtoms + 1
            self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
            self%tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
            self%tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
            self%tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)
          endif
        enddo
!       -----------------------------------------------------
      class is(Addition)
        nTotalMol = nTotalMol + 1
        do iAtom = 1, curbox%nMaxAtoms
          if(curbox%MolIndx(iAtom) == disp(1)%molIndx) then
            do iDisp = 1, size(disp)
              if(disp(iDisp)%atmIndx == iAtom) then
                nCurAtoms = nCurAtoms + 1
                self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
                self%tempcoords(1, nCurAtoms) = disp(iDisp)%x_new
                self%tempcoords(2, nCurAtoms) = disp(iDisp)%y_new
                self%tempcoords(3, nCurAtoms) = disp(iDisp)%z_new
                exit
              endif
            enddo

          else if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
            cycle
          else
            nCurAtoms = nCurAtoms + 1
            self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
            self%tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
            self%tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
            self%tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)
          endif
        enddo
!       -----------------------------------------------------
      class is(Deletion)
        nTotalMol = nTotalMol - 1
        do iAtom = 1, curbox%nMaxAtoms
          if(curbox%MolIndx(iAtom) == disp(1)%molIndx) then
            cycle
          else if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
            cycle
          else
            nCurAtoms = nCurAtoms + 1
            self%atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
            self%tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
            self%tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
            self%tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)
          endif
        enddo
!       -----------------------------------------------------
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
          self%tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom) + dx
          self%tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom) + dy
          self%tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom) + dz
        enddo
!       -----------------------------------------------------
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

     if(nCurAtoms == 0) then
       E_Diff = -curbox%Etotal
       return
     endif

     self%boxrecp = 0E0_dp
     self%boxrecp(:,:) = geo_recip_lattice(self%box)
     self%tempcoords(1:3,1:nCurAtoms) = matmul(self%boxrecp(1:3,1:3), self%tempcoords(1:3,1:nCurAtoms))/ (2E0_dp * pi)
!     do iAtom = 1, 3
!       write(*,*) self%box(1:3, iAtom)
!     enddo
!     do iAtom = 1, nCurAtoms
!       write(*,*) self%tempcoords(1:3, iAtom)
!     enddo
!     do iAtom = 1, nCurAtoms
!       self%tempcoords(1, iAtom) = self%tempcoords(1, iAtom)/self%box(1,1)
!       self%tempcoords(2, iAtom) = self%tempcoords(2, iAtom)/self%box(2,1)
!       self%tempcoords(3, iAtom) = self%tempcoords(3, iAtom)/self%box(3,1)
!     enddo
     call get_energy_lib(self%box(1:3,1:3), nCurAtoms, self%tempcoords(1:3, 1:nCurAtoms), self%atomTypes(1:nCurAtoms), pbc, Ecoh, E_T)
!     call get_energy_lib(self%box, nCurAtoms, &
!                         self%tempcoords, &
!                         self%atomTypes, pbc, &
!                         Ecoh, E_T)


!     E_T = E_T * nTotalMol / boltz
     E_T = E_T / boltz
!     E_T = Ecoh * nTotalMol / boltz
     E_Diff = E_T - curbox%ETotal


  end subroutine
!=============================================================================+
!End -DAENET safety block
#endif
!=============================================================================+
  subroutine ProcessIO_AENet(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    implicit none
    class(AENet), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, rCut, rMin
  

!    call GetXCommand(line, command, 1, lineStat)
!    select case(trim(adjustl(command)))
!      case("rcut")
!        call GetXCommand(line, command, 2, lineStat)
!        read(command, *) rCut
!        self % rCut = rCut
!        self % rCutSq = rCut * rCut
!      case default
!        param = .true.
!    end select


!    if(param) then
!    endif
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
        read(line, *) type1, type2, rMin
        self%rMinTable(type1, type2) = rMin**2
        self%rMinTable(type2, type1) = rMin**2
      case default
        lineStat = -1
    end select


  end subroutine
!=============================================================================+
  subroutine Prologue_AENet(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes, nAtomTypes
    implicit none
    class(AENet), intent(inout) :: self
    integer :: atomLimit, iBox
      atomLimit = 0
     do iBox = 1, size(BoxArray)
       if( atomLimit < boxArray(iBox) % box % GetMaxAtoms()) then
         atomLimit = boxArray(iBox) % box % GetMaxAtoms()
       endif
     enddo
     allocate(self%atomTypes(1:atomLimit))
     allocate(self%tempcoords(3, 1:atomLimit))


  end subroutine
!=============================================================================+
end module
!=============================================================================+
