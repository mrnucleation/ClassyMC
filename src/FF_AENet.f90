!=============================================================================+
! AENet forcefield interface function. The AENet files are not currently
! included with the Classy code base and must be compiled separately.
! To add AENet to the compilation apply the -DAENET flag in the Makefile
!=============================================================================+
module FF_AENet
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use SimpleSimBox, only: SimpleBox
  use VarPrecision
  use CoordinateTypes

  !AENet Library functions
#ifdef AENET
  use predict_lib, only: initialize_lib, get_energy_lib
  use geometry, only: geo_recip_lattice
  use input, only: InputData
#endif

  real(dp), parameter :: boltz = 8.6173303E-5_dp

  type, public, extends(forcefield) :: AENet
!    real(dp) :: rCut, rCutSq
    logical :: initialized = .false.
    integer, allocatable :: atomTypes(:)
    real(dp), allocatable :: tempcoords(:,:)
    real(dp)  :: box(3,3)
    real(dp)  :: boxrecp(3,3)

    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    contains
#ifdef AENET
      procedure, pass :: Constructor => Constructor_AENet
      procedure, pass :: DetailedECalc => DetailedECalc_AENet
      procedure, pass :: DiffECalc => DiffECalc_AENet
#endif
      procedure, pass :: ProcessIO => ProcessIO_AENet
!      procedure, pass :: GetCutOff
  end type

  contains
#ifdef AENET
!=============================================================================+
  subroutine Constructor_AENet(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes
    use Constants, only: pi
    implicit none
    class(AENet), intent(inout) :: self
    integer :: iBox, atomLimit

    if(.not. self%initialized) then
      atomLimit = 0
      do iBox = 1, size(BoxArray)
        if( atomLimit < BoxArray(iBox) % box % nMaxAtoms) then
          atomLimit = BoxArray(iBox) % box % nMaxAtoms
        endif
      enddo
      allocate(atomTypes(1:atomLimit))
      allocate(tempcoords(3, 1:atomLimit))
      call initialize_lib(str1, str2, indata)
!      write(*,*) str1
!      write(*,*) str2
      self%initialized = .true.
    endif


  end subroutine
!==========================================================================+
  subroutine DetailedECalc_AENet(self, curbox, E_T, accept)
    use BoxData, only: BoxArray
    implicit none
    class(AENet), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    logical :: pbc = .false.
    integer :: iAtom
    integer :: nTotalMol
    type(InputData) :: indata
    character(len=80) :: str1, str2
    real(dp) :: xoffset, yoffset, zoffset
    real(dp) :: tempdim(1:3,1:3)


    accept = .true.
    xoffset = 0E0_dp
    yoffset = 0E0_dp
    zoffset = 0E0_dp
    nCurAtoms = 0
    box = 0E0_dp

    xmax = 0E0_dp
    ymax = 0E0_dp
    zmax = 0E0_dp
    do iAtom = 1, curbox%nMaxAtoms
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
        cycle
      endif
      nCurAtoms = nCurAtoms + 1
      atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
      tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
      tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
      tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)

      select type(curbox)
        class is(SimpleBox)
          if(xoffset > tempcoords(1, nCurAtoms)) then
            xoffset = tempcoords(1, nCurAtoms)
          endif
          if(yoffset > tempcoords(2, nCurAtoms)) then
            yoffset = tempcoords(2, nCurAtoms)
          endif
          if(zoffset > tempcoords(3, nCurAtoms)) then
            zoffset = tempcoords(3, nCurAtoms)
          endif
       end select

     enddo
     !Collect the box dimensions and determine the offset. 
     select type(curbox)
       class is(CubeBox)
         pbc = .true.
         call curbox%GetDimensions(tempdim)
         xoffset = tempdim(1, 1)
         yoffset = tempdim(1, 2)
         zoffset = tempdim(1, 3)
         box(1,1) = tempdim(2, 1) - tempdim(1, 1)
         box(2,2) = tempdim(2, 2) - tempdim(1, 2)
         box(3,3) = tempdim(2, 3) - tempdim(1, 3)
     end select

     do iAtom = 1, nCurAtoms
       tempcoords(1, iAtom) = tempcoords(1, iAtom) - xoffset
       tempcoords(2, iAtom) = tempcoords(2, iAtom) - yoffset
       tempcoords(3, iAtom) = tempcoords(3, iAtom) - zoffset
       select type(curbox)
         class is(SimpleBox)
           if(xmax < tempcoords(1, iAtom)) then
             xmax = tempcoords(1, iAtom)
           endif
           if(ymax < tempcoords(2, iAtom)) then
             ymax = tempcoords(2, iAtom)
           endif
           if(zmax < tempcoords(3, iAtom)) then
             zmax = tempcoords(3, iAtom)
           endif
       end select
     enddo
     select type(curbox)
       class is(SimpleBox)
         box(1,1) = xmax + 1e-5
         box(2,2) = ymax + 1e-5
         box(3,3) = zmax + 1e-5
     end select



      boxrecp(:,:) = geo_recip_lattice(box)
      tempcoords(1:3,1:nCurAtoms) = matmul(boxrecp, tempcoords)/(2.0d0*pi)
      call get_energy_lib(box, nCurAtoms, tempcoords(1:3, 1:nCurAtoms), atomTypes(1:nCurAtoms), pbc, Ecoh, E_T)

      E_T = E_T * curbox%nTotalMol / boltz
      write(nout, *) "Total Neuro Net Energy:", E_T, Ecoh

  end subroutine
!============================================================================
  subroutine DiffECalc_AENet(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(AENet), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    logical :: pbc
    integer :: iAtom, jATom
    integer :: atmType1, atmType2
    integer :: iDisp, nTotalMol
    real(dp) :: E_T
    real(dp) :: rmin_ij
    real(dp) :: rx, ry, rz, rsq

    accept = .true.
    E_Diff = 0E0_dp

    nTotalMol = curbox%nTotalMol

    !Check the rMin criteria first to ensure there is no overlap prior to
    !passing the configuration to AENet
    select type(disp)
      class is(DisplacementNew)
        do iDisp = 1, dispLen
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
      class is(Addition)
        do iDisp = 1, dispLen
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
    end select

    ! Convert the Classy Array into an Array that AENet can read.
    nCurAtoms = 0
    select case(disp)
!       -----------------------------------------------------
      class is(DisplacementNew)
        do iAtom = 1, curbox%nMaxAtoms
          if(curbox%MolIndx(iAtom) == disp(1)%molIndx) then
            do iDisp = 1, size(disp)
              if(disp(iDisp)%atmIndx == iAtom) then
                nCurAtoms = nCurAtoms + 1
                atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
                tempcoords(1, nCurAtoms) = disp(iDisp)%x_new
                tempcoords(2, nCurAtoms) = disp(iDisp)%y_new
                tempcoords(3, nCurAtoms) = disp(iDisp)%z_new
                exit
              endif
            enddo

          else if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
            cycle
          else:
            nCurAtoms = nCurAtoms + 1
            atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
            tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
            tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
            tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)
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
                atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
                tempcoords(1, nCurAtoms) = disp(iDisp)%x_new
                tempcoords(2, nCurAtoms) = disp(iDisp)%y_new
                tempcoords(3, nCurAtoms) = disp(iDisp)%z_new
                exit
              endif
            enddo

          else if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
            cycle
          else:
            nCurAtoms = nCurAtoms + 1
            atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
            tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
            tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
            tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)
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
          else:
            nCurAtoms = nCurAtoms + 1
            atomTypes(nCurAtoms) = curbox % AtomType(iAtom)
            tempcoords(1, nCurAtoms) = curBox%atoms(1, iAtom)
            tempcoords(2, nCurAtoms) = curBox%atoms(2, iAtom)
            tempcoords(3, nCurAtoms) = curBox%atoms(3, iAtom)
          endif
        enddo
    end select

    accept = .true.
    xoffset = 0E0_dp
    yoffset = 0E0_dp
    zoffset = 0E0_dp
    nCurAtoms = 0
    box = 0E0_dp

    xmax = 0E0_dp
    ymax = 0E0_dp
    zmax = 0E0_dp

     !Collect the box dimensions and determine the offset. 
     select type(curbox)
       class is(CubeBox)
         pbc = .true.
         call curbox%GetDimensions(tempdim)
         ! Box dimensions given back as xlow, xhigh, ylow, yhigh, etc.  
         ! Need to be converted to AENet format
         xoffset = tempdim(1, 1)
         yoffset = tempdim(1, 2)
         zoffset = tempdim(1, 3)
         box(1,1) = tempdim(2, 1) - tempdim(1, 1)
         box(2,2) = tempdim(2, 2) - tempdim(1, 2)
         box(3,3) = tempdim(2, 3) - tempdim(1, 3)
     end select

     do iAtom = 1, nCurAtoms
       tempcoords(1, iAtom) = tempcoords(1, iAtom) - xoffset
       tempcoords(2, iAtom) = tempcoords(2, iAtom) - yoffset
       tempcoords(3, iAtom) = tempcoords(3, iAtom) - zoffset
       select type(curbox)
         class is(SimpleBox)
           if(xmax < tempcoords(1, iAtom)) then
             xmax = tempcoords(1, iAtom)
           endif
           if(ymax < tempcoords(2, iAtom)) then
             ymax = tempcoords(2, iAtom)
           endif
           if(zmax < tempcoords(3, iAtom)) then
             zmax = tempcoords(3, iAtom)
           endif
       end select
     enddo
     select type(curbox)
       class is(SimpleBox)
         box(1,1) = xmax + 1e-5
         box(2,2) = ymax + 1e-5
         box(3,3) = zmax + 1e-5
     end select



     boxrecp(:,:) = geo_recip_lattice(box)
     tempcoords(1:3,1:nCurAtoms) = matmul(boxrecp, tempcoords)/(2.0d0*pi)
     call get_energy_lib(box, nCurAtoms, tempcoords(1:3, 1:nCurAtoms), atomTypes(1:nCurAtoms), pbc, Ecoh, E_T)

     E_T = E_T * nTotalMol / boltz
     E_Diff = curbox%ETotal - E_T
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
!  function GetCutOff(self) result(rCut)
!    implicit none
!    class(AENet), intent(inout) :: self
!    real(dp) :: rCut
!
!    rCut = self%rCut
!  end function
!=============================================================================+
end module
!=============================================================================+
