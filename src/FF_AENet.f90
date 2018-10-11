module FF_AENet
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use VarPrecision
  use CoordinateTypes

  !AENet Library functions
  use predict_lib, only: initialize_lib, get_energy_lib
  use geometry, only: geo_recip_lattice
  use input, only: InputData



  type, public, extends(forcefield) :: AENet
!    real(dp) :: rCut, rCutSq
    logical, parameter :: pbc = .false.
    logical :: initialized = .false.
    integer, allocatable :: atomTypes(:)
    real(dp), allocatable :: tempcoords(:,:)
    real(dp)  :: box(3,3)
    real(dp)  :: boxrecp(3,3)
    real(dp), parameter  :: boltz = 8.6173303E-5_dp
    contains
      procedure, pass :: Constructor => Constructor_AENet
      procedure, pass :: DetailedECalc => DetailedECalc_AENet
      procedure, pass :: DiffECalc => DiffECalc_AENet
      procedure, pass :: ProcessIO => ProcessIO_AENet
!      procedure, pass :: GetCutOff
  end type

  contains
!=============================================================================+
  subroutine Constructor_AENet(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(AENet), intent(inout) :: self
    type(InputData) :: indata
    character(len=80) :: str1, str2

    allocate(atomTypes(1:atomLimit))
    allocate(tempcoords(3, 1:atomLimit))
    call initialize_lib(str1, str2, indata)
    write(*,*) str1
    write(*,*) str2


  end subroutine
!=============================================================================+
  subroutine DetailedECalc_AENet(self, curbox, E_T, accept)
    implicit none
    class(AENet), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept
    integer 
    real(dp) :: xoffset, yoffset, zoffset


    accept = .true.
    xoffset = 0E0_dp
    yoffset = 0E0_dp
    zoffset = 0E0_dp
    nCurAtoms = 0
    box = 0E0_dp

    xmax = 0E0_dp
    ymax = 0E0_dp
    zmax = 0E0_dp
    do iType = 1, nMolTypes
      nCurAtoms = nCurAtoms + 1
      atomTypes(nCurAtoms) = atomArray(iType, iAtom)
      tempcoords(1, nCurAtoms) = MolArray(iType)%mol(iMol)%x(iAtom)
      tempcoords(2, nCurAtoms) = MolArray(iType)%mol(iMol)%y(iAtom)
      tempcoords(3, nCurAtoms) = MolArray(iType)%mol(iMol)%z(iAtom)

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
!      write(3,*) nCurAtoms
!      write(3,*) box
!      do iAtom = 1, nCurAtoms
!        write(3,*) tempcoords(1:3, iAtom)
!      enddo
      call get_energy_lib(box, nCurAtoms, tempcoords(1:3, 1:nCurAtoms), atomTypes(1:nCurAtoms), pbc, Ecoh, E_T)

      E_T = E_T*NPART(1)/boltz
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
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

      xoffset = 0E0_dp
      yoffset = 0E0_dp
      zoffset = 0E0_dp
      nCurAtoms = 0
      box = 0E0_dp

      xmax = 0E0_dp
      ymax = 0E0_dp
      zmax = 0E0_dp
      do iType = 1, nMolTypes
        do iMol=1,NPART(iType)
          do iAtom = 1,nAtoms(iType)
            nCurAtoms = nCurAtoms + 1
!            if(atomArray(iType, iAtom)) then
!              atomTypes(nCurAtoms) = 1
!            else
!              atomTypes(nCurAtoms) = 1
!            endif
            atomTypes(nCurAtoms) = atomArray(iType, iAtom)


            tempcoords(1, nCurAtoms) = MolArray(iType)%mol(iMol)%x(iAtom)
            tempcoords(2, nCurAtoms) = MolArray(iType)%mol(iMol)%y(iAtom)
            tempcoords(3, nCurAtoms) = MolArray(iType)%mol(iMol)%z(iAtom)

            if(xoffset > tempcoords(1, nCurAtoms)) then
              xoffset = tempcoords(1, nCurAtoms)
            endif
            if(yoffset > tempcoords(2, nCurAtoms)) then
              yoffset = tempcoords(2, nCurAtoms)
            endif
            if(zoffset > tempcoords(3, nCurAtoms)) then
              zoffset = tempcoords(3, nCurAtoms)
            endif
          enddo
        enddo
      enddo



      do iAtom = 1, nCurAtoms
        tempcoords(1, iAtom) = tempcoords(1, iAtom) - xoffset
        tempcoords(2, iAtom) = tempcoords(2, iAtom) - yoffset
        tempcoords(3, iAtom) = tempcoords(3, iAtom) - zoffset
        if(xmax < tempcoords(1, iAtom)) then
          xmax = tempcoords(1, iAtom)
        endif
        if(ymax < tempcoords(2, iAtom)) then
          ymax = tempcoords(2, iAtom)
        endif
        if(zmax < tempcoords(3, iAtom)) then
          zmax = tempcoords(3, iAtom)
        endif
      enddo
      box(1,1) = xmax + 1e-5
      box(2,2) = ymax + 1e-5
      box(3,3) = zmax + 1e-5
      boxrecp(:,:) = geo_recip_lattice(box)
      tempcoords(1:3,1:nCurAtoms) = matmul(boxrecp, tempcoords)/(2.0d0*pi)
!      write(3,*) nCurAtoms
!      write(3,*) box
!      do iAtom = 1, nCurAtoms
!        write(3,*) tempcoords(1:3, iAtom)
!      enddo
      call get_energy_lib(box, nCurAtoms, tempcoords(1:3, 1:nCurAtoms), atomTypes(1:nCurAtoms), pbc, Ecoh, E_T)

      E_T = E_T*NPART(1)/boltz
      write(nout, *) "Total Neuro Net Energy:", E_T, Ecoh
  end subroutine
!=============================================================================+
  subroutine ProcessIO_AENet(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(AENet), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

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
