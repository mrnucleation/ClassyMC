!===============================
module BoxPresets
use VarPrecision

!===============================
contains
!===============================
  subroutine Preset(iLine, lineStore, boxNum, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none

    integer, intent(in) :: iLine, boxNum
    character(len=maxLineLen), intent(in) :: lineStore(:)   
    integer, intent(out) :: lineStat

    character(len=30) :: command
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(lineStore(iLine), command, 2, lineStat)

    select case(trim(adjustl(command)))
      case("fcc")
        call GetXCommand(lineStore(iLine), command, 3, lineStat)
        read(command, *) realVal
        call GetXCommand(lineStore(iLine), command, 4, lineStat)
        read(command, *) intVal
        call FCC(BoxNum, realVal, intVal)

      case default
        lineStat = -1

    end select


  end subroutine
!===============================
  subroutine FCC(BoxNum, latConst, replicate)
    use BoxData, only: BoxArray
    use CubicBoxDef, only: CubeBox
    use OrthoBoxDef, only: OrthoBox
    use ParallelVar, only: nout

    implicit none
    integer, intent(in) :: BoxNum, replicate
    real(dp), intent(in) :: latConst

    integer :: nAtoms, i, j, k, iAtom, iUnit
    real(dp) :: boxL
    real(dp) :: baseCell(1:4, 1:3)

    write(nout,*) "Box", boxNum, "Creating FCC Structure with parameters:", latConst
!    allocate( OrthoBox::BoxArray(BoxNum)%box )
    allocate( CubeBox::BoxArray(BoxNum)%box )
    BoxArray(boxNum)%box%boxStr = "ortho"


    baseCell(1, 1) = 0.0
    baseCell(1, 2) = 0.0
    baseCell(1, 3) = 0.0

    baseCell(2, 1) = 0.5*latConst
    baseCell(2, 2) = 0.5*latConst
    baseCell(2, 3) = 0.0

    baseCell(3, 1) = 0.5*latConst
    baseCell(3, 2) = 0.0
    baseCell(3, 3) = 0.5*latConst

    baseCell(4, 1) = 0.0
    baseCell(4, 2) = 0.5*latConst
    baseCell(4, 3) = 0.5*latConst


    boxL = replicate * latConst
    nAtoms = 4*replicate**3

    select type(box => BoxArray(BoxNum)%box)
      class is(OrthoBox)
        box%boxLx = boxL
        box%boxLy = boxL
        box%boxLz = boxL
        box%boxLx2 = boxL*0.5E0_dp
        box%boxLy2 = boxL*0.5E0_dp
        box%boxLz2 = boxL*0.5E0_dp

      class is(CubeBox)
        box%boxL = boxL
        box%boxL2 = boxL*0.5E0_dp

        call box%AllocateMolBound

        box%NMol = 0
        box%NMolMax = 0

        box%NMol(1) = nAtoms
        box%NMolMin(1) = nAtoms
        box%NMolMax(1) = nAtoms
        box%nAtoms = nAtoms

        call box%Constructor

        iAtom = 0
        do i = 1, replicate
          do j = 1, replicate
            do k = 1, replicate
              do iUnit = 1,4
                iAtom = iAtom + 1
                box%atoms(1, iAtom) = baseCell(iUnit, 1) + (i-1) * latConst
                box%atoms(2, iAtom) = baseCell(iUnit, 2) + (j-1) * latConst
                box%atoms(3, iAtom) = baseCell(iUnit, 3) + (k-1) * latConst
              enddo
            enddo
          enddo
        enddo

        do iAtom = 1, nAtoms
          do while(abs(box%atoms(1,iAtom)) > boxL*0.5E0_dp)  
            box%atoms(1, iAtom) = box%atoms(1, iAtom) - sign(boxL, box%atoms(1, iAtom))
          enddo
          do while(abs(box%atoms(2,iAtom)) > boxL*0.5E0_dp)  
            box%atoms(2, iAtom) = box%atoms(2, iAtom) - sign(boxL, box%atoms(2, iAtom))
          enddo
          do while(abs(box%atoms(3,iAtom)) > boxL*0.5E0_dp)  
            box%atoms(3, iAtom) = box%atoms(3, iAtom) - sign(boxL, box%atoms(3, iAtom))
          enddo
          write(*,*) "Ar", box%atoms(1, iAtom), box%atoms(2, iAtom), box%atoms(3, iAtom)

        enddo

    end select



  end subroutine
!===============================
end module
!===============================
