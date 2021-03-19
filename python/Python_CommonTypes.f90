!=========================================================================
! Module for creating Python wrappers for common Classy data strucutres.
! Used
!=========================================================================

module ClassyPyObj
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use SimpleSimBox, only: SimpleBox
  use VarPrecision
#ifdef EMBPYTHON
#define errcheck if(ierror/=0) then;call err_print;stop;endif
  use forpy_mod, only: dict, dict_create, list, list_create, ndarray, ndarray_create, err_print, ndarray_create_nocopy
  contains
!=========================================================================
  function createboxdict(boxnum) result(boxdict)
    use BoxData, only: BoxArray
    use Input_Format, only: ReplaceText


    implicit none
    integer, intent(in) :: boxnum
    type(dict) :: boxdict

    integer :: ierror

    type(ndarray) :: np_atoms, np_boxdim, np_chempot
    type(ndarray) :: np_nmol, np_ETable, np_atomtype

    character(len=30) :: dictkey, tempstr
    real(dp), allocatable :: boxDim(:,:)



    !First we need to share the atom data arrays with Python by passing
    !a dictionary with a numpy reference to the Classy side arrays

    ierror = dict_create(boxdict) 
    errcheck

    select type(box => BoxArray(boxnum)%box)
      type is(SimpleBox)
          ierror = boxdict%setitem("boxtype", "nobox")
      class is(CubeBox)
          ierror = boxdict%setitem("boxtype", "cube")
      class is(OrthoBox)
          ierror = boxdict%setitem("boxtype", "ortho")
      class default
          ierror = boxdict%setitem("boxtype", "unknown")
    end select

    ierror = boxdict%setitem("temperature", BoxArray(boxnum)%box%temperature)

    ierror = boxdict%setitem("temperature", BoxArray(boxnum)%box%temperature)
    ierror = boxdict%setitem("pressure", BoxArray(boxnum)%box%pressure)
    ierror = boxdict%setitem("volume", BoxArray(boxnum)%box%volume)
    if(.not. allocated(boxDim)) then
      allocate( boxDim(1:2, 1:BoxArray(boxNum)%box%nDimension) )
    endif
    call BoxArray(boxNum)%box%GetDimensions(boxdim)
    ierror = ndarray_create(np_boxdim, boxDim)
    errcheck
    ierror = boxdict%setitem("boxdimensions", np_boxdim)

    ierror = ndarray_create(np_chempot, BoxArray(boxnum)%box%chempot)
    errcheck
    ierror = boxdict%setitem("chemicalpotential", np_chempot)

    ierror = ndarray_create(np_atomtype, BoxArray(boxnum)%box%AtomType)
    errcheck
    ierror = boxdict%setitem("atomtype", np_atomtype)

    ierror = ndarray_create(np_atoms, BoxArray(boxnum)%box%atoms)
    errcheck
    ierror = boxdict%setitem("raw_atoms", np_atoms)

    ierror = ndarray_create(np_nmol,  BoxArray(boxnum)%box%NMol)
    errcheck
    ierror = boxdict%setitem("moleculecount", np_nmol)

    ierror = ndarray_create(np_ETable,  BoxArray(boxnum)%box%ETable)
    errcheck
    ierror = boxdict%setitem("energytable", np_ETable)


  end function
!=========================================================================
  function createboxdict_nocopy(boxnum) result(boxdict)
    use BoxData, only: BoxArray
    use Input_Format, only: ReplaceText
    implicit none
    integer, intent(in) :: boxnum
    type(dict) :: boxdict

    integer :: ierror

    type(ndarray) :: np_atoms, np_boxdim, np_chempot
    type(ndarray) :: np_nmol, np_ETable, np_atomtype

    character(len=30) :: dictkey, tempstr
    real(dp), asynchronous, allocatable :: boxDim(:,:)



    !First we need to share the atom data arrays with Python by passing
    !a dictionary with a numpy reference to the Classy side arrays

    ierror = dict_create(boxdict) 
    errcheck
    select type(box => BoxArray(boxnum)%box)
      type is(SimpleBox)
          ierror = boxdict%setitem("boxtype", "nobox")
      class is(CubeBox)
          ierror = boxdict%setitem("boxtype", "cube")
      class is(OrthoBox)
          ierror = boxdict%setitem("boxtype", "ortho")
      class default
          ierror = boxdict%setitem("boxtype", "unknown")
    end select

    ierror = boxdict%setitem("temperature", BoxArray(boxnum)%box%temperature)
    ierror = boxdict%setitem("pressure", BoxArray(boxnum)%box%pressure)
    ierror = boxdict%setitem("volume", BoxArray(boxnum)%box%volume)
    if(.not. allocated(boxDim)) then
      allocate( boxDim(1:2, 1:BoxArray(boxNum)%box%nDimension) )
    endif
    call BoxArray(boxNum)%box%GetDimensions(boxdim)
    ierror = ndarray_create(np_boxdim, boxDim)
    errcheck
    ierror = boxdict%setitem("boxdimensions", np_boxdim)

    ierror = ndarray_create_nocopy(np_chempot, BoxArray(boxnum)%box%chempot)
    errcheck
    ierror = boxdict%setitem("chemicalpotential", np_chempot)

    ierror = ndarray_create_nocopy(np_atomtype, BoxArray(boxnum)%box%AtomType)
    errcheck
    ierror = boxdict%setitem("atomtype", np_atomtype)

    ierror = ndarray_create_nocopy(np_atoms, BoxArray(boxnum)%box%atoms)
    errcheck
    ierror = boxdict%setitem("raw_atoms", np_atoms)

    ierror = ndarray_create_nocopy(np_nmol,  BoxArray(boxnum)%box%NMol)
    errcheck
    ierror = boxdict%setitem("moleculecount", np_nmol)

    ierror = ndarray_create_nocopy(np_ETable,  BoxArray(boxnum)%box%ETable)
    errcheck
    ierror = boxdict%setitem("energytable", np_ETable)


  end function

!=========================================================================
  function createdisplist(disp) result(displist)
    use BoxData, only: BoxArray
    use Input_Format, only: ReplaceText
    use CoordinateTypes, only: Perturbation, Displacement, Deletion, Addition, VolChange,&
                               OrthoVolChange
    implicit none
    type(list) :: displist
    class(Perturbation) :: disp(:)
    integer :: ierror, iDisp, nDisp
    type(dict) :: dispdict(1:size(disp))


    nDisp = size(disp)
    ierror = list_create(displist) 
    errcheck
    do iDisp = 1, nDisp
      ierror = dict_create( dispdict(iDisp) )
      errcheck
      select type(disp)
        class is(Displacement)
          ierror = dispdict(iDisp)%setitem("type", "displacement")
          ierror = dispdict(iDisp)%setitem("moltype", disp(iDisp)%molType)
          ierror = dispdict(iDisp)%setitem("molindex", disp(iDisp)%molIndx)
!          ierror = dispdict(iDisp)%setitem("atomsubindex", disp(iDisp)%atmSubIndx)
          ierror = dispdict(iDisp)%setitem("atomindex", disp(iDisp)%atmIndx)
          ierror = dispdict(iDisp)%setitem("x_new", disp(iDisp)%x_new)
          ierror = dispdict(iDisp)%setitem("y_new", disp(iDisp)%y_new)
          ierror = dispdict(iDisp)%setitem("z_new", disp(iDisp)%z_new)

        class is(Deletion)
          ierror = dispdict(iDisp)%setitem("type", "deletion")
          ierror = dispdict(iDisp)%setitem("moltype", disp(iDisp)%molType)
          ierror = dispdict(iDisp)%setitem("molindex", disp(iDisp)%molIndx)
          ierror = dispdict(iDisp)%setitem("atomindex", disp(iDisp)%atmIndx)

        class is(Addition)
          ierror = dispdict(iDisp)%setitem("type", "addition")
          ierror = dispdict(iDisp)%setitem("moltype", disp(iDisp)%molType)
          ierror = dispdict(iDisp)%setitem("molindex", disp(iDisp)%molIndx)
          ierror = dispdict(iDisp)%setitem("atomindex", disp(iDisp)%atmIndx)
          ierror = dispdict(iDisp)%setitem("x_new", disp(iDisp)%x_new)
          ierror = dispdict(iDisp)%setitem("y_new", disp(iDisp)%y_new)
          ierror = dispdict(iDisp)%setitem("z_new", disp(iDisp)%z_new)

        class is(VolChange)
          ierror = dispdict(iDisp)%setitem("type", "volchange")
          ierror = dispdict(iDisp)%setitem("volnew", disp(iDisp)%volnew)
          ierror = dispdict(iDisp)%setitem("volold", disp(iDisp)%volold)

        class is(OrthoVolChange)
          ierror = dispdict(iDisp)%setitem("type", "orthovolchange")
          ierror = dispdict(iDisp)%setitem("xscale", disp(iDisp)%xScale)
          ierror = dispdict(iDisp)%setitem("yscale", disp(iDisp)%yScale)
          ierror = dispdict(iDisp)%setitem("zscale", disp(iDisp)%zScale)
          ierror = dispdict(iDisp)%setitem("volnew", disp(iDisp)%volnew)
          ierror = dispdict(iDisp)%setitem("volold", disp(iDisp)%volold)

        class default
          error stop "Unable to cast Python Type using this perturbation type"
      end select
      ierror = displist%append( dispdict(iDisp) )
    enddo




  end function

!========================================================================
#endif
end module
!=========================================================================
