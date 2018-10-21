!==========================================================================================
! The purpose of this module is to provide the base class for the molecule
! constructor class.
!==========================================================================================
module MolCon_RidgidRegrowth
  use CoordinateTypes, only: Displacement, Perturbation
  use MasterTemplate, only: classyClass
  use VarPrecision


  type, public, extends(MolConstructor) :: RidgidRegrowth
    real(dp), allocatable :: tempcoords(:, :)
    contains
      procedure, public, pass :: Constructor => SimpleRegrwoth_Constructor
      procedure, public, pass :: GenerateConfig => SimpleRegrwoth_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine RidgidRegrowth_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(MolConstructor), intent(inout) :: self
    integer :: iType, maxAtoms

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo

    allocate( tempcoords(1:3, 1:maxAtoms) )

  end subroutine
!==========================================================================================
  subroutine RidgidRegrowth_GenerateConfig(self, trialBox, disp, probconstruct)
    use Constants, only: two_pi
    use Common_MolInfo, only: MolData
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimpleBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct = 1E0_dp
    real(dp) :: c_term, s_term
    real(dp) :: x_shift, y_shift, z_shift
    real(dp) :: x1, y1, z1
    real(dp) :: x_rid_cm, y_rid_cm, z_rid_cm

    x1 = tempcoords(1, 1)
    y1 = tempcoords(2, 1)
    z1 = tempcoords(3, 1)
    do i = 1, 
      tempcoords(1, i) = tempcoords(1, i) - x1
      tempcoords(2, i) = tempcoords(2, i) - y1
      tempcoords(3, i) = tempcoords(3, i) - z1
    enddo
    if(nAtoms(nType) .ne. 1) then
      x_rid_cm = gasConfig(nType)%x(1)
      y_rid_cm = gasConfig(nType)%y(1)
      z_rid_cm = gasConfig(nType)%z(1)    
    
!          Rotate xz axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do i = 1,nAtoms(nType)
        x_shift = tempcoords(1, i) - x_rid_cm
        z_shift = tempcoords(3, i) - z_rid_cm        
        tempcoords(1, i) = c_term*x_shift - s_term*z_shift + x_rid_cm
        tempcoords(3, i) = s_term*x_shift + c_term*z_shift + z_rid_cm
      enddo   
    
!          Rotate yz axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do i = 1,nAtoms(nType)
        y_shift = tempcoords(2, i) - y_rid_cm
        z_shift = tempcoords(3, i) - z_rid_cm        
        tempcoords(2, i) = c_term*y_shift - s_term*z_shift + y_rid_cm
        tempcoords(3, i) = s_term*y_shift + c_term*z_shift + z_rid_cm
      enddo   

!          Rotate xy axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do i = 1,nAtoms(nType)
        x_shift = tempcoords(1, i) - x_rid_cm
        y_shift = tempcoords(2, i) - y_rid_cm        
        tempcoords(1, i) = c_term*x_shift - s_term*y_shift + x_rid_cm
        tempcoords(2, i) = s_term*x_shift + c_term*y_shift + y_rid_cm
      enddo    
    endif


  end subroutine
!==========================================================================================
end module
!==========================================================================================
