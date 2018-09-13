!==========================================================================================
! The purpose of this module is to provide the base class for the molecule
! constructor class.
!==========================================================================================
module MolCon_SimpleRegrowth
  use CoordinateTypes, only: Displacement, Perturbation
  use MasterTemplate, only: classyClass
  use VarPrecision


  type, public, extends(MolConstructor) :: SimpleRegrowth
    contains
      procedure, public, pass :: Constructor => SimpleRegrwoth_Constructor
      procedure, public, pass :: GenerateConfig => SimpleRegrwoth_GenerateConfig
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine SimpleRegrowth_Constructor(self)
    implicit none
    class(MolConstructor), intent(inout) :: self
  end subroutine
!==========================================================================================
  subroutine SimpleRegrowth_GenerateConfig(self, trialBox, disp, probconstruct)
    implicit none
    class(MolConstructor), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimpleBox), intent(inout) :: trialBox
    real(dp), intent(out) :: probconstruct = 1E0_dp

      x1 = rosenTrial(iRosen)%x(1)
      y1 = rosenTrial(iRosen)%y(1)
      z1 = rosenTrial(iRosen)%z(1)
      do i = 1, nAtoms(nType)
        rosenTrial(iRosen)%x(i) = rosenTrial(iRosen)%x(i) - x1
        rosenTrial(iRosen)%y(i) = rosenTrial(iRosen)%y(i) - y1
        rosenTrial(iRosen)%z(i) = rosenTrial(iRosen)%z(i) - z1
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
          x_shift = rosenTrial(iRosen)%x(i) - x_rid_cm
          z_shift = rosenTrial(iRosen)%z(i) - z_rid_cm        
          rosenTrial(iRosen)%x(i) = c_term*x_shift - s_term*z_shift + x_rid_cm
          rosenTrial(iRosen)%z(i) = s_term*x_shift + c_term*z_shift + z_rid_cm
        enddo   
      
!          Rotate yz axis
        rotang = two_pi*grnd()
        c_term = cos(rotang)
        s_term = sin(rotang)
        do i = 1,nAtoms(nType)
          y_shift = rosenTrial(iRosen)%y(i) - y_rid_cm
          z_shift = rosenTrial(iRosen)%z(i) - z_rid_cm        
          rosenTrial(iRosen)%y(i) = c_term*y_shift - s_term*z_shift + y_rid_cm
          rosenTrial(iRosen)%z(i) = s_term*y_shift + c_term*z_shift + z_rid_cm
        enddo   

!          Rotate xy axis
        rotang = two_pi*grnd()
        c_term = cos(rotang)
        s_term = sin(rotang)
        do i = 1,nAtoms(nType)
          x_shift = rosenTrial(iRosen)%x(i) - x_rid_cm
          y_shift = rosenTrial(iRosen)%y(i) - y_rid_cm        
          rosenTrial(iRosen)%x(i) = c_term*x_shift - s_term*y_shift + x_rid_cm
          rosenTrial(iRosen)%y(i) = s_term*x_shift + c_term*y_shift + y_rid_cm
        enddo    
      endif


  end subroutine
!==========================================================================================
end module
!==========================================================================================
