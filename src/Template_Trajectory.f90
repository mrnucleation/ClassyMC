!====================================================================
module TrajectoryTemplate
  use MasterTemplate, only: classyClass
  use VarPrecision

  type, public, extends(classyClass) :: trajectory
    integer :: fileUnit = -1
    integer :: boxNum = -1
    integer :: outFreq = 5000
    character(len=50) :: fileName
    contains
       procedure, pass :: SetUnit
       procedure, pass :: SetBox
       procedure, pass :: SetFileName
       procedure, pass :: SetFreq
       procedure, pass :: OpenFile
       procedure, pass :: WriteFrame
       procedure, pass :: CloseFile
  end type
!====================================================================
  contains
!====================================================================
  subroutine SetUnit(self, fileUnit) 
    implicit none
    class(trajectory), intent(inout) :: self
    integer, intent(in) :: fileUnit

    self%fileUnit = fileUnit

  end subroutine
!====================================================================
  subroutine SetFreq(self, freq) 
    implicit none
    class(trajectory), intent(inout) :: self
    integer, intent(in) :: freq

    self%outFreq = freq

  end subroutine
!====================================================================
  subroutine SetBox(self, boxNum) 
    implicit none
    class(trajectory), intent(inout) :: self
    integer, intent(in) :: boxNum

    self%boxNum = boxNum

  end subroutine
!====================================================================
  subroutine SetFileName(self, fileName) 
    implicit none
    class(trajectory), intent(inout) :: self
    character(len=*), intent(in) :: fileName
    
    self%fileName = ""
    self%fileName(1:len(fileName)) = fileName(1:len(fileName))

  end subroutine
!====================================================================
  subroutine OpenFile(self) 
    implicit none
    class(trajectory), intent(in) :: self

    open( unit=self%fileUnit, file=trim(adjustl(self%filename)) )
  end subroutine
!====================================================================
  subroutine WriteFrame(self) 
    implicit none
    class(trajectory), intent(inout) :: self


  end subroutine
!====================================================================
  subroutine CloseFile(self) 
    implicit none
    class(trajectory), intent(in) :: self

    close(self%fileUnit)

  end subroutine
!====================================================================
end module
!====================================================================
