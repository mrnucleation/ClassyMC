!====================================================================
module SimStyleTemplate
  type, public :: SimStyle
    contains
       procedure, pass :: RunSimulation
  end type
!====================================================================
  contains
!====================================================================
  subroutine RunSimulation
    implicit none
  end subroutine
!====================================================================
end module
!====================================================================
