!===========================================================================
#define __StdErr__ 0
!===========================================================================
module SimLibrary
  use ParallelVar, only: myid, ierror, nout
use VarPrecision
!===========================================================================
contains
!===========================================================================
  subroutine Library_FullSimulation bind(C,name='Full')
    use SimMonteCarlo, only: RunMonteCarlo
    implicit none

    call RunMonteCarlo
 
  end subroutine
!===========================================================================
  subroutine Library_RunMove bind(C,name='Full')
    use SimMonteCarlo, only: RunMonteCarlo
    implicit none

 
  end subroutine
!===========================================================================
end module
!===========================================================================
