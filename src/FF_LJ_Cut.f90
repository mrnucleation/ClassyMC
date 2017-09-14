module Pair_LJ_Cut
  use ForceFieldTemplate
  use VarPrecision
  type, public :: Pair_LJ_Cut
    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: rMinTable(:,:)
    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: DetailedECalc 
      procedure, pass :: ShiftECalc
      procedure, pass :: SwapInECalc
      procedure, pass :: SwapOutECalc
      procedure, pass :: SetParameter
  end type

  contains
  !===================================================================================
  subroutine DetailedECalc(self, curbox)
    implicit none
    class(forcefield), intent(in) :: self
    type(simBox), intent(in) :: curbox
      real(dp), intent(inOut) :: E_T
      integer :: iAtom,jAtom
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: E_LJ
      real(dp) :: rmin_ij      

      E_LJ = 0E0
      PairList = 0E0      
      ETable = 0E0
      do iAtom = 1, curbox%
        atmType1 = atomArray(iType,iAtom)
        do jAtom = 1,nAtoms(jType)        
          atmType2 = atomArray(jType,jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = r_min_tab(atmType1,atmType2)          


          rx = curbox % atoms(0, iAtom)  -  curbox % atoms(0, jAtom)
          ry = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
          rz = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
          r = rx**2 + ry**2 + rz**2
          if(r .lt. rmin_ij) then
            stop "ERROR! Overlaping atoms found in the current configuration!"
          endif 
          LJ = (sig_sq/r)**3
          LJ = ep * LJ * (LJ-1E0)              
          E_LJ = E_LJ + LJ

          ETable(iIndx) = ETable(iIndx) + Ele + LJ
          ETable(jIndx) = ETable(jIndx) + Ele + LJ              
        enddo
      enddo

      
      write(nout,*) "Lennard-Jones Energy:", E_LJ
      
      E_T = E_T + E_Ele + E_LJ    
      E_Inter_T = E_Ele + E_LJ   
   end subroutine
  !=====================================================================
  subroutine ShiftCheck(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine
  !=====================================================================
  subroutine SwapInCheck(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine
  !=====================================================================
  subroutine SwapOutCheck(self)
    implicit none
    class(forcefield), intent(in) :: self
  end subroutine
  !=====================================================================
  subroutine SetParameter(self, parIndex,  parVal)
    implicit none
    class(forcefield), intent(in) :: self
    integer, intent(in) :: parIndex(:)
    real(dp), intent(in) :: parVal

    select case( parIndex(1) )
    case(1) !Epsilon
      epsTable(parIndex(2), parIndex(3)) = parVal
    case(2) !Sigma
      sigTable(parIndex(2), parIndex(3)) = parVal
    case(3) !rMin
      rMinTable(parIndex(2), parIndex(3)) = parVal
    case default
      write(*,*) "ERROR! An invalid paramter set was given to the LJ-Cut pair function."
      stop
    end select
  end subroutine
  !=====================================================================

end module
