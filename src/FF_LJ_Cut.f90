module FF_Pair_LJ_Cut
  use ForceFieldTemplate
  use VarPrecision
  use SimBoxDef, only: SimBox
  type, extends(forcefield) :: Pair_LJ_Cut
    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: rMinTable(:,:)
    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: DetailedECalc => Detailed_LJ_Cut
      procedure, pass :: ShiftECalc_Single => Shift_LJ_Cut_Single
      procedure, pass :: ShiftECalc_Multi => Shift_LJ_Cut_Multi
      procedure, pass :: SwapInECalc => SwapIn_LJ_Cut
      procedure, pass :: SwapOutECalc => SwapOut_LJ_Cut
      procedure, pass :: SetParameter => SetPar_LJ_Cut
  end type

  contains
  !===================================================================================
  subroutine Detailed_LJ_Cut(self, curbox)
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
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
          rmin_ij = self%r_min_tab(atmType1,atmType2)          

          rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
          rsq = rx**2 + ry**2 + rz**2
          if(rsq < self%rCutSq) then
            if(rsq < rmin_ij) then
              stop "ERROR! Overlaping atoms found in the current configuration!"
            endif 
            LJ = (sig_sq/rsq)**3
            LJ = ep * LJ * (LJ-1E0)              
            E_LJ = E_LJ + LJ
            ETable(iAtom) = ETable(iAtom) + LJ
            ETable(jAtom) = ETable(jAtom) + LJ 
          endif
        enddo
      enddo

      
      write(nout,*) "Lennard-Jones Energy:", E_LJ
      
      E_T = E_LJ    
   end subroutine
  !=====================================================================
  subroutine Shift_LJ_Cut_Single(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
      type(simBox), intent(in) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff
      integer :: iDisp, iAtom, jAtom, dispLen
      integer :: maxIndx, minIndx
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      dispLen = len(disp)
      E_Diff = 0E0
      dETable = 0E0
      maxIndx = 0
      minIndx = curbox % nAtoms
      do iDisp = 1, dispLen
        if(disp(iDisp)% atmindx > maxIndx) then
          maxIndx = disp(iDisp)% atmindx
        endif
        if(disp(iDisp)% atmindx < minIndx) then
          minIndx = disp(iDisp)% atmindx
        endif
      enddo      

      if(minIndx .eq. 1) then
        minIndx = 0
      endif
      if(maxIndx .eq. curbox % nAtoms) then
        maxIndx = curbox % nAtoms + 1
      endif


      do iDisp = 1, dispLen
        iAtom = disp(iDisp)%atmIndx
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, minIndx        
          atmType1 = curbox % AtomType(jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = r_min_tab(atmType1,atmType2)          

          rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
          ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
          rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < rCutSq) then
            if(rsq < rmin_ij) then
            endif 
            LJ = (sig_sq/rsq)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)              
            E_Diff = E_Diff + LJ
            dETable(iAtom) = dETable(iAtom) + LJ
            dETable(jAtom) = dETable(jAtom) + LJ
          endif

          rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < self%rCutSq) then
            if(rsq < rmin_ij) then
            endif 
            LJ = (sig_sq/rsq)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)              
            E_Diff = E_Diff - LJ
            dETable(iAtom) = dETable(iAtom) - LJ
            dETable(jAtom) = dETable(jAtom) - LJ
          endif
        enddo

        do jAtom = maxIndx, curbox % nAtoms
          atmType1 = curbox % AtomType(jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = r_min_tab(atmType1,atmType2)          

          rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
          ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
          rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
          r = rx*rx + ry*ry + rz*rz
          LJ = (sig_sq/r)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0)              
          E_Diff = E_Diff + LJ
          dETable(iAtom) = dETable(iAtom) + LJ
          dETable(jAtom) = dETable(jAtom) + LJ
!          if(r .lt. rmin_ij) then
!          endif 

          rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
          r = rx*rx + ry*ry + rz*rz
          LJ = (sig_sq/r)
          LJ = LJ * LJ * LJ
          LJ = ep * LJ * (LJ-1E0)              
          E_Diff = E_Diff - LJ
          dETable(iAtom) = dETable(iAtom) - LJ
          dETable(jAtom) = dETable(jAtom) - LJ
        enddo

      enddo
 
  end subroutine
  !=====================================================================
  subroutine SwapIn_LJ_Cut(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
      type(simBox), intent(in) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff
      integer :: iDisp, iAtom, jAtom, dispLen
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      dispLen = len(disp)
      E_Diff = 0E0
      dETable = 0E0


      do iDisp = 1, dispLen
        iAtom = disp(iDisp)%atmIndx
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, curbox % nAtoms        
          atmType1 = curbox % AtomType(jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = r_min_tab(atmType1,atmType2)          

          rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
          ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
          rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < self%rCutSq) then
            LJ = (sig_sq/r)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)              
            E_Diff = E_Diff - LJ
            dETable(iAtom) = dETable(iAtom) - LJ
            dETable(jAtom) = dETable(jAtom) - LJ
          endif
        enddo
      enddo
  end subroutine
  !=====================================================================
  subroutine SwapOut_LJ_Cut(self)
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
      type(simBox), intent(in) :: curbox
      real(dp), intent(inOut) :: E_Diff
      integer, intent(in) :: atmIndx(:)
      integer :: iIndx, iAtom, jAtom, remLen
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      remLen = len(iRemove)
      E_Diff = 0E0
      dETable = 0E0


      do iIndx = 1, remLen
        iAtom = atmIndx(iIndx)
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, curbox % nAtoms        
          atmType1 = curbox % AtomType(jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = r_min_tab(atmType1,atmType2)          

          rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          dETable(iAtom) = dETable(iAtom) + LJ
          dETable(jAtom) = dETable(jAtom) + LJ
        enddo
      enddo
  end subroutine
  !=====================================================================
  subroutine SetPar_LJ_Cut(self, parIndex,  parVal)
    implicit none
    class(Pair_LJ_Cut), intent(in) :: self
    integer, intent(in) :: parIndex(:)
    real(dp), intent(in) :: parVal

    select case( parIndex(1) )
    case(1) !Epsilon
      self%epsTable(parIndex(2), parIndex(3)) = parVal
    case(2) !Sigma
      self%sigTable(parIndex(2), parIndex(3)) = parVal
    case(3) !rMin
      self%rMinTable(parIndex(2), parIndex(3)) = parVal
    case(3) !rCut
      self%rCut = parVal
      self%rCutSq = parVal * parVal
    case default
      write(*,*) "ERROR! An invalid paramter set was given to the LJ-Cut pair function."
      stop
    end select
  end subroutine
  !=====================================================================

end module
