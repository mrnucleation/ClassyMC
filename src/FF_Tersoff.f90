!================================================================================
module FF_Pair_Tersoff
  use CoordinateTypes
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use VarPrecision

  type, extends(forcefield) :: Pair_Tersoff
    real(dp), allocatable :: epsTable(:,:)
    real(dp), allocatable :: sigTable(:,:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Fc_Func
      procedure, pass :: gik_Func
      procedure, pass :: angleCalc
      procedure, pass :: Constructor => Constructor_Tersoff
      procedure, pass :: DetailedECalc => Detailed_Tersoff
      procedure, pass :: ShiftECalc_Single => Shift_Tersoff_Single
      procedure, pass :: ShiftECalc_Multi => Shift_Tersoff_Multi
      procedure, pass :: SwapInECalc => SwapIn_Tersoff
      procedure, pass :: SwapOutECalc => SwapOut_Tersoff
      procedure, pass :: SetParameter => SetPar_Tersoff
      procedure, pass :: ReadParFile => ReadPar_Tersoff
      procedure, pass :: GetCutOff => GetCutOff_Tersoff
  end type

!====================================================================================== 
  contains
!====================================================================================== 
  pure function Fc_Func(self, r, R_eq, D) result(val)
    use Constants, only: pi
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    real(dp), intent(in) :: r, R_eq, D
    real(dp) :: val  
 
    if( r .lt. (R_eq-D) ) then
      val = 1E0_dp
    elseif( r .lt. (R_eq+D) ) then
      val = 0.5E0_dp * (1E0_dp - sin(pi*(r-R_eq)/(2E0_dp*D)) )
    else
      val = 0E0_dp
    endif

 end function
!======================================================================================      
   pure function gik_Func(self, theta, c, d, h) result(val)
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    real(dp), intent(in) :: theta, c, d, h
    real(dp) :: c_sq, d_sq
    real(dp) :: val  
 
    c_sq = c * c
    d_sq = d * d
    val = 1E0_dp + c_sq/d_sq - c_sq/(d_sq + (cos(theta) - h)**2) 

  end function
!======================================================================================
  pure function angleCalc(self, rx12, ry12, rz12, r12, rx23, ry23, rz23, r23) result(Angle)
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    real(dp), intent(in) :: rx12, ry12, rz12, r12, rx23, ry23, rz23, r23
    real(dp) :: Angle  
             
    Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
    Angle = Angle/(r12*r23)
    if(abs(Angle) .gt. 1E0_dp) then
      Angle = sign(1E0_dp, Angle)
    endif
    Angle = acos(Angle)


  end function
  !=============================================================================+
  subroutine Constructor_Tersoff(self)
    use Common_MolDef, only: nAtomTypes
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    integer :: AllocateStat

!    allocate(self%epsTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rCut = 3E0_dp
    self%rCutSq = 3E0_dp**2

!    write(*,*) 
    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_Tersoff(self, curbox, E_T)
    use ParallelVar, only: nout
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    integer :: iAtom, jAtom, kAtom
    integer :: atmType1, atmType2, atmType3
    real(dp) :: A, B, c, d, R_eq, D2 
    real(dp) :: E_Tersoff
    real(dp) :: lam1, lam2
    real(dp) :: Zeta
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik

    real(dp) :: rMax, rMax_sq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik

    E_T = 0E0_dp
    E_Tersoff = 0E0_dp
    curbox%ETable = 0E0_dp

!    A = tersoffData(atmType1)%A
!    B = tersoffData(atmType1)%B
!    c = tersoffData(atmType1)%c
!    d = tersoffData(atmType1)%d
!    R_eq = tersoffData(atmType1)%R
!    D2 = tersoffData(atmType1)%D2
!    BetaPar = tersoffData(atmType1)%beta
!    n = tersoffData(atmType1)%n
!    h = tersoffData(atmType1)%h
!    lam1 = tersoffData(atmType1)%lam1
!    lam2 = tersoffData(atmType1)%lam2

    rMax = R_eq + D2
    rMax_sq = rMax * rMax

    do iAtom = 1, curbox%nAtoms
      do jAtom = 1, curbox%nAtoms
        if(iAtom .eq. jAtom) then
          cycle
        endif
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij .gt. rMax) then
          cycle
        endif
        rij = sqrt(rij)
        Zeta = 0E0_dp
        do kAtom = 1, curbox%nAtoms
          if( (kAtom == iAtom) .or. (kAtom == jAtom) ) then
            cycle
          endif
          rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
          ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
          rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
          rik = rxij*rxij + ryij*ryij + rzij*rzij
          if(rik .lt. rMax) then
            angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
            Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, R_eq, D2)
          endif
        enddo
        if(Zeta .ne. 0E0_dp) then
          b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
        else
          b1 = 1E0_dp
        endif      
   
        V1 = 0.5E0_dp * self%Fc_Func(rij, R_eq, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
        E_Tersoff = E_Tersoff + V1
!          write(*,*) "V1:", b1, V1
        curbox%ETable(iAtom) = curbox%ETable(iAtom) + V1
        curbox%ETable(jAtom) = curbox%ETable(jAtom) + V1
      enddo
    enddo
!      E_Short = 0.5E0_dp*E_Short 
    write(nout,*) "Tersoff Energy:", E_Tersoff
    E_T = E_Tersoff

   end subroutine
  !=====================================================================
  subroutine Shift_Tersoff_Single(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, iNei, jNei, jAtom, kNei, kAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rMax, rMax_sq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik
    real(dp) :: sub
    real(dp) :: rmin_ij      
    real(dp) :: A, B, c, d, R_eq, D2 
    real(dp) :: E_Tersoff
    real(dp) :: lam1, lam2
    real(dp) :: Zeta, Zeta2
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik
    integer :: nRecalc
    integer :: recalcList(1:60)

    nRecalc = 0
    recalcList = 0

    dispLen = size(disp)
    E_Diff = 0E0_dp
    curbox%dETable = 0E0_dp
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        rxij = curbox % atoms(1, jAtom)  -  disp(iDisp) % x_New
        ryij = curbox % atoms(2, jAtom)  -  disp(iDisp) % y_New
        rzij = curbox % atoms(3, jAtom)  -  disp(iDisp) % z_New
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij .lt. rMax) then
          nRecalc = nRecalc + 1
          recalcList(nRecalc) = jAtom
          rij = sqrt(rij)
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp

          !Compute the Tersoff U_ij component
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if((kAtom .eq. iAtom) .or. (kAtom .eq. jAtom)) then
              cycle
            endif
            rxik = curbox % atoms(1, kAtom)  -  disp(iDisp) % x_New
            ryik = curbox % atoms(2, kAtom)  -  disp(iDisp) % y_New
            rzik = curbox % atoms(3, kAtom)  -  disp(iDisp) % z_New
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik .lt. rMax) then
              rik = sqrt(rik)
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, R_eq, D2)
            endif     
          enddo
          !Compute the Tersoff U_ji component
          do kNei = 1, curbox%NeighList(1)%nNeigh(jAtom)
            kAtom = curbox%NeighList(1)%list(kNei, jAtom)
            if((kAtom .eq. iAtom) .or. (kAtom .eq. jAtom)) then
              cycle
            endif
            rxjk = curbox % atoms(1, kAtom)  -  curbox % atoms(1, jAtom)
            ryjk = curbox % atoms(2, kAtom)  -  curbox % atoms(2, jAtom)
            rzjk = curbox % atoms(3, kAtom)  -  curbox % atoms(3, jAtom)
            rjk = rxij*rxij + ryij*ryij + rzij*rzij
            if(rjk .lt. rMax) then
              rjk = sqrt(rjk)
              angijk = self%angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, R_eq, D2)
            endif     
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif
          V1 = 0.5E0_dp * self%Fc_Func(rij, R_eq, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
          curbox%dETable(jAtom) = curbox%dETable(iAtom) + V1
          E_Diff = E_Diff + V1
        endif
      enddo

      !Computing the contribution from the old poisition. 
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij .lt. rMax) then
          rij = sqrt(rij)
          nRecalc = nRecalc + 1
          recalcList(nRecalc) = jAtom
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp

          !Compute the Tersoff U_ij component
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if((kAtom .eq. iAtom) .or. (kAtom .eq. jAtom)) then
              cycle
            endif
            rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
            ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
            rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik .lt. rMax) then
              rik = sqrt(rik)
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, R_eq, D2)
            endif     
          enddo
          !Compute the Tersoff U_ji component
          do kNei = 1, curbox%NeighList(1)%nNeigh(jAtom)
            kAtom = curbox%NeighList(1)%list(kNei, jAtom)
            if((kAtom .eq. iAtom) .or. (kAtom .eq. jAtom)) then
              cycle
            endif
            rxjk = curbox % atoms(1, kAtom)  -  curbox % atoms(1, jAtom)
            ryjk = curbox % atoms(2, kAtom)  -  curbox % atoms(2, jAtom)
            rzjk = curbox % atoms(3, kAtom)  -  curbox % atoms(3, jAtom)
            rjk = rxij*rxij + ryij*ryij + rzij*rzij
            if(rjk .lt. rMax) then
              rjk = sqrt(rjk)
              angijk = self%angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, R_eq, D2)
            endif     
          enddo
          if(Zeta .ne. 0E0_dp) then
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif
          V1 = 0.5E0_dp * self%Fc_Func(rij, R_eq, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          curbox%dETable(iAtom) = curbox%dETable(iAtom) - V1
          curbox%dETable(jAtom) = curbox%dETable(iAtom) - V1
          E_Diff = E_Diff - V1
        endif
      enddo
    enddo

    !Since the Tersoff is a three body potential, moving a single particle can change the bonded interaction between the first and second
    !neighbor shells.  Therefore these ineteractions must be recomputed. 
    do iNei = 1, nRecalc
      iAtom = recalcList(iNei)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        if(iAtom .eq. jAtom) then
          cycle
        endif
        if( any(disp(:)%atmIndx .eq. jAtom) ) then
          cycle
        endif
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij < rMax) then
          rij = sqrt(rij)
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(jNei, iAtom)
            if( (kAtom == iAtom) .or. (kAtom == jAtom) ) then
              cycle
            endif
            if(kAtom == disp(1)%atmIndx) then
              rxik = disp(1)%x_new  -  curbox % atoms(1, iAtom)
              ryik = disp(1)%y_new  -  curbox % atoms(2, iAtom)
              rzik = disp(1)%z_new  -  curbox % atoms(3, iAtom)
              rik = rxij*rxij + ryij*ryij + rzij*rzij
              if(rik .lt. rMax) then
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, R_eq, D2)
                Zeta = Zeta + sub
              endif
              rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
              ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
              rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
              rik = rxij*rxij + ryij*ryij + rzij*rzij
              if(rik .lt. rMax) then
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, R_eq, D2)
                Zeta2 = Zeta2 + sub
              endif
            else
              rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
              ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
              rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
              rik = rxij*rxij + ryij*ryij + rzij*rzij
              if(rik .lt. rMax) then
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, R_eq, D2)
                Zeta = Zeta + sub
                Zeta2 = Zeta2 + sub
              endif
            endif
          enddo
        endif
        if(Zeta .ne. 0E0_dp) then
          b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
        else
          b1 = 1E0_dp
        endif
        if(Zeta2 .ne. 0E0_dp) then
          b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
        else
          b2 = 1E0_dp
        endif    
   
        V1 = 0.5E0_dp * self%Fc_Func(rij, R_eq, D2) * (B*exp(-lam2*rij))*(b2 - b1)
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
        curbox%dETable(jAtom) = curbox%dETable(iAtom) + V1
        E_Diff = E_Diff + V1
      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine Shift_Tersoff_Multi(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_Tersoff), intent(in) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
   
  end subroutine
  !=====================================================================
  subroutine SwapIn_Tersoff(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_Tersoff), intent(in) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inOut) :: E_Diff
      integer :: iDisp, iAtom, jAtom, dispLen
      integer :: atmType1, atmType2
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      dispLen = size(disp)
      E_Diff = 0E0
      curbox%dETable = 0E0


      do iDisp = 1, dispLen
        iAtom = disp(iDisp)%atmIndx
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, curbox % nAtoms        
          atmType2 = curbox % AtomType(jAtom)
          ep = self%epsTable(atmType1,atmType2)
          sig_sq = self%sigTable(atmType1,atmType2)          
          rmin_ij = self%rMinTable(atmType1,atmType2)          

          rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
          ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
          rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          if(rsq < self%rCutSq) then
            LJ = (sig_sq/rsq)
            LJ = LJ * LJ * LJ
            LJ = ep * LJ * (LJ-1E0)              
            E_Diff = E_Diff - LJ
            curbox % dETable(iAtom) = curbox % dETable(iAtom) - LJ
            curbox % dETable(jAtom) = curbox % dETable(jAtom) - LJ
          endif
        enddo
      enddo
  end subroutine
  !=====================================================================
  subroutine SwapOut_Tersoff(self, curbox, atmIndx, E_Diff)
    implicit none
      class(Pair_Tersoff), intent(in) :: self
      class(SimBox), intent(inout) :: curbox
      real(dp), intent(inOut) :: E_Diff
      integer, intent(in) :: atmIndx(:)
      integer :: iIndx, iAtom, jAtom, remLen
      integer :: atmType1, atmType2
      real(dp) :: rx, ry, rz, rsq
      real(dp) :: ep, sig_sq
      real(dp) :: LJ
      real(dp) :: rmin_ij      

      remLen = size(atmIndx)
      E_Diff = 0E0
      curbox%dETable = 0E0


      do iIndx = 1, remLen
        iAtom = atmIndx(iIndx)
        atmType1 = curbox % AtomType(iAtom)
        do jAtom = 1, curbox % nAtoms        
          atmType1 = curbox % AtomType(jAtom)
          ep = self % epsTable(atmType1,atmType2)
          sig_sq = self % sigTable(atmType1,atmType2)          
          rmin_ij = self % rMinTable(atmType1,atmType2)          

          rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
          ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
          rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
          rsq = rx*rx + ry*ry + rz*rz
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + LJ
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + LJ
        enddo
      enddo
  end subroutine
  !=====================================================================
  subroutine SetPar_Tersoff(self, parIndex,  parVal)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    integer, intent(in) :: parIndex(:)
    real(dp), intent(in) :: parVal

    select case( parIndex(1) )
    case(1) !Epsilon
      self%epsTable(parIndex(2), parIndex(3)) = parVal
    case(2) !Sigma
      self%sigTable(parIndex(2), parIndex(3)) = parVal
    case(3) !rMin
      self%rMinTable(parIndex(2), parIndex(3)) = parVal
    case(4) !rCut
      self%rCut = parVal
      self%rCutSq = parVal * parVal
    case default
      write(*,*) "ERROR! An invalid paramter set was given to the LJ-Cut pair function."
      stop
    end select
  end subroutine

  !=====================================================================
  subroutine ReadPar_Tersoff(self, fileName)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    character(len=*), intent(in) :: fileName
    write(*,*) "Tersoff SAYING HELLO!!!"
  end subroutine
  !=============================================================================+
    function GetCutOff_Tersoff(self) result(rCut)
      implicit none
      class(Pair_Tersoff), intent(inout) :: self
      real(dp) :: rCut

      rCut = self%rCut
    end function
  !=====================================================================

end module