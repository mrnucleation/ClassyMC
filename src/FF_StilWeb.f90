!================================================================================
module FF_Pair_StilWeb
  use CoordinateTypes
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use VarPrecision

  type :: SW2Body
    real(dp) :: A, B
  end type

  type :: SW3Body
    real(dp) :: h, lam3, c, d, gam
  end type


  type, extends(forcefield) :: Pair_StilWeb
    type(SW2Body), allocatable :: stilwebPair(:,:)
    type(SW3Body), allocatable :: stilwebAngle(:,:,:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: angleCalc
      procedure, pass :: Constructor => Constructor_StilWeb
      procedure, pass :: DetailedECalc => Detailed_StilWeb
      procedure, pass :: ShiftECalc_Single => Shift_StilWeb_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_StilWeb_Multi
      procedure, pass :: NewECalc => New_StilWeb
!      procedure, pass :: OldECalc => Old_StilWeb
      procedure, pass :: ProcessIO => ProcessIO_StilWeb
      procedure, pass :: GetCutOff => GetCutOff_StilWeb
  end type

!================================================================================== 
  contains
!======================================================================================
   pure function angleCalc(self, rx12, ry12, rz12, r12, rx23, ry23, rz23, r23) result(Angle)
    implicit none
    class(Pair_StilWeb), intent(in) :: self
    real(dp), intent(in) :: rx12, ry12, rz12, r12, rx23, ry23, rz23, r23
    real(dp) :: Angle  
             
    Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
    Angle = Angle/(r12*r23)
!    if(abs(Angle) > 1E0_dp) then
!      Angle = sign(1E0_dp, Angle)
!    endif
!    Angle = acos(Angle)


  end function
  !=============================================================================+
  subroutine Constructor_StilWeb(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_StilWeb), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%stilwebPair(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%stilwebAngle(1:nAtomTypes,1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMinTable = 0.5E0_dp
    self%rCut = 3E0_dp
    self%rCutSq = 3E0_dp**2

    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_StilWeb(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    implicit none
    class(Pair_StilWeb), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iAtom, jAtom, kAtom
    integer :: atmType1, atmType2, atmType3
    real(dp) :: A, B, c, d, Reqij, Reqik, Dij, Dik
    real(dp) :: E_StilWeb
    real(dp) :: lam1, lam2
    real(dp) :: Zeta
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1
    real(dp) :: angijk, angjik

    real(dp) :: rMax, rMax_sq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik

    E_T = 0E0_dp
    E_StilWeb = 0E0_dp
    curbox%ETable = 0E0_dp

    write(nout,*) "StilWeb Energy:", E_StilWeb
    E_T = E_StilWeb

   end subroutine
  !=====================================================================
  subroutine Shift_StilWeb_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_StilWeb), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    integer :: iDisp, iAtom, iNei, jNei, jAtom, kNei, kAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2, atmType3
    real(dp) :: rMaxSq, rMinSq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik
    real(dp) :: sub
    real(dp) :: rmin_ij
    real(dp) :: A, B, c, d, Req, D2 
    real(dp) :: lam1, lam2
    real(dp) :: Zeta, Zeta2
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik
    real(dp) :: E_StilWeb
    integer :: nRecalc
    integer :: recalcList(1:200)

     ! The recalcList is a list of particles whose intermolecular interactions have
     ! changed, but the particles themselves did not move.
    nRecalc = 0
    recalcList = 0

    dispLen = size(disp)
    E_Diff = 0E0_dp
    curbox%dETable = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        atmType2 = curbox % AtomType(jAtom)

        rxij = curbox % atoms(1, jAtom)  -  disp(iDisp) % x_New
        ryij = curbox % atoms(2, jAtom)  -  disp(iDisp) % y_New
        rzij = curbox % atoms(3, jAtom)  -  disp(iDisp) % z_New
        call curbox%Boundary(rxij, ryij, rzij)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        rMaxSq = self%stilwebPair(atmType1, atmType2) % rMaxSq

        if(rij < rMaxSq) then
          rMinSq = self % rMinTable(atmType1, atmType2)          
          if(rij < rMinSq) then
            accept = .false.
            return
          endif
          if(all(recalcList /= jAtom) ) then
            nRecalc = nRecalc + 1
            recalcList(nRecalc) = jAtom
          endif
          rij = sqrt(rij)
          Zeta = 0E0_dp
          Zeta2 = 0E0_dp

          !Compute the StilWeb U_ij component
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if((kAtom == iAtom) .or. (kAtom == jAtom)) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%stilwebPair(atmType1, atmType3) % rMaxSq
            rxik = curbox % atoms(1, kAtom)  -  disp(iDisp) % x_New
            ryik = curbox % atoms(2, kAtom)  -  disp(iDisp) % y_New
            rzik = curbox % atoms(3, kAtom)  -  disp(iDisp) % z_New
            call curbox%Boundary(rxik, ryik, rzik)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              rik = sqrt(rik)
              D2 = self%stilwebPair(atmType1, atmType3) % D
              Req = self%stilwebPair(atmType1, atmType3) % Req
              c = self%stilwebAngle(atmType1, atmType2, atmType3) % c
              d = self%stilwebAngle(atmType1, atmType2, atmType3) % d
              h = self%stilwebAngle(atmType1, atmType2, atmType3) % h            
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
            endif     
          enddo
          !Compute the StilWeb U_ji component
          do kNei = 1, curbox%NeighList(1)%nNeigh(jAtom)
            kAtom = curbox%NeighList(1)%list(kNei, jAtom)
            if((kAtom == iAtom) .or. (kAtom == jAtom)) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%stilwebPair(atmType2, atmType3) % rMaxSq

            rxjk = curbox % atoms(1, kAtom)  -  curbox % atoms(1, jAtom)
            ryjk = curbox % atoms(2, kAtom)  -  curbox % atoms(2, jAtom)
            rzjk = curbox % atoms(3, kAtom)  -  curbox % atoms(3, jAtom)
            call curbox%Boundary(rxjk, ryjk, rzjk)
            rjk = rxjk*rxjk + ryjk*ryjk + rzjk*rzjk
            if(rjk < rMaxSq) then
              D2 = self%stilwebPair(atmType2, atmType3) % D
              Req = self%stilwebPair(atmType2, atmType3) % Req
              c = self%stilwebAngle(atmType2, atmType1, atmType3) % c
              d = self%stilwebAngle(atmType2, atmType1, atmType3) % d
              h = self%stilwebAngle(atmType2, atmType1, atmType3) % h            
              rjk = sqrt(rjk)
              angijk = self%angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta2 = Zeta2 + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rjk, Req, D2)
            endif     
          enddo
          if(Zeta .ne. 0E0_dp) then
            BetaPar = self%stilwebPair(atmType1, atmType2)%beta
            n =  self%stilwebPair(atmType1, atmType2)%n
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 /= 0E0_dp) then
            BetaPar = self%stilwebPair(atmType2, atmType1)%beta
            n =  self%stilwebPair(atmType2, atmType1)%n
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif

          if(self%symetric) then
            A = self%stilwebPair(atmType1, atmType2) % A
            B = self%stilwebPair(atmType1, atmType2) % B
            lam1 = self%stilwebPair(atmType1, atmType2) % lam1
            lam2 = self%stilwebPair(atmType1, atmType2) % lam2
            Req = self%stilwebPair(atmType1, atmType2) % REq
            D2 = self%stilwebPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          else
            A = self%stilwebPair(atmType1, atmType2) % A
            B = self%stilwebPair(atmType1, atmType2) % B
            lam1 = self%stilwebPair(atmType1, atmType2) % lam1
            lam2 = self%stilwebPair(atmType1, atmType2) % lam2
            Req = self%stilwebPair(atmType1, atmType2) % REq
            D2 = self%stilwebPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 

            A = self%stilwebPair(atmType2, atmType1) % A
            B = self%stilwebPair(atmType2, atmType1) % B
            lam1 = self%stilwebPair(atmType2, atmType1) % lam1
            lam2 = self%stilwebPair(atmType2, atmType1) % lam2
            Req = self%stilwebPair(atmType2, atmType1) % REq
            D2 = self%stilwebPair(atmType2, atmType1) % D
            V1 = V1 + 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b2*B*exp(-lam2*rij)) 
          
          endif
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
          curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
          E_Diff = E_Diff + V1

        endif
      enddo

      !Computing the contribution from the old poisition. 
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        atmType2 = curbox % AtomType(jAtom)
        rMaxSq = self%stilwebPair(atmType1, atmType2) % rMaxSq
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxij, ryij, rzij)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij < rMaxSq) then
          rij = sqrt(rij)
          if(all(recalcList /= jAtom)) then
            nRecalc = nRecalc + 1
            recalcList(nRecalc) = jAtom
          endif
          Zeta = 0E0_dp !Zeta ij
          Zeta2 = 0E0_dp !Zeta ji

          !Compute the Old StilWeb U_ij component
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if((kAtom == iAtom) .or. (kAtom == jAtom)) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
            ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
            rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
            call curbox%Boundary(rxik, ryik, rzik)
            rMaxSq = self%stilwebPair(atmType1, atmType3) % rMaxSq
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              rik = sqrt(rik)
              D2 = self%stilwebPair(atmType1, atmType3) % D
              Req = self%stilwebPair(atmType1, atmType3) % Req
              c = self%stilwebAngle(atmType1, atmType2, atmType3) % c
              d = self%stilwebAngle(atmType1, atmType2, atmType3) % d
              h = self%stilwebAngle(atmType1, atmType2, atmType3) % h            
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
            endif     
          enddo
          !Compute the Old StilWeb U_ji component
          do kNei = 1, curbox%NeighList(1)%nNeigh(jAtom)
            kAtom = curbox%NeighList(1)%list(kNei, jAtom)
            if((kAtom == iAtom) .or. (kAtom == jAtom)) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rxjk = curbox % atoms(1, kAtom)  -  curbox % atoms(1, jAtom)
            ryjk = curbox % atoms(2, kAtom)  -  curbox % atoms(2, jAtom)
            rzjk = curbox % atoms(3, kAtom)  -  curbox % atoms(3, jAtom)
            call curbox%Boundary(rxjk, ryjk, rzjk)
            rMaxSq = self%stilwebPair(atmType2, atmType3) % rMaxSq
            rjk = rxjk*rxjk + ryjk*ryjk + rzjk*rzjk
            if(rjk < rMaxSq) then
              rjk = sqrt(rjk)
              D2 = self%stilwebPair(atmType2, atmType3) % D
              Req = self%stilwebPair(atmType2, atmType3) % Req
              c = self%stilwebAngle(atmType2, atmType1, atmType3) % c
              d = self%stilwebAngle(atmType2, atmType1, atmType3) % d
              h = self%stilwebAngle(atmType2, atmType1, atmType3) % h            

              angijk = self%angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta2 = Zeta2 + self%gik_Func(angijk, c, d, h) * self%Fc_Func(rjk, Req, D2)
            endif     
          enddo
          if(Zeta .ne. 0E0_dp) then
            BetaPar = self%stilwebPair(atmType1, atmType2) % beta
            n = self%stilwebPair(atmType1, atmType2) % n
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b1 = 1E0_dp
          endif
          if(Zeta2 .ne. 0E0_dp) then
            BetaPar = self%stilwebPair(atmType2, atmType1) % beta
            n = self%stilwebPair(atmType2, atmType1) % n
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
            b2 = 1E0_dp
          endif

          if(self%symetric) then
            A = self%stilwebPair(atmType1, atmType2) % A
            B = self%stilwebPair(atmType1, atmType2) % B
            lam1 = self%stilwebPair(atmType1, atmType2) % lam1
            lam2 = self%stilwebPair(atmType1, atmType2) % lam2
            Req = self%stilwebPair(atmType1, atmType2) % REq
            D2 = self%stilwebPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          else
            A = self%stilwebPair(atmType1, atmType2) % A
            B = self%stilwebPair(atmType1, atmType2) % B
            lam1 = self%stilwebPair(atmType1, atmType2) % lam1
            lam2 = self%stilwebPair(atmType1, atmType2) % lam2
            Req = self%stilwebPair(atmType1, atmType2) % REq
            D2 = self%stilwebPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 

            A = self%stilwebPair(atmType2, atmType1) % A
            B = self%stilwebPair(atmType2, atmType1) % B
            lam1 = self%stilwebPair(atmType2, atmType1) % lam1
            lam2 = self%stilwebPair(atmType2, atmType1) % lam2
            Req = self%stilwebPair(atmType2, atmType1) % REq
            D2 = self%stilwebPair(atmType2, atmType1) % D
            V1 = V1 + 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b2*B*exp(-lam2*rij)) 
          
          endif
 
!          write(*,*) "Old:", iAtom, jAtom, -V1
          curbox%dETable(iAtom) = curbox%dETable(iAtom) - V1
          curbox%dETable(jAtom) = curbox%dETable(iAtom) - V1
          E_Diff = E_Diff - V1
        endif
      enddo
    enddo

    !Since the StilWeb is a three body potential, moving a single particle can change the bonded 
    !interaction of particles that did not move.  These interactions must also be recomputed.
    do iNei = 1, nRecalc
      iAtom = recalcList(iNei)
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        if(iAtom == jAtom) then
          cycle
        endif
        if( any(disp(:)%atmIndx == jAtom) ) then
          cycle
        endif
        atmType2 = curbox % AtomType(jAtom)
        rMaxSq = self%stilwebPair(atmType1, atmType2) % rMaxSq
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxij, ryij, rzij)

        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij < rMaxSq) then
          rij = sqrt(rij)
          Zeta = 0E0_dp  !Zeta is the angle term for the new config
          Zeta2 = 0E0_dp !Zeta is the angle term for the old config
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if( (kAtom == iAtom) .or. (kAtom == jAtom) ) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%stilwebPair(atmType1, atmType3) % rMaxSq
             ! If kAtom is the particle that moved old and new must
             ! be recomputed separately. Otherwise, old and new are the same.
            if(kAtom == disp(1)%atmIndx) then
              rxik = disp(1)%x_new  -  curbox % atoms(1, iAtom)
              ryik = disp(1)%y_new  -  curbox % atoms(2, iAtom)
              rzik = disp(1)%z_new  -  curbox % atoms(3, iAtom)
              call curbox%Boundary(rxik, ryik, rzik)
              rik = rxik*rxik + ryik*ryik + rzik*rzik
              if(rik < rMaxSq) then
                D2 = self%stilwebPair(atmType1, atmType3) % D
                Req = self%stilwebPair(atmType1, atmType3) % Req
                c = self%stilwebAngle(atmType1, atmType2, atmType3) % c
                d = self%stilwebAngle(atmType1, atmType2, atmType3) % d
                h = self%stilwebAngle(atmType1, atmType2, atmType3) % h            
                rik = sqrt(rik)
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
                Zeta = Zeta + sub
              endif
              rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
              ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
              rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
              call curbox%Boundary(rxik, ryik, rzik)
              rik = rxik*rxik + ryik*ryik + rzik*rzik
              if(rik < rMaxSq) then
                D2 = self%stilwebPair(atmType1, atmType3) % D
                Req = self%stilwebPair(atmType1, atmType3) % Req
                c = self%stilwebAngle(atmType1, atmType2, atmType3) % c
                d = self%stilwebAngle(atmType1, atmType2, atmType3) % d
                h = self%stilwebAngle(atmType1, atmType2, atmType3) % h            
                rik = sqrt(rik)
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
                Zeta2 = Zeta2 + sub
              endif
            else
              rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
              ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
              rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
              call curbox%Boundary(rxik, ryik, rzik)
              rik = rxik*rxik + ryik*ryik + rzik*rzik
              if(rik < rMaxSq) then
                D2 = self%stilwebPair(atmType1, atmType3) % D
                Req = self%stilwebPair(atmType1, atmType3) % Req
                c = self%stilwebAngle(atmType1, atmType2, atmType3) % c
                d = self%stilwebAngle(atmType1, atmType2, atmType3) % d
                h = self%stilwebAngle(atmType1, atmType2, atmType3) % h            
                rik = sqrt(rik)

                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
                Zeta = Zeta + sub
                Zeta2 = Zeta2 + sub
              endif
            endif
          enddo
        endif

        BetaPar = self%stilwebPair(atmType1, atmType2) % beta
        n = self%stilwebPair(atmType1, atmType2) % n
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
   
!        A = self%stilwebPair(atmType1, atmType2) % A
        B = self%stilwebPair(atmType1, atmType2) % B
!        lam1 = self%stilwebPair(atmType1, atmType2) % lam1
        lam2 = self%stilwebPair(atmType1, atmType2) % lam2
        D2 = self%stilwebPair(atmType1, atmType2) % D
        Req = self%stilwebPair(atmType1, atmType2) % Req
        V1 = 0.5E0_dp*self%Fc_Func(rij, Req, D2) * (B*exp(-lam2*rij))*(b2 - b1)
!        write(*,*) "Recalc", iAtom, jAtom, V1
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
        E_Diff = E_Diff + V1
      enddo
    enddo
 
!    write(*,*) E_Diff
  end subroutine  
!===================================================================== 
  subroutine New_StilWeb(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_StilWeb), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, iNei, jNei, jAtom, kNei, kAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2, atmType3
    real(dp) :: rMaxSq, rMinSq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik
    real(dp) :: sub
    real(dp) :: rmin_ij
    real(dp) :: A, B, c, d, Req, D2 
    real(dp) :: lam1, lam2
    real(dp) :: Zeta, Zeta2
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik
    real(dp) :: E_StilWeb
    integer :: nRecalc
    integer :: recalcList(1:200)

  end subroutine
!=====================================================================
  subroutine ProcessIO_StilWeb(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetAllCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: outEngUnit, outLenUnit
    implicit none
    class(Pair_StilWeb), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

    character(len=30) :: command
    logical :: logicVal
    integer :: iPar, lineStat
    integer :: type1, type2, type3
    real(dp) :: parList(1:8)
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("symetric")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) logicVal
        self%symetric = logicVal

      case("pair")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        do iPar = 1, 8
          call GetXCommand(line, command, 3+iPar, lineStat)
          read(command, *) parList(iPar)
        enddo
        self%stilwebPair(type1, type2)%A = parList(1) * outEngUnit
        self%stilwebPair(type1, type2)%B = parList(2) * outEngUnit
        self%stilwebPair(type1, type2)%lam1 = parList(3) * outLenUnit
        self%stilwebPair(type1, type2)%lam2 = parList(4) * outLenUnit
        self%stilwebPair(type1, type2)%Req = parList(5) * outLenUnit
        self%stilwebPair(type1, type2)%D = parList(6) * outLenUnit
        self%stilwebPair(type1, type2)%beta = parList(7)
        self%stilwebPair(type1, type2)%n = parList(8)
        self%stilwebPair(type1, type2)%rMax = (parList(5) + parList(6))*outLenUnit
        self%stilwebPair(type1, type2)%rMaxSq = ((parList(5) + parList(6))*outLenUnit)**2
        if(self%symetric) then
          self%stilwebPair(type2, type1)%A = parList(1) * outEngUnit
          self%stilwebPair(type2, type1)%B = parList(2) * outEngUnit
          self%stilwebPair(type2, type1)%lam1 = parList(3) * outLenUnit
          self%stilwebPair(type2, type1)%lam2 = parList(4) * outLenUnit
          self%stilwebPair(type2, type1)%Req = parList(5) * outLenUnit
          self%stilwebPair(type2, type1)%D = parList(6) * outLenUnit
          self%stilwebPair(type2, type1)%beta = parList(7)
          self%stilwebPair(type2, type1)%n = parList(8)
          self%stilwebPair(type2, type1)%rMax = (parList(5) + parList(6))*outLenUnit
          self%stilwebPair(type2, type1)%rMaxSq = ((parList(5) + parList(6))*outLenUnit)**2
        endif

      case("angle")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        call GetXCommand(line, command, 4, lineStat)
        read(command, *) type3
        do iPar = 1, 5
          call GetXCommand(line, command, 4+iPar, lineStat)
          read(command, *) parList(iPar)
        enddo
        self%stilwebAngle(type1, type2, type3)%h = parList(1)
        self%stilwebAngle(type1, type2, type3)%lam3 = parList(2)* outLenUnit
        self%stilwebAngle(type1, type2, type3)%c = parList(3)
        self%stilwebAngle(type1, type2, type3)%d = parList(4)
        self%stilwebAngle(type1, type2, type3)%gam = parList(5)
        if(self%symetric) then
          self%stilwebAngle(type2, type1, type3)%h = parList(1)
          self%stilwebAngle(type2, type1, type3)%lam3 = parList(2)* outLenUnit
          self%stilwebAngle(type2, type1, type3)%c = parList(3)
          self%stilwebAngle(type2, type1, type3)%d = parList(4)
          self%stilwebAngle(type2, type1, type3)%gam = parList(5)
        endif

      case default
        linestat = -1
    end select


  end subroutine

 !=============================================================================+
    function GetCutOff_StilWeb(self) result(rCut)
      use Common_MolInfo, only: nAtomTypes
      implicit none
      class(Pair_StilWeb), intent(inout) :: self
      real(dp) :: rCut

      integer :: type1, type2
      real(dp) :: Req, D, rMax

      rCut = -1E0_dp
      do type1 = 1, nAtomTypes
        do type2 = 1, nAtomTypes
          rMax = self%stilwebPair(type1, type2) % RMax
          if(rCut < rMax) then
            rCut = rMax
          endif
        enddo
      enddo

      rCut = self%rCut
    end function
  !=====================================================================

end module
