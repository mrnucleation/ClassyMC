!================================================================================
module FF_Pair_Tersoff
  use CoordinateTypes
  use Template_ForceField, only: ForceField
  use Template_SimBox, only: SimBox
  use VarPrecision

  type :: Tersoff2Body
    real(dp) :: A, B, lam1, lam2, Req, D, beta, n, rMax, rMaxSq
  end type

  type :: Tersoff3Body
    real(dp) :: h, lam3, c, d, gam
  end type


  type, extends(forcefield) :: Pair_Tersoff
    logical :: symetric = .true.
    type(Tersoff2Body), allocatable :: tersoffPair(:,:)
    type(Tersoff3Body), allocatable :: tersoffAngle(:,:,:)
    real(dp), allocatable :: rMinTable(:,:)
    real(dp), allocatable :: rMin(:)
    real(dp) :: dimer = 1E0_dp
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Fc_Func
      procedure, pass :: gik_Func
      procedure, pass :: angleCalc
      procedure, pass :: Constructor => Constructor_Tersoff
      procedure, pass :: DetailedECalc => Detailed_Tersoff
      procedure, pass :: DiffECalc => DiffECalc_Tersoff
      procedure, pass :: ShiftECalc_Single => Shift_Tersoff_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_Tersoff_Multi
      procedure, pass :: NewECalc => New_Tersoff
      procedure, pass :: OldECalc => Old_Tersoff
      procedure, pass :: OrthoVolECalc => OrthoVol_Tersoff
      procedure, pass :: ProcessIO => ProcessIO_Tersoff
      procedure, pass :: GetCutOff => GetCutOff_Tersoff
  end type

!================================================================================== 
  contains
!==================================================================================
  pure function Fc_Func(self, r, Req, D) result(val)
    use ClassyConstants, only: pi
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    real(dp), intent(in) :: r, Req, D
    real(dp) :: val  
 
    if( r < (Req-D) ) then
      val = 1E0_dp
    elseif( r < (Req+D) ) then
      val = 0.5E0_dp * (1E0_dp - sin(pi*(r-Req)/(2E0_dp*D)) )
    else
      val = 0E0_dp
    endif

 end function
!===============================================================================      
  pure function gik_Func(self, theta, c, d, h) result(val)
    implicit none
    class(Pair_Tersoff), intent(in) :: self
    real(dp), intent(in) :: theta, c, d, h
    real(dp) :: c_sq, d_sq
    real(dp) :: val, t_term
 
    c_sq = c * c
    d_sq = d * d
!    val = 1E0_dp + c_sq/d_sq - c_sq/(d_sq + (cos(theta) - h)**2) 
    t_term = (theta-h)
    val = 1E0_dp + c_sq/d_sq - c_sq/(d_sq + t_term*t_term) 

  end function
!======================================================================================
   pure function angleCalc(self, rx12, ry12, rz12, r12, rx23, ry23, rz23, r23) result(Angle)
    implicit none
    class(Pair_Tersoff), intent(in) :: self
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
  subroutine Constructor_Tersoff(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%tersoffPair(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)
    allocate(self%tersoffAngle(1:nAtomTypes,1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMinTable = 0.5E0_dp
    self%rCut = 3E0_dp
    self%rCutSq = 3E0_dp**2

    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_Tersoff(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iAtom, jAtom, kAtom
    integer :: atmType1, atmType2, atmType3
    real(dp) :: A, B, c, d, Reqij, Reqik, Dij, Dik
    real(dp) :: E_Tersoff
    real(dp) :: lam1, lam2, lam3, gam
    real(dp) :: Zeta
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1
    real(dp) :: angijk, angjik

    real(dp) :: rMax, rMax_sq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik

    E_T = 0E0_dp
    E_Tersoff = 0E0_dp
    curbox%ETable = 0E0_dp

    accept = .true.
    do iAtom = 1, curbox%nMaxAtoms
      atmType1 = curbox % AtomType(iAtom)
      if( .not. curbox%IsActive(iAtom) ) then
        cycle
      endif

      do jAtom = 1, curbox%nMaxAtoms
        if(iAtom == jAtom) then
          cycle
        endif
        atmType2 = curbox % AtomType(jAtom)
        if( .not. curbox%IsActive(jAtom) ) then        
          cycle
        endif
        Reqij = self%tersoffPair(atmType2, atmType1) % REq
        Dij = self%tersoffPair(atmType2, atmType1) % D
        rMax = Dij + Reqij
        rMax_sq = rMax * rMax

        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxij, ryij, rzij)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij > rMax_Sq) then
          cycle
        endif
        rij = sqrt(rij)
        Zeta = 0E0_dp
        do kAtom = 1, curbox%nMaxAtoms
          if( (kAtom == iAtom) .or. (kAtom == jAtom) ) then
            cycle
          endif
          if( .not. curbox%IsActive(kAtom) ) then        
            cycle
          endif
          atmType3 = curbox % AtomType(kAtom)
          Reqik = self%tersoffPair(atmType1, atmType3) % REq
          Dik = self%tersoffPair(atmType1, atmType3) % D
          rMax = self%tersoffPair(atmType1, atmType3) % RMax
          rMax_sq = self%tersoffPair(atmType1, atmType3) % RMaxSq
!          rMax_sq = rMax * rMax

          rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
          ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
          rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
          call curbox%Boundary(rxik, ryik, rzik)
          rik = rxik*rxik + ryik*ryik + rzik*rzik
          if(rik < rMax_sq) then
            rik = sqrt(rik)

            c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
            d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
            h = self%tersoffAngle(atmType1, atmType2, atmType3) % h
            lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
            gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
            angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
            Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Reqik, Dik)
          endif
        enddo

        if(Zeta /= 0E0_dp) then
          BetaPar = self%tersoffPair(atmType1, atmType2)%beta
          n =  self%tersoffPair(atmType1, atmType2)%n
          b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
        else
!          b1 = 1E0_dp
          b1 = self%dimer
        endif      
   
        A = self%tersoffPair(atmType1, atmType2) % A
        B = self%tersoffPair(atmType1, atmType2) % B
        lam1 = self%tersoffPair(atmType1, atmType2)%lam1
        lam2 = self%tersoffPair(atmType1, atmType2)%lam2
        V1 = 0.5E0_dp * self%Fc_Func(rij, Reqij, Dij) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
        E_Tersoff = E_Tersoff + V1
        curbox%ETable(iAtom) = curbox%ETable(iAtom) + V1
        curbox%ETable(jAtom) = curbox%ETable(jAtom) + V1
      enddo
    enddo

    write(nout,*) "ETable Sum:", 0.5E0_dp*sum(curbox%ETable)
    write(nout,*) "Tersoff Energy:", E_Tersoff
    E_T = E_Tersoff

   end subroutine

!============================================================================
  subroutine DiffECalc_Tersoff(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

    select type(disp)
      class is(Displacement)
         call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)

      class is(Addition)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

      class is(OrthoVolChange)
         call self % OrthoVolECalc( curbox, disp, E_Diff, accept)

      class default
        write(*,*) "Unknown Perturbation Type Encountered by the Tersoff Forcefield Style."
    end select


  end subroutine

  !=====================================================================
  subroutine Shift_Tersoff_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Displacement), intent(in) :: disp(:)
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
    real(dp) :: lam1, lam2, lam3, gam
    real(dp) :: Zeta, Zeta2
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik
    real(dp) :: E_Tersoff
    integer :: nRecalc
    integer :: recalcList(1:200)

     ! The recalcList is a list of particles whose intermolecular interactions have
     ! changed, but the particles themselves did not move.
    nRecalc = 0
    recalcList = 0

    dispLen = size(disp)
    E_Diff = 0E0_dp
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
        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq

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

          !Compute the Tersoff U_ij component
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if((kAtom == iAtom) .or. (kAtom == jAtom)) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
            rxik = curbox % atoms(1, kAtom)  -  disp(iDisp) % x_New
            ryik = curbox % atoms(2, kAtom)  -  disp(iDisp) % y_New
            rzik = curbox % atoms(3, kAtom)  -  disp(iDisp) % z_New
            call curbox%Boundary(rxik, ryik, rzik)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              rik = sqrt(rik)
              D2 = self%tersoffPair(atmType1, atmType3) % D
              Req = self%tersoffPair(atmType1, atmType3) % Req
              c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
              d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
              h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
            endif     
          enddo
          !Compute the Tersoff U_ji component
          do kNei = 1, curbox%NeighList(1)%nNeigh(jAtom)
            kAtom = curbox%NeighList(1)%list(kNei, jAtom)
            if((kAtom == iAtom) .or. (kAtom == jAtom)) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%tersoffPair(atmType2, atmType3) % rMaxSq

            rxjk = curbox % atoms(1, kAtom)  -  curbox % atoms(1, jAtom)
            ryjk = curbox % atoms(2, kAtom)  -  curbox % atoms(2, jAtom)
            rzjk = curbox % atoms(3, kAtom)  -  curbox % atoms(3, jAtom)
            call curbox%Boundary(rxjk, ryjk, rzjk)
            rjk = rxjk*rxjk + ryjk*ryjk + rzjk*rzjk
            if(rjk < rMaxSq) then
              D2 = self%tersoffPair(atmType2, atmType3) % D
              Req = self%tersoffPair(atmType2, atmType3) % Req
              c = self%tersoffAngle(atmType2, atmType1, atmType3) % c
              d = self%tersoffAngle(atmType2, atmType1, atmType3) % d
              h = self%tersoffAngle(atmType2, atmType1, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              rjk = sqrt(rjk)
              angijk = self%angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta2 = Zeta2 + exp(lam3*(rij-rjk))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rjk, Req, D2)
            endif     
          enddo
          if(Zeta /= 0E0_dp) then
            BetaPar = self%tersoffPair(atmType1, atmType2)%beta
            n =  self%tersoffPair(atmType1, atmType2)%n
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
!            b1 = 1E0_dp
            b1 = self%dimer
          endif
          if(Zeta2 /= 0E0_dp) then
            BetaPar = self%tersoffPair(atmType2, atmType1)%beta
            n =  self%tersoffPair(atmType2, atmType1)%n
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
!            b2 = 1E0_dp
            b2 = self%dimer
          endif

          if(self%symetric) then
            A = self%tersoffPair(atmType1, atmType2) % A
            B = self%tersoffPair(atmType1, atmType2) % B
            lam1 = self%tersoffPair(atmType1, atmType2) % lam1
            lam2 = self%tersoffPair(atmType1, atmType2) % lam2
            Req = self%tersoffPair(atmType1, atmType2) % REq
            D2 = self%tersoffPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          else
            A = self%tersoffPair(atmType1, atmType2) % A
            B = self%tersoffPair(atmType1, atmType2) % B
            lam1 = self%tersoffPair(atmType1, atmType2) % lam1
            lam2 = self%tersoffPair(atmType1, atmType2) % lam2
            Req = self%tersoffPair(atmType1, atmType2) % REq
            D2 = self%tersoffPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 

            A = self%tersoffPair(atmType2, atmType1) % A
            B = self%tersoffPair(atmType2, atmType1) % B
            lam1 = self%tersoffPair(atmType2, atmType1) % lam1
            lam2 = self%tersoffPair(atmType2, atmType1) % lam2
            Req = self%tersoffPair(atmType2, atmType1) % REq
            D2 = self%tersoffPair(atmType2, atmType1) % D
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
        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq
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

          !Compute the Old Tersoff U_ij component
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
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              rik = sqrt(rik)
              D2 = self%tersoffPair(atmType1, atmType3) % D
              Req = self%tersoffPair(atmType1, atmType3) % Req
              c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
              d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
              h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
            endif     
          enddo
          !Compute the Old Tersoff U_ji component
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
            rMaxSq = self%tersoffPair(atmType2, atmType3) % rMaxSq
            rjk = rxjk*rxjk + ryjk*ryjk + rzjk*rzjk
            if(rjk < rMaxSq) then
              rjk = sqrt(rjk)
              D2 = self%tersoffPair(atmType2, atmType3) % D
              Req = self%tersoffPair(atmType2, atmType3) % Req
              c = self%tersoffAngle(atmType2, atmType1, atmType3) % c
              d = self%tersoffAngle(atmType2, atmType1, atmType3) % d
              h = self%tersoffAngle(atmType2, atmType1, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam

              angijk = self%angleCalc(-rxij, -ryij, -rzij, rij, rxjk, ryjk, rzjk, rjk)
              Zeta2 = Zeta2 + exp(lam3*(rij-rjk))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rjk, Req, D2)
            endif     
          enddo
          if(Zeta /= 0E0_dp) then
            BetaPar = self%tersoffPair(atmType1, atmType2) % beta
            n = self%tersoffPair(atmType1, atmType2) % n
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
!            b1 = 1E0_dp
            b1 = self%dimer
          endif
          if(Zeta2 /= 0E0_dp) then
            BetaPar = self%tersoffPair(atmType2, atmType1) % beta
            n = self%tersoffPair(atmType2, atmType1) % n
            b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
!            b2 = 1E0_dp
            b2 = self%dimer
          endif

          if(self%symetric) then
            A = self%tersoffPair(atmType1, atmType2) % A
            B = self%tersoffPair(atmType1, atmType2) % B
            lam1 = self%tersoffPair(atmType1, atmType2) % lam1
            lam2 = self%tersoffPair(atmType1, atmType2) % lam2
            Req = self%tersoffPair(atmType1, atmType2) % REq
            D2 = self%tersoffPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (2d0*A*exp(-lam1*rij) - (b1+b2)*B*exp(-lam2*rij)) 
          else
            A = self%tersoffPair(atmType1, atmType2) % A
            B = self%tersoffPair(atmType1, atmType2) % B
            lam1 = self%tersoffPair(atmType1, atmType2) % lam1
            lam2 = self%tersoffPair(atmType1, atmType2) % lam2
            Req = self%tersoffPair(atmType1, atmType2) % REq
            D2 = self%tersoffPair(atmType1, atmType2) % D
            V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 

            A = self%tersoffPair(atmType2, atmType1) % A
            B = self%tersoffPair(atmType2, atmType1) % B
            lam1 = self%tersoffPair(atmType2, atmType1) % lam1
            lam2 = self%tersoffPair(atmType2, atmType1) % lam2
            Req = self%tersoffPair(atmType2, atmType1) % REq
            D2 = self%tersoffPair(atmType2, atmType1) % D
            V1 = V1 + 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b2*B*exp(-lam2*rij)) 
          
          endif
 
!          write(*,*) "Old:", iAtom, jAtom, -V1
          curbox%dETable(iAtom) = curbox%dETable(iAtom) - V1
          curbox%dETable(jAtom) = curbox%dETable(jAtom) - V1
          E_Diff = E_Diff - V1
        endif
      enddo
    enddo

    !Since the Tersoff is a three body potential, moving a single particle can change the bonded 
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
        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq
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
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
             ! If kAtom is the particle that moved old and new must
             ! be recomputed separately. Otherwise, old and new are the same.
            if(kAtom == disp(1)%atmIndx) then
              rxik = disp(1)%x_new  -  curbox % atoms(1, iAtom)
              ryik = disp(1)%y_new  -  curbox % atoms(2, iAtom)
              rzik = disp(1)%z_new  -  curbox % atoms(3, iAtom)
              call curbox%Boundary(rxik, ryik, rzik)
              rik = rxik*rxik + ryik*ryik + rzik*rzik
              if(rik < rMaxSq) then
                D2 = self%tersoffPair(atmType1, atmType3) % D
                Req = self%tersoffPair(atmType1, atmType3) % Req
                c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
                d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
                h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
                lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
                gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
                rik = sqrt(rik)
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
                Zeta = Zeta + sub
              endif
              rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
              ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
              rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
              call curbox%Boundary(rxik, ryik, rzik)
              rik = rxik*rxik + ryik*ryik + rzik*rzik
              if(rik < rMaxSq) then
                D2 = self%tersoffPair(atmType1, atmType3) % D
                Req = self%tersoffPair(atmType1, atmType3) % Req
                c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
                d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
                h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
                lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
                gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
                rik = sqrt(rik)
                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
                Zeta2 = Zeta2 + sub
              endif
            else
              rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
              ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
              rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
              call curbox%Boundary(rxik, ryik, rzik)
              rik = rxik*rxik + ryik*ryik + rzik*rzik
              if(rik < rMaxSq) then
                D2 = self%tersoffPair(atmType1, atmType3) % D
                Req = self%tersoffPair(atmType1, atmType3) % Req
                c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
                d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
                h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
                lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
                gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
                rik = sqrt(rik)

                angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
                sub = exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
                Zeta = Zeta + sub
                Zeta2 = Zeta2 + sub
              endif
            endif
          enddo
        endif

        BetaPar = self%tersoffPair(atmType1, atmType2) % beta
        n = self%tersoffPair(atmType1, atmType2) % n
        if(Zeta /= 0E0_dp) then
          b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
        else
!          b1 = 1E0_dp
          b1 = self%dimer
        endif
        if(Zeta2 /= 0E0_dp) then
          b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
        else
!          b2 = 1E0_dp
          b2 = self%dimer
        endif    
   
!        A = self%tersoffPair(atmType1, atmType2) % A
        B = self%tersoffPair(atmType1, atmType2) % B
!        lam1 = self%tersoffPair(atmType1, atmType2) % lam1
        lam2 = self%tersoffPair(atmType1, atmType2) % lam2
        D2 = self%tersoffPair(atmType1, atmType2) % D
        Req = self%tersoffPair(atmType1, atmType2) % Req
        V1 = 0.5E0_dp*self%Fc_Func(rij, Req, D2) * (B*exp(-lam2*rij))*(b2 - b1)
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
        E_Diff = E_Diff + V1
      enddo
    enddo
 
!    write(*,*) E_Diff
  end subroutine  
!====================================================================================== 
  subroutine New_Tersoff(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    integer :: iDisp, iAtom, iNei, jNei, jAtom, kNei, kAtom, dispLen, nNei
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2, atmType3
    integer :: listIndx
    real(dp) :: rMaxSq, rMinSq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik
    real(dp) :: sub
    real(dp) :: rmin_ij
    real(dp) :: A, B, c, d, Req, D2 
    real(dp) :: lam1, lam2, lam3, gam
    real(dp) :: Zeta, Zeta2
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik
    real(dp) :: E_Tersoff
    integer :: nRecalc
    integer :: recalcList(1:200)

     ! The recalcList is a list of particles whose intermolecular interactions have
     ! changed, but the particles themselves did not move.
    nRecalc = 0
    recalcList = 0

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      listIndx = disp(iDisp)%listIndex
      nNei = tempNNei(listIndx)
      do jNei = 1, nNei
        jAtom = templist(jNei, listIndx)
        atmType2 = curbox % AtomType(jAtom)

        rxij = curbox % atoms(1, jAtom)  -  disp(iDisp) % x_New
        ryij = curbox % atoms(2, jAtom)  -  disp(iDisp) % y_New
        rzij = curbox % atoms(3, jAtom)  -  disp(iDisp) % z_New
        call curbox%Boundary(rxij, ryij, rzij)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq

        if(rij < rMaxSq) then
          rMinSq = self % rMinTable(atmType1, atmType2)          
          if(rij < rMinSq) then
            accept = .false.
            return
          endif
          if(all(recalcList(1:nRecalc) /= jAtom) ) then
            nRecalc = nRecalc + 1
            recalcList(nRecalc) = jAtom
          endif
          rij = sqrt(rij)
          Zeta = 0E0_dp
          !Compute the Tersoff U_ij component
          do kNei = 1, nNei
            kAtom = templist(kNei, listIndx)
            if(kAtom == jAtom) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
            rxik = curbox % atoms(1, kAtom)  -  disp(iDisp) % x_New
            ryik = curbox % atoms(2, kAtom)  -  disp(iDisp) % y_New
            rzik = curbox % atoms(3, kAtom)  -  disp(iDisp) % z_New
            call curbox%Boundary(rxik, ryik, rzik)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              rik = sqrt(rik)
              D2 = self%tersoffPair(atmType1, atmType3) % D
              Req = self%tersoffPair(atmType1, atmType3) % Req
              c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
              d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
              h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
            endif     
          enddo
          if(Zeta /= 0E0_dp) then
            BetaPar = self%tersoffPair(atmType1, atmType2)%beta
            n =  self%tersoffPair(atmType1, atmType2)%n
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
!            b1 = 1E0_dp
            b1 = self%dimer
          endif


          A = self%tersoffPair(atmType1, atmType2) % A
          B = self%tersoffPair(atmType1, atmType2) % B
          lam1 = self%tersoffPair(atmType1, atmType2) % lam1
          lam2 = self%tersoffPair(atmType1, atmType2) % lam2
          Req = self%tersoffPair(atmType1, atmType2) % REq
          D2 = self%tersoffPair(atmType1, atmType2) % D
          V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
          curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
          E_Diff = E_Diff + V1

        endif
      enddo
    enddo


    !Since the Tersoff is a three body potential, moving a single particle can change the bonded 
    !interaction of particles that did not move.  These interactions must also be recomputed.

    !Compute the remaining U_ji interaction
    jAtom = disp(1)%atmindx
    do iNei = 1, nRecalc
      Zeta = 0E0_dp  !Zeta is the angle term for the new config
      iAtom = recalcList(iNei)
      atmType1 = curbox % AtomType(iAtom)
      atmType2 = curbox % AtomType(disp(1)%atmIndx)
      rxij = disp(1)%x_new  -  curbox % atoms(1, iAtom)
      ryij = disp(1)%y_new  -  curbox % atoms(2, iAtom)
      rzij = disp(1)%z_new  -  curbox % atoms(3, iAtom)
      call curbox%Boundary(rxij, ryij, rzij)
      rij = rxij*rxij + ryij*ryij + rzij*rzij
      rij = sqrt(rij)


      do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        kAtom = curbox%NeighList(1)%list(kNei, iAtom)
        atmType3 = curbox % AtomType(kAtom)
        rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
        rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
        ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
        rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxik, ryik, rzik)
        rik = rxik*rxik + ryik*ryik + rzik*rzik
        if(rik < rMaxSq) then
          rik = sqrt(rik)
          D2 = self%tersoffPair(atmType1, atmType3) % D
          Req = self%tersoffPair(atmType1, atmType3) % Req
          c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
          d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
          h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
          lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
          gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
          angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
          Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
        endif

      enddo
      BetaPar = self%tersoffPair(atmType1, atmType2) % beta
      n = self%tersoffPair(atmType1, atmType2) % n
      if(Zeta /= 0E0_dp) then
        b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
      else
!        b1 = 1E0_dp
        b1 = self%dimer
      endif

      A = self%tersoffPair(atmType1, atmType2) % A
      B = self%tersoffPair(atmType1, atmType2) % B
      lam1 = self%tersoffPair(atmType1, atmType2) % lam1
      lam2 = self%tersoffPair(atmType1, atmType2) % lam2
      Req = self%tersoffPair(atmType1, atmType2) % REq
      D2 = self%tersoffPair(atmType1, atmType2) % D
      V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 
      curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
      curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
      E_Diff = E_Diff + V1
    enddo
 

    if(nRecalc == 0) then
      return
    endif
    !Compute the new interactions between the particles already in the system
    do iNei = 1, nRecalc
      iAtom = recalcList(iNei)
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)

        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        atmType2 = curbox % AtomType(jAtom)
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxij, ryij, rzij)

        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij < rMaxSq) then
          Zeta = 0E0_dp  !Zeta is the angle term for the new config
          Zeta2 = 0E0_dp !Zeta is the angle term for the old config
          rij = sqrt(rij)
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if(kAtom == jAtom) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
             ! If kAtom is the particle that moved old and new must
             ! be recomputed separately. Otherwise, old and new are the same.
            rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
            ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
            rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
            call curbox%Boundary(rxik, ryik, rzik)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              D2 = self%tersoffPair(atmType1, atmType3) % D
              Req = self%tersoffPair(atmType1, atmType3) % Req
              c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
              d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
              h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              rik = sqrt(rik)

              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              sub = exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
              Zeta = Zeta + sub
              Zeta2 = Zeta2 + sub
            endif
          enddo

        !
          atmType3 = curbox % AtomType(disp(1)%atmIndx)
          rxik = disp(1)%x_New  -  curbox % atoms(1, iAtom)
          ryik = disp(1)%y_New  -  curbox % atoms(2, iAtom)
          rzik = disp(1)%z_New  -  curbox % atoms(3, iAtom)
          call curbox%Boundary(rxik, ryik, rzik)
          rik = rxik*rxik + ryik*ryik + rzik*rzik
          rik = sqrt(rik)
          D2 = self%tersoffPair(atmType1, atmType3) % D
          Req = self%tersoffPair(atmType1, atmType3) % Req
          lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
          c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
          d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
          h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
          gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
          angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
          sub = exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
          Zeta = Zeta + sub

          BetaPar = self%tersoffPair(atmType1, atmType2) % beta
          n = self%tersoffPair(atmType1, atmType2) % n
          if(Zeta /= 0E0_dp) then
              b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
!              b1 = 1E0_dp
              b1 = self%dimer
          endif
          if(Zeta2 /= 0E0_dp) then
              b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
!              b2 = 1E0_dp
              b2 = self%dimer
          endif    

          B = self%tersoffPair(atmType1, atmType2) % B
          lam2 = self%tersoffPair(atmType1, atmType2) % lam2
          D2 = self%tersoffPair(atmType1, atmType2) % D
          Req = self%tersoffPair(atmType1, atmType2) % Req
          V1 = 0.5E0_dp*self%Fc_Func(rij, Req, D2) * (B*exp(-lam2*rij))*(b2 - b1)
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
          curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
          E_Diff = E_Diff + V1
        endif
      enddo
    enddo

!    write(*,*) E_Diff
!    write(*,*) curbox%dETable(1:3)
!    write(*,*) curbox%ETable(1:3)+curbox%dETable(1:3)
  end subroutine  
!====================================================================================== 
  subroutine Old_Tersoff(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: atmType1, atmType2, atmType3
    integer :: molStart, molEnd, iAtom, iMol, jNei, jAtom
    integer :: kNei, iNei, kAtom, nAtom
    real(dp) :: rMaxSq, rMinSq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik
    real(dp) :: sub
    real(dp) :: rmin_ij
    real(dp) :: A, B, c, d, Req, D2 
    real(dp) :: lam1, lam2, lam3, gam
    real(dp) :: Zeta, Zeta2
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1, V2
    real(dp) :: angijk, angjik
    real(dp) :: E_Tersoff
    integer :: nRecalc
    integer :: recalcList(1:200)

     ! The recalcList is a list of particles whose intermolecular interactions have
     ! changed, but the particles themselves did not move.
    nRecalc = 0
    recalcList = 0
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)
    E_Diff = 0E0_dp


    iAtom = molStart
    atmType1 = curbox % AtomType(iAtom)
    do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        atmType2 = curbox % AtomType(jAtom)

        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxij, ryij, rzij)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq
        if(rij < rMaxSq) then
          if(all(recalcList(1:nRecalc) /= jAtom) ) then
            nRecalc = nRecalc + 1
            recalcList(nRecalc) = jAtom
          endif
          rij = sqrt(rij)
          Zeta = 0E0_dp
          !Compute the Tersoff U_ij component
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)

            if( kAtom == jAtom ) then
              cycle
            endif

            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
            rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
            ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
            rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
            call curbox%Boundary(rxik, ryik, rzik)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              rik = sqrt(rik)
              D2 = self%tersoffPair(atmType1, atmType3) % D
              Req = self%tersoffPair(atmType1, atmType3) % Req
              c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
              d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
              h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
              Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
            endif     
          enddo
          if(Zeta /= 0E0_dp) then
            BetaPar = self%tersoffPair(atmType1, atmType2)%beta
            n =  self%tersoffPair(atmType1, atmType2)%n
            b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
!             b1 = 1E0_dp
             b1 = self%dimer
          endif


          A = self%tersoffPair(atmType1, atmType2) % A
          B = self%tersoffPair(atmType1, atmType2) % B
          lam1 = self%tersoffPair(atmType1, atmType2) % lam1
          lam2 = self%tersoffPair(atmType1, atmType2) % lam2
          Req = self%tersoffPair(atmType1, atmType2) % REq
          D2 = self%tersoffPair(atmType1, atmType2) % D
          V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 
          curbox%dETable(iAtom) = curbox%dETable(iAtom) - V1
          curbox%dETable(jAtom) = curbox%dETable(jAtom) - V1
          E_Diff = E_Diff - V1

        endif
      enddo

    !Since the Tersoff is a three body potential, moving a single particle can change the bonded 
    !interaction of particles that did not move.  These interactions must also be recomputed.

    !Compute the remaining U_ji interaction
    jAtom = molStart
    atmType2 = curbox % AtomType(jAtom)
    do iNei = 1, nRecalc
      Zeta = 0E0_dp  !Zeta is the angle term for the new config
      iAtom = recalcList(iNei)
      atmType1 = curbox % AtomType(iAtom)
      rxij = curbox % atoms(1, jAtom) - curbox % atoms(1, iAtom)
      ryij = curbox % atoms(2, jAtom) - curbox % atoms(2, iAtom)
      rzij = curbox % atoms(3, jAtom) - curbox % atoms(3, iAtom)
      call curbox%Boundary(rxij, ryij, rzij)
      rij = rxij*rxij + ryij*ryij + rzij*rzij
      rij = sqrt(rij)


      do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        kAtom = curbox%NeighList(1)%list(kNei, iAtom)
        if( kAtom == jAtom ) then
          cycle
        endif


        atmType3 = curbox % AtomType(kAtom)
        rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
        rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
        ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
        rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxik, ryik, rzik)
        rik = rxik*rxik + ryik*ryik + rzik*rzik
        if(rik < rMaxSq) then
          rik = sqrt(rik)
          D2 = self%tersoffPair(atmType1, atmType3) % D
          Req = self%tersoffPair(atmType1, atmType3) % Req
          lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
          c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
          d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
          h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
          gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
          angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
          Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Req, D2)
        endif

      enddo
      BetaPar = self%tersoffPair(atmType1, atmType2) % beta
      n = self%tersoffPair(atmType1, atmType2) % n
      if(Zeta /= 0E0_dp) then
        b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
      else
!        b1 = 1E0_dp
         b1 = self%dimer
      endif

      A = self%tersoffPair(atmType1, atmType2) % A
      B = self%tersoffPair(atmType1, atmType2) % B
      lam1 = self%tersoffPair(atmType1, atmType2) % lam1
      lam2 = self%tersoffPair(atmType1, atmType2) % lam2
      Req = self%tersoffPair(atmType1, atmType2) % REq
      D2 = self%tersoffPair(atmType1, atmType2) % D
      V1 = 0.5E0_dp * self%Fc_Func(rij, Req, D2) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij)) 
      curbox%dETable(iAtom) = curbox%dETable(iAtom) - V1
      curbox%dETable(jAtom) = curbox%dETable(jAtom) - V1
      E_Diff = E_Diff - V1
    enddo
 

    if(nRecalc == 0) then
      return
    endif
    !Compute the new interactions between the particles already in the system
    nAtom = molStart
    do iNei = 1, nRecalc
      iAtom = recalcList(iNei)
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        

        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        if(jAtom == nAtom) then
          cycle
        endif

        atmType2 = curbox % AtomType(jAtom)
        rxij = curbox % atoms(1, jAtom)  -  curbox % atoms(1, iAtom)
        ryij = curbox % atoms(2, jAtom)  -  curbox % atoms(2, iAtom)
        rzij = curbox % atoms(3, jAtom)  -  curbox % atoms(3, iAtom)
        call curbox%Boundary(rxij, ryij, rzij)

        rMaxSq = self%tersoffPair(atmType1, atmType2) % rMaxSq
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij < rMaxSq) then
          Zeta = 0E0_dp  !Zeta is the angle term for the new config
          Zeta2 = 0E0_dp !Zeta is the angle term for the old config
          rij = sqrt(rij)
          do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
            kAtom = curbox%NeighList(1)%list(kNei, iAtom)
            if(kAtom == jAtom) then
              cycle
            endif
            atmType3 = curbox % AtomType(kAtom)
            rMaxSq = self%tersoffPair(atmType1, atmType3) % rMaxSq
             ! If kAtom is the particle that moved old and new must
             ! be recomputed separately. Otherwise, old and new are the same.
            rxik = curbox % atoms(1, kAtom)  -  curbox % atoms(1, iAtom)
            ryik = curbox % atoms(2, kAtom)  -  curbox % atoms(2, iAtom)
            rzik = curbox % atoms(3, kAtom)  -  curbox % atoms(3, iAtom)
            call curbox%Boundary(rxik, ryik, rzik)
            rik = rxik*rxik + ryik*ryik + rzik*rzik
            if(rik < rMaxSq) then
              D2 = self%tersoffPair(atmType1, atmType3) % D
              Req = self%tersoffPair(atmType1, atmType3) % Req
              lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
              c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
              d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
              h = self%tersoffAngle(atmType1, atmType2, atmType3) % h            
              gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
              rik = sqrt(rik)

              angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
!              sub = self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
              sub = exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) *  self%Fc_Func(rik, Req, D2)
              if(kAtom == nAtom) then
                Zeta2 = Zeta2 + sub
              else
                Zeta = Zeta + sub
                Zeta2 = Zeta2 + sub
              endif
            endif
          enddo

          BetaPar = self%tersoffPair(atmType1, atmType2) % beta
          n = self%tersoffPair(atmType1, atmType2) % n
          if(Zeta /= 0E0_dp) then
              b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
          else
!              b1 = 1E0_dp
              b1 = self%dimer
          endif
          if(Zeta2 /= 0E0_dp) then
              b2 = (1E0_dp + (BetaPar*Zeta2)**n)**(-1E0_dp/(2E0_dp*n))
          else
!              b2 = 1E0_dp
              b2 = self%dimer
          endif    

          B = self%tersoffPair(atmType1, atmType2) % B
          lam2 = self%tersoffPair(atmType1, atmType2) % lam2
          D2 = self%tersoffPair(atmType1, atmType2) % D
          Req = self%tersoffPair(atmType1, atmType2) % Req
          V1 = 0.5E0_dp*self%Fc_Func(rij, Req, D2) * (B*exp(-lam2*rij))*(b2 - b1)
          curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
          curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
          E_Diff = E_Diff + V1
        endif
      enddo
    enddo

  end subroutine
  !=====================================================================
  subroutine OrthoVol_Tersoff(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(OrthoVolChange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    integer :: jNei, kNei
    integer :: iAtom, jAtom, kAtom
    integer :: atmType1, atmType2, atmType3
    integer :: molIndx1, molIndx2, molIndx3
    real(dp) :: A, B, c, d, Reqij, Reqik, Dij, Dik
    real(dp) :: E_Tersoff
    real(dp) :: lam1, lam2, lam3, gam
    real(dp) :: Zeta
    real(dp) :: BetaPar, n, h
    real(dp) :: b1, b2, V1
    real(dp) :: angijk, angjik

    real(dp) :: rMax, rMax_sq
    real(dp) :: rxij, ryij, rzij, rij
    real(dp) :: rxjk, ryjk, rzjk, rjk
    real(dp) :: rxik, ryik, rzik, rik
    real(dp) :: dxi, dyi, dzi
    real(dp) :: dxj, dyj, dzj
    real(dp) :: dxk, dyk, dzk


    accept = .true.
    do iAtom = 1, curbox%nMaxAtoms
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
        cycle
      endif
      atmType1 = curbox % AtomType(iAtom)
      molIndx1 = curbox % MolIndx(iAtom)
      dxi = curbox % centerMass(1, molIndx1) * (disp(1)%xScale-1E0_dp)
      dyi = curbox % centerMass(2, molIndx1) * (disp(1)%yScale-1E0_dp)
      dzi = curbox % centerMass(3, molIndx1) * (disp(1)%zScale-1E0_dp)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)


        atmType2 = curbox % AtomType(jAtom)
        molIndx2 = curbox % MolIndx(jAtom)
        Reqij = self%tersoffPair(atmType2, atmType1) % REq
        Dij = self%tersoffPair(atmType2, atmType1) % D
        rMax = Dij + Reqij
        rMax_sq = rMax * rMax
        dxj = curbox % centerMass(1, molIndx2) * (disp(1)%xScale-1E0_dp)
        dyj = curbox % centerMass(2, molIndx2) * (disp(1)%yScale-1E0_dp)
        dzj = curbox % centerMass(3, molIndx2) * (disp(1)%zScale-1E0_dp)

        rxij = curbox % atoms(1, jAtom) + dxj - curbox % atoms(1, iAtom) - dxi
        ryij = curbox % atoms(2, jAtom) + dyj - curbox % atoms(2, iAtom) - dyi
        rzij = curbox % atoms(3, jAtom) + dzj - curbox % atoms(3, iAtom) - dzi
        call curbox%BoundaryNew(rxij, ryij, rzij, disp)
        rij = rxij*rxij + ryij*ryij + rzij*rzij
        if(rij > rMax_Sq) then
          cycle
        endif
        rij = sqrt(rij)
        Zeta = 0E0_dp
        do kNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
          kAtom = curbox%NeighList(1)%list(kNei, iAtom)
          if(kAtom == jAtom) then
            cycle
          endif

          atmType3 = curbox % AtomType(kAtom)
          molIndx3 = curbox % MolIndx(kAtom)
          Reqik = self%tersoffPair(atmType1, atmType3) % REq
          Dik = self%tersoffPair(atmType1, atmType3) % D
          rMax = Reqik + Dik
          rMax_sq = rMax * rMax
          dxk = curbox % centerMass(1, molIndx3) * (disp(1)%xScale-1E0_dp)
          dyk = curbox % centerMass(2, molIndx3) * (disp(1)%yScale-1E0_dp)
          dzk = curbox % centerMass(3, molIndx3) * (disp(1)%zScale-1E0_dp)
          rxik = curbox % atoms(1, kAtom) + dxk - curbox % atoms(1, iAtom) - dxi
          ryik = curbox % atoms(2, kAtom) + dyk - curbox % atoms(2, iAtom) - dyi
          rzik = curbox % atoms(3, kAtom) + dzk - curbox % atoms(3, iAtom) - dzi
          call curbox%BoundaryNew(rxik, ryik, rzik, disp)
          rik = rxik*rxik + ryik*ryik + rzik*rzik
          if(rik < rMax_sq) then
            rik = sqrt(rik)

            c = self%tersoffAngle(atmType1, atmType2, atmType3) % c
            d = self%tersoffAngle(atmType1, atmType2, atmType3) % d
            h = self%tersoffAngle(atmType1, atmType2, atmType3) % h
            lam3 = self%tersoffAngle(atmType1, atmType2, atmType3) % lam3
            gam = self%tersoffAngle(atmType1, atmType2, atmType3) % gam
            angijk = self%angleCalc(rxij, ryij, rzij, rij, rxik, ryik, rzik, rik)
            Zeta = Zeta + exp(lam3*(rij-rik))*gam*self%gik_Func(angijk, c, d, h) * self%Fc_Func(rik, Reqik, Dik)
          endif
        enddo

        if(Zeta /= 0E0_dp) then
          BetaPar = self%tersoffPair(atmType1, atmType2)%beta
          n =  self%tersoffPair(atmType1, atmType2)%n
          b1 = (1E0_dp + (BetaPar*Zeta)**n)**(-1E0_dp/(2E0_dp*n))
        else
!          b1 = 1E0_dp
          b1 = self%dimer
        endif      
   
        A = self%tersoffPair(atmType1, atmType2) % A
        B = self%tersoffPair(atmType1, atmType2) % B
        lam1 = self%tersoffPair(atmType1, atmType2)%lam1
        lam2 = self%tersoffPair(atmType1, atmType2)%lam2
        V1 = 0.5E0_dp * self%Fc_Func(rij, Reqij, Dij) * (A*exp(-lam1*rij) - b1*B*exp(-lam2*rij))
        E_Diff = E_Diff + V1
        curbox%dETable(iAtom) = curbox%dETable(iAtom) + V1
        curbox%dETable(jAtom) = curbox%dETable(jAtom) + V1
      enddo
    enddo

    E_Diff = E_Diff - curbox%E_Inter
    curbox%dETable = curbox%dETable - curbox%ETable

  end subroutine
!=====================================================================
  subroutine ProcessIO_Tersoff(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetAllCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(Pair_Tersoff), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

    character(len=30) :: command
    logical :: logicVal
    integer :: iPar, lineStat
    integer :: type1, type2, type3
    real(dp) :: parList(1:8), realVal
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("symetric")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) logicVal
        self%symetric = logicVal

      case("rmin")
        do iPar = 1, nAtomTypes
          call GetXCommand(line, command, iPar+1, lineStat)
          read(command, *) realVal
          self%rMin(iPar) = realVal
        enddo

        do type1 = 1, nAtomTypes
          do type2 = 1, nAtomTypes
              self%rMinTable(type1, type2) = max(self%rMin(type1), &
                                                 self%rMin(type2))**2

              self%rMinTable(type2, type1) = max(self%rMin(type1), &
                                                 self%rMin(type2))**2
          enddo
        enddo

      case("dimer")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) realVal
        self%dimer = realVal


      case("pair")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) type1
        call GetXCommand(line, command, 3, lineStat)
        read(command, *) type2
        do iPar = 1, 8
          call GetXCommand(line, command, 3+iPar, lineStat)
          read(command, *) parList(iPar)
        enddo
        self%tersoffPair(type1, type2)%A = parList(1) * inEngUnit
        self%tersoffPair(type1, type2)%B = parList(2) * inEngUnit
        self%tersoffPair(type1, type2)%lam1 = parList(3) * inLenUnit
        self%tersoffPair(type1, type2)%lam2 = parList(4) * inLenUnit
        self%tersoffPair(type1, type2)%Req = parList(5) * inLenUnit
        self%tersoffPair(type1, type2)%D = parList(6) * inLenUnit
        self%tersoffPair(type1, type2)%beta = parList(7)
        self%tersoffPair(type1, type2)%n = parList(8)
        self%tersoffPair(type1, type2)%rMax = (parList(5) + parList(6))*inLenUnit
        self%tersoffPair(type1, type2)%rMaxSq = ((parList(5) + parList(6))*inLenUnit)**2
        if(self%symetric) then
          self%tersoffPair(type2, type1)%A = parList(1) * inEngUnit
          self%tersoffPair(type2, type1)%B = parList(2) * inEngUnit
          self%tersoffPair(type2, type1)%lam1 = parList(3) * inLenUnit
          self%tersoffPair(type2, type1)%lam2 = parList(4) * inLenUnit
          self%tersoffPair(type2, type1)%Req = parList(5) * inLenUnit
          self%tersoffPair(type2, type1)%D = parList(6) * inLenUnit
          self%tersoffPair(type2, type1)%beta = parList(7)
          self%tersoffPair(type2, type1)%n = parList(8)
          self%tersoffPair(type2, type1)%rMax = (parList(5) + parList(6))*inLenUnit
          self%tersoffPair(type2, type1)%rMaxSq = ((parList(5) + parList(6))*inLenUnit)**2
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
        self%tersoffAngle(type1, type2, type3)%h = parList(1)
        self%tersoffAngle(type1, type2, type3)%lam3 = parList(2)* inLenUnit
        self%tersoffAngle(type1, type2, type3)%c = parList(3)
        self%tersoffAngle(type1, type2, type3)%d = parList(4)
        self%tersoffAngle(type1, type2, type3)%gam = parList(5)
        if(self%symetric) then
          self%tersoffAngle(type2, type1, type3)%h = parList(1)
          self%tersoffAngle(type2, type1, type3)%lam3 = parList(2)* inLenUnit
          self%tersoffAngle(type2, type1, type3)%c = parList(3)
          self%tersoffAngle(type2, type1, type3)%d = parList(4)
          self%tersoffAngle(type2, type1, type3)%gam = parList(5)
        endif

      case default
        linestat = -1
    end select


  end subroutine

 !=============================================================================+
    function GetCutOff_Tersoff(self) result(rCut)
      use Common_MolInfo, only: nAtomTypes
      implicit none
      class(Pair_Tersoff), intent(inout) :: self
      real(dp) :: rCut

      integer :: type1, type2
      real(dp) :: Req, D, rMax

      rCut = -1E0_dp
      do type1 = 1, nAtomTypes
        do type2 = 1, nAtomTypes
          rMax = self%tersoffPair(type1, type2) % RMax
          if(rCut < rMax) then
            rCut = rMax
          endif
        enddo
      enddo

      self%rCut = rCut
      self%rCutSq = rCut*rCut
    end function
  !=====================================================================

end module
