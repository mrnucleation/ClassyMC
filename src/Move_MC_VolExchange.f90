!=========================================================================
! Monte Carlo move for the Isobaric ensemble. This move scales the dimensions
! of the simulation box uniformly in all directions. 
!=========================================================================
module MCMove_VolExchange
  use CoordinateTypes, only: OrthoVolChange, TriVolChange
  use SimpleSimBox, only: SimpleBox
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use VarPrecision
!  use MoveClassDef
  use MultiBoxMoveDef, only: MCMultiBoxMove

  type, public, extends(MCMultiBoxMove) :: VolExchange
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp   
    integer :: style = 1
    real(dp) :: maxDv = 0.01E0_dp
    logical :: tuneMax = .true.
    real(dp) :: limit = 3.00E0_dp
    real(dp) :: targAccpt = 50E0_dp

    integer :: nDim = 3
    type(OrthoVolChange) :: disp1(1:1)
    type(OrthoVolChange) :: disp2(1:1)
!    type(TriVolChange) :: disptri(1:1)
    contains
      procedure, pass :: Constructor => VolExchange_Constructor
!      procedure, pass :: FullMove => VolExchange_FullMove
      procedure, pass :: MultiBox => VolExchange_MultiBox
!      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance => VolExchange_Maintenance
      procedure, pass :: Prologue => VolExchange_Prologue
      procedure, pass :: Epilogue => VolExchange_Epilogue
      procedure, pass :: ProcessIO => VolExchange_ProcessIO
  end type

 contains
!========================================================
  subroutine VolExchange_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(VolExchange), intent(inout) :: self
    integer :: nBoxes



    allocate( self%tempNNei(1) )
    allocate( self%tempList(1, 1) )
  end subroutine
!=========================================================================
  subroutine VolExchange_MultiBox(self, accept)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use RandomGen, only: grnd, ListRNG
    use CommonSampling, only: sampling

    implicit none
    class(VolExchange), intent(inout) :: self
    logical, intent(out) :: accept
    integer :: i, boxNum
    class(SimpleBox), pointer :: box1, box2
    real(dp) :: dV
    real(dp) :: Prob, Norm, extraTerms, half
    real(dp) :: E_Diff1, E_Diff2, scaleFactor
    real(dp) :: rescale(1:size(self%boxprob))



    !Randomly choose which boxes will exchange volume
    boxNum = ListRNG(self%boxProb)
    box1 => BoxArray(boxNum)%box

    !To avoid picking the same box twice, rescale the probability such that
    !box1's probability is equal to 0.
    rescale = self%boxprob
    rescale(boxNum) = 0E0_dp
    norm = 0E0_dp
    do i = 1, size(rescale)
      norm = norm + rescale(i)
    enddo
    if(norm == 0E0_dp) then
      write(0,*) "WARNING! To use VolExchange a nonzero probability must"
      write(0,*) "be set for more than one box!"
      stop 
    endif
    do i = 1, size(rescale)
      rescale(i) = rescale(i)/norm
    enddo
    boxNum = ListRNG(rescale)
    box2 => BoxArray(boxNum)%box



    !Increment attempt counter.
    self % atmps = self % atmps + 1E0_dp

    !Randomly chose the amount of volume that will be exchanged.
    select case(self%style)
      case(1) !Log Scale
        dV = self%maxDv * (2E0_dp*grnd()-1E0_dp)
        dV = exp(dV) - 1E0_dp
        self%disp1(1)%volNew = box1%volume + dV
        self%disp1(1)%volOld = box1%volume

        self%disp2(1)%volNew = box2%volume - dV
        self%disp2(1)%volOld = box2%volume
      case(2) !Linear Scale
        dV = self%maxDv * (2E0_dp*grnd()-1E0_dp)
        self%disp1(1)%volNew = box1%volume + dV
        self%disp1(1)%volOld = box1%volume

        self%disp2(1)%volNew = box2%volume - dV
        self%disp2(1)%volOld = box2%volume
      case default
        stop

    end select
    if(self%disp1(1)%volNew < 0E0_dp) then
      return
    endif
    if(self%disp2(1)%volNew < 0E0_dp) then
      return
    endif


    !Compute how much to scale the sides of the box by for box1
    select type(box1)
      class is(CubeBox)
        scaleFactor = (self%disp1(1)%volNew/self%disp1(1)%volOld)**(1E0_dp/3E0_dp)
        self%disp1(1)%xScale = scaleFactor
        self%disp1(1)%yScale = scaleFactor
        self%disp1(1)%zScale = scaleFactor

      class is(OrthoBox)
        scaleFactor = (self%disp1(1)%volNew/self%disp1(1)%volOld)**(1E0_dp/3E0_dp)
        self%disp1(1)%xScale = scaleFactor
        self%disp1(1)%yScale = scaleFactor
        self%disp1(1)%zScale = scaleFactor

      class default
        stop "This type of box is not compatible with volume change moves."
    end select

    !Compute how much to scale the sides of the box by for box2
    select type(box2)
      class is(CubeBox)
        scaleFactor = (self%disp2(1)%volNew/self%disp2(1)%volOld)**(1E0_dp/3E0_dp)
        self%disp2(1)%xScale = scaleFactor
        self%disp2(1)%yScale = scaleFactor
        self%disp2(1)%zScale = scaleFactor

      class is(OrthoBox)
        scaleFactor = (self%disp2(1)%volNew/self%disp2(1)%volOld)**(1E0_dp/3E0_dp)
        self%disp2(1)%xScale = scaleFactor
        self%disp2(1)%yScale = scaleFactor
        self%disp2(1)%zScale = scaleFactor

      class default
        stop "This type of box is not compatible with volume change moves."
    end select




    !Check Constraint of Box1
    accept = box1 % CheckConstraint( self%disp1(1:1) )
    if(.not. accept) then
      return
    endif
    !Check Constraint of Box2
    accept = box2 % CheckConstraint( self%disp2(1:1) )
    if(.not. accept) then
      return
    endif


    !Energy Calculation for Box 1
    call box1 % EFunc % Method % DiffECalc(box1, self%disp1(1:1), self%tempList, self%tempNNei, E_Diff1, accept)
    if(.not. accept) then
      return
    endif
    !Check Post Energy Constraint for Box 1
    accept = box1 % CheckPostEnergy( self%disp1(1:1), E_Diff1 )
    if(.not. accept) then
      return
    endif

    !Energy Calculation for Box 2
    call box2 % EFunc % Method % DiffECalc(box2, self%disp2(1:1), self%tempList, self%tempNNei, E_Diff2, accept)
    if(.not. accept) then
      return
    endif
    !Check Post Energy Constraint for Box 2
    accept = box2 % CheckPostEnergy( self%disp2(1:1), E_Diff2 )
    if(.not. accept) then
      return
    endif





    !Compute the Proposal probability for Box1
    select case(self%style)
      case(1) !Log Scale
        prob = (box1%nMolTotal+1) * log(self%disp1(1)%volNew / self%disp1(1)%volOld) 
      case(2) !Linear Scale
        prob = box1%nMolTotal * log(self%disp1(1)%volNew / self%disp1(1)%volOld) 
    end select

    !Compute the Proposal probability for Box2
    select case(self%style)
      case(1) !Log Scale
        prob = prob + (box2%nMolTotal+1) * log(self%disp2(1)%volNew / self%disp2(1)%volOld) 
      case(2) !Linear Scale
        prob = prob + box2%nMolTotal * log(self%disp2(1)%volNew / self%disp2(1)%volOld) 
    end select


    !Get the PV contribution to the Boltzmann weight for both boxes.
    extraTerms = sampling % GetExtraTerms(self%disp1(1:1), box1)
    half = sampling % GetExtraTerms(self%disp2(1:1), box2)
    extraTerms = extraTerms + half

!    write(*,*)  E_Diff1, E_Diff2, extraTerms, prob
!    accept = sampling % MakeDecision(box1, E_Diff, self%disp1(1:1), logProb=prob,&
!                                     extraIn=extraTerms)
    accept = sampling % MakeDecision2Box(box1,  box2, E_Diff1, E_Diff2, &
                            self%disp1(1:1), self%disp2(1:1), logProb=prob, &
                            extraIn=extraTerms )

    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call box1 % UpdateEnergy(E_Diff1)
      call box1 % UpdatePosition(self%disp1(1:1), self%tempList, self%tempNNei)

      call box2 % UpdateEnergy(E_Diff2)
      call box2 % UpdatePosition(self%disp2(1:1), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine VolExchange_Maintenance(self)
    implicit none
    class(VolExchange), intent(inout) :: self
!    real(dp), parameter :: limit = 3.0E0_dp
      
    if(self%tuneMax) then
      if(self%atmps < 0.5E0_dp) then
        return
      endif

      if(self%GetAcceptRate() > self%targAccpt) then
        if(self%maxDv*1.01E0_dp < self%limit) then
          self%maxDv = self%maxDv * 1.01E0_dp
        else 
          self%maxDv = self%limit       
        endif
      else
        self%maxDv = self%maxDv * 0.99E0_dp
      endif
    endif

  end subroutine
!=========================================================================
  subroutine VolExchange_Prologue(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(VolExchange), intent(inout) :: self
    integer :: nBoxes

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    write(nout,"(1x,A,F15.8)") "(Volume Exchange) Maximum Volume Change: ", self%maxdV

  end subroutine
!=========================================================================
  subroutine VolExchange_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(VolExchange), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Volume Exchange Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Volume Exchange Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Volume Exchange Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,F15.8)") "Final Maximum Volume Change: ", self%maxDv
    endif
 

  end subroutine
!=========================================================================
  subroutine VolExchange_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(VolExchange), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
      case("tunemax")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%tunemax = logicVal

      case("dynamiclimit")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%limit = realVal

      case("dynamictarget")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%targAccpt = realVal

      case("maxdv")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%maxDv = realVal

      case("style")
        call GetXCommand(line, command, 5, lineStat)
        select case(trim(adjustl(command)))
          case("log")
            self%style = 1
          case("linear")
            self%style = 2
        end select
      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine

!=========================================================================
end module
!=========================================================================
