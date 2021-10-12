!=========================================================================
! Monte Carlo move for the Isobaric ensemble. This move scales the dimensions
! of the simulation box uniformly in all directions. 
!=========================================================================
module MCMove_Isovol
  use CoordinateTypes, only: OrthoVolChange, TriVolChange
  use SimpleSimBox, only: SimpleBox
  use CubicBoxDef, only: CubeBox
  use OrthoBoxDef, only: OrthoBox
  use VarPrecision
  use MoveClassDef

  type, public, extends(MCMove) :: IsoVol
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp   

    
    logical :: verbose = .true.
    integer :: style = 1
    real(dp) :: maxDv = 0.01E0_dp
    logical :: tuneMax = .true.
    real(dp) :: limit = 3.00E0_dp
    real(dp) :: targAccpt = 50E0_dp

    !Rejection Counters
    integer :: volrejection = 0 
    integer :: ovlaprej = 0 
    integer :: constrainrej = 0 
    integer :: detailedrej = 0 

    integer :: nDim = 3
    type(OrthoVolChange) :: disp(1:1)
    type(TriVolChange) :: disptri(1:1)
    contains
      procedure, pass :: Constructor => IsoVol_Constructor
      procedure, pass :: FullMove => IsoVol_FullMove
!      procedure, pass :: GetAcceptRate
      procedure, pass :: Maintenance => IsoVol_Maintenance
      procedure, pass :: Prologue => IsoVol_Prologue
      procedure, pass :: Epilogue => IsoVol_Epilogue
      procedure, pass :: ProcessIO => IsoVol_ProcessIO
  end type

 contains
!========================================================
  subroutine IsoVol_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(IsoVol), intent(inout) :: self


  end subroutine
!=========================================================================
  subroutine IsoVol_FullMove(self, trialBox, accept)
    use Common_MolInfo, only: nMolTypes
    use Box_Utility, only: FindAtom, FindFirstEmptyMol
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    use CommonSampling, only: sampling

    implicit none
    class(IsoVol), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: i
    real(dp) :: dV
    real(dp) :: OldProb, NewProb, Prob, extraTerms
    real(dp) :: E_Diff, E_Inter, E_Intra, scaleFactor

    self % atmps = self % atmps + 1E0_dp
    select case(self%style)
      case(1) !Log Scale
        dV = self%maxDv * (2E0_dp*grnd()-1E0_dp)
        self%disp(1)%volNew = trialBox%volume * exp(dV)
        self%disp(1)%volOld = trialBox%volume

      case(2) !Linear Scale
        dV = self%maxDv * (2E0_dp*grnd()-1E0_dp)
        self%disp(1)%volNew = trialBox%volume + dV
        self%disp(1)%volOld = trialBox%volume
      case default
        stop

    end select

    if(self%disp(1)%volNew < 0E0_dp) then
      return
    endif

    select type(trialBox)
      class is(CubeBox)
        scaleFactor = (self%disp(1)%volNew/self%disp(1)%volOld)**(1E0_dp/3E0_dp)
        self%disp(1)%xScale = scaleFactor
        self%disp(1)%yScale = scaleFactor
        self%disp(1)%zScale = scaleFactor

      class is(OrthoBox)
        scaleFactor = (self%disp(1)%volNew/self%disp(1)%volOld)**(1E0_dp/3E0_dp)
        self%disp(1)%xScale = scaleFactor
        self%disp(1)%yScale = scaleFactor
        self%disp(1)%zScale = scaleFactor

      class default
        stop "This type of box is not compatible with volume change moves."
    end select

    if(.not. allocated(self%tempNNei)) allocate( self%tempNNei(1) )
    if(.not. allocated(self%templist)) allocate( self%templist(1,1) )
    

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:1) )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif

    !Energy Calculation
!    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:1), self%tempList, self%tempNNei, E_Diff, accept)
    call trialBox%ComputeEnergyDelta(self%disp(1:1),&
                                     self%templist, &
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.false.)
    if(.not. accept) then
      self%ovlaprej = self%ovlaprej + 1
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:1), E_Diff )
    if(.not. accept) then
      self%constrainrej = self%constrainrej + 1
      return
    endif


!    write(*,*) E_Inter

    select case(self%style)
      case(1) !Log Scale
        prob = (trialBox%nMolTotal+1) * log(self%disp(1)%volNew / self%disp(1)%volOld) 
      case(2) !Linear Scale
        prob = trialBox%nMolTotal * log(self%disp(1)%volNew / self%disp(1)%volOld) 
    end select

    !Get PV term
    extraTerms = sampling % GetExtraTerms(self%disp(1:1), trialBox)
    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:1), logProb=prob, extraIn=extraTerms)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter, E_Intra)
      call trialBox % UpdatePosition(self%disp(1:1), self%tempList, self%tempNNei)
    else
      self%detailedrej = self%detailedrej + 1
    endif

  end subroutine
!=========================================================================
  subroutine IsoVol_Maintenance(self)
    implicit none
    class(IsoVol), intent(inout) :: self
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
  subroutine IsoVol_Prologue(self)
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    implicit none
    class(IsoVol), intent(inout) :: self
    integer :: nBoxes

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    write(nout,"(1x,A,F15.8)") "(Iso-Volume) Maximum Volume Change: ", self%maxdV

  end subroutine
!=========================================================================
  subroutine IsoVol_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(IsoVol), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Iso-Volume  Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Iso-Volume  Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Iso-Volume Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,F15.8)") "Final Maximum Volume Change: ", self%maxDv
    endif

    if(self%verbose) then
      write(nout, "(1x,A,I15)") "Iso-Volume, Rejections due to overlap:", self%ovlaprej
      write(nout, "(1x,A,I15)") "Iso-Volume, Rejections due to constraint:", self%constrainrej
      write(nout, "(1x,A,I15)") "Iso-Volume, Rejections due to detailed balance:", self%detailedrej
    endif

  end subroutine
!=========================================================================
  subroutine IsoVol_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(IsoVol), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
    real(dp) :: realVal

    call GetXCommand(line, command, 4, lineStat)
    select case( trim(adjustl(command)) )
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

      case("tunemax")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%tunemax = logicVal

      case("updatefreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%maintFreq = intVal

      case("verbose")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%verbose = logicVal

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine

!=========================================================================
end module
!=========================================================================
