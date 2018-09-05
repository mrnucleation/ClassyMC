!========================================================
module MCMove_PlaneRotation
use CoordinateTypes, only: DisplacementNew
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision

  type, public, extends(MCMove) :: PlaneRotate
!    real(dp) :: atmps = 1E-30_dp
!    real(dp) :: accpt = 0E0_dp
    real(dp) :: max_rot = 0.2E0_dp
    type(DisplacementNew), allocatable :: disp(:)
!    integer, allocatable :: tempNnei(:)
!    integer, allocatable :: tempList(:, :)

    contains
      procedure, pass :: Constructor => PlaneRotate_Constructor
!      procedure, pass :: GeneratePosition => PlaneRotate_GeneratePosition
      procedure, pass :: FullMove => PlaneRotate_FullMove
      procedure, pass :: Maintenance => PlaneRotate_Maintenance
      procedure, pass :: Prologue => PlaneRotate_Prologue
      procedure, pass :: Epilogue => PlaneRotate_Epilogue
  end type
!========================================================
 contains
!========================================================
  subroutine PlaneRotate_Constructor(self)
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(PlaneRotate), intent(inout) :: self
    integer :: iType, maxAtoms

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo


    allocate( self%tempNNei(maxAtoms) )
    allocate( self%tempList(1000, maxAtoms) )
    allocate( self%disp(1:maxAtoms) )
  end subroutine
!========================================================
!  subroutine PlaneRotate_GeneratePosition(self, disp)
!    use RandomGen, only: grnd
!    implicit none
!    class(PlaneRotate), intent(in) :: self
!    type(DisplacementNew), intent(inout) :: disp
!    real(dp) :: dx, dy, dz
!      dx = self % max_rot * (2E0_dp*grnd() - 1E0_dp)
!      dy = self % max_rot * (2E0_dp*grnd() - 1E0_dp)
!      dz = self % max_rot * (2E0_dp*grnd() - 1E0_dp)
!  end subroutine
!===============================================
  subroutine PlaneRotate_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(PlaneRotate), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    integer :: nPlane
    real(dp) :: dx, dy, dz
    real(dp) :: E_Diff, biasE
    real(dp) :: s_term, c_term, cnt, angle
    real(dp) :: x_scale, y_scale, z_scale
    real(dp) :: xcm, ycm, zcm
    real(dp), parameter :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    accept = .true.

    !Propose move
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)


    nAtoms = MolData(molType)%nAtoms
    ! Select plane of rotation and how far to rotate the molecule
    nPlane = floor( 3.0E0_dp *grnd() + 1E0_dp)
    angle = self%max_rot*(2E0_dp*grnd()-1E0_dp)
    c_term = cos(angle)
    s_term = sin(angle)

!     Determine the geometric center which will act as the pivot point for the rotational motion. 
    xcm = 0E0_dp
    ycm = 0E0_dp
    zcm = 0E0_dp
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      xcm = xcm + trialBox%atoms(1, atomIndx) 
      ycm = ycm + trialBox%atoms(2, atomIndx)  
      zcm = zcm + trialBox%atoms(3, atomIndx)  
    enddo
    xcm = xcm/real(nAtoms, dp)
    ycm = ycm/real(nAtoms, dp)
    zcm = zcm/real(nAtoms, dp)



    select case(nPlane)
      case(1) !xy plane
!        Generate a random translational displacement      
        do iAtom = 1, nAtoms
          atomIndx = molStart + iAtom - 1
          self%disp(iAtom)%molType = molType
          self%disp(iAtom)%molIndx = nMove
          self%disp(iAtom)%atmIndx = atomIndx

 
          self%disp(iAtom)%newlist = .false.
          self%disp(iAtom)%listIndex = iAtom

          x_scale = trialBox%atoms(1, atomIndx) - xcm
          y_scale = trialBox%atoms(2, atomIndx) - ycm
          self%disp(iAtom)%x_new = c_term*x_scale - s_term*y_scale + xcm
          self%disp(iAtom)%y_new = s_term*x_scale + c_term*y_scale + ycm
          self%disp(iAtom)%z_new = trialBox%atoms(3, atomIndx) 
        enddo


      case(2) !xz
!        Generate a random translational displacement      
        do iAtom = 1, nAtoms
          atomIndx = molStart + iAtom - 1
          self%disp(iAtom)%molType = molType
          self%disp(iAtom)%molIndx = nMove
          self%disp(iAtom)%atmIndx = atomIndx

 
          self%disp(iAtom)%newlist = .false.
          self%disp(iAtom)%listIndex = iAtom

          x_scale = trialBox%atoms(1, atomIndx) - xcm
          z_scale = trialBox%atoms(3, atomIndx) - zcm
          self%disp(iAtom)%x_new = c_term*x_scale - s_term*z_scale + xcm
          self%disp(iAtom)%y_new = trialBox%atoms(2, atomIndx) 
          self%disp(iAtom)%z_new = s_term*x_scale + c_term*z_scale + zcm
        enddo


      case(3) !yz plane
!        Generate a random translational displacement      
        do iAtom = 1, nAtoms
          atomIndx = molStart + iAtom - 1
          self%disp(iAtom)%molType = molType
          self%disp(iAtom)%molIndx = nMove
          self%disp(iAtom)%atmIndx = atomIndx

 
          self%disp(iAtom)%newlist = .false.
          self%disp(iAtom)%listIndex = iAtom

          y_scale = trialBox%atoms(2, atomIndx) - ycm
          z_scale = trialBox%atoms(3, atomIndx) - zcm
          self%disp(iAtom)%x_new = trialBox%atoms(1, atomIndx) 
          self%disp(iAtom)%y_new = c_term*y_scale - s_term*z_scale + ycm
          self%disp(iAtom)%z_new = s_term*y_scale + c_term*z_scale + zcm
        enddo


      case default
        return
    end select

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
!    call trialbox% EFunc % Method % ShiftECalc_Single(trialBox, self%disp(1:1), E_Diff)
    call trialbox% EFunc % Method % DiffECalc(trialBox, self%disp(1:nAtoms), self%tempList, self%tempNNei, E_Diff, accept)
    if(.not. accept) then
      return
    endif


    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, Prob, self%disp(1:nAtoms))
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine PlaneRotate_Maintenance(self)
    use Constants, only: pi
    implicit none
    class(PlaneRotate), intent(inout) :: self
    real(dp), parameter :: limit = pi
      
    if(self%atmps .lt. 0.5E0_dp) then
      return
    endif

    if(self%GetAcceptRate() .gt. 50E0_dp) then
      if(self%max_rot*1.01E0_dp .lt. limit) then
        self%max_rot = self%max_rot * 1.01E0_dp
      else 
        self%max_rot = limit       
      endif
    else
      self%max_rot = self%max_rot * 0.99E0_dp
    endif

 

  end subroutine
!=========================================================================
  subroutine PlaneRotate_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(PlaneRotate), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
      


  end subroutine
!=========================================================================
  subroutine PlaneRotate_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(PlaneRotate), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,"(1x,A,I15)") "Molecule Rotation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "Molecule Rotation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "Molecule Rotation Acceptance Rate: ", accptRate
    write(nout,"(1x,A,F15.8)") "Final Maximum Displacement: ", self%max_rot
 

  end subroutine
!========================================================
end module
!========================================================
