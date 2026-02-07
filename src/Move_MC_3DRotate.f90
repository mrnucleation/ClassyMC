!========================================================
module MCMove_3DRotation
use CoordinateTypes, only: Displacement
use MoveClassDef
use SimpleSimBox, only: SimpleBox
use VarPrecision
use ClassyConstants, only: pi

  type, public, extends(MCMove) :: Rotate3D
    logical :: proportional = .false.
    logical :: tuneMax = .true.
    real(dp) :: targAccpt = 50E0_dp
    real(dp) :: max_rot = 0.02E0_dp
    real(dp) :: limit = pi
    type(Displacement), allocatable :: disp(:)

    contains
      procedure, pass :: Constructor => Rotate3D_Constructor
      procedure, pass :: FullMove => Rotate3D_FullMove
      procedure, pass :: Maintenance => Rotate3D_Maintenance
      procedure, pass :: Prologue => Rotate3D_Prologue
      procedure, pass :: Epilogue => Rotate3D_Epilogue
      procedure, pass :: ProcessIO => Rotate3D_ProcessIO
  end type
!========================================================
 contains
!========================================================
  subroutine Rotate3D_Constructor(self)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    class(Rotate3D), intent(inout) :: self
    integer :: iType, maxAtoms, nBoxes

    if(.not. allocated(self%boxProb)) then
      nBoxes = size(boxArray)
      allocate( self%boxProb(1:nBoxes) )
      self%boxProb = 1E0_dp/real(nBoxes,dp)
    endif

    maxAtoms = 0
    do iType = 1, nMolTypes
      if(MolData(iType)%nAtoms > maxAtoms) then
        maxAtoms = MolData(iType)%nAtoms 
      endif
    enddo

    allocate( self%disp(1:maxAtoms) )
    call self%CreateTempArray(maxAtoms)
  end subroutine
!========================================================
  subroutine Rotate3D_FullMove(self, trialBox, accept) 
    use Box_Utility, only: FindAtom, FindMolecule
    use CommonSampling, only: sampling
    use Common_NeighData, only: neighSkin
    use Common_MolInfo, only: MolData, nMolTypes
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: grnd
    implicit none
    class(Rotate3D), intent(inout) :: self
    class(SimpleBox), intent(inout) :: trialBox
    logical, intent(out) :: accept
    integer :: iAtom, nAtoms, atomIndx
    integer :: nMove, rawIndx, iConstrain
    integer :: CalcIndex, molStart, molEnd, molType
    real(dp) :: E_Diff, E_Inter, E_Intra, biasE
    real(dp) :: angle, c_a, s_a, omc
    real(dp) :: ux, uy, uz, theta, phi
    real(dp) :: xcm, ycm, zcm
    real(dp) :: rx, ry, rz
    real(dp) :: rx_new, ry_new, rz_new
    real(dp) :: R11, R12, R13
    real(dp) :: R21, R22, R23
    real(dp) :: R31, R32, R33
    real(dp), parameter :: Prob = 1E0_dp

    self % atmps = self % atmps + 1E0_dp
    accept = .true.
    call self%LoadBoxInfo(trialbox, self%disp)

    !Propose move - pick a random molecule
    rawIndx = floor( trialBox%nMolTotal * grnd() + 1E0_dp)
    call FindMolecule(trialbox, rawIndx, nMove)
    call trialBox % GetMolData(nMove, molStart=molStart, molEnd=molEnd, &
                               molType=molType)

    nAtoms = MolData(molType)%nAtoms

    ! Generate a random rotation axis uniformly on the unit sphere
    theta = acos(2E0_dp*grnd() - 1E0_dp)
    phi   = 2E0_dp * pi * grnd()
    ux = sin(theta) * cos(phi)
    uy = sin(theta) * sin(phi)
    uz = cos(theta)

    ! Generate a random rotation angle
    angle = self%max_rot * (2E0_dp*grnd() - 1E0_dp)
    c_a = cos(angle)
    s_a = sin(angle)
    omc = 1E0_dp - c_a

    ! Build rotation matrix (Rodrigues' rotation formula)
    R11 = c_a + ux*ux*omc
    R12 = ux*uy*omc - uz*s_a
    R13 = ux*uz*omc + uy*s_a

    R21 = uy*ux*omc + uz*s_a
    R22 = c_a + uy*uy*omc
    R23 = uy*uz*omc - ux*s_a

    R31 = uz*ux*omc - uy*s_a
    R32 = uz*uy*omc + ux*s_a
    R33 = c_a + uz*uz*omc

    ! Use the first atom as the pivot point
    atomIndx = molStart
    xcm = trialBox%atoms(1, atomIndx) 
    ycm = trialBox%atoms(2, atomIndx)  
    zcm = trialBox%atoms(3, atomIndx)  

    ! Apply the rotation to each atom in the molecule
    do iAtom = 1, nAtoms
      atomIndx = molStart + iAtom - 1
      self%disp(iAtom)%molType = molType
      self%disp(iAtom)%molIndx = nMove
      self%disp(iAtom)%atmIndx = atomIndx
      self%disp(iAtom)%newlist = .false.
      self%disp(iAtom)%listIndex = iAtom

      rx = trialBox%atoms(1, atomIndx) - xcm
      ry = trialBox%atoms(2, atomIndx) - ycm
      rz = trialBox%atoms(3, atomIndx) - zcm
      call trialBox % Boundary(rx, ry, rz)

      ! Apply rotation matrix: r_new = R * r
      rx_new = R11*rx + R12*ry + R13*rz
      ry_new = R21*rx + R22*ry + R23*rz
      rz_new = R31*rx + R32*ry + R33*rz

      self%disp(iAtom)%x_new = xcm + rx_new
      self%disp(iAtom)%y_new = ycm + ry_new
      self%disp(iAtom)%z_new = zcm + rz_new
    enddo

    !Check Constraint
    accept = trialBox % CheckConstraint( self%disp(1:nAtoms) )
    if(.not. accept) then
      return
    endif

    !Energy Calculation
    call trialBox%ComputeEnergyDelta(self%disp(1:nAtoms),&
                                     self%templist,&
                                     self%tempNNei, &
                                     E_Inter, &
                                     E_Intra, &
                                     E_Diff, &
                                     accept, &
                                     computeintra=.false.)
    if(.not. accept) then
      return
    endif

    !Check Post Energy Constraint
    accept = trialBox % CheckPostEnergy( self%disp(1:nAtoms), E_Diff )
    if(.not. accept) then
      return
    endif

    !Accept/Reject
    accept = sampling % MakeDecision(trialBox, E_Diff, self%disp(1:nAtoms), inProb=Prob)
    if(accept) then
      self % accpt = self % accpt + 1E0_dp
      call trialBox % UpdateEnergy(E_Diff, E_Inter)
      call trialBox % UpdatePosition(self%disp(1:nAtoms), self%tempList, self%tempNNei)
    endif

  end subroutine
!=========================================================================
  subroutine Rotate3D_Maintenance(self)
    use ClassyConstants, only: pi
    implicit none
    class(Rotate3D), intent(inout) :: self
      
    if(self%tunemax) then
      if(self%atmps < 0.5E0_dp) then
        return
      endif

      if(self%GetAcceptRate() > self%targAccpt) then
        if(self%max_rot*1.01E0_dp < self%limit) then
          self%max_rot = self%max_rot * 1.01E0_dp
        else 
          self%max_rot = self%limit
        endif
      else
        self%max_rot = self%max_rot * 0.99E0_dp
      endif
    endif
  end subroutine
!=========================================================================
  subroutine Rotate3D_Prologue(self)
    use ParallelVar, only: nout
    implicit none
    class(Rotate3D), intent(inout) :: self

    if(.not. allocated(self%disp)) then
      call self % Constructor
    endif
    write(nout,"(1x,A,F15.8)") "3D Rotate: Maximum Rotation Angle: ", self%max_rot
  end subroutine
!=========================================================================
  subroutine Rotate3D_Epilogue(self)
    use ParallelVar, only: nout
    implicit none
    class(Rotate3D), intent(inout) :: self
    real(dp) :: accptRate
      
    write(nout,*) 
    write(nout,"(1x,A,I15)") "3D Rotation Moves Accepted: ", nint(self%accpt)
    write(nout,"(1x,A,I15)") "3D Rotation Moves Attempted: ", nint(self%atmps)
    accptRate = self%GetAcceptRate()
    write(nout,"(1x,A,F15.8)") "3D Rotation Acceptance Rate: ", accptRate
    if(self%tunemax) then
      write(nout,"(1x,A,F15.8)") "3D Rotate: Final Maximum Angle: ", self%max_rot
    endif
 
  end subroutine
!=========================================================================
  subroutine Rotate3D_ProcessIO(self, line, lineStat)
    use Input_Format, only: GetXCommand, maxLineLen
    implicit none
    class(Rotate3D), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat
    character(len=30) :: command
    logical :: logicVal
    integer :: intVal
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
        if(realVal > pi) then
          self%limit = realVal
        else
          self%limit = pi
        endif

      case("dynamictarget")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        self%targAccpt = realVal

      case("maxangle")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) realVal
        if(realVal > pi) then
          self%max_rot = realVal
        else
          self%max_rot = pi
        endif

      case("proportional")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) logicVal
        self%proportional = logicVal

      case("updatefreq")
        call GetXCommand(line, command, 5, lineStat)
        read(command, *) intVal
        self%maintFreq = intVal

      case default
        lineStat = -1
        return

    end select
    lineStat = 0

  end subroutine
!========================================================
end module
!========================================================
