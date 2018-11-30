!==========================================================================================
! The purpose of this module is to provide the base class for the molecule
! constructor class.
!==========================================================================================
#define _IOINDEX 30
!==========================================================================================
module MolCon_RidgidRegrowth
  use CoordinateTypes, only: Perturbation, Addition
  use Template_SimBox, only: SimBox
  use Template_MolConstructor, only: MolConstructor
  use VarPrecision


  type, public, extends(MolConstructor) :: RidgidRegrowth
    integer :: nAtoms
    real(dp), allocatable :: tempcoords(:, :)
    real(dp), allocatable :: ridgidconfig(:,:)
    character(len=60) :: filename
    contains
!      procedure, public, pass :: Constructor => RidgidRegrowth_Constructor
      procedure, public, pass :: Prologue => RidgidRegrowth_Prologue
      procedure, public, pass :: GenerateConfig => RidgidRegrowth_GenerateConfig
      procedure, public, pass :: ProcessIO => RidgidRegrowth_ProcessIO
  end type
!==========================================================================================
  contains
!==========================================================================================
  subroutine RidgidRegrowth_Prologue(self)
    use Common_MolInfo, only: MolData, nMolTypes
    use ParallelVar, only: nout
    implicit none
    class(RidgidRegrowth), intent(inout) :: self
    integer :: j, iType, nAtoms, iLine, nLines
    integer :: iostatus
    character(len=5) :: symbol

    nAtoms = MolData(self%molType)%nAtoms
    if(nAtoms < 1) then
      write(0,*) "ERROR! Regrowth Type defined before the molecular topology is defined!"
      stop
    endif
    self%nAtoms = nAtoms
    allocate( self%tempcoords(1:3, 1:nAtoms) )
    allocate( self%ridgidconfig(1:3, 1:nAtoms) )

    do j = 1, len(self%filename)
      if(self%filename(j:j) == '"') then
        self%filename(j:j) = " "
      endif
    enddo

    open(unit=_IOINDEX, file=self%filename, status='old')

    nLines = 0
    do
      read(_IOINDEX, *, ioStat=iostatus)
      if(iostatus /= 0) then
        exit
      endif
      nLines = nLines + 1
    enddo
    rewind(_IOINDEX)

    if(nLines /= nAtoms+2) then
      write(0,"(A,A)") "File Name: ", trim(adjustl(self%filename))
      write(0,*) "INPUT ERROR! The an invalid configuration file has been given"
      write(0,*) "Expected number of lines:", nAtoms+2
      write(0,*) "Received number of lines:", nLines
      stop
    endif

    read(_IOINDEX, *)
    read(_IOINDEX, *)
    do iLine = 1, nLines-2
      read(_IOINDEX, *) symbol, (self%ridgidconfig(j,iLine), j=1,3)
    enddo

    write(nout,"(A,A)") trim(adjustl(self%filename)), " sucessfully read!"
    close(_IOINDEX)



  end subroutine
!==========================================================================================
  subroutine RidgidRegrowth_GenerateConfig(self, trialBox, disp, probconstruct, insPoint)
    use ClassyConstants, only: two_pi
    use Common_MolInfo, only: MolData
    use RandomGen, only: grnd
    implicit none
    class(RidgidRegrowth), intent(inout) :: self
    class(Perturbation), intent(inout) :: disp(:)
    class(SimBox), intent(inout) :: trialBox
    real(dp), intent(in), optional :: insPoint(:)

    integer :: iDisp, iAtom
    real(dp), intent(out) :: probconstruct
    real(dp) :: c_term, s_term, rotang
    real(dp) :: x_shift, y_shift, z_shift
    real(dp) :: x1, y1, z1
    real(dp) :: x_rid_cm, y_rid_cm, z_rid_cm

    probconstruct = 1E0
    self%tempcoords(1:3,1:self%nAtoms) = self%ridgidconfig(1:3, 1:self%nAtoms)

    x1 = self%tempcoords(1, 1)
    y1 = self%tempcoords(2, 1)
    z1 = self%tempcoords(3, 1)
    do iAtom = 1, self%nAtoms
      self%tempcoords(1, iAtom) = self%tempcoords(1, iAtom) - x1
      self%tempcoords(2, iAtom) = self%tempcoords(2, iAtom) - y1
      self%tempcoords(3, iAtom) = self%tempcoords(3, iAtom) - z1
    enddo
    if(self%nAtoms /= 1) then
      x_rid_cm = 0E0_dp
      y_rid_cm = 0E0_dp
      z_rid_cm = 0E0_dp
    
!          Rotate xz axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do iAtom = 1,self%nAtoms
        x_shift = self%tempcoords(1, iAtom) - x_rid_cm
        z_shift = self%tempcoords(3, iAtom) - z_rid_cm        
        self%tempcoords(1, iAtom) = c_term*x_shift - s_term*z_shift + x_rid_cm
        self%tempcoords(3, iAtom) = s_term*x_shift + c_term*z_shift + z_rid_cm
      enddo   
    
!          Rotate yz axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do iAtom = 1, self%nAtoms
        y_shift = self%tempcoords(2, iAtom) - y_rid_cm
        z_shift = self%tempcoords(3, iAtom) - z_rid_cm        
        self%tempcoords(2, iAtom) = c_term*y_shift - s_term*z_shift + y_rid_cm
        self%tempcoords(3, iAtom) = s_term*y_shift + c_term*z_shift + z_rid_cm
      enddo   

!          Rotate xy axis
      rotang = two_pi*grnd()
      c_term = cos(rotang)
      s_term = sin(rotang)
      do iAtom = 1, self%nAtoms
        x_shift = self%tempcoords(1, iAtom) - x_rid_cm
        y_shift = self%tempcoords(2, iAtom) - y_rid_cm        
        self%tempcoords(1, iAtom) = c_term*x_shift - s_term*y_shift + x_rid_cm
        self%tempcoords(2, iAtom) = s_term*x_shift + c_term*y_shift + y_rid_cm
      enddo    
    endif

    select type(disp)
      class is(Addition)
        do iDisp = 1, self%nAtoms
          if( present(insPoint) ) then
            disp(iDisp)%x_new = self%tempcoords(1, iDisp) + insPoint(1)
            disp(iDisp)%y_new = self%tempcoords(2, iDisp) + insPoint(2)
            disp(iDisp)%z_new = self%tempcoords(3, iDisp) + insPoint(3)
          else
            disp(iDisp)%x_new = self%tempcoords(1, iDisp)
            disp(iDisp)%y_new = self%tempcoords(2, iDisp)
            disp(iDisp)%z_new = self%tempcoords(3, iDisp)
          endif
        enddo
    end select


  end subroutine
!==========================================================================================
  subroutine RidgidRegrowth_ProcessIO(self, line, linestat)
    use Input_Format, only: maxLineLen, GetXCommand
    implicit none
    class(RidgidRegrowth), intent(inout) :: self
    character(len=*), intent(in) :: line
    integer, intent(out) :: linestat

    call GetXCommand(line, self%filename, 3, lineStat)
    if(lineStat /= 0) then
      return  
    endif   


  end subroutine
!==========================================================================================
end module
!==========================================================================================
