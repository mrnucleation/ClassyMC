module Template_ForceField
  use MasterTemplate, only: classyClass
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, public, extends(classyClass) :: forcefield
    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor 
      procedure, pass :: DetailedECalc 
      procedure, pass :: DiffECalc
      procedure, pass :: SinglePair
      procedure, pass :: SinglePair_Approx
      procedure, pass :: ManyBody
!      procedure, pass :: ShiftECalc_Single
!      procedure, pass :: ShiftECalc_Multi
!      procedure, pass :: NewECalc
!      procedure, pass :: OldECalc
!      procedure, pass :: VolECalc
      procedure, pass :: ProcessIO
      procedure, pass :: GetCutOff
  end type

  contains
!=============================================================================+
  subroutine Constructor(self)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(forcefield), intent(inout) :: self


  end subroutine
!=============================================================================+
  subroutine DetailedECalc(self, curbox, E_T, accept)
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    real(dp), intent(inout) :: E_T
    logical, intent(out) :: accept

    accept = .true.
  end subroutine
!============================================================================
  subroutine DiffECalc(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(forcefield), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

!    select type(disp)
!      class is(DisplacementNew)
!         call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)
!
!      class is(Displacement)
!
!        if(disp(1)%newAtom .and. disp(1)%oldAtom) then
!          if(disp(1)%oldAtmIndx == disp(1)%atmIndx) then
    !      if(.false.) then
!            call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)
!            if(.not. accept) then
!              return
!            endif
!          else
!            E_Diff = 0E0_dp
!            call self % NewECalc(curbox, disp, tempList, tempNNei, E_Half, accept)
!            if(.not. accept) then
!              return
!            endif
!
!            E_Diff = E_Diff + E_Half
!            call self % OldECalc(curbox, disp, E_Half)
!            E_Diff = E_Diff + E_Half
!          endif
!
!          return
!        endif
!
!        if(disp(1)%newAtom) then
!          call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)
!          return
!        endif
!
!        if(disp(1)%oldAtom) then
!          call self % OldECalc(curbox, disp, E_Diff)
!          return
!        endif
!      class default
!        write(*,*) "Unknown Perturbation Type."
!    end select
!
!
  end subroutine
!=============================================================================+
  function SinglePair(self, atmtype1, atmtype2, rsq) result(E_Pair)
    implicit none
    class(forcefield), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

  end function
!=============================================================================+
  function SinglePair_Approx(self, atmtype1, atmtype2, rsq) result(E_Pair)
    implicit none
    class(forcefield), intent(in) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

  end function
!=============================================================================+
  function ManyBody(self, curbox, atmtype1, pos1, atmtypes, posN  ) result(E_Atom)
    implicit none
    class(forcefield), intent(in) :: self
    class(simBox), intent(inout) :: curbox
    integer, intent(in) :: atmtype1
    integer, intent(in) :: atmtypes(:)
    real(dp), intent(in) :: pos1(:)
    real(dp), intent(in) :: posN(:,:)
    real(dp) :: E_Atom

  end function
!=============================================================================+
  subroutine ProcessIO(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(forcefield), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
!=============================================================================+
  function GetCutOff(self) result(rCut)
    implicit none
    class(forcefield), intent(inout) :: self
    real(dp) :: rCut

!    write(*,*) self%rCut
    rCut = self%rCut
  end function
!=============================================================================+
end module
!=============================================================================+
