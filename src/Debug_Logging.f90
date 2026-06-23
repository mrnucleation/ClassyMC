!========================================================================================
! Debug Logging Module for ClassyMC
! 
! This module provides debugging/visualization tools for Monte Carlo moves and constraints.
! All functionality is conditionally compiled with the LOGGING preprocessor flag.
!
! COMPILATION:
!   make LOGGING=1              # Enable debug logging
!   make LOGGING=1 DEBUG=1      # Enable both debug logging and debug flags
!
! OUTPUT FILES:
!   - debug_moves.log           # Human-readable log of all move details
!   - debug_trajectory.lammpstrj # LAMMPS dump format for Ovito visualization
!
! KEY SUBROUTINES:
!   - DebugLog_Initialize(logfile, dumpfile) - Initialize logging system
!   - DebugLog_Finalize()                    - Close all files
!   - DebugLog_DumpBox(box, label)           - Dump current box state
!   - DebugLog_DumpProposal(box, disp, label) - Dump proposed coordinates (Displacement)
!   - DebugLog_DumpProposalAddition(box, disp, label) - Dump proposed coordinates (Addition)
!   - DebugLog_MoveResult(movename, accepted, E_Diff, ...) - Log move acceptance
!   - DebugLog_ConstraintCheck(name, passed, details) - Log constraint checks
!   - DebugLog_Message(message, level)       - Write custom log message
!   - DebugLog_IsEnabled()                   - Check if logging is active
!
! USAGE EXAMPLE IN MOVE FILES:
!
!   use Debug_Logging, only: DebugLog_DumpBox, DebugLog_DumpProposal, DebugLog_MoveResult
!
!   subroutine MyMove_FullMove(self, trialBox, accept)
!     ...
!     ! After generating proposal, before energy calculation:
!     call DebugLog_DumpProposal(trialBox, self%disp(1:nAtoms), "MolTranslate Pre-Energy")
!
!     ! After energy calculation, compute delta:
!     call trialBox%ComputeEnergyDelta(...)
!
!     ! After accept/reject decision:
!     call DebugLog_MoveResult("MolTranslate", accept, E_Diff=E_Diff, E_Inter=E_Inter)
!
!     ! Optionally dump post-move state if accepted:
!     if(accept) then
!       call DebugLog_DumpBox(trialBox, "MolTranslate Post-Accept")
!     endif
!     ...
!   end subroutine
!
! OVITO VISUALIZATION:
!   The dump file uses a custom column "c_proposed" to distinguish:
!     - c_proposed = 0: Current/existing atoms
!     - c_proposed = 1: Proposed new positions (shown with type + 100)
!   
!   In Ovito, color by "c_proposed" or "Particle Type" to see proposals vs current.
!
! LOG FILE FORMAT:
!   Each move logs:
!   - Frame number (matches dump file timestep)
!   - Current and proposed coordinates
!   - Energy differences (E_Diff, E_Inter, E_Intra)
!   - Acceptance decision
!
!========================================================================================
module Debug_Logging
  use VarPrecision
  use CoordinateTypes
  implicit none

#ifdef LOGGING
  private
  public :: DebugLog_Initialize, DebugLog_Finalize
  public :: DebugLog_DumpBox, DebugLog_DumpProposal, DebugLog_DumpProposalAddition
  public :: DebugLog_DumpProposalDeletion, DebugLog_DumpVolumeChange
  public :: DebugLog_MoveResult, DebugLog_ConstraintCheck
  public :: DebugLog_Message, DebugLog_Separator
  public :: DebugLog_IsEnabled

  integer, parameter :: LOG_UNIT = 991
  integer, parameter :: DUMP_UNIT = 992
  
  logical :: logging_initialized = .false.
  integer(kind=8) :: dump_frame_counter = 0
  integer(kind=8) :: log_entry_counter = 0
  character(len=256) :: current_logfile = ""
  character(len=256) :: current_dumpfile = ""

contains

!========================================================================================
! Check if logging is enabled (always true when compiled with LOGGING)
!========================================================================================
function DebugLog_IsEnabled() result(enabled)
  logical :: enabled
  enabled = .true.
end function

!========================================================================================
! Initialize the logging system
!========================================================================================
subroutine DebugLog_Initialize(logfile, dumpfile)
  character(len=*), intent(in), optional :: logfile
  character(len=*), intent(in), optional :: dumpfile
  
  if(logging_initialized) then
    call DebugLog_Finalize()
  endif
  
  if(present(logfile)) then
    current_logfile = trim(adjustl(logfile))
  else
    current_logfile = "debug_moves.log"
  endif
  
  if(present(dumpfile)) then
    current_dumpfile = trim(adjustl(dumpfile))
  else
    current_dumpfile = "debug_trajectory.lammpstrj"
  endif
  
  open(unit=LOG_UNIT, file=trim(current_logfile), status='replace', action='write')
  open(unit=DUMP_UNIT, file=trim(current_dumpfile), status='replace', action='write')
  
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, '(A)') "  ClassyMC Debug Log"
  write(LOG_UNIT, '(A)') "  Logging Mode: ENABLED"
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, *)
  
  dump_frame_counter = 0
  log_entry_counter = 0
  logging_initialized = .true.
  
end subroutine

!========================================================================================
! Finalize and close all logging files
!========================================================================================
subroutine DebugLog_Finalize()
  
  if(.not. logging_initialized) return
  
  write(LOG_UNIT, *)
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, '(A,I12)') "  Total Log Entries: ", log_entry_counter
  write(LOG_UNIT, '(A,I12)') "  Total Dump Frames: ", dump_frame_counter
  write(LOG_UNIT, '(A)') "========================================================"
  
  close(LOG_UNIT)
  close(DUMP_UNIT)
  
  logging_initialized = .false.
  
end subroutine

!========================================================================================
! Write a separator line to the log
!========================================================================================
subroutine DebugLog_Separator(char)
  character(len=1), intent(in), optional :: char
  character(len=1) :: sep_char
  integer :: i
  
  if(.not. logging_initialized) return
  
  if(present(char)) then
    sep_char = char
  else
    sep_char = "-"
  endif
  
  write(LOG_UNIT, '(80A1)') (sep_char, i=1,80)
  
end subroutine

!========================================================================================
! Write a custom message to the log
!========================================================================================
subroutine DebugLog_Message(message, level)
  character(len=*), intent(in) :: message
  integer, intent(in), optional :: level
  integer :: msg_level
  character(len=12) :: level_str
  
  if(.not. logging_initialized) return
  
  if(present(level)) then
    msg_level = level
  else
    msg_level = 1
  endif
  
  select case(msg_level)
    case(0)
      level_str = "[DEBUG]"
    case(1)
      level_str = "[INFO]"
    case(2)
      level_str = "[WARNING]"
    case(3)
      level_str = "[ERROR]"
    case default
      level_str = "[INFO]"
  end select
  
  log_entry_counter = log_entry_counter + 1
  write(LOG_UNIT, '(A12,1X,A)') trim(level_str), trim(message)
  
end subroutine

!========================================================================================
! Dump current box state to LAMMPS dump file and log file
!========================================================================================
subroutine DebugLog_DumpBox(box, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  use Common_MolInfo, only: nMolTypes, AtomData
  
  class(SimpleBox), intent(inout) :: box
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
  
  integer :: iAtom, iType, typeStart, typeEnd
  integer :: atomType, molType, nDim
  real(dp) :: localBoxDim(2,3)
  
  if(.not. logging_initialized) return
  
  nDim = box%nDimension
  dump_frame_counter = dump_frame_counter + 1
  log_entry_counter = log_entry_counter + 1
  
  if(present(boxDim)) then
    localBoxDim = boxDim
  else
    localBoxDim(1,:) = -50.0_dp
    localBoxDim(2,:) = 50.0_dp
    call GetBoxDimensions(box, localBoxDim)
  endif
  
  write(DUMP_UNIT, '(A)') "ITEM: TIMESTEP"
  write(DUMP_UNIT, '(I12)') dump_frame_counter
  
  write(DUMP_UNIT, '(A)') "ITEM: NUMBER OF ATOMS"
  write(DUMP_UNIT, '(I12)') box%nAtoms
  
  write(DUMP_UNIT, '(A)') "ITEM: BOX BOUNDS pp pp pp"
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,1), localBoxDim(2,1)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,2), localBoxDim(2,2)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,3), localBoxDim(2,3)
  
  write(DUMP_UNIT, '(A)') "ITEM: ATOMS id type mol x y z"
  
  do iType = 1, nMolTypes
    call box%GetTypeAtoms(iType, typeStart, typeEnd)
    if(typeStart < 1) cycle
    do iAtom = typeStart, typeEnd
      atomType = box%AtomType(iAtom)
      molType = box%MolType(iAtom)
      write(DUMP_UNIT, '(I8,2I4,3E20.10)') iAtom, atomType, box%MolIndx(iAtom), &
                                           box%atoms(1,iAtom), box%atoms(2,iAtom), box%atoms(3,iAtom)
    enddo
  enddo
  
  write(LOG_UNIT, '(A)') "--------------------------------------------------------"
  write(LOG_UNIT, '(A,A)') "DUMP FRAME: ", trim(label)
  write(LOG_UNIT, '(A,I12)') "  Frame Number: ", dump_frame_counter
  write(LOG_UNIT, '(A,I8)') "  Box ID: ", box%boxID
  write(LOG_UNIT, '(A,I8)') "  Number of Atoms: ", box%nAtoms
  write(LOG_UNIT, '(A,I8)') "  Number of Molecules: ", box%nMolTotal
  write(LOG_UNIT, '(A,E20.10)') "  Total Energy: ", box%ETotal
  write(LOG_UNIT, '(A,E20.10)') "  Volume: ", box%volume
  write(LOG_UNIT, '(A)') "--------------------------------------------------------"
  
end subroutine

!========================================================================================
! Helper to get box dimensions from different box types
!========================================================================================
subroutine GetBoxDimensions(box, boxDim)
  use SimpleSimBox, only: SimpleBox
  use OrthoBoxDef, only: OrthoBox
  use CubicBoxDef, only: CubeBox
  
  class(SimpleBox), intent(in) :: box
  real(dp), intent(inout) :: boxDim(2,3)
  
  select type(box)
    class is(OrthoBox)
      boxDim(1,1) = -box%boxLx2
      boxDim(2,1) = box%boxLx2
      boxDim(1,2) = -box%boxLy2
      boxDim(2,2) = box%boxLy2
      boxDim(1,3) = -box%boxLz2
      boxDim(2,3) = box%boxLz2
    class is(CubeBox)
      boxDim(1,1) = -box%boxL2
      boxDim(2,1) = box%boxL2
      boxDim(1,2) = -box%boxL2
      boxDim(2,2) = box%boxL2
      boxDim(1,3) = -box%boxL2
      boxDim(2,3) = box%boxL2
    class default
      boxDim(1,:) = -50.0_dp
      boxDim(2,:) = 50.0_dp
  end select
  
end subroutine

!========================================================================================
! Dump proposed coordinates (Displacement moves)
! Shows current state plus proposed positions as separate "atoms" for visualization
!========================================================================================
subroutine DebugLog_DumpProposal(box, disp, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  use Common_MolInfo, only: nMolTypes, AtomData
  
  class(SimpleBox), intent(inout) :: box
  type(Displacement), intent(in) :: disp(:)
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
  
  integer :: iAtom, iDisp, iType, typeStart, typeEnd
  integer :: atomType, molType, nDim, nDisp
  integer :: totalAtoms
  real(dp) :: localBoxDim(2,3)
  
  if(.not. logging_initialized) return
  
  nDim = box%nDimension
  nDisp = size(disp)
  totalAtoms = box%nAtoms + nDisp
  
  dump_frame_counter = dump_frame_counter + 1
  log_entry_counter = log_entry_counter + 1
  
  if(present(boxDim)) then
    localBoxDim = boxDim
  else
    localBoxDim(1,:) = -50.0_dp
    localBoxDim(2,:) = 50.0_dp
    call GetBoxDimensions(box, localBoxDim)
  endif
  
  write(DUMP_UNIT, '(A)') "ITEM: TIMESTEP"
  write(DUMP_UNIT, '(I12)') dump_frame_counter
  
  write(DUMP_UNIT, '(A)') "ITEM: NUMBER OF ATOMS"
  write(DUMP_UNIT, '(I12)') totalAtoms
  
  write(DUMP_UNIT, '(A)') "ITEM: BOX BOUNDS pp pp pp"
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,1), localBoxDim(2,1)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,2), localBoxDim(2,2)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,3), localBoxDim(2,3)
  
  write(DUMP_UNIT, '(A)') "ITEM: ATOMS id type mol x y z c_proposed"
  
  do iType = 1, nMolTypes
    call box%GetTypeAtoms(iType, typeStart, typeEnd)
    if(typeStart < 1) cycle
    do iAtom = typeStart, typeEnd
      atomType = box%AtomType(iAtom)
      write(DUMP_UNIT, '(I8,2I4,3E20.10,I4)') iAtom, atomType, box%MolIndx(iAtom), &
                                              box%atoms(1,iAtom), box%atoms(2,iAtom), box%atoms(3,iAtom), 0
    enddo
  enddo
  
  do iDisp = 1, nDisp
    atomType = box%AtomType(disp(iDisp)%atmIndx)
    write(DUMP_UNIT, '(I8,2I4,3E20.10,I4)') box%nAtoms + iDisp, atomType + 100, &
                                            disp(iDisp)%molIndx, &
                                            disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new, 1
  enddo
  
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, '(A,A)') "PROPOSAL: ", trim(label)
  write(LOG_UNIT, '(A,I12)') "  Frame Number: ", dump_frame_counter
  write(LOG_UNIT, '(A,I8)') "  Box ID: ", box%boxID
  write(LOG_UNIT, '(A,I4)') "  Number of Displaced Atoms: ", nDisp
  write(LOG_UNIT, '(A)') "  ---"
  write(LOG_UNIT, '(A)') "  Displacement Details:"
  
  do iDisp = 1, nDisp
    atomType = box%AtomType(disp(iDisp)%atmIndx)
    write(LOG_UNIT, '(A,I4)') "    Atom ", iDisp
    write(LOG_UNIT, '(A,I8)') "      Global Index: ", disp(iDisp)%atmIndx
    write(LOG_UNIT, '(A,I4)') "      Mol Type: ", disp(iDisp)%molType
    write(LOG_UNIT, '(A,I8)') "      Mol Index: ", disp(iDisp)%molIndx
    write(LOG_UNIT, '(A,I4)') "      Atom Type: ", atomType
    write(LOG_UNIT, '(A,3E15.6)') "      Current:  ", box%atoms(1,disp(iDisp)%atmIndx), &
                                                      box%atoms(2,disp(iDisp)%atmIndx), &
                                                      box%atoms(3,disp(iDisp)%atmIndx)
    write(LOG_UNIT, '(A,3E15.6)') "      Proposed: ", disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new
    write(LOG_UNIT, '(A,3E15.6)') "      Delta:    ", disp(iDisp)%x_new - box%atoms(1,disp(iDisp)%atmIndx), &
                                                      disp(iDisp)%y_new - box%atoms(2,disp(iDisp)%atmIndx), &
                                                      disp(iDisp)%z_new - box%atoms(3,disp(iDisp)%atmIndx)
  enddo
  write(LOG_UNIT, '(A)') "========================================================"
  
end subroutine

!========================================================================================
! Dump proposed coordinates (Addition moves)
!========================================================================================
subroutine DebugLog_DumpProposalAddition(box, disp, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  use Common_MolInfo, only: nMolTypes, AtomData
  
  class(SimpleBox), intent(inout) :: box
  type(Addition), intent(in) :: disp(:)
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
  
  integer :: iAtom, iDisp, iType, typeStart, typeEnd
  integer :: atomType, molType, nDim, nDisp
  integer :: totalAtoms
  real(dp) :: localBoxDim(2,3)
  
  if(.not. logging_initialized) return
  
  nDim = box%nDimension
  nDisp = size(disp)
  totalAtoms = box%nAtoms + nDisp
  
  dump_frame_counter = dump_frame_counter + 1
  log_entry_counter = log_entry_counter + 1
  
  if(present(boxDim)) then
    localBoxDim = boxDim
  else
    localBoxDim(1,:) = -50.0_dp
    localBoxDim(2,:) = 50.0_dp
    call GetBoxDimensions(box, localBoxDim)
  endif
  
  write(DUMP_UNIT, '(A)') "ITEM: TIMESTEP"
  write(DUMP_UNIT, '(I12)') dump_frame_counter
  
  write(DUMP_UNIT, '(A)') "ITEM: NUMBER OF ATOMS"
  write(DUMP_UNIT, '(I12)') totalAtoms
  
  write(DUMP_UNIT, '(A)') "ITEM: BOX BOUNDS pp pp pp"
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,1), localBoxDim(2,1)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,2), localBoxDim(2,2)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,3), localBoxDim(2,3)
  
  write(DUMP_UNIT, '(A)') "ITEM: ATOMS id type mol x y z c_proposed"
  
  do iType = 1, nMolTypes
    call box%GetTypeAtoms(iType, typeStart, typeEnd)
    if(typeStart < 1) cycle
    do iAtom = typeStart, typeEnd
      atomType = box%AtomType(iAtom)
      write(DUMP_UNIT, '(I8,2I4,3E20.10,I4)') iAtom, atomType, box%MolIndx(iAtom), &
                                              box%atoms(1,iAtom), box%atoms(2,iAtom), box%atoms(3,iAtom), 0
    enddo
  enddo
  
  do iDisp = 1, nDisp
    atomType = box%AtomType(disp(iDisp)%atmIndx)
    write(DUMP_UNIT, '(I8,2I4,3E20.10,I4)') box%nAtoms + iDisp, atomType + 100, &
                                            disp(iDisp)%molIndx, &
                                            disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new, 1
  enddo
  
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, '(A,A)') "ADDITION PROPOSAL: ", trim(label)
  write(LOG_UNIT, '(A,I12)') "  Frame Number: ", dump_frame_counter
  write(LOG_UNIT, '(A,I8)') "  Box ID: ", box%boxID
  write(LOG_UNIT, '(A,I4)') "  Number of Added Atoms: ", nDisp
  write(LOG_UNIT, '(A)') "  ---"
  write(LOG_UNIT, '(A)') "  Addition Details:"
  
  do iDisp = 1, nDisp
    atomType = box%AtomType(disp(iDisp)%atmIndx)
    write(LOG_UNIT, '(A,I4)') "    Atom ", iDisp
    write(LOG_UNIT, '(A,I8)') "      Target Index: ", disp(iDisp)%atmIndx
    write(LOG_UNIT, '(A,I4)') "      Mol Type: ", disp(iDisp)%molType
    write(LOG_UNIT, '(A,I8)') "      Mol Index: ", disp(iDisp)%molIndx
    write(LOG_UNIT, '(A,I4)') "      Atom Type: ", atomType
    write(LOG_UNIT, '(A,3E15.6)') "      Position: ", disp(iDisp)%x_new, disp(iDisp)%y_new, disp(iDisp)%z_new
  enddo
  write(LOG_UNIT, '(A)') "========================================================"
  
end subroutine

!========================================================================================
! Log move acceptance/rejection result
!========================================================================================
subroutine DebugLog_MoveResult(moveName, accepted, E_Diff, E_Inter, E_Intra, biasProb)
  character(len=*), intent(in) :: moveName
  logical, intent(in) :: accepted
  real(dp), intent(in), optional :: E_Diff
  real(dp), intent(in), optional :: E_Inter
  real(dp), intent(in), optional :: E_Intra
  real(dp), intent(in), optional :: biasProb
  
  character(len=12) :: result_str
  
  if(.not. logging_initialized) return
  
  log_entry_counter = log_entry_counter + 1
  
  if(accepted) then
    result_str = "ACCEPTED"
  else
    result_str = "REJECTED"
  endif
  
  write(LOG_UNIT, '(A)') "--------------------------------------------------------"
  write(LOG_UNIT, '(A,A,A,A)') "MOVE RESULT: ", trim(moveName), " -> ", result_str
  
  if(present(E_Diff)) then
    write(LOG_UNIT, '(A,E20.10)') "  Energy Difference: ", E_Diff
  endif
  
  if(present(E_Inter)) then
    write(LOG_UNIT, '(A,E20.10)') "  Inter Energy:      ", E_Inter
  endif
  
  if(present(E_Intra)) then
    write(LOG_UNIT, '(A,E20.10)') "  Intra Energy:      ", E_Intra
  endif
  
  if(present(biasProb)) then
    write(LOG_UNIT, '(A,E20.10)') "  Bias Probability:  ", biasProb
  endif
  
  write(LOG_UNIT, '(A)') "--------------------------------------------------------"
  
end subroutine

!========================================================================================
! Log constraint check results
!========================================================================================
subroutine DebugLog_ConstraintCheck(constraintName, passed, details)
  character(len=*), intent(in) :: constraintName
  logical, intent(in) :: passed
  character(len=*), intent(in), optional :: details
  
  character(len=12) :: result_str
  
  if(.not. logging_initialized) return
  
  log_entry_counter = log_entry_counter + 1
  
  if(passed) then
    result_str = "PASSED"
  else
    result_str = "FAILED"
  endif
  
  write(LOG_UNIT, '(A,A,A,A)') "CONSTRAINT: ", trim(constraintName), " -> ", result_str
  
  if(present(details)) then
    write(LOG_UNIT, '(A,A)') "  Details: ", trim(details)
  endif
  
end subroutine

!========================================================================================
! Dump deletion proposal (marks atoms to be deleted)
!========================================================================================
subroutine DebugLog_DumpProposalDeletion(box, disp, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  use Common_MolInfo, only: nMolTypes
  
  class(SimpleBox), intent(inout) :: box
  type(Deletion), intent(in) :: disp(:)
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
  
  integer :: iAtom, iDisp, iType, typeStart, typeEnd
  integer :: atomType, nDim, nDisp
  integer :: delAtom, molStart, molEnd
  real(dp) :: localBoxDim(2,3)
  logical :: isDeleted
  
  if(.not. logging_initialized) return
  
  nDim = box%nDimension
  nDisp = size(disp)
  
  dump_frame_counter = dump_frame_counter + 1
  log_entry_counter = log_entry_counter + 1
  
  if(present(boxDim)) then
    localBoxDim = boxDim
  else
    localBoxDim(1,:) = -50.0_dp
    localBoxDim(2,:) = 50.0_dp
    call GetBoxDimensions(box, localBoxDim)
  endif
  
  write(DUMP_UNIT, '(A)') "ITEM: TIMESTEP"
  write(DUMP_UNIT, '(I12)') dump_frame_counter
  
  write(DUMP_UNIT, '(A)') "ITEM: NUMBER OF ATOMS"
  write(DUMP_UNIT, '(I12)') box%nAtoms
  
  write(DUMP_UNIT, '(A)') "ITEM: BOX BOUNDS pp pp pp"
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,1), localBoxDim(2,1)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,2), localBoxDim(2,2)
  write(DUMP_UNIT, '(2E20.10)') localBoxDim(1,3), localBoxDim(2,3)
  
  write(DUMP_UNIT, '(A)') "ITEM: ATOMS id type mol x y z c_deleted"
  
  do iType = 1, nMolTypes
    call box%GetTypeAtoms(iType, typeStart, typeEnd)
    if(typeStart < 1) cycle
    do iAtom = typeStart, typeEnd
      atomType = box%AtomType(iAtom)
      isDeleted = .false.
      do iDisp = 1, nDisp
        call box%GetMolData(disp(iDisp)%molIndx, molStart=molStart, molEnd=molEnd)
        if(iAtom >= molStart .and. iAtom <= molEnd) then
          isDeleted = .true.
          exit
        endif
      enddo
      if(isDeleted) then
        write(DUMP_UNIT, '(I8,2I4,3E20.10,I4)') iAtom, atomType + 100, box%MolIndx(iAtom), &
                                                box%atoms(1,iAtom), box%atoms(2,iAtom), box%atoms(3,iAtom), 1
      else
        write(DUMP_UNIT, '(I8,2I4,3E20.10,I4)') iAtom, atomType, box%MolIndx(iAtom), &
                                                box%atoms(1,iAtom), box%atoms(2,iAtom), box%atoms(3,iAtom), 0
      endif
    enddo
  enddo
  
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, '(A,A)') "DELETION PROPOSAL: ", trim(label)
  write(LOG_UNIT, '(A,I12)') "  Frame Number: ", dump_frame_counter
  write(LOG_UNIT, '(A,I8)') "  Box ID: ", box%boxID
  write(LOG_UNIT, '(A,I4)') "  Number of Molecules to Delete: ", nDisp
  write(LOG_UNIT, '(A)') "  ---"
  write(LOG_UNIT, '(A)') "  Deletion Details:"
  
  do iDisp = 1, nDisp
    write(LOG_UNIT, '(A,I4)') "    Molecule ", iDisp
    write(LOG_UNIT, '(A,I4)') "      Mol Type: ", disp(iDisp)%molType
    write(LOG_UNIT, '(A,I8)') "      Mol Index: ", disp(iDisp)%molIndx
  enddo
  write(LOG_UNIT, '(A)') "========================================================"
  
end subroutine

!========================================================================================
! Dump volume change proposal
!========================================================================================
subroutine DebugLog_DumpVolumeChange(box, volOld, volNew, label)
  use SimpleSimBox, only: SimpleBox
  
  class(SimpleBox), intent(inout) :: box
  real(dp), intent(in) :: volOld, volNew
  character(len=*), intent(in) :: label
  
  if(.not. logging_initialized) return
  
  log_entry_counter = log_entry_counter + 1
  
  write(LOG_UNIT, '(A)') "========================================================"
  write(LOG_UNIT, '(A,A)') "VOLUME CHANGE PROPOSAL: ", trim(label)
  write(LOG_UNIT, '(A,I8)') "  Box ID: ", box%boxID
  write(LOG_UNIT, '(A,E20.10)') "  Old Volume: ", volOld
  write(LOG_UNIT, '(A,E20.10)') "  New Volume: ", volNew
  write(LOG_UNIT, '(A,E20.10)') "  Volume Ratio: ", volNew/volOld
  write(LOG_UNIT, '(A,E20.10)') "  Scale Factor: ", (volNew/volOld)**(1.0_dp/3.0_dp)
  write(LOG_UNIT, '(A)') "========================================================"
  
end subroutine

#else

!========================================================================================
! Stub implementations when LOGGING is not defined
! These do nothing but allow the code to compile without #ifdef everywhere
!========================================================================================
contains

function DebugLog_IsEnabled() result(enabled)
  logical :: enabled
  enabled = .false.
end function

subroutine DebugLog_Initialize(logfile, dumpfile)
  character(len=*), intent(in), optional :: logfile, dumpfile
end subroutine

subroutine DebugLog_Finalize()
end subroutine

subroutine DebugLog_Separator(char)
  character(len=1), intent(in), optional :: char
end subroutine

subroutine DebugLog_Message(message, level)
  character(len=*), intent(in) :: message
  integer, intent(in), optional :: level
end subroutine

subroutine DebugLog_DumpBox(box, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  class(SimpleBox), intent(inout) :: box
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
end subroutine

subroutine DebugLog_DumpProposal(box, disp, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  class(SimpleBox), intent(inout) :: box
  type(Displacement), intent(in) :: disp(:)
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
end subroutine

subroutine DebugLog_DumpProposalAddition(box, disp, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  class(SimpleBox), intent(inout) :: box
  type(Addition), intent(in) :: disp(:)
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
end subroutine

subroutine DebugLog_DumpProposalDeletion(box, disp, label, boxDim)
  use SimpleSimBox, only: SimpleBox
  class(SimpleBox), intent(inout) :: box
  type(Deletion), intent(in) :: disp(:)
  character(len=*), intent(in) :: label
  real(dp), intent(in), optional :: boxDim(2,3)
end subroutine

subroutine DebugLog_DumpVolumeChange(box, volOld, volNew, label)
  use SimpleSimBox, only: SimpleBox
  class(SimpleBox), intent(inout) :: box
  real(dp), intent(in) :: volOld, volNew
  character(len=*), intent(in) :: label
end subroutine

subroutine DebugLog_MoveResult(moveName, accepted, E_Diff, E_Inter, E_Intra, biasProb)
  character(len=*), intent(in) :: moveName
  logical, intent(in) :: accepted
  real(dp), intent(in), optional :: E_Diff, E_Inter, E_Intra, biasProb
end subroutine

subroutine DebugLog_ConstraintCheck(constraintName, passed, details)
  character(len=*), intent(in) :: constraintName
  logical, intent(in) :: passed
  character(len=*), intent(in), optional :: details
end subroutine

#endif

end module Debug_Logging
!========================================================================================
