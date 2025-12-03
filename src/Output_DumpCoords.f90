!============================================================
module Output_DumpCoords
contains
!============================================================
  subroutine Output_DumpData
  use ParallelVar, only: myid
  use BoxData, only: BoxArray
  implicit none
  integer :: iBox
  character(len=50) :: format_string
  character(len=500) :: fileName

  if (myid .lt. 10) then
    format_string = "(A,I1,A,I1,A)"
  elseif(myid .lt. 100) then
    format_string = "(A,I2,A,I1,A)"
  elseif(myid .lt. 1000) then
    format_string = "(A,I3,A,I1,A)"
  elseif(myid .lt. 10000) then
    format_string = "(A,I4,A,I1,A)"
  else
    format_string = "(A,I5,A,I1,A)"
  endif

  do iBox = 1, size(BoxArray)
    write(fileName, format_string) "Config_", myid,"_", iBox, ".clssy"
    call BoxArray(iBox) % box % DumpData(filename)
  enddo
  

  end subroutine 
!============================================================
end module
!============================================================
