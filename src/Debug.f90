!===========================================
module Debug
!===========================================
contains
!===========================================
  subroutine Debug_DumpNeiList(iAtom, boxNum, listNum)
  use BoxData, only: BoxArray
  implicit none
  integer, intent(in) :: iAtom, listNum, boxNum
  integer :: iNei, jNei, nNei

!  if(.not. inquire(file="listdump.dat", slist="od")) then
!    open(unit=35, file="listdump.dat")
!  endif

!  do iAtom = 1, size(BoxArray(boxNum)%box%NeighList(listNum)%nNeigh)
    nNei = BoxArray(boxNum)%box%NeighList(listNum)%nNeigh(iAtom)
    write(35,"(I5, A, I5, A, 10000I5)") iAtom, "/", nNei,"/", &
        (BoxArray(boxNum)%box%NeighList(listNum)%list(jNei,iAtom), jNei=1,nNei)
!  enddo

!  close(35)

  end subroutine
!===========================================
end module
!===========================================