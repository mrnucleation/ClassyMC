!==========================================================================           
module MolSearch
!==========================================================================           
contains
!==========================================================================           
  subroutine FindBond(molType, mem1, mem2, bondType)
    use Common_MolInfo, only: MolData, BondData
    implicit none
    integer, intent(in) :: molType, mem1, mem2
    integer, intent(out) :: bondType
    integer :: iBond
  
    do iBond = 1, size(MolData(molType)%bond)
      if( (MolData(molType)%bond(iBond)%mem1 .eq. mem1) .or. &
          (MolData(molType)%bond(iBond)%mem2 .eq. mem1) ) then
          if( (MolData(molType)%bond(iBond)%mem1 .eq. mem2) .or. &
              (MolData(molType)%bond(iBond)%mem2 .eq. mem2) ) then
                 bondType = MolData(molType)%bond(iBond)%bondType
                 return
          endif
      endif
    enddo

    write(*,*) "Error! FindBond function unable to find a bond"
    write(*,*) "containing memebers:", mem1, mem2
    stop
  
  end subroutine
!==========================================================================           
  subroutine FindAngle(molType, mem1, mem2, mem3, angleType)
    use Common_MolInfo, only: MolData, AngleData
    implicit none
    integer, intent(in) :: molType, mem1, mem2, mem3
    integer, intent(out) :: angleType
    integer :: iAngle
  
    do iAngle = 1, size(MolData(molType)%angle)
      if( (MolData(molType)%angle(iAngle)%mem1 .eq. mem1) .or. &
          (MolData(molType)%angle(iAngle)%mem2 .eq. mem1) .or. &
          (MolData(molType)%angle(iAngle)%mem3 .eq. mem1) ) then
          if( (MolData(molType)%angle(iAngle)%mem1 .eq. mem2) .or. &
              (MolData(molType)%angle(iAngle)%mem2 .eq. mem2) .or. &
              (MolData(molType)%angle(iAngle)%mem3 .eq. mem2) ) then
              if( (MolData(molType)%angle(iAngle)%mem1 .eq. mem3) .or. &
                  (MolData(molType)%angle(iAngle)%mem2 .eq. mem3) .or. &
                  (MolData(molType)%angle(iAngle)%mem3 .eq. mem3) ) then
                    angleType = MolData(molType)%angle(iAngle)%angleType
                    return
              endif
          endif
      endif
    enddo

    write(*,*) "Error! FindAngle function unable to find a angle"
    write(*,*) "containing memebers:", mem1, mem2, mem3
    stop
  
  end subroutine
!==========================================================================           
end module
