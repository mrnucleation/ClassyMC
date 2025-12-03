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

    write(0,*) "Error! FindBond function unable to find a bond"
    write(0,*) "containing memebers:", mem1, mem2
    error stop
  
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

    write(0,*) "Error! FindAngle function unable to find a angle"
    write(0,*) "containing memebers:", mem1, mem2, mem3
    error stop
  
  end subroutine
!==========================================================================           
  subroutine FindTorsion(molType, mem1, mem2, mem3,mem4, torsType)
    use Common_MolInfo, only: MolData, AngleData
    implicit none
    integer, intent(in) :: molType, mem1, mem2, mem3, mem4
    integer, intent(out) :: torsType
    integer :: iTors
  
    do iTors = 1, size(MolData(molType)%torsion)
      if( (MolData(molType)%torsion(iTors)%mem1 .eq. mem1) .or. &
          (MolData(molType)%torsion(iTors)%mem2 .eq. mem1) .or. &
          (MolData(molType)%torsion(iTors)%mem3 .eq. mem1)  .or. &
          (MolData(molType)%torsion(iTors)%mem4 .eq. mem1) ) then
          if( (MolData(molType)%torsion(iTors)%mem1 .eq. mem2) .or. &
              (MolData(molType)%torsion(iTors)%mem2 .eq. mem2) .or. &
              (MolData(molType)%torsion(iTors)%mem3 .eq. mem2) .or. &
              (MolData(molType)%torsion(iTors)%mem4 .eq. mem2) ) then
              if( (MolData(molType)%torsion(iTors)%mem1 .eq. mem3) .or. &
                  (MolData(molType)%torsion(iTors)%mem2 .eq. mem3) .or. &
                  (MolData(molType)%torsion(iTors)%mem3 .eq. mem3) .or. &
                  (MolData(molType)%torsion(iTors)%mem4 .eq. mem3) ) then
                  if( (MolData(molType)%torsion(iTors)%mem1 .eq. mem4) .or. &
                      (MolData(molType)%torsion(iTors)%mem2 .eq. mem4) .or. &
                      (MolData(molType)%torsion(iTors)%mem3 .eq. mem4) .or. &
                      (MolData(molType)%torsion(iTors)%mem4 .eq. mem4) ) then
                       TorsType = MolData(molType)%torsion(iTors)%TorsType
                       return
                 endif
              endif
          endif
      endif
    enddo

    write(0,*) "Error! FindTorsion function unable to find a angle"
    write(0,*) "containing memebers:", mem1, mem2, mem3, mem4
    error stop
  
  end subroutine

!==========================================================================           
end module
