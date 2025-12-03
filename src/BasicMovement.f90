!===========================================================================================
      module SimpleMCMoves_Module
      contains
!===========================================================================================
      subroutine SingleAtom_Translation(E_T, acc_x, atmp_x)
      use CBMC_Variables      
      use Coords
!      use E_Interface
      use DistanceCriteria      
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies 
      use EnergyCriteria
      use EnergyTables
      use Forcefield
      use IndexingFunctions
      use PairStorage, only: UpdateDistArray, useDistStore
      use Pressure_LJ_Electro, only: Shift_PressCalc_Inter
      use SimParameters
      implicit none
      
      real(dp),intent(inout) :: E_T,acc_x,atmp_x      
      real(dp) :: max_distx

      logical, parameter :: useIntra(1:4) = [.true., .true., .true., .true.]
      
      logical :: rejMove      
      integer :: nType,nMol,nIndx,nMove, nAtom
      real(dp) :: grnd 
      real(dp) :: dx,dy,dz      
      real(dp) :: E_Diff,E_Inter, E_Intra
      type (displacement) :: disp(1:1)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
      prevMoveAccepted = .false.
      max_distx = 0.1E0
      rejMove = .false.
      atmp_x = atmp_x + 1E0
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)

      call Get_MolIndex(nMove, NPart, nType, nMol)
!      if(nType .ne. 1) then
!        return
!      endif
      nIndx = MolArray(nType)%mol(nMol)%indx
      if(regrowType(nType) .eq. 0) return
      
      nAtom = floor(nAtoms(nType)*grnd() + 1E0)
      

!     Generate a random translational displacement             
      dx = max_dist_single(nType) * (2E0*grnd() - 1E0)
      dy = max_dist_single(nType) * (2E0*grnd() - 1E0)
      dz = max_dist_single(nType) * (2E0*grnd() - 1E0)

!     Construct the Displacement Vectors for each atom in the molecule that was chosen.'
      disp(1)%molType = int(nType,atomIntType)
      disp(1)%molIndx = int(nMol,atomIntType)
      disp(1)%atmIndx = int(nAtom,atomIntType)
        
      disp(1)%x_old => MolArray(nType)%mol(nMol)%x(nAtom)
      disp(1)%y_old => MolArray(nType)%mol(nMol)%y(nAtom)
      disp(1)%z_old => MolArray(nType)%mol(nMol)%z(nAtom)

      disp(1)%x_new = disp(1)%x_old + dx
      disp(1)%y_new = disp(1)%y_old + dy
      disp(1)%z_new = disp(1)%z_old + dz

      
!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
      E_Inter = 0E0
      E_Intra = 0E0
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:1), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:1), PairList, dETable, useIntra, rejMove)
      if(rejMove) return
      E_Diff = E_Inter + E_Intra

!     Calculate Acceptance and determine if the move is accepted or not     
      if(E_Diff .le. 0E0) then
        if(calcPressure) then
          call Shift_PressCalc_Inter(P_Diff, disp(1:1))
          pressure = pressure + P_Diff
        endif
        disp(1)%x_old = disp(1)%x_new
        disp(1)%y_old = disp(1)%y_new
        disp(1)%z_old = disp(1)%z_new
        E_T = E_T + E_Diff
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        if(useDistStore) then
          call UpdateDistArray
        endif
        if( distCriteria ) then
!          call NeighborUpdate_Distance(PairList,nIndx)    
          call NeighborUpdate_Distance(PairList, nIndx)            
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
        call Update_SubEnergies        
        prevMoveAccepted = .true.
      elseif(-beta*E_Diff .gt. log(grnd())) then
        if(calcPressure) then
          call Shift_PressCalc_Inter(P_Diff, disp(1:1))
          pressure = pressure + P_Diff
        endif
        disp(1)%x_old = disp(1)%x_new
        disp(1)%y_old = disp(1)%y_new
        disp(1)%z_old = disp(1)%z_new
        E_T = E_T + E_Diff
        ETable = ETable + dETable
        acc_x = acc_x + 1E0
        if(useDistStore) then
          call UpdateDistArray
        endif
        if( distCriteria ) then
!          call NeighborUpdate_Distance(PairList,nIndx)    
          call NeighborUpdate_Distance(PairList, nIndx)            
        else
          call NeighborUpdate(PairList, nIndx)
        endif  
        call Update_SubEnergies        
        prevMoveAccepted = .true.
      endif

     
      end subroutine
!===========================================================================================
      subroutine Translation(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables
      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none
      
      real(dp),intent(inout) :: E_T, acc_x, atmp_x      
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: iAtom, nType, nMol, nIndx, nMove
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: grnd
      real(dp) :: dx, dy, dz      
      real(dp) :: E_Inter, E_Intra
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

      prevMoveAccepted = .false.
      
      if(NTotal .eq. 1) return

     
      rejMove = .false.
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)

      atmp_x = atmp_x + 1E0
      atmpTrans(nType) = atmpTrans(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx
!     Generate a random translational displacement             
      dx = max_dist(nType) * (2E0*grnd() - 1E0)
      dy = max_dist(nType) * (2E0*grnd() - 1E0)
      dz = max_dist(nType) * (2E0*grnd() - 1E0)

!     Construct the Displacement Vectors for each atom in the molecule that was chosen.
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%molType = int(nType, atomIntType)
        disp(iAtom)%molIndx = int(nMol, atomIntType)
        disp(iAtom)%atmIndx = int(iAtom, atomIntType)
        
        disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
        disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
        disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
        
        disp(iAtom)%x_new = disp(iAtom)%x_old + dx
        disp(iAtom)%y_new = disp(iAtom)%y_old + dy
        disp(iAtom)%z_new = disp(iAtom)%z_old + dz        
      enddo
      
!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
      E_Inter = 0E0
      E_Intra = 0E0
      !call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) then
        return
      endif
      
      biasDiff = 0E0
!      write(*,*) useUmbrella
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
        

!     Calculate Acceptance and determine if the move is accepted or not     
      if(biasEnergy .le. 0E0) then
        acptTrans(nType) = acptTrans(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif(-biasEnergy .gt. log(grnd())) then
        acptTrans(nType) = acptTrans(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      endif

      	  
      end subroutine
!===========================================================================================
      subroutine Rotation(E_T, acc_x, atmp_x)
      use SimParameters    
      implicit none
      real(dp), intent(inout) :: atmp_x,acc_x,E_T
      real(dp) :: ran_num, grnd


      if(NTotal .eq. 1) then
!       acc_x=acc_x+1E0
       return            
      endif      

      prevMoveAccepted = .false.      

      ran_num = grnd()
      if(ran_num .lt. 1E0/3E0) then
         call Rot_xy(E_T, acc_x, atmp_x)      
      elseif(ran_num .lt. 2E0/3E0) then
         call Rot_xz(E_T, acc_x, atmp_x)  
      else
         call Rot_yz(E_T, acc_x, atmp_x)  
      endif      
      
      end subroutine   

!=======================================================      
      subroutine Rot_xy(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies 
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables
!      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T,acc_x,atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: iAtom, nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd   
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term,s_term
      real(dp) :: x_scale, y_scale
      real(dp) :: xcm,ycm,angle
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if(nAtoms(nType) .eq. 1) then
        return
      endif
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array 
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%molType = int(nType, atomIntType)
        disp(iAtom)%molIndx = int(nMol, atomIntType)
        disp(iAtom)%atmIndx = int(iAtom, atomIntType)
        
        disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
        disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
        disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
      enddo
      
!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd()-1E0)
      c_term = cos(angle)
      s_term = sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion. 
      xcm=0E0
      ycm=0E0
      do iAtom = 1, nAtoms(nType)
        atmType = atomArray(nType, iAtom)
        xcm = xcm + atomData(atmType)%mass*disp(iAtom)%x_old
        ycm = ycm + atomData(atmType)%mass*disp(iAtom)%y_old
      enddo

      xcm = xcm/totalMass(nType)   
      ycm = ycm/totalMass(nType)

!     Generate a random translational displacement      
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%z_new = disp(iAtom)%z_old
        x_scale = disp(iAtom)%x_old - xcm
        y_scale = disp(iAtom)%y_old - ycm
        disp(iAtom)%x_new = c_term*x_scale - s_term*y_scale + xcm
        disp(iAtom)%y_new = s_term*x_scale + c_term*y_scale + ycm
      enddo

!     Calculate the Energy Difference Associated with the move   
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) return

      biasDiff = 0E0
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
      if(biasEnergy .le. 0E0) then
        acptRot(nType) = acptRot(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)      
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif(-biasEnergy .gt. log(grnd())) then
        acptRot(nType) = acptRot(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      endif
      end subroutine
!=======================================================      
      subroutine Rot_xz(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria      
      use EnergyTables      
!      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T,acc_x,atmp_x

      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: iAtom, nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: angle
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd   
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term,s_term
      real(dp) :: x_scale, z_scale
      real(dp) :: xcm,zcm
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

      prevMoveAccepted = .false.
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if(nAtoms(nType) .eq. 1) then
        return
      endif
      atmp_x = atmp_x+1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array 
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%molType = int(nType, atomIntType)
        disp(iAtom)%molIndx = int(nMol, atomIntType)
        disp(iAtom)%atmIndx = int(iAtom, atomIntType)
        
        disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
        disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
        disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
      enddo
      
!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType) * (2E0 * grnd() - 1E0)
      c_term = cos(angle)
      s_term = sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion. 
      xcm = 0E0
      zcm = 0E0
      do iAtom = 1, nAtoms(nType)
        atmType = atomArray(nType, iAtom)
        xcm = xcm + atomData(atmType)%mass*disp(iAtom)%x_old
        zcm = zcm + atomData(atmType)%mass*disp(iAtom)%z_old
      enddo

      xcm = xcm/totalMass(nType)    
      zcm = zcm/totalMass(nType)

!     Generate a random translational displacement      
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%y_new = disp(iAtom)%y_old
        x_scale = disp(iAtom)%x_old - xcm
        z_scale = disp(iAtom)%z_old - zcm
        disp(iAtom)%x_new = c_term*x_scale - s_term*z_scale + xcm
        disp(iAtom)%z_new = s_term*x_scale + c_term*z_scale + zcm
      enddo

!     Calculate the Energy Difference Associated with the move   
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) return

      biasDiff = 0E0
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
      if(biasEnergy .le. 0E0) then
        acptRot(nType) = acptRot(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)      
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif(-biasEnergy .gt. log(grnd())) then
        acptRot(nType) = acptRot(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      endif

      end subroutine
!=======================================================      
      subroutine Rot_yz(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
!      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none
      
      real(dp), intent(inout) :: E_T,acc_x,atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]      
      
      logical :: rejMove      
      integer :: iAtom, nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: angle
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd   
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term,s_term
      real(dp) :: y_scale, z_scale
      real(dp) :: ycm,zcm
      type (displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if(nAtoms(nType) .eq. 1) then
        return
      endif
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array 
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%molType = int(nType, atomIntType)
        disp(iAtom)%molIndx = int(nMol, atomIntType)
        disp(iAtom)%atmIndx = int(iAtom, atomIntType)
        
        disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
        disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
        disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
      enddo
      
!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd()-1E0)
      c_term=cos(angle)
      s_term=sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion. 
      ycm=0E0
      zcm=0E0
      do iAtom = 1, nAtoms(nType)
        atmType = atomArray(nType,iAtom)
        ycm = ycm + atomData(atmType)%mass * disp(iAtom)%y_old
        zcm = zcm + atomData(atmType)%mass * disp(iAtom)%z_old
      enddo

      ycm = ycm/totalMass(nType)
      zcm = zcm/totalMass(nType)

!     Generate a random translational displacement      
      do iAtom = 1, nAtoms(nType)
        disp(iAtom)%x_new = disp(iAtom)%x_old
        y_scale = disp(iAtom)%y_old - ycm
        z_scale = disp(iAtom)%z_old - zcm
        disp(iAtom)%y_new = c_term*y_scale - s_term*z_scale + ycm
        disp(iAtom)%z_new = s_term*y_scale + c_term*z_scale + zcm
      enddo

!     Calculate the Energy Difference Associated with the move   
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if(rejMove) return

      biasDiff = 0E0
      if(useUmbrella) then
        call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = beta*E_Inter - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
      if(biasEnergy .le. 0E0) then
        acptRot(nType) = acptRot(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)      
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif(-biasEnergy .gt. log(grnd())) then
        acptRot(nType) = acptRot(nType) + 1E0
        call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      endif
      end subroutine


!=======================================================      
      subroutine Update_Shift(disp, nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      use AcceptRates
      use SimParameters
      use Coords
      use Forcefield, only: nAtoms
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use PairStorage, only: UpdateDistArray, useDistStore
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      use Pressure_LJ_Electro, only: Shift_PressCalc_Inter
      implicit none
      
      real(dp), intent(inout) :: E_T,acc_x,atmp_x
      integer, intent(in) :: nIndx, nType
      type (displacement), intent(inout) :: disp(:)
      real(dp), intent(in) :: PairList(1:maxMol)
      real(dp), intent(in) :: dETable(1:maxMol)
      real(dp), intent(in) :: E_Inter

      integer :: iAtom     
    

      if(calcPressure) then
        call Shift_PressCalc_Inter(P_Diff, disp(1:nAtoms(nType) ) )
        pressure = pressure + P_Diff
      endif

      do iAtom = 1, nAtoms(nType)      
        disp(iAtom)%x_old = disp(iAtom)%x_new
        disp(iAtom)%y_old = disp(iAtom)%y_new
        disp(iAtom)%z_old = disp(iAtom)%z_new
      enddo
      E_T = E_T + E_Inter
      ETable = ETable + dETable
      acc_x = acc_x + 1E0
      if(useDistStore) then
        call UpdateDistArray
      endif
      if( distCriteria ) then
!        call NeighborUpdate_Distance(PairList,nIndx)    
        call NeighborUpdate_Distance(PairList, nIndx)            
      else
        call NeighborUpdate(PairList, nIndx)
      endif  
      call Update_SubEnergies        
      prevMoveAccepted = .true.

      end subroutine

!=======================================================      
!     Experimental Temperature Change Move.  Not guarenteed to give accurate results.
      subroutine TemperatureMove(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Temperature
      implicit none

      real(dp), intent(inout) :: E_T,acc_x,atmp_x   
      real(dp), parameter :: power = (6d0/2d0)      


      logical :: rejMove      
      integer :: iAtom, nMove 
      integer :: atmType,nMol,nType,nIndx
      real(dp) :: grnd, betaNew
      real(dp) :: biasNew, biasOld, biasDiff, biasEnergy

      if(NTotal .ne. 1) then
        return
      endif
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      atmp_x = atmp_x + 1E0_dp
      TempNew = temperature + 15E0_dp * (2E0_dp*grnd()-1E0_dp)
      betaNew = 1E0_dp/TempNew
      biasOld = 0E0_dp 
      biasNew = 0E0_dp
      if(useUmbrella) then
        call GetUmbrellaBias_Temperature(biasDiff, rejMove)
        if(rejMove) then
          return
        endif
      endif
      biasEnergy = (betaNew-beta)*E_T - biasDiff
       
!      Calculate Acceptance and determine if the move is accepted or not       
!      if(biasEnergy .le. 0E0_dp) then
!        acc_x = acc_x + 1E0_dp     
!        temperature = TempNew 
!        beta = betaNew
      if((TempNew/temperature)**power * exp(-biasEnergy) .gt. grnd()) then
        acc_x = acc_x + 1E0_dp
        temperature = TempNew
        beta = betaNew
      endif
!      write(*,*) temperature

      end subroutine
!===========================================================================================
      end module
