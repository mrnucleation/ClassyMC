!==============================================================
! This file contains physical constants and unit conversion factors
! that can be used in the simulation code. The factors contained
! within all units are scaled into units native to the simulation code.
!
! Energy = Boltzmann units (Defined so the boltzman constant = 1 K^-1)
! Length = Angstroms
! Temperature = Kelvin
! Angles = Radians
! 
! All other units are expressed in these terms.
! 
! Using this convension it is also possible to work in reduced
! units such as the Lennard-Jones units.
!==============================================================
module ClassyConstants
  use VarPrecision
  real(dp), parameter :: pi=4E0_dp*datan(1E0_dp) 
  real(dp), parameter :: two_pi=8E0_dp*datan(1E0_dp)
        
  !All units are in Boltzmann Reduced Units (kb = 1/K)
  real(dp), parameter :: coulombConst = 1.671009770E5_dp
end module
!===================================================================      
module Units
  use VarPrecision
  real(dp), parameter :: PaToBoltz = 0.724297303E-7_dp
  real(dp) :: inEngUnit = 1E0_dp
  real(dp) :: inLenUnit = 1E0_dp
  real(dp) :: inAngUnit = 1E0_dp
  real(dp) :: inPressUnit = 1E0_dp
  character(len=6) :: inengStr = "kb"
  character(len=6) :: inlenStr = "ang"
  character(len=6) :: inangStr = "rad"
  character(len=6) :: inPressStr = "atm"

  real(dp) :: outEngUnit = 1E0_dp
  real(dp) :: outLenUnit = 1E0_dp
  real(dp) :: outAngUnit = 1E0_dp
  real(dp) :: outPressUnit = 1E0_dp
  character(len=6) :: engStr = "kb"
  character(len=6) :: lenStr = "ang"
  character(len=6) :: angStr = "rad"
  character(len=6) :: pressStr = "atm"


!===================================================================      
  contains
!===================================================================      
  function FindEngUnit(unitName) result(units)
    use ParallelVar, only: nout
     implicit none 
     character(len=*), intent(in) :: unitName          
     real(dp) :: units
        
     write(nout, *) "Setting Energy Units to: ", trim(adjustl(unitName))
!     engStr = 
     select case(trim(adjustl(unitName)))
       case("j-mol")
         units = 1E0_dp/8.3144621E0_dp
       case("kj-mol")
         units = 1E0_dp/8.3144621d-3
       case("cal-mol")
         units = 1E0_dp/1.9872041E0_dp
       case("kcal-mol")
         units = 1E0_dp/1.9872041d-3
       case("ev") 
         units = 1E0_dp/8.6173303d-5
       case("kb")
         units = 1E0_dp
       case default
         write(*,*) "Error! Invalid Energy Unit Type!"
         write(*,*) unitName
         stop
       end select
        
     end function
!===================================================================      
  function FindLengthUnit(unitName) result(units)
     implicit none 
     character(len=*) :: unitName       
     real(dp) :: units
      
     select case(trim(adjustl(unitName)))
       case("nm")
         units = 1E-1_dp
       case("a")
         units = 1E0_dp
       case("ang")
         units = 1E0_dp
       case("sigma")
         units = 1E0_dp                  
       case default
         write(*,*) "Error! Invalid Length Unit Type!"
         write(*,*) unitName
         stop
      end select

      
  end function
!===================================================================      
  function FindAngularUnit(unitName) result(units)
    use ClassyConstants, only: pi
    implicit none 
    character(len=*), intent(in) :: unitName       
     real(dp) :: units
        
    select case(trim(adjustl(unitName)))
      case("deg")
        units = pi/180E0_dp
      case("degree")
        units = pi/180E0_dp
      case("degrees")
        units = pi/180E0_dp                  
      case("rad")
        units = 1E0_dp
      case("radian")
        units = 1E0_dp
      case("radians")
        units = 1E0_dp                 
      case default
        write(*,*) "Error! Invalid Angular Unit Type!"
        write(*,*) unitName
        stop
      end select
   
  end function          
!===================================================================      
  function FindPressureUnit(unitName) result(units)
     implicit none 
     character(len=*) :: unitName       
     real(dp) :: units
      
     select case(trim(adjustl(unitName)))
       case("atm")
         units = 1E0_dp/(PaToBoltz*101325)

       case("pa")
         units = 1E0_dp/PaToBoltz

       case("bar")
         units = 1E0_dp/(PaToBoltz*1E-5)

       case("mmhg")
         units = 1E0_dp/(PaToBoltz*133.322E0_dp)

       case("lj")
         units = 1E0_dp

       case default
         write(*,*) "Error! Invalid Pressure Unit Type!"
         write(*,*) unitName
         stop
      end select

      
  end function

!===================================================================      
end module    
!===================================================================
