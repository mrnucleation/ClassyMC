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
module Constants
  use VarPrecision
  real(dp), parameter :: pi=4d0*datan(1d0) 
  real(dp), parameter :: two_pi=8d0*datan(1d0)
        
end module
!===================================================================      
module Units
  use VarPrecision
  real(dp) :: outEngUnit = 1E0_dp
  real(dp) :: outLenUnit = 1E0_dp
  real(dp) :: outAngUnit = 1E0_dp
  contains
        
!     !----------------------------------------------------------
   function FindEngUnit(unitName) result(units)
     implicit none 
     character(len=*), intent(in) :: unitName          
     real(dp) :: units
        
     select case(trim(adjustl(unitName)))
       case("j-mol")
         units = 1d0/8.3144621d0
       case("kj-mol")
         units = 1d0/8.3144621d-3
       case("cal-mol")
         units = 1d0/1.9872041d0
       case("kcal-mol")
         units = 1d0/1.9872041d-3
       case("ev") 
         units = 1d0/8.6173303d-5
       case("kb")
         units = 1d0
       case default
         write(*,*) "Error! Invalid Energy Unit Type!"
         write(*,*) unitName
         stop
       end select
        
     end function
!    !----------------------------------------------------------
     real(dp) function FindLengthUnit(unitName)
       implicit none 
       character(len=*) :: unitName       
        
       select case(trim(adjustl(unitName)))
         case("nm")
           FindLengthUnit = 1d-1
         case("a")
           FindLengthUnit = 1d0
         case("ang")
           FindLengthUnit = 1d0
         case("sigma")
           FindLengthUnit = 1d0                  
         case default
           write(*,*) "Error! Invalid Length Unit Type!"
           write(*,*) unitName
           stop
        end select
 
        
      end function
!     !----------------------------------------------------------
  real(dp) function FindAngularUnit(unitName)
    use Constants, only: pi
    implicit none 
    character(len=*), intent(in) :: unitName       
        
    select case(trim(adjustl(unitName)))
      case("deg")
        FindAngularUnit = pi/180d0
      case("degree")
        FindAngularUnit = pi/180d0
      case("degrees")
        FindAngularUnit = pi/180d0                  
      case("rad")
        FindAngularUnit = 1d0
      case("radian")
        FindAngularUnit = 1d0
      case("radians")
        FindAngularUnit = 1d0                 
      case default
        write(*,*) "Error! Invalid Angular Unit Type!"
        write(*,*) unitName
        stop
      end select
   
  end function          
!      !----------------------------------------------------------
end module    
!===================================================================
