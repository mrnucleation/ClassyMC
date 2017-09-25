!==============================================================
      module Constants
      use VarPrecision
      real(dp),parameter :: pi=4d0*datan(1d0) 
      real(dp),parameter :: two_pi=8d0*datan(1d0)
        
      end module
!===================================================================      
      module Units
      use VarPrecision
      contains
        
!     !----------------------------------------------------------
        real(dp) function FindEngUnit(unitName)
        implicit none 
        character(len=*), intent(in) :: unitName       
        
        
        select case(trim(adjustl(unitName)))
          case("j-mol")
                FindEngUnit = 1d0/8.3144621d0
          case("kj-mol")
                FindEngUnit = 1d0/8.3144621d-3
          case("cal-mol")
                FindEngUnit = 1d0/1.9872041d0
          case("kcal-mol")
                FindEngUnit = 1d0/1.9872041d-3
          case("eV") 
                FindEngUnit = 1d0/8.6173303d-5
          case("ev") 
                FindEngUnit = 1d0/8.6173303d-5
          case("kB")
                FindEngUnit = 1d0
          case("kb")
                FindEngUnit = 1d0
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
        use Constants         
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
