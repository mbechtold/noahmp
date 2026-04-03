module PeatlandPhysicsMod

!!! Specifies peatland physical processes specific options
!!! Introduced by Chakraborty et al., (2025)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ApplyPeatlandPhysics(noahmp)
  
! ------------------------ Code history --------------------------------------------------
! Turning on peatland physical processes (Chakraborty & Bechtold, 2025)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                                             &
              OptSoilWaterTranspiration => noahmp%config%nmlist%OptSoilWaterTranspiration ,&
              OptRunoffSubsurface       => noahmp%config%nmlist%OptRunoffSubsurface       ,&
              OptPeatlandPhysics       => noahmp%config%nmlist%OptPeatlandPhysics         &
             )
! ----------------------------------------------------------------------


    ! Set peatland-specific physics option for runoff and transpiration
    !write(*,*) "DEBUG: OptPeatlandPhysics =", OptPeatlandPhysics

    if ( OptPeatlandPhysics == 1 ) then
        OptRunoffSubsurface = 9
        OptSoilWaterTranspiration = 4
        !write(*,*) "Peatland Physics option working with OptRunoffSubsurface option 9 and OptSoilWaterTranspiration option 4"
    endif
    
    end associate

  end subroutine ApplyPeatlandPhysics

end module PeatlandPhysicsMod
