module RunoffSubSurfacePeatlandMod

!!! Calculate subsurface runoff using Ivanov-based Peatland Runoff Scheme - Enabled by (Chakraborty & Bechtold, 2025)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  implicit none

contains

  subroutine RunoffSubSurfacePeatland(noahmp)

! ------------------------ Code history --------------------------------------------------
! Modified to include Peatland-specific Ivanov-based runoff scheme (Chakraborty & Bechtold, 2025)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! Define double precision kind parameter
    integer, parameter :: dp = kind(1.0d0)

    ! Declare local peatland-specific parameters
    real(dp) :: Ksz_zero, m_Ivanov, v_slope
    real(dp) :: Ta, BFLOW

! --------------------------------------------------------------------
    associate(                                                           &
              SoilImpervFracMax => noahmp%water%state%SoilImpervFracMax ,& ! in,    maximum soil imperviousness fraction
              WaterTableDepth   => noahmp%water%state%WaterTableDepth   ,& ! out,   water table depth [m]
              RunoffSubsurface  => noahmp%water%flux%RunoffSubsurface    & ! out,   subsurface runoff [mm/s] 
             )
! ----------------------------------------------------------------------

    ! ------------------------------------------
    ! Option 9: Ivanov-based Peatland Runoff Scheme (Chakraborty & Bechtold, 2025)
    ! ------------------------------------------

    ! Assign parameter values for peatland runoff scheme
    Ksz_zero = 3165.38_dp      ! Saturated hydraulic conductivity [m^2/s]
    m_Ivanov = 2.06_dp         ! Ivanov exponent
    v_slope = 1.5e-08_dp       ! Slope factor for runoff generation [unitless]

    ! Compute transmissivity function (Ta) [m^2/s]
    Ta = (Ksz_zero * (24.5_dp + 100.0_dp * max(-0.2449_dp, WaterTableDepth))**(1.0_dp - m_Ivanov)) / &
         (100.0_dp * (m_Ivanov - 1.0_dp))

    ! Compute baseflow (BFLOW) in mm/s
    BFLOW = v_slope * Ta * 1000.0_dp  ! Convert from m/s to mm/s

    ! Compute subsurface runoff using Peatland-specific equation
    RunoffSubsurface = (1.0_dp - SoilImpervFracMax) * BFLOW
    
    RunoffSubsurface = min(0.0002_dp,RunoffSubSurface)
    RunoffSubsurface = max(0.0_dp, RunoffSubsurface)

    end associate

  end subroutine RunoffSubSurfacePeatland

end module RunoffSubSurfacePeatlandMod
