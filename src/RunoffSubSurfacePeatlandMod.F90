module RunoffSubSurfacePeatlandMod

!!! Calculate subsurface runoff using Ivanov-based Peatland Runoff Scheme
!!! Introduced by Chakraborty & Bechtold (2025), WTD equilibrium revised by Bechtold (2026)
!!!
!!! Note: WaterTableDepth is NOT re-diagnosed here. It is assumed to be
!!! already set by the previous timestep's end (via FindWaterTableFlat in
!!! SoilWaterMainMod) or by initialization. This avoids inconsistency when
!!! SoilLiqWater represents microtopo-integrated values rather than a flat
!!! 1D profile.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine RunoffSubSurfacePeatland(noahmp)

! ------------------------ Code history --------------------------------------------------
! Peatland-specific Ivanov-based runoff scheme (Chakraborty & Bechtold, 2025; revised Bechtold, 2026)
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
              WaterTableDepth   => noahmp%water%state%WaterTableDepth   ,& ! in,    water table depth [m] (already diagnosed)
              FSW_change         => noahmp%water%state%FSW_change        ,& ! inout, 
              RunoffSubsurface  => noahmp%water%flux%RunoffSubsurface    & ! out,   subsurface runoff [mm/s] 
             )
! ----------------------------------------------------------------------

    ! WaterTableDepth is assumed valid from the previous timestep end or
    ! from model initialization.  No re-diagnosis here — SoilLiqWater may
    ! hold microtopo-integrated values that are incompatible with the flat
    ! equilibrium deficit functions.

    ! ------------------------------------------
    ! Option 9: Ivanov-based Peatland Runoff Scheme (Chakraborty & Bechtold, 2025)
    ! ------------------------------------------

    ! Assign parameter values for peatland runoff scheme
    Ksz_zero = 3165.38_dp      ! Saturated hydraulic conductivity [m^2/s]
    m_Ivanov = 2.06_dp         ! Ivanov exponent
    v_slope = 1.5e-08_dp       ! Slope factor for runoff generation [unitless]

    ! Compute transmissivity function (Ta) [m^2/s]
    ! Clamp to -1.0 m (maximum allowed water above surface with microtopography)
    Ta = (Ksz_zero * (24.5_dp + 100.0_dp * max(-0.2449_dp, WaterTableDepth))**(1.0_dp - m_Ivanov)) / &
         (100.0_dp * (m_Ivanov - 1.0_dp))

    ! Compute baseflow (BFLOW) in mm/s
    BFLOW = v_slope * Ta * 1000.0_dp  ! Convert from m/s to mm/s

    ! Compute subsurface runoff using Peatland-specific equation
    RunoffSubsurface = (1.0_dp - SoilImpervFracMax) * BFLOW
    
    RunoffSubsurface = min(0.0002,RunoffSubSurface)
    
    ! Set FSW_change to zero for following calculations in SoilWaterMain and WaterBalanceError Check
    FSW_change = 0.0

    end associate

  end subroutine RunoffSubSurfacePeatland

end module RunoffSubSurfacePeatlandMod
