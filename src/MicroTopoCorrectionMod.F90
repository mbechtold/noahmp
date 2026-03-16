module MicroTopoCorrectionMod

!!! Microtopography correction for peatlands based on Dettmann & Bechtold (2015)
!!! Replaces empirical PEATCLSM equations with physically-based Gaussian
!!! microtopography distribution.
!!!
!!! Computes:
!!!   - FloodedFraction: fraction of surface below water table (from Gaussian CDF)
!!!   - f_soil: fraction of fluxes going to soil (from Sy decomposition)
!!!   - SySoilLocal / SySurfLocal: local specific yield components
!!!
!!! Introduced by Chakraborty & Bechtold (2025), revised with Dettmann & Bechtold (2015)
!!! microtopography theory by Bechtold (2026)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use PeatMicroTopoMod, only : FloodedFrac, FsoilMicroTopo, &
                                SurfaceWaterStorage,          &
                                InitGaussLegendre, gl_initialized

  implicit none

contains

  subroutine MicroTopoCorrection(noahmp)

! ------------------------ Code history --------------------------------------------------
! Original: Specific yield from WaterTableDepth (Chakraborty & Bechtold, 2025)
! Revised:  Dettmann & Bechtold (2015) Gaussian microtopography physics (Bechtold, 2026)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! Local variables
    real(kind=kind_noahmp) :: z_wt             ! water table in D&B convention (positive up)
    real(kind=kind_noahmp) :: ae, bb, thetas   ! Campbell soil parameters

! --------------------------------------------------------------------
    associate(                                                           &
              SoilMoistureSat     => noahmp%water%param%SoilMoistureSat    ,& ! in, saturated water content [m3/m3]
              SoilMatPotentialSat => noahmp%water%param%SoilMatPotentialSat,& ! in, air-entry potential [m]
              SoilExpCoeffB       => noahmp%water%param%SoilExpCoeffB      ,& ! in, Campbell b exponent
              f_soil              => noahmp%water%state%f_soil             ,& ! out, fraction of flux to soil [-]
              FloodedFraction     => noahmp%water%state%FloodedFraction    ,& ! out, fraction of surface flooded [-]
              WaterTableDepth     => noahmp%water%state%WaterTableDepth     & ! in, water table depth [m] positive downward
             )
! ----------------------------------------------------------------------

    ! Ensure Gauss-Legendre quadrature is initialized
    if (.not. gl_initialized) call InitGaussLegendre()

    ! Campbell soil parameters (use top layer, uniform for peat)
    thetas = SoilMoistureSat(1)
    ae     = abs(SoilMatPotentialSat(1))   ! positive air-entry head [m]
    bb     = SoilExpCoeffB(1)

    ! Convert to D&B convention: z positive upward from mean surface
    z_wt = -WaterTableDepth

    ! Compute flooded fraction from Gaussian CDF of surface elevations
    FloodedFraction = FloodedFrac(z_wt)

    ! Compute f_soil from Sy decomposition (Dettmann & Bechtold 2015, Eq. 2)
    ! f_soil = Sy_soil / (Sy_soil + Sy_surface)
    ! This replaces the old f_soil empirical equation from PEATCLSM
    f_soil = FsoilMicroTopo(z_wt, thetas, ae, bb)

    end associate

  end subroutine MicroTopoCorrection

  ! --------------------------------------------------------------------
  ! Surface water storage [mm] from Dettmann & Bechtold (2015)
  ! Gaussian microtopography theory.
  ! Convenience wrapper for external callers (e.g. LIS TWS diagnostic).
  !
  ! Input:  wtd         - water table depth [m], NoahMP convention (positive downward)
  ! Output: storage_mm  - surface water storage [mm]
  ! --------------------------------------------------------------------
  subroutine CalcSurfaceWaterStorage_mm(wtd, storage_mm)
    implicit none
    real, intent(in)  :: wtd            ! water table depth [m], positive downward
    real, intent(out) :: storage_mm     ! surface water storage [mm]

    ! Convert to D&B convention (z positive upward) and call physics routine
    ! SurfaceWaterStorage returns meters; multiply by 1000 for mm
    if (.not. gl_initialized) call InitGaussLegendre()
    storage_mm = real(SurfaceWaterStorage(-real(wtd, kind_noahmp))) * 1000.0

  end subroutine CalcSurfaceWaterStorage_mm

end module MicroTopoCorrectionMod
