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
  ! Helper: flooded area fraction AR1 as function of water table depth 
  ! --------------------------------------------------------------------

  pure real function FloodedAreaFrac_AR1(wtd, bf1, bf2, ars1, ars2,ars3) result(ar)
  implicit none

  real, intent(in) :: wtd, bf1, bf2, ars1, ars2, ars3
  real :: catdef

  catdef = max(((wtd + bf2)**2 - 1.0e-20), 0.0) * bf1
  ar = (1.0 + ars1*catdef) / (1.0 + ars2*catdef + ars3*catdef**2)
  ! clamp to [0,1] for safety
  ar = max(0.0, min(1.0, ar))

  end function FloodedAreaFrac_AR1

  ! --------------------------------------------------------------------
  ! Surface water storage [mm] from integrating AR1 from z_ref=1m to current WT depth
  ! Assumptions: 
  ! - WaterTableDepth (wtd) is positive downward [m]
  ! - At z_ref = 1.0 m, flooded area is ~0 and storage is 0
  ! - Storage increases as WT becomes shallower than 1 m
  ! If wtd >= z_ref: storage = 0
  ! If wtd < 0 (WT above surface): adds extra ponded depth (-wtd)*1000 mm
  ! --------------------------------------------------------------------

  subroutine CalcSurfaceWaterStorage_mm(wtd, storage_mm, nsteps)
    implicit none

    real, intent(in) :: wtd  ! water table depth [m], positive downward
    real, intent(out) :: storage_mm  ! surface water storage [mm]
    integer, intent(in), optional :: nsteps  ! number of steps for numerical integration (default: 100)

    ! peatland parameters (same as in MicroTopoCorrection)
    real, parameter :: bf1 = 1.8064707e+02
    real, parameter :: bf2 = 2.4298242e-01
    real, parameter :: ars1 = -8.9250673e-03
    real, parameter :: ars2 = 5.7296452e-02
    real, parameter :: ars3 = 2.2656331e-03

    real, parameter :: zref = 1.0 ! [m] where flooded area ~ 0
    integer :: n, i
    real :: z0, z1, dz
    real :: ar0, ar1, integral_m
    storage_mm = 0.0

    ! number of integration sub-intervals
    n = 200
    if (present(nsteps)) n = max(10, nsteps)

    ! Case 1: WT is deeper than reference => no surface storage
    if (wtd >= zref) then
      storage_mm = 0.0
      return
    end if

    ! We integrate from z0 up to z1=zref with trapezoids, where z is depth [m]
    ! If wtd < 0, integrate from 0 to zref and then add ponded (-wtd) 
    if (wtd < 0.0) then
      z0 = 0.0
    else
      z0 = wtd
    end if
    z1 = zref

    dz = (z1 - z0) / real(n)
    integral_m = 0.0

    ar0 = FloodedAreaFrac_AR1(z0, bf1, bf2, ars1, ars2, ars3)
    do i = 1, n
      ar1 = FloodedAreaFrac_AR1(z0 + real(i)*dz, bf1, bf2, ars1, ars2, ars3)
      integral_m = integral_m + 0.5 * (ar0 + ar1) * dz
      ar0 = ar1
    end do

    storage_mm = integral_m * 1000.0  ! convert m to mm

    ! If WT is above surface, add ponded depth
    if (wtd < 0.0) then
      storage_mm = storage_mm + (-wtd)*1000.0
    end if

    ! clamp to non-negative storage
    storage_mm = max(0.0, storage_mm)

  end subroutine CalcSurfaceWaterStorage_mm


end module MicroTopoCorrectionMod
