module PeatMicroTopoMod

!!! Peatland microtopography module based on Dettmann & Bechtold (2015)
!!! Hydrological Processes, DOI: 10.1002/hyp.10475
!!!
!!! Provides:
!!!   - Gaussian CDF of surface elevations (flooded fraction)
!!!   - Surface water storage in hollows
!!!   - Effective soil water storage accounting for microtopography
!!!   - Specific yield decomposition into soil + surface components
!!!   - Gauss-Legendre quadrature for microtopography integrals
!!!
!!! Sign convention:
!!!   z is positive upward, z=0 at mean ground surface
!!!   z_wt = -WaterTableDepth (NoahMP WaterTableDepth is positive downward)
!!!   z_wt > 0 means flooding above mean surface
!!!   z_wt < 0 means water table below mean surface
!!!
!!! Microtopography:
!!!   Surface elevation ~ N(0, sigma_elev^2), truncated at [-z_trunc, +z_trunc]
!!!   sigma_elev = 0.16 m (standard deviation)
!!!   z_trunc    = 1.0  m (truncation limit)
!!!
!!! Soil hydraulic model: Campbell (as in NoahMP)
!!!   theta(h) = theta_s                          for h >= -h_e
!!!   theta(h) = theta_s * (-h/h_e)^(-1/b)       for h <  -h_e
!!!   where h = pressure head (negative in unsaturated zone, h = z_wt - z for hydrostatic)
!!!         h_e = air-entry suction head (positive value, = abs(SoilMatPotentialSat))
!!!         b   = Campbell exponent (SoilExpCoeffB)
!!!
!!! Reference: Dettmann & Bechtold (2015), Hydrological Processes

  use Machine
  use ConstantDefineMod

  implicit none

  ! ------ Module parameters ------
  real(kind=kind_noahmp), parameter :: sigma_elev  = 0.16_kind_noahmp   ! std dev of microtopography [m]
  real(kind=kind_noahmp), parameter :: z_trunc     = 1.0_kind_noahmp    ! truncation limit [m]
  real(kind=kind_noahmp), parameter :: pi_noahmp   = 3.14159265358979323846_kind_noahmp

  ! ------ Gauss-Legendre quadrature (20 points on [-1,1]) ------
  integer, parameter :: n_gl = 20
  real(kind=kind_noahmp), dimension(n_gl) :: gl_nodes, gl_weights

  ! Flag for initialization
  logical, save :: gl_initialized = .false.

contains

  !========================================================================
  ! Initialize Gauss-Legendre nodes and weights for 20-point quadrature
  !========================================================================
  subroutine InitGaussLegendre()
    implicit none

    ! 20-point Gauss-Legendre nodes (on [-1,1]) and weights
    ! Tabulated to high precision
    gl_nodes( 1) = -0.9931285991850949_kind_noahmp
    gl_nodes( 2) = -0.9639719272779138_kind_noahmp
    gl_nodes( 3) = -0.9122344282513259_kind_noahmp
    gl_nodes( 4) = -0.8391169718222188_kind_noahmp
    gl_nodes( 5) = -0.7463319064601508_kind_noahmp
    gl_nodes( 6) = -0.6360536807265150_kind_noahmp
    gl_nodes( 7) = -0.5108670019508271_kind_noahmp
    gl_nodes( 8) = -0.3737060887154195_kind_noahmp
    gl_nodes( 9) = -0.2277858511416451_kind_noahmp
    gl_nodes(10) = -0.0765265211334973_kind_noahmp
    gl_nodes(11) =  0.0765265211334973_kind_noahmp
    gl_nodes(12) =  0.2277858511416451_kind_noahmp
    gl_nodes(13) =  0.3737060887154195_kind_noahmp
    gl_nodes(14) =  0.5108670019508271_kind_noahmp
    gl_nodes(15) =  0.6360536807265150_kind_noahmp
    gl_nodes(16) =  0.7463319064601508_kind_noahmp
    gl_nodes(17) =  0.8391169718222188_kind_noahmp
    gl_nodes(18) =  0.9122344282513259_kind_noahmp
    gl_nodes(19) =  0.9639719272779138_kind_noahmp
    gl_nodes(20) =  0.9931285991850949_kind_noahmp

    gl_weights( 1) = 0.0176140071391521_kind_noahmp
    gl_weights( 2) = 0.0406014298003869_kind_noahmp
    gl_weights( 3) = 0.0626720483341091_kind_noahmp
    gl_weights( 4) = 0.0832767415767048_kind_noahmp
    gl_weights( 5) = 0.1019301198172404_kind_noahmp
    gl_weights( 6) = 0.1181945319615184_kind_noahmp
    gl_weights( 7) = 0.1316886384491766_kind_noahmp
    gl_weights( 8) = 0.1420961093183820_kind_noahmp
    gl_weights( 9) = 0.1491729864726037_kind_noahmp
    gl_weights(10) = 0.1527533871307258_kind_noahmp
    gl_weights(11) = 0.1527533871307258_kind_noahmp
    gl_weights(12) = 0.1491729864726037_kind_noahmp
    gl_weights(13) = 0.1420961093183820_kind_noahmp
    gl_weights(14) = 0.1316886384491766_kind_noahmp
    gl_weights(15) = 0.1181945319615184_kind_noahmp
    gl_weights(16) = 0.1019301198172404_kind_noahmp
    gl_weights(17) = 0.0832767415767048_kind_noahmp
    gl_weights(18) = 0.0626720483341091_kind_noahmp
    gl_weights(19) = 0.0406014298003869_kind_noahmp
    gl_weights(20) = 0.0176140071391521_kind_noahmp

    gl_initialized = .true.

  end subroutine InitGaussLegendre

  !========================================================================
  ! Approximate error function using Abramowitz & Stegun (1964) 7.1.26
  ! Maximum error: 1.5e-7
  !========================================================================
  pure function erf_approx(x) result(res)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x
    real(kind=kind_noahmp) :: res
    real(kind=kind_noahmp) :: t, ax, poly
    real(kind=kind_noahmp), parameter :: p  =  0.3275911_kind_noahmp
    real(kind=kind_noahmp), parameter :: a1 =  0.254829592_kind_noahmp
    real(kind=kind_noahmp), parameter :: a2 = -0.284496736_kind_noahmp
    real(kind=kind_noahmp), parameter :: a3 =  1.421413741_kind_noahmp
    real(kind=kind_noahmp), parameter :: a4 = -1.453152027_kind_noahmp
    real(kind=kind_noahmp), parameter :: a5 =  1.061405429_kind_noahmp

    ax = abs(x)
    t  = 1.0_kind_noahmp / (1.0_kind_noahmp + p * ax)
    poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))
    res = 1.0_kind_noahmp - poly * exp(-ax * ax)
    if (x < 0.0_kind_noahmp) res = -res
  end function erf_approx

  !========================================================================
  ! Standard normal CDF: Phi(x) = 0.5 * (1 + erf(x/sqrt(2)))
  !========================================================================
  pure function phi_normal(x) result(res)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x
    real(kind=kind_noahmp) :: res
    res = 0.5_kind_noahmp * (1.0_kind_noahmp + erf_approx(x / sqrt(2.0_kind_noahmp)))
  end function phi_normal

  !========================================================================
  ! Standard normal PDF: phi_pdf(x) = exp(-x^2/2) / sqrt(2*pi)
  !========================================================================
  pure function phi_pdf(x) result(res)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x
    real(kind=kind_noahmp) :: res
    res = exp(-0.5_kind_noahmp * x * x) / sqrt(2.0_kind_noahmp * pi_noahmp)
  end function phi_pdf

  !========================================================================
  ! Cumulative distribution function of surface elevations: F_s(z)
  ! F_s(z) = CDF of N(0, sigma_elev^2), clamped at z_trunc
  ! Returns fraction of ground surface BELOW elevation z
  !========================================================================
  pure function Fs_cdf(z) result(res)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z
    real(kind=kind_noahmp) :: res

    if (z <= -z_trunc) then
       res = 0.0_kind_noahmp
    else if (z >= z_trunc) then
       res = 1.0_kind_noahmp
    else
       res = phi_normal(z / sigma_elev)
    endif
  end function Fs_cdf

  !========================================================================
  ! Flooded fraction: fraction of surface below water table
  ! z_wt: water table position (positive upward from mean surface)
  !========================================================================
  pure function FloodedFrac(z_wt) result(res)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt
    real(kind=kind_noahmp) :: res
    res = Fs_cdf(z_wt)
  end function FloodedFrac

  !========================================================================
  ! Campbell soil water content at given pressure head h
  ! h < 0: unsaturated; h >= 0: saturated
  ! h_e: air-entry suction (positive), b: Campbell exponent
  ! theta_s: saturated water content
  !========================================================================
  pure function theta_campbell(h, theta_s, h_e, b) result(theta)
    implicit none
    real(kind=kind_noahmp), intent(in) :: h, theta_s, h_e, b
    real(kind=kind_noahmp) :: theta

    if (h >= -h_e) then
       ! Saturated (pressure head > air entry)
       theta = theta_s
    else
       ! Unsaturated: theta = theta_s * (|h|/h_e)^(-1/b)
       theta = theta_s * (abs(h) / h_e) ** (-1.0_kind_noahmp / b)
       theta = max(0.0_kind_noahmp, min(theta_s, theta))
    endif
  end function theta_campbell

  !========================================================================
  ! Surface water storage volume [m] for water table at z_wt
  ! Integral of F_s(z) dz from z_min to z_wt  (Eq. 3 in D&B 2015)
  ! For z_wt <= -z_trunc: no surface water
  !========================================================================
  function SurfaceWaterStorage(z_wt) result(vol)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt
    real(kind=kind_noahmp) :: vol
    real(kind=kind_noahmp) :: z_lo, z_hi, z_mid, z_half, z_pt
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    if (z_wt <= -z_trunc) then
       vol = 0.0_kind_noahmp
       return
    endif

    ! Integration limits
    z_lo = -z_trunc
    z_hi = min(z_wt, z_trunc)

    ! Gauss-Legendre integration of F_s(z) from z_lo to z_hi
    z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
    z_half = 0.5_kind_noahmp * (z_hi - z_lo)

    vol = 0.0_kind_noahmp
    do k = 1, n_gl
       z_pt = z_mid + z_half * gl_nodes(k)
       vol  = vol + gl_weights(k) * Fs_cdf(z_pt)
    enddo
    vol = vol * z_half

    ! If z_wt > z_trunc, all surface above z_trunc is open water (F_s=1)
    if (z_wt > z_trunc) then
       vol = vol + (z_wt - z_trunc)
    endif

  end function SurfaceWaterStorage

  !========================================================================
  ! Soil water storage [m] for water table at z_wt, integrated over
  ! microtopography (Eq. 4 in D&B 2015)
  !
  ! A_soil(z_wt) = integral from -z_trunc to +z_trunc of
  !                (1 - F_s(z)) * theta(z_wt - z) dz
  !
  ! where theta is the Campbell retention with hydrostatic assumption
  ! h = z_wt - z (pressure head at elevation z for water table at z_wt)
  !
  ! This gives the total water stored in the soil column per unit area
  ! from the lowest point of the microtopography upward.
  !========================================================================
  function SoilWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp) result(A_soil)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: A_soil
    real(kind=kind_noahmp) :: z_lo, z_hi, z_mid, z_half, z_pt, h_pt
    real(kind=kind_noahmp) :: soil_frac, theta_val
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    ! Integration from -z_trunc to +z_trunc
    ! (1-F_s(z))*theta(z_wt-z) is zero for z > z_trunc and F_s=1 beyond z_trunc
    z_lo = -z_trunc
    z_hi =  z_trunc

    z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
    z_half = 0.5_kind_noahmp * (z_hi - z_lo)

    A_soil = 0.0_kind_noahmp
    do k = 1, n_gl
       z_pt = z_mid + z_half * gl_nodes(k)
       soil_frac = 1.0_kind_noahmp - Fs_cdf(z_pt)
       h_pt = z_wt - z_pt   ! pressure head at elevation z_pt
       theta_val = theta_campbell(h_pt, theta_s, h_e, b_camp)
       A_soil = A_soil + gl_weights(k) * soil_frac * theta_val
    enddo
    A_soil = A_soil * z_half

  end function SoilWaterStorageMicroTopo

  !========================================================================
  ! Total water storage [m] = soil + surface for water table at z_wt
  !========================================================================
  function TotalWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp) result(W_total)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: W_total

    W_total = SoilWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp) + &
              SurfaceWaterStorage(z_wt)

  end function TotalWaterStorageMicroTopo

  !========================================================================
  ! Effective soil moisture for a layer [z_bot, z_top] (both relative to
  ! mean surface, z positive upward) accounting for microtopography
  !
  ! theta_eff = (1 / layer_thickness) * integral from z_bot to z_top of
  !             (1 - F_s(z)) * theta(z_wt - z) dz
  !            / (1 / layer_thickness) * integral from z_bot to z_top of
  !             (1 - F_s(z)) dz
  !
  ! Actually, for comparison with flat-surface layer moisture, we want:
  ! theta_eff = integral of (1-Fs(z))*theta(zwt-z) dz / integral of (1-Fs(z)) dz
  ! over the layer range, so it represents average moisture of the soil
  ! fraction only.
  !
  ! But for total water [m] in the layer (soil + surface), we return:
  ! W_layer = integral of [(1-Fs(z))*theta(zwt-z) + Fs(z)*I(z<zwt)] dz
  ! where I(z<zwt) = 1 when z < zwt (surface water present)
  !
  ! For NoahMP compatibility, we return:
  ! theta_eff_flat = W_layer / layer_thickness
  ! This is what should be compared with satellite top-10cm data.
  !========================================================================
  function EffectiveSoilMoistureLayer(z_bot, z_top, z_wt, theta_s, h_e, b_camp) result(theta_eff)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_bot, z_top, z_wt, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: theta_eff
    real(kind=kind_noahmp) :: layer_thick, z_lo, z_hi, z_mid, z_half, z_pt
    real(kind=kind_noahmp) :: h_pt, soil_frac, theta_val, w_soil, w_surface
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    layer_thick = z_top - z_bot
    if (layer_thick <= 0.0_kind_noahmp) then
       theta_eff = 0.0_kind_noahmp
       return
    endif

    ! Clamp integration limits to microtopography range
    z_lo = max(z_bot, -z_trunc)
    z_hi = min(z_top,  z_trunc)

    if (z_lo >= z_hi) then
       ! Layer entirely outside microtopography range
       if (z_bot >= z_trunc) then
          ! Above all hummocks: no soil, just surface water if z_wt > z_top
          if (z_wt > z_bot) then
             theta_eff = 1.0_kind_noahmp  ! open water
          else
             theta_eff = 0.0_kind_noahmp
          endif
       else
          ! Below all hollows: all soil
          h_pt = z_wt - 0.5_kind_noahmp*(z_bot+z_top)
          theta_eff = theta_campbell(h_pt, theta_s, h_e, b_camp)
       endif
       return
    endif

    z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
    z_half = 0.5_kind_noahmp * (z_hi - z_lo)

    w_soil    = 0.0_kind_noahmp
    w_surface = 0.0_kind_noahmp

    do k = 1, n_gl
       z_pt = z_mid + z_half * gl_nodes(k)
       soil_frac = 1.0_kind_noahmp - Fs_cdf(z_pt)
       h_pt = z_wt - z_pt
       theta_val = theta_campbell(h_pt, theta_s, h_e, b_camp)

       ! Soil water contribution
       w_soil = w_soil + gl_weights(k) * soil_frac * theta_val

       ! Surface water contribution: where z_pt < z_wt and no soil
       if (z_pt < z_wt) then
          w_surface = w_surface + gl_weights(k) * Fs_cdf(z_pt)
       endif
    enddo
    w_soil    = w_soil * z_half
    w_surface = w_surface * z_half

    ! Add contributions from parts of the layer below -z_trunc (all soil)
    if (z_bot < -z_trunc) then
       ! From z_bot to -z_trunc: all soil (F_s = 0)
       z_mid  = 0.5_kind_noahmp * (-z_trunc + z_bot)
       z_half = 0.5_kind_noahmp * (-z_trunc - z_bot)
       do k = 1, n_gl
          z_pt = z_mid + z_half * gl_nodes(k)
          h_pt = z_wt - z_pt
          theta_val = theta_campbell(h_pt, theta_s, h_e, b_camp)
          w_soil = w_soil + gl_weights(k) * theta_val * z_half
       enddo
    endif

    ! Add contributions from parts of the layer above +z_trunc (no soil, F_s=1)
    if (z_top > z_trunc) then
       ! From z_trunc to z_top: no soil, all surface water if z_wt > z
       if (z_wt > z_trunc) then
          w_surface = w_surface + min(z_wt, z_top) - z_trunc
       endif
    endif

    ! Total water in layer [m] divided by layer thickness -> effective theta
    theta_eff = (w_soil + w_surface) / layer_thick

    ! Clamp for safety
    theta_eff = max(0.0_kind_noahmp, min(1.0_kind_noahmp, theta_eff))

  end function EffectiveSoilMoistureLayer

  !========================================================================
  ! Specific yield of soil component (Eq. 5 / Eq. 6 in D&B 2015)
  ! Sy_soil for water level change from z_l to z_u
  !
  ! Sy_soil = (1/dz) * integral from -z_trunc to +z_trunc of
  !           (1-Fs(z)) * [theta(z_u-z) - theta(z_l-z)] dz
  !========================================================================
  function SysoilMicroTopo(z_l, z_u, theta_s, h_e, b_camp) result(sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_l, z_u, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: sy
    real(kind=kind_noahmp) :: dz, z_lo, z_hi, z_mid, z_half, z_pt
    real(kind=kind_noahmp) :: soil_frac, theta_u, theta_l
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    dz = z_u - z_l
    if (abs(dz) < 1.0e-10_kind_noahmp) then
       sy = 0.0_kind_noahmp
       return
    endif

    z_lo = -z_trunc
    z_hi =  z_trunc

    z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
    z_half = 0.5_kind_noahmp * (z_hi - z_lo)

    sy = 0.0_kind_noahmp
    do k = 1, n_gl
       z_pt = z_mid + z_half * gl_nodes(k)
       soil_frac = 1.0_kind_noahmp - Fs_cdf(z_pt)
       theta_u = theta_campbell(z_u - z_pt, theta_s, h_e, b_camp)
       theta_l = theta_campbell(z_l - z_pt, theta_s, h_e, b_camp)
       sy = sy + gl_weights(k) * soil_frac * (theta_u - theta_l)
    enddo
    sy = sy * z_half / dz

    sy = max(0.0_kind_noahmp, sy)

  end function SysoilMicroTopo

  !========================================================================
  ! Specific yield of surface component (Eq. 3 in D&B 2015)
  ! Sy_surface for water level change from z_l to z_u
  !
  ! Sy_surface = (1/dz) * integral from z_l to z_u of F_s(z) dz
  !========================================================================
  function SysurfaceMicroTopo(z_l, z_u) result(sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_l, z_u
    real(kind=kind_noahmp) :: sy
    real(kind=kind_noahmp) :: dz, z_lo, z_hi, z_mid, z_half, z_pt
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    dz = z_u - z_l
    if (abs(dz) < 1.0e-10_kind_noahmp) then
       sy = 0.0_kind_noahmp
       return
    endif

    ! Integration limits clamped to relevant range
    z_lo = max(z_l, -z_trunc)
    z_hi = min(z_u,  z_trunc)

    sy = 0.0_kind_noahmp

    ! Below -z_trunc: Fs=0, no contribution
    ! Within [-z_trunc, +z_trunc]: integrate Fs
    if (z_lo < z_hi) then
       z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
       z_half = 0.5_kind_noahmp * (z_hi - z_lo)
       do k = 1, n_gl
          z_pt = z_mid + z_half * gl_nodes(k)
          sy = sy + gl_weights(k) * Fs_cdf(z_pt)
       enddo
       sy = sy * z_half
    endif

    ! Above +z_trunc: Fs=1, contribution = min(z_u,.) - z_trunc
    if (z_u > z_trunc) then
       sy = sy + (z_u - max(z_trunc, z_l))
    endif

    sy = sy / dz
    sy = max(0.0_kind_noahmp, min(1.0_kind_noahmp, sy))

  end function SysurfaceMicroTopo

  !========================================================================
  ! Total specific yield Sy = Sy_soil + Sy_surface (Eq. 2 in D&B 2015)
  !========================================================================
  function SytotalMicroTopo(z_l, z_u, theta_s, h_e, b_camp) result(sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_l, z_u, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: sy

    sy = SysoilMicroTopo(z_l, z_u, theta_s, h_e, b_camp) + &
         SysurfaceMicroTopo(z_l, z_u)

    sy = max(1.0e-6_kind_noahmp, min(1.0_kind_noahmp, sy))

  end function SytotalMicroTopo

  !========================================================================
  ! Fraction of flux going to soil (f_soil) based on Sy decomposition
  ! f_soil = Sy_soil / Sy_total
  ! Uses a small dz increment around current water table to compute
  ! local specific yield
  !========================================================================
  function FsoilMicroTopo(z_wt, theta_s, h_e, b_camp) result(f_s)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: f_s
    real(kind=kind_noahmp) :: dz_local, z_l, z_u
    real(kind=kind_noahmp) :: sy_soil, sy_surf, sy_tot

    ! Use a small increment to compute local Sy
    dz_local = 0.005_kind_noahmp  ! 5 mm increment
    z_l = z_wt - dz_local
    z_u = z_wt + dz_local

    sy_soil = SysoilMicroTopo(z_l, z_u, theta_s, h_e, b_camp)
    sy_surf = SysurfaceMicroTopo(z_l, z_u)
    sy_tot  = sy_soil + sy_surf

    if (sy_tot > 1.0e-10_kind_noahmp) then
       f_s = sy_soil / sy_tot
    else
       f_s = 1.0_kind_noahmp  ! deep water table: all flux to soil
    endif

    f_s = max(0.0_kind_noahmp, min(1.0_kind_noahmp, f_s))

  end function FsoilMicroTopo

  !========================================================================
  ! Water table depth from total water deficit
  ! Given current soil moisture profile and microtopography, find the
  ! equilibrium water table position z_wt such that the total water
  ! storage (soil + surface) matches the given total water content.
  !
  ! Uses bisection method on the monotonically increasing relationship
  ! between z_wt and total water storage.
  !
  ! total_water_content: total water stored [m] in soil + surface
  ! theta_s, h_e, b_camp: Campbell soil parameters
  !========================================================================
  function FindWaterTable(total_water_content, theta_s, h_e, b_camp) result(z_wt)
    implicit none
    real(kind=kind_noahmp), intent(in) :: total_water_content, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: z_wt
    real(kind=kind_noahmp) :: z_lo, z_hi, z_mid, W_lo, W_hi, W_mid
    integer :: iter
    integer, parameter :: max_iter = 60
    real(kind=kind_noahmp), parameter :: tol = 1.0e-5_kind_noahmp

    if (.not. gl_initialized) call InitGaussLegendre()

    ! Bracket: z_wt can range from very deep to above surface
    z_lo = -3.0_kind_noahmp   ! 3 m below surface
    z_hi =  z_trunc           ! at truncation limit (1 m above)

    W_lo = TotalWaterStorageMicroTopo(z_lo, theta_s, h_e, b_camp)
    W_hi = TotalWaterStorageMicroTopo(z_hi, theta_s, h_e, b_camp)

    ! Check if target is within range
    if (total_water_content <= W_lo) then
       z_wt = z_lo
       return
    endif
    if (total_water_content >= W_hi) then
       z_wt = z_hi
       return
    endif

    ! Bisection
    do iter = 1, max_iter
       z_mid = 0.5_kind_noahmp * (z_lo + z_hi)
       W_mid = TotalWaterStorageMicroTopo(z_mid, theta_s, h_e, b_camp)

       if (abs(W_mid - total_water_content) < tol .or. (z_hi - z_lo) < tol) then
          z_wt = z_mid
          return
       endif

       if (W_mid < total_water_content) then
          z_lo = z_mid
       else
          z_hi = z_mid
       endif
    enddo

    z_wt = 0.5_kind_noahmp * (z_lo + z_hi)

  end function FindWaterTable

  !========================================================================
  ! Compute the soil-only water storage (deficit relative to full 
  ! saturation) for a given water table, including microtopography effect.
  ! This is used for the WaterTableEquilibrium calculation.
  !
  ! deficit = integral of (theta_s - theta(z)) * (1-Fs(z)) dz
  !           over the full soil column
  ! We compute it as: theta_s * total_soil_volume - A_soil(z_wt)
  !========================================================================
  function SoilDeficitMicroTopo(z_wt, theta_s, h_e, b_camp) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: deficit
    real(kind=kind_noahmp) :: total_soil_vol, A_soil
    real(kind=kind_noahmp) :: z_lo, z_hi, z_mid, z_half, z_pt, soil_frac
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    ! Total soil volume per unit area = integral of (1-Fs(z)) dz from -z_trunc to z_trunc
    ! For N(0,sigma^2): integral of (1-Phi(z/sigma)) dz from -L to L
    !   = integral of Phi(-z/sigma) dz from -L to L  (by symmetry)
    !   = L  (exactly, by symmetry of the standard normal)
    ! Actually: integral from -L to L of (1-Fs(z)) dz = L (half the range, by symmetry)
    ! More precisely: integral of (1-Phi(z/sigma)) dz from -infinity to infinity = 0 
    ! (since the mean of the distribution is at z=0)
    ! No -- let me compute it numerically for correctness.
    
    z_lo = -z_trunc
    z_hi =  z_trunc
    z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
    z_half = 0.5_kind_noahmp * (z_hi - z_lo)

    total_soil_vol = 0.0_kind_noahmp
    do k = 1, n_gl
       z_pt = z_mid + z_half * gl_nodes(k)
       soil_frac = 1.0_kind_noahmp - Fs_cdf(z_pt)
       total_soil_vol = total_soil_vol + gl_weights(k) * soil_frac
    enddo
    total_soil_vol = total_soil_vol * z_half

    A_soil = SoilWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp)

    deficit = theta_s * total_soil_vol - A_soil
    deficit = max(0.0_kind_noahmp, deficit)

  end function SoilDeficitMicroTopo

end module PeatMicroTopoMod
