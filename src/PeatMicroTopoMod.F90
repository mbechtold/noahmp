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
!!!   Surface elevation ~ N(0, sigma_elev^2), truncated at +z_trunc only
!!!   sigma_elev = 0.16 m (standard deviation)
!!!   z_trunc    = 1.0  m (upper truncation limit, top of hummocks)
!!!   No lower truncation: Fs_cdf evaluates the Gaussian CDF at any depth.
!!!   At depth −2 m (z/σ ≈ −12.5), Fs ≈ 10^−35, so soil_frac ≈ 1 naturally.
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

  ! ------ Gauss-Legendre quadrature ------
  ! Standard 40-point rule for functions with smooth integrands
  ! (surface water, specific yield, soil thickness, FindWaterTable, etc.)
  integer, parameter :: n_gl = 40
  real(kind=kind_noahmp), dimension(n_gl) :: gl_nodes, gl_weights

  ! 20-point rule for equilibrium SM integrals
  ! (used in EquilibriumSMMicroTopo; resolves the Campbell air-entry kink)
  integer, parameter :: n_gl_lite = 20
  real(kind=kind_noahmp), dimension(n_gl_lite) :: gl_nodes_lite, gl_weights_lite

  ! Fine 400-point rule for the microtopo zone in SoilWaterStorageMicroTopo
  ! where the Campbell air-entry kink demands high resolution
  integer, parameter :: n_gl_fine = 400
  real(kind=kind_noahmp), dimension(n_gl_fine) :: gl_nodes_fine, gl_weights_fine

  ! Precomputed quadrature coefficients for microtopo-zone integrals
  ! wt_soil_micro(k) = gl_weights(k) * (1 - Fs_cdf(z_trunc * gl_nodes(k)))
  ! Eliminates repeated Fs_cdf/erf_approx calls in hot GL loops
  real(kind=kind_noahmp), dimension(n_gl) :: wt_soil_micro

  ! Flag for initialization
  logical, save :: gl_initialized = .false.

contains

  !========================================================================
  ! Initialize Gauss-Legendre nodes and weights
  ! All rules computed via Newton iteration on Legendre polynomials.
  ! - 40-point rule: general use (Sy, surface water, FindWaterTable, etc.)
  ! - 20-point rule: equilibrium SM integrals (EquilibriumSMMicroTopo)
  ! - 400-point rule: fine diagnostic integrals (SoilWaterStorageMicroTopo)
  !========================================================================
  subroutine InitGaussLegendre()
    implicit none
    integer :: k

    ! Compute all quadrature rules via Newton iteration
    call ComputeGaussLegendre(n_gl,      gl_nodes,      gl_weights)
    call ComputeGaussLegendre(n_gl_lite, gl_nodes_lite, gl_weights_lite)
    call ComputeGaussLegendre(n_gl_fine, gl_nodes_fine, gl_weights_fine)

    ! Precompute weighted soil fractions at microtopo-zone GL nodes
    ! Physical z = z_trunc * gl_nodes(k) for domain [-z_trunc, +z_trunc]
    do k = 1, n_gl
       wt_soil_micro(k) = gl_weights(k) * &
           (1.0_kind_noahmp - Fs_cdf(z_trunc * gl_nodes(k)))
    enddo

    gl_initialized = .true.

  end subroutine InitGaussLegendre

  !========================================================================
  ! Compute n-point Gauss-Legendre nodes and weights on [-1,1]
  ! using Newton iteration on the Legendre polynomial P_n(x).
  ! Exploits symmetry: only computes half the roots.
  !========================================================================
  subroutine ComputeGaussLegendre(n, nodes, weights)
    implicit none
    integer, intent(in) :: n
    real(kind=kind_noahmp), intent(out) :: nodes(n), weights(n)
    integer :: i, j, m, iter
    real(kind=kind_noahmp) :: x, x_old, p0, p1, p2, dp
    integer, parameter :: max_iter = 100
    real(kind=kind_noahmp), parameter :: tol = 1.0e-15_kind_noahmp

    m = (n + 1) / 2  ! number of roots to compute (positive half + centre)

    do i = 1, m
       ! Initial guess (Tricomi approximation)
       x = cos(pi_noahmp * (real(i, kind_noahmp) - 0.25_kind_noahmp) / &
                            (real(n, kind_noahmp) + 0.5_kind_noahmp))

       ! Newton iteration to find root of P_n(x)
       do iter = 1, max_iter
          ! Evaluate P_n(x) and P_n'(x) via recurrence
          p0 = 1.0_kind_noahmp          ! P_0(x)
          p1 = x                         ! P_1(x)
          do j = 2, n
             p2 = ((2.0_kind_noahmp * real(j, kind_noahmp) - 1.0_kind_noahmp) * x * p1 &
                  - (real(j, kind_noahmp) - 1.0_kind_noahmp) * p0) / real(j, kind_noahmp)
             p0 = p1
             p1 = p2
          enddo
          ! p1 = P_n(x)
          ! Derivative: P_n'(x) = n*(x*P_n - P_{n-1}) / (x^2 - 1)
          dp = real(n, kind_noahmp) * (x * p1 - p0) / (x * x - 1.0_kind_noahmp)
          x_old = x
          x = x - p1 / dp
          if (abs(x - x_old) < tol) exit
       enddo

       ! Assign symmetric pairs
       nodes(i)     = -x
       nodes(n+1-i) =  x
       weights(i)       = 2.0_kind_noahmp / ((1.0_kind_noahmp - x*x) * dp*dp)
       weights(n+1-i)   = weights(i)
    enddo

  end subroutine ComputeGaussLegendre

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
  ! F_s(z) = CDF of N(0, sigma_elev^2), truncated at +z_trunc only.
  ! Returns fraction of ground surface BELOW elevation z.
  !
  ! No lower truncation: the Gaussian CDF is evaluated for all z < z_trunc.
  ! At depth −2 m (z/σ ≈ −12.5), Φ ≈ 10^−35, so soil_frac = 1−Fs ≈ 1
  ! naturally without requiring an artificial clamp.
  !========================================================================
  pure function Fs_cdf(z) result(res)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z
    real(kind=kind_noahmp) :: res

    if (z >= z_trunc) then
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
  ! For z_wt <= -z_trunc: Fs is negligible, so surface water ≈ 0.
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
  ! microtopography (Eq. 4 in D&B 2015, extended to full soil column)
  !
  ! A_soil(z_wt) = integral from -z_col_bot to +z_trunc of
  !                (1 - F_s(z)) * theta(z_wt - z) dz
  !
  ! where theta is the Campbell retention with hydrostatic assumption
  ! h = z_wt - z (pressure head at elevation z for water table at z_wt)
  !
  ! z_col_bot:  soil column bottom depth [m], positive downward from
  !             mean surface.  The integration extends from elevation
  !             -z_col_bot up to +z_trunc (top of hummocks).
  !             Deep below the surface, Fs_cdf ≈ 0 so soil_frac ≈ 1
  !             naturally from the Gaussian CDF tail.
  !
  ! Composite Gauss-Legendre quadrature: the domain is split at
  ! z = -z_trunc into two sub-intervals to maintain resolution.
  ! The Campbell air-entry transition (width h_e ≈ 0.024 m) would
  ! create derivative discontinuities at each quadrature node when a
  ! single GL rule spans the full 3 m domain.  Splitting concentrates
  ! 20 nodes on the 2 m microtopo zone and 20 on the 1 m deep zone,
  ! yielding much smoother integrals as z_wt varies.
  !========================================================================
  function SoilWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp, z_col_bot) result(A_soil)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: A_soil
    real(kind=kind_noahmp) :: z_lo, z_hi, z_mid, z_half, z_pt, h_pt
    real(kind=kind_noahmp) :: soil_frac, theta_val, A_sub
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    A_soil = 0.0_kind_noahmp

    ! --- Sub-interval 1: deep zone [-z_col_bot, -z_trunc] ---
    ! soil_frac ≈ 1 here (Fs negligible at depth).
    ! Included for generality; typically fully saturated for peatland WTDs.
    if (z_col_bot > z_trunc) then
      z_lo = -z_col_bot
      z_hi = -z_trunc
      z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
      z_half = 0.5_kind_noahmp * (z_hi - z_lo)

      A_sub = 0.0_kind_noahmp
      do k = 1, n_gl
         z_pt = z_mid + z_half * gl_nodes(k)
         soil_frac = 1.0_kind_noahmp - Fs_cdf(z_pt)
         h_pt = z_wt - z_pt
         theta_val = theta_campbell(h_pt, theta_s, h_e, b_camp)
         A_sub = A_sub + gl_weights(k) * soil_frac * theta_val
      enddo
      A_soil = A_soil + A_sub * z_half
    endif

    ! --- Sub-interval 2: microtopo zone [-z_trunc, +z_trunc] ---
    ! Full Gaussian CDF variation of soil_frac occurs here.
    ! Uses 100-point GL rule for smooth deficit-vs-WTD relationship.
    z_lo = -z_trunc
    z_hi =  z_trunc
    z_mid  = 0.5_kind_noahmp * (z_hi + z_lo)
    z_half = 0.5_kind_noahmp * (z_hi - z_lo)

    A_sub = 0.0_kind_noahmp
    do k = 1, n_gl_fine
       z_pt = z_mid + z_half * gl_nodes_fine(k)
       soil_frac = 1.0_kind_noahmp - Fs_cdf(z_pt)
       h_pt = z_wt - z_pt
       theta_val = theta_campbell(h_pt, theta_s, h_e, b_camp)
       A_sub = A_sub + gl_weights_fine(k) * soil_frac * theta_val
    enddo
    A_soil = A_soil + A_sub * z_half

  end function SoilWaterStorageMicroTopo

  !========================================================================
  ! Lightweight soil water storage [m] for water table at z_wt.
  ! Same physics as SoilWaterStorageMicroTopo but uses n_gl (10-pt)
  ! instead of n_gl_fine (100-pt) for the microtopo zone.
  !
  ! Designed for fast repeated evaluation inside FindWaterTable bisection.
  ! Accuracy is sufficient for WTD diagnosis (tolerance 1e-5 m).
  !========================================================================
  function SoilWaterStorageMicroTopoLite(z_wt, theta_s, h_e, b_camp, z_col_bot) result(A_soil)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: A_soil
    real(kind=kind_noahmp) :: h_pt, inv_b, A_sub, z_mid_d, z_half_d
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    A_soil = 0.0_kind_noahmp
    inv_b  = -1.0_kind_noahmp / b_camp

    ! --- Deep zone [-z_col_bot, -z_trunc]: soil_frac ≈ 1 (Fs negligible) ---
    if (z_col_bot > z_trunc) then
       if (z_wt >= -(z_trunc + h_e)) then
          ! Entire deep zone saturated (common case for peatland WTDs)
          A_soil = theta_s * (z_col_bot - z_trunc)
       else
          ! Partially unsaturated — GL without Fs_cdf (soil_frac = 1)
          z_mid_d  = -0.5_kind_noahmp * (z_col_bot + z_trunc)
          z_half_d =  0.5_kind_noahmp * (z_col_bot - z_trunc)
          A_sub = 0.0_kind_noahmp
          do k = 1, n_gl
             h_pt = z_wt - (z_mid_d + z_half_d * gl_nodes(k))
             if (h_pt >= -h_e) then
                A_sub = A_sub + gl_weights(k) * theta_s
             else
                A_sub = A_sub + gl_weights(k) * theta_s * (abs(h_pt) / h_e) ** inv_b
             endif
          enddo
          A_soil = A_sub * z_half_d
       endif
    endif

    ! --- Microtopo zone [-z_trunc, +z_trunc]: precomputed soil fractions ---
    A_sub = 0.0_kind_noahmp
    do k = 1, n_gl
       h_pt = z_wt - z_trunc * gl_nodes(k)
       if (h_pt >= -h_e) then
          A_sub = A_sub + wt_soil_micro(k) * theta_s
       else
          A_sub = A_sub + wt_soil_micro(k) * theta_s * (abs(h_pt) / h_e) ** inv_b
       endif
    enddo
    A_soil = A_soil + A_sub * z_trunc

  end function SoilWaterStorageMicroTopoLite

  !========================================================================
  ! Total water storage [m] = soil + surface for water table at z_wt
  ! z_col_bot: soil column bottom depth [m], positive downward
  !========================================================================
  function TotalWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp, z_col_bot) result(W_total)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: W_total

    W_total = SoilWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp, z_col_bot) + &
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

    ! Add contributions from parts of the layer below -z_trunc (soil_frac ≈ 1)
    if (z_bot < -z_trunc) then
       ! From z_bot to -z_trunc: soil_frac ≈ 1 (Fs negligible)
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
       ! From z_trunc to z_top: no soil (all surface), surface water if z_wt > z
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
  ! Sy_soil = (1/dz) * integral from -z_col_bot to +z_trunc of
  !           (1-Fs(z)) * [theta(z_u-z) - theta(z_l-z)] dz
  !
  ! Note: currently integrates over [-z_trunc, +z_trunc] since Sy
  ! differences are negligible in the deep zone where soil_frac ≈ 1
  ! and both theta(z_u-z) and theta(z_l-z) ≈ theta_s.
  !========================================================================
  function SysoilMicroTopo(z_l, z_u, theta_s, h_e, b_camp) result(sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_l, z_u, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: sy
    real(kind=kind_noahmp) :: dz, theta_u, theta_l
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    dz = z_u - z_l
    if (abs(dz) < 1.0e-10_kind_noahmp) then
       sy = 0.0_kind_noahmp
       return
    endif

    ! Use precomputed wt_soil_micro (= gl_weights * soil_frac at GL nodes)
    sy = 0.0_kind_noahmp
    do k = 1, n_gl
       theta_u = theta_campbell(z_u - z_trunc * gl_nodes(k), theta_s, h_e, b_camp)
       theta_l = theta_campbell(z_l - z_trunc * gl_nodes(k), theta_s, h_e, b_camp)
       sy = sy + wt_soil_micro(k) * (theta_u - theta_l)
    enddo
    sy = sy * z_trunc / dz

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

    ! Below -z_trunc: Fs ≈ 0, negligible contribution
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
  ! Find water table from soil water storage (microtopo-aware)
  !
  ! Given total soil water content [m] (sum of SoilLiqWater × dz),
  ! bisects on SoilWaterStorageMicroTopoLite to find z_wt such that
  ! the horizontally-integrated hydrostatic soil water equals the
  ! observed amount.  Surface water is NOT included — it is a
  ! diagnostic computed from SurfaceWaterStorage(z_wt) after the
  ! water table is found.
  !
  ! Uses warm-start bracketing when z_wt_prev is provided: starts
  ! with a ±0.5 m bracket around the hint, expanding to the full
  ! range [-z_col_bot, z_trunc] if the target is not contained.
  !
  ! soil_water_content: total soil water [m], sum of SoilLiqWater(i)*dz(i)
  ! theta_s, h_e, b_camp: Campbell soil parameters
  ! z_col_bot: soil column bottom depth [m], positive downward
  ! z_wt_prev: optional warm-start hint (z_wt from previous timestep)
  ! Returns:   z_wt [m], positive upward (D&B convention)
  !========================================================================
  function FindWaterTable(soil_water_content, theta_s, h_e, b_camp, z_col_bot, z_wt_prev) result(z_wt)
    implicit none
    real(kind=kind_noahmp), intent(in) :: soil_water_content, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp), intent(in), optional :: z_wt_prev
    real(kind=kind_noahmp) :: z_wt
    real(kind=kind_noahmp) :: z_lo, z_hi, z_new, W_lo, W_hi
    real(kind=kind_noahmp) :: W_val, dW_val, h_pt, theta_val, inv_b, A_deep
    integer :: iter, k
    integer, parameter :: max_iter = 20
    real(kind=kind_noahmp), parameter :: tol = 1.0e-5_kind_noahmp
    real(kind=kind_noahmp), parameter :: warm_margin = 0.5_kind_noahmp

    if (.not. gl_initialized) call InitGaussLegendre()

    inv_b = -1.0_kind_noahmp / b_camp

    ! Deep-zone contribution (constant when fully saturated = common case)
    A_deep = 0.0_kind_noahmp
    if (z_col_bot > z_trunc) A_deep = theta_s * (z_col_bot - z_trunc)

    ! --- Set up bracket with optional warm-start ---
    if (present(z_wt_prev)) then
       z_lo = max(-z_col_bot, z_wt_prev - warm_margin)
       z_hi = min( z_trunc,   z_wt_prev + warm_margin)
       W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot)
       W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot)
       if (soil_water_content < W_lo .or. soil_water_content > W_hi) then
          z_lo = -z_col_bot
          z_hi =  z_trunc
          W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot)
          W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot)
       endif
    else
       z_lo = -z_col_bot
       z_hi =  z_trunc
       W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot)
       W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot)
    endif

    if (soil_water_content <= W_lo) then
       z_wt = z_lo;  return
    endif
    if (soil_water_content >= W_hi) then
       z_wt = z_hi;  return
    endif

    ! --- Safeguarded Newton-Raphson ---
    ! Start from warm-start hint (clamped to bracket) or midpoint
    if (present(z_wt_prev)) then
       z_wt = max(z_lo, min(z_hi, z_wt_prev))
    else
       z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
    endif

    do iter = 1, max_iter
       ! Compute W(z_wt) and dW/dz_wt in a single pass
       ! Deep zone: fully saturated → constant, derivative = 0
       W_val  = A_deep
       dW_val = 0.0_kind_noahmp
       ! Microtopo zone: use precomputed wt_soil_micro
       do k = 1, n_gl
          h_pt = z_wt - z_trunc * gl_nodes(k)
          if (h_pt >= -h_e) then
             W_val = W_val + wt_soil_micro(k) * theta_s * z_trunc
          else
             theta_val = theta_s * (abs(h_pt) / h_e) ** inv_b
             W_val  = W_val  + wt_soil_micro(k) * theta_val * z_trunc
             dW_val = dW_val + wt_soil_micro(k) * theta_val / &
                      (b_camp * abs(h_pt)) * z_trunc
          endif
       enddo

       if (abs(W_val - soil_water_content) < tol) return

       ! Update bracket
       if (W_val < soil_water_content) then
          z_lo = z_wt
       else
          z_hi = z_wt
       endif
       if ((z_hi - z_lo) < tol) then
          z_wt = 0.5_kind_noahmp * (z_lo + z_hi);  return
       endif

       ! Newton step with bracket safeguard
       if (dW_val > 1.0e-12_kind_noahmp) then
          z_new = z_wt + (soil_water_content - W_val) / dW_val
          if (z_new > z_lo .and. z_new < z_hi) then
             z_wt = z_new
          else
             z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
          endif
       else
          z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
       endif
    enddo

  end function FindWaterTable

  !========================================================================
  ! Compute the soil-only water deficit (relative to full saturation)
  ! for a given water table, including microtopography effect.
  ! Used for the WaterTableEquilibrium calculation.
  !
  ! deficit = integral of (theta_s - theta(z)) * (1-Fs(z)) dz
  !           from elevation -z_col_bot to +z_trunc
  !
  ! Computed as: theta_s * total_soil_volume - A_soil(z_wt)
  !
  ! z_wt:      water table position [m], positive upward (D&B convention)
  ! z_col_bot: soil column bottom depth [m], positive downward from mean
  !            surface.  Deep below the surface, Fs ≈ 0 so soil_frac ≈ 1
  !            naturally from the Gaussian CDF tail.
  !========================================================================
  function SoilDeficitMicroTopo(z_wt, theta_s, h_e, b_camp, z_col_bot) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: deficit
    real(kind=kind_noahmp) :: total_soil_vol, A_soil

    ! Total soil volume = EffSoilThickMicroTopo from hummock top to column bottom
    ! (EffSoilThickMicroTopo takes depth coordinates: positive downward)
    total_soil_vol = EffSoilThickMicroTopo(-z_trunc, z_col_bot)

    A_soil = SoilWaterStorageMicroTopo(z_wt, theta_s, h_e, b_camp, z_col_bot)

    deficit = theta_s * total_soil_vol - A_soil
    deficit = max(0.0_kind_noahmp, deficit)

  end function SoilDeficitMicroTopo

  !========================================================================
  ! Single-column hydrostatic equilibrium deficit [m of water]
  ! for a standard 1D soil column (no microtopography weighting).
  !
  ! This computes the analytical integral:
  !   deficit = integral from 0 to WTD of (theta_s - theta_campbell(h)) dh
  ! where h = WTD - z is the suction at depth z above the water table.
  !
  ! For the Campbell retention curve, the integral has a closed form:
  !   deficit = theta_s * [(WTD - h_e) - h_e*b/(b-1) * ((WTD/h_e)^((b-1)/b) - 1)]
  ! valid for WTD > h_e.  For WTD <= h_e, the entire column above WT
  ! is within the air-entry zone and fully saturated: deficit = 0.
  !
  ! This is used for WTD diagnosis from the NoahMP 4-layer soil moisture
  ! profile, which represents a single vertical column (not an areally-
  ! averaged microtopographic distribution).
  !========================================================================
  pure function SingleColumnDeficit(WTD, theta_s, h_e, b_camp) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD      ! water table depth [m], positive downward
    real(kind=kind_noahmp), intent(in) :: theta_s  ! saturated water content [m3/m3]
    real(kind=kind_noahmp), intent(in) :: h_e      ! air-entry suction head [m], positive
    real(kind=kind_noahmp), intent(in) :: b_camp   ! Campbell b exponent
    real(kind=kind_noahmp) :: deficit
    real(kind=kind_noahmp) :: bm1_over_b, b_over_bm1

    ! No deficit when water table is at or above the surface,
    ! or when WTD is within the air-entry zone
    if (WTD <= h_e) then
       deficit = 0.0_kind_noahmp
       return
    endif

    bm1_over_b = (b_camp - 1.0_kind_noahmp) / b_camp
    b_over_bm1 = b_camp / (b_camp - 1.0_kind_noahmp)

    ! Analytical integral of (theta_s - theta_campbell(h)) from h_e to WTD
    deficit = theta_s * ( (WTD - h_e) - h_e * b_over_bm1 * &
              ((WTD / h_e)**bm1_over_b - 1.0_kind_noahmp) )
    deficit = max(0.0_kind_noahmp, deficit)

  end function SingleColumnDeficit

  !========================================================================
  ! Effective soil thickness [m] for a depth range, accounting for
  ! microtopography.
  !
  ! z_top, z_bot: depth limits [m], positive downward from mean surface.
  !   z_top < z_bot.  z_top may be negative (above mean surface, hummocks).
  !
  ! At depth z (positive downward), the soil fraction is:
  !   z in [-z_trunc, z_trunc]: Phi(z / sigma_elev)
  !     = CDF of the surface elevation distribution
  !     (fraction of area that has soil at this depth)
  !   z > z_trunc:  1  (below all surface points, entire area is soil)
  !   z < -z_trunc: ≈ 0  (above virtually all hummock peaks)
  !
  ! Returns the integral of soil_fraction over [z_top, z_bot].
  !========================================================================
  function EffSoilThickMicroTopo(z_top, z_bot) result(eff_thick)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_top, z_bot
    real(kind=kind_noahmp) :: eff_thick
    real(kind=kind_noahmp) :: z_a, z_b, z_mid, z_half, z_pt
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    eff_thick = 0.0_kind_noahmp

    ! Part within the microtopo zone: clamp to [-z_trunc, z_trunc]
    z_a = max(z_top, -z_trunc)
    z_b = min(z_bot,  z_trunc)

    if (z_b > z_a) then
      z_mid  = 0.5_kind_noahmp * (z_b + z_a)
      z_half = 0.5_kind_noahmp * (z_b - z_a)
      do k = 1, n_gl
        z_pt = z_mid + z_half * gl_nodes(k)
        ! Soil fraction at depth z_pt = Phi(z_pt / sigma_elev)
        eff_thick = eff_thick + gl_weights(k) * phi_normal(z_pt / sigma_elev)
      enddo
      eff_thick = eff_thick * z_half
    endif

    ! Part below the microtopo zone (z > z_trunc): soil_frac = 1
    if (z_bot > z_trunc) then
      eff_thick = eff_thick + (z_bot - max(z_top, z_trunc))
    endif

  end function EffSoilThickMicroTopo

  !========================================================================
  ! Hydrostatic equilibrium soil moisture for a FLAT surface layer
  !
  ! Analytical closed-form average of the Campbell retention curve over
  ! a layer at hydrostatic equilibrium.  No numerical quadrature needed.
  !
  ! The soil column is split at d_sat = WTD - h_e:
  !   d >= d_sat: saturated (within air-entry zone or below WT), theta = theta_s
  !   d <  d_sat: unsaturated, theta = theta_s * ((WTD-d)/h_e)^(-1/b)
  ! The unsaturated integral has the closed form:
  !   theta_s * h_e^(1/b) * b/(b-1) * [u_top^((b-1)/b) - u_bot^((b-1)/b)]
  ! where u = WTD - d is the suction at depth d.
  !
  ! d_top, d_bot: layer depth bounds [m], positive downward (d_top < d_bot)
  ! WTD:          water table depth [m], positive downward
  !========================================================================
  pure function EquilibriumSMFlat(d_top, d_bot, WTD, theta_s, h_e, b_camp) result(theta_eq)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d_top, d_bot, WTD, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: theta_eq
    real(kind=kind_noahmp) :: layer_thick, d_sat
    real(kind=kind_noahmp) :: bm1_over_b, b_over_bm1
    real(kind=kind_noahmp) :: u_top, u_bot, integral_unsat

    layer_thick = d_bot - d_top

    ! d_sat = depth above which soil is unsaturated (suction > h_e)
    d_sat = WTD - h_e

    ! Case 1: entire layer at or below saturation boundary
    if (d_top >= d_sat) then
       theta_eq = theta_s
       return
    endif

    bm1_over_b = (b_camp - 1.0_kind_noahmp) / b_camp
    b_over_bm1 = b_camp / (b_camp - 1.0_kind_noahmp)

    if (d_bot <= d_sat) then
       ! Case 2: entire layer unsaturated
       u_top = WTD - d_top   ! suction at layer top (> h_e)
       u_bot = WTD - d_bot   ! suction at layer bottom (>= h_e)
       integral_unsat = theta_s * h_e**(1.0_kind_noahmp / b_camp) * b_over_bm1 * &
                        (u_top**bm1_over_b - u_bot**bm1_over_b)
       theta_eq = integral_unsat / layer_thick
    else
       ! Case 3: layer spans saturation boundary at d_sat
       u_top = WTD - d_top   ! suction at layer top
       u_bot = h_e            ! suction at d_sat boundary
       integral_unsat = theta_s * h_e**(1.0_kind_noahmp / b_camp) * b_over_bm1 * &
                        (u_top**bm1_over_b - u_bot**bm1_over_b)
       theta_eq = (integral_unsat + theta_s * (d_bot - d_sat)) / layer_thick
    endif

  end function EquilibriumSMFlat

  !========================================================================
  ! Hydrostatic equilibrium soil moisture integrated over microtopography
  ! (SOIL WATER ONLY — no surface water in hollows)
  !
  ! Computes the average volumetric soil moisture in a layer, weighted by
  ! the fraction of the microtopographic surface that has soil at each
  ! depth:
  !
  !   theta_eq = (1 / layer_thick) *
  !              integral_{d_top}^{d_bot} (1 - Fs(-d)) * theta(d - WTD) dd
  !
  ! where (1 - Fs(-d)) = fraction of area with soil at depth d, and
  ! theta(d - WTD) is the Campbell retention at hydrostatic pressure.
  !
  ! This gives total soil water in the layer divided by layer thickness.
  ! Surface water in hollows is NOT included (tracked via FSW_change).
  !
  ! The result is always <= theta_s because (1 - Fs(-d)) <= 1.
  ! For deep layers (d >> z_trunc), soil_frac -> 1 and the result
  ! converges to the flat-surface equilibrium.
  !
  ! d_top, d_bot: layer depth bounds [m], positive downward (d_top < d_bot)
  ! WTD:          water table depth [m], positive downward
  !========================================================================
  function EquilibriumSMMicroTopo(d_top, d_bot, WTD, theta_s, h_e, b_camp) result(theta_eq)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d_top, d_bot, WTD, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: theta_eq
    real(kind=kind_noahmp) :: d_mid, d_half, d_pt, h_pt, soil_frac
    integer :: k

    ! Fast path: for layers well below the microtopo zone (d_top > ~2.3*sigma),
    ! soil_frac ≈ 1 everywhere and the result equals EquilibriumSMFlat.
    ! Avoids 5-pt GL quadrature + Fs_cdf calls for layers 3 and 4.
    if (d_top > 2.326_kind_noahmp * sigma_elev) then
       theta_eq = EquilibriumSMFlat(d_top, d_bot, WTD, theta_s, h_e, b_camp)
       return
    endif

    if (.not. gl_initialized) call InitGaussLegendre()

    d_mid  = 0.5_kind_noahmp * (d_top + d_bot)
    d_half = 0.5_kind_noahmp * (d_bot - d_top)

    theta_eq = 0.0_kind_noahmp
    do k = 1, n_gl_lite
       d_pt = d_mid + d_half * gl_nodes_lite(k)
       ! Soil fraction at depth d_pt: fraction of surface ABOVE elevation -d_pt
       soil_frac = 1.0_kind_noahmp - Fs_cdf(-d_pt)
       ! Pressure head at depth d_pt for WT at WTD
       h_pt = d_pt - WTD
       theta_eq = theta_eq + gl_weights_lite(k) * soil_frac * theta_campbell(h_pt, theta_s, h_e, b_camp)
    enddo
    theta_eq = theta_eq * d_half / (d_bot - d_top)

  end function EquilibriumSMMicroTopo

  !========================================================================
  ! Find water table depth from flat-surface soil moisture deficit.
  !
  ! Given the total column deficit [m] computed from a flat 1D soil
  ! moisture profile:
  !   deficit = sum_layers( (theta_s - SM(i)) * dz(i) )
  !
  ! finds WTD such that SingleColumnDeficit(WTD) = deficit.
  !
  ! Uses bisection on the monotonically increasing SingleColumnDeficit.
  !
  ! deficit:   target deficit [m], must be >= 0
  ! z_col_bot: maximum column depth [m], positive downward
  ! Returns:   WTD [m], positive downward (0 = surface, z_col_bot = bottom)
  !========================================================================
  function FindWaterTableFlat(deficit, theta_s, h_e, b_camp, z_col_bot) result(WTD)
    implicit none
    real(kind=kind_noahmp), intent(in) :: deficit, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: WTD
    real(kind=kind_noahmp) :: zwt_lo, zwt_hi, zwt_mid, def_mid
    integer :: iter
    integer, parameter :: max_iter = 60
    real(kind=kind_noahmp), parameter :: tol = 1.0e-5_kind_noahmp

    ! If deficit is negligible, water table is at the surface
    if (deficit <= tol) then
       WTD = 0.0_kind_noahmp
       return
    endif

    ! Bisection range: WTD from 0 (surface) to z_col_bot (column bottom)
    zwt_lo = 0.0_kind_noahmp
    zwt_hi = z_col_bot

    ! Check if deficit exceeds column capacity
    if (SingleColumnDeficit(zwt_hi, theta_s, h_e, b_camp) < deficit) then
       WTD = zwt_hi
       return
    endif

    ! Bisection
    do iter = 1, max_iter
       zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
       def_mid = SingleColumnDeficit(zwt_mid, theta_s, h_e, b_camp)

       if (abs(def_mid - deficit) < tol .or. (zwt_hi - zwt_lo) < tol) then
          WTD = zwt_mid
          return
       endif

       if (def_mid < deficit) then
          zwt_lo = zwt_mid   ! deficit too small → go deeper
       else
          zwt_hi = zwt_mid   ! deficit too large → go shallower
       endif
    enddo

    WTD = 0.5_kind_noahmp * (zwt_lo + zwt_hi)

  end function FindWaterTableFlat

end module PeatMicroTopoMod
