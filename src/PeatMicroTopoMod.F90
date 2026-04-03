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

  ! Precomputed quadrature coefficients for column-averaging over the
  ! truncated Gaussian surface elevation distribution.
  ! GL nodes map [-1,1] → elevation [-z_trunc, +z_trunc].
  ! wt_pdf(k)    = gl_weights(k) * z_trunc * phi_pdf(z_s_nodes(k)/sigma) / (sigma * norm_factor)
  ! z_s_nodes(k) = z_trunc * gl_nodes(k)
  ! All microtopo-aware functions become: sum_k wt_pdf(k) * FlatColumnFunction(z_s_nodes(k) - z_wt)
  real(kind=kind_noahmp), dimension(n_gl) :: wt_pdf
  real(kind=kind_noahmp), dimension(n_gl) :: z_s_nodes

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

    ! Precompute column-averaging quadrature weights.
    ! GL nodes map [-1,1] → surface elevation [-z_trunc, +z_trunc].
    ! wt_pdf(k) encodes the truncated Gaussian PDF weight at each node.
    block
      real(kind=kind_noahmp) :: norm_factor
      norm_factor = phi_normal(z_trunc / sigma_elev) - phi_normal(-z_trunc / sigma_elev)
      do k = 1, n_gl
         z_s_nodes(k) = z_trunc * gl_nodes(k)
         wt_pdf(k)    = gl_weights(k) * z_trunc * &
                         phi_pdf(z_s_nodes(k) / sigma_elev) / (sigma_elev * norm_factor)
      enddo
    end block

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
  ! microtopography within the NoahMP layer domain (elevation -z_col_bot
  ! to 0, i.e. depth 0 to z_col_bot).
  !
  ! A_soil(z_wt) = integral from -z_col_bot to 0 of
  !                (1 - F_s(z)) * theta(z_wt - z) dz
  !
  ! where theta is the Campbell retention with hydrostatic assumption
  ! h = z_wt - z (pressure head at elevation z for water table at z_wt)
  !
  ! The integration stops at z = 0 (mean surface).  Hummock soil above
  ! the mean surface (z > 0) is NOT tracked by any NoahMP soil layer
  ! and is therefore excluded.
  !
  ! Composite Gauss-Legendre quadrature: the domain is split at
  ! z = -z_trunc into two sub-intervals to maintain resolution.
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

    ! --- Sub-interval 2: microtopo zone [-z_trunc, 0] ---
    ! Full Gaussian CDF variation of soil_frac occurs here.
    ! Uses 400-point GL rule for smooth deficit-vs-WTD relationship.
    ! Upper limit = 0 (mean surface): soil above hummock peaks (z > 0)
    ! is not tracked by NoahMP layers.
    z_lo = -z_trunc
    z_hi =  0.0_kind_noahmp
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
  ! Column-averaged soil water storage [m] for water table at z_wt.
  !
  ! W(z_wt) = sum_k wt_pdf(k) * W_flat(z_s_nodes(k) - z_wt, z_col_bot)
  !
  ! where W_flat(WTD) = theta_s * z_cb - SingleColumnDeficit(WTD)
  ! for WTD > 0, or W_flat = theta_s * z_cb for WTD <= 0.
  !
  ! This function is monotonically increasing over
  ! [-z_col_bot - z_trunc, +z_trunc], allowing unique inversion.
  !
  ! Designed for fast repeated evaluation inside FindWaterTable.
  !========================================================================
  function SoilWaterStorageMicroTopoLite(z_wt, theta_s, h_e, b_camp, z_col_bot) result(W_col_avg)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: W_col_avg
    real(kind=kind_noahmp) :: W_sat, WTD_local
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    W_sat = theta_s * z_col_bot  ! storage of a fully saturated column

    W_col_avg = 0.0_kind_noahmp
    do k = 1, n_gl
       WTD_local = z_s_nodes(k) - z_wt   ! local WTD for this column (positive = water table below surface)
       if (WTD_local <= 0.0_kind_noahmp) then
          ! Column is fully saturated (or flooded)
          W_col_avg = W_col_avg + wt_pdf(k) * W_sat
       else
          W_col_avg = W_col_avg + wt_pdf(k) * (W_sat - SingleColumnDeficit(WTD_local, theta_s, h_e, b_camp))
       endif
    enddo

  end function SoilWaterStorageMicroTopoLite

  !========================================================================
  ! Specific yield of soil component (Eq. 5 / Eq. 6 in D&B 2015)
  ! Sy_soil for water level change from z_l to z_u, column-averaged.
  !
  ! Sy_soil = (1/dz) * sum_k wt_pdf(k) * [W_flat(z_s(k)-z_u) - W_flat(z_s(k)-z_l)]
  !
  ! Only columns with WTD > 0 contribute non-trivially.
  !========================================================================
  function SysoilMicroTopo(z_l, z_u, theta_s, h_e, b_camp) result(sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_l, z_u, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: sy
    real(kind=kind_noahmp) :: dz, WTD_u, WTD_l, D_u, D_l
    integer :: k

    if (.not. gl_initialized) call InitGaussLegendre()

    dz = z_u - z_l
    if (abs(dz) < 1.0e-10_kind_noahmp) then
       sy = 0.0_kind_noahmp
       return
    endif

    ! Column-averaging: integrate over surface elevation distribution
    sy = 0.0_kind_noahmp
    do k = 1, n_gl
       WTD_u = z_s_nodes(k) - z_u   ! local WTD at upper water table position
       WTD_l = z_s_nodes(k) - z_l   ! local WTD at lower water table position
       ! W_flat = theta_s * z_cb - deficit; difference = deficit_l - deficit_u
       D_u = SingleColumnDeficit(max(0.0_kind_noahmp, WTD_u), theta_s, h_e, b_camp)
       D_l = SingleColumnDeficit(max(0.0_kind_noahmp, WTD_l), theta_s, h_e, b_camp)
       sy = sy + wt_pdf(k) * (D_l - D_u)
    enddo
    sy = sy / dz

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
  ! Find water table from soil water storage (column-averaged)
  !
  ! Given total soil water content [m] (sum of SoilLiqWater × dz),
  ! finds z_wt such that the column-averaged hydrostatic soil water
  ! storage equals the observed amount.
  !
  ! Uses the column-averaged SoilWaterStorageMicroTopoLite which
  ! integrates W_flat(z_s - z_wt) over the Gaussian surface distribution.
  !
  ! Bracket: [-z_col_bot - z_trunc, +z_trunc]
  !   - At upper bound: all columns saturated, W = theta_s * z_cb
  !   - At lower bound: driest possible state
  !
  ! Derivative: dW/dz_wt = sum_k wt_pdf(k) * Sy_flat(z_s(k) - z_wt)
  ! where Sy_flat(WTD) = theta_s * [1 - (WTD/h_e)^(-1/b)] for WTD > h_e
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
    real(kind=kind_noahmp) :: W_val, dW_val, WTD_local, Sy_flat_val, inv_b
    integer :: iter, k
    integer, parameter :: max_iter = 30
    real(kind=kind_noahmp), parameter :: tol = 1.0e-5_kind_noahmp
    real(kind=kind_noahmp), parameter :: warm_margin = 0.5_kind_noahmp

    if (.not. gl_initialized) call InitGaussLegendre()

    inv_b = -1.0_kind_noahmp / b_camp

    ! --- Full bracket: [-z_col_bot - z_trunc, +z_trunc] ---
    ! z_wt = +z_trunc: all columns saturated, W = theta_s * z_cb
    ! z_wt = -(z_col_bot + z_trunc): driest possible state

    ! --- Set up bracket with optional warm-start ---
    if (present(z_wt_prev)) then
       z_lo = max(-(z_col_bot + z_trunc), z_wt_prev - warm_margin)
       z_hi = min( z_trunc,               z_wt_prev + warm_margin)
       W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot)
       W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot)
       if (soil_water_content < W_lo .or. soil_water_content > W_hi) then
          z_lo = -(z_col_bot + z_trunc)
          z_hi =  z_trunc
          W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot)
          W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot)
       endif
    else
       z_lo = -(z_col_bot + z_trunc)
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
    if (present(z_wt_prev)) then
       z_wt = max(z_lo, min(z_hi, z_wt_prev))
    else
       z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
    endif

    do iter = 1, max_iter
       ! Compute W(z_wt) and dW/dz_wt in a single pass using column-averaging
       W_val  = 0.0_kind_noahmp
       dW_val = 0.0_kind_noahmp
       do k = 1, n_gl
          WTD_local = z_s_nodes(k) - z_wt
          if (WTD_local <= 0.0_kind_noahmp) then
             ! Saturated column: W = theta_s * z_cb, dW/dz_wt = 0
             W_val = W_val + wt_pdf(k) * theta_s * z_col_bot
          else
             W_val = W_val + wt_pdf(k) * (theta_s * z_col_bot - &
                     SingleColumnDeficit(WTD_local, theta_s, h_e, b_camp))
             ! Sy_flat = theta_s * [1 - (WTD/h_e)^(-1/b)] for WTD > h_e, else 0
             if (WTD_local > h_e) then
                Sy_flat_val = theta_s * (1.0_kind_noahmp - (WTD_local / h_e) ** inv_b)
             else
                Sy_flat_val = 0.0_kind_noahmp
             endif
             dW_val = dW_val + wt_pdf(k) * Sy_flat_val
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
  ! Find z_wt from TOTAL water storage (soil + surface ponding).
  !
  ! W_total(z_wt) = SoilWaterStorageMicroTopoLite(z_wt) + SurfaceWaterStorage(z_wt)
  !
  ! This is monotonically increasing and unbounded above (ponded water
  ! grows linearly above z_trunc).  Safeguarded Newton-Raphson.
  ! Derivative: dW_total/dz_wt = Sy_soil + FloodedFrac(z_wt)
  !========================================================================
  function FindWaterTableTotal(total_water, theta_s, h_e, b_camp, z_col_bot, z_wt_prev) result(z_wt)
    implicit none
    real(kind=kind_noahmp), intent(in) :: total_water  ! W_soil + W_surface [m]
    real(kind=kind_noahmp), intent(in) :: theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp), intent(in), optional :: z_wt_prev
    real(kind=kind_noahmp) :: z_wt
    real(kind=kind_noahmp) :: z_lo, z_hi, z_new, W_lo, W_hi
    real(kind=kind_noahmp) :: W_val, dW_val, WTD_local, Sy_flat_val, inv_b
    real(kind=kind_noahmp) :: flood_frac, sqrt2_inv, norm_den, Phi_neg
    integer :: iter, k
    integer, parameter :: max_iter = 40
    real(kind=kind_noahmp), parameter :: tol = 1.0e-6_kind_noahmp
    real(kind=kind_noahmp), parameter :: warm_margin = 0.5_kind_noahmp

    if (.not. gl_initialized) call InitGaussLegendre()

    inv_b = -1.0_kind_noahmp / b_camp
    sqrt2_inv = 1.0_kind_noahmp / sqrt(2.0_kind_noahmp)

    ! Precompute truncated-Gaussian CDF helpers
    norm_den  = erf(z_trunc * sqrt2_inv / sigma_elev)
    Phi_neg   = 0.5_kind_noahmp * (1.0_kind_noahmp - norm_den)

    ! --- Bracket ---
    z_lo = -(z_col_bot + z_trunc)
    ! Above z_trunc FSW grows linearly with slope 1
    z_hi = z_trunc + max(0.0_kind_noahmp, total_water - theta_s * z_col_bot)

    ! Warm start
    if (present(z_wt_prev)) then
       z_lo = max(z_lo, z_wt_prev - warm_margin)
       z_hi = min(z_hi, z_wt_prev + warm_margin)
       W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot) + &
              SurfaceWaterStorage(z_lo)
       W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot) + &
              SurfaceWaterStorage(z_hi)
       if (total_water < W_lo .or. total_water > W_hi) then
          z_lo = -(z_col_bot + z_trunc)
          z_hi = z_trunc + max(0.0_kind_noahmp, total_water - theta_s * z_col_bot)
       endif
    endif

    W_lo = SoilWaterStorageMicroTopoLite(z_lo, theta_s, h_e, b_camp, z_col_bot) + &
           SurfaceWaterStorage(z_lo)
    W_hi = SoilWaterStorageMicroTopoLite(z_hi, theta_s, h_e, b_camp, z_col_bot) + &
           SurfaceWaterStorage(z_hi)

    if (total_water <= W_lo) then
       z_wt = z_lo;  return
    endif
    if (total_water >= W_hi) then
       z_wt = z_hi;  return
    endif

    ! Starting point
    if (present(z_wt_prev)) then
       z_wt = max(z_lo, min(z_hi, z_wt_prev))
    else
       z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
    endif

    do iter = 1, max_iter
       ! --- W_total(z_wt) and dW_total/dz_wt in one pass ---
       ! Soil part
       W_val  = 0.0_kind_noahmp
       dW_val = 0.0_kind_noahmp
       do k = 1, n_gl
          WTD_local = z_s_nodes(k) - z_wt
          if (WTD_local <= 0.0_kind_noahmp) then
             W_val = W_val + wt_pdf(k) * theta_s * z_col_bot
          else
             W_val = W_val + wt_pdf(k) * (theta_s * z_col_bot - &
                     SingleColumnDeficit(WTD_local, theta_s, h_e, b_camp))
             if (WTD_local > h_e) then
                Sy_flat_val = theta_s * (1.0_kind_noahmp - (WTD_local / h_e) ** inv_b)
             else
                Sy_flat_val = 0.0_kind_noahmp
             endif
             dW_val = dW_val + wt_pdf(k) * Sy_flat_val
          endif
       enddo

       ! Surface water part
       W_val = W_val + SurfaceWaterStorage(z_wt)

       ! dFSW/dz_wt = FloodedFrac(z_wt), computed analytically
       if (z_wt <= -z_trunc) then
          flood_frac = 0.0_kind_noahmp
       elseif (z_wt >= z_trunc) then
          flood_frac = 1.0_kind_noahmp
       else
          flood_frac = (0.5_kind_noahmp * (1.0_kind_noahmp + &
                        erf(z_wt * sqrt2_inv / sigma_elev)) - Phi_neg) / norm_den
       endif
       dW_val = dW_val + flood_frac

       if (abs(W_val - total_water) < tol) return

       ! Update bracket
       if (W_val < total_water) then
          z_lo = z_wt
       else
          z_hi = z_wt
       endif
       if ((z_hi - z_lo) < tol) then
          z_wt = 0.5_kind_noahmp * (z_lo + z_hi);  return
       endif

       ! Newton step with safeguard
       if (dW_val > 1.0e-12_kind_noahmp) then
          z_new = z_wt + (total_water - W_val) / dW_val
          if (z_new > z_lo .and. z_new < z_hi) then
             z_wt = z_new
          else
             z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
          endif
       else
          z_wt = 0.5_kind_noahmp * (z_lo + z_hi)
       endif
    enddo

  end function FindWaterTableTotal

  !========================================================================
  ! Compute the soil-only water deficit (relative to full saturation)
  ! for a given water table, including microtopography effect.
  ! Used for the WaterTableEquilibrium calculation.
  !
  ! deficit = integral of (theta_s - theta(z)) * (1-Fs(z)) dz
  !           from elevation -z_col_bot to 0  (NoahMP soil domain)
  !
  ! Computed as: theta_s * total_soil_volume - A_soil(z_wt)
  !
  ! z_wt:      water table position [m], positive upward (D&B convention)
  ! z_col_bot: soil column bottom depth [m], positive downward from mean
  !            surface.
  !========================================================================
  function SoilDeficitMicroTopo(z_wt, theta_s, h_e, b_camp, z_col_bot) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_wt, theta_s, h_e, b_camp, z_col_bot
    real(kind=kind_noahmp) :: deficit
    real(kind=kind_noahmp) :: total_soil_vol, A_soil

    ! Total soil volume = EffSoilThickMicroTopo from surface to column bottom
    ! (depth 0 to z_col_bot; hummock soil above surface excluded)
    total_soil_vol = EffSoilThickMicroTopo(0.0_kind_noahmp, z_col_bot)

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
  ! Column-averaged hydrostatic equilibrium soil moisture
  !
  ! theta_eq(d_top, d_bot, z_wt) = sum_k wt_pdf(k) *
  !     EquilibriumSMFlat(d_top, d_bot, z_s_nodes(k) - z_wt, theta_s, h_e, b)
  !
  ! Each column has its surface at z_s, so its local WTD = z_s - z_wt.
  ! Columns with z_s <= z_wt are saturated (WTD <= 0) and contribute theta_s.
  ! Columns with z_s > z_wt have genuinely unsaturated soil in upper layers.
  !
  ! Fast path for deep layers (d_top > 2.3*sigma): all columns have 
  ! essentially the same unsaturated zone structure at this depth, so
  ! the result equals EquilibriumSMFlat directly.
  !
  ! d_top, d_bot: layer depth bounds [m], positive downward (d_top < d_bot)
  ! WTD:          water table depth [m], positive downward (= -z_wt)
  !========================================================================
  function EquilibriumSMMicroTopo(d_top, d_bot, WTD, theta_s, h_e, b_camp) result(theta_eq)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d_top, d_bot, WTD, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: theta_eq
    real(kind=kind_noahmp) :: z_wt_local, WTD_local
    integer :: k

    ! Fast path: layer is fully saturated in ALL columns (including the
    ! highest at z_trunc).  Highest column's local WTD = z_trunc + WTD.
    ! Saturation boundary: d_sat = WTD_local - h_e = z_trunc + WTD - h_e.
    ! Layer top must be below that for all columns to be saturated.
    if (d_top > z_trunc + WTD - h_e) then
       theta_eq = theta_s
       return
    endif

    if (.not. gl_initialized) call InitGaussLegendre()

    z_wt_local = -WTD  ! convert to D&B convention

    theta_eq = 0.0_kind_noahmp
    do k = 1, n_gl
       WTD_local = z_s_nodes(k) - z_wt_local   ! local WTD for this column
       ! EquilibriumSMFlat handles WTD <= 0 (returns theta_s)
       theta_eq = theta_eq + wt_pdf(k) * &
                  EquilibriumSMFlat(d_top, d_bot, WTD_local, theta_s, h_e, b_camp)
    enddo

  end function EquilibriumSMMicroTopo

  !========================================================================
  ! Column-averaged soil moisture for the microtopography-aware column
  ! after applying a uniform pressure-head anomaly to the hydrostatic
  ! equilibrium profile.
  !========================================================================
  function ThetaFromHeadShiftMicro(d_top, d_bot, WTD_ref, head_shift, theta_s, h_e, b_camp) result(theta_shift)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d_top, d_bot, WTD_ref, head_shift, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: theta_shift

    theta_shift = EquilibriumSMMicroTopo(d_top, d_bot, WTD_ref - head_shift, theta_s, h_e, b_camp)

  end function ThetaFromHeadShiftMicro

  !========================================================================
  ! Diagnose the uniform pressure-head anomaly that reproduces a target
  ! microtopography-averaged layer-mean soil moisture.
  !========================================================================
  function HeadShiftFromThetaMicro(theta_target, d_top, d_bot, WTD_ref, theta_s, h_e, b_camp) result(head_shift)
    implicit none
    real(kind=kind_noahmp), intent(in) :: theta_target, d_top, d_bot, WTD_ref, theta_s, h_e, b_camp
    real(kind=kind_noahmp) :: head_shift
    real(kind=kind_noahmp) :: shift_lo, shift_hi, shift_mid
    real(kind=kind_noahmp) :: theta_lo, theta_hi, theta_mid, theta_tgt
    integer :: iter
    integer, parameter :: max_iter = 50
    real(kind=kind_noahmp), parameter :: tol = 1.0e-8_kind_noahmp

    shift_lo = -(d_bot + z_trunc + max(WTD_ref, 0.0_kind_noahmp) + 10.0_kind_noahmp * h_e)
    shift_hi =   d_bot + z_trunc + max(WTD_ref, 0.0_kind_noahmp) + 10.0_kind_noahmp * h_e

    theta_lo = ThetaFromHeadShiftMicro(d_top, d_bot, WTD_ref, shift_lo, theta_s, h_e, b_camp)
    theta_hi = ThetaFromHeadShiftMicro(d_top, d_bot, WTD_ref, shift_hi, theta_s, h_e, b_camp)
    theta_tgt = max(theta_lo, min(theta_hi, theta_target))

    if (theta_tgt <= theta_lo + tol) then
       head_shift = shift_lo
       return
    endif
    if (theta_tgt >= theta_hi - tol) then
       head_shift = shift_hi
       return
    endif

    do iter = 1, max_iter
       shift_mid = 0.5_kind_noahmp * (shift_lo + shift_hi)
       theta_mid = ThetaFromHeadShiftMicro(d_top, d_bot, WTD_ref, shift_mid, theta_s, h_e, b_camp)
       if (abs(theta_mid - theta_tgt) < tol) exit
       if (theta_mid < theta_tgt) then
          shift_lo = shift_mid
       else
          shift_hi = shift_mid
       endif
    enddo

    head_shift = 0.5_kind_noahmp * (shift_lo + shift_hi)

  end function HeadShiftFromThetaMicro

end module PeatMicroTopoMod
