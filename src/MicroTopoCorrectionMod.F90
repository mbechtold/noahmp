module MicroTopoCorrectionMod

!!! Gaussian microtopography for peatland soil water retention
!!! Based on Dettmann & Bechtold (2015, Hydrological Processes)
!!! Replaces PEATCLSM empirical formulation (Chakraborty et al., 2026)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

  ! Microtopography parameters
  real(kind=kind_noahmp), parameter :: sigma_z    = 0.16_kind_noahmp    ! std dev of surface elevations [m]
  real(kind=kind_noahmp), parameter :: z_elev_max = 1.0_kind_noahmp     ! truncation height above mean [m]
  real(kind=kind_noahmp), parameter :: pi_val     = 3.14159265358979_kind_noahmp
  real(kind=kind_noahmp), parameter :: sqrt2      = 1.41421356237310_kind_noahmp
  real(kind=kind_noahmp), parameter :: sqrt2pi    = 2.50662827463100_kind_noahmp

  ! Numerical integration
  integer, parameter :: NumQuadPoints = 50   ! quadrature points for ensemble integration
  real(kind=kind_noahmp), parameter :: delta_Sy = 0.001_kind_noahmp  ! finite difference for Sy [m]

contains

  !=======================================================================
  ! GaussianCDF: Phi(x) = P(Z <= x) for standard N(0,1) scaled to N(0,sigma^2)
  !=======================================================================
  pure function GaussianCDF(x, sigma) result(phi)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x, sigma
    real(kind=kind_noahmp) :: phi

    phi = 0.5_kind_noahmp * erfc(-x / (sigma * sqrt2))
  end function GaussianCDF

  !=======================================================================
  ! GaussianPDF: phi(x) for N(0, sigma^2)
  !=======================================================================
  pure function GaussianPDF(x, sigma) result(pdf)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x, sigma
    real(kind=kind_noahmp) :: pdf

    pdf = exp(-0.5_kind_noahmp * (x/sigma)**2) / (sigma * sqrt2pi)
  end function GaussianPDF

  !=======================================================================
  ! FloodedFraction: fraction of surface that is flooded at given WTD
  ! WTD positive down => water table at elevation -WTD relative to mean
  ! A column at z_s is flooded if z_s < -WTD (WT above its surface)
  !=======================================================================
  pure function FloodedFraction(WTD, sigma) result(ff)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD, sigma
    real(kind=kind_noahmp) :: ff

    ff = GaussianCDF(-WTD, sigma)
    ff = max(0.0_kind_noahmp, min(1.0_kind_noahmp, ff))
  end function FloodedFraction

  !=======================================================================
  ! SurfaceWaterStorage_mm: total surface water [mm] from Gaussian
  ! microtopography using partial expectation formula.
  ! Surface water = integral over all flooded columns of (water depth) × pdf
  ! = 1000 * [sigma * phi_std(-WTD/sigma) + WTD * Phi(-WTD/sigma)]
  ! where phi_std is standard Gaussian PDF, Phi is standard Gaussian CDF
  ! WTD positive down: negative WTD means ponding above mean surface
  !=======================================================================
  pure function SurfaceWaterStorage_mm(WTD, sigma) result(storage)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD, sigma
    real(kind=kind_noahmp) :: storage
    real(kind=kind_noahmp) :: u  ! standardized variable

    u = -WTD / sigma
    ! Partial expectation: E[max(-WTD - z_s, 0)] = sigma*phi(u) + (-WTD)*Phi(u)
    ! where phi is std normal pdf, Phi is std normal CDF
    storage = 1000.0_kind_noahmp * (sigma * GaussianPDF(0.0_kind_noahmp, 1.0_kind_noahmp) * &
              exp(-0.5_kind_noahmp * u * u) / GaussianPDF(0.0_kind_noahmp, 1.0_kind_noahmp) + &
              (-WTD) * GaussianCDF(-WTD, sigma))

    ! Simplify: sigma * phi_std(u) = sigma * exp(-u^2/2)/sqrt(2*pi)
    !           (-WTD) * Phi_std(u)  = (-WTD) * 0.5*erfc(-u/sqrt(2))
    storage = 1000.0_kind_noahmp * ( &
              sigma * exp(-0.5_kind_noahmp * u * u) / sqrt2pi + &
              (-WTD) * 0.5_kind_noahmp * erfc(-u / sqrt2) )

    storage = max(0.0_kind_noahmp, storage)
  end function SurfaceWaterStorage_mm

  !=======================================================================
  ! CampbellTheta: equilibrium soil moisture from Campbell retention
  ! h = height above water table [m] (positive upward)
  ! Returns theta in [0, theta_s]
  !=======================================================================
  pure function CampbellTheta(h, theta_s, psi_ae, bcoeff) result(theta)
    implicit none
    real(kind=kind_noahmp), intent(in) :: h, theta_s, psi_ae, bcoeff
    real(kind=kind_noahmp) :: theta

    if (h <= psi_ae) then
       ! Below or at air-entry: saturated
       theta = theta_s
    else
       ! Campbell: theta = theta_s * (h/psi_ae)^(-1/b)
       theta = theta_s * (h / psi_ae) ** (-1.0_kind_noahmp / bcoeff)
       theta = max(0.0_kind_noahmp, min(theta_s, theta))
    endif
  end function CampbellTheta

  !=======================================================================
  ! SingleColumnDeficit: deficit for one column at surface elevation z_s
  ! with NumSoilLayer layers. WTD is positive down from mean surface.
  ! Deficit = integral of (theta_s - theta_equil) over the soil column,
  ! but only where soil exists (z < z_s, i.e., below column surface).
  !=======================================================================
  pure function SingleColumnDeficit(z_s, WTD, NumSoilLayer, DepthSoilLayer, &
                                     theta_s, psi_ae, bcoeff) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: z_s          ! surface elevation of this column [m] above mean
    real(kind=kind_noahmp), intent(in) :: WTD          ! water table depth from mean surface [m], positive down
    integer,                intent(in) :: NumSoilLayer
    real(kind=kind_noahmp), intent(in) :: DepthSoilLayer(NumSoilLayer)  ! cumulative depths, negative downward [m]
    real(kind=kind_noahmp), intent(in) :: theta_s, psi_ae, bcoeff
    real(kind=kind_noahmp)             :: deficit

    integer :: k, j
    integer, parameter :: Nsub = 10   ! sub-intervals per layer for integration
    real(kind=kind_noahmp) :: z_top, z_bot, z_top_eff, dz_sub, z_mid
    real(kind=kind_noahmp) :: h_above_wt, theta_eq, layer_deficit

    deficit = 0.0_kind_noahmp

    do k = 1, NumSoilLayer
       ! Layer boundaries relative to mean surface (z positive up):
       ! z_top = top of layer, z_bot = bottom of layer
       if (k == 1) then
          z_top = 0.0_kind_noahmp                     ! top of layer 1 = mean surface
       else
          z_top = DepthSoilLayer(k-1)                  ! negative
       endif
       z_bot = DepthSoilLayer(k)                       ! negative

       ! Effective top: soil only exists below z_s
       ! If z_s < z_top, this column doesn't reach above z_bot for this layer
       if (z_s <= z_bot) cycle                        ! column surface is below this layer entirely
       z_top_eff = min(z_top, z_s)                    ! clip layer top to column surface

       ! Integrate deficit over [z_bot, z_top_eff] using sub-intervals
       dz_sub = (z_top_eff - z_bot) / real(Nsub, kind_noahmp)
       if (dz_sub <= 0.0_kind_noahmp) cycle

       layer_deficit = 0.0_kind_noahmp
       do j = 1, Nsub
          z_mid = z_bot + (real(j, kind_noahmp) - 0.5_kind_noahmp) * dz_sub
          ! Height above water table: WT is at elevation -WTD
          h_above_wt = z_mid - (-WTD)         ! = z_mid + WTD
          if (h_above_wt < 0.0_kind_noahmp) then
             ! Below water table: saturated, no deficit
             theta_eq = theta_s
          else
             theta_eq = CampbellTheta(h_above_wt, theta_s, psi_ae, bcoeff)
          endif
          layer_deficit = layer_deficit + (theta_s - theta_eq) * dz_sub
       enddo

       deficit = deficit + layer_deficit
    enddo
  end function SingleColumnDeficit

  !=======================================================================
  ! MicroTopoStorageDeficit: ensemble-averaged deficit [m] across all
  ! columns weighted by Gaussian pdf of surface elevations.
  ! Numerically integrates over z_s from -4*sigma to z_elev_max.
  !=======================================================================
  pure function MicroTopoStorageDeficit(WTD, NumSoilLayer, DepthSoilLayer, &
                                         theta_s, psi_ae, bcoeff, sigma) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD
    integer,                intent(in) :: NumSoilLayer
    real(kind=kind_noahmp), intent(in) :: DepthSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp), intent(in) :: theta_s, psi_ae, bcoeff, sigma
    real(kind=kind_noahmp)             :: deficit

    integer :: i
    real(kind=kind_noahmp) :: z_lo, z_hi, dz, z_s, weight, col_def
    real(kind=kind_noahmp) :: total_weight

    deficit = 0.0_kind_noahmp
    total_weight = 0.0_kind_noahmp

    z_lo = -4.0_kind_noahmp * sigma    ! lower bound of integration
    z_hi = z_elev_max                   ! upper bound (truncation)
    dz = (z_hi - z_lo) / real(NumQuadPoints, kind_noahmp)

    do i = 1, NumQuadPoints
       z_s = z_lo + (real(i, kind_noahmp) - 0.5_kind_noahmp) * dz
       weight = GaussianPDF(z_s, sigma) * dz

       col_def = SingleColumnDeficit(z_s, WTD, NumSoilLayer, DepthSoilLayer, &
                                      theta_s, psi_ae, bcoeff)
       deficit = deficit + col_def * weight
       total_weight = total_weight + weight
    enddo

    ! Normalize by total weight to account for truncation
    if (total_weight > 0.0_kind_noahmp) then
       deficit = deficit / total_weight
    endif

    deficit = max(0.0_kind_noahmp, deficit)
  end function MicroTopoStorageDeficit

  !=======================================================================
  ! MicroTopoStorageDeficitExtended: deficit including below-model extension
  ! When WTD > column_depth, extends with flat-soil Campbell integral below.
  !=======================================================================
  pure function MicroTopoStorageDeficitExtended(WTD, NumSoilLayer, DepthSoilLayer, &
                                                 theta_s, psi_ae, bcoeff, sigma) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD
    integer,                intent(in) :: NumSoilLayer
    real(kind=kind_noahmp), intent(in) :: DepthSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp), intent(in) :: theta_s, psi_ae, bcoeff, sigma
    real(kind=kind_noahmp)             :: deficit

    real(kind=kind_noahmp) :: column_depth, extra_depth, h_top, h_bot
    real(kind=kind_noahmp) :: dz_ext, z_mid, h_above_wt, theta_eq
    integer :: j
    integer, parameter :: Nsub_ext = 50   ! sub-intervals for below-model integration

    ! Deficit within the model domain (multi-column)
    deficit = MicroTopoStorageDeficit(WTD, NumSoilLayer, DepthSoilLayer, &
                                      theta_s, psi_ae, bcoeff, sigma)

    ! Below-model extension: compute extra deficit below the domain bottom
    column_depth = -DepthSoilLayer(NumSoilLayer)   ! positive [m]
    if (WTD > column_depth) then
       ! Below model bottom, assume flat soil (f_soil = 1.0)
       ! Integrate from column_depth down to WTD
       extra_depth = WTD - column_depth
       dz_ext = extra_depth / real(Nsub_ext, kind_noahmp)
       do j = 1, Nsub_ext
          z_mid = -column_depth - (real(j, kind_noahmp) - 0.5_kind_noahmp) * dz_ext
          h_above_wt = z_mid + WTD   ! height above WT
          if (h_above_wt > 0.0_kind_noahmp) then
             theta_eq = CampbellTheta(h_above_wt, theta_s, psi_ae, bcoeff)
             deficit = deficit + (theta_s - theta_eq) * dz_ext
          endif
       enddo
    endif
  end function MicroTopoStorageDeficitExtended

  !=======================================================================
  ! MicroTopoSpecificYield: numerical derivative of deficit function
  !=======================================================================
  pure function MicroTopoSpecificYield(WTD, NumSoilLayer, DepthSoilLayer, &
                                        theta_s, psi_ae, bcoeff, sigma) result(Sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD
    integer,                intent(in) :: NumSoilLayer
    real(kind=kind_noahmp), intent(in) :: DepthSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp), intent(in) :: theta_s, psi_ae, bcoeff, sigma
    real(kind=kind_noahmp)             :: Sy

    real(kind=kind_noahmp) :: D_plus, D_minus

    D_plus  = MicroTopoStorageDeficitExtended(WTD + delta_Sy, NumSoilLayer, DepthSoilLayer, &
                                               theta_s, psi_ae, bcoeff, sigma)
    D_minus = MicroTopoStorageDeficitExtended(WTD - delta_Sy, NumSoilLayer, DepthSoilLayer, &
                                               theta_s, psi_ae, bcoeff, sigma)

    Sy = (D_plus - D_minus) / (2.0_kind_noahmp * delta_Sy)
    Sy = max(1.0e-6_kind_noahmp, Sy)   ! clamp to avoid division by zero
  end function MicroTopoSpecificYield

  !=======================================================================
  ! MicroTopoCorrection: main entry point called from SoilWaterMainMod
  ! Computes FloodedFrac, f_soil (= f_part), and SurfaceWaterStorage
  !=======================================================================
  subroutine MicroTopoCorrection(noahmp)
    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    real(kind=kind_noahmp) :: Sy_soil, ff

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSoilLayer      => noahmp%config%domain%NumSoilLayer              ,& ! in
              DepthSoilLayer    => noahmp%config%domain%DepthSoilLayer            ,& ! in
              SoilMoistureSat   => noahmp%water%param%SoilMoistureSat             ,& ! in
              SoilMatPotentialSat => noahmp%water%param%SoilMatPotentialSat       ,& ! in
              SoilExpCoeffB     => noahmp%water%param%SoilExpCoeffB               ,& ! in
              WaterTableDepth   => noahmp%water%state%WaterTableDepth             ,& ! in
              f_soil            => noahmp%water%state%f_soil                      ,& ! out, Sy-weighted partition fraction
              FloodedFrac       => noahmp%water%state%FloodedFrac                 & ! out, flooded fraction
             )
! ----------------------------------------------------------------------

    ! Flooded fraction from Gaussian CDF
    ff = FloodedFraction(WaterTableDepth, sigma_z)
    FloodedFrac = ff

    ! Sy-weighted flux partition
    Sy_soil = MicroTopoSpecificYield(WaterTableDepth, NumSoilLayer, DepthSoilLayer, &
                                      SoilMoistureSat(1), abs(SoilMatPotentialSat(1)), &
                                      SoilExpCoeffB(1), sigma_z)

    ! f_part = ((1-FF)*Sy_soil) / ((1-FF)*Sy_soil + FF)
    if (ff >= (1.0_kind_noahmp - 1.0e-10_kind_noahmp)) then
       f_soil = 0.0_kind_noahmp
    else
       f_soil = ((1.0_kind_noahmp - ff) * Sy_soil) / &
                ((1.0_kind_noahmp - ff) * Sy_soil + ff)
       f_soil = max(0.0_kind_noahmp, min(1.0_kind_noahmp, f_soil))
    endif

    end associate
  end subroutine MicroTopoCorrection

end module MicroTopoCorrectionMod
