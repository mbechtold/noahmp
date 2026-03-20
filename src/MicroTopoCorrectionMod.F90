module MicroTopoCorrectionMod

!!! Gaussian microtopography correction for peatland soil water (Dettmann & Bechtold 2015)
!!! Replaces PEATCLSM empirical formulation with physically consistent approach.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

  ! Module-level parameters
  real(kind=kind_noahmp), parameter :: sigma_z       = 0.16_kind_noahmp
  real(kind=kind_noahmp), parameter :: z_elev_max    = 1.0_kind_noahmp
  real(kind=kind_noahmp), parameter :: sqrt2         = 1.4142135623730951_kind_noahmp
  real(kind=kind_noahmp), parameter :: inv_sqrt2pi   = 0.3989422804014327_kind_noahmp

contains

  ! ====================================================================
  ! GaussianCDF: Phi(x) = 0.5 * erfc(-x / sqrt(2))
  ! ====================================================================
  pure function GaussianCDF(x) result(phi)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x
    real(kind=kind_noahmp) :: phi
    phi = 0.5_kind_noahmp * erfc(-x / sqrt2)
  end function GaussianCDF

  ! ====================================================================
  ! GaussianPDF: phi(x) = exp(-x^2/2) / sqrt(2*pi)
  ! ====================================================================
  pure function GaussianPDF(x) result(phi)
    implicit none
    real(kind=kind_noahmp), intent(in) :: x
    real(kind=kind_noahmp) :: phi
    phi = inv_sqrt2pi * exp(-0.5_kind_noahmp * x * x)
  end function GaussianPDF

  ! ====================================================================
  ! FloodedFracFromWTD: fraction of landscape flooded
  ! ====================================================================
  pure function FloodedFracFromWTD(WTD) result(ff)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD
    real(kind=kind_noahmp) :: ff
    ff = GaussianCDF(-WTD / sigma_z)
    ff = max(0.0_kind_noahmp, min(1.0_kind_noahmp, ff))
  end function FloodedFracFromWTD

  ! ====================================================================
  ! SoilFractionAtDepth: fraction of landscape with soil at depth d
  ! ====================================================================
  pure function SoilFractionAtDepth(d) result(fs)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d
    real(kind=kind_noahmp) :: fs
    fs = GaussianCDF(d / sigma_z)
  end function SoilFractionAtDepth

  ! ====================================================================
  ! SurfaceWaterStorage_mm: total surface water storage [mm]
  !   S = 1000 * [sigma * phi_std(-WTD/sigma) - WTD * Phi(-WTD/sigma)]
  ! ====================================================================
  pure function SurfaceWaterStorage_mm(WTD) result(Smm)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD
    real(kind=kind_noahmp) :: Smm, arg
    arg = -WTD / sigma_z
    Smm = 1000.0_kind_noahmp * (sigma_z * GaussianPDF(arg) - WTD * GaussianCDF(arg))
    Smm = max(0.0_kind_noahmp, Smm)
  end function SurfaceWaterStorage_mm

  ! ====================================================================
  ! CampbellEquilibriumTheta: equilibrium moisture at depth d
  ! ====================================================================
  pure function CampbellEquilibriumTheta(d, WTD, thetas, ae, bb) result(theta)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d, WTD, thetas, ae, bb
    real(kind=kind_noahmp) :: theta, psi_abs
    psi_abs = WTD - d
    if (psi_abs <= ae) then
       theta = thetas
    else
       theta = thetas * (psi_abs / ae) ** (-1.0_kind_noahmp / bb)
       theta = max(0.0_kind_noahmp, min(thetas, theta))
    endif
  end function CampbellEquilibriumTheta

  ! ====================================================================
  ! MicroTopoEquilibriumMoisture: area-averaged equilibrium moisture
  !   for a layer from d_top to d_bot, integrating over microtopography
  !   theta_eq_k = (1/dz) * integral Phi(d/sigma) * theta_eq(d,WTD) dd
  ! ====================================================================
  function MicroTopoEquilibriumMoisture(WTD, d_top, d_bot, thetas, ae, bb) result(theta_eq)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD, d_top, d_bot, thetas, ae, bb
    real(kind=kind_noahmp) :: theta_eq
    integer :: j, Nsub
    real(kind=kind_noahmp) :: dz_sub, d_mid, dz_layer

    dz_layer = d_bot - d_top
    if (dz_layer <= 0.0_kind_noahmp) then
       theta_eq = 0.0_kind_noahmp
       return
    endif
    Nsub = max(20, nint(dz_layer / 0.005_kind_noahmp))
    dz_sub = dz_layer / real(Nsub, kind_noahmp)

    theta_eq = 0.0_kind_noahmp
    do j = 1, Nsub
       d_mid = d_top + (real(j, kind_noahmp) - 0.5_kind_noahmp) * dz_sub
       theta_eq = theta_eq + SoilFractionAtDepth(d_mid) * &
                             CampbellEquilibriumTheta(d_mid, WTD, thetas, ae, bb) * dz_sub
    enddo
    theta_eq = theta_eq / dz_layer
  end function MicroTopoEquilibriumMoisture

  ! ====================================================================
  ! MicroTopoEquilibriumMoisture1D: flat-column equilibrium moisture
  !   (no microtopography weighting, for 1D Richards)
  ! ====================================================================
  function MicroTopoEquilibriumMoisture1D(WTD, d_top, d_bot, thetas, ae, bb) result(theta_eq)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD, d_top, d_bot, thetas, ae, bb
    real(kind=kind_noahmp) :: theta_eq
    integer :: j, Nsub
    real(kind=kind_noahmp) :: dz_sub, d_mid, dz_layer

    dz_layer = d_bot - d_top
    if (dz_layer <= 0.0_kind_noahmp) then
       theta_eq = 0.0_kind_noahmp
       return
    endif
    Nsub = max(20, nint(dz_layer / 0.005_kind_noahmp))
    dz_sub = dz_layer / real(Nsub, kind_noahmp)

    theta_eq = 0.0_kind_noahmp
    do j = 1, Nsub
       d_mid = d_top + (real(j, kind_noahmp) - 0.5_kind_noahmp) * dz_sub
       theta_eq = theta_eq + CampbellEquilibriumTheta(d_mid, WTD, thetas, ae, bb) * dz_sub
    enddo
    theta_eq = theta_eq / dz_layer
  end function MicroTopoEquilibriumMoisture1D

  ! ====================================================================
  ! LayerAverageSoilFraction: depth-integrated soil fraction for a layer
  ! ====================================================================
  function LayerAverageSoilFraction(d_top, d_bot) result(f_soil)
    implicit none
    real(kind=kind_noahmp), intent(in) :: d_top, d_bot
    real(kind=kind_noahmp) :: f_soil
    integer :: j, Nsub
    real(kind=kind_noahmp) :: dz_sub, d_mid, dz_layer

    dz_layer = d_bot - d_top
    if (dz_layer <= 0.0_kind_noahmp) then
       f_soil = 0.0_kind_noahmp
       return
    endif
    Nsub = max(20, nint(dz_layer / 0.005_kind_noahmp))
    dz_sub = dz_layer / real(Nsub, kind_noahmp)

    f_soil = 0.0_kind_noahmp
    do j = 1, Nsub
       d_mid = d_top + (real(j, kind_noahmp) - 0.5_kind_noahmp) * dz_sub
       f_soil = f_soil + SoilFractionAtDepth(d_mid) * dz_sub
    enddo
    f_soil = f_soil / dz_layer
  end function LayerAverageSoilFraction

  ! ====================================================================
  ! MicroTopoStorageDeficit: total moisture deficit [m] for given WTD
  !   With below-model extension for deep WTD
  ! ====================================================================
  function MicroTopoStorageDeficit(WTD, NumSoilLayer, DepthSoilLayer, &
                                    ThicknessSoilLayer, thetas, ae, bb) result(deficit)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD, thetas, ae, bb
    integer, intent(in) :: NumSoilLayer
    real(kind=kind_noahmp), intent(in) :: DepthSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp), intent(in) :: ThicknessSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp) :: deficit
    integer :: k, j, Nsub_ext
    real(kind=kind_noahmp) :: d_top, d_bot, theta_eq_k, theta_sat_k
    real(kind=kind_noahmp) :: depth_bot, dz_ext, dz_sub, d_mid, theta_ext

    deficit = 0.0_kind_noahmp

    do k = 1, NumSoilLayer
       if (k == 1) then
          d_top = 0.0_kind_noahmp
       else
          d_top = abs(DepthSoilLayer(k-1))
       endif
       d_bot = abs(DepthSoilLayer(k))
       theta_eq_k = MicroTopoEquilibriumMoisture(WTD, d_top, d_bot, thetas, ae, bb)
       theta_sat_k = thetas * LayerAverageSoilFraction(d_top, d_bot)
       deficit = deficit + (theta_sat_k - theta_eq_k) * ThicknessSoilLayer(k)
    enddo

    ! Below-model extension for deep WTD (flat soil, f_soil ~ 1.0)
    depth_bot = abs(DepthSoilLayer(NumSoilLayer))
    if (WTD > depth_bot) then
       dz_ext = WTD - depth_bot
       Nsub_ext = max(20, nint(dz_ext / 0.01_kind_noahmp))
       dz_sub = dz_ext / real(Nsub_ext, kind_noahmp)
       do j = 1, Nsub_ext
          d_mid = depth_bot + (real(j, kind_noahmp) - 0.5_kind_noahmp) * dz_sub
          theta_ext = CampbellEquilibriumTheta(d_mid, WTD, thetas, ae, bb)
          deficit = deficit + (thetas - theta_ext) * dz_sub
       enddo
    endif
  end function MicroTopoStorageDeficit

  ! ====================================================================
  ! MicroTopoSpecificYield: Sy = dDeficit/dWTD (central difference)
  ! ====================================================================
  function MicroTopoSpecificYield(WTD, NumSoilLayer, DepthSoilLayer, &
                                   ThicknessSoilLayer, thetas, ae, bb) result(Sy)
    implicit none
    real(kind=kind_noahmp), intent(in) :: WTD, thetas, ae, bb
    integer, intent(in) :: NumSoilLayer
    real(kind=kind_noahmp), intent(in) :: DepthSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp), intent(in) :: ThicknessSoilLayer(NumSoilLayer)
    real(kind=kind_noahmp) :: Sy
    real(kind=kind_noahmp), parameter :: delta = 0.001_kind_noahmp

    Sy = (MicroTopoStorageDeficit(WTD + delta, NumSoilLayer, DepthSoilLayer, &
                                   ThicknessSoilLayer, thetas, ae, bb) - &
          MicroTopoStorageDeficit(WTD - delta, NumSoilLayer, DepthSoilLayer, &
                                   ThicknessSoilLayer, thetas, ae, bb)) &
         / (2.0_kind_noahmp * delta)
    Sy = max(0.001_kind_noahmp, Sy)
  end function MicroTopoSpecificYield

  ! ====================================================================
  ! MicroTopoCorrection: main subroutine
  !   Computes FloodedFrac, f_part, f_soil_k
  ! ====================================================================
  subroutine MicroTopoCorrection(noahmp)
    implicit none
    type(noahmp_type), intent(inout) :: noahmp

    real(kind=kind_noahmp) :: Sy_soil, ff, d_top, d_bot
    integer :: k

    associate(                                                                             &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer              ,&
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer            ,&
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer    ,&
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat             ,&
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat         ,&
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB               ,&
              WaterTableDepth        => noahmp%water%state%WaterTableDepth             ,&
              FloodedFrac            => noahmp%water%state%FloodedFrac                 ,&
              f_part                 => noahmp%water%state%f_part                      ,&
              f_soil_k               => noahmp%water%state%f_soil_k                     &
             )

    ff = FloodedFracFromWTD(WaterTableDepth)
    FloodedFrac = ff

    do k = 1, NumSoilLayer
       if (k == 1) then
          d_top = 0.0_kind_noahmp
       else
          d_top = abs(DepthSoilLayer(k-1))
       endif
       d_bot = abs(DepthSoilLayer(k))
       f_soil_k(k) = LayerAverageSoilFraction(d_top, d_bot)
    enddo

    Sy_soil = MicroTopoSpecificYield(WaterTableDepth, NumSoilLayer, DepthSoilLayer, &
                  ThicknessSnowSoilLayer, SoilMoistureSat(1), &
                  abs(SoilMatPotentialSat(1)), SoilExpCoeffB(1))

    if (ff > (1.0_kind_noahmp - 1.0e-6_kind_noahmp)) then
       f_part = 0.0_kind_noahmp
    else
       f_part = ((1.0_kind_noahmp - ff) * Sy_soil) / &
                ((1.0_kind_noahmp - ff) * Sy_soil + ff)
       f_part = max(0.0_kind_noahmp, min(1.0_kind_noahmp, f_part))
    endif

    end associate
  end subroutine MicroTopoCorrection

end module MicroTopoCorrectionMod
