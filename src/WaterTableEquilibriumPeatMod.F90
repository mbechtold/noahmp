module WaterTableEquilibriumPeatMod

!!! Calculate equilibrium water table depth for peatlands
!!!
!!! Both the target deficit (from NoahMP soil moisture) and the equilibrium
!!! deficit function (SoilDeficitMicroTopo) are computed over the same
!!! microtopography-consistent domain, extending from the top of the
!!! hummocks (+z_trunc elevation) to the soil column bottom (-z_col_bot):
!!!
!!!   - Hummock zone (0 to +z_trunc elevation):
!!!     soil fraction = 1 - Phi(z/sigma), weighted by Gaussian CDF.
!!!     Assigned the top-layer soil moisture for the deficit target.
!!!   - Below mean surface (0 to z_col_bot depth):
!!!     soil fraction varies via Gaussian CDF in the microtopo zone and
!!!     equals 1 below z_trunc depth.  Integrated via Gauss-Legendre.
!!!
!!! This ensures domain-consistency: deficit_target and SoilDeficitMicroTopo
!!! both integrate over the same spatial domain with identical weighting.
!!!
!!! WaterTableDepth convention: positive downward (NoahMP standard)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use PeatMicroTopoMod, only : SoilDeficitMicroTopo, EffSoilThickMicroTopo, &
                                z_trunc, InitGaussLegendre, gl_initialized
  
  implicit none

contains

  subroutine WaterTableEquilibriumPeat(noahmp)
  
! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: ZWTEQ
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactored:    C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! This version:  Microtopography-consistent deficit bisection (Bechtold, 2026)
!                Both deficit_target and equilibrium function use the same
!                spatial domain: hummock zone + microtopo zone + deep zone
! ----------------------------------------------------------------------------------------

    implicit none
    
    type(noahmp_type), intent(inout) :: noahmp

    integer                          :: i, iter
    real(kind=kind_noahmp)           :: ae, bb, thetas
    real(kind=kind_noahmp)           :: deficit_target, deficit_mid
    real(kind=kind_noahmp)           :: zwt_lo, zwt_hi, zwt_mid
    real(kind=kind_noahmp)           :: z_col_bot               ! column bottom depth [m]
    real(kind=kind_noahmp)           :: z_depth_top, z_depth_bot ! layer depth bounds [m]
    real(kind=kind_noahmp)           :: eff_thick               ! effective soil thickness [m]
    real(kind=kind_noahmp)           :: V_hummock               ! hummock soil thickness [m]
    real(kind=kind_noahmp), parameter:: tol_def  = 1.0e-6_kind_noahmp
    real(kind=kind_noahmp), parameter:: tol_zwt  = 1.0e-4_kind_noahmp
    integer,          parameter      :: max_iter = 60
! -----------------------------------------------------------------------------------------------------------------------------
    associate(                                                                        &
      NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in
      DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in  layer-bottom depths [m], negative downward
      ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in
      SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! in  [m3/m3]
      SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in  theta_s [m3/m3]
      SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in  psi_e (air-entry), typically negative [m]
      SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in  Campbell b
      WaterTableDepth        => noahmp%water%state%WaterTableDepth            & ! out z_wt from surface [m], positive downward
    )
! ------------------------------------------------------------------------------------------------------------------------------

      ! Ensure Gauss-Legendre quadrature is initialized
      if (.not. gl_initialized) call InitGaussLegendre()

      ! Constants for the single soil type used (peat)
      thetas = SoilMoistureSat(1)
      ae     = abs(SoilMatPotentialSat(1))   ! air-entry suction head [m], positive
      bb     = SoilExpCoeffB(1)
      z_col_bot = abs(DepthSoilLayer(NumSoilLayer))  ! e.g. 2.0 m

      ! ================================================================
      ! Compute microtopography-consistent deficit target from NoahMP state
      ! ================================================================
      deficit_target = 0.0_kind_noahmp

      ! 1. Hummock zone (above mean surface, depth from -z_trunc to 0)
      !    Assign top-layer soil moisture (hummock peaks extend only a few cm;
      !    most of the hummock soil volume is within the top-layer depth range).
      V_hummock = EffSoilThickMicroTopo(-z_trunc, 0.0_kind_noahmp)
      deficit_target = deficit_target + (thetas - SoilLiqWater(1)) * V_hummock

      ! 2. NoahMP layers (each weighted by microtopo soil fraction)
      do i = 1, NumSoilLayer
        if (i == 1) then
          z_depth_top = 0.0_kind_noahmp
        else
          z_depth_top = abs(DepthSoilLayer(i-1))
        endif
        z_depth_bot = abs(DepthSoilLayer(i))
        eff_thick = EffSoilThickMicroTopo(z_depth_top, z_depth_bot)
        deficit_target = deficit_target + (thetas - SoilLiqWater(i)) * eff_thick
      end do

      ! ================================================================
      ! Bisect on WaterTableDepth using SoilDeficitMicroTopo
      ! which integrates over the same domain as deficit_target:
      !   hummock zone + microtopo zone + deep zone down to z_col_bot
      ! SoilDeficitMicroTopo takes z_wt in D&B convention (positive up),
      ! so we pass -WTD.
      ! ================================================================
      ! Bisection range:
      !   zwt_lo = -z_trunc (WT at top of hummocks, deficit ~ 0)
      !   zwt_hi = z_col_bot (WT at column bottom, maximum deficit)
      zwt_lo = -z_trunc
      zwt_hi =  z_col_bot

      ! If deficit_target is effectively zero, WT is at or above all hummocks
      if (deficit_target <= tol_def) then
        WaterTableDepth = zwt_lo
      else
        ! Check upper bracket
        deficit_mid = SoilDeficitMicroTopo(-zwt_hi, thetas, ae, bb, z_col_bot)
        if (deficit_mid < deficit_target) then
          ! Target exceeds capacity: deepest water table
          WaterTableDepth = zwt_hi
        else
          ! Bisection
          do iter = 1, max_iter
            zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
            deficit_mid = SoilDeficitMicroTopo(-zwt_mid, thetas, ae, bb, z_col_bot)

            if (abs(deficit_mid - deficit_target) <= tol_def .or. &
                (zwt_hi - zwt_lo) <= tol_zwt) then
              WaterTableDepth = zwt_mid
              exit
            end if

            if (deficit_mid > deficit_target) then
              ! Deficit too large (too deep) -> move shallower
              zwt_hi = zwt_mid
            else
              ! Deficit too small (too shallow) -> move deeper
              zwt_lo = zwt_mid
            end if

            if (iter == max_iter) WaterTableDepth = zwt_mid
          end do
        end if
      end if

    end associate

  end subroutine WaterTableEquilibriumPeat

end module WaterTableEquilibriumPeatMod
