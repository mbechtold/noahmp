module WaterTableEquilibriumPeatMod

!!! Calculate equilibrium water table depth for peatlands
!!! Now includes microtopography effects following Dettmann & Bechtold (2015)
!!! Hydrological Processes, DOI: 10.1002/hyp.10475
!!!
!!! The soil deficit calculation accounts for the Gaussian distribution of
!!! surface elevations: hummocks add soil above mean surface, hollows remove
!!! soil below mean surface. Surface water in hollows is also tracked.
!!!
!!! WaterTableDepth convention: positive downward (NoahMP standard)
!!! Internally: z_wt = -WaterTableDepth (positive upward, D&B convention)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use PeatMicroTopoMod, only : SoilDeficitMicroTopo, &
                                InitGaussLegendre, gl_initialized
  
  implicit none

contains

  subroutine WaterTableEquilibriumPeat(noahmp)
  
! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: ZWTEQ
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactored:    C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! This version:  Microtopography-aware equilibrium using Dettmann & Bechtold (2015)
!                Bisection on soil deficit (Bechtold, 2026)
! ----------------------------------------------------------------------------------------

    implicit none
    
    type(noahmp_type), intent(inout) :: noahmp

    integer                          :: i, iter
    real(kind=kind_noahmp)           :: ae, bb, thetas
    real(kind=kind_noahmp)           :: deficit_target, deficit_mid
    real(kind=kind_noahmp)           :: zwt_lo, zwt_hi, zwt_mid
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

      ! Compute target deficit from current coarse NoahMP profile
      ! This is the amount of water below full saturation in the modeled layers.
      ! deficit_target >= 0 always because SoilLiqWater <= theta_s (clamped
      ! in SoilWaterMainMod). With microtopography extending to +1 m,
      ! there is always unsaturated soil in the hummocks, so deficit_target > 0
      ! under any realistic condition.
      deficit_target = 0.0_kind_noahmp
      do i = 1, NumSoilLayer
        deficit_target = deficit_target + (thetas - SoilLiqWater(i)) * ThicknessSnowSoilLayer(i)
      end do

      ! Find water table using microtopography-aware deficit bisection.
      ! SoilDeficitMicroTopo accounts for:
      !   - Reduced soil volume in hollows (below mean surface)
      !   - Extra soil volume in hummocks (above mean surface)
      !   - Hydrostatic equilibrium with Campbell retention
      ! Bisect on WaterTableDepth (positive downward):
      !   zwt_lo = -1.0 (z_wt = +1.0 m, top of microtopography, deficit ~ 0)
      !   zwt_hi = +3.0 (z_wt = -3.0 m, deep water table, maximum deficit)
      zwt_lo = -1.0_kind_noahmp
      zwt_hi =  3.0_kind_noahmp

      ! Check bracket
      deficit_mid = SoilDeficitMicroTopo(-zwt_hi, thetas, ae, bb)
      if (deficit_mid < deficit_target) then
        ! Target exceeds capacity: deepest water table
        WaterTableDepth = zwt_hi
      else
        deficit_mid = SoilDeficitMicroTopo(-zwt_lo, thetas, ae, bb)
        if (deficit_mid > deficit_target) then
          ! Very little deficit: water table near top of microtopography
          WaterTableDepth = zwt_lo
        else
          ! Bisection
          do iter = 1, max_iter
            zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
            deficit_mid = SoilDeficitMicroTopo(-zwt_mid, thetas, ae, bb)

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
