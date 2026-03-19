module WaterTableEquilibriumPeatMod

!!! Diagnose equilibrium water table depth from soil moisture state
!!! using multi-column Gaussian microtopography deficit (Dettmann & Bechtold 2015)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use MicroTopoCorrectionMod, only : MicroTopoStorageDeficitExtended
  
  implicit none

contains

  subroutine WaterTableEquilibriumPeat(noahmp)
  
! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: ZWTEQ
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactored:    C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! This version:  multi-column Gaussian microtopo deficit bisection (2025)
! ----------------------------------------------------------------------------------------

    implicit none
    
    type(noahmp_type), intent(inout) :: noahmp

    integer                          :: i, iter
    real(kind=kind_noahmp)           :: deficit_target, deficit_mid
    real(kind=kind_noahmp)           :: zwt_lo, zwt_hi, zwt_mid, zmax
    real(kind=kind_noahmp)           :: ae, bb, thetas
    real(kind=kind_noahmp), parameter:: tol_def  = 1.0e-6_kind_noahmp
    real(kind=kind_noahmp), parameter:: tol_zwt  = 1.0e-4_kind_noahmp
    integer,          parameter      :: max_iter = 60
    real(kind=kind_noahmp), parameter:: WaterTableDepthMinPeat = -0.2449_kind_noahmp
! -----------------------------------------------------------------------------------------------------------------------------
    associate(                                                                        &
      NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in
      DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in  layer-bottom depths [m], negative downward
      ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in
      SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! in  [m3/m3]
      SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in  θ_s [m3/m3]
      SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in  ψ_e (air-entry), typically negative [m]
      SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in  Campbell b
      WaterTableDepth        => noahmp%water%state%WaterTableDepth            & ! out z_wt from surface [m], positive downward
    )
! ------------------------------------------------------------------------------------------------------------------------------

      ! constants for the single peat soil type
      thetas = SoilMoistureSat(1)
      ae     = abs(SoilMatPotentialSat(1))   ! air-entry suction head [m], positive
      bb     = SoilExpCoeffB(1)

      ! target deficit from current soil moisture profile
      deficit_target = 0.0_kind_noahmp
      do i = 1, NumSoilLayer
        deficit_target = deficit_target + (thetas - SoilLiqWater(i)) * ThicknessSnowSoilLayer(i)
      end do

      if (deficit_target <= 0.0_kind_noahmp) then
        ! Fully saturated or oversaturated → WT at or above surface
        WaterTableDepth = WaterTableDepthMinPeat
      else
        ! Bisection bracket
        zwt_lo = WaterTableDepthMinPeat
        zmax   = 3.0_kind_noahmp * (-DepthSoilLayer(NumSoilLayer))
        zwt_hi = zmax

        ! Check if deficit at upper bracket is already sufficient
        deficit_mid = MicroTopoStorageDeficitExtended(zwt_hi, NumSoilLayer, DepthSoilLayer, &
                                                       thetas, ae, bb, 0.16_kind_noahmp)
        if (deficit_mid < deficit_target) then
          ! Target exceeds capacity: deepest WT
          WaterTableDepth = zwt_hi
        else
          ! Bisection on multi-column deficit
          do iter = 1, max_iter
            zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
            deficit_mid = MicroTopoStorageDeficitExtended(zwt_mid, NumSoilLayer, DepthSoilLayer, &
                                                           thetas, ae, bb, 0.16_kind_noahmp)

            if (abs(deficit_mid - deficit_target) <= tol_def .or. (zwt_hi - zwt_lo) <= tol_zwt) then
              WaterTableDepth = zwt_mid
              exit
            end if

            if (deficit_mid > deficit_target) then
              zwt_hi = zwt_mid   ! too deep → move shallower
            else
              zwt_lo = zwt_mid   ! too shallow → move deeper
            end if

            if (iter == max_iter) WaterTableDepth = zwt_mid
          end do
        end if
      end if

    end associate

  end subroutine WaterTableEquilibriumPeat

end module WaterTableEquilibriumPeatMod
