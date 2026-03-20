module WaterTableEquilibriumPeatMod

!!! Diagnose water table depth from soil moisture via microtopographic
!!! deficit bisection (replaces fine-layer approach).

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use MicroTopoCorrectionMod, only: MicroTopoStorageDeficit, &
                                     MicroTopoEquilibriumMoisture, &
                                     LayerAverageSoilFraction

  implicit none

contains

  subroutine WaterTableEquilibriumPeat(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    real(kind=kind_noahmp), parameter :: tol_def  = 1.0e-6_kind_noahmp
    real(kind=kind_noahmp), parameter :: tol_zwt  = 1.0e-4_kind_noahmp
    integer,                parameter :: max_iter = 60

    integer              :: k, iter
    real(kind=kind_noahmp) :: deficit_target, deficit_mid
    real(kind=kind_noahmp) :: zwt_lo, zwt_hi, zwt_mid, zmax
    real(kind=kind_noahmp) :: ae, bb, thetas
    real(kind=kind_noahmp) :: d_top, d_bot, dz_k, f_soil_k_local, theta_sat_k

    associate(                                                                        &
      NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,&
      DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,&
      SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,&
      SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,&
      SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,&
      SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,&
      WaterTableDepth        => noahmp%water%state%WaterTableDepth           &
    )

    thetas = SoilMoistureSat(1)
    ae     = abs(SoilMatPotentialSat(1))
    bb     = SoilExpCoeffB(1)

    ! Compute target deficit from current area-averaged SoilLiqWater
    ! deficit_target = sum over layers of (theta_s * f_soil_k - theta_k) * dz_k
    deficit_target = 0.0_kind_noahmp
    do k = 1, NumSoilLayer
       if (k == 1) then
          d_top = 0.0_kind_noahmp
       else
          d_top = abs(DepthSoilLayer(k-1))
       endif
       d_bot = abs(DepthSoilLayer(k))
       dz_k = d_bot - d_top
       f_soil_k_local = LayerAverageSoilFraction(d_top, d_bot)
       theta_sat_k = thetas * f_soil_k_local
       deficit_target = deficit_target + &
                        (theta_sat_k - SoilLiqWater(k)) * dz_k
    enddo

    if (deficit_target <= 0.0_kind_noahmp) then
       WaterTableDepth = 0.0_kind_noahmp
    else
       ! Bracket: WTD between 0 and 3x soil bottom
       zwt_lo = 0.0_kind_noahmp
       zmax = 3.0_kind_noahmp * abs(DepthSoilLayer(NumSoilLayer))
       zwt_hi = max(zmax, 1.0_kind_noahmp)

       if (MicroTopoStorageDeficit(zwt_hi, NumSoilLayer, DepthSoilLayer, &
           thetas, ae, bb) < deficit_target) then
          WaterTableDepth = zwt_hi
       else
          ! Bisection
          do iter = 1, max_iter
             zwt_mid = 0.5_kind_noahmp * (zwt_lo + zwt_hi)
             deficit_mid = MicroTopoStorageDeficit(zwt_mid, NumSoilLayer, &
                           DepthSoilLayer, thetas, ae, bb)

             if (abs(deficit_mid - deficit_target) <= tol_def .or. &
                 (zwt_hi - zwt_lo) <= tol_zwt) then
                WaterTableDepth = zwt_mid
                exit
             endif

             if (deficit_mid > deficit_target) then
                zwt_hi = zwt_mid
             else
                zwt_lo = zwt_mid
             endif

             if (iter == max_iter) WaterTableDepth = zwt_mid
          enddo
       endif
    endif

    end associate
  end subroutine WaterTableEquilibriumPeat

end module WaterTableEquilibriumPeatMod
