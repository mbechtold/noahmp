module SoilWaterTranspirationMod
 
!!! compute soil water transpiration factor that will be used for 
!!! stomata resistance and evapotranspiration calculations
!!! Includes peatland-specific drought and waterlogging stress conditions
 
  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use PeatlandPhysicsMod,                only : ApplyPeatlandPhysics
 
  implicit none
 
contains
 
  subroutine SoilWaterTranspiration(noahmp)
 
! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Modified for Peatland transpiration and waterlogging stress (Chakraborty et al. 2026)
! -------------------------------------------------------------------------
 
    implicit none
 
! in & out variables
    type(noahmp_type), intent(inout) :: noahmp
 
! local variables
    integer                          :: IndSoil       ! loop index
    real(kind=kind_noahmp)           :: SoilWetFac    ! temporary variable
    real(kind=kind_noahmp)           :: MinThr        ! minimum threshold to prevent divided by zero
    real(kind=kind_noahmp)           :: F_wilt        ! PEAT drought stress severity (0 to 1)
    real(kind=kind_noahmp)           :: F_log         ! PEAT waterlogging stress factor
!New variables
    real(kind=kind_noahmp)           :: BetaDrought   ! transpiration reduction due to drought (deep water table)
    real(kind=kind_noahmp)           :: BetaWaterlog  ! transpiration reduction due to waterlogging (shallow WTD)
    real(kind=kind_noahmp)           :: BetaPeat      ! final combined peatland stress factor
    real(kind=kind_noahmp)           :: RootFracSum   ! temporary variable used for normalization
    real(kind=kind_noahmp)           :: RootFracBase  ! baseline root distribution weight for each soil layer
 
! --------------------------------------------------------------------
    associate(                                                                             &
              SurfaceType               => noahmp%config%domain%SurfaceType               ,& ! in,  surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer    => noahmp%config%domain%ThicknessSnowSoilLayer    ,& ! in,  thickness of snow/soil layers [m]
              DepthSoilLayer            => noahmp%config%domain%DepthSoilLayer            ,& ! in,  depth [m] of layer-bottom from soil surface
              OptSoilWaterTranspiration => noahmp%config%nmlist%OptSoilWaterTranspiration ,& ! in,  option for soil moisture factor for stomatal resistance & ET
              OptPeatlandPhysics        => noahmp%config%nmlist%OptPeatlandPhysics        ,& ! in,  option for peatland physics
              NumSoilLayerRoot          => noahmp%water%param%NumSoilLayerRoot            ,& ! in,  number of soil layers with root present
              SoilMoistureWilt          => noahmp%water%param%SoilMoistureWilt            ,& ! in,  wilting point soil moisture [m3/m3]
              SoilMoistureFieldCap      => noahmp%water%param%SoilMoistureFieldCap        ,& ! in,  reference soil moisture (field capacity) [m3/m3]
              SoilMatPotentialWilt      => noahmp%water%param%SoilMatPotentialWilt        ,& ! in,  soil metric potential for wilting point [m]
              SoilMatPotentialSat       => noahmp%water%param%SoilMatPotentialSat         ,& ! in,  saturated soil matric potential [m]
              SoilMoistureSat           => noahmp%water%param%SoilMoistureSat             ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilExpCoeffB             => noahmp%water%param%SoilExpCoeffB               ,& ! in,  soil B parameter
              SoilLiqWater              => noahmp%water%state%SoilLiqWater                ,& ! in,  soil water content [m3/m3]
              SoilTranspFac             => noahmp%water%state%SoilTranspFac               ,& ! out, soil water transpiration factor (0 to 1)
              SoilTranspFacAcc          => noahmp%water%state%SoilTranspFacAcc            ,& ! out, accumulated soil water transpiration factor (0 to 1)
              SoilMatPotential          => noahmp%water%state%SoilMatPotential            ,& ! out, soil matrix potential [m]
              WaterTableDepth           => noahmp%water%state%WaterTableDepth              & ! in,  water table depth [m]
             )
! ----------------------------------------------------------------------
 
    ! soil moisture factor controlling stomatal resistance and evapotranspiration
    MinThr           = 1.0e-6
    SoilTranspFacAcc = 0.0
    RootFracSum      = 0.0
    BetaPeat         = 1.0
 
    ! Set peatland-specific transpiration option if peatland physics is enabled
    if ( OptPeatlandPhysics == 1 ) then
       OptSoilWaterTranspiration = 4
    endif
 
    ! only for soil point
    if ( SurfaceType ==1 ) then
       do IndSoil = 1, NumSoilLayerRoot
          if ( OptSoilWaterTranspiration == 1 ) then  ! Noah
             SoilWetFac                = (SoilLiqWater(IndSoil) - SoilMoistureWilt(IndSoil)) / &
                                         (SoilMoistureFieldCap(IndSoil) - SoilMoistureWilt(IndSoil))
          endif
          if ( OptSoilWaterTranspiration == 2 ) then  ! CLM
             SoilMatPotential(IndSoil) = max(SoilMatPotentialWilt, -SoilMatPotentialSat(IndSoil) * &
                                            (max(0.01,SoilLiqWater(IndSoil))/SoilMoistureSat(IndSoil)) ** &
                                            (-SoilExpCoeffB(IndSoil)))
             SoilWetFac                = (1.0 - SoilMatPotential(IndSoil)/SoilMatPotentialWilt) / &
                                         (1.0 + SoilMatPotentialSat(IndSoil)/SoilMatPotentialWilt)
          endif
          if ( OptSoilWaterTranspiration == 3 ) then  ! SSiB
             SoilMatPotential(IndSoil) = max(SoilMatPotentialWilt, -SoilMatPotentialSat(IndSoil) * &
                                            (max(0.01,SoilLiqWater(IndSoil))/SoilMoistureSat(IndSoil)) ** &
                                            (-SoilExpCoeffB(IndSoil)))
             SoilWetFac                = 1.0 - exp(-5.8*(log(SoilMatPotentialWilt/SoilMatPotential(IndSoil))))
          endif
          if ( OptSoilWaterTranspiration == 4 ) then  ! PEAT
 
             ! Drought stress severity: 0 = no stress, 1 = full stress
             if (WaterTableDepth < 0.3) then
                F_wilt = 0.0
             else if (WaterTableDepth >= 0.3 .and. WaterTableDepth < 1.15) then
                F_wilt = 1.18 * WaterTableDepth - 0.35
             else
                F_wilt = 1.0
             end if
             F_wilt = max(0.0, min(1.0, F_wilt))
 
             !------------------------Separate block added for option 4------------------------------
             ! Convert drought stress severity to transpiration reduction factor
             BetaDrought = 1.0 - F_wilt
             BetaDrought = max(0.0, min(1.0, BetaDrought))
             !----------------------------------End Block---------------------------------------------
 
             ! Waterlogging reduction factor: 1 = no stress, 0 = full stress
             if (WaterTableDepth >= 0.29) then
                F_log = 1.0
             else if (WaterTableDepth >= -0.35 .and. WaterTableDepth < 0.29) then
                F_log = 1.0 - max(0.0, min(0.95, (0.29 - WaterTableDepth) / 0.64))
             else
                F_log = 0.0
             end if
             F_log = max(0.0, min(1.0, F_log))
             ! Combine both drought and waterlogging stress factors
             !SoilWetFac = (1.0 - F_wilt) * F_log
         !endif
         !SoilWetFac = min(1.0, max(0.0, SoilWetFac))
 
         !SoilTranspFac(IndSoil) = max(MinThr, ThicknessSnowSoilLayer(IndSoil) / &
         !                              (-DepthSoilLayer(NumSoilLayerRoot)) * SoilWetFac)
         !SoilTranspFacAcc = SoilTranspFacAcc + SoilTranspFac(IndSoil)
 
             !------------------------Separate block added for option 4------------------------------
             BetaWaterlog = F_log
 
             ! Final combined peatland stress factor
             BetaPeat = BetaDrought * BetaWaterlog
             BetaPeat = max(MinThr, min(1.0, BetaPeat))
 
             ! Use only root distribution inside SoilTranspFac for PEAT.
             ! The bulk WTD stress is stored later in SoilTranspFacAcc.
             RootFracBase              = max(MinThr, ThicknessSnowSoilLayer(IndSoil) / &
                                            (-DepthSoilLayer(NumSoilLayerRoot)))
             SoilTranspFac(IndSoil)    = RootFracBase
             RootFracSum               = RootFracSum + SoilTranspFac(IndSoil)
          else
             SoilWetFac                   = min(1.0, max(0.0,SoilWetFac))
 
             SoilTranspFac(IndSoil)       = max(MinThr, ThicknessSnowSoilLayer(IndSoil) / &
                                               (-DepthSoilLayer(NumSoilLayerRoot)) * SoilWetFac)
             SoilTranspFacAcc             = SoilTranspFacAcc + SoilTranspFac(IndSoil)
          endif
         !----------------------------------End Block---------------------------------------------
       enddo
       !------------------------Separate block added for option 4------------------------------
       if ( OptSoilWaterTranspiration == 4 ) then
          SoilTranspFacAcc = BetaPeat
          RootFracSum = max(MinThr, RootFracSum)
          SoilTranspFac(1:NumSoilLayerRoot) = SoilTranspFac(1:NumSoilLayerRoot) / RootFracSum
       else
       !----------------------------------End Block---------------------------------------------
          SoilTranspFacAcc = max(MinThr, SoilTranspFacAcc)
          SoilTranspFac(1:NumSoilLayerRoot) = SoilTranspFac(1:NumSoilLayerRoot) / SoilTranspFacAcc
       endif
    endif
 
    end associate
 
  end subroutine SoilWaterTranspiration
 
end module SoilWaterTranspirationMod