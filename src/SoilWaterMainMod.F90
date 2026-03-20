module SoilWaterMainMod

!!! Main soil water module including all soil water processes & update soil moisture
!!! surface runoff, infiltration, soil water diffusion, subsurface runoff, tile drainage
!!! added peatland specific options (Chakraborty et al., 2025)
!!! Peatland hybrid two-regime microtopo (Bechtold, 2025)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use RunoffSurfaceTopModelGrdMod,       only : RunoffSurfaceTopModelGrd
  use RunoffSurfaceTopModelEquiMod,      only : RunoffSurfaceTopModelEqui
  use RunoffSurfaceFreeDrainMod,         only : RunoffSurfaceFreeDrain
  use RunoffSurfaceBatsMod,              only : RunoffSurfaceBATS
  use RunoffSurfaceTopModelMmfMod,       only : RunoffSurfaceTopModelMMF
  use RunoffSurfaceVicMod,               only : RunoffSurfaceVIC
  use RunoffSurfaceXinAnJiangMod,        only : RunoffSurfaceXinAnJiang
  use RunoffSurfaceDynamicVicMod,        only : RunoffSurfaceDynamicVic
  use RunoffSubSurfaceEquiWaterTableMod, only : RunoffSubSurfaceEquiWaterTable
  use RunoffSubSurfaceGroundWaterMod,    only : RunoffSubSurfaceGroundWater
  use RunoffSubSurfaceDrainageMod,       only : RunoffSubSurfaceDrainage
  use RunoffSubSurfaceShallowMmfMod,     only : RunoffSubSurfaceShallowWaterMMF
  use RunoffSubSurfacePeatlandMod,       only : RunoffSubSurfacePeatland
  use MicroTopoCorrectionMod,            only : MicroTopoCorrection, &
                                                 SurfaceWaterStorage_mm, &
                                                 MicroTopoEquilibriumMoisture, &
                                                 MicroTopoEquilibriumMoisture1D, &
                                                 MicroTopoSpecificYield
  use PeatlandPhysicsMod,                only : ApplyPeatlandPhysics
  use SoilWaterDiffusionRichardsMod,     only : SoilWaterDiffusionRichards
  use SoilMoistureSolverMod,             only : SoilMoistureSolver
  use TileDrainageSimpleMod,             only : TileDrainageSimple
  use TileDrainageHooghoudtMod,          only : TileDrainageHooghoudt
  use WaterTableEquilibriumMod,          only : WaterTableEquilibrium
  use WaterTableEquilibriumPeatMod,      only : WaterTableEquilibriumPeat

  implicit none

contains

  subroutine SoilWaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SOILWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Option=9; RunoffSubSurfacePeatlandMod added by Chakraborty & Bechtold (2025)
! Two-regime microtopo: Bechtold (2025)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout)  :: noahmp

! local variables
    integer                           :: LoopInd1, LoopInd2           ! loop index
    integer                           :: IndIter                      ! iteration index
    integer                           :: NumIterSoilWat               ! iteration times soil moisture
    real(kind=kind_noahmp)            :: TimeStepFine                 ! fine time step [s]
    real(kind=kind_noahmp)            :: SoilSatExcAcc                ! accumulation of soil saturation excess [m]
    real(kind=kind_noahmp)            :: SoilWatConductAcc            ! sum of SoilWatConductivity*ThicknessSnowSoilLayer
    real(kind=kind_noahmp)            :: WaterRemove                  ! water mass removal [mm]
    real(kind=kind_noahmp)            :: SoilWatRem                   ! temporary remaining soil water [mm]
    real(kind=kind_noahmp)            :: SoilWaterMin                 ! minimum soil water [mm]
    real(kind=kind_noahmp)            :: DrainSoilBotAcc              ! accumulated drainage water [mm] at fine time step
    real(kind=kind_noahmp)            :: RunoffSurfaceAcc             ! accumulated surface runoff [mm] at fine time step
    real(kind=kind_noahmp)            :: InfilSfcAcc                  ! accumulated infiltration rate [m/s]
    real(kind=kind_noahmp)            :: WaterTableDepthBegin         ! WTD at start of timestep [m]
    real(kind=kind_noahmp)            :: Q_net_soil                   ! net water flux into soil [m/s]
    real(kind=kind_noahmp)            :: Sy_local                     ! specific yield [-]
    real(kind=kind_noahmp)            :: d_top_k, d_bot_k             ! layer depth bounds [m]
    real(kind=kind_noahmp)            :: thetas_loc, ae_loc, bb_loc   ! soil params
    real(kind=kind_noahmp)            :: threshold_k                  ! per-layer disequilib threshold
    logical                           :: AnyDisequilibrium            ! any layer in disequilibrium?
    real(kind=kind_noahmp), parameter :: WaterTableDepthMinPeat = -0.2449_kind_noahmp
    real(kind=kind_noahmp), parameter :: DisequilMargin = 0.2_kind_noahmp
    real(kind=kind_noahmp), parameter :: HysteresisHalf = 0.05_kind_noahmp
    real(kind=kind_noahmp), parameter :: SoilImpPara = 4.0            ! soil impervious fraction parameter
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft1     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft2     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft3     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqTmp   ! temporary soil liquid water [mm]
    real(kind=kind_noahmp), allocatable, dimension(:) :: theta_eq_micro ! microtopo equilibrium moisture
    real(kind=kind_noahmp), allocatable, dimension(:) :: theta_eq_1D   ! flat-column equilibrium moisture

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    noahmp soil time step [s]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth of layer bottom [m]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              FlagUrban              => noahmp%config%domain%FlagUrban              ,& ! in,    logical flag for urban grid
              OptRunoffSurface       => noahmp%config%nmlist%OptRunoffSurface       ,& ! in,    options for surface runoff
              OptRunoffSubsurface    => noahmp%config%nmlist%OptRunoffSubsurface    ,& ! in,    options for subsurface runoff
              OptPeatlandPhysics     => noahmp%config%nmlist%OptPeatlandPhysics     ,& ! in,    options for peatland physics
              OptTileDrainage        => noahmp%config%nmlist%OptTileDrainage        ,& ! in,    options for tile drainage
              SoilIce                => noahmp%water%state%SoilIce                  ,& ! in,    soil ice content [m3/m3]
              TileDrainFrac          => noahmp%water%state%TileDrainFrac            ,& ! in,    tile drainage map (fraction)
              SoilSfcInflowMean      => noahmp%water%flux%SoilSfcInflowMean         ,& ! in,    mean water input on soil surface [m/s]
              SoilMoistureSat        => noahmp%water%param%SoilMoistureSat          ,& ! in,    saturated value of soil moisture [m3/m3]
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in,    saturated matric potential [m]
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in,    Campbell exponent B
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil water content [m3/m3]
              RechargeGwDeepWT       => noahmp%water%state%RechargeGwDeepWT         ,& ! inout, recharge to or from the water table when deep [m]
              DrainSoilBot           => noahmp%water%flux%DrainSoilBot              ,& ! out,   soil bottom drainage [m/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface             ,& ! out,   surface runoff [mm per soil timestep]
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface          ,& ! out,   subsurface runoff [mm per soil timestep] 
              InfilRateSfc           => noahmp%water%flux%InfilRateSfc              ,& ! out,   infiltration rate at surface [m/s]
              TileDrain              => noahmp%water%flux%TileDrain                 ,& ! out,   tile drainage [mm per soil timestep]
              Transpiration          => noahmp%water%flux%Transpiration             ,& ! in,    transpiration rate [mm/s]
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet             ,& ! in,    net ground evaporation [mm/s]
              EvapSoilSfcLiqMean     => noahmp%water%flux%EvapSoilSfcLiqMean       ,& ! in,    mean soil surface evaporation [m/s]
              TranspWatLossSoilMean  => noahmp%water%flux%TranspWatLossSoilMean     ,& ! in,    mean transpiration per layer [m/s]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! inout, water table depth [m]
              SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax        ,& ! out,   maximum soil imperviousness fraction
              SoilWatConductivity    => noahmp%water%state%SoilWatConductivity      ,& ! out,   soil hydraulic conductivity [m/s]
              SoilEffPorosity        => noahmp%water%state%SoilEffPorosity          ,& ! out,   soil effective porosity [m3/m3]
              SoilImpervFrac         => noahmp%water%state%SoilImpervFrac           ,& ! out,   impervious fraction due to frozen soil
              SoilIceFrac            => noahmp%water%state%SoilIceFrac              ,& ! out,   ice fraction in frozen soil
              SoilSaturationExcess   => noahmp%water%state%SoilSaturationExcess     ,& ! out,   saturation excess of the total soil [m]
              SoilIceMax             => noahmp%water%state%SoilIceMax               ,& ! out,   maximum soil ice content [m3/m3]
              FSW_change             => noahmp%water%state%FSW_change               ,& ! inout, surface storage change [mm]
              FloodedFrac            => noahmp%water%state%FloodedFrac              ,& ! inout, flooded fraction [-]
              f_part                 => noahmp%water%state%f_part                   ,& ! inout, soil flux partition fraction [-]
              f_soil_k               => noahmp%water%state%f_soil_k                ,& ! inout, per-layer soil fraction
              InDisequilibrium       => noahmp%water%state%InDisequilibrium         ,& ! inout, per-layer disequilibrium flag
              SoilMoistureDeficit    => noahmp%water%state%SoilMoistureDeficit      ,& ! inout, per-layer deficit
              SoilLiqWaterMin        => noahmp%water%state%SoilLiqWaterMin          & ! out,   minimum soil liquid water content [m3/m3]
             )
! ----------------------------------------------------------------------

    ! initialization
    if (.not. allocated(MatRight)  ) allocate(MatRight  (1:NumSoilLayer))
    if (.not. allocated(MatLeft1)  ) allocate(MatLeft1  (1:NumSoilLayer))
    if (.not. allocated(MatLeft2)  ) allocate(MatLeft2  (1:NumSoilLayer))
    if (.not. allocated(MatLeft3)  ) allocate(MatLeft3  (1:NumSoilLayer))
    if (.not. allocated(SoilLiqTmp)) allocate(SoilLiqTmp(1:NumSoilLayer))
    MatRight         = 0.0
    MatLeft1         = 0.0
    MatLeft2         = 0.0
    MatLeft3         = 0.0
    SoilLiqTmp       = 0.0
    RunoffSurface    = 0.0
    RunoffSubsurface = 0.0
    InfilRateSfc     = 0.0
    SoilSatExcAcc    = 0.0
    InfilSfcAcc      = 1.0e-06

    if ( OptPeatlandPhysics == 1 ) then
       call ApplyPeatlandPhysics(noahmp)
    endif

    ! for the case when snowmelt water is too large
    do LoopInd1 = 1, NumSoilLayer
       SoilEffPorosity(LoopInd1) = max(1.0e-4, (SoilMoistureSat(LoopInd1) - SoilIce(LoopInd1)))
       SoilSatExcAcc             = SoilSatExcAcc + max(0.0, SoilLiqWater(LoopInd1) - SoilEffPorosity(LoopInd1)) * &
                                                   ThicknessSnowSoilLayer(LoopInd1)
       SoilLiqWater(LoopInd1)    = min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1))
    enddo

    ! impermeable fraction due to frozen soil
    do LoopInd1 = 1, NumSoilLayer
       SoilIceFrac(LoopInd1)    = min(1.0, SoilIce(LoopInd1) / SoilMoistureSat(LoopInd1))
       SoilImpervFrac(LoopInd1) = max(0.0, exp(-SoilImpPara*(1.0-SoilIceFrac(LoopInd1))) - exp(-SoilImpPara)) / &
                                  (1.0 - exp(-SoilImpPara))
    enddo

    ! maximum soil ice content and minimum liquid water of all layers
    SoilIceMax        = 0.0
    SoilImpervFracMax = 0.0
    SoilLiqWaterMin   = SoilMoistureSat(1)
    do LoopInd1 = 1, NumSoilLayer
       if ( SoilIce(LoopInd1) > SoilIceMax )               SoilIceMax        = SoilIce(LoopInd1)
       if ( SoilImpervFrac(LoopInd1) > SoilImpervFracMax ) SoilImpervFracMax = SoilImpervFrac(LoopInd1)
       if ( SoilLiqWater(LoopInd1) < SoilLiqWaterMin )     SoilLiqWaterMin   = SoilLiqWater(LoopInd1)
    enddo

    ! subsurface runoff for runoff scheme option 2
    if ( OptRunoffSubsurface == 2 ) call RunoffSubSurfaceEquiWaterTable(noahmp)

    ! Peatland: diagnose WTD and compute subsurface runoff
    if ( OptRunoffSubsurface == 9 ) then
       call RunoffSubSurfacePeatland(noahmp)
    endif

    ! jref impermable surface at urban
    if ( FlagUrban .eqv. .true. ) SoilImpervFrac(1) = 0.95

    WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)

    ! surface runoff and infiltration rate using different schemes
    if ( OptRunoffSurface == 1 ) call RunoffSurfaceTopModelGrd(noahmp)
    if ( OptRunoffSurface == 2 ) call RunoffSurfaceTopModelEqui(noahmp)
    if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 4 ) call RunoffSurfaceBATS(noahmp)
    if ( OptRunoffSurface == 5 ) call RunoffSurfaceTopModelMMF(noahmp)
    if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,SoilTimeStep,InfilSfcAcc)

    ! For peatland: add surface runoff to infiltration (no ponding)
    if ( OptPeatlandPhysics == 1 ) then
       InfilRateSfc = InfilRateSfc + RunoffSurface
       RunoffSurface = 0.0
    endif

    ! =====================================================================
    !  PEATLAND TWO-REGIME APPROACH
    ! =====================================================================
    peatland_block: if ( OptPeatlandPhysics == 1 ) then

       ! Allocate temporaries
       if (.not. allocated(theta_eq_micro)) allocate(theta_eq_micro(1:NumSoilLayer))
       if (.not. allocated(theta_eq_1D))    allocate(theta_eq_1D(1:NumSoilLayer))

       ! Diagnose WTD from current soil moisture
       call WaterTableEquilibriumPeat(noahmp)
       WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)
       WaterTableDepthBegin = WaterTableDepth

       ! Compute microtopography variables (FloodedFrac, f_part, f_soil_k)
       call MicroTopoCorrection(noahmp)

       ! Split infiltration: soil vs surface water
       InfilRateSfc = f_part * InfilRateSfc   ! soil portion (already overwritten)

       ! Soil hydraulic parameters (single type from layer 1)
       thetas_loc = SoilMoistureSat(1)
       ae_loc     = abs(SoilMatPotentialSat(1))
       bb_loc     = SoilExpCoeffB(1)

       ! Determine per-layer disequilibrium (with hysteresis)
       AnyDisequilibrium = .false.
       do LoopInd1 = 1, NumSoilLayer
          d_bot_k = abs(DepthSoilLayer(LoopInd1))
          threshold_k = d_bot_k + DisequilMargin
          if (InDisequilibrium(LoopInd1)) then
             ! Currently in disequilibrium: exit if WTD drops below threshold
             if (WaterTableDepth < threshold_k - HysteresisHalf) then
                InDisequilibrium(LoopInd1) = .false.
                SoilMoistureDeficit(LoopInd1) = 0.0_kind_noahmp
             endif
          else
             ! Currently in equilibrium: enter disequilibrium if WTD rises above threshold
             if (WaterTableDepth > threshold_k + HysteresisHalf) then
                InDisequilibrium(LoopInd1) = .true.
             endif
          endif
          if (InDisequilibrium(LoopInd1)) AnyDisequilibrium = .true.
       enddo

       ! Compute equilibrium moisture profiles at current WTD
       do LoopInd1 = 1, NumSoilLayer
          if (LoopInd1 == 1) then
             d_top_k = 0.0_kind_noahmp
          else
             d_top_k = abs(DepthSoilLayer(LoopInd1-1))
          endif
          d_bot_k = abs(DepthSoilLayer(LoopInd1))
          theta_eq_micro(LoopInd1) = MicroTopoEquilibriumMoisture(WaterTableDepth, &
                                     d_top_k, d_bot_k, thetas_loc, ae_loc, bb_loc)
          theta_eq_1D(LoopInd1) = MicroTopoEquilibriumMoisture1D(WaterTableDepth, &
                                  d_top_k, d_bot_k, thetas_loc, ae_loc, bb_loc)
       enddo

       regime_select: if (.not. AnyDisequilibrium) then
          ! =============================================================
          !  REGIME 1: Full Equilibrium — bypass Richards
          ! =============================================================
          ! Net water flux into the soil [m/s]
          Q_net_soil = InfilRateSfc
          do LoopInd1 = 1, NumSoilLayer
             Q_net_soil = Q_net_soil - f_soil_k(LoopInd1) * TranspWatLossSoilMean(LoopInd1)
          enddo
          Q_net_soil = Q_net_soil - f_soil_k(1) * EvapSoilSfcLiqMean
          Q_net_soil = Q_net_soil - RunoffSubsurface / 1000.0_kind_noahmp

          ! Specific yield at current WTD
          Sy_local = MicroTopoSpecificYield(WaterTableDepth, NumSoilLayer, &
                     DepthSoilLayer, thetas_loc, ae_loc, bb_loc)

          ! Update WTD: soil gains water -> WTD decreases (water table rises)
          WaterTableDepth = WaterTableDepth - Q_net_soil * SoilTimeStep / Sy_local
          WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)

          ! Set soil moisture from equilibrium at new WTD
          do LoopInd1 = 1, NumSoilLayer
             if (LoopInd1 == 1) then
                d_top_k = 0.0_kind_noahmp
             else
                d_top_k = abs(DepthSoilLayer(LoopInd1-1))
             endif
             d_bot_k = abs(DepthSoilLayer(LoopInd1))
             SoilLiqWater(LoopInd1) = MicroTopoEquilibriumMoisture(WaterTableDepth, &
                                      d_top_k, d_bot_k, thetas_loc, ae_loc, bb_loc)
             SoilLiqWater(LoopInd1) = max(1.0e-4_kind_noahmp, &
                                     min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1)))
             SoilMoistureDeficit(LoopInd1) = 0.0_kind_noahmp
          enddo

          DrainSoilBot = 0.0

       else regime_select
          ! =============================================================
          !  REGIME 2: Deficit-Based Richards
          ! =============================================================
          ! Compute deficits and set up 1D profile for Richards
          do LoopInd1 = 1, NumSoilLayer
             SoilMoistureDeficit(LoopInd1) = theta_eq_micro(LoopInd1) - SoilLiqWater(LoopInd1)
             ! 1D profile: flat-column equilibrium minus deficit
             SoilLiqWater(LoopInd1) = theta_eq_1D(LoopInd1) - SoilMoistureDeficit(LoopInd1)
             SoilLiqWater(LoopInd1) = max(1.0e-4_kind_noahmp, &
                                     min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1)))
          enddo
          SoilMoisture(1:NumSoilLayer) = SoilLiqWater(1:NumSoilLayer) + SoilIce(1:NumSoilLayer)

          ! Iteration loop for Richards
          NumIterSoilWat = 3
          if ( (InfilRateSfc*SoilTimeStep) > (ThicknessSnowSoilLayer(1)*SoilMoistureSat(1)) ) then
             NumIterSoilWat = NumIterSoilWat*2
          endif
          TimeStepFine = SoilTimeStep / NumIterSoilWat

          DrainSoilBotAcc  = 0.0
          do IndIter = 1, NumIterSoilWat
             call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
             ! Pin equilibrium layers: Dirichlet boundary (no change)
             do LoopInd1 = 1, NumSoilLayer
                if (.not. InDisequilibrium(LoopInd1)) then
                   MatLeft1(LoopInd1) = 0.0_kind_noahmp
                   MatLeft2(LoopInd1) = 1.0_kind_noahmp
                   MatLeft3(LoopInd1) = 0.0_kind_noahmp
                   MatRight(LoopInd1) = 0.0_kind_noahmp
                endif
             enddo
             call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)
             SoilSatExcAcc   = SoilSatExcAcc + SoilSaturationExcess
             DrainSoilBotAcc = DrainSoilBotAcc + DrainSoilBot
          enddo
          DrainSoilBot = DrainSoilBotAcc / NumIterSoilWat

          ! Remove subsurface runoff from the 1D profile (proportional to K*dz)
          SoilWatConductAcc = 0.0_kind_noahmp
          do LoopInd1 = 1, NumSoilLayer
             SoilWatConductAcc = SoilWatConductAcc + &
                                 SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)
          enddo
          if (SoilWatConductAcc > 0.0_kind_noahmp) then
             do LoopInd1 = 1, NumSoilLayer
                WaterRemove = RunoffSubsurface * SoilTimeStep * &
                              (SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc
                SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1) * 1000.0)
             enddo
          endif

          ! Convert back from 1D profile to area-averaged moisture
          do LoopInd1 = 1, NumSoilLayer
             if (InDisequilibrium(LoopInd1)) then
                SoilMoistureDeficit(LoopInd1) = theta_eq_1D(LoopInd1) - SoilLiqWater(LoopInd1)
                SoilLiqWater(LoopInd1) = theta_eq_micro(LoopInd1) - SoilMoistureDeficit(LoopInd1)
             else
                SoilLiqWater(LoopInd1) = theta_eq_micro(LoopInd1)
             endif
             SoilLiqWater(LoopInd1) = max(1.0e-4_kind_noahmp, &
                                     min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1)))
          enddo
       endif regime_select

       ! Update SoilMoisture
       do LoopInd1 = 1, NumSoilLayer
          SoilMoisture(LoopInd1) = SoilLiqWater(LoopInd1) + SoilIce(LoopInd1)
       enddo

       ! Diagnose WTD from updated moisture
       call WaterTableEquilibriumPeat(noahmp)
       WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)

       ! Definitional FSW_change from surface water storage function
       FSW_change = SurfaceWaterStorage_mm(WaterTableDepth) - SurfaceWaterStorage_mm(WaterTableDepthBegin)

       ! DrainSoilBot already set (0 for Regime 1, accumulated for Regime 2)
       DrainSoilBot = DrainSoilBot * 1000.0  ! m/s -> mm/s
       RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep  ! mm/s

       ! Deallocate peatland temporaries
       if (allocated(theta_eq_micro)) deallocate(theta_eq_micro)
       if (allocated(theta_eq_1D))    deallocate(theta_eq_1D)

    else peatland_block
    ! =====================================================================
    !  STANDARD (NON-PEATLAND) FLOW
    ! =====================================================================

    ! determine iteration times
    NumIterSoilWat = 3
    if ( (InfilRateSfc*SoilTimeStep) > (ThicknessSnowSoilLayer(1)*SoilMoistureSat(1)) ) then
       NumIterSoilWat = NumIterSoilWat*2
    endif
    TimeStepFine = SoilTimeStep / NumIterSoilWat

    ! solve soil moisture
    InfilSfcAcc      = 1.0e-06
    DrainSoilBotAcc  = 0.0
    RunoffSurfaceAcc = 0.0

    do IndIter = 1, NumIterSoilWat
       if ( SoilSfcInflowMean > 0.0 ) then
          if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,TimeStepFine,InfilSfcAcc)
       endif
       call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
       call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)
       SoilSatExcAcc    = SoilSatExcAcc + SoilSaturationExcess
       DrainSoilBotAcc  = DrainSoilBotAcc + DrainSoilBot
       RunoffSurfaceAcc = RunoffSurfaceAcc + RunoffSurface
    enddo

    DrainSoilBot  = DrainSoilBotAcc / NumIterSoilWat
    RunoffSurface = RunoffSurfaceAcc / NumIterSoilWat
    RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep  ! m/s -> mm/s
    DrainSoilBot  = DrainSoilBot * 1000.0  ! m/s -> mm/s

    ! compute tile drainage
    if ( (OptTileDrainage == 1) .and. (TileDrainFrac > 0.3) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageSimple(noahmp)
    endif
    if ( (OptTileDrainage == 2) .and. (TileDrainFrac > 0.1) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageHooghoudt(noahmp)
    END IF

    ! removal of soil water due to subsurface runoff (option 2)
    if ( OptRunoffSubsurface == 2 ) then
       SoilWatConductAcc = 0.0
       do LoopInd1 = 1, NumSoilLayer
          SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)
       enddo
       do LoopInd1 = 1, NumSoilLayer
          WaterRemove            = RunoffSubsurface * SoilTimeStep * &
                                  (SoilWatConductivity(LoopInd1)*ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc
          SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1)*1000.0)
       enddo
    endif

    ! Limit SoilLiqTmp to be greater than or equal to watmin.
    if ( OptRunoffSubsurface /= 1 ) then
       do LoopInd2 = 1, NumSoilLayer
          SoilLiqTmp(LoopInd2) = SoilLiqWater(LoopInd2) * ThicknessSnowSoilLayer(LoopInd2) * 1000.0
       enddo

       SoilWaterMin = 0.01   ! mm
       do LoopInd2 = 1, NumSoilLayer-1
          if ( SoilLiqTmp(LoopInd2) < 0.0 ) then
             SoilWatRem = SoilWaterMin - SoilLiqTmp(LoopInd2)
          else
             SoilWatRem = 0.0
          endif
          SoilLiqTmp(LoopInd2  ) = SoilLiqTmp(LoopInd2  ) + SoilWatRem
          SoilLiqTmp(LoopInd2+1) = SoilLiqTmp(LoopInd2+1) - SoilWatRem
       enddo
       LoopInd2 = NumSoilLayer
       if ( SoilLiqTmp(LoopInd2) < SoilWaterMin ) then
           SoilWatRem = SoilWaterMin - SoilLiqTmp(LoopInd2)
       else
           SoilWatRem = 0.0
       endif
       SoilLiqTmp(LoopInd2) = SoilLiqTmp(LoopInd2) + SoilWatRem
       RunoffSubsurface     = RunoffSubsurface - SoilWatRem/SoilTimeStep

       if ( OptRunoffSubsurface == 5 ) RechargeGwDeepWT = RechargeGwDeepWT - SoilWatRem * 1.0e-3

       do LoopInd2 = 1, NumSoilLayer
          SoilLiqWater(LoopInd2) = SoilLiqTmp(LoopInd2) / (ThicknessSnowSoilLayer(LoopInd2)*1000.0)
       enddo
    endif

    ! compute groundwater and subsurface runoff
    if ( OptRunoffSubsurface == 1 ) call RunoffSubSurfaceGroundWater(noahmp)

    ! compute subsurface runoff based on drainage rate
    if ( (OptRunoffSubsurface == 3) .or. (OptRunoffSubsurface == 4) .or. (OptRunoffSubsurface == 6) .or. &
         (OptRunoffSubsurface == 7) .or. (OptRunoffSubsurface == 8) ) then
         call RunoffSubSurfaceDrainage(noahmp)
    endif

    ! update soil moisture
    do LoopInd2 = 1, NumSoilLayer
        SoilMoisture(LoopInd2) = SoilLiqWater(LoopInd2) + SoilIce(LoopInd2)
    enddo

    ! compute subsurface runoff and shallow water table for MMF scheme
    if ( OptRunoffSubsurface == 5 ) call RunoffSubSurfaceShallowWaterMMF(noahmp)

    endif peatland_block

    ! accumulated water flux over soil timestep [mm]
    RunoffSurface    = RunoffSurface    * SoilTimeStep
    RunoffSubsurface = RunoffSubsurface * SoilTimeStep
    TileDrain        = TileDrain        * SoilTimeStep

    ! deallocate local arrays
    deallocate(MatRight  )
    deallocate(MatLeft1  )
    deallocate(MatLeft2  )
    deallocate(MatLeft3  )
    deallocate(SoilLiqTmp)

    end associate

  end subroutine SoilWaterMain

end module SoilWaterMainMod
