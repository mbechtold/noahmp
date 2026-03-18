module SoilWaterMainMod

!!! Main soil water module including all soil water processes & update soil moisture
!!! surface runoff, infiltration, soil water diffusion, subsurface runoff, tile drainage
!!! Peatland-specific options (Chakraborty et al., 2025; revised Bechtold, 2026)

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
  use MicroTopoCorrectionMod,            only : MicroTopoCorrection
  use PeatlandPhysicsMod,                only : ApplyPeatlandPhysics
  use PeatMicroTopoMod,                  only : FloodedFrac, FsoilMicroTopo, &
                                                 EquilibriumSMFlat,           &
                                                 EquilibriumSMMicroTopo,      &
                                                 FindWaterTable
  use SoilWaterDiffusionRichardsMod,     only : SoilWaterDiffusionRichards
  use SoilMoistureSolverMod,             only : SoilMoistureSolver
  use TileDrainageSimpleMod,             only : TileDrainageSimple
  use TileDrainageHooghoudtMod,          only : TileDrainageHooghoudt
  use WaterTableEquilibriumMod,          only : WaterTableEquilibrium

  implicit none

contains

  subroutine SoilWaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SOILWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Option=9; Peatland runoff scheme (Chakraborty & Bechtold, 2025; revised Bechtold, 2026)
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
    real(kind=kind_noahmp)            :: InfilRateSfc_FSW_change      ! surface water infiltration flux
    real(kind=kind_noahmp)            :: W_soil_peat                  ! total soil water content [m]
    real(kind=kind_noahmp)            :: z_wt_peat                    ! water table elevation (z positive up) [m]
    real(kind=kind_noahmp), parameter :: WaterTableDepthMinPeat = -1.0
    real(kind=kind_noahmp), parameter :: SoilImpPara = 4.0            ! soil impervious fraction parameter
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatRight     ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft1     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft2     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft3     ! left-hand side term
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqTmp   ! temporary soil liquid water [mm]
    real(kind=kind_noahmp)            :: d_top_peat, d_bot_peat       ! layer depth bounds for peat [m]
    real(kind=kind_noahmp)            :: z_col_bot_peat               ! column bottom depth [m]
    real(kind=kind_noahmp)            :: thetas_peat, ae_peat, bb_peat! Campbell peat parameters
    real(kind=kind_noahmp)            :: SM_eq_flat_tmp                ! temporary flat equilibrium SM
    real(kind=kind_noahmp)            :: SM_eq_micro_tmp               ! temporary microtopo equilibrium SM
    real(kind=kind_noahmp)            :: SM_old_peat                   ! SoilLiqWater before transfer [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilLiqGap    ! per-layer gap added in forward transfer [m3/m3]

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SoilTimeStep           => noahmp%config%domain%SoilTimeStep           ,& ! in,    noahmp soil time step [s]
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
              SoilMatPotentialSat    => noahmp%water%param%SoilMatPotentialSat      ,& ! in,    saturated soil matric potential [m]
              SoilExpCoeffB          => noahmp%water%param%SoilExpCoeffB            ,& ! in,    Campbell b exponent [-]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil water content [m3/m3]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! inout, total soil water content [m3/m3]
              RechargeGwDeepWT       => noahmp%water%state%RechargeGwDeepWT         ,& ! inout, recharge to or from the water table when deep [m]
              DrainSoilBot           => noahmp%water%flux%DrainSoilBot              ,& ! out,   soil bottom drainage [m/s]
              RunoffSurface          => noahmp%water%flux%RunoffSurface             ,& ! out,   surface runoff [mm per soil timestep]
              RunoffSubsurface       => noahmp%water%flux%RunoffSubsurface          ,& ! out,   subsurface runoff [mm per soil timestep] 
              InfilRateSfc           => noahmp%water%flux%InfilRateSfc              ,& ! out,   infiltration rate at surface [m/s]
              TileDrain              => noahmp%water%flux%TileDrain                 ,& ! out,   tile drainage [mm per soil timestep]
              Transpiration          => noahmp%water%flux%Transpiration             ,& ! in,    transpiration rate [mm/s]
              EvapGroundNet          => noahmp%water%flux%EvapGroundNet             ,& ! in,    net ground (soil/snow) evaporation [mm/s]
              WaterTableDepth        => noahmp%water%state%WaterTableDepth          ,& ! inout,   water table depth [m]
              SoilImpervFracMax      => noahmp%water%state%SoilImpervFracMax        ,& ! out,   maximum soil imperviousness fraction
              SoilWatConductivity    => noahmp%water%state%SoilWatConductivity      ,& ! out,   soil hydraulic conductivity [m/s]
              SoilEffPorosity        => noahmp%water%state%SoilEffPorosity          ,& ! out,   soil effective porosity [m3/m3]
              SoilImpervFrac         => noahmp%water%state%SoilImpervFrac           ,& ! out,   impervious fraction due to frozen soil
              SoilIceFrac            => noahmp%water%state%SoilIceFrac              ,& ! out,   ice fraction in frozen soil
              SoilSaturationExcess   => noahmp%water%state%SoilSaturationExcess     ,& ! out,   saturation excess of the total soil [m]
              SoilIceMax             => noahmp%water%state%SoilIceMax               ,& ! out,   maximum soil ice content [m3/m3]
              FSW_change             => noahmp%water%state%FSW_change               ,& ! inout,   surface storage change [mm]
              FloodedFraction        => noahmp%water%state%FloodedFraction          ,& ! inout,   flooded fraction [-]
              f_soil                 => noahmp%water%state%f_soil                   ,& ! inout, fraction of flux in and out of soil [-]
              SoilLiqWaterMin        => noahmp%water%state%SoilLiqWaterMin         ,& ! out,   minimum soil liquid water content [m3/m3]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer          & ! in,    depth [m] of layer-bottom from soil surface
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
    
    ! ================================================================
    ! Peatland: Diagnose WTD from soil water (microtopo-aware)
    ! Then compute Ivanov runoff using the diagnosed WTD.
    ! (Bechtold, 2026)
    ! ================================================================
    if ( OptPeatlandPhysics == 1 ) then
       thetas_peat    = SoilMoistureSat(1)
       ae_peat        = abs(SoilMatPotentialSat(1))
       bb_peat        = SoilExpCoeffB(1)
       z_col_bot_peat = abs(DepthSoilLayer(NumSoilLayer))

       ! Compute total soil water [m]
       W_soil_peat = 0.0_kind_noahmp
       do LoopInd1 = 1, NumSoilLayer
          W_soil_peat = W_soil_peat + &
              SoilLiqWater(LoopInd1) * abs(ThicknessSnowSoilLayer(LoopInd1))
       enddo

       ! Diagnose microtopo-aware water table from soil storage
       z_wt_peat = FindWaterTable(W_soil_peat, thetas_peat, ae_peat, &
           bb_peat, z_col_bot_peat, -WaterTableDepth)
       WaterTableDepth = max(-z_wt_peat, WaterTableDepthMinPeat)
    endif

    ! Peatland: Ivanov runoff (uses diagnosed WTD)
    if ( OptRunoffSubsurface == 9 ) then
            call RunoffSubSurfacePeatland(noahmp)
    endif


    ! jref impermable surface at urban
    if ( FlagUrban .eqv. .true. ) SoilImpervFrac(1) = 0.95
    
    WaterTableDepth = max(WaterTableDepth, WaterTableDepthMinPeat)

    ! ================================================================
    ! Peatland forward excess transfer (Bechtold, 2026)
    ! SoilLiqWater at this point = microtopo-integrated soil moisture.
    ! Compute departure from microtopo equilibrium and map to flat 1D
    ! profile for the Richards equation solver.
    ! ================================================================
    if ( OptPeatlandPhysics == 1 ) then
       if (.not. allocated(SoilLiqGap)) allocate(SoilLiqGap(1:NumSoilLayer))

       thetas_peat    = SoilMoistureSat(1)
       ae_peat        = abs(SoilMatPotentialSat(1))
       bb_peat        = SoilExpCoeffB(1)
       z_col_bot_peat = abs(DepthSoilLayer(NumSoilLayer))

       do LoopInd1 = 1, NumSoilLayer
          if (LoopInd1 == 1) then
             d_top_peat = 0.0_kind_noahmp
          else
             d_top_peat = abs(DepthSoilLayer(LoopInd1 - 1))
          endif
          d_bot_peat = abs(DepthSoilLayer(LoopInd1))

          SM_eq_micro_tmp = EquilibriumSMMicroTopo(d_top_peat, d_bot_peat, &
              WaterTableDepth, thetas_peat, ae_peat, bb_peat)
          SM_eq_flat_tmp  = EquilibriumSMFlat(d_top_peat, d_bot_peat, &
              WaterTableDepth, thetas_peat, ae_peat, bb_peat)

          ! Store the actual gap applied per layer for exact reversal
          SM_old_peat = SoilLiqWater(LoopInd1)
          SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) + (SM_eq_flat_tmp - SM_eq_micro_tmp)
          SoilLiqWater(LoopInd1) = max(0.001_kind_noahmp, &
              min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1)))
          ! Record actual change (accounts for clamping)
          SoilLiqGap(LoopInd1) = SoilLiqWater(LoopInd1) - SM_old_peat
       enddo
    endif

    ! surface runoff and infiltration rate using different schemes
    ! MB: Alternative idea which will produce the same output as PEATCLSM: 
    if ( OptRunoffSurface == 1 ) call RunoffSurfaceTopModelGrd(noahmp)
    if ( OptRunoffSurface == 2 ) call RunoffSurfaceTopModelEqui(noahmp)
    if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 4 ) call RunoffSurfaceBATS(noahmp)
    if ( OptRunoffSurface == 5 ) call RunoffSurfaceTopModelMMF(noahmp)
    if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,SoilTimeStep)
    if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,SoilTimeStep,InfilSfcAcc)
    
    ! MB: We add RunoffSurface to the water that needs to infiltrate:
    ! MB: ToDo, test with peat soil properties, currently leading to too high water balance errors
    if ( OptPeatlandPhysics == 1 ) then
       ! NOTE: from RunoffSurface routines, the RunoffSurface is still in m/s same as InfilRateSfc
       InfilRateSfc = InfilRateSfc + RunoffSurface
       RunoffSurface = 0.0
    endif

    ! determine iteration times  to solve soil water diffusion and moisture
    NumIterSoilWat = 3
    if ( (InfilRateSfc*SoilTimeStep) > (ThicknessSnowSoilLayer(1)*SoilMoistureSat(1)) ) then
       NumIterSoilWat = NumIterSoilWat*2
    endif
    TimeStepFine = SoilTimeStep / NumIterSoilWat

    ! solve soil moisture
    InfilSfcAcc      = 1.0e-06
    DrainSoilBotAcc  = 0.0
    RunoffSurfaceAcc = 0.0

    !if ( OptRunoffSubsurface == 9 ) then
    if ( OptPeatlandPhysics == 1 ) then
               call MicroTopoCorrection(noahmp) ! Compute f_soil and FloodedFraction
               ! Partition fluxes: f_soil fraction goes to soil Richards solver,
               ! (1-f_soil) fraction goes to surface water storage
               InfilRateSfc_FSW_change = (1.0_kind_noahmp - f_soil) * InfilRateSfc
               InfilRateSfc = f_soil * InfilRateSfc
    endif

    do IndIter = 1, NumIterSoilWat
       if ( SoilSfcInflowMean > 0.0 ) then
          if ( OptRunoffSurface == 3 ) call RunoffSurfaceFreeDrain(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 6 ) call RunoffSurfaceVIC(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 7 ) call RunoffSurfaceXinAnJiang(noahmp,TimeStepFine)
          if ( OptRunoffSurface == 8 ) call RunoffSurfaceDynamicVic(noahmp,TimeStepFine,InfilSfcAcc)
          ! MB: During iteration, again we add RunoffSurface to the water that needs to infiltrate:
          ! MB: ToDo, test with peat soil properties, currently leading to too high water balance errors
          !if ( OptRunoffSubsurface == 9 ) then
          !   InfilRateSfc = InfilRateSfc + RunoffSurface / SoilTimeStep / 1000.0 
          ! division because units between InfilRateSfc and RunoffSurface differ
          !   RunoffSurface = 0.0
          !endif
          !if ( OptRunoffSubsurface == 9 ) then
          !if ( OptPeatlandPhysics == 1 ) then
          !     InfilRateSfc_FSW_change = (1-f_soil)*InfilRateSfc
          !     InfilRateSfc = f_soil*InfilRateSfc
          !endif
       endif
       
       ! MB: Here the reduced InfilRateSfc will be redistributed as usual
       call SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight)
       call SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight)
       SoilSatExcAcc    = SoilSatExcAcc + SoilSaturationExcess
       DrainSoilBotAcc  = DrainSoilBotAcc + DrainSoilBot
       RunoffSurfaceAcc = RunoffSurfaceAcc + RunoffSurface
    enddo

    !if ( OptRunoffSubsurface == 9 ) then
    if ( OptPeatlandPhysics == 1 ) then
               ! Compute surface water storage change using flux partitioning
               ! (1-f_soil) fraction of each flux goes to/from surface water
               FSW_change = 0.0
               ! Surface water receives (1-f_soil) fraction of infiltration
               FSW_change = FSW_change + InfilRateSfc_FSW_change * SoilTimeStep * 1000.0_kind_noahmp
               ! Surface water loses (1-f_soil) fraction of evaporation
               FSW_change = FSW_change - (1.0_kind_noahmp - f_soil) * EvapGroundNet * SoilTimeStep
               ! Surface water loses (1-f_soil) fraction of transpiration  
               FSW_change = FSW_change - (1.0_kind_noahmp - f_soil) * Transpiration * SoilTimeStep
    endif

    DrainSoilBot  = DrainSoilBotAcc / NumIterSoilWat
    RunoffSurface = RunoffSurfaceAcc / NumIterSoilWat
    RunoffSurface = RunoffSurface * 1000.0 + SoilSatExcAcc * 1000.0 / SoilTimeStep  ! m/s -> mm/s
    DrainSoilBot  = DrainSoilBot * 1000.0  ! m/s -> mm/s

    ! compute tile drainage ! pvk
    if ( (OptTileDrainage == 1) .and. (TileDrainFrac > 0.3) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageSimple(noahmp)  ! simple tile drainage
    endif
    if ( (OptTileDrainage == 2) .and. (TileDrainFrac > 0.1) .and. (OptRunoffSurface == 3) ) then
       call TileDrainageHooghoudt(noahmp)  ! Hooghoudt tile drain
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

    ! Peatland: remove runoff from soil and surface water (Bechtold, 2026)
    !if ( OptRunoffSubsurface == 9 ) then
    if ( OptPeatlandPhysics == 1 ) then
       ! Subsurface runoff is removed proportionally from soil layers
       ! and from surface water according to f_soil / (1-f_soil) split
       SoilWatConductAcc = 0.0
       do LoopInd1 = 1, NumSoilLayer
          SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)
       enddo
       if (SoilWatConductAcc > 0.0) then
          do LoopInd1 = 1, NumSoilLayer
             WaterRemove = RunoffSubsurface * SoilTimeStep * &
                          (SoilWatConductivity(LoopInd1)*ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc
             ! Surface water portion of runoff removal
             FSW_change = FSW_change - (1.0_kind_noahmp - f_soil) * WaterRemove
             ! Soil portion
             WaterRemove = f_soil * WaterRemove
             SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1)*1000.0)
          enddo
       endif
    endif

    ! Limit SoilLiqTmp to be greater than or equal to watmin.
    ! Get water needed to bring SoilLiqTmp equal SoilWaterMin from lower layer.
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
    endif ! OptRunoffSubsurface /= 1

    ! compute groundwater and subsurface runoff
    if ( OptRunoffSubsurface == 1 ) call RunoffSubSurfaceGroundWater(noahmp)

    ! compute subsurface runoff based on drainage rate
    if ( (OptRunoffSubsurface == 3) .or. (OptRunoffSubsurface == 4) .or. (OptRunoffSubsurface == 6) .or. &
         (OptRunoffSubsurface == 7) .or. (OptRunoffSubsurface == 8) ) then
         call RunoffSubSurfaceDrainage(noahmp)
    endif
    
    ! ================================================================
    ! Peatland: backward excess transfer + final WTD (Bechtold, 2026)
    !
    ! After Richards, SoilLiqWater is a flat 1D profile.
    ! 1. Subtract the forward gap to restore microtopo-integrated values.
    ! 2. Diagnose final WTD from soil water via FindWaterTable.
    ! ================================================================
    if ( OptPeatlandPhysics == 1 ) then
        ! --- 1. Backward excess transfer: flat 1D → microtopo ---
        ! Subtract the SAME per-layer gap that was added in the forward
        ! transfer.  This is inherently water-conserving: the round-trip
        ! adds and removes exactly the same amount per layer.
        do LoopInd2 = 1, NumSoilLayer
           SoilLiqWater(LoopInd2) = SoilLiqWater(LoopInd2) - SoilLiqGap(LoopInd2)
           SoilLiqWater(LoopInd2) = max(0.001_kind_noahmp, &
               min(SoilEffPorosity(LoopInd2), SoilLiqWater(LoopInd2)))
        enddo

        ! --- 2. Diagnose final WTD (microtopo-aware) ---
        W_soil_peat = 0.0_kind_noahmp
        do LoopInd2 = 1, NumSoilLayer
           W_soil_peat = W_soil_peat + &
               SoilLiqWater(LoopInd2) * abs(ThicknessSnowSoilLayer(LoopInd2))
        enddo
        z_wt_peat = FindWaterTable(W_soil_peat, thetas_peat, ae_peat, &
            bb_peat, z_col_bot_peat, z_wt_peat)
        WaterTableDepth = max(-z_wt_peat, WaterTableDepthMinPeat)

        ! Update FloodedFraction from final WTD
        FloodedFraction = FloodedFrac(-WaterTableDepth)
    endif

    ! update soil moisture
    do LoopInd2 = 1, NumSoilLayer
        SoilMoisture(LoopInd2) = SoilLiqWater(LoopInd2) + SoilIce(LoopInd2)
    enddo

    ! compute subsurface runoff and shallow water table for MMF scheme
    if ( OptRunoffSubsurface == 5 ) call RunoffSubSurfaceShallowWaterMMF(noahmp)

    ! accumulated water flux over soil timestep [mm]
    RunoffSurface    = RunoffSurface    * SoilTimeStep
    RunoffSubsurface = RunoffSubsurface * SoilTimeStep
    TileDrain        = TileDrain        * SoilTimeStep

    ! deallocate local arrays to avoid memory leaks
    deallocate(MatRight  )
    deallocate(MatLeft1  )
    deallocate(MatLeft2  )
    deallocate(MatLeft3  )
    deallocate(SoilLiqTmp)
    if (allocated(SoilLiqGap)) deallocate(SoilLiqGap)

    end associate

  end subroutine SoilWaterMain

end module SoilWaterMainMod
